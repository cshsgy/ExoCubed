// athena
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// Parameters
Real om_earth = 1.75E-4;
Real phi0, ros_fac, rand_lim;


std::string dfile(std::string fname) {
  std::stringstream msg;

  std::ifstream file(fname.c_str(), std::ios::in);
  std::string ss;
  char c;
  while (file) {
    file.get(c);
    if (c == '#') {
      while (c != '\n' && file) file.get(c);
      continue;
    }
    ss += c;
  }
  return ss;
}

void read_data(AthenaArray<Real> *data, std::string fname, char c=char(32)) {
  // remove comment
  std::string str_file = dfile(fname);

  // read first time to determine dimension
  std::stringstream inp(str_file);
  // std::ifstream inp(fname.c_str(), std::ios::in);
  std::string line;
  std::getline(inp, line);
  int rows = 0, cols = 0;
  if (!line.empty()) {
    rows = 1;
    cols = line[0] == c ? 0 : 1;
    for (int i = 1; i < line.length(); ++i)
      if (line[i - 1] == c && line[i] != c) cols++;
  }
  while (std::getline(inp, line)) ++rows;
  rows--;
  // inp.close();

  // str_file = DecommentFile(fname);

  // read second time
  data->NewAthenaArray(rows, cols);
  std::stringstream inp2(str_file);
  // inp.open(fname.c_str(), std::ios::in);

  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) inp2 >> (*data)(i, j);
}


void Forcing2(MeshBlock *pmb, Real const time, Real const dt,
              AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
              AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
              AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        // coriolis force
        Real f = 2. * om_earth * sin(lat);
        Real U, V;
        pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);
        Real ll_acc_U = f * V;
        Real ll_acc_V = -f * U;
        Real acc1, acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        u(IM2, k, j, i) += dt * w(IDN, k, j, i) * acc2;
        u(IM3, k, j, i) += dt * w(IDN, k, j, i) * acc3;
      }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Application::Logger app("main");
  app->Log("ProblemGenerator: Cubed Sphere");

  AthenaArray<Real> vpos;
  read_data(&vpos, pin->GetString("problem", "vfile"));
  // By default, lat goes from -90 to 90, lon goes from 0 to 360, both linearly
  int lon_num = vpos.GetDim1();  // input file allows positions
  int lat_num = vpos.GetDim2();  // input file allows positions

  auto pexo3 = pimpl->pexo3;


  // setup initial height and velocity field
  for (int k = ks-NGHOST; k <= ke+NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        // Make sure lon is in [0, 2pi] and lat is in [0, pi]
        lon = lon < 0 ? lon + 2 * PI : lon;
        lat = lat + PI/2.0;
        // Find the left points along lon and lat in the input file
        Real dlat = PI / (lat_num - 1);
        Real dlon = 2 * PI / (lon_num - 1);
        int i1 = (int) (lon / dlon);
        int i2 = (int) (lat / dlat);

        // Interpolate the data
        Real lat1 = i2 * dlat;
        Real lat2 = (i2 + 1) * dlat;
        Real lon1 = i1 * dlon;
        Real lon2 = (i1 + 1) * dlon;
        Real data11 = vpos(i1, i2);
        Real data12 = vpos(i1, i2 + 1);
        Real data21 = vpos(i1 + 1, i2);
        Real data22 = vpos(i1 + 1, i2 + 1);
        Real interpolated_data = (data11 * (lat2 - lat) * (lon2 - lon) +
                                  data12 * (lat - lat1) * (lon2 - lon) +
                                  data21 * (lat2 - lat) * (lon - lon1) +
                                  data22 * (lat - lat1) * (lon - lon1)) /
                                 ((lat2 - lat1) * (lon2 - lon1));
        // Calculate the velocities
        phydro->w(IDN,k,j,i) = phi0;
        // Random noise
        float tmp = 1 - 2 * ((float) rand()/RAND_MAX - 0.5) * rand_lim;
        phydro->w(IDN,k,j,i) += interpolated_data * ros_fac * tmp;
      }

  // convert to conserved variables
  app->Log("Done with problem generator");
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);

}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // Calculate lat, lon, U, V for output and visualization
  for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        Real x = tan(pcoord->x2v(j));
        Real y = tan(pcoord->x3v(k));
        Real C = sqrt(1.0 + x * x);
        Real D = sqrt(1.0 + y * y);
        Real delta = 1.0 / (1.0 + x * x + y * y);

        Real Vy = phydro->w(IVY, k, j, i);
        Real Vz = phydro->w(IVZ, k, j, i);
        Real lat, lon;
        pimpl->pexo3->GetLatLon(&lat, &lon, k, j, i);
        user_out_var(0, k, j, i) = lat / PI * 180.0;
        user_out_var(1, k, j, i) = lon / PI * 180.0;
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1 * x1 + x2 * x2);
        Real f = 0.0;
        Real U, V;
        pimpl->pexo3->GetUV(&U, &V, Vy, Vz, k, j, i);
        user_out_var(2, k, j, i) = U;
        user_out_var(3, k, j, i) = V;
        user_out_var(4, k, j, i) = delta / (C * D);
      }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
  SetUserOutputVariableName(4, "sqrtg");
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Read parameters
  phi0 = pin->GetReal("problem", "phi0");
  ros_fac = pin->GetReal("problem", "ros_fac");
  rand_lim = pin->GetReal("problem", "rand_lim");
  EnrollUserExplicitSourceFunction(Forcing2);
}
