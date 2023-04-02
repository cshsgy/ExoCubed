#include <parameter_input.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>
#include <configure.hpp>
#include <random>

#include <cubed_sphere.hpp>
#include "../src/eos/eos.hpp"
#include "../src/field/field.hpp"

Real phi0, dphi, interval, tseed;
Real omega, vphi, vrad, polarity, skewness;
Real vis;

std::default_random_engine gen;
std::uniform_real_distribution<double> uniform;


void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Input field is a uniform 10 m/s zonal wind
  Real V = 0.0;
  Real U = 0.0;

  std::cout << "======Start of Problem generator======" << loc.lx2 << "||" << loc.lx3 <<std::endl;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN,k,j,i) = phi0; // / (R*R);
        Real Vy, Vz;
        GetVyVz(&Vy, &Vz, pcoord, U, V, k, j, i);
        phydro->w(IVY,k,j,i) = Vy;
        phydro->w(IVZ,k,j,i) = Vz;
      }
  std::cout << "Done with problem generator..." << std::endl;
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  // Calculate lat, lon, U, V for output and visualization
  for (int k = ks-NGHOST; k <= ke+NGHOST; ++k)
    for (int j = js-NGHOST; j <= je+NGHOST; ++j)
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
        Real Vy = phydro->w(IVY,k,j,i);
        Real Vz = phydro->w(IVZ,k,j,i);
        Real lat, lon;
        GetLatLon(&lat, &lon, pcoord, k, j, i);
        user_out_var(0,k,j,i) = lat/PI*180.0;
        user_out_var(1,k,j,i) = lon/PI*180.0;
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1*x1 + x2*x2);
        Real f = omega * sin(lat);
        Real U, V, U1, V1;
        GetUV(&U, &V, pcoord, Vy, Vz, k, j, i);
        user_out_var(2,k,j,i) = U;
        user_out_var(3,k,j,i) = V;
        // Calculate the vorticity
        Vy = phydro->w(IVY,k+1,j,i);
        Vz = phydro->w(IVZ,k+1,j,i);
        GetUV(&U, &V, pcoord, Vy, Vz, k+1, j, i);
        Vy = phydro->w(IVY,k-1,j,i);
        Vz = phydro->w(IVZ,k-1,j,i);
        GetUV(&U1, &V1, pcoord, Vy, Vz, k-1, j, i);
        Real dz = pcoord->dx3v(k) + pcoord->dx3v(k-1);
        Real dUdz = (U - U1)/dz;
        Real dVdz = (V - V1)/dz;

        Vy = phydro->w(IVY,k,j+1,i);
        Vz = phydro->w(IVZ,k,j+1,i);
        GetUV(&U, &V, pcoord, Vy, Vz, k, j+1, i);
        Vy = phydro->w(IVY,k,j-1,i);
        Vz = phydro->w(IVZ,k,j-1,i);
        GetUV(&U1, &V1, pcoord, Vy, Vz, k, j-1, i);
        Real dy = pcoord->dx2v(j) + pcoord->dx2v(j-1);
        Real dUdy = (U - U1)/dy;
        Real dVdy = (V - V1)/dy;

        // Decompose the gradients
        Real dUdlat, dUdlon, dVdlat, dVdlon;
        GetUV(&dUdlat, &dUdlon, pcoord, dUdy, dUdz, k, j, i);
        GetUV(&dVdlat, &dVdlon, pcoord, dVdy, dVdz, k, j, i);

        user_out_var(4,k,j,i) = dVdlon - dUdlat;
        user_out_var(5,k,j,i) = (user_out_var(4,k,j,i) + f) / phydro->w(IDN,k,j,i);
      }
}

void MeshBlock::UserWorkInLoop()
{
  Real time = pmy_mesh->time;
  Real dt = pmy_mesh->dt;
    // check if need to seed a new vortex
    if (time > tseed) {
      Real x1, x2, range;
      Real min_lat = 0.0;
      Real vpol = uniform(gen) > polarity ? 1 : -1;

      if (vpol < 0) // Anticyclone
        vpol *= skewness;
      Real vlat = -PI/2.0+uniform(gen)*PI;
      Real vlon = uniform(gen)*2.*M_PI;
    for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
      for (int j = js - NGHOST; j <= je + NGHOST; ++j)
        for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
          Real lat, lon, vel;
          GetLatLon(&lat, &lon, pcoord, k, j, i);
          Real radius = pcoord -> x1v(i); // Earth's radius in km
          Real dlat = (lat - vlat) * M_PI / 180.0;
          Real dlon = (lon - vlon) * M_PI / 180.0;
          Real a = sin(dlat / 2.0) * sin(dlat / 2.0) +
                cos(vlat * M_PI / 180.0) * cos(lat * M_PI / 180.0) *
                sin(dlon / 2.0) * sin(dlon / 2.0);
          Real c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
          Real dist = radius * c;
          Real phi = vpol*vphi*exp(-0.5*(dist * dist)/(vrad * vrad));
          Real fcor = 2.*omega*sin(lat);

          if (vpol > 0)  // cyclone
            vel = -dist*fcor/2. + dist*fcor/2.*sqrt(1. - 4.*phi/((fcor*vrad)*(fcor*vrad)));
          else  // anticyclone
            vel = -phi*dist/(fcor*vrad*vrad);
          phydro->w(IDN,j,i) += phi;
          // Geostrophic balance still need to be implemented
        //   phydro->w(IVY,j,i) += vel*(x2 - pcoord->x2v(j))/dist;
        //   phydro->w(IVZ,j,i) += -vel*(x1 - pcoord->x1v(i))/dist;
        }
      peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord,
        is - NGHOST, ie + NGHOST, js - NGHOST, je + NGHOST, ks, ke);

      tseed = time - interval*log(uniform(gen));
  }

}

void CoriolisForcing(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &w, const AthenaArray<Real> &prim_scalar, AthenaArray<Real> const &bcc, 
  AthenaArray<Real> &u, AthenaArray<Real> &cons_scalar)
{
    for (int k = pmb->ks; k<=pmb->ke; ++k)
        for (int j = pmb->js; j <= pmb->je; ++j)
            for (int i = pmb->is; i <= pmb->ie; ++i) {
                Real R = pmb->pcoord->x1v(i);
                Real x2 = tan(pmb->pcoord->x2v(j));
                Real x3 = tan(pmb->pcoord->x3v(k));
                Real lat, lon;
                GetLatLon(&lat, &lon, pmb->pcoord, k, j, i);
                // coriolis force
                Real f = 2.*omega*sin(lat);
                u(IM1,j,i) += dt*f*w(IDN,j,i)*w(IVY,j,i);
                u(IM2,j,i) += -dt*f*w(IDN,j,i)*w(IVX,j,i);

      // topography and viscosity
    //   Real phib = Phib(pmb->pcoord,j,i);
    //   u(IM1,j,i) += dt*w(IDN,j,i)*phib*x1*pow(dist/lambda,alpha)/(dist * dist);
    //   u(IM2,j,i) += dt*w(IDN,j,i)*phib*x2*pow(dist/lambda,alpha)/(dist * dist);
    //   Real lat_now = 90. - dist/radius/M_PI*180.;
    //   Real dx = pmb->pcoord->dx1f(i);
    //   Real dy = pmb->pcoord->dx2f(j);
    //   u(IM1,j,i) += vis*dt/(dx*dy)*w(IDN,j,i)*
    //       (w(IM1,j+1,i)+w(IM1,j-1,i)+w(IM1,j,i+1)+w(IM1,j,i-1)-4*w(IM1,j,i));
    //   u(IM2,j,i) += vis*dt/(dx*dy)*w(IDN,j,i)*
    //       (w(IM2,j+1,i)+w(IM2,j-1,i)+w(IM2,j,i+1)+w(IM2,j,i-1)-4*w(IM2,j,i));
    }
}


void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // forcing parameters
  phi0 = pin->GetReal("problem", "phi0");
  omega = pin->GetReal("problem", "omega");
  interval = pin->GetReal("problem", "interval");
  vrad = pin->GetReal("problem", "vrad");
  vphi = pin->GetReal("problem", "vphi");
  vis = pin->GetOrAddReal("problem", "vis", 0.);
  polarity = pin->GetOrAddReal("problem", "polarity", 0.);
  skewness = pin->GetOrAddReal("problem", "skewness", 1.);

  // next seeding time
  tseed = time;

  // forcing function
  EnrollUserExplicitSourceFunction(CoriolisForcing);
}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
  // Vort and PV: the derivative is in an unstructured grid, implemented but correctness not verified
  SetUserOutputVariableName(4, "Vort");
  SetUserOutputVariableName(5, "PV");
}
