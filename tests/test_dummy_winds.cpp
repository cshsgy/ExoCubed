#include <parameter_input.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>
#include <configure.hpp>
#include <cubed_sphere.hpp>
#include "../src/eos/eos.hpp"
#include "../src/field/field.hpp"

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Input field is a uniform 10 m/s zonal wind
  Real R = 6371000.0;
  Real V = 0.0/R;
  Real U = 0.0/R;

  std::cout << "======Start of Problem generator======" << loc.lx2 << "||" << loc.lx3 <<std::endl;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real rad = sqrt(pcoord->x3v(k) * pcoord->x3v(k) + pcoord->x2v(j) * pcoord->x2v(j));
        if ((rad < 0.4) && (i == is)) // Circular dam breaking
          phydro->w(IDN,k,j,i) = 1.2; // / (R*R);
        else
          phydro->w(IDN,k,j,i) = 1.0; // / (R*R);
        // (*Following used for comm test*) (k-ks)*(je-js)*(ie-is) + (j-js)*(ie-is) + (i-is) + 0.1*(loc.lx2+1) + 0.01*(loc.lx3+1);
        //phydro->w(IPR,k,j,i) = 1.0;
        Real Vy, Vz;
        GetVyVz(&Vy, &Vz, pcoord, U, V, k, j, i);
        phydro->w(IVY,k,j,i) = Vy;
        phydro->w(IVZ,k,j,i) = Vz;
      }
  std::cout << "Done with problem generator..." << std::endl;
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

// void MeshBlock::UserWorkInLoop()
// {
//   // Print density if not equal to 1.0
//   int flag = 0;
//   for (int k = ks; k <= ke; ++k)
//     for (int j = js; j <= je; ++j)
//       for (int i = is; i <= ie; ++i) {
//         if (phydro->w(IDN,k,j,i) != 1.0) {
//           std::cout << "Density not equal to 1.0 at (" << i << "," << j << "," << k << "): " << phydro->w(IDN,k,j,i) << std::endl;
//           // print the velocities at the location
//           std::cout << "Velocities: " << phydro->w(IVX,k,j,i) << "," << phydro->w(IVY,k,j,i) << "," << phydro->w(IVZ,k,j,i) << std::endl;
//           std::cout << "Location: " << loc.lx2 << "," << loc.lx3 << std::endl;
//           std::cout << "Time: " << pmy_mesh->time << std::endl;
//           // Terminate the program after loops done
//           flag = 1;
//         }
//       }
//   // if (flag == 1){
//   //   std::cout << "Terminating the program..." << std::endl;
//   //   exit(EXIT_FAILURE);
//   // }
// }

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
        Real f = 0.0;
        Real U, V;
        GetUV(&U, &V, pcoord, Vy, Vz, k, j, i);
        user_out_var(2,k,j,i) = U*6371000.0;
        user_out_var(3,k,j,i) = V*6371000.0;
      }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(4);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
}
