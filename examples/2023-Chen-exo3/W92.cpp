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

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Application::Logger app("main");
  app->Log("ProblemGenerator: Cubed Sphere");

  auto pexo3 = pimpl->pexo3;

  // Parameters
  Real h0 = 8000.;
  Real g = 9.80616;
  Real omega = 7.848E-6;
  Real a = 6.37122E6;
  Real K = 7.848E-6;
  Real R = 4.0;
  Real om_earth = 7.292E-5;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x = tan(pcoord->x2v(j));
        Real y = tan(pcoord->x3v(k));
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);

        // h terms, formatted to avoid singularity at poles
        Real A = 0.5 * omega * (2 * om_earth + omega) * cos(lat) * cos(lat) +
                 0.25 * K * K *
                     ((R + 1) * pow(cos(lat), 2 * R + 2) +
                      (2 * R * R - R - 2) * pow(cos(lat), 2 * R) -
                      2 * R * R * pow(cos(lat), 2 * R - 2));
        Real B -
            2 * (om_earth + omega) * K / (R + 1) / (R + 2) * pow(cos(lat), R) *
                ((R * R + 2 * R + 2) - (R + 1) * (R + 1) * cos(lat) * cos(lat));
        Real C = 0.25 * K * K * pow(cos(lat), 2 * R) *
                 ((R + 1) * cos(lat) * cos(lat) - (R + 2));

        phydro->w(IDN, k, j, i) = 500.0;
        Real Vy, Vz;
        pexo3->GetVyVz(&Vy, &Vz, U, V, k, j, i);
        phydro->w(IVY, k, j, i) = Vy;
        phydro->w(IVZ, k, j, i) = Vz;
      }
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
