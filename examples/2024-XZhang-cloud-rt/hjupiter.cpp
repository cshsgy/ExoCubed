// C++ headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <impl.hpp>

// climath
#include <climath/core.h>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// astro
#include <astro/celestrial_body.hpp>

// harp
#include "harp/radiation.hpp"

Real Ps, Ts, Omega, grav, sponge_tau, bflux, gammad, Tmin;
int sponge_layer;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(7);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "lat");
  SetUserOutputVariableName(3, "lon");
  SetUserOutputVariableName(4, "vlat");
  SetUserOutputVariableName(5, "vlon");
  SetUserOutputVariableName(6, "zenith");
  SetUserOutputVariableName(7, "AMZ");
  SetUserOutputVariableName(8, "Energy");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pexo3 = pimpl->pexo3;
  auto pthermo = Thermodynamics::GetInstance();

  Real lat, lon;
  Real U, V;
  Direction ray = pimpl->prad->GetRayInput(0);
  Real zenith;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, Ps, k, j, i);

        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, phydro->w(IVY, k, j, i), phydro->w(IVZ, k, j, i),
                     k, j, i);
        user_out_var(2, k, j, i) = lat;
        user_out_var(3, k, j, i) = lon;
        user_out_var(4, k, j, i) = U;
        user_out_var(5, k, j, i) = V;

        ray = pimpl->planet->ParentZenithAngle(pmy_mesh->time, M_PI / 2. - lat,
                                               lon);
        zenith = std::acos(ray.mu) / M_PI * 180.0;
        user_out_var(6, k, j, i) = zenith;

        /* AM
        Real x1l = pcoord->x1f(i);
        Real x1u = pcoord->x1f(i + 1);
        Real xt = tan(pcoord->x2v(j));
        Real yt = tan(pcoord->x3v(k));
        Real sin_theta =
            sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

        Real x1 = tan(pcoord->x2f(j));
        Real x2 = tan(pcoord->x2f(j + 1));
        Real y = tan(pcoord->x3v(k));
        Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
        Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
        Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

        Real x = tan(pcoord->x2v(j));
        Real y1 = tan(pcoord->x3f(k));
        Real y2 = tan(pcoord->x3f(k + 1));
        delta1 = sqrt(1.0 + x * x + y1 * y1);
        delta2 = sqrt(1.0 + x * x + y2 * y2);
        Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

        Real vol = pcoord->dx1f(i) * dx2_ang * dx3_ang * sin_theta;

        user_out_var(7, k, j, i) = phydro->w(IDN, k, j, i) * vol *
               sqrt((sqr(x1l) + sqr(x1u)) / 2.) * cos(lat) *
               (Omega * sqrt(0.5 * (sqr(x1l) + sqr(x1u))) * cos(lat) + U);

        user_out_var(8, k, j, i) = phydro->w(IEN, k, j, i) * vol;*/
      }
    }
}

Real AngularMomentum(MeshBlock *pmb, int iout) {
  auto pexo3 = pmb->pimpl->pexo3;
  Real AMz = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks,
      ke = pmb->ke;

  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        Real x1l = pmb->pcoord->x1f(i);
        Real x1u = pmb->pcoord->x1f(i + 1);
        Real U, V;
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, pmb->phydro->w(IVY, k, j, i),
                     pmb->phydro->w(IVZ, k, j, i), k, j, i);

        Real xt = tan(pmb->pcoord->x2v(j));
        Real yt = tan(pmb->pcoord->x3v(k));
        Real sin_theta =
            sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

        Real x1 = tan(pmb->pcoord->x2f(j));
        Real x2 = tan(pmb->pcoord->x2f(j + 1));
        Real y = tan(pmb->pcoord->x3v(k));
        Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
        Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
        Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

        Real x = tan(pmb->pcoord->x2v(j));
        Real y1 = tan(pmb->pcoord->x3f(k));
        Real y2 = tan(pmb->pcoord->x3f(k + 1));
        delta1 = sqrt(1.0 + x * x + y1 * y1);
        delta2 = sqrt(1.0 + x * x + y2 * y2);
        Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

        Real vol = pmb->pcoord->dx1f(i) * dx2_ang * dx3_ang * sin_theta;

        AMz += pmb->phydro->w(IDN, k, j, i) * vol *
               sqrt((sqr(x1l) + sqr(x1u)) / 2.) * cos(lat) *
               (Omega * sqrt(0.5 * (sqr(x1l) + sqr(x1u))) * cos(lat) + U);
      }
    }
  }

  return AMz;
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  auto pthermo = Thermodynamics::GetInstance();
  auto prad = pmb->pimpl->prad;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);

        // coriolis force
        Real f = 2. * Omega * sin(lat);
        Real f2 = 2. * Omega * cos(lat);
        Real U, V;

        pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);

        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * U;
        Real m3 = w(IDN, k, j, i) * V;

        Real ll_acc_U = f * m3;
        Real ll_acc_V = -f * m2;
        Real acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        du(IM2, k, j, i) += dt * acc2;
        du(IM3, k, j, i) += dt * acc3;
      }

  Real Rd = pthermo->GetRd();
  Real cv = 1. / (gammad - 1.) * Rd;

  // Sponge Layer
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->ie - sponge_layer; i <= pmb->ie; ++i) {
        Real scale = cos(M_PI / 2. * (pmb->ie - i) / sponge_layer);
        du(IVX, k, j, i) -=
            w(IVX, k, j, i) * (dt / sponge_tau) * w(IDN, k, j, i) * scale;
        du(IVY, k, j, i) -=
            w(IVY, k, j, i) * (dt / sponge_tau) * w(IDN, k, j, i) * scale;
        du(IVZ, k, j, i) -=
            w(IVZ, k, j, i) * (dt / sponge_tau) * w(IDN, k, j, i) * scale;
      }
    }

  // bottom flux
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      int i = pmb->is;
      du(IEN, k, j, i) += bflux / pmb->pcoord->dx1f(i) * dt;
    }
}

Real TimeStep(MeshBlock *pmb) {
  auto prad = pmb->pimpl->prad;
  Real time = pmb->pmy_mesh->time;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      prad->CalTimeStep(pmb, k, j, pmb->is, pmb->ie);
    }

  return prad->GetTimeStep() * (1. + time / prad->GetRelaxTime());
}

//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // forcing parameters
  Omega = pin->GetReal("problem", "Omega");
  grav = -pin->GetReal("hydro", "grav_acc1");
  Ts = pin->GetReal("problem", "Ts");
  Ps = pin->GetReal("problem", "Ps");
  bflux = pin->GetOrAddReal("problem", "bflux", 0.);
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  sponge_layer = pin->GetInteger("problem", "sponge_layer");
  gammad = pin->GetReal("hydro", "gamma");
  Tmin = pin->GetReal("problem", "Tmin");

  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);
  // EnrollUserTimeStepFunction(TimeStep);

  // Z-Angular Momentum
  // AllocateUserHistoryOutput(2);
  // EnrollUserHistoryOutput(0, AngularMomentum, "z-angular-mom");
}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));

  auto pexo3 = pimpl->pexo3;
  auto pthermo = Thermodynamics::GetInstance();

  AirParcel air(AirParcel::Type::MoleFrac);
  // construct atmosphere from bottom up
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.SetZero();
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2.,
                           Thermodynamics::Method::ReversibleAdiabat, grav);

      int i = is;
      for (; i <= ie; ++i) {
        if (air.w[IDN] < Tmin) break;
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::PseudoAdiabat, grav);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::Isothermal, grav);
      }
    }
}
