//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details
//========================================================================================
//! \file hs94.cpp
//  \brief Problem generator for Held-Suarez-94 GCM bench mark.
//
// REFERENCE: I.M Held & M.J Suarez, "A Proposal for the Intercomparison of the
// Dynamical Cores of Atmospheric General Circulation Models"

// C++ headers
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include <athena.hpp>
#include <athena_arrays.hpp>
#include <coordinates/coordinates.hpp>
#include <cubed_sphere.hpp>
#include <eos/eos.hpp>
#include <field/field.hpp>
#include <globals.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>
#include <parameter_input.hpp>
// #include "../math/interpolation.h"
#include <utils/utils.hpp>  // replaceChar
// #include "../thermodynamics/thermodynamics.hpp"
// #include "../thermodynamics/thermodynamic_funcs.hpp"
// #include "../radiation/radiation.hpp"
// #include "../radiation/hydrogen_cia.hpp"
// #include "../radiation/freedman_mean.hpp"
// #include "../radiation/freedman_simple.hpp"
// #include "../radiation/correlatedk_absorber.hpp"
// #include "../physics/physics.hpp"

#define _sqr(x) ((x) * (x))
#define _qur(x) ((x) * (x) * (x) * (x))

using namespace std;

static Real p0, Omega, Rd, cp, sigmab, Kf, Ts, dT, dtheta, Ka, Ks, Rp, scaled_z,
    z_iso, sponge_tau, sponge_width, grav;

// \brief Held-Suarez atmosphere benchmark test. Refernce: Held & Suarez,
// (1994). Forcing parameters are given in the paper.

//! \fn void Damping(MeshBlock *pmb, Real const time, Real const dt,
//  AthenaArray<Real> const& w, AthenaArray<Real> const& bcc, AthenaArray<Real>
//  &u) \brief Pseudo radiative damping of Earth atmosphere for HS94 test.
void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &cons_scalar) {
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real omega1, omega2;
        // Load Latitude
        Real lat, lon;
        GetLatLon(&lat, &lon, pmb->pcoord, k, j, i);
        Real theta = lat;

        omega1 = cos(theta) * Omega;
        omega2 = sin(theta) * Omega;

        Real U, V;
        GetUV(&U, &V, pmb->pcoord, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);

        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * U;
        Real m3 = w(IDN, k, j, i) * V;

        u(IM1, k, j, i) -= -2. * dt * (omega1 * m3);
        Real acc2, acc3;
        Real tmp_acc2, tmp_acc3;
        tmp_acc2 = -2. * dt * (omega1 * m3);
        tmp_acc3 = -2. * dt * (omega1 * m1 - omega2 * m2);
        GetVyVz(&acc2, &acc3, pmb->pcoord, tmp_acc2, tmp_acc3, k, j, i);
        u(IM2, k, j, i) += acc2;
        u(IM3, k, j, i) += acc3;
      }
    }
  }

  Real kappa;  // Rd/Cp
  kappa = Rd / cp;
  Real iso_temp = Ts + grav * z_iso / cp;

  // Newtonian cooling and Rayleigh drag
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // Load latitude
        Real lat, lon;
        GetLatLon(&lat, &lon, pmb->pcoord, k, j, i);
        // Momentum damping coefficient, Kv
        Real scaled_z = w(IPR, k, j, i) / p0;
        Real sigma = w(IPR, k, j, i) / p0;  // pmb->phydro->pbot(k,j);
        Real sigma_p = (sigma - sigmab) / (1. - sigmab);
        Real Kv = (sigma_p < 0.0) ? 0.0 : sigma_p * Kf;

        // Temperature (energy) damping coefficient
        // Temperature difference, T - Teq
        Real Teq_p =
            Ts - dT * _sqr(sin(lat)) - dtheta * log(scaled_z) * _sqr(cos(lat));
        Teq_p *= pow(scaled_z, kappa);
        Real Teq = (Teq_p > 200.) ? Teq_p : 200.;
        Real temp =
            pmb->phydro->w(IPR, k, j, i) / pmb->phydro->w(IDN, k, j, i) / Rd;
        // Temperature damping coefficient, Kt
        sigma_p = (sigma_p < 0.0) ? 0.0 : sigma_p * _qur(cos(lat));
        Real Kt = Ka + (Ks - Ka) * sigma_p;

        // Momentum and energy damping
        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * w(IVY, k, j, i);
        Real m3 = w(IDN, k, j, i) * w(IVZ, k, j, i);

        u(IM1, k, j, i) += -dt * Kv * m1;
        u(IM2, k, j, i) += -dt * Kv * m2;
        u(IM3, k, j, i) += -dt * Kv * m3;
        u(IEN, k, j, i) +=
            -dt * (cp - Rd) * w(IDN, k, j, i) * Kt * (temp - Teq);
      }
}

Real AngularMomentum(MeshBlock *pmb, int iout) {
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
        GetLatLon(&lat, &lon, pmb->pcoord, k, j, i);
        GetUV(&U, &V, pmb->pcoord, pmb->phydro->w(IVY, k, j, i),
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

        // Originally here used cos(lat), which is x2v-pi, strange
        AMz += pmb->phydro->w(IDN, k, j, i) * vol *
               sqrt((_sqr(x1l) + _sqr(x1u)) / 2.) * cos(lat) *
               (Omega * sqrt(0.5 * (_sqr(x1l) + _sqr(x1u))) * cos(lat) + U);
      }
    }
  }

  return AMz;
}

//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real day_to_s = 8.64E4;
  // forcing parameters
  Omega = pin->GetReal("problem", "Omega");
  // thermodynamic parameters
  Real gamma = pin->GetReal("hydro", "gamma");
  grav = pin->GetReal("hydro", "grav_acc1");
  Ts = pin->GetReal("problem", "Ts");
  p0 = pin->GetReal("problem", "p0");
  dtheta = pin->GetReal("thermodynamics", "dtheta");
  dT = pin->GetReal("thermodynamics", "dT");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;
  // damping parameters
  Kf = pin->GetReal("problem", "Kf");
  Kf /= day_to_s;
  Ka = pin->GetReal("problem", "Ka");
  Ka /= day_to_s;
  Ks = pin->GetReal("problem", "Ks");
  Ks /= day_to_s;
  sigmab = pin->GetReal("problem", "sigmab");
  z_iso = pin->GetReal("problem", "z_iso");
  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);

  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, AngularMomentum, "z-angular-mom");
  return;
}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Held-Suarez problem generator
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real grav = -pin->GetReal("hydro", "grav_acc1");
  // Real grav = -phydro->psrc->GetG1();
  Real gamma = pin->GetReal("hydro", "gamma");
  p0 = pin->GetReal("problem", "p0");
  Ts = pin->GetReal("problem", "Ts");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;
  Rp = pin->GetReal("problem", "Rp");
  z_iso = pin->GetReal("problem", "z_iso");

  // Set up perturbation as the initial condition to break the symmetry of the
  // original initial condition.
  long int iseed = -1;
  Real k3 =
      10. * M_PI / (pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        if ((pcoord->x1v(i) - Rp) < z_iso) {
          Real lat, lon;
          GetLatLon(&lat, &lon, pcoord, k, j, i);
          Real x1 = pcoord->x1v(i) - Rp;
          Real temp = Ts - grav * x1 / cp;
          phydro->w(IPR, k, j, i) = p0 * pow(temp / Ts, cp / Rd);
          phydro->w(IDN, k, j, i) =
              phydro->w(IPR, k, j, i) /
              (Rd * (temp + 20. * (ran2(&iseed) - 0.5) * (1. + cos(k3 * lon))));
          phydro->w(IVX, k, j, i) = 0.;
          phydro->w(IVY, k, j, i) = 0.;
          phydro->w(IVZ, k, j, i) = 0.;
        } else {
          Real x1 = pcoord->x1v(i) - Rp;
          Real temp = Ts - grav * x1 / cp;
          Real iso_temp = Ts - grav * z_iso / cp;
          phydro->w(IPR, k, j, i) = p0 * pow(temp / Ts, cp / Rd);
          phydro->w(IDN, k, j, i) = phydro->w(IPR, k, j, i) / (Rd * iso_temp);
          phydro->w(IVX, k, j, i) = 0.;
          phydro->w(IVY, k, j, i) = 0.;
          phydro->w(IVZ, k, j, i) = 0.;
        }
      }
    }
  }

  // transfer to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  UserWorkInLoop();
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

// \brif Output distributions of temperature and potential temperature.
void MeshBlock::UserWorkInLoop() {
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real prim[NHYDRO];
        for (int n = 0; n < NHYDRO; ++n) prim[n] = phydro->w(n, j, i);
        Real temp = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) / Rd;
        user_out_var(0, k, j, i) = temp;
        user_out_var(1, k, j, i) =
            temp * pow(p0 / phydro->w(IPR, k, j, i), Rd / cp);
      }
    }
  }
}
