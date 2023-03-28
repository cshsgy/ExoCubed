//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//! \brief Calculate hydro/MHD fluxes

// C headers

// C++ headers
#include <algorithm>   // min,max
#include <sstream>     // string stream

// Athena++ headers
#include <athena.hpp>
#include <athena_arrays.hpp>
#include <coordinates/coordinates.hpp>
#include <eos/eos.hpp>   // reapply floors to face-centered reconstructed states
#include <field/field.hpp>
#include <field/field_diffusion/field_diffusion.hpp>
#include <gravity/gravity.hpp>
#include <reconstruct/reconstruction.hpp>
#include <scalars/scalars.hpp>
#include <hydro/hydro.hpp>
#include <hydro/hydro_diffusion/hydro_diffusion.hpp>
#include <cubed_sphere.hpp>

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//! \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                            AthenaArray<Real> &bcc, const int order) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  // b and bcc only used in Magnetic fields for now, should be able to do it quick
  // b,bcc are passed as fn parameters because clients may want to pass different bcc1,
  // b1, b2, etc., but the remaining members of the Field class are accessed directly via
  // pointers because they are unique. NOTE: b, bcc are nullptrs if no MHD.

  AthenaArray<Real> &flux_fc = scr1_nkji_;
  AthenaArray<Real> &laplacian_all_fc = scr2_nkji_;

  //--------------------------------------------------------------------------------------
  // i-direction

  AthenaArray<Real> &x1flux = flux[X1DIR];
  // set the loop limits

#ifndef HYDROSTATIC // sw not do X1

  jl = js, ju = je, kl = ks, ku = ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      }
      pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#ifdef CUBED_SPHERE // Rieman solver run later
      SaveLR3DValues(wl_, wr_, X1DIR, k, j, is, ie+1); // is to ie+1 is what the RiemannSolver below uses...
#else
      RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, x1flux, dxw_);

#endif
    }
  }
#endif  // end HYDROSTATIC

  //--------------------------------------------------------------------------------------
  // j-direction
  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux = flux[X2DIR];
    // set the loop limits
    il = is-1, iu = ie+1, kl = ks, ku = ke;

    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      }
      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, w, bcc, wlb_, wr_);
        }
        pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
#ifdef CUBED_SPHERE // Rieman solver run later
        SaveLR3DValues(wl_, wr_, X2DIR, k, j, il, iu); // il to iu is what the RiemannSolver below uses...
        wl_.SwapAthenaArray(wlb_);
#else
        RiemannSolver(k, j, il, iu, IVY, wl_, wr_, x2flux, dxw_);
        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);

#endif
      }
    }
  }

  //--------------------------------------------------------------------------------------
  // k-direction
  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux = flux[X3DIR];
    // set the loop limits
    il = is, iu = ie, jl = js, ju = je;

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      }
      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, w, bcc, wlb_, wr_);
        }
        pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#ifdef CUBED_SPHERE // Rieman solver run later
        SaveLR3DValues(wl_, wr_, X3DIR, k, j, il, iu); // il to iu is what the RiemannSolver below uses...
        wl_.SwapAthenaArray(wlb_);
#else

        RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, x3flux, dxw_);
        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);

#endif
      }
    }
  }

// Cubed Sphere: recover the stored values, run riemann solvers
#ifdef CUBED_SPHERE
  // Temporarily comment out the following 2 lines to avoid the error
  SynchronizeFluxesSend();
  SynchronizeFluxesRecv();
  //--------------------------------------------------------------------------------------
  // i-direction
#ifndef HYDROSTATIC
  jl = js, ju = je, kl = ks, ku = ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      LoadLR3DValues(wl_, wr_, X1DIR, k, j, is, ie+1);
      pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
      RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, x1flux, dxw_);
    }
  }
#endif // HYDROSTATIC

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux = flux[X2DIR];
    // set the loop limits
    il = is-1, iu = ie+1, kl = ks, ku = ke;

    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
        LoadLR3DValues(wl_, wr_, X2DIR, k, j, il, iu); // il to iu is what the RiemannSolver below uses...
        RiemannSolver(k, j, il, iu, IVY, wl_, wr_, x2flux, dxw_);
      }
    }
  }


  //--------------------------------------------------------------------------------------
  // k-direction

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux = flux[X3DIR];
    // set the loop limits
    il = is, iu = ie, jl = js, ju = je;

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      for (int k=ks; k<=ke+1; ++k) {
        pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
        LoadLR3DValues(wl_, wr_, X3DIR, k, j, il, iu);

        RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, x3flux, dxw_);
      }
    }
  }
#endif

  if (!STS_ENABLED)
    AddDiffusionFluxes();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes_STS
//! \brief Calculate Hydrodynamic Diffusion Fluxes for STS

void Hydro::CalculateFluxes_STS() {
  AddDiffusionFluxes();
}

void Hydro::AddDiffusionFluxes() {
  Field *pf = pmy_block->pfield;
  // add diffusion fluxes
  if (hdif.hydro_diffusion_defined) {
    if (hdif.nu_iso > 0.0 || hdif.nu_aniso > 0.0)
      hdif.AddDiffusionFlux(hdif.visflx,flux);
    if (NON_BAROTROPIC_EOS) {
      if (hdif.kappa_iso > 0.0 || hdif.kappa_aniso > 0.0)
        hdif.AddDiffusionEnergyFlux(hdif.cndflx,flux);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED && NON_BAROTROPIC_EOS) {
    if (pf->fdif.field_diffusion_defined)
      pf->fdif.AddPoyntingFlux(pf->fdif.pflux);
  }
  return;
}
