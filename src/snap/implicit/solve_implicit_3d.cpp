// athena
#include <athena/mesh/mesh.hpp>
#include <athena/hydro/hydro.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>

// snap
#include "implicit_solver.hpp"

#ifdef ENABLE_GLOG
#include <glog/logging.h>
#endif

void ImplicitSolver::SolveImplicit3D(AthenaArray<Real> &du, AthenaArray<Real> const &w,
      Real dt)
{
  auto pmb = pmy_block_;
  auto ph = pmb->phydro;
  int is, ie, js, je, ks, ke;

  // X3DIR
  if ((implicit_flag_ & (1 << 2)) && (pmb->ncells3 > 1)) {
    SetDirection(X3DIR);
    FindNeighbors();

    ks = pmb->js, js = pmb->is, is = pmb->ks;
    ke = pmb->je, je = pmb->ie, ie = pmb->ke;

    // shuffle dimension
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) du_(n, k, j, i) = du(n, i, k, j);

    // do implicit 
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        if (implicit_flag_ & (1 << 3))
          FullCorrection(du_, w, dt, k, j, is, ie);
        else
          PartialCorrection(du_, w, dt, k, j, is, ie);
      }

    // shuffle back
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN, i, k, j) = du_(IDN, k, j, i);
          du(IVZ, i, k, j) = du_(IVZ, k, j, i);
          du(IEN, i, k, j) = du_(IEN, k, j, i);
        }

  }

  // X2DIR
  if ((implicit_flag_ & (1 << 1)) && (pmb->ncells2 > 1)) {
    SetDirection(X2DIR);
    FindNeighbors();

    ks = pmb->is, js = pmb->ks, is = pmb->js;
    ke = pmb->ie, je = pmb->ke, ie = pmb->je;

    // shuffle dimension
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) du_(n, k, j, i) = du(n, j, i, k);

    // do implicit
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        if (implicit_flag_ & (1 << 3))
          FullCorrection(du_, w, dt, k, j, is, ie);
        else
          PartialCorrection(du_, w, dt, k, j, is, ie);
      }

    // shuffle back
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN, j, i, k) = du_(IDN, k, j, i);
          du(IVY, j, i, k) = du_(IVY, k, j, i);
          du(IEN, j, i, k) = du_(IEN, k, j, i);
        }
  }

  // X1DIR
  if (implicit_flag_ & 1) {
    SetDirection(X1DIR);
    FindNeighbors();

    ks = pmb->ks, js = pmb->js, is = pmb->is;
    ke = pmb->ke, je = pmb->je, ie = pmb->ie;

    // swap in
    du_.SwapAthenaArray(du);

    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        // save a copy of du for redo
        for (int n = 0; n < NHYDRO; ++n)
          for (int i = 0; i < du_.GetDim1(); ++i)
            du(n, k, j, i) = du_(n, k, j, i);

        bool redo;
        int factor = 1, iter = 0, max_iter = 5;
        do {
          redo = false;
          iter++;
          // reset du
          if (factor > 1) {
            for (int n = 0; n < NHYDRO; ++n)
              for (int i = 0; i < du_.GetDim1(); ++i)
                du_(n, k, j, i) = du(n, k, j, i);
          }

          // do implicit
          if (implicit_flag_ & (1 << 3))
            FullCorrection(du_, w, dt * factor, k, j, is, ie);
          else
            PartialCorrection(du_, w, dt * factor, k, j, is, ie);

          // check for negative density and internal energy
          for (int i = is; i <= ie; i++) {
#ifdef ENABLE_GLOG
            LOG_IF(ERROR, ph->u(IEN,k,j,i) + du_(IEN,k,j,i) < 0.)
                << "rank = " << Globals::my_rank << ", (k,j,i) = "
                << "(" << k << "," << j << "," << i << ")"
                << ", u[IEN] = " << ph->u(IEN,k,j,i) + du(IEN,k,j,i)
                << ", u[IVX] = " << ph->u(IVX,k,j,i) << std::endl;

            LOG_IF(ERROR, ph->u(IDN,k,j,i) + du_(IDN,k,j,i) < 0.)
                << "rank = " << Globals::my_rank << ", (k,j,i) = "
                << "(" << k << "," << j << "," << i << ")"
                << ", u[IDN] = " << ph->u(IDN,k,j,i) + du_(IDN,k,j,i)
                << ", u[IVX] = " << ph->u(IVX,k,j,i) << std::endl;
#endif
            if ((ph->u(IEN,k,j,i) + du_(IEN,k,j,i) < 0.) ||
                (ph->u(IDN,k,j,i) + du_(IDN,k,j,i) < 0.)) {
              factor *= 2;
              redo = true;
              break;
            }
          }
          
          if (iter > max_iter) {
            throw RuntimeError("solver_implicit_3d", "negative density or internal energy");
          }
        } while (redo);
      }

    // swap out
    du_.SwapAthenaArray(du);
  }
}
