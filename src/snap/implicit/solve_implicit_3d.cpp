// athena
#include <athena/mesh/mesh.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/stride_iterator.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>

// exo3
#include <exo3/gnomonic_equiangle.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// snap
#include "implicit_solver.hpp"

#ifdef ENABLE_GLOG
#include <glog/logging.h>
#endif

#ifdef CUBED_SPHERE
namespace cs = CubedSphereUtility;
#endif

void ImplicitSolver::SolveImplicit3D(AthenaArray<Real> &du, AthenaArray<Real> &w,
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
        // project velocity
#ifdef CUBED_SPHERE
        auto pco = static_cast<GnomonicEquiangle*>(pmb->pcoord);
        Real cos_theta = pco->GetCosineCell(k, j);
        Real sin_theta = pco->GetSineCell(k, j);

        for (int i = 0; i < w.GetDim1(); ++i) {
          w(IVY, k, j, i) += w(IVZ, k, j, i) * cos_theta;
          w(IVZ, k, j, i) *= sin_theta;

          //cs::CovariantToContravariant(du_.at(k,j,i), cos_theta);
          //du_(IVY, k, j, i) += du_(IVZ, k, j, i) * cos_theta;
          //du_(IVZ, k, j, i) *= sin_theta;
        }
#endif  // CUBED_SPHERE
        
        // do implicit
        if (implicit_flag_ & (1 << 3)) {
          FullCorrection(du_, w, dt, k, j, is, ie);
        } else {
          PartialCorrection(du_, w, dt, k, j, is, ie);
        }

#ifdef ENABLE_GLOG
        // check for negative density and internal energy
        for (int i = is; i <= ie; i++) {
          LOG_IF(WARNING, ph->u(IEN,k,j,i) + du_(IEN,k,j,i) < 0.)
              << "rank = " << Globals::my_rank << ", (k,j,i) = "
              << "(" << k << "," << j << "," << i << ")" << std::endl
              << "(before) u[IDN] = " << ph->u(IDN,k,j,i) + du(IDN,k,j,i)
              << ", u[IVX] = " << ph->u(IVX,k,j,i) + du(IVX,k,j,i)
              << ", u[IEN] = " << ph->u(IEN,k,j,i) + du(IEN,k,j,i) << std::endl

              << "(after) u[IDN] = " << ph->u(IDN,k,j,i) + du_(IDN,k,j,i)
              << ", u[IVX] = " << ph->u(IVX,k,j,i) + du_(IVX,k,j,i)
              << ", u[IEN] = " << ph->u(IEN,k,j,i) + du_(IEN,k,j,i) << std::endl;

          LOG_IF(WARNING, ph->u(IDN,k,j,i) + du_(IDN,k,j,i) < 0.)
              << "rank = " << Globals::my_rank << ", (k,j,i) = "
              << "(" << k << "," << j << "," << i << ")"
              << "(before) u[IDN] = " << ph->u(IDN,k,j,i) + du(IDN,k,j,i)
              << ", u[IVX] = " << ph->u(IVX,k,j,i) + du(IVX,k,j,i)
              << ", u[IEN] = " << ph->u(IEN,k,j,i) + du(IEN,k,j,i) << std::endl

              << "(after) u[IDN] = " << ph->u(IDN,k,j,i) + du_(IDN,k,j,i)
              << ", u[IVX] = " << ph->u(IVX,k,j,i) + du_(IVX,k,j,i)
              << ", u[IEN] = " << ph->u(IEN,k,j,i) + du_(IEN,k,j,i) << std::endl;
        }
#endif  // ENABLE_GLOG
          
        // de-project velocity
#ifdef CUBED_SPHERE
        for (int i = 0; i < w.GetDim1(); ++i) {
          w(IVZ, k, j, i) /= sin_theta;
          w(IVY, k, j, i) -= w(IVZ, k, j, i) * cos_theta;

          //du_(IVZ, k, j, i) /= sin_theta;
          //du_(IVY, k, j, i) -= du_(IVZ, k, j) * cos_theta;
          //cs::ContravariantToCovariant(du_.at(k, j, i), cos_theta);
        }
#endif  // CUBED_SPHERE
      }

    // swap out
    du_.SwapAthenaArray(du);
  }
}
