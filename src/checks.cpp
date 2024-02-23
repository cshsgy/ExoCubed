// C/C++
#include <iomanip>
#include <sstream>
#include <string>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

#ifdef ENABLE_GLOG
#include <glog/logging.h>
#endif

std::string print_column_table(std::string name, AthenaArray<Real> const& var,
                               int n, int k, int j, int il, int iu,
                               int width = 1) {
  std::stringstream msg;

  msg << "rank = " << Globals::my_rank << ", " << name << " dump:" << std::endl;

  for (int i = il; i <= iu; ++i) {
    msg << "i = " << i << " |";
    msg << std::setw(8) << var(n, k, j, i) << ", ";
    msg << std::endl;
  }

  return msg.str();
}

std::string print_column_table(std::string name, AthenaArray<Real> const& var,
                               int n, int il, int iu) {
  std::stringstream msg;

  msg << "rank = " << Globals::my_rank << ", " << name << " dump:" << std::endl;

  for (int i = il; i <= iu; ++i) {
    msg << "i = " << i << " |";
    msg << std::setw(8) << var(n, i) << ", ";
    msg << std::endl;
  }

  return msg.str();
}

void check_eos_cons2prim(AthenaArray<Real> const& prim, int k, int j, int il,
                         int iu) {
#ifdef ENABLE_GLOG
  for (int i = il; i <= iu; ++i) {
    Real const& w_d = prim(IDN, k, j, i);
    Real const& w_p = prim(IPR, k, j, i);

    LOG_IF(FATAL, std::isnan(w_d) || (w_d < 0.))
        << print_column_table("prim", prim, IDN, k, j, il, iu);

    LOG_IF(FATAL, std::isnan(w_p) || (w_p < 0.))
        << print_column_table("prim", prim, IPR, k, j, il, iu);
  }
#endif  // ENABLE_GLOG
}

// dirty fix for negative pressure and density
void fix_eos_cons2prim(MeshBlock* pmb, AthenaArray<Real>& prim, int k, int j,
                       int il, int iu) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;
  Real Rd = pthermo->GetRd();
  Real grav = pmb->phydro->hsrc.GetG1();

  for (int i = il; i <= iu; ++i) {
    int ifix = il;
    for (; ifix <= iu; ++ifix) {
      if ((prim(IDN, k, j, ifix) < 0.) || std::isnan(prim(IDN, k, j, ifix)) ||
          (prim(IPR, k, j, ifix) < 0.) || std::isnan(prim(IPR, k, j, ifix))) {
        break;
      }
    }

    Real temp = pthermo->GetTemp(pmb, k, j, ifix - 1);
    for (int i = ifix; i <= iu; ++i) {
      Real z = pcoord->x1v(i) - pcoord->x1v(ifix - 1);
      prim(IDN, k, j, i) =
          prim(IDN, k, j, ifix - 1) * exp(grav * z / (Rd * temp));
      prim(IPR, k, j, i) =
          prim(IPR, k, j, ifix - 1) * exp(grav * z / (Rd * temp));
    }
  }
}

void fix_reconstruct_x2(AthenaArray<Real>& wl, AthenaArray<Real>& wr,
                        AthenaArray<Real> const& w, int k, int j, int il,
                        int iu) {
  for (int i = il; i <= iu; ++i) {
    if (wl(IDN, i) < 0.) wl(IDN, i) = w(IDN, k, j - 1, i);
    if (wr(IDN, i) < 0.) wr(IDN, i) = w(IDN, k, j, i);
    if (wl(IPR, i) < 0.) wl(IPR, i) = w(IPR, k, j - 1, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = w(IPR, k, j, i);
  }
}

void fix_reconstruct_x3(AthenaArray<Real>& wl, AthenaArray<Real>& wr,
                        AthenaArray<Real> const& w, int k, int j, int il,
                        int iu) {
  for (int i = il; i <= iu; ++i) {
    if (wl(IDN, i) < 0.) wl(IDN, i) = w(IDN, k - 1, j, i);
    if (wr(IDN, i) < 0.) wr(IDN, i) = w(IDN, k, j, i);
    if (wl(IPR, i) < 0.) wl(IPR, i) = w(IPR, k - 1, j, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = w(IPR, k, j, i);
  }
}

void check_reconstruct(AthenaArray<Real> const& wl, AthenaArray<Real> const& wr,
                       int dir, int k, int j, int il, int iu) {
#ifdef ENABLE_GLOG
  for (int i = il; i <= iu; ++i) {
    char name[80];
    snprintf(name, 80, "wl-den-%d", dir + 1);
    LOG_IF(FATAL, wl(IDN, i) < 0.) << print_column_table(name, wl, IDN, il, iu);

    snprintf(name, 80, "wr-den-%d", dir + 1);
    LOG_IF(FATAL, wr(IDN, i) < 0.) << print_column_table(name, wr, IDN, il, iu);

    snprintf(name, 80, "wl-pre-%d", dir + 1);
    LOG_IF(FATAL, wl(IPR, i) < 0.) << print_column_table(name, wl, IPR, il, iu);

    snprintf(name, 80, "wr-pre-%d", dir + 1);
    LOG_IF(FATAL, wr(IPR, i) < 0.) << print_column_table(name, wr, IPR, il, iu);
  }
#endif  // ENABLE_GLOG
}

void check_hydro_riemann_solver_flux(AthenaArray<Real> const& flux, int ivx,
                                     int k, int j, int il, int iu) {
#ifdef ENABLE_GLOG
  for (int i = il; i <= iu; ++i) {
    for (int n = 0; n < NHYDRO; ++n) {
      char name[80];
      snprintf(name, 80, "flux%d-%d", ivx, n);
      LOG_IF(FATAL, std::isnan(flux(n, k, j, i)))
          << print_column_table(name, flux, n, k, j, il, iu);
    }
  }
#endif  // ENABLE_GLOG
}

void check_decomposition(AthenaArray<Real> const& wl,
                         AthenaArray<Real> const& wr, int k, int j, int il,
                         int iu) {
#ifdef ENABLE_GLOG
  for (int i = il; i <= iu; ++i) {
    LOG_IF(FATAL, wl(IDN, k, j, i) < 0.)
        << print_column_table("wl-den", wl, IDN, k, j, il, iu);

    LOG_IF(FATAL, wr(IDN, k, j, i) < 0.)
        << print_column_table("wr-den", wr, IDN, k, j, il, iu);

    LOG_IF(ERROR, wl(IPR, k, j, i) < 0.)
        << print_column_table("wl-pre", wl, IPR, k, j, il, iu);

    LOG_IF(ERROR, wr(IPR, k, j, i) < 0.)
        << print_column_table("wl-pre", wr, IPR, k, j, il, iu);
  }
#endif  // ENABLE_GLOG
}

void check_implicit_cons(AthenaArray<Real> const& cons, int il, int iu, int jl,
                         int ju, int kl, int ku) {
#ifdef ENABLE_GLOG
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        Real const& u_d = cons(IDN, k, j, i);
        Real const& u_e = cons(IEN, k, j, i);

        LOG_IF(FATAL, std::isnan(u_d) || (u_d < 0.))
            << print_column_table("cons", cons, IDN, k, j, il, iu);

        LOG_IF(FATAL, std::isnan(u_e) || (u_e < 0.))
            << print_column_table("cons", cons, IEN, k, j, il, iu);
      }
#endif  // ENABLE_GLOG
}

void fix_implicit_cons(MeshBlock* pmb, AthenaArray<Real>& cons, int il, int iu,
                       int jl, int ju, int kl, int ku) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;
  Real Rd = pthermo->GetRd();
  Real grav = pmb->phydro->hsrc.GetG1();

  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      int ifix = il;

      for (; ifix <= iu; ++ifix) {
        if ((cons(IDN, k, j, ifix) < 0.) || std::isnan(cons(IDN, k, j, ifix)) ||
            (cons(IEN, k, j, ifix) < 0.) || std::isnan(cons(IEN, k, j, ifix))) {
          break;
        }
      }

      AirParcel air0 =
          AirParcelHelper::gather_from_conserved(pmb, k, j, ifix - 1);
      air0.ToMoleFraction();
      Real temp = air0.w[IDN];

      for (int i = ifix; i <= iu; ++i) {
        AirParcel air = AirParcelHelper::gather_from_conserved(pmb, k, j, i);
        air.ToMassFraction();

        Real z = pcoord->x1v(i) - pcoord->x1v(ifix - 1);
        air.w[IDN] = cons(IDN, k, j, ifix - 1) * exp(grav * z / (Rd * temp));
        air.w[IPR] = cons(IPR, k, j, ifix - 1) * exp(grav * z / (Rd * temp));

        AirParcelHelper::distribute_to_conserved(pmb, k, j, i, air);
      }
    }
}
