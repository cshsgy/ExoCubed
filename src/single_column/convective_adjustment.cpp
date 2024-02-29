// C/C++
#include <array>
#include <cmath>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/stride_iterator.hpp>

// canoe
#include <air_parcel.hpp>

// climath
#include <climath/broyden_root.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// scm
#include "single_column.hpp"

struct TPBottomSolver {
  MeshBlock *pmb;
  SingleColumn *pscm;
  int k, j, il, iu;
  Real mass0, enthalpy0, vol0;
};

//! wrapper function passing to broyden_root
//! over SingleColumn::findTPBottom
void find_tp_bottom(int n, double *x, double *f, void *arg) {
  TPBottomSolver *psolver = static_cast<TPBottomSolver *>(arg);

  auto result = psolver->pscm->findTPBottom(
      x[0], x[1], psolver->mass0, psolver->enthalpy0, psolver->pmb, psolver->k,
      psolver->j, psolver->il, psolver->iu);

  std::cout << "error0 = " << result[0] << std::endl;
  std::cout << "error1 = " << result[1] << std::endl;

  f[0] = result[0];
  f[1] = result[1];
}

std::array<Real, 2> SingleColumn::findTPBottom(Real tempf, Real presf,
                                               Real mass0, Real enthalpy0,
                                               MeshBlock *pmb, int k, int j,
                                               int il, int iu) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;
  auto phydro = pmb->phydro;
  Real grav = -phydro->hsrc.GetG1();
  Real Tmin = GetPar<Real>("Tmin");

  // adjust pressure and density, and examine mass and energy conservation
  Real mass1 = 0., enthalpy1 = 0., vol = 0.;
  AirParcel air = AirParcelHelper::gather_from_conserved(pmb, k, j, il);

  // half a grid to cell center
  air.ToMoleFraction();
  air.w[IDN] = tempf;
  air.w[IPR] = presf;

  std::cout << "tempf = " << tempf << std::endl;
  std::cout << "presf = " << presf << std::endl;

  pthermo->Extrapolate(&air, pcoord->dx1f(il) / 2.,
                       Thermodynamics::Method::ReversibleAdiabat, grav);

  air.ToMassFraction();
  mass1 += air.w[IDN] * vol_(il);
  enthalpy1 += (air.w[IEN] + air.w[IDN] * grav * pcoord->x1v(il)) * vol_(il);
  vol += vol_(il);

  int i = il;
  for (; i <= iu; ++i) {
    if (air.w[IDN] < Tmin) break;
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    air.ToMassFraction();
    mass1 += air.w[IDN] * vol_(il);
    enthalpy1 += (air.w[IEN] + air.w[IDN] * grav * pcoord->x1v(il)) * vol_(il);
    vol += vol_(il);
  }

  // Replace adiabatic atmosphere with isothermal atmosphere
  for (; i <= iu; ++i) {
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::Isothermal, grav);
    air.ToMassFraction();
    mass1 += air.w[IDN] * vol_(il);
    enthalpy1 += (air.w[IEN] + air.w[IDN] * grav * pcoord->x1v(il)) * vol_(il);
    vol += vol_(il);
  }

  std::cout << std::endl << "mass1 = " << mass1 << std::endl;
  std::cout << "enthalpy1 = " << enthalpy1 << std::endl;

  // return std::array<Real, 2>({(mass1 - mass0) / mass0,
  //     (enthalpy1 - enthalpy0) / enthalpy0});
  return std::array<Real, 2>(
      {(mass1 - mass0) / vol, (enthalpy1 - enthalpy0) / vol});
}

// type of the function
std::array<int, 2> SingleColumn::findUnstableRange(MeshBlock *pmb, int k,
                                                   int j) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;
  auto phydro = pmb->phydro;

  Real den_tol = GetPar<Real>("den_tol");
  Real grav = -pmb->phydro->hsrc.GetG1();

  // determine il where convective adjustment starts
  int il = pmb->is;
  for (; il <= pmb->ie; ++il) {
    AirParcel air = AirParcelHelper::gather_from_conserved(pmb, k, j, il);
    pthermo->Extrapolate(&air, pcoord->dx1f(il),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    air.ToMassFraction();
    Real density = phydro->u(k, j, il + 1);
    if (air.w[IDN] - density < -den_tol) break;
  }

  // determine iu where convective adjustment stops
  int iu = il + 1;
  for (; iu < pmb->ie; ++iu) {
    AirParcel air = AirParcelHelper::gather_from_conserved(pmb, k, j, iu);
    pthermo->Extrapolate(&air, pcoord->dx1f(iu),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    air.ToMassFraction();
    Real density = phydro->u(k, j, iu + 1);
    if (air.w[IDN] - density > den_tol) break;
  }

  return std::array<int, 2>({il, iu});
}

void SingleColumn::ConvectiveAdjustment(MeshBlock *pmb, int k, int j) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;
  auto phydro = pmb->phydro;

  Real tp_tol = GetPar<Real>("rel_tol");
  int max_iter_broyden = GetPar<int>("max_iter_broyden");
  int max_iter = GetPar<int>("max_iter");

  Real Tmin = GetPar<Real>("Tmin");
  Real grav = -pmb->phydro->hsrc.GetG1();

  TPBottomSolver solver;
  solver.pmb = pmb;
  solver.pscm = this;
  solver.k = k;
  solver.j = j;

  auto range = findUnstableRange(pmb, k, j);
  AthenaArray<Real> &u = pmb->phydro->u;

  std::cout << "Find unstable range" << std::endl;
  std::cout << range[0] << " " << range[1] << std::endl;

  int iter = 0;
  while ((range[0] != -1) && (iter++ < max_iter)) {
    solver.il = range[0];
    solver.iu = range[1];

    pcoord->CellVolume(k, j, solver.il, solver.iu, vol_);

    // sum the energy and mass of all air parcels that to be adjusted
    // (conservation of energy and mass)
    solver.mass0 = 0.;
    solver.enthalpy0 = 0.;
    for (int i = solver.il; i <= solver.iu; ++i) {
      Real density = phydro->u(k, j, i);
      solver.mass0 += density * vol_(i);
      solver.enthalpy0 +=
          (phydro->u(IEN, k, j, i) + density * grav * pcoord->x1v(i)) * vol_(i);
    }

    std::cout << "mass0 = " << solver.mass0 << std::endl;
    std::cout << "enthalpy0 = " << solver.enthalpy0 << std::endl;

    // take a guess that the temperature at the bottom of the
    // convective adjustment column is the same as the cell averaged
    // temperature (slightly colder than an adiabatic profile)
    AirParcel air =
        AirParcelHelper::gather_from_conserved(pmb, k, j, solver.il);
    air.ToMoleFraction();
    Real tempf = air.w[IDN];

    // take a guess that the pressure at the bottom of the
    // convective adjustment column is the integrate mass
    Real presf = air.w[IPR] * 1.2;

    std::cout << "guess tempf = " << tempf << std::endl;
    std::cout << "guess presf = " << presf << std::endl;

    Real tp_bot[2] = {tempf, presf};
    int status = broyden_root(2, tp_bot, find_tp_bottom, tp_tol,
                              max_iter_broyden, &solver);
    if (status != 0) {
      throw std::runtime_error(
          "ConvectiveAdjustment: "
          "Broyden's method failed to converge");
    }

    std::cout << "final tempf = " << tp_bot[0] << std::endl;
    std::cout << "final presf = " << tp_bot[1] << std::endl;

    air.w[IDN] = tp_bot[0];
    air.w[IPR] = tp_bot[1];

    // half a grid to cell center
    pthermo->Extrapolate(&air, pcoord->dx1f(solver.il) / 2.,
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    AirParcelHelper::distribute_to_conserved(pmb, k, j, solver.il, air);

    int i = solver.il;
    for (; i <= solver.iu; ++i) {
      if (air.w[IDN] < Tmin) break;
      pthermo->Extrapolate(&air, pcoord->dx1f(i),
                           Thermodynamics::Method::ReversibleAdiabat, grav);
      AirParcelHelper::distribute_to_conserved(pmb, k, j, i, air);
    }

    // Replace adiabatic atmosphere with isothermal atmosphere
    for (; i <= solver.iu; ++i) {
      pthermo->Extrapolate(&air, pcoord->dx1f(i),
                           Thermodynamics::Method::Isothermal, grav);
      AirParcelHelper::distribute_to_conserved(pmb, k, j, i, air);
    }

    range = findUnstableRange(pmb, k, j);
  }
}
