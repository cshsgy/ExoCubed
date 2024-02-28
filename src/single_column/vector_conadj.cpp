#include <cmath>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>

// canoe
#include <air_parcel.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

std::array<Real, 2> SingleColumn::findTPBottom(AthenaArray<Real> const& u,
                                               int k, int j, int il, int iu) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmy_block_->pcoord;
  Real grav = -phydro->hsrc.GetG1();

  pcoord->CellVolumn(k, j, il, iu, vol_);

  Real mass0 = 0., enthalpy0 = 0.;
  Real presf = sqrt(pthermo->GetPres(u.at(k, j, iu)) *
                    pthermo->GetPres(u.at(k, j, iu + 1)));
  Real tempf = pthermo->GetTemp(u.at(k, j, i));

  // sum the energy and mass of all air parcels that to be adjusted
  // (conservation of energy and mass)
  for (int i = il; i <= iu; ++i) {
    mass0 += u(IDN, k, j, i) * vol_(i);
    enthalpy0 +=
        (u(IEN, k, j, i) + u(IDN, k, j, i) * grav * pcoord->x1v(i)) * vol_(i);
    presf += u(IDN, k, j, i) * pcoord->dx1f(i) * grav;
  }

  // adjust pressure and density, and examine mass and energy conservation
  Real mass1 = 0., enthalpy = 0.;
  AirParcel air(u.at(k, j, il), AirParcel::Type::MoleConc);

  // half a grid to cell center
  air.ToMoleFrac();
  air.w[IDN] = tempf;
  air.w[IPR] = presf;

  pthermo->Extrapolate(&air, pcoord->dx1f(il) / 2.,
                       Thermodynamics::Method::ReversibleAdiabat, grav);

  air.ToMassConcentration();
  mass1 += air.w[IDN] * vol_(il);
  enthalpy1 += (air.w[IEN] + air.[IDN] * grav * pcoord -> x1v(il)) * vol_(il);

  int i = is;
  for (; i <= iu; ++i) {
    if (air.w[IDN] < Tmin) break;
    air.ToMoleFraction();
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::ReversibleAdiabat, grav);

    air.ToMassConcentration();
    mass1 += air.w[IDN] * vol_(il);
    enthalpy1 += (air.w[IEN] + air.[IDN] * grav * pcoord -> x1v(il)) * vol_(il);
  }

  // Replace adiabatic atmosphere with isothermal atmosphere
  for (; i <= ie; ++i) {
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::Isothermal, grav);

    mass1 += air.w[IDN] * vol_(il);
    enthalpy1 += (air.w[IEN] + air.[IDN] * grav * pcoord -> x1v(il)) * vol_(il);
  }

  return std::array<Real, 2>({mass1 - mass0, enthalpy1 - enthalpy0});
}

// type of the function
std::array<int, 2> SingleColumn::findUnstableRange(AthenaArray<Real> const& u,
                                                   int k, int j) {
  Real den_tol = GetPar<Real>("den_tol");

  // determine from il to iu where convection adjustment is performed
  int il = pmy_block_->is;
  for (; il <= pmy_block_->ie; ++il) {
    AirParcel air(u.at(k, j, il), AirParcel::Type::MassConc);
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    air.ToMassConcentration()
        : if (air.w[IDN] - u(k, j, il + 1) < -den_tol) break;
  }

  // determine iu where theta starts to increase with height
  int iu = il + 1;

  for (; iu < pmy_block_->ie; ++iu) {
    AirParcel air(u.at(k, j, il), AirParcel::Type::MassConc);
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    air.ToMassConcentration()
        : if (air.w[IDN] - u(k, j, il + 1) > den_tol) break;
  }

  return std::array<int, 2>({il, iu});
}

void SingleColumn::ConvectiveAdjustment(AthenaArray<Real>& u, int k, int j) {
  auto range = findUnstableRange(u, k, j);

  while (range.first != -1) {
    auto tp_bot = findTPBottom(u, k, j, range.first, range.second);

    Real temp = tp_bot.first;
    Real pres = tp_bot.second;

    range = findUnstableRange(u, k, j)
  }
}
