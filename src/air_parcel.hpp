#ifndef SRC_AIR_PARCEL_HPP_
#define SRC_AIR_PARCEL_HPP_

// C/C++
#include <array>
#include <iostream>
#include <vector>

// athena
#include <athena/athena.hpp>  // Real

// canoe
#include <configure.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

//! \class AirParcel
//! \brief a collection of all physical data in a computational cell
class AirParcel {
 public:
  enum { Size = NHYDRO + NCLOUD + NCHEMISTRY + NTRACER + NTURBULENCE };

  enum class Type { MassFrac = 0, MassConc = 1, MoleFrac = 2, MoleConc = 3 };

  friend std::ostream &operator<<(std::ostream &os, Type const &type);
  friend std::ostream &operator<<(std::ostream &os, AirParcel const &var);

  //! reference polytropic index of dry air
  static Real const gammad_ref;

  //! adiaibatic index of dry air [1]
  //! \return $\gamma_d = \frac{c_p}{c_v}$
  Real Gammad() const;

  //! \return $\chi = \frac{R}{c_p}$
  Real Chi() const {
    auto pthermo = Thermodynamics::GetInstance();
    Real gammad = Gammad();

    Real qsig = 1., feps = 1.;
#pragma omp simd reduction(+ : qsig)
    for (int n = 1; n <= NVAPOR; ++n) {
      qsig += w[n] * (pthermo->GetCpRatioMole(n) - 1.);
    }

#pragma omp simd reduction(+ : qsig, feps)
    for (int n = 0; n < NCLOUD; ++n) {
      feps += -c[n];
      qsig += c[n] * (pthermo->GetCpRatioMole(n + 1 + NVAPOR) - 1.);
    }

    return (gammad - 1.) / gammad / qsig;
  }

  //! Scaled inverse of the mean molecular weight (with cloud)
  //! Eq.94 in Li2019
  Real RovRd() const {
    auto pthermo = Thermodynamics::GetInstance();
    Real fgas = 1., feps = 1.;

#pragma omp simd reduction(+ : feps)
    for (int n = 1; n <= NVAPOR; ++n) {
      feps += w[n] * (pthermo->GetMuRatio(n) - 1.);
    }

#pragma omp simd reduction(+ : fgas, feps)
    for (int n = 0; n < NCLOUD; ++n) {
      fgas += -c[n];
      feps += c[n] * (pthermo->GetMuRatio(n) - 1.);
    }

    return fgas / feps;
  }

  //! Mean molecular weight [kg/mol]
  Real Mu() const {
    auto pthermo = Thermodynamics::GetInstance();
    Real mud = Constants::Rgas / pthermo->GetRd();
    return mud / RovRd();
  }

  //! Specific heat capacity [J/(kg K)] at constant volume
  //! $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
  //! $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
  //! \return $c_{v,i}$
  Real CvMass(int n) const {
    auto pthermo = Thermodynamics::GetInstance();
    Real cvd = Rd_ / (Gammad() - 1.);
    return pthermo->GetCvRatioMass(n) * cvd;
  }

  //! Specific heat capacity [J/(mol K)] of the air parcel at constant volume
  //! \return $\hat{c}_v$
  Real CvMole(int n) const {
    auto pthermo = Thermodynamics::GetInstance();
    return CvMass(n) * pthermo->GetMu(n);
  }

  //! Specific heat capacity [J/(kg K)] at constant pressure
  //! $c_{p,d} = \frac{\gamma_d}{\gamma_d - 1}R_d$ \n
  //! $c_{p,i} = \frac{c_{p,i}}{c_{p,d}}\times c_{p,d}$
  //! \return $c_p$
  Real CpMass(int n) const {
    Real gammad = Gammad();
    Real cpd = Rd_ * gammad / (gammad - 1.);
    auto pthermo = Thermodynamics::GetInstance();
    return pthermo->GetCpRatioMass(n) * cpd;
  }

  //! Specific heat capacity [J/(mol K)] of the air parcel at constant pressure
  //! \return $\hat{c}_v$
  Real CpMole(int n) const {
    auto pthermo = Thermodynamics::GetInstance();
    return CpMass(n) * pthermo->GetMu(n);
  }

  //! Molar density [mol/m^3]
  Real DensityMole() const {
    Real qgas = 1.;
#pragma omp simd reduction(+ : qgas)
    for (int n = 0; n < NCLOUD; ++n) qgas += -c[n];
    return w[IPR] / (Constants::Rgas * w[IDN] * qgas);
  }

  //! Molar internal energy [J/mol]
  Real InternalEnergyMole() const {
    auto pthermo = Thermodynamics::GetInstance();

    Real cvd = Constants::Rgas / (Gammad() - 1.);
    Real fsig = 1., LE = 0.;

    for (int i = 1; i <= NVAPOR; ++i) {
      // vapor
      fsig += (pthermo->GetCvRatioMole(i) - 1.) * w[i];

      // clouds
      for (auto j : pthermo->GetCloudIndices(i)) {
        int n = j + 1 + NVAPOR;
        Real qc = c[j];

        fsig += (pthermo->GetCvRatioMole(n) - 1.) * qc;
        LE += latent_energy_mole_[n] * qc;
      }
    }

    return cvd * w[IDN] * fsig - LE;
  }

  //! \brief moist adiabatic temperature gradient
  //!
  //! $\Gamma_m = (\frac{d\ln T}{d\ln P})_m$
  //! \return $\Gamma_m$
  Real DlnTDlnP(Real latent[]) const;

  Real RelativeHumidity(int n) const {
    auto pthermo = Thermodynamics::GetInstance();
    auto rates = pthermo->TryEquilibriumTP_VaporCloud(*this, n, 0., true);
    return w[n] / (w[n] + rates[0]);
  }

  //! Thermodnamic equilibrium at current TP
  void EquilibrateTP() {
    setTotalEquivalentVapor();

    // vapor <=> cloud
    for (int i = 1; i <= NVAPOR; ++i) {
      auto rates = pthermo->TryEquilibriumTP_VaporCloud(*this, i);

      // vapor condensation rate
      w[i] += rates[0];

      // cloud concentration rates
      auto &cloud_index = pthermo->GetCloudIndices(i);
      for (int n = 1; n < rates.size(); ++n) c[cloud_index[n - 1]] += rates[n];
    }

    // vapor + vapor <=> cloud
    for (auto const &[ij, info] : cloud_reaction_map_) {
      auto rates = TryEquilibriumTP_VaporVaporCloud(*this, ij);
      auto &indx = info.first;

      // vapor condensation rate
      w[indx[0]] += rates[0];
      w[indx[1]] += rates[1];

      // cloud concentration rates
      c[indx[2]] += rates[3];
    }
  }

 protected:
  // data holder
  std::array<Real, Size> data_;

  // type
  Type mytype_;

 public:
  //! data pointers
  //! hydro data
  Real *const w;

  //! cloud data
  Real *const c;

  //! chemistry data
  Real *const q;

  //! tracer data
  Real *const x;

  //! turbulence data
  Real *const t;

  //! particle data
  Real const *d;

  // constructor
  explicit AirParcel(Type type = Type::MoleFrac)
      : mytype_(type),
        w(data_.data()),
        c(w + NHYDRO),
        q(c + NCLOUD),
        x(q + NCHEMISTRY),
        t(x + NTRACER),
        d(t + NTURBULENCE) {
    std::fill(data_.begin(), data_.end(), 0.0);
  }

  // copy constructor
  AirParcel(AirParcel const &other)
      : mytype_(other.mytype_),
        data_(other.data_),
        w(data_.data()),
        c(w + NHYDRO),
        q(c + NCLOUD),
        x(q + NCHEMISTRY),
        t(x + NTRACER),
        d(t + NTURBULENCE) {}

  // Assignment operator
  AirParcel &operator=(const AirParcel &other) {
    // Check for self-assignment
    if (this == &other) {
      return *this;
    }

    data_ = other.data_;
    mytype_ = other.mytype_;
    return *this;
  }

  void SetType(Type type) { mytype_ = type; }

  Type GetType() const { return mytype_; }

  void SetZero() { std::fill(data_.begin(), data_.end(), 0.0); }

  AirParcel &ConvertTo(AirParcel::Type type);

  AirParcel &ToMassFraction();
  AirParcel &ToMassConcentration();
  AirParcel &ToMoleFraction();
  AirParcel &ToMoleConcentration();

 protected:
  void massFractionToMassConcentration();
  void massConcentrationToMassFraction();

  void massFractionToMoleFraction();
  void moleFractionToMassFraction();

  void massConcentrationToMoleFraction();
  void moleFractionToMassConcentration();

  void moleFractionToMoleConcentration();
  void moleConcentrationToMoleFraction();

  void massConcentrationToMoleConcentration();
  void moleConcentrationToMassConcentration();

  void massFractionToMoleConcentration();
  void moleConcentrationToMassFraction();

  void setTotalEquivalentVapor() {
    auto pthermo = Thermodynamics::GetInstance();

    // vpaor <=> cloud
    for (int i = 1; i <= NVAPOR; ++i)
      for (auto &j : pthermo->GetCloudIndices(i)) {
        w[i] += c[j];
        c[j] = 0.;
      }

    // vapor + vapor <=> cloud
    for (auto const &[ij, info] : pthermo->GetCloudReactions()) {
      auto &indx = info.first;
      auto &stoi = info.second;

      w[indx[0]] += stoi[0] / stoi[2] * c[indx[2]];
      c[indx[2]] = 0.;
      w[indx[1]] += stoi[1] / stoi[2] * c[indx[2]];
      c[indx[2]] = 0.;
    }
  }
};

using AirColumn = std::vector<AirParcel>;

namespace AirParcelHelper {
AirParcel gather_from_primitive(MeshBlock const *pmb, int k, int j, int i);
AirColumn gather_from_primitive(MeshBlock const *pmb, int k, int j);

inline AirColumn gather_from_primitive(MeshBlock const *pmb, int k, int j,
                                       int il, int iu) {
  AirColumn ac(iu - il + 1);
  for (int i = il; i <= iu; ++i) {
    ac[i - il] = gather_from_primitive(pmb, k, j, i);
  }

  return ac;
}

AirParcel gather_from_conserved(MeshBlock const *pmb, int k, int j, int i);
AirColumn gather_from_conserved(MeshBlock const *pmb, int k, int j);

inline AirColumn gather_from_conserved(MeshBlock const *pmb, int k, int j,
                                       int il, int iu) {
  AirColumn ac(iu - il + 1);
  for (int i = il; i <= iu; ++i) {
    ac[i - il] = gather_from_conserved(pmb, k, j, i);
  }

  return ac;
}

void distribute_to_primitive(MeshBlock *pmb, int k, int j, int i,
                             AirParcel const &ac);
void distribute_to_primitive(MeshBlock *pmb, int k, int j, AirColumn const &ac);

inline void distribute_to_primitive(MeshBlock *pmb, int k, int j, int il,
                                    int iu, AirColumn const &ac) {
  for (int i = il; i <= iu; ++i) {
    distribute_to_primitive(pmb, k, j, i, ac[i - il]);
  }
}

void distribute_to_conserved(MeshBlock *pmb, int k, int j, int i,
                             AirParcel const &ac);
void distribute_to_conserved(MeshBlock *pmb, int k, int j, AirColumn const &ac);

inline void distribute_to_conserved(MeshBlock *pmb, int k, int j, int il,
                                    int iu, AirColumn const &ac) {
  for (int i = il; i <= iu; ++i) {
    distribute_to_conserved(pmb, k, j, i, ac[i - il]);
  }
}

}  // namespace AirParcelHelper

#endif  // SRC_AIR_PARCEL_HPP_
