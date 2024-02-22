//thermodynamics/saturation_adjustment.cpp  as reference

// variable names with a _ following remain to be obtained from other places

void convective_adjustment(std::vector<AirParcel>& air, Coordinates *pcoord, Real grav,
													 int k, int j, int il, int iu) {
  auto pthermo = Thermodynamics::GetInstance();
  Real Rd = pthermo->GetRd();
	Real gammad = pthermo->GetGammad();
	Real cp = 
  
  Real total_energy = 0.;
  Real total_mass = 0.;
  
  AirParcel* parcel;

  // sum the energy and mass of all air parcels that to be adjusted (conservation of energy and mass)
  for (int i = il; i <= iu; ++i) {
    parcel = &air[i];
    //auto&& air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);
    parcel.ToMassConcentration();    // if this is necessary?  parcel->ToMassConcentration() ?
    
		Real volume = pcoord->GetCellVolume(k, j, i);
    Real temp = parcel->w[IPR] / Rd / parcel->w[IDN]; // IDN is density, IPR is pressure
		                                                  // primitive
    total_mass += parcel->w[IDN] * volume;
    total_energy += parcel->w[IDN] * volume * (cv_ * temp + grav * pcoord->x1v(i));

    if (i == il) {
      Real pres_l = parcel->w[IPR];
      Real rho_l = parcel->w[IDN];
    }
    if (i == iu) {
      Real pres_u = parcel->w[IPR];
      Real rho_u = parcel->w[IDN];
    }
  }
  
  Real guess_pres_0 = 0.5 * (pres_l + pres_u); 
  Real guess_rho_0 = 0.5 * (rho_l + rho_u);
  
  bool energy_conservd = false;
  bool mass_conversed = false;
  
  while (energy_conservd == false || mass_conversed == false) {
		Real new_total_energy = 0.;
		Real new_total_mass = 0.;
		for (int i = il; i <= iu; ++i) {
			parcel = &air[i];
			parcel.w[IDN] = guess_temp_0 - grav / cp_ * (pcoord->x1v(i) - pcoord->x1v(il));
			rho = guess_rho_0 * pow(parcel.w[IDN] / guess_temp_0, cp_ / Rd -1); 
			parcel.w[IPR] = rho * Rd * parcel.w[IDN];
			
			new_total_mass += rho * volume_; 
			new_total_energy += rho * volume_ * (cv_ * parcel.w[IDN] + grav * pcoord->x1v(i));
		}

		if (new_total_energy - total_energy > 5.) {
			guess_temp_0 -= 0.5;  // need control guess_temp_0 positive
														// how much delta P -> how much delta E
		} else if (new_total_energy - total_energy < -5.) {
			guess_temp_0 += 0.5;
		} else {
			energy_conservd = true;
		}

		if (new_total_mass - total_mass > 5.) {
			guess_rho_0 -= 0.5;   // need control guess_rho_0 positive
														// how much delta rho -> how much delta mass
		} else if (new_total_mass - total_mass < -5.) {
			guess_rho_0 += 0.5;
		} else {
			energy_conservd = true;
		}
  }

#pragma omp simd reduction(+ : IE, rho)
  for (int n = 0; n <= NVAPOR; ++n) {
    IE += air.w[n] * GetCpMassRef(n) * temp;
    rho += air.w[n];
  }

#pragma omp simd reduction(+ : IE, LE, rho)
  for (int n = 0; n < NCLOUD; ++n) {
    IE += air.c[n] * GetCpMassRef(1 + NVAPOR + n) * temp;
    LE += -latent_energy_mass_[1 + NVAPOR + n] * air.c[n];
    rho += air.c[n];
  }

  return (IE + LE) / rho + gz;
}
