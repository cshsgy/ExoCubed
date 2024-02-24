#include "air_parcel.hpp"
// snap
#include "snap/thermodynamics/thermodynamics.hpp"

//thermodynamics/saturation_adjustment.cpp  as reference

// variable names with a _ following remain to be obtained from other places

auto pthermo = Thermodynamics::GetInstance();
Real Rd = pthermo->GetRd();
Real gammad = pthermo->GetGammad();
Real cp = Rd * gammad / (gammad-1.);

// type of the function
void search_convective_adjustment(std::vector<AirParcel>& air_column, Coordinates *pcoord, Real grav,
													 int k, int j) {
	
	size_t length = air_column.size();

	// initialize il and iu, the lower and upper level between which the airparcels need
	// convective adjustment 
	int il = -1;
	int iu = ie;   //REVISE ie AS LENGTH
	
	// determine il where theta starts to decrease with height	
	for (int i = is; i <= ie-1; ++i) {   // where to get is and ie?
		auto& air2 = air_column[i+1].ToMassFraction();
		auto& air1 = air_column[i].ToMassFraction();
		if (GetTheta(air2) - GetTheta(air1) < -1e-3) {
			il = i;
			break;
		}	
	}

	// if there is nowhere theta decreases with height, end the function
	if (il == -1) {return;}

	// determine iu where theta starts to increase with height
	for (int i = il + 1; i <= ie; ++i) {
		if (i == ie) {
			iu = i;
			break;
		}

		auto& air2 = air_column[i+1].ToMassFraction();
		auto& air1 = air_column[i].ToMassFraction();
		if (GetTheta(air2) - GetTheta(air1) > 1e-3) {
			iu = i;
			break;
		}
	}

	// execute convective adjustment to the aircolumn segment 'il-iu' 
	if (il != -1) {
		convective_adjustment(air_column, pcoord, grav, k, j, il, iu);
	}

	// return the current adjusted air column 
	return search_convective_adjustment(air_column, pcoord, grav, k,j);
}


void convective_adjustment(std::vector<AirParcel>& air, Coordinates *pcoord, Real grav,
													 int k, int j, int il, int iu) {
  
  Real total_energy = 0.;
  Real total_mass = 0.;
  
  AirParcel* parcel;

  // sum the energy and mass of all air parcels that to be adjusted (conservation of energy and mass)
  for (int i = il; i <= iu; ++i) {
    parcel = &air[i];
    parcel->ToMassFraction();

		Real volume = pcoord->GetCellVolume(k, j, i);
    Real temp = parcel->w[IPR] / Rd / parcel->w[IDN]; // IDN is density, IPR is pressure
		                                                  // primitive
    total_mass += parcel->w[IDN] * volume;
    total_energy += parcel->w[IDN] * volume * (cp / gammad * temp + grav * pcoord->x1v(i));

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

  bool energy_conserved = false;
  bool mass_conserved = false;
  
  while (energy_conserved == false || mass_conserved == false) {
  	Real guess_temp_0 = guess_pres_0 / guess_rho_0 / Rd;
		
		Real new_total_energy = 0.;
		Real new_total_mass = 0.;
		for (int i = il; i <= iu; ++i) {
			parcel = &air[i];
			temp = guess_temp_0 - grav / cp * (pcoord->x1v(i) - pcoord->x1v(il));
			parcel->w[IPR] = guess_pres_0 * pow(temp / guess_temp_0, cp / Rd); 
			parcel->w[IDN] = parcel->w[IPR] / Rd / temp;
			
			new_total_mass += parcel->w[IDN] * volume; 
			new_total_energy += parcel->w[IDN] * volume * (cp / gammad * temp + grav * pcoord->x1v(i));
		}

		if (new_total_energy - total_energy > 5.) {
			guess_pres_0 -= 0.5;  // need control guess_temp_0 positive
														// how much delta P -> how much delta E
		} else if (new_total_energy - total_energy < -5.) {
			guess_pres_0 += 0.5;
		} else {
			energy_conserved = true;
		}

		if (new_total_mass - total_mass > 5.) {
			guess_rho_0 -= 0.5;   // need control guess_rho_0 positive
														// how much delta rho -> how much delta mass
		} else if (new_total_mass - total_mass < -5.) {
			guess_rho_0 += 0.5;
		} else {
			mass_conserved = true;
		}
  }
}



Real GetTheta(AirParcel const& air) {
	Real pres_ref = 1.0e5;
	Real pres = air.w[IPR];
	Real rho = air.w[IDN];
	Real temp = pres / Rd / rho; 
	
	Real theta = temp * pow(pres_ref / pres, Rd / cp); 
	return theta;	
}
