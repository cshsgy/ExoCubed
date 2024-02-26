#include <cmath>

// athena
#include <athena/coordinates/coordinates.hpp>

// canoe
#include <air_parcel.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

//thermodynamics/saturation_adjustment.cpp  as reference



void convective_adjustment(std::vector<AirParcel>& air_column, Coordinates *pcoord, Real grav,
													 int k, int j, int il, int iu) {
	auto pthermo = Thermodynamics::GetInstance();
	Real Rd = pthermo->GetRd();
	Real gammad;
	Real cp;
	Real cv;
  
	Real temp, volume;
	
	// the pressure and density of the lower and upper parcel of the air column section
	Real pres_l = air_column[il].w[IPR];
	Real rho_l = air_column[il].w[IDN];
	Real pres_u = air_column[iu].w[IPR];
	Real rho_u = air_column[iu].w[IDN];


  Real total_energy = 0.;
  Real total_mass = 0.;
  
  AirParcel* parcel;

  // sum the energy and mass of all air parcels that to be adjusted (conservation of energy and mass)
  for (int i = il; i <= iu; ++i) {
    parcel = &air_column[i];
    parcel->ToMassFraction();

		gammad = pthermo->GetGammad(*parcel);
		cv = Rd / (gammad-1.);
		volume = pcoord->GetCellVolume(k, j, i);
    temp = parcel->w[IPR] / Rd / parcel->w[IDN]; // IDN is density, IPR is pressure
		                                                  // primitive
    total_mass += parcel->w[IDN] * volume;
    total_energy += parcel->w[IDN] * volume * (cv * temp + grav * pcoord->x1v(i));
		
  }
 
  // the initial guess of pressure and density of lower parcel
  Real guess_pres_0 = 0.5 * (pres_l + pres_u); 
  Real guess_rho_0 = 0.5 * (rho_l + rho_u);
  Real guess_temp_0;

	Real new_total_energy;
	Real new_total_mass;

	Real total_energy_bias;
  Real total_mass_bias;

	bool energy_conserved = false;
  bool mass_conserved = false;
  
	// adjust pressure and density, and examine mass and energy conservation
  while (energy_conserved == false || mass_conserved == false) {
  	guess_temp_0 = guess_pres_0 / guess_rho_0 / Rd;

		new_total_energy = 0.;
		new_total_mass = 0.;
		for (int i = il; i <= iu; ++i) {
			parcel = &air_column[i];
    	parcel->ToMassFraction();
			
			gammad = pthermo->GetGammad(*parcel);
			cp = Rd * gammad / (gammad-1.);
			volume = pcoord->GetCellVolume(k, j, i);
			temp = guess_temp_0 - grav / cp * (pcoord->x1v(i) - pcoord->x1v(il));
			parcel->w[IPR] = guess_pres_0 * pow(temp / guess_temp_0, cp / Rd); 
			parcel->w[IDN] = parcel->w[IPR] / Rd / temp;
			
			new_total_mass += parcel->w[IDN] * volume; 
			new_total_energy += parcel->w[IDN] * volume * (cp / gammad * temp + grav * pcoord->x1v(i));
		}

		Real volume_0 = pcoord->GetCellVolume(k, j, il); 
		
		// examine if energy conservation is satisfied
		total_energy_bias = new_total_energy - total_energy;
		if (fabs(total_energy_bias) > 5.) {
			Real dE_dP = (cp-Rd+grav*pcoord->x1v(il)/guess_temp_0) * volume_0 / Rd;
			Real pres_change = total_energy_bias / dE_dP;
			while (pres_change >= guess_pres_0) {pres_change *= 0.4;} // 
			guess_pres_0 -= pres_change;
		} else {energy_conserved = true;}

	  // examine if mass conservation is satisfied
		total_mass_bias = new_total_mass - total_mass;
		if (fabs(total_mass_bias) > 5.) {
			Real dm_drho = volume_0;
			Real rho_change = total_mass_bias / dm_drho;
			while (rho_change >= guess_rho_0) {rho_change *= 0.4;}
			guess_rho_0 -= rho_change;
		} else {mass_conserved = true;}
  }
}


Real GetTheta(AirParcel const& air) {
	auto pthermo = Thermodynamics::GetInstance();
	Real Rd = pthermo->GetRd();
	Real pres_ref = 1.0e5;
	Real pres = air.w[IPR];
	Real rho = air.w[IDN];
	Real temp = pres / Rd / rho; 
	
	Real gammad = pthermo->GetGammad(air);
	Real cp = Rd * gammad / (gammad-1.);
	Real theta = temp * pow(pres_ref / pres, Rd /cp); 
	return theta;	
}

// type of the function
void recursively_search_convective_adjustment(std::vector<AirParcel>& air_column, 
																							Coordinates *pcoord, Real grav, int k, int j) {
	
	size_t nlayers = air_column.size();

	// initialize il and iu, the lower and upper level between which the airparcels need
	// convective adjustment 
	int il = -1;
	int iu = nlayers-1;
	
	// determine il where theta starts to decrease with height	
	for (int i = 0; i < nlayers-1; ++i) {   // where to get is and ie?
		auto& air2 = air_column[i+1].ToMassFraction();
		auto& air1 = air_column[i].ToMassFraction();
		if (GetTheta(air2) - GetTheta(air1) < -1e-2) {
			il = i;
			break;
		}	
	}

	// if there is nowhere theta decreases with height, end the function
	if (il == -1) {
		std::cout << "no adjustment" << std::endl;
		return;}

	// determine iu where theta starts to increase with height
	for (int i = il + 1; i <= iu; ++i) {
		if (i == iu) {break;}

		auto& air2 = air_column[i+1].ToMassFraction();
		auto& air1 = air_column[i].ToMassFraction();
		if (GetTheta(air2) - GetTheta(air1) > 1e-2) {
			iu = i;
			break;
		}
	}

	// execute convective adjustment to the aircolumn segment 'il-iu' 
	if (il != -1) {
		convective_adjustment(air_column, pcoord, grav, k, j, il, iu);
	}

	// Recursive call to search further 
	return recursively_search_convective_adjustment(air_column, pcoord, grav, k,j);
}


