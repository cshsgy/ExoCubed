// utils
#include <cmath>
#include <iostream>
#include <fstream>
#include <others/vector_conadj.hpp>
#include <utils/fileio.hpp>  // read_data_vector
#include <air_parcel.hpp>
#include <cstdlib>   // for exit program
#include <random>

int main() {
	Real grav = 23.1;
	Real Rd = 3714.;
	Real Ps = 1.e7;
	Real Ts = 8000.;
	Real cp = 3.5*Rd;
	Real e = 2.718281828459045;
	Real temp;
	Real T_fluct;
	Real theta;

	AirParcel air(AirParcel::Type::MassFrac);
	// read air column height
	std::string input_atm_path =
	"/home/linfel/ExoCubedlinfel/data/convadj_test_atm2.txt";
	DataVector atm = read_data_vector(input_atm_path);

	Real* x1_ptr;
	x1_ptr = atm["HGT"].data();
	size_t nx1 = atm["HGT"].size();
  
	// prepare random real number
	std::random_device rd;
	std::mt19937 gen(rd()); // Mersenne Twister generator
	// Define the distribution to be uniform for real numbers between 0.0 and 1.0
	std::uniform_real_distribution<> dis(-20., 20.0);

	// output file
	std::ofstream outFile1("ac_before.csv");
	outFile1 << "i,T_fluct,x1,pres,temp,rho,theta" << std::endl;

	// construt the air column
	std::cout << "Constructing air column" << std::endl;
	std::vector<AirParcel> vector_ac;
	for (int i = 0; i < nx1; ++i) {
		air.SetZero();
		
		T_fluct = dis(gen);
		temp = Ts - grav / cp * (x1_ptr[i]-x1_ptr[0]) + T_fluct;
		air.w[IPR] = Ps * pow(e, -grav/Rd/temp*(x1_ptr[i]-x1_ptr[0]));
		air.w[IDN] = air.w[IPR] / Rd / temp;
		theta = GetTheta(air);
		std::cout << i << "   " << T_fluct << "  " << x1_ptr[i] << "   " << air.w[IPR] << "   " <<  temp << "   " << air.w[IDN] << std::endl;
		
		outFile1 << i << "," << T_fluct << "," << x1_ptr[i] << "," << air.w[IPR] << "," <<
						temp << "," << air.w[IDN] << "," << theta << std::endl;
			

		vector_ac.push_back(air);
	}
	
	outFile1.close();

	std::cout << "=======================" << std::endl;
	
	recursively_search_convective_adjustment(vector_ac, x1_ptr, grav, 0, 0);

	// output air column after convective adjustment
	std::ofstream outFile2("ac_after.csv");
	outFile2 << "i,x1,pres,rho,theta" << std::endl;
	for (int i = 0; i < nx1; ++i) {
		outFile2 << i << "," << x1_ptr[i] << "," << vector_ac[i].w[IPR] << "," <<
				 vector_ac[i].w[IDN] << "," << GetTheta(vector_ac[i]) << std::endl;
	}
	outFile2.close();

	return 0;
}





