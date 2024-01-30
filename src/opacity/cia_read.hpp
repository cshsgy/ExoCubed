#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <utility>

// athena
#include <athena/athena.hpp>

//read cia reform format file on temperature vs. spectral wavelength 2D data set. 
AthenaArray<Real> reform_read (std::string filename);

//read cia ff format file on temperature vs. spectral wavelength 2D data set. 
AthenaArray<Real> ff_read (std::string filename);