#include <cmath>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// athena
#include <athena/athena.hpp>

// The first array is the 2D data sheet of cross sections, the second vector is
// the temperature axis, and the third is spectral axis

// read cia reform format file on temperature vs. spectral wavelength 2D data
// set.
std::tuple<AthenaArray<Real>, std::vector<double>, std::vector<double>>
reform_read(std::string filename);

// read cia ff format file on temperature vs. spectral wavelength 2D data set.
std::tuple<AthenaArray<Real>, std::vector<double>, std::vector<double>> ff_read(
    std::string filename);
