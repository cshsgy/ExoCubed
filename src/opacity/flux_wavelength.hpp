#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <utility>

//read two textfiles and combine the list of data into a pair of vector of double specifically sw_band_flux_HD189_11.txt and wavelengths_GCM_11.txt
std::pair <std::vector <double>, std::vector <double>> read_file (std::string file1, std::string file2);