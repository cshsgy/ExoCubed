// C/C++
#include <algorithm>
#include <cmath>
#include <iostream>

#include <fstream>
#include <vector>
#include <string>


// external
#include <gtest/gtest.h>

// athena
#include <athena/reconstruct/interpolation.hpp>

TEST(interp_weno3, test_case1) {

    std::string fileName = "../../data/CK_data/sw_band_flux_HD189_11.txt";
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << fileName << std::endl;
    }

    std::vector<double> numbers;
    std::string line;
    while (std::getline(file, line)) {
        try {
            double num = std::stod(line);
            numbers.push_back(num);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid number found: " << line << std::endl;
        }
    }

    file.close();
}


TEST(interp_weno3, test_case2) {

    std::string fileName = "../../data/CK_data/wavelengths_GCM_11.txt";
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << fileName << std::endl;
    }

    std::vector<double> numbers;
    std::string line;
    while (std::getline(file, line)) {
        try {
            double num = std::stod(line);
            numbers.push_back(num);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid number found: " << line << std::endl;
        }
    }

    file.close();
}



int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
