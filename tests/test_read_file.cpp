//C/C++
#include <algorithm>
#include <vector>
#include <cmath>

// external
#include <gtest/gtest.h>

// ExoCubed
#include <opacity/read_file.hpp>

TEST(read_file, test_case1) {
  std::vector <double> output1 = (1059434.8118773764, 6035941.569908597, 7520131.630087908, 
                                  8488784.244992392, 5729875.721026806, 1400979.114668621,
                                  1083717.0198599927, 331262.2156410236, 296581.7003886011,
                                  44557.440115296864, 4298.3275481968085);
  std::vector <double> output2 = (12, 0.260, 0.420, 0.610, 0.850, 1.320, 2.020, 
                                  2.500, 3.500, 4.400, 8.70, 20.00, 324.68);
  std::pair <vector <double>, vector <double>> expected_result (output1, output2);
  std::pair <vector <double>, vector <double>> result = read_file (sw_band_flux_HD189_11.txt, wavelengths_GCM_11.txt);
  EXPECT_EQ(result, expected_result);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
