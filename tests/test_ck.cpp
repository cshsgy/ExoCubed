// C/C++
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>

// external
#include <gtest/gtest.h>
#include <application/application.hpp>

// opacity
#include <opacity/absorber_ck.hpp>

std::string data_folder = "ck_data_01242024/cia/";

TEST(LoadCoefficient, bid_is_1) {
  std::string correct_output = 
  "37 34 11 8 200.0 300.0 400.0 500.0 600.0 700.0
  800.0 900.0 1000.0 1100.0 1200.0 1300.0 1400.0 1500.0
  1700.0 1900.0 2100.0 2300.0 2500.0 2700.0 2900.0 3100.0
  3300.0 3500.0 3700.0 3900.0 4100.0 4300.0 4500.0 4700.0
  4900.0 5100.0 5300.0 5500.0 5700.0 5900.0 6100.0 1e-08
  2.1877616239495518e-08 4.6773514128719815e-08 1e-07
  2.187761623949552e-07 4.677351412871981e-07 1e-06
  2.1877616239495517e-06 4.677351412871981e-06 1e-05
  2.187761623949552e-05 4.677351412871981e-05 0.0001
  0.00021877616239495518 0.00046773514128719813 0.001
  0.002187761623949552 0.004677351412871981 0.01
  0.02187761623949553 0.046773514128719815 0.1 
  0.218776162394955230.4677351412871982 1.0
  2.137962089502232 4.570881896148751 10.0
  21.379620895022324 45.708818961487495 100.0
  213.79620895022325 457.0881896148752 1000.0
  0.26 0.42";

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "PM_ck_HELIOSK_cond_11_nOPT_wcia.txt");
  HeliosCKPremix PM_ck_HELIOSK_cond_11_nOPT_wcia("PM_ck_HELIOSK_cond_11_nOPT_wcia");

  std::stringstream output;
  PM_ck_HELIOSK_cond_11_nOPT_wcia.LoadCoefficient(file, 1, output);

  EXPECT_EQ(output, correct_output);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
