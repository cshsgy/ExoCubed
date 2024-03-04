// C/C++
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

// external
#include <gtest/gtest.h>
#include <application/application.hpp>

// opacity
#include <opacity/absorber_ck.hpp>

std::string data_folder = "ck_data_01242024/ck/";

TEST(LoadCoefficient, bid_is_10) {
  auto app = Application::GetInstance();
  std::string fname = "PM_ck_HELIOSK_cond_11_nOPT_wcia.txt";
  auto file = app->FindResource(data_folder + fname);
  HeliosCKPremix PM_ck_HELIOSK_cond_11_nOPT_wcia("PM_ck_HELIOSK_cond_11_nOPT_wcia");
  PM_ck_HELIOSK_cond_11_nOPT_wcia.LoadCoefficient(file, 10, std::cout);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
