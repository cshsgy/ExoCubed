//C/C++
#include <algorithm>
#include <vector>
#include <cmath>

// external
#include <gtest/gtest.h>

// ExoCubed/src
#include <opacity/cia_read.hpp>

TEST(reform_read, test_case1) {
  std::vector <double> expected_result = {4.47996926605873e-45, 1.868712588877408e-44, 1.2799236924204295e-44, 4.467454212770098e-45, 1.2948567646502106e-45, 2.802009368955379e-46, 2.308620846855808e-47, 3.3960726642220515e-49, 1e-99, 1e-99, 1e-99};
  AthenaArray<Real> data = reform_read ("He-H_reform.cia");
  std::vector <double> result;
  for (int i = 0; i < 11; ++i) {
    result.push_back(data(3,i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(reform_read, test_case2) {
  std::vector <double> expected_result = {2.0649974391687292e-45, 1.504861456948198e-45, 9.0902488843771e-47, 5.624908034565523e-48, 4.523257975006434e-47, 3.004769764226695e-45, 6.423470305957112e-48, 2.647376247478312e-49, 1.9496282953781007e-51, 1e-99, 1e-99};
  AthenaArray<Real> data = reform_read ("H2-He_reform.cia");
  std::vector <double> result;
  for (int i = 0; i < 11; ++i) {
    result.push_back(data(9,i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(ff_read, test_case1) {
  std::vector <double> expected_result = {7.74e-2, 1.07e-1, 1.21e-1, 1.34e-1, 1.57e-1, 1.81e-1, 2.29e-1, 2.79e-1};
  AthenaArray<Real> data = ff_read ("He-_ff.txt");
  std::vector <double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(4,i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(ff_read, test_case2) {
  std::vector <double> expected_result = {8.70e-2, 1.24e-1, 1.46e-1, 1.67e-1, 2.10e-1, 2.53e-1, 3.39e-1, 4.27e-1};
  AthenaArray<Real> data = ff_read ("H2-_ff.txt");
  std::vector <double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(2,i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(ff_read, test_case3) {
  std::vector <double> expected_result = {2.80e1, 3.62e1, 4.08e1, 4.49e1, 5.26e1, 5.98e1, 7.27e1, 8.40e1};
  AthenaArray<Real> data = ff_read ("He-_ff.txt");
  std::vector <double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(15,i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(ff_read, test_case4) {
  std::vector <double> expected_result = {7.16e1, 9.23e1, 1.01e2, 1.08e2, 1.18e2, 1.26e2, 1.38e2, 1.47e2};
  AthenaArray<Real> data = ff_read ("H2-_ff.txt");
  std::vector <double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(17,i));
  }
  EXPECT_EQ(result, expected_result);
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
