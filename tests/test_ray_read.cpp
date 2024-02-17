// C/C++
#include <algorithm>
#include <cmath>
#include <vector>

// external
#include <gtest/gtest.h>

// ExoCubed/src
#include <opacity/ray_read.hpp>

TEST(read_file, test_case1) {
  std::vector<double> expected_result = {
      4.76343376712079e-28,   1.749309404604686e-28,  6.539309299693454e-29,
      1.5810255315582998e-29, 3.63245689802763e-30,   1.6897747786802582e-30,
      4.900547268018036e-31,  2.0761625566011635e-31, 1.5181456645496456e-32,
      5.806445837106188e-34,  8.774789618085861e-39};
  std::vector<double> result = read_ray("Ray_He_11.txt");
  EXPECT_EQ(result, expected_result);
}

TEST(read_file, test_case2) {
  std::vector<double> expected_result = {
      6.652458732213571e-25, 6.652458732213571e-25, 6.652458732213571e-25,
      6.652458732213571e-25, 6.652458732213571e-25, 6.652458732213571e-25,
      6.652458732213571e-25, 6.652458732213571e-25, 6.652458732213571e-25,
      6.652458732213571e-25, 6.652458732213571e-25};
  std::vector<double> result = read_ray("Ray_e-_11.txt");
  EXPECT_EQ(result, expected_result);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
