// C/C++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// opacity
#include "absorber_ck.hpp"

void HeliosCKPremix::LoadCoefficient(std::string fname, size_t bid) 
{
  std::ifstream file(fname);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + fname);
  }

  size_t num_bands;

  // temperature, pressure, band, g-points
  file >> len_[1] >> len_[0] >> num_bands >> len_[2];

  if (bid >= num_bands) {
    throw std::runtime_error("Band index out of range: " + std::to_string(bid));
  }

  axis_.resize(len_[0] + len_[1] + len_[2]);
  kcoeff_.resize(len_[0] * len_[1] * len_[2]);

  // temperature grid
  for (int i = 0; i < len_[1]; ++i) {
    file >> axis_[len_[0] + i];
  }

  // pressure grid
  for (int j = 0; j < len_[0]; ++j) {
    file >> axis_[j];
  }

  Real dummy, wmin, wmax;
  // band wavelength
  for (int b = 0; b <= num_bands; ++b) {
    if (b != bid) {
      file >> dummy;
    } else {
      file >> wmin >> wmax;
      break;
    }
  }

  // g-points and weights
  for (int g = 0; g < len_[2]; ++g) {
    Real gpoint;
    file >> gpoint >> dummy;
    axis_[len_[0] + len_[1] + g] = wmin + (wmax - wmin) * gpoint;
  }

  int n = 0;
  for (int i = 0; i < len_[0]; ++i)
    for (int j = 0; j < len_[1]; ++j)
      for (int b = num_bands - 1; b >= 0; --b) {
        if (b != bid) {
          file >> dummy;
        } else {
          for (int g = 0; g < len_[2]; ++g, ++n) {
            file >> kcoeff_[n];
            kcoeff_[n] = log(std::max(kcoeff_[n], 1.0e-99));
          }
        }
      }

  file.close();
}
