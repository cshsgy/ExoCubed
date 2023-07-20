// C/C++
#include <iostream>
#include <ostream>
#include <sstream>

// athena
#include <athena/coordinates/coordinates.hpp>

// exo3
#include "cubed_sphere.hpp"

namespace CubedSphere {

Real MeshGeneratorX2(Real x, LogicalLocation const& loc) {
  Real x_l, x_u;
  int lx2_lv2 = loc.lx2 >> (loc.level - 2);
  switch (lx2_lv2) {
    case 0:
      x_l = -0.5;
      x_u = 0.0;
      break;
    case 1:
      x_l = 0.0;
      x_u = 0.5;
  }
  return (0.5 * (x - x_l) / (x_u - x_l) - 0.25) * PI;  // Add Pi back later!!
}

Real MeshGeneratorX3(Real x, LogicalLocation const& loc) {
  Real x_l, x_u;
  int lx3_lv2 = loc.lx3 >> (loc.level - 2);
  switch (lx3_lv2) {
    case 0:
      x_l = -0.5;
      x_u = -1.0 / 6.0;
      break;
    case 1:
      x_l = -1.0 / 6.0;
      x_u = 1.0 / 6.0;
      break;
    case 2:
      x_l = 1.0 / 6.0;
      x_u = 0.5;
  }
  return (0.5 * (x - x_l) / (x_u - x_l) - 0.25) * PI;  // Add Pi back later!!
}

}  // namespace CubedSphere
