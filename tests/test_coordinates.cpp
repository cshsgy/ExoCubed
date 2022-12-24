#include <parameter_input.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>
#include <configure.hpp>
#include <cubed_sphere.hpp>

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::cout << "my rank  = " << Globals::my_rank << std::endl;
  std::cout << "my panel = " << FindBlockID(loc) << std::endl;

  std::cout << "x1v = " << std::endl;
  for (int i = is - GhostZoneSize; i < is; ++i) {
    std::cout << pcoord->x1v(i)  << " " << std::endl;
  }
  std::cout << "--------------" << std::endl;

  for (int i = is; i <= ie; ++i) {
    std::cout << pcoord->x1v(i)  << " " << std::endl;
  }
  std::cout << "--------------" << std::endl;
  for (int i = ie + 1; i <= ie + GhostZoneSize; ++i) {
    std::cout << pcoord->x1v(i)  << " " << std::endl;
  }
  std::cout << std::endl;

  std::cout << "x2f = " << std::endl;
  for (int j = js - GhostZoneSize; j < js; ++j) {
    std::cout << pcoord->x2f(j)/M_PI  << " " << std::endl;
  }
  std::cout << "--------------" << std::endl;

  for (int j = js; j <= je; ++j) {
    std::cout << pcoord->x2f(j)/M_PI  << " " << std::endl;
  }
  std::cout << "--------------" << std::endl;

  for (int j = je + 1; j <= je + GhostZoneSize; ++j) {
    std::cout << pcoord->x2f(j)/M_PI  << " " << std::endl;
  }
  std::cout << std::endl;

  std::cout << "x3f = " << std::endl;
  for (int k = ks - GhostZoneSize; k < ks; ++k) {
    std::cout << pcoord->x3f(k)/M_PI  << " " << std::endl;
  }
  std::cout << "--------------" << std::endl;

  for (int k = ks; k <= ke; ++k) {
    std::cout << pcoord->x3f(k)/M_PI  << " " << std::endl;
  }
  std::cout << "--------------" << std::endl;

  for (int k = ke + 1; k <= ke + GhostZoneSize; ++k) {
    std::cout << pcoord->x3f(k)/M_PI  << " " << std::endl;
  }
}
