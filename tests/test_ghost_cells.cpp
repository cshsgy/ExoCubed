
#include <parameter_input.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>
#include <configure.hpp>

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN,k,j,i) = (k-ks)*(je-js)*(ie-is) + (j-js)*(ie-is) + (i-is);
      }
}
