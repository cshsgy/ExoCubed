#include <parameter_input.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>
#include <configure.hpp>
#include "../src/eos/eos.hpp"
#include "../src/field/field.hpp"

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::cout << "======Start of Problem generator======" << loc.lx2 <<std::endl;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN,k,j,i) = (k-ks)*(je-js)*(ie-is) + (j-js)*(ie-is) + (i-is) + 0.1*(loc.lx2+1) + 0.01*(loc.lx3+1);
        phydro->w(IPR,k,j,i) = 1.0;
      }
  std::cout << "Done with problem generator..." << std::endl;
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
