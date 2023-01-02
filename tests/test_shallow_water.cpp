#include <parameter_input.hpp>
#include <hydro/hydro.hpp>
#include <coordinates/coordinates.hpp>
#include <mesh/mesh.hpp>
#include <configure.hpp>

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real phi = pin->GetReal("problem", "phi");
  Real uphi = pin->GetReal("problem", "uphi");
  Real dphi = pin->GetReal("problem", "dphi");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      phydro->u(IDN,k,j,is) = phi;
      if (pcoord->x3v(k) > 0.)
        phydro->u(IM2,k,j,is) = -uphi;
      else
        phydro->u(IM2,k,j,is) = uphi;

      if (pcoord->x2v(j) > 0. && pcoord->x2v(j) < 5.)
        phydro->u(IDN,k,j,is) += dphi;
    }
}
