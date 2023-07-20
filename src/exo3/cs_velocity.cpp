// C/C++
#include <cmath>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

namespace CubedSphere {

void GetUV(Real *U, Real *V, Coordinates *pcoord, Real V2, Real V3, int k,
           int j, int i) {
  // Obtain U and V (Lat-Lon) from V2 and V3 (Gnomonic Equiangle)
  // U is V_lam, V is V_phi.

  // Find the block number
  LogicalLocation loc = pcoord->pmy_block->loc;
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
  Real X = tan(pcoord->x2v(j));
  Real Y = tan(pcoord->x3v(k));
  VecTransRLLFromABP(X, Y, blockID, V2, V3, U, V);
}

void GetVyVz(Real *V2, Real *V3, Coordinates *pcoord, Real U, Real V, int k,
             int j, int i) {
  // Convert U and V (Lat-Lon) to V2 and V3 (Gnomonic Equiangle)
  // U is V_lam, V is V_phi.

  // Find the block number
  LogicalLocation loc = pcoord->pmy_block->loc;
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
  // Calculate the needed parameters
  Real X = tan(pcoord->x2v(j));
  Real Y = tan(pcoord->x3v(k));
  VecTransABPFromRLL(X, Y, blockID, U, V, V2, V3);
}

}  // namespace CubedSphere
