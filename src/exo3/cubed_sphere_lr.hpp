#ifndef SRC_EXO3_CUBED_SPHERE_LR_HPP_
#define SRC_EXO3_CUBED_SPHERE_LR_HPP_

// athena
#include <athena/athena.hpp>

class MeshBlock;

class CubedSphereLR {
 public:
  void SetMeshBlock(MeshBlock *pmb_in);
  void InitializeSizes(int nc3, int nc2, int nc1);
  void SaveLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);
  void LoadLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);
  void SynchronizeFluxes();
  void SendNeighborBlocks(LogicalLocation const &loc, int ox2, int ox3,
                          int tg_rank, int tg_gid);
  void RecvNeighborBlocks(LogicalLocation const &loc, int ox2, int ox3,
                          int tg_rank, int tg_gid);

  AthenaArray<Real> L3DValues, R3DValues;
  MeshBlock *pmb;
};

#endif  // SRC_EXO3_CUBED_SPHERE_LR_HPP_
