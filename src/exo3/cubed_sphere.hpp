#ifndef SRC_EXO3_CUBED_SPHERE_HPP_
#define SRC_EXO3_CUBED_SPHERE_HPP_

// C/C++
#include <memory>

// athena
#include <athena/athena.hpp>

// canoe
#include <configure.hpp>

class MeshBlock;

class CubedSphere {
 public:
  CubedSphereLR(MeshBlock *pmb);

  void GetLatLon(Real *lat, Real *lon, int k, int j, int i) const;
  void GetLatLonFace2(Real *lat, Real *lon, int k, int j, int i) const;
  void GetLatLonFace3(Real *lat, Real *lon, int k, int j, int i) const;

  void GetUV(Real *U, Real *V, Real V2, Real V3, int k, int j, int i) const;
  void GetVyVz(Real *V2, Real *V3, Real U, Real V, int k, int j, int i) const;

  void TransformOx(int *ox2, int *ox3, int *tox2, int *tox3) const;

  void SaveLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);

  void LoadLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);

  void SynchronizeFluxesSend();

  void SynchronizeFluxesRecv();

 protected:
  void sendNeighborBlocks(int ox2, int ox3, int tg_rank, int tg_gid);
  void recvNeighborBlocks(int ox2, int ox3, int tg_rank, int tg_gid);

  AthenaArray<Real> L3DValues[3], R3DValues[3];

#ifdef MPI_PARALLEL
  MPI_Request send_request[4];
  MPI_Request recv_request[4];
#endif

  std::vector<Real> LRDataBuffer[4];

  MeshBlock *pmy_block_;
};

using CubedSpherePtr = std::shared_ptr<CubedSphere>;

#endif  // SRC_EXO3_CUBED_SPHERE_HPP_
