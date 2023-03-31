#ifndef MESHBLOCK_IMPL_HPP
#define MESHBLOCK_IMPL_HPP

// C/C++ header
#include <memory>

// Athena++ header
#include <mesh/mesh.hpp>

class BlockIndex;
class ParameterInput;
class Decomposition;
class ImplicitSolver;
class FaceReconstruct;
class Forcing;

//! \class CellVariables
//  \brief a collection of all physical data in a computational cell
class MeshBlock::Impl {
public:
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl();

  AthenaArray<Real> du;    // stores tendency

  BlockIndex            *pblock;
  Decomposition         *pdec;
  ImplicitSolver        *phevi;
  FaceReconstruct       *precon;
  Forcing               *pforce;

private:
  MeshBlock const *pmy_block_;
};

#endif
