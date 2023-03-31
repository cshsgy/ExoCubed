// Athena++ headers
#include <athena.hpp>
#include <parameter_input.hpp>
#include <coordinates/coordinates.hpp>
#include <hydro/hydro.hpp>
#include <field/field.hpp>
#include <scalars/scalars.hpp>
#include <orbital_advection/orbital_advection.hpp>
#include <gravity/gravity.hpp>
#include <eos/eos.hpp>

// exdacs headers
#include <configure.hpp>
#include "block_index.hpp"
#include "meshblock_impl.hpp"
#include "../hydro/implicit/implicit_solver.hpp"
#include "../hydro/decomposition/decomposition.hpp"
#include "../reconstruct/face_reconstruct.hpp"

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin):
  pmy_block_(pmb)
{
  du.NewAthenaArray(NumHydros, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // block index
  pblock = new BlockIndex(pmb);

  // hydro decomposition
  pdec = new Decomposition(pmb->phydro);

  // implicit methods
  phevi = new ImplicitSolver(pmb, pin);

  // reconstruction
  precon = new FaceReconstruct(pmb, pin);
}

// Athena++ demands destruct pbval AFTER all boundary values
// But in this mod, boundary values are destructed BEFORE pbval
// TODO, check if this is OK
MeshBlock::Impl::~Impl()
{
  //std::cout << "Impl desctructor" << std::endl;
  delete pblock;
  delete pdec;
  delete phevi;
  delete precon;
}
