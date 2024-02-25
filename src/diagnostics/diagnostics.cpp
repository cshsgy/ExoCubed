// C/C++ headers
#include <cstring>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/debugger/debugger.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

#include "diagnostics.hpp"

const std::string Diagnostics::input_key = "diagnostics";

Diagnostics::Diagnostics(MeshBlock *pmb, std::string short_name,
                         std::string long_name)
    : NamedGroup(short_name, short_name, long_name),
      pmy_block_(pmb),
      ncycle_(0) {
  auto app = App

      std::stringstream msg;
  ncells1_ = pmb->block_size.nx1 + 2 * (NGHOST);
  ncells2_ = 1;
  ncells3_ = 1;
  if (pmb->pmy_mesh->f2) ncells2_ = pmb->block_size.nx2 + 2 * (NGHOST);
  if (pmb->pmy_mesh->f3) ncells3_ = pmb->block_size.nx3 + 2 * (NGHOST);

  x1edge_.NewAthenaArray(ncells1_ + 1);
  x1edge_p1_.NewAthenaArray(ncells1_);
  x2edge_.NewAthenaArray(ncells1_ + 1);
  x2edge_p1_.NewAthenaArray(ncells1_);
  x3edge_.NewAthenaArray(ncells1_ + 1);
  x3edge_p1_.NewAthenaArray(ncells1_);

  x1area_.NewAthenaArray(ncells1_ + 1);
  x2area_.NewAthenaArray(ncells1_);
  x2area_p1_.NewAthenaArray(ncells1_);
  x3area_.NewAthenaArray(ncells1_);
  x3area_p1_.NewAthenaArray(ncells1_);

  vol_.NewAthenaArray(ncells1_);
  total_vol_.NewAthenaArray(ncells1_);
  brank_.resize(Globals::nranks);
  color_.resize(Globals::nranks);

  for (int i = 0; i < Globals::nranks; ++i) {
    brank_[i] = -1;
    color_[i] = -1;
  }
}
