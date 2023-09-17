// athena
#include <athena/mesh/mesh.hpp>

// canoe
#include <impl.hpp>

void Mesh::SaveAllStates()
{
  for (int b = 0; b < nblocal; ++b) {
    my_blocks(b)->pimpl->SaveAllStates();
  }
}

void Mesh::LoadAllStates()
{
  for (int b = 0; b < nblocal; ++b) {
    my_blocks(b)->pimpl->LoadAllStates();
  }
}

void Mesh::DecreaseTimeStep()
{
  dt /= 2;
}

bool Mesh::CheckAllValid() const
{
  bool valid = true;
  for (int b = 0; b < nblocal; ++b)
    if (!my_blocks(b)->pimpl->IsStateValid())
      valid = false;

  return true;
}
