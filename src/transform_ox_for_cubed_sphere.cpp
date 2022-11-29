#include <sstream>
#include <ostream>
#include <bvals/bvals.hpp>
#include <bvals/cubed_sphere.hpp>

// Add the transformation function into the boundary base class, for two reasons:
// 1. We would anyway need to change the BoundaryBase class
// 2. We need to be in the friend class for meshblock tree (well maybe I am wrong...)
void TransformOxForCubedSphere(int *ox1, int *ox2, int *tox1, int *tox2,
  LogicalLocation const& loc)
{
  // Find the location in level 2
  int lv2_lx1 = loc.lx1 >> (loc.level - 2);
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  
  // Determine the block number, works for 6 blocks only
  int block_id = ((lv2_lx1>>1)*4+(lv2_lx1-((lv2_lx1>>1)<<1))+lv2_lx2*(2-(lv2_lx1>>1))+1);
  
  // Find relative location within block
  int local_lx1 = loc.lx1 - (lv2_lx1<<(loc.level - 2));
  int local_lx2 = loc.lx2 - (lv2_lx2<<(loc.level - 2));
  int bound_lim = (1<<(loc.level - 2)) - 1;

  // Hard code the cases...
  // No need to consider the corner cases, abandon in reading buffers.
  int target_block = -1; // Block id of target
  int target_loc_x, target_loc_y; // local x and y in target block

  switch (block_id)
  {
    case 1: // Special: left, top, right
      if ((local_lx1==0) && (*ox1==-1)){ // Left Boundary
        target_block = 2;
        target_loc_x = local_lx2;
        target_loc_y = 0;
        *tox1 = 0;
        *tox2 = -1;
      }
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 6;
        target_loc_x = bound_lim - local_lx1;
        target_loc_y = 0;
        *tox1 = 0;
        *tox2 = -1;
      }
      if ((local_lx1==bound_lim) && (*ox1==1)){ // Right Boundary
        target_block = 4;
        target_loc_x = bound_lim - local_lx2;
        target_loc_y = 0;
        *tox1 = 0;
        *tox2 = -1;
      }
      break;
    case 2:
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 1;
        target_loc_x = 0;
        target_loc_y = local_lx1;
        *tox1 = -1;
        *tox2 = 0;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 5;
        target_loc_x = 0;
        target_loc_y = bound_lim - local_lx1;
        *tox1 = -1;
        *tox2 = 0;
      }
      if ((local_lx1==bound_lim) && (*ox1==1)){ // Right Boundary
        target_block = 3;
        target_loc_x = 0;
        target_loc_y = local_lx2;
        *tox1 = 1;
        *tox2 = 0;
      }
      break;
    case 3:
      if ((local_lx1==0) && (*ox1==-1)){ // Left Boundary
        target_block = 2;
        target_loc_x = bound_lim;
        target_loc_y = local_lx2;
        *tox1 = 1;
        *tox2 = 0;
      }
      break;
    case 4:
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 1;
        target_loc_x = bound_lim;
        target_loc_y = bound_lim - local_lx1;
        *tox1 = 1;
        *tox2 = 0;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 5;
        target_loc_x = bound_lim;
        target_loc_y = local_lx1;
        *tox1 = 1;
        *tox2 = 0;
      }
      break;
    case 5:
      if ((local_lx1==0) && (*ox1==-1)){ // Left Boundary
        target_block = 2;
        target_loc_x = bound_lim - local_lx2;
        target_loc_y = bound_lim;
        *tox1 = 0;
        *tox2 = 1;
      }
      if ((local_lx1==bound_lim) && (*ox1==1)){ // Right Boundary
        target_block = 4;
        target_loc_x = local_lx2;
        target_loc_y = bound_lim;
        *tox1 = 0;
        *tox2 = 1;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 6;
        target_loc_x = bound_lim - local_lx1;
        target_loc_y = bound_lim;
        *tox1 = 0;
        *tox2 = 1;
      }
      break;
    case 6:
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 1;
        target_loc_x = bound_lim - local_lx1;
        target_loc_y = 0;
        *tox1 = 0;
        *tox2 = -1;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 5;
        target_loc_x = bound_lim - local_lx1;
        target_loc_y = bound_lim;
        *tox1 = 0;
        *tox2 = 1;
      }
      break;
    default:
      std::stringstream msg;
      msg << "Error: something wrong, check the geometry setup of the cubed sphere. \n";
      msg << "----------------------------------" << std::endl;
      ATHENA_ERROR(msg);
  }

  // Calculate ox1 and ox2
  if (target_block>0){ // Need to change only when a special boundary is crossed
    // Calculate the lx1 and lx2 positions of the neighbor block
    // First calculate the top left corner position of the block
    int lx1_0, lx2_0;
    switch (target_block)
    {
    case 1:
      lx1_0 = 0;
      lx2_0 = 0;
      break;
    case 2:
      lx1_0 = 1;
      lx2_0 = 0;
      break;
    case 3:
      lx1_0 = 0;
      lx2_0 = 1;
      break;
    case 4:
      lx1_0 = 1;
      lx2_0 = 1;
      break;
    case 5:
      lx1_0 = 2;
      lx2_0 = 0;
      break;
    case 6:
      lx1_0 = 2;
      lx2_0 = 1;
      break;
    default:
      std::stringstream msg;
      msg << "Error: something wrong, check the geometry setup of the cubed sphere. \n";
      msg << "----------------------------------" << std::endl;
      ATHENA_ERROR(msg);
    }
    lx1_0 = (lx1_0 << (loc.level - 2));
    lx2_0 = (lx2_0 << (loc.level - 2));
    // Add up first block and local positions
    int lx1_t = lx1_0 + target_loc_x;
    int lx2_t = lx2_0 + target_loc_y;
    // Calculate and pass the differences
    *ox1 = lx1_t - loc.lx1;
    *ox2 = lx2_t - loc.lx2;
  }
  return;
}
