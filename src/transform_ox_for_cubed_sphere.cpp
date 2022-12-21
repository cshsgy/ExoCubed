#include <sstream>
#include <ostream>
#include <bvals/bvals.hpp>
#include <bvals/cubed_sphere.hpp>

int FindBlockID(LogicalLocation const& loc){
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  
  // Determine the block number
  int block_id;
  switch(lv2_lx3)
  {
    case 0:
      switch(lv2_lx2)
      {
        case 0:
          block_id = 1;
          break;
        case 1:
          block_id = 3;
          break;
        default:
          std::stringstream msg;
          msg << "Error: something wrong, check the geometry setup of the cubed sphere. \n";
          msg << "----------------------------------" << std::endl;
          ATHENA_ERROR(msg);
      }
      break;
    case 1:
      switch(lv2_lx2)
      {
        case 0:
          block_id = 2;
          break;
        case 1:
          block_id = 4;
          break;
        default:
          std::stringstream msg;
          msg << "Error: something wrong, check the geometry setup of the cubed sphere. \n";
          msg << "----------------------------------" << std::endl;
          ATHENA_ERROR(msg);
      }
      break;
    case 2:
      switch(lv2_lx2)
      {
        case 0:
          block_id = 5;
          break;
        case 1:
          block_id = 6;
          break;
        default:
          std::stringstream msg;
          msg << "Error: something wrong, check the geometry setup of the cubed sphere. \n";
          msg << "----------------------------------" << std::endl;
          ATHENA_ERROR(msg);
      }
    break;
    default:
      std::stringstream msg;
      msg << "Error: something wrong, check the geometry setup of the cubed sphere. \n";
      msg << "----------------------------------" << std::endl;
      ATHENA_ERROR(msg);
  }

  return block_id;
}

void TransformOxForCubedSphere(int *ox2, int *ox3, int *tox2, int *tox3,
  LogicalLocation const& loc)
{
  // Find the block ID
  int block_id = FindBlockID(loc);

  // Find relative location within block
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx2 = loc.lx2 - (lv2_lx2<<(loc.level - 2));
  int local_lx3 = loc.lx3 - (lv2_lx3<<(loc.level - 2));
  int bound_lim = (1<<(loc.level - 2)) - 1;

  // Hard code the cases...
  // No need to consider the corner cases, abandon in reading buffers.
  int target_block = -1; // Block id of target
  int target_loc_2, target_loc_3; // local x2 and x3 in target block

  switch (block_id)
  {
    case 1: // Special: left, top, right
      if ((local_lx3==0) && (*ox3==-1)){ // Left Boundary
        target_block = 2;
        target_loc_2 = 0;
        target_loc_3 = local_lx2;
        *tox2 = -1;
        *tox3 = 0;
      }
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 6;
        target_loc_2 = 0;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = -1;
        *tox3 = 0;
      }
      if ((local_lx3==bound_lim) && (*ox3==1)){ // Right Boundary
        target_block = 4;
        target_loc_2 = 0;
        target_loc_3 = bound_lim - local_lx2;
        *tox2 = -1;
        *tox3 = 0;
      }
      break;
    case 2:
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 1;
        target_loc_2 = local_lx3;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 5;
        target_loc_2 = bound_lim - local_lx3;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
      }
      if ((local_lx3==bound_lim) && (*ox3==1)){ // Right Boundary
        target_block = 3;
        target_loc_2 = local_lx2;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
      }
      break;
    case 3:
      if ((local_lx3==0) && (*ox3==-1)){ // Left Boundary
        target_block = 2;
        target_loc_2 = local_lx2;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      break;
    case 4:
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 1;
        target_loc_2 = bound_lim - local_lx3;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 5;
        target_loc_2 = local_lx3;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      break;
    case 5:
      if ((local_lx3==0) && (*ox3==-1)){ // Left Boundary
        target_block = 2;
        target_loc_2 = bound_lim;
        target_loc_3 = bound_lim - local_lx2;
        *tox2 = 1;
        *tox3 = 0;
      }
      if ((local_lx3==bound_lim) && (*ox3==1)){ // Right Boundary
        target_block = 4;
        target_loc_2 = bound_lim;
        target_loc_3 = local_lx2;
        *tox2 = 1;
        *tox3 = 0;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 6;
        target_loc_2 = bound_lim;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = 1;
        *tox3 = 0;
      }
      break;
    case 6:
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 1;
        target_loc_2 = 0;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = -1;
        *tox3 = 0;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 5;
        target_loc_2 = bound_lim;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = 1;
        *tox3 = 0;
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
    int lx3_0, lx2_0;
    switch (target_block)
    {
    case 1:
      lx3_0 = 0;
      lx2_0 = 0;
      break;
    case 2:
      lx3_0 = 1;
      lx2_0 = 0;
      break;
    case 3:
      lx3_0 = 0;
      lx2_0 = 1;
      break;
    case 4:
      lx3_0 = 1;
      lx2_0 = 1;
      break;
    case 5:
      lx3_0 = 2;
      lx2_0 = 0;
      break;
    case 6:
      lx3_0 = 2;
      lx2_0 = 1;
      break;
    default:
      std::stringstream msg;
      msg << "Error: something wrong, check the geometry setup of the cubed sphere. \n";
      msg << "----------------------------------" << std::endl;
      ATHENA_ERROR(msg);
    }
    lx3_0 = (lx3_0 << (loc.level - 2));
    lx2_0 = (lx2_0 << (loc.level - 2));
    // Add up first block and local positions
    int lx3_t = lx3_0 + target_loc_3;
    int lx2_t = lx2_0 + target_loc_2;
    // Calculate and pass the differences
    *ox3 = lx3_t - loc.lx1;
    *ox2 = lx2_t - loc.lx2;
  }
  return;
}

void PackDataCubedSphereR3(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; n++){
    for (int j=sj; j<=ej; j++) {
      for (int k=ek; k>=sk; k--) {
#pragma omp simd
        for (int i=si; i<=ei; i++)
          buf[offset++] = src(k, j, i);
      }
    }
  }
  return;
}

void PackDataCubedSphereR2(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; n++){
    for (int k=ek; k>=sk; k--) {
      for (int j=ej; j>=sj; j--) {
#pragma omp simd
        for (int i=si; i<=ei; i++)
          buf[offset++] = src(k, j, i);
      }
    }
  }
  return;
}

void PackDataCubedSphereR1(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; n++){
    for (int j=ej; j>=sj; j--) {
      for (int k=sk; k<=ek; k++) {
#pragma omp simd
        for (int i=si; i<=ei; i++)
          buf[offset++] = src(k, j, i);
      }
    }
  }
  return;
}


void PackDataCubedSphere(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset,
         int ox2, int ox3,LogicalLocation const& loc){
// Find the block ID
int blockID = FindBlockID(loc);
switch(blockID){
  case 1:
    if(ox3==1){PackDataCubedSphereR1(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    if(ox2==-1){PackDataCubedSphereR2(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    if(ox3==-1){PackDataCubedSphereR3(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    break;
  case 2:
    if(ox2==1){PackDataCubedSphereR3(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    if(ox2==-1){PackDataCubedSphereR1(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    break;
  case 3:
    break;
  case 4:
    if(ox2==1){PackDataCubedSphereR1(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    if(ox2==-1){PackDataCubedSphereR3(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    break;
  case 5:
    if(ox3==1){PackDataCubedSphereR3(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    if(ox2==1){PackDataCubedSphereR2(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    if(ox3==-1){PackDataCubedSphereR1(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    break;
  case 6:
    if(ox2==1){PackDataCubedSphereR2(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
    if(ox2==-1){PackDataCubedSphereR2(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);}
}
return;
}