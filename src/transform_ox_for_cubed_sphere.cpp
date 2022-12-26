#include <sstream>
#include <ostream>
#include <iostream>
#include <bvals/bvals.hpp>
#include <cubed_sphere.hpp>

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
          block_id = 2;
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
          block_id = 3;
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
  // std::cout << "|Block ID: " << block_id << "||ox2: " << *ox2 << "|ox3: " << *ox3 << std::endl;

  switch (block_id)
  {
    case 1: // Special: left, top, right
      if ((local_lx3==0) && (*ox3==-1)){ // Left Boundary
        target_block = 3;
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
      if ((local_lx3==0) && (*ox3==-1)){ // Left Boundary
        target_block = 3;
        target_loc_2 = local_lx2;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      if ((local_lx2==bound_lim) && (*ox2==1)){ // Bottom Boundary
        target_block = 5;
        target_loc_2 = 0;
        target_loc_3 = local_lx3;
        *tox2 = -1;
        *tox3 = 0;
      }
      break;
    case 3:
      if ((local_lx3==0) && (*ox3==-1)){ // Left Boundary
        target_block = 6;
        target_loc_2 = local_lx2;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
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
        target_block = 2;
        target_loc_2 = local_lx2;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
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
        target_block = 3;
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
      if ((local_lx2==0) && (*ox2==-1)){ // Top Boundary
        target_block = 2;
        target_loc_2 = bound_lim;
        target_loc_3 = local_lx3;
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
      if ((local_lx3==bound_lim) && (*ox3==1)){ // Right Boundary
        target_block = 3;
        target_loc_2 = local_lx2;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
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
      lx3_0 = 0;
      lx2_0 = 1;
      break;
    case 3:
      lx3_0 = 1;
      lx2_0 = 0;
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
    *ox3 = lx3_t - loc.lx3;
    *ox2 = lx2_t - loc.lx2;
  }else{
    *tox2 = -*ox2;
    *tox3 = -*ox3;
  }
  // std::cout << "|Block ID: " << block_id << "|Target Block: " << target_block << "||ox2: " << *ox2 << "|ox3: " << *ox3 << "||tox2:" << *tox2 << "|tox3:" << *tox3 << std::endl;
  return;
}

void PackDataCubedSphereR3(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; n++){
    for (int j=sj; j<=ej; j++) {
      for (int k=ek; k>=sk; k--) {
#pragma omp simd
        for (int i=si; i<=ei; i++)
          buf[offset++] = src(n, k, j, i);
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
          buf[offset++] = src(n, k, j, i);
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
          buf[offset++] = src(n, k, j, i);
      }
    }
  }
  return;
}

void PackDataCubedSphereR0(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; n++){
    for (int k=sk; k<=ek; k++) {
      for (int j=sj; j<=ej; j++) {
#pragma omp simd
        for (int i=si; i<=ei; i++)
          buf[offset++] = src(n, k, j, i);
      }
    }
  }
  return;
}

Real CalculateInterpLocationsCubedSphere(int loc_n, int N_blk, int k, bool GhostZone){
  Real temp0 = tan((Real)1.0*(N_blk-1-2*k)/(N_blk*4)*PI);
  if(GhostZone){
    temp0 = 1.0/temp0;
    temp0 = 1.0 + temp0 * temp0;
  }else{
    temp0 = 1.0 + temp0*temp0;
  }
  if (loc_n<0)
    return -acos(sqrt(temp0/(temp0 + pow(tan((Real)1.0*(2*loc_n+1)/(N_blk*4)*PI),2.0))));
  else
    return acos(sqrt(temp0/(temp0 + pow(tan((Real)1.0*(2*loc_n+1)/(N_blk*4)*PI),2.0))));
}

void InteprolateX2CubedSphere(const AthenaArray<Real> &src, AthenaArray<Real> &tgt, LogicalLocation const& loc, int sn, int en, int si, int ei, int sj, int ej, int sk, int ek){
  // Interpolation along X2 (j) axis, used before sending data to X3 (k) axis
  // Get the local indices
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx2 = loc.lx2 - (lv2_lx2<<(loc.level - 2));
  int bound_lim = (1<<(loc.level - 2)) - 1;
  int N_blk = (ej - sj + 1) * (bound_lim+1); // N in X2 direction for each panel. This value is required to be an even number
  int n_start = local_lx2 * bound_lim - N_blk / 2;
  Real *src_x2 = new Real[ej-sj+1];
  Real *tgt_x2 = new Real[ej-sj+1];

  for (int n=sn; n<=en; n++){
    for (int k=sk; k<=ek; k++) {
      int k_now;
      if(sk<ej-NGHOST) // Calculate the location into the ghost zones
        k_now = k;
      else
        k_now = ek-k;
#pragma omp simd
        for (int i=si; i<=ei; i++){
          for (int j=sj; j<=ej; j++) {
            // Calculate coefficients
            src_x2[j-sj] = CalculateInterpLocationsCubedSphere(n_start+j-sj, N_blk, k_now, false);
            tgt_x2[j-sj] = CalculateInterpLocationsCubedSphere(n_start+j-sj, N_blk, k_now, true);
          }
          int src_pointer = 0;
          for (int j=sj; j<=ej; j++) {
            // Interpolate to target array, linear interpolation used here
            while(tgt_x2[j]>src_x2[src_pointer+1]) src_pointer++; // Find the left location of src
            Real y1 = src_x2[src_pointer];
            Real y2 = src_x2[src_pointer+1];
            Real v1 = src(n, k, sj+src_pointer, i);
            Real v2 = src(n, k, sj+src_pointer+1, i);
            Real yq = tgt_x2[j];
            tgt(n-sn, k-sk, j-sj, i-si) = ((y2-yq)*v1 + (yq-y1)*v2) / (y2-y1);
          }
        }
    }
  }
  delete[] src_x2, tgt_x2; // Release memory
}

void InteprolateX3CubedSphere(const AthenaArray<Real> &src, AthenaArray<Real> &tgt, LogicalLocation const& loc, int sn, int en, int si, int ei, int sj, int ej, int sk, int ek){
  // Interpolation along X2 (j) axis, used before sending data to X3 (k) axis
  // Get the local indices
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx3 = loc.lx3 - (lv2_lx3<<(loc.level - 2));
  int bound_lim = (1<<(loc.level - 2)) - 1;
  int N_blk = (ek - sk + 1) * (bound_lim+1); // N in X2 direction for each panel. This value is required to be an even number
  int n_start = local_lx3 * bound_lim - N_blk / 2;
  Real *src_x3 = new Real[ek-sk+1];
  Real *tgt_x3 = new Real[ek-sk+1];

  for (int n=sn; n<=en; n++){
    for (int j=sj; j<=ej; j++) {
      int k_now;
      if(sj<ek-NGHOST) // Calculate the location into the ghost zones
        k_now = j;
      else
        k_now = ej-j;
#pragma omp simd
        for (int i=si; i<=ei; i++){
          for (int k=sk; k<=ek; k++) {
            // Calculate coefficients
            src_x3[k-sk] = CalculateInterpLocationsCubedSphere(n_start+k-sk, N_blk, k_now, false);
            tgt_x3[k-sk] = CalculateInterpLocationsCubedSphere(n_start+k-sk, N_blk, k_now, true);
            // std::cout << "|" << src_x3[k-sk] << "|" << tgt_x3[k-sk] << std::endl;
          }
          int src_pointer = 0;
          for (int k=sk; k<=ek; k++) {
            // Interpolate to target array, linear interpolation used here
            while(tgt_x3[k]>src_x3[src_pointer+1]) src_pointer++; // Find the left location of src
            Real y1 = src_x3[src_pointer];
            Real y2 = src_x3[src_pointer+1];
            Real v1 = src(n, sk+src_pointer, j, i);
            Real v2 = src(n, sk+src_pointer+1, j, i);
            Real yq = tgt_x3[k];
            tgt(n-sn, k-sk, j-sj, i-si) = ((y2-yq)*v1 + (yq-y1)*v2) / (y2-y1);
          }
        }
    }
  }
  delete[] src_x3, tgt_x3; // Release memory
}

void PackDataCubedSphere(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset,
         int ox1, int ox2, int ox3, LogicalLocation const& loc){
// Find the block ID
int blockID = FindBlockID(loc);
// Bypass Corner cases
if((ox2+ox3==0)||(ox2+ox3==2)||(ox2+ox3==-2)||(ox1!=0)){
  PackDataCubedSphereR0(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);
  return;
}
// Get the local indices
int lv2_lx2 = loc.lx2 >> (loc.level - 2);
int lv2_lx3 = loc.lx3 >> (loc.level - 2);
int local_lx2 = loc.lx2 - (lv2_lx2<<(loc.level - 2));
int local_lx3 = loc.lx3 - (lv2_lx3<<(loc.level - 2));
int bound_lim = (1<<(loc.level - 2)) - 1;

// Work on interpolation
AthenaArray<Real> interpolatedSrc;
interpolatedSrc.NewAthenaArray(en-sn+1, ek-sk+1, ej-sj+1, ei-si+1);

switch(blockID){
  case 1:
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 2:
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 3:
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 4:
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 5:
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 6:
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
}
PackDataCubedSphereR0(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);
return;
}