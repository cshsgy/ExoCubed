#include <sstream>
#include <ostream>
#include <iostream>
#include <bvals/bvals.hpp>
#include <cubed_sphere.hpp>
#include <coordinates/coordinates.hpp>
#include <field/field.hpp>
#include <cmath>
#include <algorithm>

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
  // Calculate the angular locations of the ghost zone interpolations
  // N_blk: the total number of points along the boundary
  // k: levels into the ghost zone (0, 1, 2 in interior = -1, -2, -3 in ghost)
  // loc_n: local index in the panel
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

void InteprolateX2CubedSphere(const AthenaArray<Real> &src, AthenaArray<Real> &tgt, LogicalLocation const& loc, int NRot, int DirInv, int TgtSide, int sn, int en, int si, int ei, int sj, int ej, int sk, int ek){
  // Interpolation along X2 (j) axis, used before sending data to X3 (k) axis
  // Get the local indices
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx2 = loc.lx2 - (lv2_lx2<<(loc.level - 2));
  int bound_lim = (1<<(loc.level - 2)) - 1;
  int N_blk = (ej - sj + 1) * (bound_lim+1); // N in X2 direction for each panel. This value is required to be an even number
  int n_start = local_lx2 * bound_lim - N_blk / 2;
  Real *src_x2 = new Real[N_blk]; // Need to calculate source indices along the whole panel boundary
  Real *src_coord = new Real[N_blk]; // The src coordinate along panel boundary
  Real *tgt_coord = new Real[ej-sj+1]; // The tgt coordinate along panel boundary
  Real *tgt_x2 = new Real[ej-sj+1];
  int SrcSide;

  for (int n=sn; n<=en; n++){
    for (int k=sk; k<=ek; k++) {
      int k_now;
      if(sk<ej-NGHOST){ // Calculate the location into the ghost zones
        k_now = k;
        SrcSide = -1;
      }else{
        k_now = ek-k;
        SrcSide = 1;
      }

#pragma omp simd
        for (int i=si; i<=ei; i++){
          for (int j=0; j<N_blk; j++){
            // Calculate coefficients for src, need to go from 0 to N_blk to cover all
            src_x2[j] = CalculateInterpLocationsCubedSphere(j-N_blk/2, N_blk, k_now, false);
            // Calculate the coordinate locations, going from -pi/4 to pi/4
            src_coord[j] = PI/2.0/N_blk*(j+0.5)-PI/4.0;
          }
          for (int j=sj; j<=ej; j++) {
            // Calculate coefficients for tgt
            tgt_x2[j-sj] = CalculateInterpLocationsCubedSphere(n_start+j-sj, N_blk, k_now, true);
            if(DirInv==1){
              tgt_coord[j-sj] = -src_coord[n_start+j-sj];
            }else{
              tgt_coord[j-sj] = src_coord[n_start+j-sj];
            }
          }
          int src_pointer = 0;
          for (int j=sj; j<=ej; j++) {
            // Interpolate to target array, linear interpolation used here
            while(tgt_x2[j-sj]>src_x2[src_pointer+1]) src_pointer++; // Find the left location of src

            Real y1 = src_x2[src_pointer];
            Real y2 = src_x2[src_pointer+1];
            Real yq = tgt_x2[j-sj];
            if (n==IVY || n==IVZ){
              // Projection needed, find the tgt locations first
              Real v1y = src(IVY, k, src_pointer-n_start+sj, i);
              Real v1z = src(IVZ, k, src_pointer-n_start+sj, i);
              Real v2y = src(IVY, k, src_pointer+1-n_start+sj, i);
              Real v2z = src(IVZ, k, src_pointer+1-n_start+sj, i);
              Real tgt_cy = tgt_coord[j-sj];
              Real tgt_cz = TgtSide * (PI/2.0/(N_blk)*(NGHOST-k_now-0.5)+PI/4.0);
              Real src_cy1 = src_coord[src_pointer];
              Real src_cy2 = src_coord[src_pointer+1];
              Real src_cz = SrcSide * (-PI/2.0/(N_blk)*(NGHOST-k_now-0.5)+PI/4.0);
              Real vy = ((y2-yq)*v1y + (yq-y1)*v2y) / (y2-y1);
              Real vz = ((y2-yq)*v1z + (yq-y1)*v2z) / (y2-y1);
              Real src_cy = ((y2-yq)*src_cy1 + (yq-y1)*src_cy2) / (y2-y1);
              Real s_sc = sqrt(1+src_cy*src_cy+src_cz*src_cz)/sqrt(1+src_cy*src_cy)/sqrt(1+src_cz*src_cz);
              Real c_sc = -src_cy*src_cz/sqrt(1+src_cy*src_cy)/sqrt(1+src_cz*src_cz);
              Real s_tg = sqrt(1+tgt_cy*tgt_cy+tgt_cz*tgt_cz)/sqrt(1+tgt_cy*tgt_cy)/sqrt(1+tgt_cz*tgt_cz);
              Real c_tg = -tgt_cy*tgt_cz/sqrt(1+tgt_cy*tgt_cy)/sqrt(1+tgt_cz*tgt_cz);
              Real o11 = 1.0;
              Real o12 = c_sc-c_tg/s_tg*s_sc;
              Real o21 = 0.0;
              Real o22 = s_sc/s_tg;
              for (int rt=0; rt<NRot; rt++){
                  Real o11t = -o21;
                  Real o12t = -o22;
                  Real o21t = o11;
                  Real o22t = o12;
                  o11 = o11t;
                  o12 = o12t;
                  o21 = o21t;
                  o22 = o22t;
              }
              if (n==IVY){
                tgt(n-sn, k-sk, j-sj, i-si) = o11*vy+o12*vz;
              }else{ // n==IVZ
                tgt(n-sn, k-sk, j-sj, i-si) = o21*vy+o22*vz;
              }
            }else{
              Real v1 = src(n, k, src_pointer-n_start+sj, i);
              Real v2 = src(n, k, src_pointer+1-n_start+sj, i);
              tgt(n-sn, k-sk, j-sj, i-si) = ((y2-yq)*v1 + (yq-y1)*v2) / (y2-y1);
            }
          }
        }
    }
  }
  delete[] src_x2, tgt_x2, src_coord, tgt_coord; // Release memory
}

void InteprolateX3CubedSphere(const AthenaArray<Real> &src, AthenaArray<Real> &tgt, LogicalLocation const& loc, int NRot, int DirInv, int TgtSide, int sn, int en, int si, int ei, int sj, int ej, int sk, int ek){
  // Interpolation along X3 (k) axis, used before sending data to ghost zone in X2 (j) direction
  // Get the local indices
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx3 = loc.lx3 - (lv2_lx3<<(loc.level - 2));
  int bound_lim = (1<<(loc.level - 2)) - 1;
  int N_blk = (ek - sk + 1) * (bound_lim+1); // N in X2 direction for each panel. This value is required to be an even number
  int n_start = local_lx3 * bound_lim - N_blk / 2;
  Real *src_x3 = new Real[N_blk]; // Need to calculate source indices along the whole panel boundary
  Real *src_coord = new Real[N_blk]; // The src coordinate along panel boundary
  Real *tgt_coord = new Real[ek-sk+1]; // The tgt coordinate along panel boundary
  Real *tgt_x3 = new Real[ek-sk+1];
  int SrcSide;

  for (int n=sn; n<=en; n++){
    for (int j=sj; j<=ej; j++) {
      int k_now;
      if(sj<ek-NGHOST) // Calculate the location into the ghost zones
        k_now = j;
      else
        k_now = ej-j;
#pragma omp simd
        for (int i=si; i<=ei; i++){
          for (int k=0; k<N_blk; k++){
            // Calculate coefficients for src, need to go from 0 to N_blk to cover all
            src_x3[k] = CalculateInterpLocationsCubedSphere(k-N_blk/2, N_blk, k_now, false);
            // Calculate the coordinate locations, going from -pi/4 to pi/4
            src_coord[k] = PI/2.0/N_blk*(k+0.5)-PI/4.0;
          }
          for (int k=sk; k<=ek; k++) {
            // Calculate coefficients for tgt
            tgt_x3[k-sk] = CalculateInterpLocationsCubedSphere(n_start+k-sk, N_blk, k_now, true);
            if(DirInv==1){
              tgt_coord[k-sk] = -src_coord[n_start+k-sk];
            }else{
              tgt_coord[k-sk] = src_coord[n_start+k-sk];
            }
          }
          int src_pointer = 0;
          for (int k=sk; k<=ek; k++) {
            // Interpolate to target array, linear interpolation used here
            while(tgt_x3[k-sk]>src_x3[src_pointer+1]) src_pointer++; // Find the left location of src

            Real y1 = src_x3[src_pointer];
            Real y2 = src_x3[src_pointer+1];
            Real yq = tgt_x3[k-sk];
            if (n==IVY || n==IVZ){
              // Projection needed, find the tgt locations first
              Real v1y = src(IVY, src_pointer-n_start+sk, j, i);
              Real v1z = src(IVZ, src_pointer-n_start+sk, j, i);
              Real v2y = src(IVY, src_pointer+1-n_start+sk, j, i);
              Real v2z = src(IVZ, src_pointer+1-n_start+sk, j, i);
              Real tgt_cz = tgt_coord[k-sk];
              Real tgt_cy = TgtSide * (PI/2.0/(N_blk)*(NGHOST-k_now-0.5)+PI/4.0);
              Real src_cz1 = src_coord[src_pointer];
              Real src_cz2 = src_coord[src_pointer+1];
              Real src_cy = SrcSide * (-PI/2.0/(N_blk)*(NGHOST-k_now-0.5)+PI/4.0);
              Real vy = ((y2-yq)*v1y + (yq-y1)*v2y) / (y2-y1);
              Real vz = ((y2-yq)*v1z + (yq-y1)*v2z) / (y2-y1);
              Real src_cz = ((y2-yq)*src_cz1 + (yq-y1)*src_cz2) / (y2-y1);
              Real s_sc = sqrt(1+src_cy*src_cy+src_cz*src_cz)/sqrt(1+src_cy*src_cy)/sqrt(1+src_cz*src_cz);
              Real c_sc = -src_cy*src_cz/sqrt(1+src_cy*src_cy)/sqrt(1+src_cz*src_cz);
              Real s_tg = sqrt(1+tgt_cy*tgt_cy+tgt_cz*tgt_cz)/sqrt(1+tgt_cy*tgt_cy)/sqrt(1+tgt_cz*tgt_cz);
              Real c_tg = -tgt_cy*tgt_cz/sqrt(1+tgt_cy*tgt_cy)/sqrt(1+tgt_cz*tgt_cz);
              Real o11 = 1.0;
              Real o12 = c_sc-c_tg/s_tg*s_sc;
              Real o21 = 0.0;
              Real o22 = s_sc/s_tg;
              for (int rt=0; rt<NRot; rt++){
                  Real o11t = -o21;
                  Real o12t = -o22;
                  Real o21t = o11;
                  Real o22t = o12;
                  o11 = o11t;
                  o12 = o12t;
                  o21 = o21t;
                  o22 = o22t;
              }
              if (n==IVY){
                tgt(n-sn, k-sk, j-sj, i-si) = o11*vy+o12*vz;
              }else{ // n==IVZ
                tgt(n-sn, k-sk, j-sj, i-si) = o21*vy+o22*vz;
              }
            }else{
              if((src_pointer-n_start+sk>ek)||(src_pointer-n_start+sk<sk)){
                std::cout<< "N_blk=" << N_blk << "|n_start=" << n_start << "|src_pt=" << src_pointer << "|start_k=" << sk <<std::endl;
              }
              Real v1 = src(n, src_pointer-n_start+sk, j, i);
              Real v2 = src(n, src_pointer+1-n_start+sk, j, i);
              tgt(n-sn, k-sk, j-sj, i-si) = ((y2-yq)*v1 + (yq-y1)*v2) / (y2-y1);
            }
          }
        }
    }
  }
  delete[] src_x3, tgt_x3, src_coord, tgt_coord; // Release memory
}

void PackDataCubedSphere(const AthenaArray<Real> &src, Real *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset,
         int ox1, int ox2, int ox3, LogicalLocation const& loc){
// Find the block ID
int blockID = FindBlockID(loc);

// Table of #rot needed
const int rot[6][4] = { // To access: rot[source_id][target_dir]
    {2,0,1,3},
    {0,0,0,0},
    {3,1,0,0},
    {1,3,0,0},
    {0,2,3,1},
    {2,2,0,0}
} ;

// Table of wether direction inversion is needed
const int dinv[6][4] = { // To access: dinv[source_id][target_dir]
    {1,0,0,1},
    {0,0,0,0},
    {0,1,0,0},
    {1,0,0,0},
    {0,1,1,0},
    {1,1,0,0}
} ;

// Table of which side (+pi/4 or -pi/4) is in touch with this panel
const int tgside[6][4] = { // To access: dinv[source_id][target_dir]
    {-1,-1,-1,-1},
    {1,-1,1,-1},
    {-1,-1,1,-1},
    {1,1,1,-1},
    {1,1,1,1},
    {-1,1,1,-1}
} ;

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
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][3], dinv[blockID-1][3], tgside[blockID-1][3], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][1], dinv[blockID-1][1], tgside[blockID-1][1], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][0], dinv[blockID-1][0], tgside[blockID-1][0], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][2], dinv[blockID-1][2], tgside[blockID-1][2], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 2:
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][3], dinv[blockID-1][3], tgside[blockID-1][3], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][1], dinv[blockID-1][1], tgside[blockID-1][1], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][0], dinv[blockID-1][0], tgside[blockID-1][0], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][2], dinv[blockID-1][2], tgside[blockID-1][2], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 3:
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][1], dinv[blockID-1][1], tgside[blockID-1][1], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);
      return; 
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][0], dinv[blockID-1][0], tgside[blockID-1][0], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][3], dinv[blockID-1][3], tgside[blockID-1][3], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][2], dinv[blockID-1][2], tgside[blockID-1][2], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 4:
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][1], dinv[blockID-1][1], tgside[blockID-1][1], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][0], dinv[blockID-1][0], tgside[blockID-1][0], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][3], dinv[blockID-1][3], tgside[blockID-1][3], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][2], dinv[blockID-1][2], tgside[blockID-1][2], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 5:
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][3], dinv[blockID-1][3], tgside[blockID-1][3], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR3(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][1], dinv[blockID-1][1], tgside[blockID-1][1], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][2], dinv[blockID-1][2], tgside[blockID-1][2], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR1(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][0], dinv[blockID-1][0], tgside[blockID-1][0], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    break;
  case 6:
    if(ox2==1 && local_lx2==bound_lim){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][1], dinv[blockID-1][1], tgside[blockID-1][1], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox2==-1 && local_lx2==0){
      InteprolateX3CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][0], dinv[blockID-1][0], tgside[blockID-1][0], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR2(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
    if(ox3==1 && local_lx3==bound_lim){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][3], dinv[blockID-1][3], tgside[blockID-1][3], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset); 
      return;
    }
    if(ox3==-1 && local_lx3==0){
      InteprolateX2CubedSphere(src, interpolatedSrc, loc, rot[blockID-1][2], dinv[blockID-1][2], tgside[blockID-1][2], sn, en, si, ei, sj, ej, sk, ek);
      PackDataCubedSphereR0(interpolatedSrc, buf, 0, en-sn, 0, ei-si, 0, ej-sj, 0, ek-sk, offset);  
      return;
    }
}
PackDataCubedSphereR0(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);
return;
}











// #include <sstream>
// #include <ostream>
// #include <iostream>
// #include <athena.hpp>
// #include <athena_arrays.hpp>
// #include <coordinates/coordinates.hpp>
// #include <field/field.hpp>
// #include <cubed_sphere.hpp>
// #include <cmath>
// #include <algorithm>

#define DBL_EPSILON 1.0e-10

void GetLatLon(Real *lat, Real *lon, Coordinates *pcoord, int k, int j, int i){
    // Obtain Lat and Lon (radians) from x2 and x3
    // k is not used for now
    // Find the block number
    LogicalLocation loc = pcoord->pmy_block->loc;
    int lv2_lx2 = loc.lx2 >> (loc.level - 2);
    int lv2_lx3 = loc.lx3 >> (loc.level - 2);
    int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
    // Calculate the needed parameters
    Real dX = pcoord->x2v(j);
    Real dY = pcoord->x3v(k);
    RLLFromXYP(dY, -dX, blockID-1, *lon, *lat);
}

void GetUV(Real *U, Real *V, Coordinates *pcoord, Real V2, Real V3, int k, int j, int i){
    // Obtain U and V (Lat-Lon) from V2 and V3 (Gnomonic Equiangle)
    // U is V_lam, V is V_phi.

    // Find the block number
    LogicalLocation loc = pcoord->pmy_block->loc;
    int lv2_lx2 = loc.lx2 >> (loc.level - 2);
    int lv2_lx3 = loc.lx3 >> (loc.level - 2);
    int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
    // Calculate the needed parameters
    Real dX = pcoord->x2v(j);
    Real dY = pcoord->x3v(k);
    Real tmp_U, tmp_V;
    // Calls Paul Ullrich's code. Note how the values are transformed here.
    VecTransRLLFromABP(dY, -dX, blockID-1, V3, -V2, tmp_V, tmp_U);
    *U = tmp_U;
    *V = tmp_V;
}

void GetVyVz(Real *V2, Real *V3, Coordinates *pcoord, Real U, Real V, int k, int j, int i){
    // Convert U and V (Lat-Lon) to V2 and V3 (Gnomonic Equiangle)
    // U is V_lam, V is V_phi.

    // Find the block number
    LogicalLocation loc = pcoord->pmy_block->loc;
    int lv2_lx2 = loc.lx2 >> (loc.level - 2);
    int lv2_lx3 = loc.lx3 >> (loc.level - 2);
    int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
    // Calculate the needed parameters
    Real dX = pcoord->x2v(j);
    Real dY = pcoord->x3v(k);
    Real tmp_V2, tmp_V3;
    // Calls Paul Ullrich's code. Note how the values are transformed here.
    VecTransABPFromRLL(dY, -dX, blockID-1, V, U, tmp_V2, tmp_V3);
    *V2 = -tmp_V3;
    *V3 = tmp_V2;
}


////////////////////////////////////////////////////////////////////////////////
// Following are codes adapted from P.U. tempest model 
// https://github.com/paullric/tempestmodel/blob/master/src/atm/CubedSphereTrans.cpp
////////////////////////////////////////////////////////////////////////////////

void VecTransABPFromRLL(
	Real dX,
	Real dY,
	int nP,
	Real dUlon,
	Real dUlat,
	Real & dUalpha,
	Real & dUbeta
) {
	Real dDelta2 = 1.0 + dX * dX + dY * dY;
	Real dRadius;

	Real lat;

	if (((nP==0)||(nP==4)) && (fabs(dX) < 1.0e-13) && (fabs(dY) < 1.0e-13)) {
		if (nP == 0) {
			dUalpha = dUlon;
		} else {
			dUalpha = - dUlon;
		}
		dUbeta = dUlat;
		return;
	}

	switch (nP) {
		// Equatorial panels
		case 1:
		case 2:
		case 3:
		case 5:
			// Convert spherical coords to geometric basis
			lat = atan(dY / sqrt(1.0 + dX * dX));
			dUlon = dUlon / cos(lat);

			// Calculate new vector components
			dUalpha = dUlon;
			dUbeta =
				dX * dY / (1.0 + dY * dY) * dUlon
				+ dDelta2 / ((1.0 + dY * dY) * sqrt(1.0 + dX * dX)) * dUlat;
			break;

		// North polar panel
		case 0:
			// Convert spherical coords to geometric basis
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon / cos(lat);

			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUalpha =
				- dY / (1.0 + dX * dX) * dUlon
				- dDelta2 * dX / ((1.0 + dX * dX) * dRadius) * dUlat;

			dUbeta =
				dX / (1.0 + dY * dY) * dUlon
				- dDelta2 * dY / ((1.0 + dY * dY) * dRadius) * dUlat;
			break;

		// South polar panel
		case 4:
			// Convert spherical coords to geometric basis
			lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon / cos(lat);

			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUalpha =
				dY / (1.0 + dX * dX) * dUlon
				+ dDelta2 * dX / ((1.0 + dX * dX) * dRadius) * dUlat;

			dUbeta =
				- dX / (1.0 + dY * dY) * dUlon
				+ dDelta2 * dY / ((1.0 + dY * dY) * dRadius) * dUlat;
			break;

	}
}

////////////////////////////////////////////////////////////////////////////////

void VecTransRLLFromABP(
	Real dX,
	Real dY,
	int nP,
	Real dUalpha,
	Real dUbeta,
	Real & dUlon,
	Real & dUlat
) {
	Real dDelta2 = 1.0 + dX * dX + dY * dY;
	Real dRadius;

	Real lat;

	switch (nP) {
		// Equatorial panels
		case 1:
		case 2:
		case 3:
		case 5:
			// Calculate new vector components
			dUlon = dUalpha;
			dUlat = 
				- dX * dY * sqrt(1.0 + dX * dX) / dDelta2 * dUalpha
				+ (1.0 + dY * dY) * sqrt(1.0 + dX * dX) / dDelta2 * dUbeta;

			// Convert spherical coords to unit basis
			lat = atan(dY / sqrt(1.0 + dX * dX));
			dUlon = dUlon * cos(lat);
			break;

		// North polar panel
		case 0:
			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUlon = 
				- dY * (1.0 + dX * dX) / (dRadius * dRadius) * dUalpha
				+ dX * (1.0 + dY * dY) / (dRadius * dRadius) * dUbeta;

			dUlat =
				- dX * (1.0 + dX * dX) / (dDelta2 * dRadius) * dUalpha
				- dY * (1.0 + dY * dY) / (dDelta2 * dRadius) * dUbeta;

			// Convert spherical coords to unit basis
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon * cos(lat);

			break;

		// South polar panel
		case 4:
			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUlon =
				dY * (1.0 + dX * dX) / (dRadius * dRadius) * dUalpha
				- dX * (1.0 + dY * dY) / (dRadius * dRadius) * dUbeta;

			dUlat =
				dX * (1.0 + dX * dX) / (dDelta2 * dRadius) * dUalpha
				+ dY * (1.0 + dY * dY) / (dDelta2 * dRadius) * dUbeta;
		
			// Convert spherical coords to unit basis
			lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon * cos(lat);

			break;
	}
}

void RLLFromXYP(
	Real dX,
	Real dY,
	int nP,
	Real &lon,
	Real &lat
) {
	switch (nP) {
		// Equatorial panel 2
		case 1:
			lon = atan(dX);
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// Equatorial panel 4
		case 3:
			lon = atan(dX) + 0.5 * M_PI;
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// Equatorial panel 6
		case 5:
			lon = atan(dX) + M_PI;
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// Equatorial panel 3
		case 2:
			lon = atan(dX) + 1.5 * M_PI;
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// North polar panel
		case 0:
			if (fabs(dX) > DBL_EPSILON) {
				lon = atan2(dX, -dY);
			} else if (dY <= 0.0) {
				lon = 0.0;
			} else {
				lon = M_PI;
			}
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			break;

		// South polar panel
		case 4:
			if (fabs(dX) > DBL_EPSILON) {
				lon = atan2(dX, dY);
			} else if (dY > 0.0) {
				lon = 0.0;
			} else {
				lon = M_PI;
			}
			lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
			break;
	}

	// Map to the interval [0, 2 pi]
	if (lon < 0.0) {
		lon += 2.0 * M_PI;
	}
}

////////////////////////////////////////////////////////////////////////////////
// May be used later... Note that mo modification of this has been done yet.
void XYPFromRLL(
	Real lon,
	Real lat,
	Real &dX,
	Real &dY,
	int &nP
) {
	// Default panel to unattainable value
	nP = 6;

	// Translate from RLL coordinates to XYZ space
	Real xx, yy, zz, pm;

	xx = cos(lon) * cos(lat);
	yy = sin(lon) * cos(lat);
	zz = sin(lat);

	pm = std::max(fabs(xx), std::max(fabs(yy), fabs(zz)));

	// Check maxmality of the x coordinate
	if (pm == fabs(xx)) {
		if (xx > 0) {
			nP = 0;
		} else {
			nP = 2;
		}
	}

	// Check maximality of the y coordinate
	if (pm == fabs(yy)) {
		if (yy > 0) {
			nP = 1;
		} else {
			nP = 3;
		}
	}

	// Check maximality of the z coordinate
	if (pm == fabs(zz)) {
		if (zz > 0) {
			nP = 4;
		} else {
			nP = 5;
		}
	}

	// Panel assignments
	Real sx, sy, sz;
	if (nP == 0) {
		sx = yy;
		sy = zz;
		sz = xx;

	} else if (nP == 1) {
		sx = -xx;
		sy = zz;
		sz = yy;

	} else if (nP == 2) {
		sx = -yy;
		sy = zz;
		sz = -xx;

	} else if (nP == 3) {
		sx = xx;
		sy = zz;
		sz = -yy;

	} else if (nP == 4) {
		sx = yy;
		sy = -xx;
		sz = zz;

	} else if (nP == 5) {
		sx = yy;
		sy = xx;
		sz = -zz;

	}

	// Convert to gnomonic coordinates
	dX = sx / sz;
	dY = sy / sz;
}
