//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_fc.cpp
//  \brief functions that perform flux correction for FACE_CENTERED variables

// C headers

// C++ headers
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "bvals_fc.hpp"

// this is not added in flux_correction_cc.cpp:
// #include "../bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadEMFBoundaryBufferSameLevel(Real *buf,
//                                                   const NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the same level

int FaceCenteredBoundaryVariable::LoadEMFBoundaryBufferSameLevel(
    Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;

  // KGF: shearing box:
  // Real qomL = qshear_*Omega_0_*x1size_;
  // AthenaArray<Real> &bx1=pmb->pfield->b.x1f;

  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        // pack e2
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i);
        }
        // pack e3

        // KGF: shearing box
        // shift azmuthal velocity if shearing boundary blocks
        // if (nb.shear && nb.fid==INNER_X1) {
        //   for (int k=pmb->ks; k<=pmb->ke; k++) {
        //     for (int j=pmb->js; j<=pmb->je+1; j++)
        //       buf[p++]=e3(k,j,i)-0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
        //   }
        // } else if (nb.shear && nb.fid==OUTER_X1) {
        //   for (int k=pmb->ks; k<=pmb->ke; k++) {
        //     for (int j=pmb->js; j<=pmb->je+1; j++)
        //       buf[p++]=e3(k,j,i)+0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
        //   }
        // } else {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              buf[p++]=e3(k,j,i);
          }
          // } // KGF: shearing box
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        // pack e1
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            buf[p++]=e1(k,j,i);
        }
        // pack e3
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++]=e3(k,j,i);
        }
        // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if (nb.fid==INNER_X3) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        // pack e1
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            buf[p++]=e1(k,j,i);
        }
        // pack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++]=e2(k,j,i);
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        // pack e2
        // KGF: shearing box
        // shift azimuthal velocity for x-z shearing
        // if (SHEARING_BOX) {
        //   if (ShBoxCoord_==2 && (pmb->loc.lx1==0) && (nb.ox1==-1)) {
        //     for (int j=pmb->js; j<=pmb->je; j++)
        //       buf[p++]=e2(k,j,i)+qomL*bx1(k,j,i);
        //   } else if (ShBoxCoord_==2 && (pmb->loc.lx1==(pmb->pmy_mesh->nrbx1-1))
        //              && nb.ox1==1) {
        //     for (int j=pmb->js; j<=pmb->je; j++)
        //       buf[p++]=e2(k,j,i)-qomL*bx1(k,j,i);
        //   } else {
        //     for (int j=pmb->js; j<=pmb->je; j++)
        //       buf[p++]=e2(k,j,i);
        //   }
        // } else {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++]=e2(k,j,i);
          // } // KGF: shearing box
        // pack e3
        for (int j=pmb->js; j<=pmb->je+1; j++)
          buf[p++]=e3(k,j,i);
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        // pack e1
        for (int i=pmb->is; i<=pmb->ie; i++)
          buf[p++]=e1(k,j,i);
        // pack e3
        for (int i=pmb->is; i<=pmb->ie+1; i++)
          buf[p++]=e3(k,j,i);
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      // pack e2 and e3
      buf[p++]=e2(k,j,i);
      buf[p++]=e3(k,j,i);
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    // x1x2 edge (both 2D and 3D)
    if (nb.eid>=0 && nb.eid<4) {
      int i, j;
      if ((nb.eid & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((nb.eid & 2)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      // KGF: shearing box
      // shift azmuthal velocity if shearing boundary blocks
      // if (nb.shear && nb.ox1==-1) {
      //   for (int k=pmb->ks; k<=pmb->ke; k++)
      //     buf[p++]=e3(k,j,i)-0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
      // } else if (nb.shear && nb.ox1==1) {
      //   for (int k=pmb->ks; k<=pmb->ke; k++)
      //     buf[p++]=e3(k,j,i)+0.5*qomL*(bx1(k,j,i)+bx1(k,j-1,i));
      // } else {
        // pack e3
        for (int k=pmb->ks; k<=pmb->ke; k++)
          buf[p++]=e3(k,j,i);
        //      }         // KGF: shearing box
      // x1x3 edge
    } else if (nb.eid>=4 && nb.eid<8) {
      int i, k;
      if ((nb.eid & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((nb.eid & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      // pack e2
      // KGF: shearing box
      // shift azimuthal velocity for x-z shearing
      // if (SHEARING_BOX) {
      //   if (ShBoxCoord_==2 && (pmb->loc.lx1==0) && (nb.ox1==-1))   {
      //     for (int j=pmb->js; j<=pmb->je; j++)
      //       buf[p++]=e2(k,j,i)+qomL*bx1(k,j,i);
      //   } else if (ShBoxCoord_==2 && (pmb->loc.lx1==(pmb->pmy_mesh->nrbx1-1))
      //              && nb.ox1==1) {
      //     for (int j=pmb->js; j<=pmb->je; j++)
      //       buf[p++]=e2(k,j,i)-qomL*bx1(k,j,i);
      //   } else {
      //     for (int j=pmb->js; j<=pmb->je; j++)
      //       buf[p++]=e2(k,j,i);
      //   }
      // } else {
        for (int j=pmb->js; j<=pmb->je; j++)
          buf[p++]=e2(k,j,i);
        // } // KGF: shearing box
      // x2x3 edge
    } else if (nb.eid>=8 && nb.eid<12) {
      int j, k;
      if ((nb.eid & 1)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      if ((nb.eid & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      // pack e1
      for (int i=pmb->is; i<=pmb->ie; i++)
        buf[p++]=e1(k,j,i);
    }
  }
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadEMFBoundaryBufferToCoarser(Real *buf,
//                                                         const NeighborBlock& nb)
//  \brief Set EMF correction buffers for sending to a block on the coarser level

int FaceCenteredBoundaryVariable::LoadEMFBoundaryBufferToCoarser(
    Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  Coordinates *pco=pmb->pcoord;
  // use the surface area aray as the edge length array
  AthenaArray<Real> &le1=pbval_->sarea_[0];
  AthenaArray<Real> &le2=pbval_->sarea_[1];
  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        // restrict and pack e2
        for (int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          for (int j=pmb->js; j<=pmb->je; j+=2) {
            Real el1=pco->GetEdge2Length(k,j,i);
            Real el2=pco->GetEdge2Length(k,j+1,i);
            buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
          }
        }
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          for (int j=pmb->js; j<=pmb->je+1; j+=2) {
            bool pole = pco->IsPole(j);
            Real el1, el2;
            if (!pole) {
              el1 = pco->GetEdge3Length(k,j,i);
              el2 = pco->GetEdge3Length(k+1,j,i);
            } else {
              el1 = pco->dx3f(k);
              el2 = pco->dx3f(k+1);
            }
            buf[p++]=(e3(k,j,i)*el1+e3(k+1,j,i)*el2)/(el1+el2);
          }
        }
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        bool pole = pco->IsPole(j);
        // restrict and pack e1
        for (int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          if (!pole || !GENERAL_RELATIVITY) {
            pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          } else {
            for (int i = pmb->is; i <= pmb->ie+1; i+=2) {
              le1(i) = pco->dx1f(i);
              le1(i+1) = pco->dx1f(i+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        }
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          if (!pole) {
            pco->Edge3Length(k,   j, pmb->is, pmb->ie+1, le1);
            pco->Edge3Length(k+1, j, pmb->is, pmb->ie+1, le2);
          } else {
            for (int i = pmb->is; i <= pmb->ie+1; i+=2) {
              le1(i) = pco->dx3f(k);
              le2(i) = pco->dx3f(k+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++]=(e3(k,j,i)*le1(i)+e3(k+1,j,i)*le2(i))/(le1(i)+le2(i));
        }
        // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if (nb.fid==INNER_X3) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        // restrict and pack e1
        for (int j=pmb->js; j<=pmb->je+1; j+=2) {
          bool pole = pco->IsPole(j);
          if (!pole || !GENERAL_RELATIVITY) {
            pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          } else {
            for (int i = pmb->is; i <= pmb->ie; i+=2) {
              le1(i) = pco->dx1f(i);
              le1(i+1) = pco->dx1f(i+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        }
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          pco->Edge2Length(k,   j, pmb->is, pmb->ie+1, le1);
          pco->Edge2Length(k, j+1, pmb->is, pmb->ie+1, le2);
          for (int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++]=(e2(k,j,i)*le1(i)+e2(k,j+1,i)*le2(i))/(le1(i)+le2(i));
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1=pco->GetEdge2Length(k,j,i);
          Real el2=pco->GetEdge2Length(k,j+1,i);
          buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
        }
        // pack e3
        for (int j=pmb->js; j<=pmb->je+1; j+=2)
          buf[p++]=e3(k,j,i);
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        bool pole = pco->IsPole(j);
        // restrict and pack e1
        if (!pole || !GENERAL_RELATIVITY) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        } else {
          for (int i = pmb->is; i <= pmb->ie; i+=2) {
            le1(i) = pco->dx1f(i);
            le1(i+1) = pco->dx1f(i+1);
          }
        }
        for (int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
        // pack e3
        for (int i=pmb->is; i<=pmb->ie+1; i+=2)
          buf[p++]=e3(k,j,i);
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      // pack e2 and e3
      buf[p++]=e2(k,j,i);
      buf[p++]=e3(k,j,i);
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if (nb.eid>=0 && nb.eid<4) {
        int i, j;
        if ((nb.eid & 1)==0) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        if ((nb.eid & 2)==0) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        bool pole = pco->IsPole(j);
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          Real el1, el2;
          if (!pole) {
            el1 = pco->GetEdge3Length(k,j,i);
            el2 = pco->GetEdge3Length(k+1,j,i);
          } else {
            el1 = pco->dx3f(k);
            el2 = pco->dx3f(k+1);
          }
          buf[p++]=(e3(k,j,i)*el1+e3(k+1,j,i)*el2)/(el1+el2);
        }
        // x1x3 edge
      } else if (nb.eid>=4 && nb.eid<8) {
        int i, k;
        if ((nb.eid & 1)==0) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        if ((nb.eid & 2)==0) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1=pco->GetEdge2Length(k,j,i);
          Real el2=pco->GetEdge2Length(k,j+1,i);
          buf[p++]=(e2(k,j,i)*el1+e2(k,j+1,i)*el2)/(el1+el2);
        }
        // x2x3 edge
      } else if (nb.eid>=8 && nb.eid<12) {
        int j, k;
        if ((nb.eid & 1)==0) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        bool pole = pco->IsPole(j);
        if ((nb.eid & 2)==0) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        // restrict and pack e1
        if (!pole || !GENERAL_RELATIVITY) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        } else {
          for (int i = pmb->is; i <= pmb->ie; i+=2) {
            le1(i) = pco->dx1f(i);
            le1(i+1) = pco->dx1f(i+1);
          }
        }
        for (int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++]=(e1(k,j,i)*le1(i)+e1(k,j,i+1)*le1(i+1))/(le1(i)+le1(i+1));
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      // x1x2 edge
      int i, j;
      if ((nb.eid & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((nb.eid & 2)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      // pack e3
      buf[p++]=e3(pmb->ks,j,i);
    }
  }
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadEMFBoundaryPolarBuffer(Real *buf,
//                                                           const PolarNeighborBlock &nb)
//  \brief Load EMF values along polar axis into send buffers

int FaceCenteredBoundaryVariable::LoadEMFBoundaryPolarBuffer(
    Real *buf, const PolarNeighborBlock &nb) {
  MeshBlock *pmb = pmy_block_;
  int count = 0;
  int j = nb.north ? pmb->js : pmb->je+1;
  for (int i = pmb->is; i <= pmb->ie; ++i) {
    Real val = 0.0;
    for (int k = pmb->ks; k <= pmb->ke; ++k) {  // avoid double counting right ends
      val += pmb->pfield->e.x1e(k, j, i);
    }
    buf[count++] = val / (pmb->ke - pmb->ks + 1);
  }
  return count;
}

// KGF: helper function for below SendFluxCorrection()
void FaceCenteredBoundaryVariable::CopyPolarBufferSameProcess(
    const PolarNeighborBlock& nb, int ssize, int polar_block_index, bool is_north) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block=pmy_mesh_->FindMeshBlock(nb.gid);
  // 2) which element in vector of BoundaryVariable *?
  // KGF: do we need to typecast from generic "BoundaryVariable *" ?
  FaceCenteredBoundaryVariable *ptarget_pfbval =
      static_cast<FaceCenteredBoundaryVariable *>(
          ptarget_block->pbval->bvars[bvar_index]);
  Real *target_buf, *send_buf;
  enum BoundaryStatus target_flag;
  if (is_north) {
    target_buf= ptarget_pfbval->emf_north_recv_[pmy_block_->loc.lx3];
    send_buf = emf_north_send_[polar_block_index];
    target_flag = ptarget_pfbval->emf_north_flag_[pmy_block_->loc.lx3];
  } else {
    target_buf= ptarget_pfbval->emf_south_recv_[pmy_block_->loc.lx3];
    send_buf = emf_south_send_[polar_block_index];
    target_flag = ptarget_pfbval->emf_south_flag_[pmy_block_->loc.lx3];
  }
  std::memcpy(target_buf, send_buf, ssize*sizeof(Real));
  target_flag = BNDRY_ARRIVED;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendEMFCorrection()
//  \brief Restrict, pack and send the surface EMF to the coarse neighbor(s) if
//  needed
void FaceCenteredBoundaryVariable::SendFluxCorrection() {
  MeshBlock *pmb=pmy_block_;

  // Send non-polar EMF values
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if ((nb.type!=NEIGHBOR_FACE) && (nb.type!=NEIGHBOR_EDGE)) break;
    int p=0;
    if (nb.level==pmb->loc.level) {
      if ((nb.type==NEIGHBOR_FACE)
          || ((nb.type==NEIGHBOR_EDGE) && (pbval_->edge_flag_[nb.eid]==true))) {
        p=LoadEMFBoundaryBufferSameLevel(bd_var_flcor_.send[nb.bufid], nb);
      } else {
        continue;
      }
    } else if (nb.level==pmb->loc.level-1) {
      p=LoadEMFBoundaryBufferToCoarser(bd_var_flcor_.send[nb.bufid], nb);
    } else {
      continue;
    }
    if (nb.rank==Globals::my_rank) { // on the same MPI rank
      CopyFluxCorrectionBufferSameProcess(nb, p);
      // KGF: double check
      // std::memcpy(pbl->pbval->bd_emfcor_.recv[nb.targetid],
      //             bd_var_flcor_.send[nb.bufid], p*sizeof(Real));
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&(bd_var_flcor_.req_send[nb.bufid]));
#endif
  }

  // Send polar EMF values
  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pbval_->polar_neighbor_north[n];
    int count = LoadEMFBoundaryPolarBuffer(emf_north_send_[n], nb);
    if (nb.rank == Globals::my_rank) { // on the same MPI rank
      CopyPolarBufferSameProcess(nb, count, n, true);
      // KGF: check this
      // MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      // std::memcpy(pbl->pbval->emf_north_recv_[pmb->loc.lx3],
      //             emf_north_send_[n], count * sizeof(Real));
      // pbl->pbval->emf_north_flag_[pmb->loc.lx3] = BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_emf_north_send_[n]);
#endif
  }
  for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pbval_->polar_neighbor_south[n];
    int count = LoadEMFBoundaryPolarBuffer(emf_south_send_[n], nb);
    if (nb.rank == Globals::my_rank) { // on the same node
      CopyPolarBufferSameProcess(nb, count, n, false);
      // KGF: check this
      // MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      // std::memcpy(pbl->pbval->emf_south_recv_[pmb->loc.lx3],
      //             emf_south_send_[n], count * sizeof(Real));
      // pbl->pbval->emf_south_flag_[pmb->loc.lx3] = BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_emf_south_send_[n]);
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetEMFBoundarySameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the same level
//         Later they will be divided in the AverageEMFBoundary function

void FaceCenteredBoundaryVariable::SetEMFBoundarySameLevel(Real *buf,
                                                           const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        // KGF: shearing box
        // if (nb.shear && nb.fid==INNER_X1) {
        //   // store e2 for shearing periodic bcs
        //   for (int k=pmb->ks; k<=pmb->ke+1; k++) {
        //     for (int j=pmb->js; j<=pmb->je; j++)
        //       shboxvar_inner_emf_.x2e(k,j) += buf[p++];
        //   }
        //   // store e3 for shearing periodic bcs
        //   for (int k=pmb->ks; k<=pmb->ke; k++) {
        //     for (int j=pmb->js; j<=pmb->je+1; j++)
        //       shboxvar_inner_emf_.x3e(k,j) += buf[p++];
        //   }
        // } else if (nb.shear && nb.fid==OUTER_X1) {
        //   // store e2 for shearing periodic bcs
        //   for (int k=pmb->ks; k<=pmb->ke+1; k++) {
        //     for (int j=pmb->js; j<=pmb->je; j++)
        //       shboxvar_outer_emf_.x2e(k,j) += buf[p++];
        //   }
        //   // store e3 for shearing periodic bcs
        //   for (int k=pmb->ks; k<=pmb->ke; k++) {
        //     for (int j=pmb->js; j<=pmb->je+1; j++)
        //       shboxvar_outer_emf_.x3e(k,j) += buf[p++];
        //   }
        // } else {
          // unpack e2
          for (int k=pmb->ks; k<=pmb->ke+1; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)+=buf[p++];
          }
          // unpack e3
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              e3(k,j,i)+=buf[p++];
          }
          //} // KGF: shearing box
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        // unpack e1
        Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)+=sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            e3(k,j,i)+=sign*buf[p++];
        }
        // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k;
        if (nb.fid==INNER_X3) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        // unpack e1
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            e2(k,j,i)+=buf[p++];
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          e2(k+1,j,i)+=buf[p];
          e2(k,  j,i)+=buf[p++];
        }
        // unpack e3
        for (int j=pmb->js; j<=pmb->je+1; j++)
          e3(k,j,i)+=buf[p++];
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        // unpack e1
        for (int i=pmb->is; i<=pmb->ie; i++) {
          e1(k+1,j,i)+=buf[p];
          e1(k  ,j,i)+=buf[p++];
        }
        // unpack e3
        for (int i=pmb->is; i<=pmb->ie+1; i++)
          e3(k,j,i)+=buf[p++];
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      // unpack e2
      e2(k+1,j,i)+=buf[p];
      e2(k  ,j,i)+=buf[p++];
      // unpack e3
      e3(k,j+1,i)+=buf[p];
      e3(k  ,j,i)+=buf[p++];
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    // x1x2 edge (2D and 3D)
    if (nb.eid>=0 && nb.eid<4) {
      int i, j;
      if ((nb.eid & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((nb.eid & 2)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      // KGF: shearing box
      // if (nb.shear && nb.ox1==-1) {
      //   // store e3 for shearing periodic bcs
      //   for (int k=pmb->ks; k<=pmb->ke; k++)
      //     shboxvar_inner_emf_.x3e(k,j) += buf[p++];
      // } else if (nb.shear && nb.ox1==1) {
      //   // store e3 for shearing periodic bcs
      //   for (int k=pmb->ks; k<=pmb->ke; k++)
      //     shboxvar_outer_emf_.x3e(k,j) += buf[p++];
      // } else {
        // unpack e3
        Real sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          e3(k,j,i)+=sign*buf[p++];
        }
        //} // KGF: shearing box
      // x1x3 edge
    } else if (nb.eid>=4 && nb.eid<8) {
      int i, k;
      if ((nb.eid & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((nb.eid & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      // KGF: shearing box
      // if (nb.shear && nb.ox1==-1) {
      //   // store e2 for shearing periodic bcs
      //   for (int j=pmb->js; j<=pmb->je; j++)
      //     shboxvar_inner_emf_.x2e(k,j) += buf[p++];
      // } else if (nb.shear && nb.ox1==1) {
      //   // store e2 for shearing periodic bcs
      //   for (int j=pmb->js; j<=pmb->je; j++)
      //     shboxvar_outer_emf_.x2e(k,j) += buf[p++];
      // } else {
        // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i)+=buf[p++];
        // } // KGF: shearing box
      // x2x3 edge
    } else if (nb.eid>=8 && nb.eid<12) {
      int j, k;
      if ((nb.eid & 1)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      if ((nb.eid & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      // unpack e1
      Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)+=sign*buf[p++];
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetEMFBoundaryFromFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Add up the EMF received from a block on the finer level
//         Later they will be divided in the AverageEMFBoundary function

void FaceCenteredBoundaryVariable::SetEMFBoundaryFromFiner(Real *buf,
                                                           const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int p=0;
  if (nb.type==NEIGHBOR_FACE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        if (nb.fi1==0) {
          ju=pmb->js+pmb->block_size.nx2/2-1;
        } else {
          jl=pmb->js+pmb->block_size.nx2/2;
        }
        if (nb.fi2==0) {
          ku=pmb->ks+pmb->block_size.nx3/2-1;
        } else {
          kl=pmb->ks+pmb->block_size.nx3/2;
        }
        // unpack e2
        for (int k=kl; k<=ku+1; k++) {
          for (int j=jl; j<=ju; j++)
            e2(k,j,i)+=buf[p++];
        }
        // unpack e3
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju+1; j++)
            e3(k,j,i)+=buf[p++];
        }
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j, il=pmb->is, iu=pmb->ie, kl=pmb->ks, ku=pmb->ke;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        if (nb.fi1==0) {
          iu=pmb->is+pmb->block_size.nx1/2-1;
        } else {
          il=pmb->is+pmb->block_size.nx1/2;
        }
        if (nb.fi2==0) {
          ku=pmb->ks+pmb->block_size.nx3/2-1;
        } else {
          kl=pmb->ks+pmb->block_size.nx3/2;
        }
        // unpack e1
        Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku+1; k++) {
          for (int i=il; i<=iu; i++)
            e1(k,j,i)+=sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku; k++) {
          for (int i=il; i<=iu+1; i++)
            e3(k,j,i)+=sign*buf[p++];
        }
        // x3 direction
      } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
        int k, il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je;
        if (nb.fid==INNER_X3) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        if (nb.fi1==0) {
          iu=pmb->is+pmb->block_size.nx1/2-1;
        } else {
          il=pmb->is+pmb->block_size.nx1/2;
        }
        if (nb.fi2==0) {
          ju=pmb->js+pmb->block_size.nx2/2-1;
        } else {
          jl=pmb->js+pmb->block_size.nx2/2;
        }
        // unpack e1
        for (int j=jl; j<=ju+1; j++) {
          for (int i=il; i<=iu; i++)
            e1(k,j,i)+=buf[p++];
        }
        // unpack e2
        for (int j=jl; j<=ju; j++) {
          for (int i=il; i<=iu+1; i++)
            e2(k,j,i)+=buf[p++];
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k=pmb->ks;
      // x1 direction
      if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
        int i, jl=pmb->js, ju=pmb->je;
        if (nb.fid==INNER_X1) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        if (nb.fi1==0) {
          ju=pmb->js+pmb->block_size.nx2/2-1;
        } else {
          jl=pmb->js+pmb->block_size.nx2/2;
        }
        // unpack e2
        for (int j=jl; j<=ju; j++) {
          e2(k+1,j,i)+=buf[p];
          e2(k,  j,i)+=buf[p++];
        }
        // unpack e3
        for (int j=jl; j<=ju+1; j++)
          e3(k,j,i)+=buf[p++];
        // x2 direction
      } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
        int j, il=pmb->is, iu=pmb->ie;
        if (nb.fid==INNER_X2) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        if (nb.fi1==0) {
          iu=pmb->is+pmb->block_size.nx1/2-1;
        } else {
          il=pmb->is+pmb->block_size.nx1/2;
        }
        // unpack e1
        for (int i=il; i<=iu; i++) {
          e1(k+1,j,i)+=buf[p];
          e1(k  ,j,i)+=buf[p++];
        }
        // unpack e3
        for (int i=il; i<=iu+1; i++)
          e3(k,j,i)+=buf[p++];
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid==INNER_X1) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      // unpack e2
      e2(k+1,j,i)+=buf[p];
      e2(k  ,j,i)+=buf[p++];
      // unpack e3
      e3(k,j+1,i)+=buf[p];
      e3(k  ,j,i)+=buf[p++];
    }
  } else if (nb.type==NEIGHBOR_EDGE) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if (nb.eid>=0 && nb.eid<4) {
        int i, j, kl=pmb->ks, ku=pmb->ke;
        if ((nb.eid & 1)==0) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        if ((nb.eid & 2)==0) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        if (nb.fi1==0) {
          ku=pmb->ks+pmb->block_size.nx3/2-1;
        } else {
          kl=pmb->ks+pmb->block_size.nx3/2;
        }
        // unpack e3
        Real sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku; k++)
          e3(k,j,i)+=sign*buf[p++];
        // x1x3 edge
      } else if (nb.eid>=4 && nb.eid<8) {
        int i, k, jl=pmb->js, ju=pmb->je;
        if ((nb.eid & 1)==0) {
          i=pmb->is;
        } else {
          i=pmb->ie+1;
        }
        if ((nb.eid & 2)==0) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        if (nb.fi1==0) {
          ju=pmb->js+pmb->block_size.nx2/2-1;
        } else {
          jl=pmb->js+pmb->block_size.nx2/2;
        }
        // unpack e2
        for (int j=jl; j<=ju; j++)
          e2(k,j,i)+=buf[p++];
        // x2x3 edge
      } else if (nb.eid>=8 && nb.eid<12) {
        int j, k, il=pmb->is, iu=pmb->ie;
        if ((nb.eid & 1)==0) {
          j=pmb->js;
        } else {
          j=pmb->je+1;
        }
        if ((nb.eid & 2)==0) {
          k=pmb->ks;
        } else {
          k=pmb->ke+1;
        }
        if (nb.fi1==0) {
          iu=pmb->is+pmb->block_size.nx1/2-1;
        } else {
          il=pmb->is+pmb->block_size.nx1/2;
        }
        // unpack e1
        Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int i=il; i<=iu; i++)
          e1(k,j,i)+=sign*buf[p++];
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int i, j, k=pmb->ks;
      if ((nb.eid & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((nb.eid & 2)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      // unpack e3
      e3(k,j,i)+=buf[p++];
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetEMFBoundaryPolar(Real **buf_list,
//                                                             int num_bufs, bool north)
//  \brief Overwrite EMF values along polar axis with azimuthal averages

void FaceCenteredBoundaryVariable::SetEMFBoundaryPolar(Real **buf_list, int num_bufs,
                                                       bool north) {
  MeshBlock *pmb = pmy_block_;
  if (pmb->block_size.nx3 > 1) {
    int j = north ? pmb->js : pmb->je+1;
    int count = 0;
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      Real val = 0.0;
      for (int n = 0; n < num_bufs; ++n)
        val += buf_list[n][count];
      for (int k = pmb->ks-NGHOST; k <= pmb->ke+NGHOST+1; ++k)
        pmb->pfield->e.x1e(k, j, i) = val / num_bufs;
      ++count;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ClearCoarseEMFBoundary()
//  \brief Clear the EMFs on the surface/edge contacting with a finer block

void FaceCenteredBoundaryVariable::ClearCoarseEMFBoundary() {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int nl;
  // face
  for (int n=0; n<pbval_->nface_; n++) {
    if (n==INNER_X1 || n==OUTER_X1) {
      int i;
      if (n==INNER_X1) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      nl=pbval_->nblevel[1][1][2*n];
      if (nl>pmb->loc.level) { // finer
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)=0.0;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i)=0.0;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i)=e2(pmb->ks+1,j,i)=0.0;
          for (int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i)=0.0;
        } else { // 1D
          e2(pmb->ks,pmb->js,i)=e2(pmb->ks+1,pmb->js,i)=0.0;
          e3(pmb->ks,pmb->js,i)=e3(pmb->ks,pmb->js+1,i)=0.0;
        }
      }
    }
    if (n==INNER_X2 || n==OUTER_X2) {
      int j;
      if (n==INNER_X2) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      nl=pbval_->nblevel[1][2*n-4][1];
      if (nl>pmb->loc.level) { // finer
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i)=0.0;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i)=0.0;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i)=e1(pmb->ks+1,j,i)=0.0;
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i)=0.0;
        }
      }
    }
    if (n==INNER_X3 || n==OUTER_X3) {
      int k;
      if (n==INNER_X3) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      nl=pbval_->nblevel[2*n-8][1][1];
      if (nl>pmb->loc.level) { // finer
        // this is always 3D
        for (int j=pmb->js+1; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)=0.0;
        }
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i)=0.0;
        }
      }
    }
  }
  // edge
  for (int n=0; n<pbval_->nedge_; n++) {
    if (pbval_->edge_flag_[n]==true) continue;
    if (n>=0 && n<4) {
      int i, j;
      if ((n & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((n & 2)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      for (int k=pmb->ks; k<=pmb->ke; k++)
        e3(k,j,i)=0.0;
      // x1x3 edge
    } else if (n>=4 && n<8) {
      int i, k;
      if ((n & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((n & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      for (int j=pmb->js; j<=pmb->je; j++)
        e2(k,j,i)=0.0;
      // x2x3 edge
    } else if (n>=8 && n<12) {
      int k, j;
      if ((n & 1)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      if ((n & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)=0.0;
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::AverageEMFBoundary()
// \brief Set EMF boundary received from a block on the finer level

void FaceCenteredBoundaryVariable::AverageEMFBoundary() {
  MeshBlock *pmb=pmy_block_;
  AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  AthenaArray<Real> &e3=pmb->pfield->e.x3e;
  int nl;
  // face
  for (int n=0; n<pbval_->nface_; n++) {
    if ((pbval_->block_bcs[n] != BLOCK_BNDRY) && (pbval_->block_bcs[n] != PERIODIC_BNDRY)
        && (pbval_->block_bcs[n] != POLAR_BNDRY)) continue;
    if (n==INNER_X1 || n==OUTER_X1) {
      int i;
      if (n==INNER_X1) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      nl=pbval_->nblevel[1][1][2*n];
      if (nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i)*=0.5;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i)*=0.5;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i)*=0.5, e2(pmb->ks+1,j,i)*=0.5;
          for (int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i)*=0.5;
        } else { // 1D
          e2(pmb->ks,pmb->js,i)*=0.5, e2(pmb->ks+1,pmb->js,i)*=0.5;
          e3(pmb->ks,pmb->js,i)*=0.5, e3(pmb->ks,pmb->js+1,i)*=0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          int k=pmb->ks+pmb->block_size.nx3/2;
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(k,j,i)*=0.5;
        }
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          int j=pmb->js+pmb->block_size.nx2/2;
          for (int k=pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i)*=0.5;
        }
      }
    }
    if (n==INNER_X2 || n==OUTER_X2) {
      int j;
      if (n==INNER_X2) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      nl=pbval_->nblevel[1][2*n-4][1];
      if (nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        if (pmb->block_size.nx3 > 1) {
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i)*=0.5;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i)*=0.5;
          }
        } else if (pmb->block_size.nx2 > 1) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i)*=0.5, e1(pmb->ks+1,j,i)*=0.5;
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i)*=0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          int k=pmb->ks+pmb->block_size.nx3/2;
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)*=0.5;
        }
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          int i=pmb->is+pmb->block_size.nx1/2;
          for (int k=pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i)*=0.5;
        }
      }
    }
    if (n==INNER_X3 || n==OUTER_X3) {
      int k;
      if (n==INNER_X3) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      nl=pbval_->nblevel[2*n-8][1][1];
      if (nl==pmb->loc.level) { // same ; divide all the face EMFs by 2
        for (int j=pmb->js+1; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i)*=0.5;
        }
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i)*=0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        // this is always 3D
        int j=pmb->js+pmb->block_size.nx2/2;
        for (int i=pmb->is; i<=pmb->ie; i++)
          e1(k,j,i)*=0.5;
        int i=pmb->is+pmb->block_size.nx1/2;
        for (int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i)*=0.5;
      }
    }
  }
  // edge
  for (int n=0; n<pbval_->nedge_; n++) {
    if (pbval_->nedge_fine_[n]==1) continue;
    Real div=1.0/static_cast<Real>(pbval_->nedge_fine_[n]);
    NeighborBlock& nb=pbval_->neighbor[n+6];
    Real half_div=div;
    if (nb.shear) half_div=0.5;
    // x1x2 edge (both 2D and 3D)
    if (n>=0 && n<4) {
      int i, j;
      if ((n & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((n & 2)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      for (int k=pmb->ks; k<=pmb->ke; k++)
        e3(k,j,i)*=half_div;
      // x1x3 edge
    } else if (n>=4 && n<8) {
      int i, k;
      if ((n & 1)==0) {
        i=pmb->is;
      } else {
        i=pmb->ie+1;
      }
      if ((n & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      for (int j=pmb->js; j<=pmb->je; j++)
        e2(k,j,i)*=half_div;
      // x2x3 edge
    } else if (n>=8 && n<12) {
      int j, k;
      if ((n & 1)==0) {
        j=pmb->js;
      } else {
        j=pmb->je+1;
      }
      if ((n & 2)==0) {
        k=pmb->ks;
      } else {
        k=pmb->ke+1;
      }
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)*=div;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlockEMF()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void FaceCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlockEMF() {
  MeshBlock *pmb=pmy_block_;
  if (pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1
      && pmb->block_size.nx3 > 1) {
    AthenaArray<Real> &e1=pmb->pfield->e.x1e;
    AthenaArray<Real> &e3=pmb->pfield->e.x3e;
    if (pbval_->block_bcs[INNER_X2]==POLAR_BNDRY
        || pbval_->block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE) {
      int j=pmb->js;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real tote1=0.0;
        for (int k=pmb->ks; k<=pmb->ke; k++)
          tote1 += e1(k,j,i);
        Real e1a=tote1/static_cast<double>(pmb->ke-pmb->ks+1);
        for (int k=pmb->ks; k<=pmb->ke+1; k++)
          e1(k,j,i) = e1a;
      }
      for (int i=pmb->is; i<=pmb->ie+1; i++) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          pbval_->azimuthal_shift_(k) = e3(k,j,i);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          int k_shift = k;
          k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
          e3(k,j,i) = pbval_->azimuthal_shift_(k_shift);
        }
      }
    }

    if (pbval_->block_bcs[OUTER_X2]==POLAR_BNDRY
        || pbval_->block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE) {
      int j=pmb->je+1;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real tote1=0.0;
        for (int k=pmb->ks; k<=pmb->ke; ++k)
          tote1 += e1(k,j,i);
        Real e1a=tote1/static_cast<double>(pmb->ke-pmb->ks+1);
        for (int k=pmb->ks; k<=pmb->ke+1; ++k)
          e1(k,j,i) = e1a;
      }
      for (int i=pmb->is; i<=pmb->ie+1; i++) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          pbval_->azimuthal_shift_(k) = e3(k,j,i);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          int k_shift = k;
          k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
          e3(k,j,i) = pbval_->azimuthal_shift_(k_shift);
        }
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReceiveFluxCorrection()
//  \brief Receive and Apply the surface EMF to the coarse neighbor(s) if needed

bool FaceCenteredBoundaryVariable::ReceiveFluxCorrection() {
  MeshBlock *pmb=pmy_block_;
  bool flag=true;

  // Receive same-level non-polar EMF values
  if (pbval_->firsttime_==true) {
    for (int n=0; n<pbval_->nneighbor; n++) { // first correct the same level
      NeighborBlock& nb = pbval_->neighbor[n];
      if (nb.type!=NEIGHBOR_FACE && nb.type!=NEIGHBOR_EDGE) break;
      if (nb.level!=pmb->loc.level) continue;
      if ((nb.type==NEIGHBOR_FACE) || ((nb.type==NEIGHBOR_EDGE) &&
                                       (pbval_->edge_flag_[nb.eid]==true))) {
        if (bd_var_flcor_.flag[nb.bufid]==BNDRY_COMPLETED) continue;
        if (bd_var_flcor_.flag[nb.bufid]==BNDRY_WAITING) {
          if (nb.rank==Globals::my_rank) { // on the same process
            flag=false;
            continue;
          }
#ifdef MPI_PARALLEL
          else { // NOLINT
            int test;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                       MPI_STATUS_IGNORE);
            MPI_Test(&(bd_var_flcor_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
            if (static_cast<bool>(test)==false) {
              flag=false;
              continue;
            }
            bd_var_flcor_.flag[nb.bufid] = BNDRY_ARRIVED;
          }
#endif
        }
        // boundary arrived; apply EMF correction
        SetEMFBoundarySameLevel(bd_var_flcor_.recv[nb.bufid], nb);
        bd_var_flcor_.flag[nb.bufid] = BNDRY_COMPLETED;
      }
    }
    if (flag==false) return flag;
    if (pmb->pmy_mesh->multilevel==true)
      ClearCoarseEMFBoundary();
    pbval_->firsttime_=false;
  }

  // Receive finer non-polar EMF values
  if (pmb->pmy_mesh->multilevel==true) {
    for (int n=0; n<pbval_->nneighbor; n++) { // then from finer
      NeighborBlock& nb = pbval_->neighbor[n];
      if (nb.type!=NEIGHBOR_FACE && nb.type!=NEIGHBOR_EDGE) break;
      if (nb.level!=pmb->loc.level+1) continue;
      if (bd_var_flcor_.flag[nb.bufid]==BNDRY_COMPLETED) continue;
      if (bd_var_flcor_.flag[nb.bufid]==BNDRY_WAITING) {
        if (nb.rank==Globals::my_rank) {// on the same process
          flag=false;
          continue;
        }
#ifdef MPI_PARALLEL
        else { // NOLINT
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&(bd_var_flcor_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
          if (static_cast<bool>(test)==false) {
            flag=false;
            continue;
          }
          bd_var_flcor_.flag[nb.bufid] = BNDRY_ARRIVED;
        }
#endif
      }
      // boundary arrived; apply EMF correction
      SetEMFBoundaryFromFiner(bd_var_flcor_.recv[nb.bufid], nb);
      bd_var_flcor_.flag[nb.bufid] = BNDRY_COMPLETED;
    }
  }

  // Receive polar EMF values
  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pbval_->polar_neighbor_north[n];
    if (emf_north_flag_[n] == BNDRY_WAITING) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int recv_flag;
        MPI_Test(&req_emf_north_recv_[n], &recv_flag, MPI_STATUS_IGNORE);
        if (!recv_flag) {
          flag = false;
          continue;
        }
        emf_north_flag_[n] = BNDRY_ARRIVED;
      }
#endif
    }
  }
  for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pbval_->polar_neighbor_south[n];
    if (emf_south_flag_[n] == BNDRY_WAITING) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int recv_flag;
        MPI_Test(&req_emf_south_recv_[n], &recv_flag, MPI_STATUS_IGNORE);
        if (!recv_flag) {
          flag = false;
          continue;
        }
        emf_south_flag_[n] = BNDRY_ARRIVED;
      }
#endif
    }
  }

  if (flag==true) {
    AverageEMFBoundary();
    if (pbval_->num_north_polar_blocks_ > 0)
      SetEMFBoundaryPolar(emf_north_recv_, pbval_->num_north_polar_blocks_, true);
    for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n)
      emf_north_flag_[n] = BNDRY_COMPLETED;
    if (pbval_->num_south_polar_blocks_ > 0)
      SetEMFBoundaryPolar(emf_south_recv_, pbval_->num_south_polar_blocks_, false);
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n)
      emf_south_flag_[n] = BNDRY_COMPLETED;
    if (pbval_->block_bcs[INNER_X2]==POLAR_BNDRY
        || pbval_->block_bcs[OUTER_X2]==POLAR_BNDRY
        || pbval_->block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE
        || pbval_->block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE)
      PolarBoundarySingleAzimuthalBlockEMF();
  }
  return flag;
}
