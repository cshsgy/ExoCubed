#include <sstream>
#include <iostream>

#include <globals.hpp>
#include <hydro/hydro.hpp>
#include <cubed_sphere.hpp>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

void Hydro::SaveLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
  int direction,  int k, int j, int il, int iu) {
    for (int n=0; n<NWAVE; n++){
        for (int i=il; i<=iu; i++){
            L3DValues[direction](n, k, j, i) = L_in(n, i);
            R3DValues[direction](n, k, j, i) = R_in(n, i);
        }
    }
}

void Hydro::LoadLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
  int direction,  int k, int j, int il, int iu) {
    for (int n=0; n<NWAVE; n++){
        for (int i=il; i<=iu; i++){
            L_in(n, i) = L3DValues[direction](n, k, j, i);
            R_in(n, i) = R3DValues[direction](n, k, j, i);
        }
    }
}

void Hydro::SynchronizeFluxesSend(){
    MeshBlock *pmb = pmy_block;
    for (int n=0; n<pmb->pbval->nneighbor; n++){
        NeighborBlock &nb = pmb->pbval->neighbor[n];
        if(nb.ni.ox1==0 && nb.ni.ox2*nb.ni.ox3==0){ // On x2 and x3 face boundaries only
            SendNeighborBlocks(pmb->loc, nb.ni.ox2, nb.ni.ox3, nb.snb.rank, nb.snb.gid);
        }
    }
}

void Hydro::SynchronizeFluxesRecv(){
    MeshBlock *pmb = pmy_block;
    for (int n=0; n<pmb->pbval->nneighbor; n++){
        NeighborBlock &nb = pmb->pbval->neighbor[n];
        if(nb.ni.ox1==0 && nb.ni.ox2*nb.ni.ox3==0){ // On x2 and x3 face boundaries only
            RecvNeighborBlocks(pmb->loc, nb.ni.ox2, nb.ni.ox3, nb.snb.rank, nb.snb.gid);
        }
    }
}

void Hydro::SendNeighborBlocks(LogicalLocation const& loc, int ox2, int ox3, int tg_rank, int tg_gid){
    MeshBlock *pmb = pmy_block;

    int lv2_lx2 = loc.lx2 >> (loc.level - 2);
    int lv2_lx3 = loc.lx3 >> (loc.level - 2);
    int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
    // Calculate local ID
    int local_lx2 = loc.lx2 - (lv2_lx2<<(loc.level - 2));
    int local_lx3 = loc.lx3 - (lv2_lx3<<(loc.level - 2));
    int bound_lim = (1<<(loc.level - 2)) - 1;
    int tox2, tox3;
    int DirTag, DirNum, ownTag; // Tag for target, and the numbering of axis
    int ox2_bkp = ox2;
    int ox3_bkp = ox3;
    bool invDir, Left; // Marking whether need to reverse direction in packing or unpacking; Left marking left or right.
    // Hard code the boundary cases
    if(local_lx2==bound_lim && ox2==1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = true;
        // Determine whether need to inverse direction
        switch(blockID){
            case 1: invDir = false; break;
            case 2: invDir = false; break;
            case 3: invDir = true; break;
            case 4: invDir = false; break;
            case 5: invDir = true; break;
            case 6: invDir = true; break;
        }
    }else if(local_lx2==0 && ox2==-1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = false;
        switch(blockID){
            case 1: invDir = true; break;
            case 2: invDir = false; break;
            case 3: invDir = false; break;
            case 4: invDir = true; break;
            case 5: invDir = false; break;
            case 6: invDir = true; break;
        }
    }else if(local_lx3==bound_lim && ox3==1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = true;
        switch(blockID){
            case 1: invDir = true; break;
            case 2: invDir = false; break;
            case 3: invDir = false; break;
            case 4: invDir = false; break;
            case 5: invDir = false; break;
            case 6: invDir = false; break;
        }
    }else if(local_lx3==0 && ox3==-1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = false;
        switch(blockID){
            case 1: invDir = false; break;
            case 2: invDir = false; break;
            case 3: invDir = false; break;
            case 4: invDir = false; break;
            case 5: invDir = true; break;
            case 6: invDir = false; break;
        }
    }else{ // No need to communicate fluxes, return
        return;
    }
    // Pack the data
    int kb1, kb2, jb1, jb2, ib1, ib2;
    if (ox2==1){
        DirNum = X2DIR;
        jb1 = pmb->je + 1;
        jb2 = pmb->je + 1;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ks;
        kb2 = pmb->ke + 1;
    }
    if (ox2==-1){
        DirNum = X2DIR;
        jb1 = pmb->js;
        jb2 = pmb->js;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ks;
        kb2 = pmb->ke + 1;
    }
    if (ox3==1){
        DirNum = X3DIR;
        jb1 = pmb->js;
        jb2 = pmb->je + 1;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ke + 1;
        kb2 = pmb->ke + 1;
    }
    if (ox3==-1){
        DirNum = X3DIR;
        jb1 = pmb->js;
        jb2 = pmb->je + 1;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ks;
        kb2 = pmb->ks;
    }
    int dsize = ((kb2 - kb1 + 1) * (jb2 - jb1 + 1) * (ib2 - ib1 + 1) * NWAVE);
    Real *data = new Real[dsize];
    int offset = 0;
    if (invDir){
        for (int n=0; n<NWAVE; n++)
            for (int k=kb2; k>=kb1; k--)
                for (int j=jb2; j>=jb1; j--)
                    for (int i=ib1; i<=ib2; i++)
                        if (Left)
                            data[offset++] = L3DValues[DirNum](n,k,j,i);
                        else
                            data[offset++] = R3DValues[DirNum](k,j,i);
    }else{
        for (int n=0; n<NWAVE; n++)
            for (int k=kb1; k<=kb2; k++)
                for (int j=jb1; j<=jb2; j++)
                    for (int i=ib1; i<=ib2; i++)
                        if (Left)
                            data[offset++] = L3DValues[DirNum](n,k,j,i);
                        else
                            data[offset++] = R3DValues[DirNum](n,k,j,i);
    }

    // Calculate the tag of destination
    if (tox2==-1) DirTag = 0 + 4 * pmb->gid + 24*(1 << (loc.level - 2)) * tg_gid;
    if (tox2==1) DirTag = 1 + 4 * pmb->gid + 24*(1 << (loc.level - 2)) * tg_gid;
    if (tox3==-1) DirTag = 2 + 4 * pmb->gid + 24*(1 << (loc.level - 2)) * tg_gid;
    if (tox3==1) DirTag = 3 + 4 * pmb->gid + 24*(1 << (loc.level - 2)) * tg_gid;
    // Send by MPI: we don't care whether it is in the same process for now
    if (ox2==-1) ownTag = 0;
    if (ox2==1) ownTag = 1;
    if (ox3==-1) ownTag = 2;
    if (ox3==1) ownTag = 3;
    MPI_Isend(data, dsize, MPI_DOUBLE, tg_rank, DirTag, MPI_COMM_WORLD, &send_request[ownTag]);
    std::cout << "===============================" << std::endl;
    std::cout << "MPI Message: Sent data with size " << dsize << " from rank " << Globals::my_rank << " to " << tg_rank << " on tag number " << DirTag << std::endl;
    std::cout << "This is a message from block gid " << pmb->gid << " to " << tg_gid << std::endl;
    std::cout << "The direction is ox2=" << ox2 << ", ox3=" << ox3 << " and inversion is " << invDir << std::endl;
    delete[] data;
}

void Hydro::RecvNeighborBlocks(LogicalLocation const& loc, int ox2, int ox3, int tg_rank, int tg_gid){
    MeshBlock *pmb = pmy_block;

    int lv2_lx2 = loc.lx2 >> (loc.level - 2);
    int lv2_lx3 = loc.lx3 >> (loc.level - 2);
    int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
    // Calculate local ID
    int local_lx2 = loc.lx2 - (lv2_lx2<<(loc.level - 2));
    int local_lx3 = loc.lx3 - (lv2_lx3<<(loc.level - 2));
    int bound_lim = (1<<(loc.level - 2)) - 1;
    int tox2, tox3;
    int DirTag, DirNum, ownTag; // Tag for receiving, and the numbering of axis to place the values
    int ox2_bkp = ox2;
    int ox3_bkp = ox3;
    bool invDir, Left; // Marking whether need to reverse direction in packing or unpacking; Left marking left or right.
    // Hard code the boundary cases
    if(local_lx2==bound_lim && ox2==1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = false;
        // Determine whether need to inverse direction
    }else if(local_lx2==0 && ox2==-1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = true;
    }else if(local_lx3==bound_lim && ox3==1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = false;
    }else if(local_lx3==0 && ox3==-1){
        TransformOxForCubedSphere(&ox2_bkp, &ox3_bkp, &tox2, &tox3, loc);
        Left = true;
    }else{ // No need to communicate fluxes, return
        return;
    }
    // Pack the data
    int kb1, kb2, jb1, jb2, ib1, ib2;
    if (ox2==1){
        DirNum = X2DIR;
        jb1 = pmb->je + 1;
        jb2 = pmb->je + 1;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ks;
        kb2 = pmb->ke + 1;
    }
    if (ox2==-1){
        DirNum = X2DIR;
        jb1 = pmb->js;
        jb2 = pmb->js;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ks;
        kb2 = pmb->ke + 1;
    }
    if (ox3==1){
        DirNum = X3DIR;
        jb1 = pmb->js;
        jb2 = pmb->je + 1;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ke + 1;
        kb2 = pmb->ke + 1;
    }
    if (ox3==-1){
        DirNum = X3DIR;
        jb1 = pmb->js;
        jb2 = pmb->je + 1;
        ib1 = pmb->is;
        ib2 = pmb->ie + 1;
        kb1 = pmb->ks;
        kb2 = pmb->ks;
    }
    int dsize = ((kb2 - kb1 + 1) * (jb2 - jb1 + 1) * (ib2 - ib1 + 1) * NWAVE);
    Real *data = new Real[dsize];
    // Calculate the tag for receiving
    if (ox2==-1) DirTag = 0 + 4 * tg_gid + 24*(1 << (loc.level - 2)) * pmb->gid;
    if (ox2==1) DirTag = 1 + 4 * tg_gid + 24*(1 << (loc.level - 2)) * pmb->gid;
    if (ox3==-1) DirTag = 2 + 4 * tg_gid + 24*(1 << (loc.level - 2)) * pmb->gid;
    if (ox3==1) DirTag = 3 + 4 * tg_gid + 24*(1 << (loc.level - 2)) * pmb->gid;
    MPI_Recv(data, dsize, MPI_DOUBLE, tg_rank, DirTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
    std::cout << "===============================" << std::endl;
    std::cout << "MPI Message: Received data from " << tg_gid << " to " << pmb->gid << " on tag number " << DirTag << std::endl;

    // =======
    // int test;
    // probe MPI communications.  This is a bit of black magic that seems to promote
    // communications to top of stack and gets them to complete more quickly
    // MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
    //             MPI_STATUS_IGNORE);
    // =======
    // MPI_Test(&recv_request[ownTag], &test, MPI_STATUS_IGNORE);
    

    int offset = 0;
    for (int n=0; n<NWAVE; n++)
        for (int k=kb1; k<=kb2; k++)
            for (int j=jb1; j<=jb2; j++)
                for (int i=ib1; i<=ib2; i++)
                    if (Left)
                        L3DValues[DirNum](n,k,j,i) = data[offset++];
                    else
                        R3DValues[DirNum](n,k,j,i) = data[offset++];
    delete[] data;
}
