
template <typename T>
void PlanarExchanger::setColor_(int *color, CoordinateDirection dir) {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock bblock, tblock;
  pmb->FindNeighbors(dir, bblock, tblock);

  if (dir == X1DIR) {
    if (pmb->block_size.x1min <= pmb->pmy_mesh->mesh_size.x1min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x1max >= pmb->pmy_mesh->mesh_size.x1max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  } else if (dir == X2DIR) {
    if (pmb->block_size.x2min <= pmb->pmy_mesh->mesh_size.x2min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x2max >= pmb->pmy_mesh->mesh_size.x2max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  } else {  // X3DIR
    if (pmb->block_size.x3min <= pmb->pmy_mesh->mesh_size.x3min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x3max >= pmb->pmy_mesh->mesh_size.x3max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allgather(&bblock.snb.rank, 1, MPI_INT, brank_.data(), 1, MPI_INT,
                MPI_COMM_WORLD);
#else
  brank_[0] = -1;
#endif

  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    // color[i] = brank_[i] == -1 ? color[i] : color[brank_[i]];
    if (brank_[i] == -1) {
      if (color[i] == -1) color[i] = c++;
    } else
      color[i] = color[brank_[i]];
  }
}

void Diagnostics::gatherAllData23_(AthenaArray<Real> &total_vol,
                                   AthenaArray<Real> &total_data) {
  MeshBlock *pmb = pmy_block_;

  // calculate total volume
  total_vol.ZeroClear();
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
      for (int i = pmb->is; i <= pmb->ie; ++i) total_vol(i) += vol_(i);
    }

    // sum over all ranks
#ifdef MPI_PARALLEL
  MPI_Comm comm;
  std::fill(color_.data(), color_.data() + Globals::nranks, -1);
  setColor_(color_.data(), X2DIR);
  setColor_(color_.data(), X3DIR);
  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank,
                 &comm);
  // int size;
  // MPI_Comm_size(comm, &size);
  // std::cout << size << std::endl;
  MPI_Allreduce(MPI_IN_PLACE, total_vol.data(), total_vol.GetSize(),
                MPI_ATHENA_REAL, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, total_data.data(), total_data.GetSize(),
                MPI_ATHENA_REAL, MPI_SUM, comm);
  MPI_Comm_free(&comm);
#endif
}
