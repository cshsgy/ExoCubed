// Athena++ headers
#include <parameter_input.hpp>
#include <mesh/mesh.hpp>
//#include <globals.hpp>

// debugger headers
#include <debugger.hpp>

// cliutils header
#include <cliutils/command_line.hpp>
#include <cliutils/program_setup.hpp>

// canoe headers
#include <configure.hpp>
#include "mesh/meshblock_impl.hpp"

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

using CL = CommandLine;

void setup_input_and_init_mesh(ParameterInput* &pinput, Mesh* &pmesh)
{
  IOWrapper infile, restartfile;

  try {
    pinput = new ParameterInput;

    if (CL::res_flag == 1) {
      restartfile.Open(CL::restart_filename, IOWrapper::FileMode::read);
      pinput->LoadFromFile(restartfile);
      // If both -r and -i are specified, make sure next_time gets corrected.
      // This needs to be corrected on the restart file because we need the old dt.
      if (CL::iarg_flag == 1) pinput->RollbackNextTime();
      // leave the restart file open for later use
    }

    if (CL::iarg_flag == 1) {
      // if both -r and -i are specified, override the parameters using the input file
      infile.Open(CL::input_filename, IOWrapper::FileMode::read);
      pinput->LoadFromFile(infile);
      infile.Close();
    }
    pinput->ModifyFromCmdline(CL::argc ,CL::argv);
  } catch (std::bad_alloc& ba) {
    if (CL::res_flag == 1) restartfile.Close();
    Debugger::Fatal("main", "memory allocation failed initializing class ParameterInput:", ba.what());
  } catch(std::exception const& ex) {
    if (CL::res_flag == 1) restartfile.Close();
    Debugger::Fatal("main", ex.what());
  }

  try {
    if (CL::res_flag == 0) {
      pmesh = new Mesh(pinput, CL::mesh_flag);
    } else {
      pmesh = new Mesh(pinput, restartfile, CL::mesh_flag);
    }
  } catch (std::bad_alloc& ba) {
    if (CL::res_flag == 1) restartfile.Close();
    Debugger::Fatal("main", "memory allocation failed initializing class Mesh:", ba.what());
  } catch (std::exception const& ex) {
    if (CL::res_flag == 1) restartfile.Close();
    Debugger::Fatal("main", ex.what());
  }

  // With current mesh time possibly read from restart file, correct next_time for outputs
  if (CL::iarg_flag == 1 && CL::res_flag == 1) {
    // if both -r and -i are specified, ensure that next_time  >= mesh_time - dt
    pinput->ForwardNextTime(pmesh->time);
  }

  // Dump input parameters and quit if code was run with -n option.
  if (CL::narg_flag) {
    if (Globals::my_rank == 0) pinput->ParameterDump(std::cout);
    if (CL::res_flag == 1) restartfile.Close();
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    exit(0);
  }

  // everything works fine here, close the restart file
  if (CL::res_flag == 1) restartfile.Close();

  // Quit if -m was on cmdline.  This option builds and outputs mesh structure.
  if (CL::mesh_flag > 0) {
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    exit(0);
  }

  // set up additional components
  for (int b = 0; b < pmesh->nblocal; ++b) {
    MeshBlock *pmb = pmesh->my_blocks(b);
    pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
  }

  // initialize mesh
  pmesh->Initialize(CL::res_flag, pinput);
}
