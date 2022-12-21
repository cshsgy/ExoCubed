// Athena++ header
#include <parameter_input.hpp>
#include <mesh/mesh.hpp>
#include <outputs/outputs.hpp>
#include <utils/utils.hpp>
#include <globals.hpp>

// debugger header
#include <debugger.hpp>

// cliutils header
#include <cliutils/command_line.hpp>
#include <cliutils/program_setup.hpp>

// canoe headers
#include <configure.hpp>
//#include "task_list/implicit_hydro_tasks.hpp"

// MPI/OpenMP headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

void setup_input_and_init_mesh(ParameterInput* &pinput, Mesh* &pmesh);

int main(int argc, char *argv[]) {

  program_start(argc, argv);

  ParameterInput *pinput;
  Mesh *pmesh;

  setup_input_and_init_mesh(pinput, pmesh);

  // Construct and initialize TaskList

  TimeIntegratorTaskList  *ptlist;
  try {
    ptlist = new TimeIntegratorTaskList(pinput, pmesh);
  }
  catch(std::bad_alloc& ba) {
    Debugger::Fatal("canoe", "memory allocation failed");
    return(0);
  }

  // Change to run directory, initialize outputs object, and make output of ICs

  Outputs *pouts;
  try {
    ChangeRunDir(CommandLine::prundir);
    pouts = new Outputs(pmesh, pinput);
    if (CommandLine::res_flag == 0)
      pouts->MakeOutputs(pmesh, pinput);
  }
  catch(std::bad_alloc& ba) {
    Debugger::Fatal("canoe", "memory allocation failed");
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    return(0);
  }

  // START OF MAIN INTEGRATION LOOP 
  // For performance, there is no error handler protecting this step (except outputs)

  if (Globals::my_rank == 0) {
    std::cout << "\nSetup complete, entering main loop...\n" << std::endl;
  }

  while ((pmesh->time < pmesh->tlim) &&
         (pmesh->nlim < 0 || pmesh->ncycle < pmesh->nlim)) {
    if (Globals::my_rank == 0)
      pmesh->OutputCycleDiagnostics();

    for (int stage=1; stage<=ptlist->nstages; ++stage) {
      ptlist->DoTaskListOneStage(pmesh, stage);
    }

    pmesh->UserWorkInLoop();

    pmesh->ncycle++;
    pmesh->time += pmesh->dt;
    Globals::mbcnt += pmesh->nbtotal;
    pmesh->step_since_lb++;

    pmesh->LoadBalancingAndAdaptiveMeshRefinement(pinput);

    pmesh->NewTimeStep();

    try {
      if (pmesh->time < pmesh->tlim) // skip the final output as it happens later
        pouts->MakeOutputs(pmesh,pinput);
    }
    catch(std::bad_alloc& ba) {
      Debugger::Fatal("canoe", "memory allocation failed");
      return(0);
    }
    catch(std::exception const& ex) {
      std::cout << ex.what() << std::endl;  // prints diagnostic message
#ifdef MPI_PARALLEL
      MPI_Finalize();
#endif
      return(0);
    }

    // check for signals
    if (SignalHandler::CheckSignalFlags() != 0) {
      break;
    }
  } 

  // END OF MAIN INTEGRATION LOOP 
  // Make final outputs, print diagnostics, clean up and terminate

  // Output the final cycle diagnostics and make the final outputs

  if (Globals::my_rank == 0)
    pmesh->OutputCycleDiagnostics();

  pmesh->UserWorkAfterLoop(pinput);

  try {
    pouts->MakeOutputs(pmesh,pinput,true);
  }
  catch(std::bad_alloc& ba) {
    Debugger::Fatal("canoe", "memory allocation failed");
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }

  // Print diagnostic messages related to the end of the simulation
  program_end(pmesh);

  delete pinput;
  delete pmesh;
  delete ptlist;
  delete pouts;

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif
  return(0);
}
