// C/C++ headers
#include <cstring>

// Athena++ headers
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <debugger/debugger.hpp>
#include <snap/meshblock_impl.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>
#include <utils/vectorize.hpp>

#include "profile_inversion.hpp"

void new_inversion_queue(std::vector<Inversion *> &fitq, MeshBlock *pmb,
                         ParameterInput *pin) {
  std::string str = pin->GetOrAddString("inversion", "tasks", "");
  std::vector<std::string> task_names =
      Vectorize<std::string>(str.c_str(), " ,");

  Inversion *pfit;
  for (auto p : task_names) {
    if (p == "VLAProfileInversion") {
      pfit = new VLAProfileInversion(pmb, pin);
    } else if (p == "JunoProfileInversion") {
      pfit = new JunoProfileInversion(pmb, pin);
    } else if (p == "VLACompositionInversion") {
    } else if (p == "JunoCompositionInversion") {
    } else {
      Debugger::Fatal("new_inversion_queue", "task::" + p, "unrecognized");
    }
    fitq.push_back(pfit);
  }

  int jl = pmb->js;
  for (auto q : fitq) {
    q->InitializePositions();
    q->setX2Indices(jl);
    jl += q->getX2Span();
  }
}