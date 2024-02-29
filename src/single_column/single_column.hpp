// C/C++
#include <array>
#include <vector>

// athena
#include <athena/athena.hpp>

// canoe
#include <virtual_groups.hpp>

class MeshBlock;

class SingleColumn : public ParameterGroup {
 public:
  friend void find_tp_bottom(int n, double *x, double *f, void *arg);

  SingleColumn(MeshBlock *pmb, ParameterInput *pin);
  virtual ~SingleColumn();

  void ConvectiveAdjustment(MeshBlock *pmb, int k, int j);

 protected:  // convective adjustment functions
  std::array<Real, 2> findTPBottom(Real tempf, Real presf, Real mass0,
                                   Real enthalpy0, MeshBlock *pmb, int k, int j,
                                   int il, int iu);
  std::array<int, 2> findUnstableRange(MeshBlock *pmb, int k, int j);

 protected:
  // scrach arrays
  AthenaArray<Real> vol_;
};
