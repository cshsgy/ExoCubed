// C/C++
#include <array>
#include <vector>

// athena
#include <athena/athena.hpp>

// canoe
#include <virtual_groups.hpp>

class SingleColumn : public ParameterGroup {
 public:
  SingleColumn(MeshBlock *pmb, ParameterInput *pin);
  virtual ~SingleColumn() {}

  void ConvectiveAdjustment(AthenaArray<Real> &u, int k, int j);

 protected:
  std::array<Real, 2> findTPBottom(AthenaArray<Real> const &u, int k, int j,
                                   int il, int iu);
  std::array<int, 2> findUnstableRange(AthenaArray<Real> const &u, int k,
                                       int j);

 protected:
  MeshBlock *pmy_block_;

  // scrach arrays
  AthenaArray<Real> vol_;
};

Real GetTheta(AirParcel const &air);

void recursively_search_convective_adjustment(
    std::vector<AirParcel> &air_column, Real *x1_ptr, Real grav, int k, int j);
