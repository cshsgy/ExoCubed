// athena
#include <athena/coordinates/coordinates.hpp>

// canoe
#include "air_parcel.hpp"

void convective_adjustment(std::vector<AirParcel>& air_column,
                           Coordinates* pcoord, Real grav, int k, int j, int il,
                           int iu);

Real GetTheta(AirParcel const& air);

void recursively_search_convective_adjustment(
    std::vector<AirParcel>& air_column, Coordinates* pcoord, Real grav, int k,
    int j);
