
namespace CubedSphere {

#define DBL_EPSILON 1.0e-10

void GetLatLon(Real *lat, Real *lon, Coordinates *pcoord, int k, int j, int i) {
  // Obtain Lat and Lon (radians) from x2 and x3
  // k is not used for now
  // Find the block number
  LogicalLocation loc = pcoord->pmy_block->loc;
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
  // Calculate the needed parameters
  Real dX = tan(pcoord->x2v(j));
  Real dY = tan(pcoord->x3v(k));
  RLLFromXYP(dY, -dX, blockID - 1, *lon, *lat);
}

void GetLatLonFace2(Real *lat, Real *lon, Coordinates *pcoord, int k, int j,
                    int i) {
  // Obtain Lat and Lon (radians) from x2 and x3
  // k is not used for now
  // Find the block number
  LogicalLocation loc = pcoord->pmy_block->loc;
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
  // Calculate the needed parameters
  Real dX = tan(pcoord->x2f(j));
  Real dY = tan(pcoord->x3v(k));
  RLLFromXYP(dY, -dX, blockID - 1, *lon, *lat);
}

void GetLatLonFace3(Real *lat, Real *lon, Coordinates *pcoord, int k, int j,
                    int i) {
  // Obtain Lat and Lon (radians) from x2 and x3
  // k is not used for now
  // Find the block number
  LogicalLocation loc = pcoord->pmy_block->loc;
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
  // Calculate the needed parameters
  Real dX = tan(pcoord->x2v(j));
  Real dY = tan(pcoord->x3f(k));
  RLLFromXYP(dY, -dX, blockID - 1, *lon, *lat);
}

}  // namespace CubedSphere
