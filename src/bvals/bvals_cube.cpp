// Utilities for finding locations for cubed sphere 

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"
#include "../mesh/meshblock_tree.hpp"

void TransformOxForCubedSphere(LogicalLocation loc, int *ox1, int *ox2, int *ox3, MeshBlockTree &tree){

}