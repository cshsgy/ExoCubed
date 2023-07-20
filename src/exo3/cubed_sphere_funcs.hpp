#ifndef SRC_EXO3_CUBED_SPHERE_FUNCS_HPP_
#define SRC_EXO3_CUBED_SPHERE_FUNCS_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

int FindBlockID(LogicalLocation const &loc);

void TransformOx(int *ox2, int *ox3, int *tox2, int *tox3,
                 LogicalLocation const &loc);

void PackData(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
              int ei, int sj, int ej, int sk, int ek, int &offset, int ox1,
              int ox2, int ox3, LogicalLocation const &loc);

Real MeshGeneratorX2(Real x, LogicalLocation const &loc);
Real MeshGeneratorX3(Real x, LogicalLocation const &loc);

// Helper functions
void GetLatLon(Real *lat, Real *lon, Coordinates *pcoord, int k, int j, int i);
void GetLatLonFace2(Real *lat, Real *lon, Coordinates *pcoord, int k, int j,
                    int i);
void GetLatLonFace3(Real *lat, Real *lon, Coordinates *pcoord, int k, int j,
                    int i);

void GetUV(Real *U, Real *V, Coordinates *pcoord, Real V2, Real V3, int k,
           int j, int i);
void GetVyVz(Real *V2, Real *V3, Coordinates *pcoord, Real U, Real V, int k,
             int j, int i);

// Helper functions adapted from Paul
void VecTransABPFromRLL(Real X, Real Y, int blockID, Real U, Real V, Real *V2,
                        Real *V3);
void VecTransRLLFromABP(Real X, Real Y, int blockID, Real V2, Real V3, Real *U,
                        Real *V);
void RLLFromXYP(Real dX, Real dY, int nP, Real &lon, Real &lat);
void XYPFromRLL(Real lon, Real lat, Real &dX, Real &dY, int &nP);

#endif  // SRC_EXO3_CUBED_SPHERE_FUNCS_HPP_
