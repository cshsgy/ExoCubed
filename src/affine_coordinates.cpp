// Athena++ headers
#include <athena.hpp>
#include <athena_arrays.hpp>
#include <mesh/mesh.hpp>
#include <parameter_input.hpp>
#include <coordinates/coordinates.hpp>
#include <cubed_sphere.hpp>

//----------------------------------------------------------------------------------------
//! Cartesian coordinates constructor

AffineCoordinate::AffineCoordinate(MeshBlock *pmb, ParameterInput *pin, bool flag)
    : Coordinates(pmb, pin, flag)
{
  // initialize volume-averaged coordinates and spacing
  // x1-direction: x1v = dx/2
  for (int i=il-ng; i<=iu+ng; ++i) {
    x1v(i) = 0.5*(x1f(i+1) + x1f(i));
  }
  for (int i=il-ng; i<=iu+ng-1; ++i) {
    if (pmb->block_size.x1rat != 1.0) {
      dx1v(i) = x1v(i+1) - x1v(i);
    } else {
      // dx1v = dx1f constant for uniform mesh; may disagree with x1v(i+1) - x1v(i)
      dx1v(i) = dx1f(i);
    }
  }

  // x2-direction: x2v = dy/2
  if (pmb->block_size.nx2 == 1) {
    x2v(jl) = 0.5*(x2f(jl+1) + x2f(jl));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      x2v(j) = 0.5*(x2f(j+1) + x2f(j));
    }
    for (int j=jl-ng; j<=ju+ng-1; ++j) {
      if (pmb->block_size.x2rat != 1.0) {
        dx2v(j) = x2v(j+1) - x2v(j);
      } else {
        // dx2v = dx2f constant for uniform mesh; may disagree with x2v(j+1) - x2v(j)
        dx2v(j) = dx2f(j);
      }
    }
  }

  // x3-direction: x3v = dz/2
  if (pmb->block_size.nx3 == 1) {
    x3v(kl) = 0.5*(x3f(kl+1) + x3f(kl));
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k=kl-ng; k<=ku+ng; ++k) {
      x3v(k) = 0.5*(x3f(k+1) + x3f(k));
    }
    for (int k=kl-ng; k<=ku+ng-1; ++k) {
      if (pmb->block_size.x3rat != 1.0) {
        dx3v(k) = x3v(k+1) - x3v(k);
      } else {
        // dxkv = dx3f constant for uniform mesh; may disagree with x3v(k+1) - x3v(k)
        dx3v(k) = dx3f(k);
      }
    }
  }
  // initialize geometry coefficients
  // x1-direction
  for (int i=il-ng; i<=iu+ng; ++i) {
    h2v(i) = 1.0;
    h2f(i) = 1.0;
    h31v(i) = 1.0;
    h31f(i) = 1.0;
    dh2vd1(i) = 0.0;
    dh2fd1(i) = 0.0;
    dh31vd1(i) = 0.0;
    dh31fd1(i) = 0.0;
  }

  // x2-direction
  if (pmb->block_size.nx2 == 1) {
    h32v(jl) = 1.0;
    h32f(jl) = 1.0;
    dh32vd2(jl) = 0.0;
    dh32fd2(jl) = 0.0;
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      h32v(j) = 1.0;
      h32f(j) = 1.0;
      dh32vd2(j) = 0.0;
      dh32fd2(j) = 0.0;
    }
  }

  // initialize area-averaged coordinates used with MHD AMR
  if ((pmb->pmy_mesh->multilevel) && MAGNETIC_FIELDS_ENABLED) {
    for (int i=il-ng; i<=iu+ng; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j=jl-ng; j<=ju+ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k=kl-ng; k<=ku+ng; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }
}


void ProjectLocalCartesianAffine(const AthenaArray<Real> &src, 
        AthenaArray<Real> &tgt, Real affine_angle, int sn, int en, int si, int ei, int sj, 
        int ej, int sk, int ek, int Dir)
{
  // int ivy = IVX + ((ivx-IVX)+1)%3;
  // int ivz = IVX + ((ivx-IVX)+2)%3;
  for(int n=sn; n<=en; n++)
    for(int k=sk; k<=ek; k++)
      for(int j=sj; j<=ej; j++)
        for(int i=si; i<=ei; i++)
          tgt(n,k,j,i) = src(n,k,j,i); // To ensure that the sources are in the right location

  for(int k=sk; k<=ek; k++)
    for(int j=sj; j<=ej; j++)
      for(int i=si; i<=ei; i++){
        Real vx = src(IVX,k,j,i);
        Real vy = src(IVY,k,j,i);
        Real vz = src(IVZ,k,j,i);
        switch(Dir)
        {
          case X1DIR: // Projecting to X Direction
            tgt(IVX,k,j,i) = vx*sin(affine_angle);
            tgt(IVY,k,j,i) = vy + vx*cos(affine_angle);
            tgt(IVZ,k,j,i) = vz;
            break;
          case X2DIR: // Projecting to Y Direction
            tgt(IVX,k,j,i) = vx+vy*cos(affine_angle);
            tgt(IVY,k,j,i) = vy*sin(affine_angle);
            tgt(IVZ,k,j,i) = vz;
            break;
          case X3DIR: // Projecting to Z Direction, no need to project... (Both affine and cubed sphere not needed)
            break;
        }
      }
}

void DeProjectLocalCartesianAffine(AthenaArray<Real> &flux, Real affine_angle, int sn, int en, int si, int ei, int sj, 
        int ej, int sk, int ek, int Dir)
{
  for(int k=sk; k<=ek; k++)
    for(int j=sj; j<=ej; j++)
      for(int i=si; i<=ei; i++){
        Real fx = flux(IVX,k,j,i);
        Real fy = flux(IVY,k,j,i);
        Real fz = flux(IVZ,k,j,i);
        switch(Dir)
        {
          case X1DIR: // Projecting to X Direction
            flux(IVX,k,j,i) = fx/sin(affine_angle);
            flux(IVY,k,j,i) = fy - fx*cos(affine_angle)/sin(affine_angle);
            flux(IVZ,k,j,i) = fz;
            break;
          case X2DIR: // Projecting to Y Direction
            flux(IVX,k,j,i) = fx-fy*cos(affine_angle)/sin(affine_angle);
            flux(IVY,k,j,i) = fy/sin(affine_angle);
            flux(IVZ,k,j,i) = fz;
            break;
          case X3DIR: // Projecting to Z Direction, no need to project... (Both affine and cubed sphere not needed)
            break;
        }
      }
}

// Put in the changes in face2area etc, similar to cylindrical.cpp
// In affine coordinates, the differences lie in face1area and face2area. face3area cancels the sqrt(g) term exactly.

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

void AffineCoordinate::Face1Area(const int k, const int j, const int il, const int iu,
                            AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real& area_i = area(i);
    area_i = dx1f(i)*dx3f(k) / (sin(PI/3.0) * sin(PI/3.0));
  }
  return;
}

void AffineCoordinate::Face2Area(const int k, const int j, const int il, const int iu,
                            AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real& area_i = area(i);
    area_i = dx1f(i)*dx3f(k) / (sin(PI/3.0) * sin(PI/3.0));
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real AffineCoordinate::GetFace1Area(const int k, const int j, const int i) {
  return dx2f(j)*dx3f(k) / (sin(PI/3.0) * sin(PI/3.0));
}

Real AffineCoordinate::GetFace2Area(const int k, const int j, const int i) {
  return dx1f(i)*dx3f(k) / (sin(PI/3.0) * sin(PI/3.0));
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// VolCenterFaceXArea functions: compute area of face with normal in X-dir as vector
// where the faces are joined by cell centers (for non-ideal MHD)

void AffineCoordinate::VolCenterFace1Area(const int k, const int j, const int il, const int iu,
                                     AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real& area_i = area(i);
    area_i = dx2v(j)*dx3v(k) / (sin(PI/3.0) * sin(PI/3.0));
  }
  return;
}

void AffineCoordinate::VolCenterFace2Area(const int k, const int j, const int il, const int iu,
                                     AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real& area_i = area(i);
    area_i = dx1v(i)*dx3v(k) / (sin(PI/3.0) * sin(PI/3.0));
  }
  return;
}