#include <parameter_input.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>
#include <configure.hpp>
#include <cubed_sphere.hpp>
#include "../src/eos/eos.hpp"
#include "../src/field/field.hpp"


void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Input field is a uniform 10 m/s zonal wind
  Real U, V;
  Real Vy, Vz;

  // Constants of the test
  Real a = 6.37122e6;
  Real R = 4.0;
  Real omega = 7.292e-5;
  Real g = 9.80616;
  Real K = 7.848e-6;
  Real h0 = 8000.0;

  std::cout << "======Start of Problem generator======" << loc.lx2 << "||" << loc.lx3 <<std::endl;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x = tan(pcoord->x2v(j));
        Real y = tan(pcoord->x3v(k));
        Real lat, lon;
        GetLatLon(&lat, &lon, pcoord, k, j, i);
        Real rad = (PI/2.0-lat)*R;
        Real f = omega * sin(lat);
        
        U = a*K*cos(lat)+a*K*pow(cos(lat),R-1.0)*(R*sin(lat)*sin(lat)-cos(lat)*cos(lat))*cos(R*lon);
        V = -a*K*R*pow(cos(lat),R-1.0)*sin(R*lon)*sin(lat);
        Real A = K/2.0*(2.0*omega+K)*cos(lat)*cos(lat)+0.25*K*K*pow(cos(lat),R*2.0)*((R+1.0)*cos(lat)*cos(lat)+(2.0*R*R-R-2.0)-2.0*R*R*pow(cos(lat),-2.0));
        Real B = 2.0*(omega+K)*K/(R+1.0)/(R+2.0)*pow(cos(lat),R)*((R*R+2.0*R+2.0)-(R*R+2.0*R+1.0)*cos(lat)*cos(lat));
        Real C = 0.25*K*K*pow(cos(lat),R*2.0)*((R+1.0)*cos(lat)*cos(lat)-(R+2.0));
        Real dphi = a*a*A+a*a*B*cos(R*lon)+a*a*C*cos(2.0*R*lon);

        GetVyVz(&Vy, &Vz, pcoord, U, V, k, j, i);

        phydro->w(IDN,k,j,i) = g*h0+dphi;
        phydro->w(IVY,k,j,i) = Vy;
        phydro->w(IVZ,k,j,i) = Vz;
      }
  std::cout << "Done with problem generator..." << std::endl;
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  // Calculate lat, lon, U, V for output and visualization
  for (int k = ks-NGHOST; k <= ke+NGHOST; ++k)
    for (int j = js-NGHOST; j <= je+NGHOST; ++j)
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
        Real Vy = phydro->w(IVY,k,j,i);
        Real Vz = phydro->w(IVZ,k,j,i);
        Real lat, lon;
        GetLatLon(&lat, &lon, pcoord, k, j, i);
        user_out_var(0,k,j,i) = lat/PI*180.0;
        user_out_var(1,k,j,i) = lon/PI*180.0;
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1*x1 + x2*x2);
        Real f = 0.0;
        Real U, V;
        GetUV(&U, &V, pcoord, Vy, Vz, k, j, i);
        user_out_var(2,k,j,i) = U;
        user_out_var(3,k,j,i) = V;
      }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(4);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
}


void CoriolisForcing(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &w, const AthenaArray<Real> &prim_scalar, AthenaArray<Real> const &bcc, 
  AthenaArray<Real> &u, AthenaArray<Real> &cons_scalar)
{
    int is = pmb->is;
    Real omega = 7.292e-5;
    for (int k = pmb->ks; k<=pmb->ke; ++k)
        for (int j = pmb->js; j <= pmb->je; ++j){
          Real lat, lon;
          GetLatLon(&lat, &lon, pmb->pcoord, k, j, is);
          // coriolis force
          Real f = 2.*omega*sin(lat);
          Real U, V;
          GetUV(&U, &V, pmb->pcoord, w(IVY,k,j,is), w(IVZ,k,j,is), k, j, is);
          Real ll_acc_U = -f*V;
          Real ll_acc_V = f*U;
          Real acc2, acc3;
          GetVyVz(&acc2, &acc3, pmb->pcoord, ll_acc_U, ll_acc_V, k, j, is);
          u(IM2,k,j,is) += dt*w(IDN,k,j,is)*acc2;
          u(IM3,k,j,is) += dt*w(IDN,k,j,is)*acc3;
    }
}
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitSourceFunction(CoriolisForcing);
}
