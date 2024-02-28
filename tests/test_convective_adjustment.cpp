// C/C++
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// scm
#include <single_column/single_column.hpp>

// canoe
#include <air_parcel.hpp>

class TestConvectiveAdjustment : public testing::Test {
 protected:
  Mesh* pmesh;
  ParameterInput* pinput;
  Real Ps, Ts, grav;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_convective_adjustment.inp", IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    Thermodynamics::InitFromAthenaInput(pinput);

    // set up mesh
    int restart = false;
    int mesh_only = false;
    pmesh = new Mesh(pinput, mesh_only);

    // set up components
    for (int b = 0; b < pmesh->nblocal; ++b) {
      MeshBlock* pmb = pmesh->my_blocks(b);
      pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
    }

    Ps = pinput->GetReal("problem", "Ps");
    Ts = pinput->GetReal("problem", "Ts");
    grav = -pinput->GetReal("hydro", "grav_acc1");
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    Thermodynamics::Destroy();
    IndexMap::Destroy();
    delete pinput;
    delete pmesh;
  }
};

TEST_F(TestConvectiveAdjustment, RandomProfile) {
  auto pmb = pmesh->my_blocks(0);
  auto pcoord = pmb->pcoord;

  // prepare random real number
  std::random_device rd;
  std::mt19937 gen(rd());  // Mersenne Twister generator

  // Define the distribution to be uniform for real numbers
  // between 0.0 and 1.0
  std::uniform_real_distribution<> dis(-20., 20.0);

  // output file
  std::ofstream outFile1("ac_before.csv");
  outFile1 << "i,T_fluct,x1,pres,temp,theta" << std::endl;

  // construt the air column
  std::vector<AirParcel> vector_ac;

  AirParcel air(AirParcel::Type::MoleFrac);
  air.SetZero();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;

  auto pthermo = Thermodynamics::GetInstance();

  // half a grid to cell center
  pthermo->Extrapolate(&air, pcoord->dx1f(pmb->is) / 2.,
                       Thermodynamics::Method::ReversibleAdiabat, grav);

  int js = pmb->js, ks = pmb->ks;
  for (int i = pmb->is; i <= pmb->ie; ++i) {
    Real T_fluct = dis(gen);
    air.w[IDN] += T_fluct;

    AirParcelHelper::distribute_to_primitive(pmb, ks, js, i, air);
    AirParcelHelper::distribute_to_conserved(pmb, ks, js, i, air);
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::PseudoAdiabat, grav);
    Real theta = pthermo->PotentialTemp(pmb, Ps, ks, js, i);

    outFile1 << i << "," << T_fluct << "," << pcoord->x1v(i) << ","
             << air.w[IPR] << "," << air.w[IDN] << "," << theta << std::endl;
  }
  outFile1.close();

  pmb->pimpl->pscm->ConvectiveAdjustment(pmb, pmb->ks, pmb->js);

  // output air column after convective adjustment
  std::ofstream outFile2("ac_after.csv");
  outFile2 << "i,x1,pres,temp,theta" << std::endl;
  for (int i = pmb->is; i <= pmb->ie; ++i) {
    air = AirParcelHelper::gather_from_conserved(pmb, ks, js, i);
    AirParcelHelper::distribute_to_primitive(pmb, ks, js, i, air);

    Real theta = pthermo->PotentialTemp(pmb, Ps, ks, js, i);
    air.ToMoleFraction();
    outFile2 << i << "," << pcoord->x1v(i) << "," << air.w[IPR] << ","
             << air.w[IDN] << "," << theta << std::endl;
  }
  outFile2.close();
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
