// C++ headers
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <climath/core.h>  // sqr

// athena
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// Parameters
Real h0, g, omega_planet, radius;
Real sponge_lat, sponge_tau;
Real vrad, vphi;
Real vis, b; 
Real max_cyclone, max_anticyclone;

bool thermal_tide;
/* Thickening effect of the fluid layer due to the movement of substellar point. 
 * 'tide_amp' is the amplitude of the tide thickness, and 'tau_t' is the relaxa-
 * -tion time of layer cooling.
*/
Real obliq, tide_amp, orbital_period, rotation_period, tau_t;

bool istopowave;
/* IS there some topographical waves
 * The topographical wave forcing is set to a sinusoidal function wrapped by a 
 * Gaussian envelope.
 * 'wave_amp' is the amplitude of topographical waves.
 * 'wave_pos_x' is the center position of wave packet in longitudal degree.
 * 'wave_wid_x' is the width of wave packet in longitudal degree.
 * Venus' plateaus are at 60-80N, 50W-10E.
 * 'time_interval' is the wave-generating time interval, which is 7 hr in Navarro & Schubert 2024
*/
Real wave_amp,wave_wid_x,wave_pos_x,wave_wid_y,wave_pos_y;
Real wave_wid_x_rad,wave_pos_x_rad,wave_wid_y_rad,wave_pos_y_rad;
Real phi_add, tau_w, time_interval;

Real wave_wid_x2,wave_pos_x2,wave_wid_y2,wave_pos_y2;
Real wave_wid_x_rad2,wave_pos_x_rad2,wave_wid_y_rad2,wave_pos_y_rad2;

bool iscircuwave;
/* IS there a circumpolar wave forcing
 * The wave forcing belt is 'wavec_wid_y' width in north-south in degree
 * and at 'wavec_pos_y' latitude in degree
*/
Real wavec_amp,wavec_wid_y,wavec_pos_y;
Real wavec_wid_y_rad,wavec_pos_y_rad;

bool iscircurrent;
/* IS there a circumpolar current
 * with 'curr_wid_y' width at 'curr_pos_y' latitude in degree
*/
Real curr_amp,curr_wid_y,curr_pos_y;
Real curr_wid_y_rad,curr_pos_y_rad;

Real sgn(Real x){
if (x==0){
        return 0.0;
} else if (x<0){
        return -1.0;
} else if (x>0){
        return 1.0;
} else{
        return 0;
}   
}   
    
Real abs_real(Real x){
if (x==0){
        return 0.0;
} else if (x<0){
        return -x;
} else if (x>0){
        return x;
} else{
        return 0;
}
}

Real pow_real(Real x, Real y){
return sgn(x)*pow(abs_real(x),y);
}

Real get_mu(Real time, Real lat, Real lon){
    Real omega_orbital = 2.*M_PI/orbital_period;
    Real lonss_deg = std::fmod(180.*(omega_orbital-omega_planet)*time/M_PI,360.)-180.;
    Real latss_deg = obliq*sin(2.*M_PI*time/orbital_period);
    Real lonss = lonss_deg*M_PI/180.;
    Real latss = latss_deg*M_PI/180.;
    return cos(lon)*cos(lat)*cos(lonss)*cos(latss)
       +sin(lon)*cos(lat)*sin(lonss)*cos(latss)
       +sin(lat)*sin(latss);
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  int is = pmb->is;
  Real omega = 2. * PI / (24. * 3600.);
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real cF2, cF3;
        pexo3->CalculateCoriolisForce2(j, k, w(IVY, k, j, i), w(IVZ, k, j, i),
                                       omega_planet, w(IDN, k, j, i), &cF2, &cF3);
        u(IVY, k, j, i) += dt * cF2;
        u(IVZ, k, j, i) += dt * cF3;
      }
    }
}

void Forcing2(MeshBlock *pmb, Real const time, Real const dt,
              AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
              AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
              AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        // coriolis force
        Real f = 2. * omega_planet * sin(lat);
        Real U, V; //zonal and meridional velocity in lat-lon grid

	pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);
	
	Real dist = radius * (PI/2. - lat);

	Real s = (90. - dist/radius/M_PI*180.)/sponge_lat;
        Real uprim, vprim, ucon, vcon;
	if (s < 1. && s > 0.) {  // sponge layer
            uprim = dt*w(IDN,k,j,i)*w(IM2,k,j,i)*pow((1.-s), 2)/sponge_tau;
            vprim = dt*w(IDN,k,j,i)*w(IM3,k,j,i)*pow((1.-s), 2)/sponge_tau;
            pexo3->ContravariantVectorToCovariant(j, k, uprim, vprim, &ucon, &vcon);
	    u(IM2,k,j,i) -= ucon;
            u(IM3,k,j,i) -= vcon;		    
        } else {  // viscosity
            u(IM2,k,j,i) -= 0.0;
	    u(IM3,k,j,i) -= 0.0;
        }
	Real ll_acc_U = f * V; //Coriolis acceleration in lat-lon grid
        Real ll_acc_V = -f * U;
        Real acc1, acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i); //convert to cubed sphere
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        u(IM2, k, j, i) += dt * w(IDN, k, j, i) * acc2;
        u(IM3, k, j, i) += dt * w(IDN, k, j, i) * acc3;
      }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Application::Logger app("main");
  app->Log("ProblemGenerator: Cubed Sphere");

  auto pexo3 = pimpl->pexo3;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN, k, j, i) = g * h0;}

       // initial polar vortex
      /*
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, ie);
   
        Real fcor = 2. * omega_planet * sin(lat);
        Real dist = radius * (PI/2. - lat);
	Real vsize = 1.0*vrad;

	Real gam = boost::math::tgamma(2./b, (1./b)*pow(dist/vsize, b));
	Real phi = -fcor*vphi*vsize*exp(1./b)*pow(b, 2./b - 1.)*gam;
        
	Real vel = (vphi/vsize)*dist*exp((1./b)*(1-pow((dist/vsize), b)));
	Real U, V;
        U = vel; V = 0.0;
        Real Vy, Vz;
        pexo3->GetVyVz(&Vy, &Vz, U, V, k, j, ie);

	phydro->w(IDN, k, j, ie) += phi;
        phydro->w(IVY, k, j, ie) = Vy;
        phydro->w(IVZ, k, j, ie) = Vz;*/
      }
  app->Log("Done with problem generator");
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // Calculate lat, lon, U, V for output and visualization
  for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        Real x = tan(pcoord->x2v(j));
        Real y = tan(pcoord->x3v(k));
        Real C = sqrt(1.0 + x * x);
        Real D = sqrt(1.0 + y * y);
        Real delta = 1.0 / (1.0 + x * x + y * y);

        Real Vy = phydro->w(IVY, k, j, i);
        Real Vz = phydro->w(IVZ, k, j, i);
        Real lat, lon;
        pimpl->pexo3->GetLatLon(&lat, &lon, k, j, i);
        user_out_var(0, k, j, i) = lat / PI * 180.0;
        user_out_var(1, k, j, i) = lon / PI * 180.0;
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1 * x1 + x2 * x2);
        Real f = 0.0;
        Real U, V;
        pimpl->pexo3->GetUV(&U, &V, Vy, Vz, k, j, i);
        user_out_var(2, k, j, i) = U;
        user_out_var(3, k, j, i) = V;
        user_out_var(4, k, j, i) = delta / (C * D);
	user_out_var(5,k,j,i) = phydro->w(IDN,k,j,i);
      }
}

void MeshBlock::UserWorkInLoop()        //called at the end of every time step
{
    Real time = pmy_mesh->time;
    Real dt = pmy_mesh->dt;

   if (thermal_tide){
   for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j) {
        Real lat, lon, phi_eq, heating;
        pimpl->pexo3->GetLatLon(&lat, &lon, k, j, ie);
           phi_eq = h0 * g + tide_amp
               *abs_real(pow_real(get_mu(time, lat, lon),1.5));
           heating = (phi_eq - phydro->w(IDN,k,j,ie))/tau_t;
        phydro->w(IDN,k,j,ie) += dt * heating;
    }
    }

   if (istopowave){
   wave_wid_x_rad = wave_wid_x/180.0*M_PI;wave_pos_x_rad = wave_pos_x/180.0*M_PI;
   wave_wid_y_rad = wave_wid_y/180.0*M_PI;wave_pos_y_rad = wave_pos_y/180.0*M_PI;
   wave_wid_x_rad2 = wave_wid_x2/180.0*M_PI;wave_pos_x_rad2 = wave_pos_x2/180.0*M_PI;
   wave_wid_y_rad2 = wave_wid_y2/180.0*M_PI;wave_pos_y_rad2 = wave_pos_y2/180.0*M_PI;

   for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j) {
       Real lat, lon, wave_sine;
       Real wave_sine2;
       Real envelope_x,envelope_y,topo_wave;
       Real envelope_x2,envelope_y2,topo_wave2;
       pimpl->pexo3->GetLatLon(&lat, &lon, k, j, ie);
       wave_sine = sin(M_PI*(lon-wave_pos_x_rad)/4/wave_wid_x_rad);
       envelope_x = exp(-abs_real(lon-wave_pos_x_rad)/wave_wid_x_rad);
       envelope_y = exp(-abs_real(lat-wave_pos_y_rad)/wave_wid_y_rad);
       topo_wave = wave_amp*wave_sine*envelope_x*envelope_y;

       wave_sine2 = sin(M_PI*(lon-wave_pos_x_rad2)/4/wave_wid_x_rad2);
       envelope_x2 = exp(-abs_real(lon-wave_pos_x_rad2)/wave_wid_x_rad2);
       envelope_y2 = exp(-abs_real(lat-wave_pos_y_rad2)/wave_wid_y_rad2);
       topo_wave2 = wave_amp*wave_sine2*envelope_x2*envelope_y2;

       phi_add = h0 * g + topo_wave + topo_wave2;
       if (std::fmod(time,time_interval)<400.0) {
	  phydro->w(IDN,k,j,ie) = phydro->w(IDN,k,j,ie) + phi_add;}
    }
   }

   if (iscircuwave){
   wavec_pos_y_rad=wavec_pos_y/180.0*M_PI;
   wavec_wid_y_rad=wavec_wid_y/180.0*M_PI;
   for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j) {
       Real lat, lon, envelopec_y, circu_wave;
       pimpl->pexo3->GetLatLon(&lat, &lon, k, j, ie);
       envelopec_y = exp(-abs_real(lat-wavec_pos_y_rad)/wavec_wid_y_rad);
       circu_wave = wavec_amp*envelopec_y;
       phi_add = h0 * g + circu_wave;
              if (std::fmod(time,time_interval)<400.0) {
          phydro->w(IDN,k,j,ie) = phydro->w(IDN,k,j,ie) + phi_add;}
   }
   }

   if (iscircurrent){
     curr_pos_y_rad=curr_pos_y/180.0*M_PI;
     curr_wid_y_rad=curr_wid_y/180.0*M_PI;
     for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
       for (int j = js - NGHOST; j <= je + NGHOST; ++j) {
          Real lat, lon;
          pimpl->pexo3->GetLatLon(&lat, &lon, k, j, ie);

	  Real add_U = curr_amp*exp(-abs_real(lat-curr_pos_y_rad)/curr_wid_y_rad); 
          Real add_V = 0.0;

	  Real addy, addz;
          pimpl->pexo3->GetVyVz(&addy, &addz, add_U, add_V, k, j, ie); //convert to cubed sphere

	  if (std::fmod(time,time_interval)<400.0) {
             phydro->w(IVY, k, j, ie) += addy;
             phydro->w(IVZ, k, j, ie) += addz;}
     }
   }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
  SetUserOutputVariableName(4, "sqrtg");
  SetUserOutputVariableName(5, "phi");
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
    h0 = pin->GetReal("problem", "h0");
    g = pin->GetReal("problem", "g");
    omega_planet = pin->GetReal("problem", "omega_planet");
    radius = pin->GetReal("problem", "radius");
    sponge_lat = pin->GetReal("problem", "sponge_lat");
    sponge_tau = pin->GetReal("problem", "sponge_tau");
    vrad = pin->GetReal("problem", "vrad");
    vphi = pin->GetReal("problem", "vphi");
    vis = pin->GetOrAddReal("problem", "vis", 0.);
    max_cyclone = pin->GetOrAddReal("problem", "max_cyclone", 90.);
    max_anticyclone = pin->GetOrAddReal("problem", "max_anticyclone", 90.);
    b = pin->GetOrAddReal("problem", "b", 2.);

    // thermal tide
    thermal_tide = pin->GetBoolean("thermal", "thermal_tide");
    obliq = pin->GetOrAddReal("thermal", "obliq", 3.39458);
    tide_amp = pin->GetOrAddReal("thermal", "tide_amp", 1.E4);
    orbital_period = pin->GetOrAddReal("thermal", "orbital_period", 224.701*86400.0);
    rotation_period = pin->GetOrAddReal("thermal", "rotation_period", 243.0*86400.0);
    tau_t = pin->GetReal("thermal", "tau_t");
 
    // topography waves
    istopowave = pin->GetBoolean("topo", "istopowave");
    wave_amp = pin->GetReal("topo", "wave_amp");
    wave_wid_x = pin->GetReal("topo", "wave_wid_x");
    wave_pos_x = pin->GetReal("topo", "wave_pos_x");
    wave_wid_y = pin->GetReal("topo", "wave_wid_y");
    wave_pos_y = pin->GetReal("topo", "wave_pos_y");
    tau_w = pin->GetReal("topo", "tau_w");
    time_interval = pin->GetReal("topo","time_interval");

    wave_wid_x2 = pin->GetReal("topo", "wave_wid_x2");
    wave_pos_x2 = pin->GetReal("topo", "wave_pos_x2");
    wave_wid_y2 = pin->GetReal("topo", "wave_wid_y2");
    wave_pos_y2 = pin->GetReal("topo", "wave_pos_y2");

    // circu waves
    iscircuwave = pin->GetBoolean("circu", "iscircuwave");
    wavec_wid_y = pin->GetReal("circu", "wavec_wid_y");
    wavec_pos_y = pin->GetReal("circu", "wavec_pos_y");
    wavec_amp = pin->GetReal("circu", "wavec_amp");

    // circu current
    iscircurrent = pin->GetBoolean("circurr", "iscircurrent");
    curr_wid_y = pin->GetReal("circurr", "curr_wid_y");
    curr_pos_y = pin->GetReal("circurr", "curr_pos_y");
    curr_amp = pin->GetReal("circurr", "curr_amp");
  EnrollUserExplicitSourceFunction(Forcing2);
}
