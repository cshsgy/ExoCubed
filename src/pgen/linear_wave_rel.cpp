// Linear wave problem generator for special relativity

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>      // sin(), sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/eos/eos.hpp"            // FluidEqnOfState

// Declarations
static Real quadratic_root(Real a1, Real a0, bool greater_root);
static Real cubic_root_real(Real a2, Real a1, Real a0);
static void quartic_roots(Real a3, Real a2, Real a1, Real a0,
    Real *px1, Real *px2, Real *px3, Real *px4);

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field
//   pin: parameters
// Outputs: (none)
// Notes:
//   initializes linear wave with sinusoidal variation
//     sets both primitive and conserved variables
//   references Anton et al. 2010, ApJS 188 1 (A, MHD)
//              Falle & Komissarov 1996, MNRAS 278 586 (FK, hydro)
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  // Prepare index bounds
  MeshBlock *pmb = pfl->pmy_block;
  int il = pmb->is - NGHOST;
  int iu = pmb->ie + NGHOST;
  int jl = pmb->js;
  int ju = pmb->je;
  if (pmb->block_size.nx2 > 1)
  {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmb->block_size.nx3 > 1)
  {
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Read information regarding desired wave and check input
  int wave_flag = pin->GetInteger("problem", "wave_flag");
  Real amp = pin->GetReal("problem", "amp");
  if (wave_flag < 0 or wave_flag > NWAVE-1)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "wave_flag=" << wave_flag << " must be between 0 and " << NWAVE-1
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Get ratio of specific heats
  Real gamma_adi = pfl->pf_eos->GetGamma();  
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read background state
  Real rho = pin->GetReal("problem", "rho");
  Real pgas = pin->GetReal("problem", "pgas");
  Real vx = pin->GetReal("problem", "vx");
  Real vy = pin->GetReal("problem", "vy");
  Real vz = pin->GetReal("problem", "vz");
  Real bx = 0.0, by = 0.0, bz = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    bx = pin->GetReal("problem", "Bx");
    by = pin->GetReal("problem", "By");
    bz = pin->GetReal("problem", "Bz");
  }

  // Calculate background 4-vectors
  Real v_sq = SQR(vx) + SQR(vy) + SQR(vz);
  Real u[4], b[4];
  u[0] = 1.0 / std::sqrt(1.0 - v_sq);
  u[1] = u[0]*vx;
  u[2] = u[0]*vy;
  u[3] = u[0]*vz;
  b[0] = bx*u[1] + by*u[2] + bz*u[3];
  b[1] = 1.0/u[0] * (bx + b[0]*u[1]);
  b[2] = 1.0/u[0] * (by + b[0]*u[2]);
  b[3] = 1.0/u[0] * (bz + b[0]*u[3]);

  // Calculate useful background scalars
  Real b_sq = -SQR(b[0]) + SQR(b[1]) + SQR(b[2]) + SQR(b[3]);
  Real ptot = pgas + 0.5*b_sq;
  Real wgas = rho + gamma_adi_red * pgas;
  Real wtot = wgas + b_sq;
  Real cs_sq = gamma_adi * pgas / wgas;
  Real cs = std::sqrt(cs_sq);

  // Calculate desired perturbation
  Real delta_rho, delta_pgas;
  Real delta_u[4], delta_b[4], delta_v[4];
  Real lambda;
  if (MAGNETIC_FIELDS_ENABLED)  // MHD
  {
    switch (wave_flag)
    {
      case 3:  // entropy (A 46)
      {
        lambda = vx;
        delta_rho = 1.0;
        delta_pgas = 0.0;
        for (int mu = 0; mu < 4; ++mu)
        {
          delta_u[mu] = 0.0;
          delta_b[mu] = 0.0;
        }
        break;
      }
      case 1: case 5:  // Alfven (A 65)
      {
        // Calculate wavespeed
        Real lambda_ap = (b[1] + std::sqrt(wtot) * u[1])
            / (b[0] + std::sqrt(wtot) * u[0]);            // (A 38)
        Real lambda_am = (b[1] - std::sqrt(wtot) * u[1])
            / (b[0] - std::sqrt(wtot) * u[0]);            // (A 38)
        Real sign = 1.0;
        if (lambda_ap > lambda_am)  // \lambda_{a,\pm} = \lambda_a^\pm
        {
          if (wave_flag == 1)  // leftgoing
            sign = -1.0;
        }
        else  // lambda_{a,\pm} = \lambda_a^\mp
          if (wave_flag == 5)  // rightgoing
            sign = -1.0;
        if (sign > 0)  // want \lambda_{a,+}
          lambda = lambda_ap;
        else  // want \lambda_{a,-} instead
          lambda = lambda_am;

        // Prepare auxiliary quantities
        Real alpha_1[4], alpha_2[4];
        alpha_1[0] = u[3];                                              // (A 58)
        alpha_1[1] = lambda * u[3];                                     // (A 58)
        alpha_1[2] = 0.0;                                               // (A 58)
        alpha_1[3] = u[0] - lambda * u[1];                              // (A 58)
        alpha_2[0] = -u[2];                                             // (A 59)
        alpha_2[1] = -lambda * u[2];                                    // (A 59)
        alpha_2[2] = lambda * u[1] - u[0];                              // (A 59)
        alpha_2[3] = 0.0;                                               // (A 59)
        Real g_1 = 1.0/u[0] * (by + lambda*vy / (1.0-lambda*vx) * bx);  // (A 60)
        Real g_2 = 1.0/u[0] * (bz + lambda*vz / (1.0-lambda*vx) * bx);  // (A 61)
        Real f_1, f_2;
        if (g_1 == 0.0 and g_2 == 0.0)
          f_1 = f_2 = ONE_OVER_SQRT2;  // (A 67)
        else
        {
          f_1 = g_1 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
          f_2 = g_2 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
        }

        // Set perturbation
        delta_rho = 0.0;
        delta_pgas = 0.0;
        for (int mu = 0; mu < 4; ++mu)
        {
          delta_u[mu] = f_1 * alpha_1[mu] + f_2 * alpha_2[mu];
          delta_b[mu] = -sign * std::sqrt(wtot) * delta_u[mu];
        }
        break;
      }
      default:  // magnetosonic (A 71)
      {
        // Calculate wavespeed
        Real factor_a = wgas * (1.0/cs_sq - 1.0);
        Real factor_b = -(wgas + b_sq/cs_sq);
        Real gamma_2 = SQR(u[0]);
        Real gamma_4 = SQR(gamma_2);
        Real coeff_4 = factor_a * gamma_4
                     - factor_b * gamma_2
                     - SQR(b[0]);
        Real coeff_3 = -factor_a * 4.0 * gamma_4 * vx
                     + factor_b * 2.0 * gamma_2 * vx
                     + 2.0 * b[0] * b[1];
        Real coeff_2 = factor_a * 6.0 * gamma_4 * SQR(vx)
                     + factor_b * gamma_2 * (1.0-SQR(vx))
                     + SQR(b[0]) - SQR(b[1]);
        Real coeff_1 = -factor_a * 4.0 * gamma_4 * vx*SQR(vx)
                     - factor_b * 2.0 * gamma_2 * vx
                     - 2.0 * b[0] * b[1];
        Real coeff_0 = factor_a * gamma_4 * SQR(SQR(vx))
                     + factor_b * gamma_2 * SQR(vx)
                     + SQR(b[1]);
        Real lambda_fl, lambda_sl, lambda_sr, lambda_fr;
        quartic_roots(
            coeff_3/coeff_4, coeff_2/coeff_4, coeff_1/coeff_4, coeff_0/coeff_4,
            &lambda_fl, &lambda_sl, &lambda_sr, &lambda_fr);
        Real lambda_other_ms;
        if (wave_flag == 0)
        {
          lambda = lambda_fl;
          lambda_other_ms = lambda_sl;
        }
        if (wave_flag == 2)
        {
          lambda = lambda_sl;
          lambda_other_ms = lambda_fl;
        }
        if (wave_flag == 4)
        {
          lambda = lambda_sr;
          lambda_other_ms = lambda_fr;
        }
        if (wave_flag == 6)
        {
          lambda = lambda_fr;
          lambda_other_ms = lambda_sr;
        }

        // Determine which sign to use
        Real lambda_ap = (b[1] + std::sqrt(wtot) * u[1])
            / (b[0] + std::sqrt(wtot) * u[0]);            // (A 38)
        Real lambda_am = (b[1] - std::sqrt(wtot) * u[1])
            / (b[0] - std::sqrt(wtot) * u[0]);            // (A 38)
        Real lambda_a = lambda_ap;
        Real sign = 1.0;
        if (lambda_ap > lambda_am)  // \lambda_{a,\pm} = \lambda_a^\pm
        {
          if (wave_flag < 3)  // leftgoing
          {
            lambda_a = lambda_am;
            sign = -1.0;
          }
        }
        else  // lambda_{a,\pm} = \lambda_a^\mp
          if (wave_flag > 3)  // rightgoing
          {
            lambda_a = lambda_am;
            sign = -1.0;
          }

        // Prepare auxiliary quantities
        Real a = u[0] * (vx - lambda);                                       // (A 39)
        Real g = 1.0 - SQR(lambda);                                          // (A 41)
        Real b_over_a = -sign * std::sqrt(-factor_b - factor_a * SQR(a)/g);  // (A 68)
        Real alpha_1[4], alpha_2[4];
        alpha_1[0] = u[3];                                                   // (A 58)
        alpha_1[1] = lambda * u[3];                                          // (A 58)
        alpha_1[2] = 0.0;                                                    // (A 58)
        alpha_1[3] = u[0] - lambda * u[1];                                   // (A 58)
        alpha_2[0] = -u[2];                                                  // (A 59)
        alpha_2[1] = -lambda * u[2];                                         // (A 59)
        alpha_2[2] = lambda * u[1] - u[0];                                   // (A 59)
        alpha_2[3] = 0.0;                                                    // (A 59)
        Real alpha_11 = -SQR(alpha_1[0]);
        Real alpha_12 = -alpha_1[0] * alpha_2[0];
        Real alpha_22 = -SQR(alpha_2[0]);
        for (int i = 1; i < 4; ++i)
        {
          alpha_11 += SQR(alpha_1[i]);
          alpha_12 += alpha_1[i] * alpha_2[i];
          alpha_22 += SQR(alpha_2[i]);
        }
        Real g_1 = 1.0/u[0] * (by + lambda*vy / (1.0-lambda*vx) * bx);       // (A 60)
        Real g_2 = 1.0/u[0] * (bz + lambda*vz / (1.0-lambda*vx) * bx);       // (A 61)
        Real c_1 = (g_1*alpha_12 + g_2*alpha_22)
            / (alpha_11*alpha_22 - SQR(alpha_12)) * u[0] * (1.0-lambda*vx);  // (A 63)
        Real c_2 = -(g_1*alpha_11 + g_2*alpha_12)
            / (alpha_11*alpha_22 - SQR(alpha_12)) * u[0] * (1.0-lambda*vx);  // (A 63)
        Real b_t[4];
        for (int mu = 0; mu < 4; ++mu)
          b_t[mu] = c_1 * alpha_1[mu] + c_2 * alpha_2[mu];  // (A 62)
        Real f_1, f_2;
        if (g_1 == 0.0 and g_2 == 0.0)
          f_1 = f_2 = ONE_OVER_SQRT2;  // (A 67)
        else
        {
          f_1 = g_1 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
          f_2 = g_2 / std::sqrt(SQR(g_1) + SQR(g_2));  // (A 66)
        }
        Real phi_plus_a_u[4];
        for (int mu = 0; mu < 4; ++mu)
          phi_plus_a_u[mu] = a * u[mu];
        phi_plus_a_u[0] += lambda;
        phi_plus_a_u[1] += 1.0;

        // Set perturbation
        if (std::abs(lambda-lambda_a)               // using closer magnetosonic wave...
            <= std::abs(lambda_other_ms-lambda_a))  // ...to the associated Alfven wave
        {
          Real b_t_normalized[4];
          Real denom = std::sqrt((alpha_11*alpha_22 - SQR(alpha_12))
              * (SQR(f_1)*alpha_11 + 2.0*f_1*f_2*alpha_12 + SQR(f_2)*alpha_22));
          for (int mu = 0; mu < 4; ++mu)
            b_t_normalized[mu] = ((f_1*alpha_12+f_2*alpha_22) * alpha_1[mu]
                - (f_1*alpha_11+f_2*alpha_12) * alpha_2[mu]) / denom;        // (A 75)
          Real b_t_norm = -SQR(b_t[0]);
          for (int i = 1; i < 4; ++i)
            b_t_norm += SQR(b_t[i]);
          b_t_norm = std::sqrt(b_t_norm);
          denom = SQR(a) - (g+SQR(a)) * cs_sq;
          if (denom == 0.0)
            delta_pgas = 0.0;
          else
            delta_pgas = -(g+SQR(a)) * cs_sq / denom * b_t_norm;  // (A 74)
          delta_rho = rho / (gamma_adi*pgas) * delta_pgas;
          for (int mu = 0; mu < 4; ++mu)
          {
            delta_u[mu] =
                -a*delta_pgas / (wgas*cs_sq*(g+SQR(a))) * phi_plus_a_u[mu]
                - b_over_a / wgas * b_t_normalized[mu];                     // (A 72)
            delta_b[mu] = -b_over_a * delta_pgas/wgas * u[mu]
                - (1.0+SQR(a)/g) * b_t_normalized[mu];                      // (A 73)
          }
        }
        else  // using more distant magnetosonic wave
        {
          delta_pgas = -1.0;                                // (A 78)
          delta_rho = rho / (gamma_adi*pgas) * delta_pgas;
          Real b_t_reduced[4] = {0.0};                      // (A 79)
          Real denom = wgas * SQR(a) - b_sq * g;
          if (denom != 0.0)
            for (int mu = 0; mu < 4; ++mu)
              b_t_reduced[mu] = b_t[mu] / denom;
          for (int mu = 0; mu < 4; ++mu)
          {
            delta_u[mu] = a / (wgas*cs_sq*(g+SQR(a))) * phi_plus_a_u[mu]
                - b_over_a * g/wgas * b_t_reduced[mu];                    // (A 76)
            delta_b[mu] = b_over_a / wgas * u[mu]
                - (1.0+SQR(a)/g) * g * b_t_reduced[mu];                   // (A 77)
          }
        }
      }
    }
  }
  else  // hydro
  {
    // Calculate perturbation in 4-velocity components (Q of FK)
    switch (wave_flag)
    {
      case 1:  // entropy 1/3
      {
        lambda = vx;
        delta_rho = 1.0;
        delta_pgas = 0.0;
        delta_u[1] = delta_u[2] = delta_u[3] = 0.0;
        break;
      }
      case 2:  // entropy 2/3
      {
        lambda = vx;
        delta_rho = 0.0;
        delta_pgas = 0.0;
        delta_u[1] = vx * vy / (1.0 - SQR(vx));
        delta_u[2] = 1.0;
        delta_u[3] = 0.0;
        break;
      }
      case 3:  // entropy 3/3
      {
        lambda = vx;
        delta_rho = 0.0;
        delta_pgas = 0.0;
        delta_u[1] = vx * vz / (1.0 - SQR(vx));
        delta_u[2] = 0.0;
        delta_u[3] = 1.0;
        break;
      }
      default:  // sound
      {
        Real delta = SQR(u[0]) * (1.0-cs_sq) + cs_sq;
        Real v_minus_lambda_a = vx * cs_sq;
        Real v_minus_lambda_b =
            cs * std::sqrt(SQR(u[0]) * (1.0-cs_sq) * (1.0-SQR(vx)) + cs_sq);
        Real v_minus_lambda;
        if (wave_flag == 0)  // leftgoing
          v_minus_lambda = (v_minus_lambda_a + v_minus_lambda_b) / delta;  // (FK A1)
        else  // rightgoing
          v_minus_lambda = (v_minus_lambda_a - v_minus_lambda_b) / delta;  // (FK A1)
        lambda = vx - v_minus_lambda;
        delta_rho = rho;
        delta_pgas = wgas * cs_sq;
        delta_u[1] = -cs_sq * u[1] - cs_sq / u[0] / v_minus_lambda;
        delta_u[2] = -cs_sq * u[2];
        delta_u[3] = -cs_sq * u[3];
      }
    }

    // Calculate perturbation in 3-velocity components (P of FK)
    delta_v[1] = (1.0-SQR(vx)) * delta_u[1] - vx*vy * delta_u[2] - vx*vz * delta_u[3];
    delta_v[2] = -vx*vy * delta_u[1] + (1.0-SQR(vy)) * delta_u[2] - vy*vz * delta_u[3];
    delta_v[3] = -vx*vz * delta_u[1] - vy*vz * delta_u[2] + (1.0-SQR(vz)) * delta_u[3];
    for (int i = 1; i < 4; ++i)
      delta_v[i] /= u[0];
  }

  // Renormalize perturbation to unit L^2 norm
  if (MAGNETIC_FIELDS_ENABLED)  // MHD
  {
    Real perturbation_size = SQR(delta_rho) + SQR(delta_pgas);
    for (int mu = 0; mu < 4; ++mu)
      perturbation_size += SQR(delta_u[mu]) + SQR(delta_b[mu]);
    perturbation_size = std::sqrt(perturbation_size);
    delta_rho /= perturbation_size;
    delta_pgas /= perturbation_size;
    for (int mu = 0; mu < 4; ++mu)
    {
      delta_u[mu] /= perturbation_size;
      delta_b[mu] /= perturbation_size;
    }
  }
  else  // hydro
  {
    Real perturbation_size = SQR(delta_rho) + SQR(delta_pgas);
    for (int i = 1; i < 4; ++i)
      perturbation_size += SQR(delta_v[i]);
    perturbation_size = std::sqrt(perturbation_size);
    delta_rho /= perturbation_size;
    delta_pgas /= perturbation_size;
    for (int i = 1; i < 4; ++i)
      delta_v[i] /= perturbation_size;
  }

  // Calculate wavenumber such that wave has single period over domain
  Real x1_min = pmb->pmy_mesh->mesh_size.x1min;
  Real x1_max = pmb->pmy_mesh->mesh_size.x1max;
  Real x2_min = pmb->pmy_mesh->mesh_size.x2min;
  Real x3_min = pmb->pmy_mesh->mesh_size.x3min;
  Real arg_min, arg_max;
  if (GENERAL_RELATIVITY)
  {
    Real t_left, x_left, y_left, z_left;
    Real t_right, x_right, y_right, z_right;
    pmb->pcoord->MinkowskiCoordinates(0.0, x1_min, x2_min, x3_min,
        &t_left, &x_left, &y_left, &z_left);
    pmb->pcoord->MinkowskiCoordinates(0.0, x1_max, x2_min, x3_min,
        &t_right, &x_right, &y_right, &z_right);
    arg_min = x_left - lambda * t_left;
    arg_max = x_right - lambda * t_right;
  }
  else
  {
    arg_min = x1_min;
    arg_max = x1_max;
  }
  Real wavenumber = 2.0*PI / (arg_max - arg_min);

  // Initialize hydro variables
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i)
      {
        // Find location of cell in spacetime
        Real t, x, y, z;
        if (GENERAL_RELATIVITY)
          pmb->pcoord->MinkowskiCoordinates(0.0, pmb->x1v(i), pmb->x2v(j), pmb->x3v(k),
              &t, &x, &y, &z);
        else
        {
          t = 0.0;
          x = pmb->x1v(i);
        }

        // Calculate scalar perturbations
        Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
        Real rho_local = rho + local_amp * delta_rho;
        Real pgas_local = pgas + local_amp * delta_pgas;

        // Calculate vector perturbations
        Real u_mink[4];
        Real b_mink[4] = {0.0};
        if (MAGNETIC_FIELDS_ENABLED)  // MHD
        {
          for (int mu = 0; mu < 4; ++mu)
          {
            u_mink[mu] = u[mu] + local_amp * delta_u[mu];
            b_mink[mu] = b[mu] + local_amp * delta_b[mu];
          }
        }
        else  // hydro
        {
          Real vx_mink = vx + local_amp * delta_v[1];
          Real vy_mink = vy + local_amp * delta_v[2];
          Real vz_mink = vz + local_amp * delta_v[3];
          u_mink[0] =
              1.0 / std::sqrt(1.0 - SQR(vx_mink) - SQR(vy_mink) - SQR(vz_mink));
          u_mink[1] = u_mink[0] * vx_mink;
          u_mink[2] = u_mink[0] * vy_mink;
          u_mink[3] = u_mink[0] * vz_mink;
        }

        // Transform vector perturbations
        Real u_local[4], b_local[4], u_local_low[4], b_local_low[4];
        if (GENERAL_RELATIVITY)
        {
          pmb->pcoord->TransformVectorCell(u_mink[0], u_mink[1], u_mink[2], u_mink[3],
              k, j, i, &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
          pmb->pcoord->TransformVectorCell(b_mink[0], b_mink[1], b_mink[2], b_mink[3],
              k, j, i, &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
          pmb->pcoord->LowerVectorCell(
              u_local[0], u_local[1], u_local[2], u_local[3], k, j, i,
              &u_local_low[0], &u_local_low[1], &u_local_low[2], &u_local_low[3]);
          pmb->pcoord->LowerVectorCell(
              b_local[0], b_local[1], b_local[2], b_local[3], k, j, i,
              &b_local_low[0], &b_local_low[1], &b_local_low[2], &b_local_low[3]);
        }
        else
          for (int mu = 0; mu < 4; ++mu)
          {
            u_local[mu] = u_mink[mu];
            b_local[mu] = b_mink[mu];
            u_local_low[mu] = (mu == 0 ? -1.0 : 1.0) * u_local[mu];
            b_local_low[mu] = (mu == 0 ? -1.0 : 1.0) * b_local[mu];
          }

        // Calculate useful local scalars
        Real b_sq_local = 0.0;
        for (int mu = 0; mu < 4; ++mu)
          b_sq_local += b_local[mu] * b_local_low[mu];
        Real wtot_local = rho_local + gamma_adi_red * pgas_local + b_sq_local;
        Real ptot_local = pgas_local + 0.5*b_sq_local;

        // Set primitive hydro variables
        pfl->w(IDN,k,j,i) = pfl->w1(IDN,k,j,i) = rho_local;
        pfl->w(IEN,k,j,i) = pfl->w1(IEN,k,j,i) = pgas_local;
        pfl->w(IVX,k,j,i) = pfl->w1(IVX,k,j,i) = u_local[1] / u_local[0];
        pfl->w(IVY,k,j,i) = pfl->w1(IVY,k,j,i) = u_local[2] / u_local[0];
        pfl->w(IVZ,k,j,i) = pfl->w1(IVZ,k,j,i) = u_local[3] / u_local[0];

        // Set conserved hydro variables
        pfl->u(IDN,k,j,i) = u_local[0] * rho_local;
        if (GENERAL_RELATIVITY)
        {
          pfl->u(IEN,k,j,i) = wtot_local*u_local[0]*u_local_low[0]
              - b_local[0]*b_local_low[0] + ptot_local;
          pfl->u(IM1,k,j,i) = wtot_local*u_local[0]*u_local_low[1]
              - b_local[0]*b_local_low[1];
          pfl->u(IM2,k,j,i) = wtot_local*u_local[0]*u_local_low[2]
              - b_local[0]*b_local_low[2];
          pfl->u(IM3,k,j,i) = wtot_local*u_local[0]*u_local_low[3]
              - b_local[0]*b_local_low[3];
        }
        else
        {
          pfl->u(IEN,k,j,i) = wtot_local*u_local[0]*u_local[0] - b_local[0]*b_local[0]
              - ptot_local;
          pfl->u(IM1,k,j,i) = wtot_local*u_local[0]*u_local[1] - b_local[0]*b_local[1];
          pfl->u(IM2,k,j,i) = wtot_local*u_local[0]*u_local[2] - b_local[0]*b_local[2];
          pfl->u(IM3,k,j,i) = wtot_local*u_local[0]*u_local[3] - b_local[0]*b_local[3];
        }
      }

  // Initialize magnetic fields
  if (MAGNETIC_FIELDS_ENABLED)
    for (int k = kl; k <= ku+1; ++k)
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
          if (GENERAL_RELATIVITY)
          {
            // Set B^1 if needed
            if (j != ju+1 and k != ku+1)
            {
              Real t, x, y, z;
              pmb->pcoord->MinkowskiCoordinates(
                  0.0, pmb->x1f(i), pmb->x2v(j), pmb->x3v(k), &t, &x, &y, &z);
              Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
              Real u_mink[4], b_mink[4];
              for (int mu = 0; mu < 4; ++mu)
              {
                u_mink[mu] = u[mu] + local_amp * delta_u[mu];
                b_mink[mu] = b[mu] + local_amp * delta_b[mu];
              }
              Real u_local[4], b_local[4];
              pmb->pcoord->TransformVectorFace1(
                  u_mink[0], u_mink[1], u_mink[2], u_mink[3], k, j, i,
                  &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
              pmb->pcoord->TransformVectorFace1(
                  b_mink[0], b_mink[1], b_mink[2], b_mink[3], k, j, i,
                  &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
              pfd->b.x1f(k,j,i) = b_local[1] * u_local[0] - b_local[0] * u_local[1];
            }

            // Set B^2 if needed
            if (i != iu+1 and k != ku+1)
            {
              Real t, x, y, z;
              pmb->pcoord->MinkowskiCoordinates(
                  0.0, pmb->x1v(i), pmb->x2f(j), pmb->x3v(k), &t, &x, &y, &z);
              Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
              Real u_mink[4], b_mink[4];
              for (int mu = 0; mu < 4; ++mu)
              {
                u_mink[mu] = u[mu] + local_amp * delta_u[mu];
                b_mink[mu] = b[mu] + local_amp * delta_b[mu];
              }
              Real u_local[4], b_local[4];
              pmb->pcoord->TransformVectorFace2(
                  u_mink[0], u_mink[1], u_mink[2], u_mink[3], k, j, i,
                  &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
              pmb->pcoord->TransformVectorFace2(
                  b_mink[0], b_mink[1], b_mink[2], b_mink[3], k, j, i,
                  &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
              pfd->b.x2f(k,j,i) = b_local[2] * u_local[0] - b_local[0] * u_local[2];
            }

            // Set B^3 if needed
            if (i != iu+1 and j != ju+1)
            {
              Real t, x, y, z;
              pmb->pcoord->MinkowskiCoordinates(
                  0.0, pmb->x1v(i), pmb->x2v(j), pmb->x3f(k), &t, &x, &y, &z);
              Real local_amp = amp * std::sin(wavenumber * (x - lambda * t));
              Real u_mink[4], b_mink[4];
              for (int mu = 0; mu < 4; ++mu)
              {
                u_mink[mu] = u[mu] + local_amp * delta_u[mu];
                b_mink[mu] = b[mu] + local_amp * delta_b[mu];
              }
              Real u_local[4], b_local[4];
              pmb->pcoord->TransformVectorFace3(
                  u_mink[0], u_mink[1], u_mink[2], u_mink[3], k, j, i,
                  &u_local[0], &u_local[1], &u_local[2], &u_local[3]);
              pmb->pcoord->TransformVectorFace3(
                  b_mink[0], b_mink[1], b_mink[2], b_mink[3], k, j, i,
                  &b_local[0], &b_local[1], &b_local[2], &b_local[3]);
              pfd->b.x3f(k,j,i) = b_local[3] * u_local[0] - b_local[0] * u_local[3];
            }
          }
          else
          {
            Real local_amp = amp * std::sin(wavenumber * pmb->x1v(i));
            Real u_local[4], b_local[4];
            for (int mu = 0; mu < 4; ++mu)
            {
              u_local[mu] = u[mu] + local_amp * delta_u[mu];
              b_local[mu] = b[mu] + local_amp * delta_b[mu];
            }
            Real by_local = b_local[2]*u_local[0] - b_local[0]*u_local[2];
            Real bz_local = b_local[3]*u_local[0] - b_local[0]*u_local[3];
            if (j != ju+1 and k != ku+1)
              pfd->b.x1f(k,j,i) = bx;
            if (i != iu+1 and k != ku+1)
              pfd->b.x2f(k,j,i) = by_local;
            if (i != iu+1 and j != ju+1)
              pfd->b.x3f(k,j,i) = bz_local;
          }
  return;
}

// Function for finding root of monic quadratic equation
// Inputs:
//   a1: linear coefficient
//   a0: constant coefficient
//   greater_root: flag indicating that larger root is to be returned
//     "larger" does not mean absolute value
// Outputs:
//   returned value: desired root
// Notes:
//   same function as in adiabatic_mhd_sr.cpp, adiabatic_mhd_gr.cpp, and hlld_rel.cpp
//   solves x^2 + a_1 x + a_0 = 0 for x
//   returns abscissa of vertex if there are no real roots
//   follows advice in Numerical Recipes, 3rd ed. (5.6) for avoiding large cancellations
static Real quadratic_root(Real a1, Real a0, bool greater_root)
{
  if (a1*a1 < 4.0*a0)  // no real roots
    return -a1/2.0;
  if (greater_root)
  {
    if (a1 >= 0.0)
      return -2.0*a0 / (a1 + std::sqrt(a1*a1 - 4.0*a0));
    else
      return (-a1 + std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
  }
  else
  {
    if (a1 >= 0.0)
      return (-a1 - std::sqrt(a1*a1 - 4.0*a0)) / 2.0;
    else
      return -2.0*a0 / (a1 - std::sqrt(a1*a1 - 4.0*a0));
  }
}

// Function for finding real root of monic cubic equation
// Inputs:
//   a2: quadratic coefficient
//   a1: linear coefficient
//   a0: constant coefficient
// Outputs:
//   returned value: a real root
// Notes:
//   solves x^3 + a_2 x^2 + a_1 x + a_0 = 0 for x
//   same function as in adiabatic_mhd_sr.cpp and adiabatic_mhd_gr.cpp
//   references Numerical Recipes, 3rd ed. (NR)
static Real cubic_root_real(Real a2, Real a1, Real a0)
{
  Real q = (a2*a2 - 3.0*a1) / 9.0;                       // (NR 5.6.10)
  Real r = (2.0*a2*a2*a2 - 9.0*a1*a2 + 27.0*a0) / 54.0;  // (NR 5.6.10)
  if (r*r - q*q*q < 0.0)
  {
    Real theta = acos(r/std::sqrt(q*q*q));                 // (NR 5.6.11)
    return -2.0 * std::sqrt(q) * cos(theta/3.0) - a2/3.0;  // (NR 5.6.12)
  }
  else
  {
    Real a = -copysign(1.0, r)
        * cbrt(std::abs(r) + std::sqrt(r*r - q*q*q));  // (NR 5.6.15)
    Real b = (a != 0.0) ? q/a : 0.0;                   // (NR 5.6.16)
    return a + b - a2/3.0;
  }
}

// Function for finding extremal real roots of monic quartic equation
// Inputs:
//   a3: cubic coefficient
//   a2: quadratic coefficient
//   a1: linear coefficient
//   a0: constant coefficient
// Outputs:
//   px1: value set to least real root
//   px2: value set to second least real root
//   px3: value set to second greatest real root
//   px4: value set to greatest real root
// Notes:
//   solves x^4 + a3 x^3 + a2 x^2 + a1 x + a0 = 0 for x
//   uses following procedure:
//     1) eliminate cubic term y^4 + b2 y^2 + b1 y + b0
//     2) construct resolvent cubic z^3 + c2 z^2 + c1 z + c0
//     3) find real root z0 of cubic
//     4) construct quadratics:
//          y^2 + d1 y + d0
//          y^2 + e1 y + e0
//     5) find roots of quadratics
//   similar function to those in adiabatic_mhd_sr.cpp and adiabatic_mhd_gr.cpp
static void quartic_roots(Real a3, Real a2, Real a1, Real a0,
    Real *px1, Real *px2, Real *px3, Real *px4)
{
  // Step 1: Find reduced quartic coefficients
  Real b2 = a2 - 3.0/8.0*SQR(a3);
  Real b1 = a1 - 1.0/2.0*a2*a3 + 1.0/8.0*a3*SQR(a3);
  Real b0 = a0 - 1.0/4.0*a1*a3 + 1.0/16.0*a2*SQR(a3) - 3.0/256.0*SQR(SQR(a3));

  // Step 2: Find resolvent cubic coefficients
  Real c2 = -b2;
  Real c1 = -4.0*b0;
  Real c0 = 4.0*b0*b2 - SQR(b1);

  // Step 3: Solve cubic
  Real z0 = cubic_root_real(c2, c1, c0);

  // Step 4: Find quadratic coefficients
  Real d1 = (z0 - b2 > 0.0) ? std::sqrt(z0 - b2) : 0.0;
  Real e1 = -d1;
  Real d0, e0;
  if (b1 < 0)
  {
    d0 = z0/2.0 + std::sqrt(SQR(z0)/4.0 - b0);
    e0 = z0/2.0 - std::sqrt(SQR(z0)/4.0 - b0);
  }
  else
  {
    d0 = z0/2.0 - std::sqrt(SQR(z0)/4.0 - b0);
    e0 = z0/2.0 + std::sqrt(SQR(z0)/4.0 - b0);
  }

  // Step 5: Solve quadratics
  Real y1 = quadratic_root(d1, d0, false);
  Real y2 = quadratic_root(d1, d0, true);
  Real y3 = quadratic_root(e1, e0, false);
  Real y4 = quadratic_root(e1, e0, true);

  // Step 6: Set original quartic roots
  *px1 = std::min(y1, y3) - a3/4.0;
  Real mid_1 = std::max(y1, y3) - a3/4.0;
  *px4 = std::max(y2, y4) - a3/4.0;
  Real mid_2 = std::min(y2, y4) - a3/4.0;
  *px2 = std::min(mid_1, mid_2);
  *px3 = std::max(mid_1, mid_2);
  return;
}