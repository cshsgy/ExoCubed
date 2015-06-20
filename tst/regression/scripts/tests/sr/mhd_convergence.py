# Test script for relativistic MHD linear wave convergence

# Modules
import numpy as np
import math
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison

# Parameters
wave_flags = range(7)
res_low = 64
res_high = 512
cutoff = 1.8
amp = 1.0e-6
gamma_adi = 4.0/3.0
rho = 4.0
pgas = 1.0
vx = 0.1
vy = 0.3
vz = -0.05
bx = 2.5
by = 1.8
bz = -1.2

# Prepare Athena++
def prepare():
  athena.configure('sb',
      prob='linear_wave_rel',
      coord='cartesian',
      flux='hlld')
  athena.make()

# Run Athena++
def run():
  wavespeeds = wavespeeds_mhd()
  for wave_flag in wave_flags:
    time = 1.0 / abs(wavespeeds[wave_flag])
    arguments = [
        'job/problem_id=sr_mhd_wave_{0}_low'.format(wave_flag),
        'mesh/nx1=' + repr(res_low),
        'time/tlim=' + repr(time),
        'output1/dt=' + repr(time),
        'fluid/gamma=' + repr(gamma_adi),
        'problem/rho=' + repr(rho), 'problem/pgas=' + repr(pgas),
        'problem/vx=' + repr(vx), 'problem/vy=' + repr(vy), 'problem/vz=' + repr(vz),
        'problem/Bx=' + repr(bx), 'problem/By=' + repr(by), 'problem/Bz=' + repr(bz),
        'problem/wave_flag=' + repr(wave_flag), 'problem/amp=' + repr(amp)]
    athena.run('mhd_sr/athinput.linear_wave', arguments)
    arguments[0] = 'job/problem_id=sr_mhd_wave_{0}_high'.format(wave_flag)
    arguments[1] = 'mesh/nx1=' + repr(res_high)
    athena.run('mhd_sr/athinput.linear_wave', arguments)

# Analyze outputs
def analyze():

  # Specify tab file columns
  headings = ('x', 'rho', 'pgas', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz')

  # Check that convergence is attained for each wave
  for wave_flag in wave_flags:

    # Read low and high resolution initial and final states
    prim_initial_low = athena.read_tab(
        'bin/sr_mhd_wave_{0}_low.block0.out1.00000.tab'.format(wave_flag),
        headings=headings, dimensions=1)
    prim_initial_high = athena.read_tab(
        'bin/sr_mhd_wave_{0}_high.block0.out1.00000.tab'.format(wave_flag),
        headings=headings, dimensions=1)
    prim_final_low = athena.read_tab(
        'bin/sr_mhd_wave_{0}_low.block0.out1.00001.tab'.format(wave_flag),
        headings=headings, dimensions=1)
    prim_final_high = athena.read_tab(
        'bin/sr_mhd_wave_{0}_high.block0.out1.00001.tab'.format(wave_flag),
        headings=headings, dimensions=1)

    # Calculate overall errors for low and high resolution runs
    epsilons_low = []
    epsilons_high = []
    for quantity in headings[1:]:
      qi = prim_initial_low[quantity][0,0,:]
      qf = prim_final_low[quantity][0,0,:]
      epsilons_low.append(math.fsum(abs(qf-qi)) / res_low)
      qi = prim_initial_high[quantity][0,0,:]
      qf = prim_final_high[quantity][0,0,:]
      epsilons_high.append(math.fsum(abs(qf-qi)) / res_high)
    epsilons_low = np.array(epsilons_low)
    epsilons_high = np.array(epsilons_high)
    epsilon_low = (math.fsum(epsilons_low**2) / len(epsilons_low))**0.5 / amp
    epsilon_high = (math.fsum(epsilons_high**2) / len(epsilons_high))**0.5 / amp

    # Test fails if convergence is not at least that specified by cutoff
    if epsilon_high / epsilon_low > (float(res_low) / float(res_high))**cutoff:
      return False

  # All waves must have converged
  return True

# MHD wavespeed calculator
def wavespeeds_mhd():

  # Handle simple entropy case
  wavespeeds = np.empty(7)
  wavespeeds[3] = vx

  # Calculate 4-vectors
  v_sq = vx**2 + vy**2 + vz**2
  u = np.empty(4)
  u[0] = 1.0 / (1.0 - v_sq)**0.5
  u[1] = u[0]*vx
  u[2] = u[0]*vy
  u[3] = u[0]*vz
  b = np.empty(4)
  b[0] = bx*u[1] + by*u[2] + bz*u[3]
  b[1] = 1.0/u[0] * (bx + b[0]*u[1])
  b[2] = 1.0/u[0] * (by + b[0]*u[2])
  b[3] = 1.0/u[0] * (bz + b[0]*u[3])

  # Calculate useful scalars
  gamma_adi_red = gamma_adi / (gamma_adi - 1.0)
  b_sq = -b[0]**2 + sum(b[1:]**2)
  wgas = rho + gamma_adi_red * pgas
  wtot = wgas + b_sq
  cs_sq = gamma_adi * pgas / wgas

  # Calculate Alfven speeds
  lambda_ap = (b[1] + wtot**0.5 * u[1]) / (b[0] + wtot**0.5 * u[0])
  lambda_am = (b[1] - wtot**0.5 * u[1]) / (b[0] - wtot**0.5 * u[0])
  wavespeeds[1] = min(lambda_ap, lambda_am)
  wavespeeds[5] = max(lambda_ap, lambda_am)

  # Calculate magnetosonic speeds
  factor_a = wgas * (1.0/cs_sq - 1.0)
  factor_b = -(wgas + b_sq/cs_sq)
  a4 = factor_a * u[0]**4 - factor_b * u[0]**2 - b[0]**2
  a3 = -factor_a * 4.0 * u[0]**4 * vx \
      + factor_b * 2.0 * u[0]**2 * vx + 2.0 * b[0] * b[1]
  a2 = factor_a * 6.0 * u[0]**4 * vx**2 \
      + factor_b * u[0]**2 * (1.0-vx**2) + b[0]**2 - b[1]**2
  a1 = -factor_a * 4.0 * u[0]**4 * vx**3 \
      - factor_b * 2.0 * u[0]**2 * vx - 2.0 * b[0] * b[1]
  a0 = factor_a * u[0]**4 * vx**4 + factor_b * u[0]**2 * vx**2 + b[1]**2
  roots = np.roots([a4, a3, a2, a1, a0])
  roots.sort()
  wavespeeds[0] = roots[0]
  wavespeeds[2] = roots[1]
  wavespeeds[4] = roots[2]
  wavespeeds[6] = roots[3]
  return wavespeeds