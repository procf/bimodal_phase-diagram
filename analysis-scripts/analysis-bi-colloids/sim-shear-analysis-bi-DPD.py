## Analyze the results of a DPD shear simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 solvent type (typeid=0)
##                           and 1 colloid type (typeid=1)
## NOTE: to select specific analyses, scroll to the bottom
##       of this file and comment out unwanted analyses in the
##       RUN CHECKS ON A SIMULATION section
"""
## This code performs the following shear analyses: 
##   - extracts temperature and pressure (the negative of stress) from all frames of a GSD file 
##   - calculates:
##        * corrected (kinetic) pressure and temperature 
##           - tracking only the effect of colloid particles, not solvents
##        * solvent velocity profile
##        * colloid velocity profile
##        * non-affinity (amount of affine vs. non-affine motion)
##        * the fabric tensor
##
## For other analyses, use the standard analysis module and analysis scripts
"""
## (Dr. Deepak Mangal, Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import pandas as pd
import gsd.hoomd
import math
import fortranmod as module
import os
import re
import glob
import sys

##########################
""" INPUT PARAMETERS """
##########################

from_path = True

if from_path == True:
  # call variables from path
  import config_from_path
  if len(sys.argv) > 3:
    data_path_name = sys.argv[1]
    gsd_sim_path = sys.argv[2]
    shear_sim_path = sys.argv[3]
    shear_filepath, gel_filepath, data_outpath, shear_style, shear_rate_const, shear_rate, Amp, omega, kT, eta0, frames_per_strain, dt_Integration, theta, Iter_cycle, period, MDtime, colloid1_typeid, colloid2_typeid = config_from_path.call_shear_variables(data_path_name, gsd_sim_path, shear_sim_path)
  # if you provide data_path_name and gel_sim_path and run this where shear_sim_path == '..'
  elif len(sys.argv) > 2:
    data_path_name = sys.argv[1]
    gsd_sim_path = sys.argv[2] 
    shear_filepath, gel_filepath, data_outpath, shear_style, shear_rate_const, shear_rate, Amp, omega, kT, eta0, frames_per_strain, dt_Integration, theta, Iter_cycle, period, MDtime, colloid1_typeid, colloid2_typeid = config_from_path.call_shear_variables(data_path_name, gsd_sim_path)
  # if you only provide data_path_name and run this where shear_sim_path == '..' and gel_sim_path == '../..'
  elif len(sys.argv) > 1:
    data_path_name = sys.argv[1]
    shear_filepath, gel_filepath, data_outpath, shear_style, shear_rate_const, shear_rate, Amp, omega, kT, eta0, frames_per_strain, dt_Integration, theta, Iter_cycle, period, MDtime, colloid1_typeid, colloid2_typeid = config_from_path.call_shear_variables(data_path_name) 
  # if you use all the defaults
  else:
    shear_filepath, gel_filepath, data_outpath, shear_style, shear_rate_const, shear_rate, Amp, omega, kT, eta0, frames_per_strain, dt_Integration, theta, Iter_cycle, period, MDtime, colloid1_typeid, colloid2_typeid = config_from_path.call_shear_variables() 

elif from_path == False:
  # manual simulation parameters
  shear_filepath = 'Shear.gsd'
  gel_filepath = 'Gelation.gsd'
  data_outpath = 'data'
  shear_style = 'constant' 
  shear_rate_const = 0.5
  Amp = 0.0
  omega = 0.0
  kT = 0.1
  eta0 = 0.3
  frames_per_strain = 500 
  dt_Integration = 1e-4 
  theta = 1.0
  if shear_style == 'const':
    shear_rate = shear_rate_const
    #vinf = shear_rate * L_Y
    t_ramp = theta/shear_rate
  else:
    omega_timesteps = omega*dt_Integration # frequency (1/N_timesteps)
    #vinf = Amp * omega * L_Y
    shear_rate = Amp * omega
    t_ramp = Amp/shear_rate
  Iter_cycle = int(t_ramp / dt_Integration)
  period= int(Iter_cycle/frames_per_strain)
  MDtime = period*dt_Integration
  colloid1_typeid = 1
  colloid2_typeid = 2


## Get nstrains from the GSD file
traj = gsd.hoomd.open(shear_filepath, 'r')
nframes = len(traj)
nstrains = (nframes-1) / frames_per_strain

## For the noaffinity calculation, be ready to compare to pre-shear gel
affine_ref_filepath = gel_filepath 

# record simulation parameters
print("\nSIMULATION TIME PARAMETERS")
print(" * shear rate:",shear_rate)
print(" * Iter_cycle (timesteps per strain):", Iter_cycle) 
print(" * period (timesteps per frame):", period)
print(" * MDtime per frame (MDtime):",MDtime) 
print(" * total number of frames:",nframes)
print(" * frames per strain:",frames_per_strain)
print(" * total number of simulated strains:",nstrains,"\n")






###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)


#######
## extract thermodynamic properties from GSD for all frames
"""
# get temperature and pressure (AKA negative of stress), including all virial components
"""
def extract_properties_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)

  # create a file to save the thermodynamic properties
  f=open(data_outpath+'/gsd-properties.txt', 'w')
  f.write('simframe Virial-Pressure Vr_CONS Vr_DISS Vr_RAND Vr_SQUE Vr_CONT PE kT tps\n')

  # for each frame
  for i in range(0, nframes-1):
    simframe = i
    # extract the total "force virial contribution" of the pressure tensor
    Pr=float(traj[i].log['md/compute/ThermodynamicQuantities/pressure_tensor'][1])
    # extract the decomposed virial compoenents of the pressure tensor
    Vr=np.asarray(traj[i].log['md/compute/ThermodynamicQuantities/virial_ind_tensor'],dtype=float)
    # extract the potential energy and scale to kilo-units (1/1000)
    Pe=float(traj[i].log['md/compute/ThermodynamicQuantities/potential_energy'][0])*0.001
    # extract the kinetic temperature (kT)
    kT=float(traj[i].log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0])
    # extract the transactios per secont (tps) for assessing speed/efficiency of the simulation
    tps=int(traj[i].log['Simulation/tps'][0])

    # write these values to the output file:
    # raw values
    f.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.format(simframe, Pr, 
      Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, kT, tps))

    # rounded values
    #f.write('%f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %d %0.2f %d\n'%((i+1)*t1, Pr, 
    #  Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, kT, tps))
  print("Data extracted from GSD file")
#######


#######
## corrected (kinetic) pressure and temperature
"""
## 1. extract ONLY the effect of colloid motion on pressure and temperature
## 2. normalize pressure by simulation volume

#  1. - divide the simulation into layers/bins along the y-axis (-L_Y/2 to L_Y/2)
#     - in each bin:
#         - find the average SOLVENT velocity in the x-direction
#         - subtract this bin/average solvent velocity from the x-velocity of EVERY 
#             particle in the bin (vel_x is now relative to solvent velocity)  
#         - use this corrected velocity to correct pressure (pressure is now ONLY
#             the effect of colloid particles)
#     - use the corrected pressure to calculate kinetic temperature (kT)
#       (normalized by dimension (3) and number of particles)
#  2. - normalize the pressure by the simulation box volume
"""

# set the parameters for binning the system
nlayers_TP = 50 # n_layers in y direction (should be at least particle size)

def corrected_temperature_py(filename):
  traj = gsd.hoomd.open(filename, 'r')
  nframe=len(traj)
  outname = "T_corrected.txt"
  f=open(data_outpath+'/'+outname,'w')
  f.write('frame P_xx P_xy P_xz kT\n')
  for i in range(0,nframe):
    nparticles=traj[i].particles.N
    typeid=traj[i].particles.typeid
    Lbox=traj[i].configuration.box[0:3]
    vel=traj[i].particles.velocity
    pos=traj[i].particles.position
    mass=traj[i].particles.mass
    pxx=0.0;
    pxy=0.0;
    pxz=0.0;
    kT=0.0
    pxx,pxy,pxz,kT=module.corrected_temperature(nparticles,nlayers_TP,pos,vel,mass,typeid,Lbox)
    f.write("%f,%f,%f,%f,%f\n"%(i,pxx,pxy,pxz,kT))

  f.close()
#######

#######
# solvent particle velocity profile
"""
# bin the velocities to create a velocity profile within your flow geometry
# NOTE: this calculates the SOLVENT velocity profile
"""

# set the parameters for binning the system
nlayers_VP = 100 # n_layers in y direction (should be at least particle size)

# choose which frames to analyze (defaults to last frame)
framechoice = [-1] # ex: [-1], [10,11,12], ['all']

def solvent_vel_profile_py(filename):
  traj = gsd.hoomd.open(filename, 'r')
  nframes = len(traj)
  solvents = np.where(traj[-1].particles.typeid == [0])[0]
  nsolvents = len(solvents)
  Lbox=traj[-1].configuration.box[0:3]

  if framechoice[0] == 'all':
    framechoice[0] = len(traj)

  df_list = []
  for frame in framechoice:
    if frame == -1:
      frame = nframes-1
    frame_df = pd.DataFrame()
    vel_s=traj[frame].particles.velocity[solvents]
    pos_s=traj[frame].particles.position[solvents]
    vel_bin,bin_array = module.vel_profile(nsolvents,nlayers_VP,pos_s,vel_s,Lbox)
    frame_df['y_pos'] = bin_array
    frame_df['velocity'] = vel_bin
    frame_df.insert(loc=0, column='frame', value=frame)
    df_list.append(frame_df)
  velprofile_df = pd.concat(df_list, ignore_index=True) 

  velprofile_df.to_csv(data_outpath+'/vel-profile-solvent.csv', index=False) 
  print("Velocity profile calculated")

#######


#######
# colloid particle velocity profile
"""
# bin the velocities to create a velocity profile within your flow geometry
# NOTE: this calculates the COLLOID velocity profile
"""

# set the parameters for binning the system
nlayers_VP = 100 # n_layers in y direction (should be at least particle size)

# choose which frames to analyze (defaults to last frame)
framechoice = [-1] # ex: [-1], [10,11,12], ['all']

def colloid_vel_profile_py(filename):
  traj = gsd.hoomd.open(filename, 'r')
  nframes = len(traj)
  #colloids = np.where(traj[-1].particles.typeid == [colloid_typeid])[0]
  colloids = np.where((traj[-1].particles.typeid == [colloid1_typeid]) | (traj[-1].particles.typeid == [colloid2_typeid]))[0]
  ncolloids = len(colloids)
  Lbox=traj[-1].configuration.box[0:3]

  if framechoice[0] == 'all':
    framechoice[0] = len(traj)

  df_list = []
  for frame in framechoice:
    if frame == -1:
      frame = nframes-1
    frame_df = pd.DataFrame()
    vel=traj[frame].particles.velocity[colloids]
    pos=traj[frame].particles.position[colloids]
    vel_bin,bin_array = module.vel_profile(ncolloids,nlayers_VP,pos,vel,Lbox)
    frame_df['y_pos'] = bin_array
    frame_df['velocity'] = vel_bin
    frame_df.insert(loc=0, column='frame', value=frame)
    df_list.append(frame_df)
  velprofile_df = pd.concat(df_list, ignore_index=True) 

  velprofile_df.to_csv(data_outpath+'/vel-profile-colloids.csv', index=False) 
  print("Velocity profile calculated")

#######


#######
# non-affinity
"""
# calculate if the motion/transformation is affine or non-affine
# by comparing the gelation with the last frame
"""

affin_frames = 100 # the number of frames to check movement in

def non_affinity_py(filename, MDtime):
  # get the gel trajectory
  traj_gel = gsd.hoomd.open(affine_ref_filepath, 'r')
  Lbox=traj_gel[-1].configuration.box[:3]
  inv_Lbox = 1.0/Lbox
  index=[]
  #index.extend(np.where(traj_gel[-1].particles.typeid == [colloid_typeid])[0])
  index.extend(np.where((traj_gel[-1].particles.typeid == [colloid1_typeid]) | (traj_gel[-1].particles.typeid == [colloid2_typeid]))[0])
  init_pos=traj_gel[-1].particles.position[index]
  pos_prev_unwrap=init_pos.copy()

  # get the shear trajectory
  traj = gsd.hoomd.open(filename, 'r')
  nframe=len(traj)

  outname='non_affinity.txt'
  f=open(data_outpath+'/'+outname,'w')
  f.write("MDtime nonaff\n")

  if nframe < affin_frames:
    frames = nframe
  else:
    frames = affin_frames

  # find how each particle moved
  for i in range(frames):
    if(i==0):
      dr = traj[i].particles.position[index,:]-pos_prev_unwrap
    else:
      dr = traj[i].particles.position[index,:]-traj[i-1].particles.position[index,:]
    dr[:,2] -= Lbox[2]*np.rint(dr[:,2]*inv_Lbox[2])
    img = Lbox[1] * np.rint(dr[:,1]*inv_Lbox[1])
    dr[:,1] -= img
    dr[:,0] -= img*traj[i].configuration.box[3]
    dr[:,0] -= Lbox[0]*np.rint(dr[:,0]*inv_Lbox[0])
    pos_prev_unwrap += dr
    pos = pos_prev_unwrap - init_pos
    pos[:,0] -= i*MDtime*init_pos[:,1]

    # find if the motion is nonaffine
    nonaff = np.mean(np.linalg.norm(pos,axis=1))
    f.write('%f,%f\n'%(i*MDtime,nonaff))
  f.close()
  print("Calculated the degree of non-affinity in the particle motion")
#######


#######
# fabric tensor
"""
# calculate the fabric tensor for all frames
# see this paper for more info:
#     Jamali, McKinley, Armstrong (2017) PRL
#     "Microstructural Rearrangements and their Rheological Implications in a Model Thixotropic Elastoviscoplastic Fluid"
#     DOI: https://doi.org/10.1103/PhysRevLett.118.048003
"""
def fabric_tensor_py(filename):
  cut_off=0.3
  traj = gsd.hoomd.open(filename, 'r')
  nframe=len(traj)
  Lbox=traj[-1].configuration.box[:3]
  inv_Lbox = 1.0/Lbox
  index=[]
  #index.extend(np.where(traj[-1].particles.typeid == [colloid_typeid])[0])
  index.extend(np.where((traj[-1].particles.typeid == [colloid1_typeid]) | (traj[-1].particles.typeid == [colloid2_typeid]))[0])
  nparticles=int(len(index))
  radius=0.5*traj[-1].particles.diameter[index]
  outname='fabric_tensor.txt'
  f=open(data_outpath+'/'+outname,'w')
  f.write("frame Rxy_avg\n")
  for i in range(0,nframe):
    pos=traj[i].particles.position[index]
    m_xy=traj[i].configuration.box[3]
    Rxy_avg=module.fabric_tensor(nparticles,pos,radius,Lbox,cut_off,m_xy)
    f.write("%f %f \n"%(i,Rxy_avg))
  f.close()
  print("Calculated the fabric tensor")
#######

if __name__ == '__main__':
  extract_properties_py(shear_filepath)
  if 'DPD' in shear_filepath:
    corrected_temperature_py(shear_filepath)
    solvent_vel_profile_py(shear_filepath)
  else:
    print('Not a DPD sim, corrected temp and solvent velocity profile NOT calculated')
  colloid_vel_profile_py(shear_filepath)
  non_affinity_py(shear_filepath,MDtime)
  fabric_tensor_py(shear_filepath)

