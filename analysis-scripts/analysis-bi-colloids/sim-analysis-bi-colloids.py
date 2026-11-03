## Analyze the results of a BD colloid simulation
## NOTE: requires matching Fortran module and solvopt module
## NOTE: this code assumes 2 colloid types (typeid=0, typeid=1)
## NOTE: to select specific analyses, scroll to the bottom
##       of this file and comment out unwanted analyses in the
##       RUN CHECKS ON A SIMULATION section
"""
## This code performs the following analyses: 
##    - extracts temperature and pressure (the negative of stress) from all frames of a GSD file 
##   - calculates:
##        * colloid coordination number:
##           - coordination number distribution (Z counts) for each frame
##           - average coordination number (<Z>) for each frame
##        * mean squared displacement (MSD) of colloids for each frame
##        * radial distribution function g(r):
##           - colloid-colloid distribution in the final frame
##        * qualitative approximation of the pair correlation function h(r1,r2) 
##           for the colloid-colloid distribution of the final frame OR an 
##           average of the last few frames
##        * static structure factor S(q) (and associated g(r)) 
##           for colloids in the final frame OR a selection of frames
##        * number density fluctuation (NDF) of colloids for each frame 
##          (AKA "variance to mean ratio" or "index of dispersion") 
##          with the option to test a range of bin-sizes
##        * pore size distribution for the final frame OR a selection of frames, 
##          using two methods: 
##           - Torquato’s Pore Size Distribution
##           - Gubbins’s Pore Size Distribution 
##           (requires the solvopt algorithm, included as a separate f90 module)
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import pandas as pd
import gsd.hoomd
import math
import fortranmod as module
from fortranmod import pore_size_calculation
import os
import re
import glob
import sys

##########################
""" INPUT PARAMETERS """
##########################

# call variables from path, or set manually
from_path = True

if from_path == True:
  # call variables from path
  import config_from_path
  if len(sys.argv) > 2:
    data_path_name = sys.argv[1]
    gsd_sim_path = sys.argv[2]
    filepath, data_outpath, period, kT, eta0, Lbox_shortest, R_C1, R_C2, R_min, R_max, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_gel_variables(data_path_name, gsd_sim_path)
  elif len(sys.argv) > 1:
    data_path_name = sys.argv[1]
    filepath, data_outpath, period, kT, eta0, Lbox_shortest, R_C1, R_C2, R_min, R_max, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_gel_variables(data_path_name)
  else:
    filepath, data_outpath, period, kT, eta0, Lbox_shortest, R_C1, R_C2, R_min, R_max, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_gel_variables()

elif from_path == False:
  # manual simulation parameters
  filepath = 'Gelation.gsd' 
  data_outpath = 'data'
  period = 10000
  kT = 0.1
  eta0 = 0.3
  L_X = 30
  Lbox_shortest = L_X # the shortest side of the simulation box
  R_C1 = 1
  R_C2 = 2
  colloid1_typeid = 1
  colloid2_typeid = 2
  kappa = 30

# set min and max particle size
R_min = min(R_C1, R_C2)
R_max = max(R_C1, R_C2)

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
  log = traj[-1].log

  # Check if any log exists
  if not log:
      print('WARNING: No log found, skipping extract_properties_py')
      return

  # Check required keys
  required_keys = [
      'md/compute/ThermodynamicQuantities/pressure_tensor',
      'md/compute/ThermodynamicQuantities/virial_ind_tensor',
      'md/compute/ThermodynamicQuantities/potential_energy',
      'md/compute/ThermodynamicQuantities/kinetic_temperature',
      'Simulation/tps'
  ]

  missing = [k for k in required_keys if k not in log]

  if missing:
      print('WARNING: Missing log keys:', missing)
      print('Skipping extract_properties_py')
      return

  else:
    # create a file to save the thermodynamic properties
    f=open(data_outpath+'/gsd-properties.txt', 'w')
    f.write('simframe Virial-Pressure Vr_CONS Vr_DISS Vr_RAND Vr_SQUE Vr_CONT PE kT tps\n')

    # for each frame
    for i in range(0, nframes):
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
## coordination number AKA contact number 
"""
# calculate the coordination number distribution (Z counts) and 
# the average coordination number (<Z>) for each frame
#
# for a given particle, Z is the number of other particles touching it
# (we define "contact" from the "attraction range" set with kappa; 
#   for a gel the average Z should plateau as the network is formed, 
#    and the distribution is usually centered around an average of Z=6)
"""

# use kappa to calculate attraction range (cut-off distance)
cut_off = round(3/kappa,2)

def coordination_number_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filepath,'r')
  # get the number of frames
  nframes = len(traj)
  # use the last frame to get simulation box size [L_X, L_Y, L_Z]
  Lbox = traj[-1].configuration.box[:3]
  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  # use last frame to get the radius of each colloid
  R_C = 0.5*traj[-1].particles.diameter[colloids]
  # use last frame to get the typeid of each colloid
  typeid = traj[-1].particles.typeid[colloids]

  # generate list of all types of interaction pairs from the data
  ntypes = int(len(np.unique(typeid)))
  pairs = []
  pair_names = []
  for colloid_type_a in range(ntypes):
    for colloid_type_b in range(ntypes):
      pairs.append([colloid_type_a+1, colloid_type_b+1])
  for pair in pairs:
    pair_names.append('Z_'+str(pair[0])+str(pair[1]))
  npairs = len(pairs)
  # and the number of total combos, Z_any, Z_11, Z_12..., Z_1any, ...
  n_allcombos = int(1+npairs+np.sqrt(npairs))

  ## 1. CALCULATE Z DISTRIBUTION

  # gather data for the whole simulation
  allpos_allframe = np.zeros((nframes,ncolloids,3))
  m_xys = np.zeros(nframes)
  for frame in range(nframes):
    # get all particle positions
    allpos_allframe[frame] = traj[frame].particles.position[colloids]
    # get the xy tilt factor (square=0.0, sheared-right=0.45)
    m_xys[frame] = traj[frame].configuration.box[3]

  print('Calculating Z-distribution data for '+str(len(pairs)+1)+' types of colloid-colloid pairs in '+str(nframes)+' times...')

  # "pair type" labels start at 1 (i.e. Z_11, Z_12, etc), but GSD typeids may not; account for this in Fortran
  if 0 in np.unique(typeid):
    fortran_pairs = []
    for pair in pairs:
      new_pair = []
      for i in pair:
        BD_i = i-1
        new_pair.append(BD_i)
      fortran_pairs.append(new_pair)
  else:
    fortran_pairs = pairs.copy()

  # run the fortran module
  Zs_array = module.coordination_number(nframes,Lbox,ncolloids,typeid,R_C,m_xys,allpos_allframe,cut_off,npairs,n_allcombos,fortran_pairs)

  """
  # convert array to data frame for easy saving
  allframes_Zs_df = []
  for frame in range(nframes):
    frame_df = pd.DataFrame(Zs_array[frame])
    frame_df.reset_index(inplace=True)
    frame_df = frame_df.rename(columns = {'index':'colloidID'})
    frame_df.insert(loc=1, column='typeid', value=typeid)
    frame_df.insert(loc=0, column='frame', value=frame)
    frame_df = frame_df.rename(columns = {0:'Z_any'})
    for pair_index in range(len(pairs)):
      name = pair_names[pair_index]
      frame_df = frame_df.rename(columns = {pair_index+1:name})
    # and the Z_1any, Z_2any, etc...
    for type_index in range(ntypes):
        colloid_type = type_index+1
        name = 'Z_'+str(colloid_type)+'any'
        frame_df = frame_df.rename(columns = {1+npairs+type_index:name})
    allframes_Zs_df.append(frame_df)

  total_Zsdf = pd.concat(allframes_Zs_df, ignore_index=True)

  total_Zsdf.to_csv(data_outpath+'/Z-counts.csv', index=False)
  """
  # vectorize
  nframes, ncolloids, _ = Zs_array.shape

  z_labels = ['Z_any'] + pair_names + [f'Z_{t+1}any' for t in range(ntypes)]

  flat = Zs_array.reshape(nframes * ncolloids, -1)

  total_Zsdf = pd.DataFrame(flat, columns=z_labels)

  total_Zsdf.insert(0, 'colloidID', np.tile(np.arange(ncolloids), nframes))
  total_Zsdf.insert(0, 'frame', np.repeat(np.arange(nframes), ncolloids))
  total_Zsdf.insert(2, 'typeid', np.tile(typeid, nframes))

  total_Zsdf.to_csv(data_outpath + '/Z-counts.csv', index=False)

  print("...Coordination number calculated for all colloids in "+str(nframes)+" frames")

  ## 2. CALCULATE THE AVERAGE Z

  print("Calculating the average coordination number for "+str(len(pairs)+1)+" types of colloid-colloid pairs in "+str(nframes)+" frames...")

  # get the tags and number of colloids for each colloid type
  colloids_subpops = []
  ncolloids_subpops = []
  for colloid_type in np.unique(typeid):
    tags = np.where(traj[-1].particles.typeid == [colloid_type])[0]
    colloids_subpops.append(tags)
    number = len(tags)
    ncolloids_subpops.append(number)

  if sum(ncolloids_subpops) != ncolloids:
    print("ERROR: total number of colloids in <Z> calculation does not equal sum of subpopulations; <Z> was NOT calculated")
    return

  # create an array to save the averages
  Zavgs_array = np.zeros((nframes,n_allcombos))
  #Zavgs_array = np.zeros((nframes,len(pairs)+1))

  #caluclate the average coordination numbers
  for frame in range(nframes):
    any_sum = sum(Zs_array[frame,:,0])
    any_nc = ncolloids
    Zavgs_array[frame][0] = any_sum / any_nc

    for pair_index in range(len(pairs)):
      sub_sum = sum(Zs_array[frame,:,pair_index+1])
      primary_colloid_typeid = pairs[pair_index][0]
      sub_nc = ncolloids_subpops[primary_colloid_typeid-1]
      Zavgs_array[frame][pair_index+1] = sub_sum/sub_nc

    # and the Z_1any, Z_2any, etc...
    for type_index in range(ntypes):
      subany_sum = sum(Zs_array[frame,:,1+npairs+type_index])
      subany_nc = ncolloids_subpops[type_index]
      Zavgs_array[frame][1+npairs+type_index] = subany_sum/subany_nc 

  # convert array to data frame for easy saving
  Zavgs_df = pd.DataFrame(Zavgs_array)
  #Zavgs_df.insert(loc=0,column='simframe', value=list(range(0,nframes)))

  #for frame in range(nframes):
  #  Zavgs_df = Zavgs_df.rename(columns = {0:'Z_any'})
  #  for pair_index in range(len(pairs)):
  #    name = pair_names[pair_index]
  #    Zavgs_df = Zavgs_df.rename(columns = {pair_index+1:name})
  # and the Z_1any, Z_2any, etc...
  #for type_index in range(ntypes):
  #    colloid_type = type_index+1
  #    name = 'Z_'+str(colloid_type)+'any'
  #    Zavgs_df = Zavgs_df.rename(columns = {1+npairs+type_index:name})

  z_labels = (
      ['Z_any'] +
      pair_names +
      [f'Z_{t+1}any' for t in range(ntypes)]
  )
  Zavgs_df = pd.DataFrame(Zavgs_array, columns=z_labels)
  Zavgs_df.insert(0, 'simframe', np.arange(nframes))


  Zavgs_df.to_csv(data_outpath+'/Zavg.csv', index=False)

  print("...Average coordination number calculated for "+str(nframes)+" frames")
#######


#######
## mean squared displacement (MSD)
"""
# calculate the mean squared displacement (MSD) for colloids in all frames,
# as well as the sample standard deviation of MSD. 
#
# MSD measures how far a particle has moved from it's original position over a period of time 
# (if the data forms a diagonal line with a fixed slope, this indicates the particle is 
#   moving steadily (fluid); but if the data forms a flat/horizontal line, the particle 
#   has stopped moving (solid) -- MSD tells you about particle dynamics) 
"""

# calculate the time it takes the smallest (fastest) particle to diffuse half the shortest box length
D = kT/(6*math.pi*eta0*R_min) # diffusion coefficient = r^2 / tau
tau_to_half = (Lbox_shortest/2)**2 / D # diffusion time to half-box (L/2)

def msd_py(filename):

  if period >= tau_to_half:
    print('ERROR: the GSD file\'s recording timestep is too large for accurate MSD calculations.\n\n' +
          'To create data that you can use to accurately calculate the MSD, you should set the'
          ' period/trigger to a small enough value that a particle will not passively move'
          ' a distance equal to 1/2 the shortest box length in between frames.\n\n' + 
          'Rerun the simulation with a smaller trigger/period before calculating MSD.')
  else:
    # open the simulation GSD file as "traj" (trajectory)
    traj = gsd.hoomd.open(filename, 'r')
    # get the simulation box size [L_X, L_Y, L_Z] from the last frame 
    Lbox = traj[-1].configuration.box[:3]
    inv_Lbox = 1.0/Lbox
    # get the number of frames
    nframes = len(traj)	
    # get all the colloid particles in the last frame, and use this to count the total number of colloids
    colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
    ncolloids = len(colloids)
    # do the same for all the colloid subpopulations
    colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
    ncolloid1 = len(colloid1)
    colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
    ncolloid2 = len(colloid2)
    if (ncolloid1 + ncolloid2) != ncolloids:
      print("ERROR: ncolloid1 + ncolloid2 != ncolloids in MSD calc; MSD was NOT calculated")
      return

    # create an empty array for xyz positon of all colloids in all frames    
    allcollpos = np.zeros((nframes,ncolloids,3))
    # get the initial colloid positions from the first frame
    allcollpos[0,:,:] = traj[0].particles.position[colloids]
    # correct the change in position for colloids crossing a box boundary 
    for i in range(1,nframes):
      # calculate the change in position since the previous frame
      delpos = traj[i].particles.position[colloids] - traj[i-1].particles.position[colloids]
      # if it is more than one box length (i.e. the particle crossed a 
      # boundary) correct the position to wrap inside the sim box
      # (matching our periodic boundaries)
      delpos -= Lbox*np.rint(delpos*inv_Lbox)
      # update and record the position of the colloids in this frame
      allcollpos[i,:,:] = allcollpos[i-1,:,:] + delpos

    # create an empty array for xyz positon of all type 1 colloids in all frames    
    allcollpos1 = np.zeros((nframes,ncolloid1,3))
    # get the initial type 1 colloid positions from the first frame
    allcollpos1[0,:,:] = traj[0].particles.position[colloid1]
    # correct the change in position for type 1 colloids crossing a box boundary 
    for i in range(1,nframes):
      # calculate the change in position since the previous frame
      delpos1 = traj[i].particles.position[colloid1] - traj[i-1].particles.position[colloid1]
      # if it is more than one box length (i.e. the particle crossed a 
      # boundary) correct value the position to be inside the sim box
      delpos1 -= Lbox*np.rint(delpos1*inv_Lbox)
      # update and record the position of the colloids in this frame
      allcollpos1[i,:,:] = allcollpos1[i-1,:,:] + delpos1

    # create an empty array for xyz positon of all type 2 colloids in all frames    
    allcollpos2 = np.zeros((nframes,ncolloid2,3))
    # get the initial type 2 colloid positions from the first frame
    allcollpos2[0,:,:] = traj[0].particles.position[colloid2]
    # correct the change in position for colloids crossing a box boundary 
    for i in range(1,nframes):
      # calculate the change in position since the previous frame
      delpos2 = traj[i].particles.position[colloid2] - traj[i-1].particles.position[colloid2]
      # if it is more than one box length (i.e. the particle crossed a 
      # boundary) correct value the position to be inside the sim box
      delpos2 -= Lbox*np.rint(delpos2*inv_Lbox)
      # update and record the position of the colloids in this frame
      allcollpos2[i,:,:] = allcollpos2[i-1,:,:] + delpos2

    # run the module to calculate MSD
    # (output file "msd.txt" is created in the Fortran module)      
    module.msd_calculation(data_outpath,nframes,ncolloids,ncolloid1,ncolloid2,allcollpos,allcollpos1,allcollpos2)

    print("MSD calculation complete")
#######


#######
## radial distribution function g(r)  
"""
# calculate the probability density per unit volume, g(r), for the LAST FRAME
# of a simulation. 
# This code calculates g(r) for both colloid-solvent and colloid-colloid 
# interactions; colloid-solvent g(r) is only used to double check simulation physics
#
# g(r) is used to measure how much the structure deviates from complete randomness;
# choose a reference particle, and at a distance r away from its center, count
# how many other particles exist, then calculate the density at that distance and 
# compare it to a random distribution -- this g(r) is independent of direction 
# (r is a SCALAR, not a vector), it is "radial" because we consider all directions
# surounding the particle
#
# g(r) can be calculated from densities:
# particle_density_at_r = bulk_density * g(r)
#
# (in a gelled system we expect g(r) to show one large peak at the distance equal to 
#   particle-particle contact (AKA the particle diameter); this is the most likely 
#   position to find another particle, and it is followed by smaller peaks that decay 
#   to the ideal-gas limit: (1 − 1/N) ≈ 1)
"""

# set the parameters binning the system
gofr_rmax = 6.0*R_max # recommended to use rmax > half Lbox_shortest
if gofr_rmax > Lbox_shortest/2:
  print("WARNING: gofr_rmax is too large, the calculation will double-count across periodic boundaries")
gofr_bin_width = 0.1 # the side-length of each bin

def gofr_py(filename):	
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the simulation box size [L_X, L_Y, L_Z] from the last frame 
  Lbox = traj[-1].configuration.box[:3]
  # get the xy tilt factor (box deformation in x direction) from the last frame
  m_xy = traj[-1].configuration.box[3]
  # get the total number of particles from the last frame
  nparticle = traj[-1].particles.N
  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  # do the same for all the colloid subpopulations
  colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
  ncolloid1 = len(colloid1)
  colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
  ncolloid2 = len(colloid2)
  if (ncolloid1 + ncolloid2) != ncolloids:
    print("ERROR: ncolloid1 + ncolloid2 != ncolloids in gofr calc; gofr was NOT calculated")
    return
  # get the typeid of all particles from the last frame
  typeid = traj[-1].particles.typeid
  # get the xyz position of all particles from the last frame
  pos=traj[-1].particles.position
  # set the maximum distance considered from the center of a particle
  # (AKA size of the search box)
  rmax = gofr_rmax
  # set the size of each bin/layer/slab
  bin_width = gofr_bin_width
  # use rmax and bin_width to calculate the number of layers
  nlayers = int(round(rmax/bin_width))

  # run the module to calculate g(r)
  # (output file "gofr.txt" is created in the Fortran module)	
  module.gofr_calc(data_outpath,Lbox,m_xy,nparticle,ncolloids,ncolloid1,ncolloid2,typeid,colloid1_typeid,colloid2_typeid,pos,bin_width,nlayers)
  print("g(r) calculation complete")
#######


#######
## pair correlation function h(r)
"""
# approximate the pair correlation function, a variant of the
# radial distribution function that includes VECTOR information
# about position in 3D space, for the LAST FRAME (or, optionally, 
# from a larger dataset created by summing data from the last few frames)
#
# the pair correlation function measures the probability of finding
# particles at both position r1 and position r2; it is calculated from the
# pair distribution function g(r1,r2):
#   h(r1,r2) = g(r1,r2) - 1
# which can be calculated from particle densities:
#   density_in_2D(r1,r2) = density_at_r1 * density_at_r2 * g(r1,r2)
# For a fluid of spherically symmetric particle in bulk, where the 
# external potential is zero, translation invariance simplifies this to
#   density_in_2D(r12) = bulk_density**2 * g(r1,r2)
# 
# for our purposes, this can be approximated by JUST COUNTING particles
# in 2D slices of 3D space and storing these counts in 2D arrays. 
# (for a disordered or weakly ordered system, each 2D plot shows uniform
#   rings at distances that match the peaks in the radial distribution
#   function; but if the system has directional order -- such as in a
#   crystal -- you will see uneven patterns that are not rings)
"""

# set the parameters for binning the system
pcf_bin_size = 0.25 # the side-length of each cubic bin
pcf_max_dims = 10*R_max # recommended to use pcf_max_dims < half Lbox_shortest
if pcf_max_dims > Lbox_shortest/2:
  print("WARNING: pcf_max_dims is too large, the calculation will double-count across periodic boundaries")
pcf_slice_width = 0.5*R_min # the width of each 2D slice

# choose how many frames to sum together from the end of the sim
# ex: '3' chooses the last 3 frames of the sim
frames_to_sum = 1

def pair_correlation_py(filename):

  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)
  # get the simulation box size [L_X, L_Y, L_Z] from the last frame 
  Lbox = traj[-1].configuration.box[:3]
  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)

  # do the same for all the colloid subpopulations
  colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
  ncolloid1 = len(colloid1)
  colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
  ncolloid2 = len(colloid2)
  if (ncolloid1 + ncolloid2) != ncolloids:
    print("ERROR: ncolloid1 + ncolloid2 != ncolloids in PCF calc; PCF was NOT calculated")
    return 

  # get the typeid of all colloids from the last frame
  typeid = traj[-1].particles.typeid[colloids]

  # set the maximum distance considered from the center of a particle
  # (AKA size of the search box)
  bin_max = pcf_max_dims
  # set the size of each bin/layer/slab
  bin_size = pcf_bin_size
  # set the size of each 2D slice
  slice_width = pcf_slice_width
  # setting array dimensions based on bin size and max/min distance
  array_dim = int((2*bin_max)/bin_size)
  # the maximum distance concidered
  dims_max = pcf_max_dims

  # include the option to sum over multiple frames
  start_frame = nframes - frames_to_sum
  end_frame = nframes
  frame_list = list(range(start_frame, end_frame))

  # gather data for the frames to be averaged
  pos_allframe = np.zeros((frames_to_sum,ncolloids,3))
  pos_c1_allframe = np.zeros((frames_to_sum,ncolloid1,3))
  pos_c2_allframe = np.zeros((frames_to_sum,ncolloid2,3))
  m_xys = np.zeros(frames_to_sum)
  for i in range(len(frame_list)):
    frame = frame_list[i] 
    # get all colloid positions in the frame
    pos_allframe[i] = traj[frame].particles.position[colloids]
    pos_c1_allframe[i] = traj[frame].particles.position[colloid1]
    pos_c2_allframe[i] = traj[frame].particles.position[colloid2]
    # get the xy tilt factor (square=0.0, sheared-right=0.45)
    m_xys[i] = traj[frame].configuration.box[3]

  # create arrays to hold the final sums
  counts_XY_CC_sum = np.zeros((array_dim, array_dim))
  counts_XZ_CC_sum = np.zeros((array_dim, array_dim))
  counts_YZ_CC_sum = np.zeros((array_dim, array_dim))

  counts_XY_C1C1_sum = np.zeros((array_dim, array_dim))
  counts_XZ_C1C1_sum = np.zeros((array_dim, array_dim))
  counts_YZ_C1C1_sum = np.zeros((array_dim, array_dim))

  counts_XY_C2C2_sum = np.zeros((array_dim, array_dim))
  counts_XZ_C2C2_sum = np.zeros((array_dim, array_dim))
  counts_YZ_C2C2_sum = np.zeros((array_dim, array_dim))

  counts_XY_C1C2_sum = np.zeros((array_dim, array_dim))
  counts_XZ_C1C2_sum = np.zeros((array_dim, array_dim))
  counts_YZ_C1C2_sum = np.zeros((array_dim, array_dim))

  # loop through each frame
  for i in range(len(frame_list)):
    m_xy = m_xys[i]
    pos_all = pos_allframe[i]
    pos_c1 = pos_c1_allframe[i]
    pos_c2 = pos_c2_allframe[i]

    # run the module and get 2D counts to approximate PCF
    #print("...binning particle loaction for PCF for all colloids")
    counts_XY_CC,counts_XZ_CC,counts_YZ_CC,counts_XY_C1C1,counts_XZ_C1C1,counts_YZ_C1C1,counts_XY_C1C2,counts_XZ_C1C2,counts_YZ_C1C2,counts_XY_C2C2,counts_XZ_C2C2,counts_YZ_C2C2 = module.pcf_calc(Lbox,m_xy,ncolloids,pos_all,bin_size,dims_max,array_dim,slice_width,typeid,colloid1_typeid,colloid2_typeid)

    counts_XY_CC_sum = counts_XY_CC_sum + counts_XY_CC
    counts_XZ_CC_sum = counts_XZ_CC_sum + counts_XZ_CC
    counts_YZ_CC_sum = counts_YZ_CC_sum + counts_YZ_CC

    counts_XY_C1C1_sum = counts_XY_C1C1_sum + counts_XY_C1C1
    counts_XZ_C1C1_sum = counts_XZ_C1C1_sum + counts_XZ_C1C1
    counts_YZ_C1C1_sum = counts_YZ_C1C1_sum + counts_YZ_C1C1

    counts_XY_C1C2_sum = counts_XY_C1C2_sum + counts_XY_C1C2
    counts_XZ_C1C2_sum = counts_XZ_C1C2_sum + counts_XZ_C1C2
    counts_YZ_C1C2_sum = counts_YZ_C1C2_sum + counts_YZ_C1C2

    counts_XY_C2C2_sum = counts_XY_C2C2_sum + counts_XY_C2C2
    counts_XZ_C2C2_sum = counts_XZ_C2C2_sum + counts_XZ_C2C2
    counts_YZ_C2C2_sum = counts_YZ_C2C2_sum + counts_YZ_C2C2


  df_counts_xy_cc = pd.DataFrame(data = counts_XY_CC_sum)
  df_counts_xz_cc = pd.DataFrame(data = counts_XZ_CC_sum)
  df_counts_yz_cc = pd.DataFrame(data = counts_YZ_CC_sum)

  df_counts_xy_cc.to_csv(data_outpath+'/PCF_counts_XY.csv', header=False, index=False)
  df_counts_xz_cc.to_csv(data_outpath+'/PCF_counts_XZ.csv', header=False, index=False)
  df_counts_yz_cc.to_csv(data_outpath+'/PCF_counts_YZ.csv', header=False, index=False)

  df_counts_xy_c1c1 = pd.DataFrame(data = counts_XY_C1C1_sum)
  df_counts_xz_c1c1 = pd.DataFrame(data = counts_XZ_C1C1_sum)
  df_counts_yz_c1c1 = pd.DataFrame(data = counts_YZ_C1C1_sum)

  df_counts_xy_c1c1.to_csv(data_outpath+'/PCF_counts_c1c1_XY.csv', header=False, index=False)
  df_counts_xz_c1c1.to_csv(data_outpath+'/PCF_counts_c1c1_XZ.csv', header=False, index=False)
  df_counts_yz_c1c1.to_csv(data_outpath+'/PCF_counts_c1c1_YZ.csv', header=False, index=False)

  df_counts_xy_c1c2 = pd.DataFrame(data = counts_XY_C1C2_sum)
  df_counts_xz_c1c2 = pd.DataFrame(data = counts_XZ_C1C2_sum)
  df_counts_yz_c1c2 = pd.DataFrame(data = counts_YZ_C1C2_sum)

  df_counts_xy_c1c2.to_csv(data_outpath+'/PCF_counts_c1c2_XY.csv', header=False, index=False)
  df_counts_xz_c1c2.to_csv(data_outpath+'/PCF_counts_c1c2_XZ.csv', header=False, index=False)
  df_counts_yz_c1c2.to_csv(data_outpath+'/PCF_counts_c1c2_YZ.csv', header=False, index=False)

  df_counts_xy_c2c2 = pd.DataFrame(data = counts_XY_C2C2_sum)
  df_counts_xz_c2c2 = pd.DataFrame(data = counts_XZ_C2C2_sum)
  df_counts_yz_c2c2 = pd.DataFrame(data = counts_YZ_C2C2_sum)

  df_counts_xy_c2c2.to_csv(data_outpath+'/PCF_counts_c2c2_XY.csv', header=False, index=False)
  df_counts_xz_c2c2.to_csv(data_outpath+'/PCF_counts_c2c2_XZ.csv', header=False, index=False)
  df_counts_yz_c2c2.to_csv(data_outpath+'/PCF_counts_c2c2_YZ.csv', header=False, index=False)

  print("PCF data collection complete, summed for "+str(frames_to_sum)+" frames")
#######


#######
## static structure factor S(q)
"""
# calculates the static structure factor S(q)
# (and the associated g(r) values used to caluclate S(q))
#
# S(q) gives us information about short-range order in our colloidal structures, 
# it is a mathematical description of how a material scatters incident radiation 
# (i.e. comparable to experimental data from X-ray diffraction or other diffraction methods)
#
# S(q) measures the fluctuation of the **Fourier component** of the fluctuations 
# of pair density in space. The q is a wavenumber or scattering vector associated 
# with the type of lightwaves used in scattering experiments.
#
# We do NOT match specific q values or calculate S(q) from q itself, instead we can
# caluclate it from h(r)=g(r)-1 (the radial correlation function)
# The result approximates S(q) over a range of high and low q values.
#
# In a log-log plot, the slope of S(q) corresponds to lengthscales in the system
# (you will often see S(q) reported as S(q*a), where a is particle size)
#
# In a linear plot, the peaks of S(q) should correspond to different types of order
# (i.e. different crystal structures); however we usually use a log-log plot
#
# Citation: Filipponi 1994, DOI 10.1088/0953-8984/6/41/006
# Citation using S(q) for colloids: Shah et.al. 2003, DOI 10.1021/la020982g
# Citation re: mean nearest-neighbor cluster separation: Huang et.al. 1998, DOI 10.1103/PhysRevE.57.875
# Citation for log-log plots: Cerdà et.al. 2004, DOI 10.1103/PhysRevE.70.051405
#
# (in a gelled system we expect S(q) to have a high peak at LOW q values, and then decay to
#   1 at HIGH q values. The more distinct a peak is at low q, the more homogeneous the structures are.
#   In depletion colloidal gels, we expect to see two slopes: one corresponding with cluster size, and
#   one with cluster separation; additionally, the plateau that intersects the y-axis is a measure of
#   the compressibility of the system (based purely on particle position NOT including bond properties, etc.)
"""

# set the min distance to integrate from (to reduce noise at low-q values) 
# (usually the smallest particle size AKA point of center-center contact AKA expected g(r) peak)
sq_gofr_rmin = 2*R_min-0.01

# calculate S(q) for the last frame only (or select a different frame)
framechoice_sq = [-1]

def sofq_py(filename):	
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)

  # convert negative framechoice values into specific frames 
  for i in range(len(framechoice_sq)):
    if framechoice_sq[i] < 0:
      framechoice_sq[i] = nframes - abs(framechoice_sq[i])
  # replace the total nframes with the desired nframes to analyze
  nframes = len(framechoice_sq)

  # get the simulation box size [L_X, L_Y, L_Z] from the last frame 
  Lbox = traj[-1].configuration.box[:3]
  Lbox_shortest = min(Lbox[0],Lbox[1],Lbox[2])

  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)

  # do the same for all the colloid subpopulations
  colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
  ncolloid1 = len(colloid1)
  colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
  ncolloid2 = len(colloid2)
  if (ncolloid1 + ncolloid2) != ncolloids:
    print("ERROR: ncolloid1 + ncolloid2 != ncolloids in gofr calc; gofr was NOT calculated")
    return

  # get the typeid of all colloids from the last frame
  typeid = traj[-1].particles.typeid[colloids]

  # gather data for the frames to be averaged
  pos_allframe = np.zeros((nframes,ncolloids,3))
  pos_c1_allframe = np.zeros((nframes,ncolloid1,3))
  pos_c2_allframe = np.zeros((nframes,ncolloid2,3))
  m_xys = np.zeros(nframes)
  for i in range(nframes):
    frame = framechoice_sq[i]
    # get all colloid positions in the frame
    pos_allframe[i] = traj[frame].particles.position[colloids]
    pos_c1_allframe[i] = traj[frame].particles.position[colloid1]
    pos_c2_allframe[i] = traj[frame].particles.position[colloid2]
    # get the xy tilt factor (square=0.0, sheared-right=0.45)
    m_xys[i] = traj[frame].configuration.box[3]


  # set the maximum distance considered from the center of a particle
  # (AKA size of the search box, recommended 10 < rmax < Lbox_shortest/2)
  rmax = Lbox_shortest/2.0
  # set the minimum distance considered from the center of a particle
  rmin = sq_gofr_rmin
  #rmin_c1 = 2*R_min-0.01
  #rmin_c2 = 2*R_max-0.01

  # set the q values (likely 0-pi or 0-10, but avoid 0 for log-scale plot)
  qmin = 4*np.pi/Lbox[0]
  qmax = np.pi
  nq = 500
  q = np.logspace(np.log10(qmin),np.log10(qmax),nq)
  #q = q*R_C # scale q -> qa

  # run the module to calculate g(r) and S(q)
  # (output files "gofr_sq.csv" and "sofq.csv" created in the Fortran module)
  #print("...calculating structure factor for all colloids")
  pop_type = "-allc"
  module.structure_factor(data_outpath,pop_type,nframes,framechoice_sq,Lbox,m_xys,ncolloids,ncolloid1,ncolloid2,typeid,colloid1_typeid,colloid2_typeid,pos_allframe,rmax,rmin,nq,q)

  print("S(q), and associated g(r), calculations complete for all partial populations")
#######


#######
## number density fluctuation as a function of bin size
## NOTE: this analysis does not use the Fortran module
"""
# a normalized measure of the dispersion of the number density: it is a measure 
# used to quantify whether a set of observed occurrences (e.g. particles) are 
# clustered or dispersed compared to a standard statistical model (random distribution)
#
# NDF may also be called "mean to variance ratio" or "index of dispersion"
# and is calculated as: (variance / mean) = (<N^2> - <N>^2) / <N>
#
# (We expect the NDF to be zero when the system has uniform, random, unclustered 
#    structure and 1 when the system is heavily ordered/clustered)
#
# [OPTIONAL: explore a range of bin sizes]
# In some situations (such as under flow), the NDF may change with bin-size,
# allowing you to use the bin-size that maximizes NDF as an measure of cluster size
#
# Citation for NDF with colloids: Rajaram and Mohraz 2012, DOI 10.1039/C2SM25936B
"""

# choose the range of bins you want to explore (default 1 size if min==max)
# NOTE: this assumes your simulation is a cube (L_X = L_Y = L_Z = size)
#       for other geometries you will need to assign binwidth in x, y, and z separately
bin_min = 10 # minimum ~ R_C
if bin_min < R_min:
  print("WARNING: NDF bin_min is smaller than the smallest particle size, will produce null data")
bin_max = 10 # maximum ~ the size of the simulation box Lbox_shortest
if bin_max > Lbox_shortest:
  print("WARNING: NDF bin_max is larger than the smallest dimension of the simulation (Lbox_shortest), will produce null data")

def ndfluc_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames 
  nframes = len(traj)
  # get the simulation box size [L_X, L_Y, L_Z] from the last frame 
  Lbox = traj[-1].configuration.box[:3]
  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  # do the same for all the colloid subpopulations
  colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
  ncolloid1 = len(colloid1)
  colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
  ncolloid2 = len(colloid2)
  if (ncolloid1 + ncolloid2) != ncolloids:
    print("ERROR: ncolloid1 + ncolloid2 != ncolloids in gofr calc; gofr was NOT calculated")
    return
  # use last frame to get the typeid of each colloid
  typeid = traj[-1].particles.typeid[colloids]

  # get the xyz position of all all colloids in all frames
  # create an empty array for xyz positon of all colloids in all frames    
  allpos = np.zeros((nframes,ncolloids,3))
  # fill the array        
  for i in range(0,nframes):
    allpos[i,:,:] = traj[i].particles.position[colloids]

  # create the data TXT file
  #f=open('ndfluc-allbins.txt','w')
  f=open(data_outpath+'/ndfluc.txt','w')
  # write first line (labels)
  f.write('frame bin-size ndfluc_C ndfluc_C1 ndfluc_C2\n')  

  print('Calculating number density fluctuation for selected bin-sizes...')

  for i in range(bin_min, bin_max+1):
    binwidth_x = i
    binwidth_y = i
    binwidth_z = i
    binvol = binwidth_x * binwidth_y * binwidth_z

    # use Lbox and bin_width to calculate the number of layers
    nlayers_x = int(np.rint(Lbox[0]/binwidth_x))
    nlayers_y = int(np.rint(Lbox[1]/binwidth_y))
    nlayers_z = int(np.rint(Lbox[2]/binwidth_z))

    # calculate the total number of bins  
    n_bins = nlayers_x * nlayers_y * nlayers_z

    # get the inverse binwidths for faster calculations
    inv_binwidth_x = 1.0/binwidth_x
    inv_binwidth_y = 1.0/binwidth_y
    inv_binwidth_z = 1.0/binwidth_z

    # loop over all frames
    for frame in range(nframes):

      # set up the history bins to track each particle type
      histbin_c = np.zeros((nlayers_x,nlayers_y,nlayers_z))
      histbin_c1 = np.zeros((nlayers_x,nlayers_y,nlayers_z))
      histbin_c2 = np.zeros((nlayers_x,nlayers_y,nlayers_z))

      for colloid in range(ncolloids): 
        # convert position from range(-L/2, L/2) to range(0, L)
        #   and use that to assign each particle to a bin
        x_bin = math.floor((allpos[frame][colloid][0]+0.5*Lbox[0])*inv_binwidth_x)
        y_bin = math.floor((allpos[frame][colloid][1]+0.5*Lbox[1])*inv_binwidth_y)
        z_bin = math.floor((allpos[frame][colloid][2]+0.5*Lbox[2])*inv_binwidth_z)

        # for particles outside the layers, assign them to the last bin
        if x_bin >= nlayers_x:
          x_bin = nlayers_x - 1
        if y_bin >= nlayers_y:
          y_bin = nlayers_y - 1
        if z_bin >= nlayers_z:
          z_bin = nlayers_z - 1

        histbin_c[x_bin][y_bin][z_bin] = histbin_c[x_bin][y_bin][z_bin] + 1

        # and also bin sub-populations separately
        if typeid[colloid] == colloid1_typeid:
          histbin_c1[x_bin][y_bin][z_bin] = histbin_c1[x_bin][y_bin][z_bin] + 1 
        elif typeid[colloid] == colloid2_typeid:
          histbin_c2[x_bin][y_bin][z_bin] = histbin_c2[x_bin][y_bin][z_bin] + 1 


      # check sub-population hist_bins
      if not ((histbin_c1 + histbin_c2) == histbin_c).all(): 
        print("ERROR: the histbin counts of each sub-population do not sum together to equal the total colloid population; NDF was NOT calculated")
        return
      
      # calculate the mean, AKA the average number density per bin, <N>
      sum_histbin_c = np.sum(histbin_c/binvol) # N_C
      mean_c = sum_histbin_c / n_bins # <N_C>
      sum_histbin_c1 = np.sum(histbin_c1/binvol) # N_C1
      mean_c1 = sum_histbin_c1 / n_bins # <N_C1>
      sum_histbin_c2 = np.sum(histbin_c2/binvol) # N_C2
      mean_c2 = sum_histbin_c2 / n_bins # <N_C2>

      # calculate the variance squared (<N^2> - <N>^2)
      sq_histbin_c = np.sum((histbin_c/binvol)**2) # N_C^2
      mean_c_sq = sq_histbin_c / n_bins # <N_C^2>
      variance_c_sq = mean_c_sq - (mean_c*mean_c) # <N_C^2> - <N_C>^2
      sq_histbin_c1 = np.sum((histbin_c1/binvol)**2) # N_C1^2
      mean_c1_sq = sq_histbin_c1 / n_bins # <N_C1^2>
      variance_c1_sq = mean_c1_sq - (mean_c1*mean_c1) # <N_C1^2> - <N_C1>^2
      sq_histbin_c2 = np.sum((histbin_c2/binvol)**2) # N_C2^2
      mean_c2_sq = sq_histbin_c2 / n_bins # <N_C2^2>
      variance_c2_sq = mean_c2_sq - (mean_c2*mean_c2) # <N_C2^2> - <N_C2>^2

      # calculate the number density fluctuation
      NDF_C = variance_c_sq / mean_c
      NDF_C1 = variance_c1_sq / mean_c1
      NDF_C2 = variance_c2_sq / mean_c2

      # - write NDF and current frame as DPD times
      f.write('%f %f %f %f %f\n'%(frame, i, NDF_C, NDF_C1, NDF_C2))

    print("...number density fluctuation calculated for bin dimensions "+str(binwidth_x)+"x"+str(binwidth_y)+"y"+str(binwidth_z)+"z")
  print("Number density fluctuation calculated for all frames and all bin sizes") 
#######


#######
## pore size calculation
"""
# calculates a distribution of the void space (approximated as spheres) in between particle clusters
# 
# There are two methods that are used:
#   - Torquato’s Pore Size Distribution (volume where the center of a particle can fit in between clusters)
#   - Gubbins’s Pore Size Distribution (volume occupied by a whole particlein in between clusters) 
#       (Gubbin's PSD requires solvopt, the Solver For Local Nonlinear Optimization Problems, 
#        incuded as a secondary fortran module) 
# 
# These methods are described in the Section 4.2 and the Appendix of
# Sorichetti, Hugouvieux, and Kob 2020, DOI: 10.1021/acs.macromol.9b02166
#
# This code assumes you are analyzing a porous medium made of uniform particles (like a colloidal gel),
#  and uses particle trajectories in gsd format; It takes particle and probe size as inputs.
#  Then it uses Linked-list method to compute minimum distance of a point in the void space from 
#  nearby porous medium particles. And then uses solvopt non-linear optimization code to compute Gubbin's pore size.
#
# (we usually use Gubbin's PSD and expect a plot of probability vs void diameter to peak at the size of the most common voids)
"""

# use kappa to calculate attraction range (cut-off distance)
pore_cut_off = round(3/kappa,2) 

# calculate poresize for the last frame only (or select a different frame)
framechoice_ps = [-1]

def pore_size_calc_py(filename):
  # set the cut-off distance for interaction
  cut_off = pore_cut_off
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)

  # convert negative framechoice values into specific frames 
  for i in range(len(framechoice_ps)):
    if framechoice_ps[i] < 0:
      framechoice_ps[i] = nframes - abs(framechoice_ps[i])
  # replace the total nframes with the desired nframes to analyze
  nframes = len(framechoice_ps)

  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  # find the radii of every colloid
  radii = 0.5*traj[-1].particles.diameter[colloids]

  # do the same for all the colloid subpopulations
  colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
  ncolloid1 = len(colloid1)
  colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
  ncolloid2 = len(colloid2)
  if (ncolloid1 + ncolloid2) != ncolloids:
    print("ERROR: ncolloid1 + ncolloid2 != ncolloids in gofr calc; gofr was NOT calculated")
    return
  radii_c1 = 0.5*traj[-1].particles.diameter[colloid1]
  radii_c2 = 0.5*traj[-1].particles.diameter[colloid2]


  ## 1. CALCULATE THE VOID SIZE DISTRIBUTION

  print('Calculating void size distribution data for '+str(len(framechoice_ps))+' frame(s)...')

  # set the number of random points used to explore void size
  nprobe = 10000 # can test quickly at 1
  # set the size of the cells used in the linked list
  dcell_init = 5.0

  # create empty arrays for holding data for all frames
  box_length = np.zeros((nframes,3))
  rxi = np.zeros((nframes,ncolloids))
  ryi = np.zeros((nframes,ncolloids))
  rzi = np.zeros((nframes,ncolloids))

  # fill arrays with data from each frame
  for i in range(nframes):
    frame = framechoice_ps[i]
    box_length[i,:]=traj[frame].configuration.box[:3]
    rxi[i,:] = traj[frame].particles.position[colloids,0]
    ryi[i,:] = traj[frame].particles.position[colloids,1]
    rzi[i,:] = traj[frame].particles.position[colloids,2]


  # do the same for all the colloid subpopulations
  rxi_c1 = np.zeros((nframes,ncolloid1))
  ryi_c1 = np.zeros((nframes,ncolloid1))
  rzi_c1 = np.zeros((nframes,ncolloid1)) 
  for i in range(nframes):
    frame = framechoice_ps[i]
    rxi_c1[i,:] = traj[frame].particles.position[colloid1,0]
    ryi_c1[i,:] = traj[frame].particles.position[colloid1,1]
    rzi_c1[i,:] = traj[frame].particles.position[colloid1,2]
  rxi_c2 = np.zeros((nframes,ncolloid2))
  ryi_c2 = np.zeros((nframes,ncolloid2))
  rzi_c2 = np.zeros((nframes,ncolloid2)) 
  for i in range(nframes):
    frame = framechoice_ps[i]
    rxi_c2[i,:] = traj[frame].particles.position[colloid2,0]
    ryi_c2[i,:] = traj[frame].particles.position[colloid2,1]
    rzi_c2[i,:] = traj[frame].particles.position[colloid2,2]


  # calculate the pore_size for all the selected data
  print("...calculating pore size for all colloids")
  pop_type = '-allc'
  pore_size_calculation.pore_size_calc(data_outpath,pop_type,ncolloids,nframes,framechoice_ps,nprobe,radii,dcell_init,rxi,ryi,rzi,box_length)

  pore_size_calculation.deallocate_arrays()
  
  # do the same for all the colloid subpopulations
  print("...calculating pore size for colloid type 1 only")
  pop_type = '-c1only'
  pore_size_calculation.pore_size_calc(data_outpath,pop_type,ncolloid1,nframes,framechoice_ps,nprobe,radii_c1,dcell_init,rxi_c1,ryi_c1,rzi_c1,box_length)

  pore_size_calculation.deallocate_arrays()

  print("...calculating pore size for colloid type 2 only")
  pop_type = '-c2only'
  pore_size_calculation.pore_size_calc(data_outpath,pop_type,ncolloid2,nframes,framechoice_ps,nprobe,radii_c2,dcell_init,rxi_c2,ryi_c2,rzi_c2,box_length)

  pore_size_calculation.deallocate_arrays()

  print("...Void size distribution calculated for "+str(nprobe)+" probe points in each of "+str(nframes)+" frames")


  ## 2. CALCULATE THE AVERAGE VOID SIZE

  print("Calculating the average void size for "+str(len(framechoice_ps))+" frame(s)...")

  allc_filename = data_outpath+'/poresize-allc.csv'
  c1_filename = data_outpath+'/poresize-c1only.csv'
  c2_filename = data_outpath+'/poresize-c2only.csv'

  filenames = [allc_filename, c1_filename, c2_filename]
  filetags = ['allc','c1only','c2only']

  frame_avg_list = []
  for f in range(len(filenames)):
    filepath = filenames[f]
    data_allframes_df = pd.read_csv(filepath)
    filetag = filetags[f]

    npores_list = []
    volumes_list = []
    T_avg_list = []
    G_avg_list = []
    frames_list = []
    for i in range(nframes):
      frame = framechoice_ps[i]

      data_df = data_allframes_df.loc[data_allframes_df['frame'] == frame]

      # calculate system volume and add to the dataframe
      L_X = int(np.ceil(max(data_df['probe_posx']))-np.floor(min(data_df['probe_posx'])))
      L_Y = int(np.ceil(max(data_df['probe_posy']))-np.floor(min(data_df['probe_posy'])))
      L_Z = int(np.ceil(max(data_df['probe_posz']))-np.floor(min(data_df['probe_posz'])))
      volumes_list.append(L_X*L_Y*L_Z)
      #print(f" - System volume: {[L_X, L_Y, L_Z]}")

      for poresize_type in ['porediameter_T','porediameter_G']:
        #print(f' - {poresize_type}')

        all_pore_diameters = data_df[poresize_type].to_numpy()

        #################################
        """ AVG PORE SIZE CALCULATION """
        #################################

        # get average
        pore_diameters = all_pore_diameters[np.isfinite(all_pore_diameters)] # clean any NaNs and infs
        pore_avg = np.mean(pore_diameters)
        if poresize_type == 'porediameter_T':
          T_avg_list.append(pore_avg)
        elif poresize_type == 'porediameter_G':
          G_avg_list.append(pore_avg)
        else:
          print(f"ERROR: sizetype must be 'Torquato' or 'Gubbin', not: {sizetype}")
          exit(1)

        # calc npores for future normalization
        npores_dataset = len(pore_diameters)

      # only add npores once, it's the same for Torquato and Gubbin
      npores_list.append(npores_dataset)

      #print("Pore size average calculated for frame "+str(frame)+'...')
      frames_list.append(frame)

    avgs_df = pd.DataFrame(T_avg_list, columns=['T_avg'])
    avgs_df['G_avg'] = G_avg_list
    avgs_df['npores'] = npores_list
    avgs_df['volume'] = volumes_list
    avgs_df['population'] = filetag
    avgs_df.insert(loc=0,column='frame', value=frames_list)

    frame_avg_list.append(avgs_df)

  all_avgs_df = pd.concat(frame_avg_list, ignore_index=True)
  all_avgs_df.to_csv(data_outpath+'/poresize-avg.csv', index=False)

  print("...Average void size calculated for "+str(len(filenames))+" populations in "+str(len(framechoice_ps))+" frame(s)")

#######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':	
  # record params:
  print("Data Parameters")
  print(" - period:",period)
  print(" - Lbox_shortest:",Lbox_shortest)
  print(" - kT:",kT)
  print(" - R_C1:",R_C1)
  print(" - R_C2:",R_C2)
  print(" - eta0:",eta0)
  print(" - kappa:",kappa)
  print(" - colloid1_typeid:",colloid1_typeid)
  print(" - colloid2_typeid:",colloid2_typeid,'\n')

  # run analyses
  extract_properties_py(filepath)	
  coordination_number_py(filepath)
  msd_py(filepath)
  gofr_py(filepath)
  pair_correlation_py(filepath)
  sofq_py(filepath)
  ndfluc_py(filepath)
  pore_size_calc_py(filepath)

  print("\nAll analyses complete")
