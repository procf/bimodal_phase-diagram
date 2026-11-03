## Extract data to CSV and run Voronoi volume analysis on all colloids
## in a BD simulation
## NOTE: this code assumes 2 colloid types (typeid=0, typeid=1)
##
"""
## This code performs the following analyses:
##    - extracts colloid position data to from GSD to CSV files for all frames
##    - calculates (for all colloids and for each colloid sub-population):
##         * Voronoi volume around each colloid and the "free volume" 
##            (Voronoi volume - particle volume) in the simulation for each frame
##         * statistical measures of the Voronoi volume distribution and free
##            volume distribution for each frame
##            - mean (1st mode)
##            - standard deviation (2nd mode)
##            - skew (3rd mode)
##            - kurtosis and excess kurtosis (4th mode)
##         * bins the distribution of Voronoi volumes and free volume for each frame 
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import gsd.hoomd
import math
import pyvoro
import pandas as pd
from scipy.stats import skew
from scipy.stats import kurtosis
import re
import os
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
    filepath, data_outpath, colloid1_typeid, colloid2_typeid = config_from_path.call_voronoi_variables(data_path_name, gsd_sim_path)
  elif len(sys.argv) > 1:
    data_path_name = sys.argv[1]
    filepath, data_outpath, colloid1_typeid, colloid2_typeid = config_from_path.call_voronoi_variables(data_path_name)
  else:
    filepath, data_outpath, colloid1_typeid, colloid2_typeid = config_from_path.call_voronoi_variables()

elif from_path == False:
  # manual simulation parameters
  filepath = 'Gelation.gsd'
  data_outpath = 'data'
  colloid1_typeid = 1
  colloid2_typeid = 2

# extra distance to make sure particles at the edge of the box are included
# range becomes (min+Lbuffer) to (max+Lbuffer) --> size+(2*Lbuffer) in each direction
Lbuffer = 0.001 

###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)

#######
## posCSV
"""
# create a CSV file with particle (i.e. node) position information for each frame
# (tag, x, y, z, typeID, radius)
"""

def posCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')

  # set path and filename  
  dir_path = data_outpath+'/frame-pos'
  pos_output = dir_path+'/positions_frame' # + <#>.csv in python loop

  # check for existing CSV data 
  if os.path.exists(dir_path) == False:
    os.mkdir(dir_path)
  if os.path.exists(dir_path) == True: 
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'positions_frame\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1

    if nframes == len(traj):  
      print('CSV files already seem to exist for all frames. Not creating new CSV files.')
      return 

  # if the CSV files didn't already exist, create them 
  nframes = len(traj)
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  typeid = traj[-1].particles.typeid[colloids]
  radii = 0.5*traj[-1].particles.diameter[colloids]

  for f in range(nframes):
    rpos = traj[f].particles.position[colloids]

    # make a data frame to export to CSV
    df_pos = pd.DataFrame()
    df_pos['tag'] = colloids
    df_pos['x'] = rpos[:,0]
    df_pos['y'] = rpos[:,1]
    df_pos['z'] = rpos[:,2] 
    df_pos['typeID'] = typeid
    df_pos['radius'] = radii

    df_pos.to_csv(pos_output+str(f)+'.csv', index=False)
  print("Position data saved to CSV for "+str(nframes)+" frames")
#######


#######
## Voronoi volumes
"""
# use pyvoro to calculate the volume surounding each colloid particle
# in each simulation frame:
# - for each particle pair, find the mid-point between them particle
# - use these midpoints to build polyhedra surrounding each particle
# - calculate the volume of each polyhedra and subtract particle volume
"""

def voronoi_volume_py(filename):
  traj = gsd.hoomd.open(filename, 'r')
  Lbox = traj[-1].configuration.box[:3]
  # set individual dimensions
  L_X = Lbox[0]
  L_Y = Lbox[1]
  L_Z = Lbox[2]

  # folder path
  dir_path = data_outpath+'/frame-pos'

  # get length of sim from the CSV data (match existing data not desired data!)
  # NOTE: expects data to be saved as "positions_frame#.csv" inside dir_path
  count = 0
  # Iterate directory
  for path in os.listdir(dir_path):
      # check if current path is a file
      if os.path.isfile(os.path.join(dir_path, path)):
          count += 1
  nframes = count

  # import all data into one dataframe
  all_frame_dfs = []
  for frame in range(nframes):
    #filepath = 'positions_frame99_plusframetag.csv'
    filepath = dir_path+'/positions_frame'+str(frame)+'.csv'
    # import CSV data
    df = pd.read_csv(filepath)
    # rename colums as needed
    df = df.rename(columns={"typeID":"typeid"})
    df.insert(loc=0, column='frame', value=frame)
    all_frame_dfs.append(df)
  allpos_df = pd.concat(all_frame_dfs, ignore_index=True)

  ###### GET SYSTEM PARAMS FROM LAST FRAME

  lastframe_df = allpos_df.loc[allpos_df['frame'] == nframes-1]

  # get the number of colloid types
  ntypes = len(np.unique(lastframe_df['typeid'].to_numpy()))

  # get the particle radii
  R_C = lastframe_df['radius']
  radii = np.unique(R_C)
  a1 = min(radii)
  a2 = max(radii)

  # number of particles
  ncolloids = len(np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0])
  ncolloid1 = len(lastframe_df.loc[lastframe_df['typeid'] == colloid1_typeid])
  ncolloid2 = len(lastframe_df.loc[lastframe_df['typeid'] == colloid2_typeid])
  if ncolloids != ncolloid1 + ncolloid2:
    print("ERROR: total number of colloids does not equal sum of subpopulations in voronoi volumes calc; Voronoi volumes NOT calculated")
    return

  # index numbers for all particles (extracted from first column of the csv file)
  colloid1 = lastframe_df.loc[lastframe_df['typeid'] == colloid1_typeid]["tag"]
  colloid2 = lastframe_df.loc[lastframe_df['typeid'] == colloid2_typeid]["tag"]
  colloids = lastframe_df["tag"] 

  # positions
  allpos1 = allpos_df.loc[allpos_df['typeid'] == colloid1_typeid][['frame','tag','typeid','x', 'y', 'z', 'radius']] # all colloid1 positions
  allpos2 = allpos_df.loc[allpos_df['typeid'] == colloid2_typeid][['frame','tag','typeid','x', 'y', 'z', 'radius']] # all colloid2 positions
  allpos = allpos_df[['frame','tag','typeid','x', 'y', 'z', 'radius']] # all colloid positions

  # create files to save the voronoi volumes
  f_all=open(data_outpath+'/voronoi-volumes.txt', 'w')
  f_all.write('frame particle typeid radius total-volume free-volume\n')

  f_c1=open(data_outpath+'/voronoi-volumes-c1only.txt', 'w')
  f_c1.write('frame particle typeid radius total-volume free-volume\n')

  f_c2=open(data_outpath+'/voronoi-volumes-c2only.txt', 'w')
  f_c2.write('frame particle typeid radius total-volume free-volume\n')

  #vorvols = np.zeros((nframes,ncolloids))

  for frame in range(0,nframes):

    frameradii = allpos.loc[allpos['frame'] == frame]['radius']
    frameradii_array = frameradii.to_numpy()
    frameradii_c1 = allpos1.loc[allpos1['frame'] == frame]['radius']
    frameradii_array_c1 = frameradii_c1.to_numpy()
    frameradii_c2 = allpos2.loc[allpos2['frame'] == frame]['radius']
    frameradii_array_c2 = frameradii_c2.to_numpy()

    frame_particle_volumes = (4/3)*np.pi*(frameradii_array**3)
    frame_particle_volumes_c1 = (4/3)*np.pi*(frameradii_array_c1**3)
    frame_particle_volumes_c2 = (4/3)*np.pi*(frameradii_array_c2**3)

    pos_c = allpos.loc[allpos['frame'] == frame][['x', 'y', 'z']]
    framepos_c = pos_c.to_numpy()
    pos_c1 = allpos1.loc[allpos1['frame'] == frame][['x', 'y', 'z']]
    framepos_c1 = pos_c1.to_numpy()
    pos_c2 = allpos2.loc[allpos2['frame'] == frame][['x', 'y', 'z']]
    framepos_c2 = pos_c2.to_numpy()  
 
    types_all = allpos.loc[allpos['frame'] == frame]['typeid'].to_numpy()
    types_1 = allpos1.loc[allpos['frame'] == frame]['typeid'].to_numpy()
    types_2 = allpos2.loc[allpos['frame'] == frame]['typeid'].to_numpy()

    cells_all = pyvoro.compute_voronoi(
      framepos_c, # point positions
      [[(-L_X/2)-Lbuffer, (L_X/2)+Lbuffer], [(-L_Y/2)-Lbuffer, (L_Y/2)+Lbuffer], [(-L_Z/2)-Lbuffer, (L_Z/2)+Lbuffer]], # limits 
      2.0, # block size
      radii=frameradii_array # matching list of particle radii -- required for polydisperse / radical tessellation 
    )

    cells_c1 = pyvoro.compute_voronoi(
      framepos_c1, # point positions
      [[(-L_X/2)-Lbuffer, (L_X/2)+Lbuffer], [(-L_Y/2)-Lbuffer, (L_Y/2)+Lbuffer], [(-L_Z/2)-Lbuffer, (L_Z/2)+Lbuffer]], # limits 
      2.0, # block size
      radii=frameradii_array_c1 # matching list of particle radii -- required for polydisperse / radical tessellation 
    )

    cells_c2 = pyvoro.compute_voronoi(
      framepos_c2, # point positions
      [[(-L_X/2)-Lbuffer, (L_X/2)+Lbuffer], [(-L_Y/2)-Lbuffer, (L_Y/2)+Lbuffer], [(-L_Z/2)-Lbuffer, (L_Z/2)+Lbuffer]], # limits 
      2.0, # block size
      radii=frameradii_array_c2 # matching list of particle radii -- required for polydisperse / radical tessellation 
    )
    
    volumes_all = np.zeros(len(cells_all)) 
    for j in range(len(volumes_all)):
      volumes_all[j] = cells_all[j].get("volume")

    free_volumes_all = np.zeros(len(volumes_all))
    for j in range(len(free_volumes_all)):
      free_volumes_all[j] = volumes_all[j] - frame_particle_volumes[j] 

    volumes_c1 = np.zeros(len(cells_c1))
    for j in range(len(volumes_c1)):
      volumes_c1[j] = cells_c1[j].get("volume")

    free_volumes_c1 = np.zeros(len(volumes_c1))
    for j in range(len(free_volumes_c1)):
      free_volumes_c1[j] = volumes_c1[j] - frame_particle_volumes_c1[j] 

    volumes_c2 = np.zeros(len(cells_c2))
    for j in range(len(volumes_c2)):
      volumes_c2[j] = cells_c2[j].get("volume")

    free_volumes_c2 = np.zeros(len(volumes_c2))
    for j in range(len(free_volumes_c2)):
      free_volumes_c2[j] = volumes_c2[j] - frame_particle_volumes_c2[j] 

 
    for particle in range(ncolloids):
      f_all.write('{0} {1} {2} {3} {4} {5}\n'.format(int(frame), particle, types_all[particle], frameradii_array[particle], volumes_all[particle], free_volumes_all[particle]))	

    for particle in range(ncolloid1):
      f_c1.write('{0} {1} {2} {3} {4} {5}\n'.format(int(frame), particle, types_1[particle], frameradii_array_c1[particle], volumes_c1[particle], free_volumes_c1[particle]))

    for particle in range(ncolloid2):
      f_c2.write('{0} {1} {2} {3} {4} {5}\n'.format(int(frame), particle, types_2[particle], frameradii_array_c2[particle], volumes_c2[particle], free_volumes_c2[particle]))

    #print(vorvols[1])

  print('Voronoi volume calculation complete for '+str(nframes)+' frames and '+str(ncolloids)+' colloid particles ('+str(ntypes)+'colloid types).')
#######


######
## Voronoi stats
"""
# calculate mean (1st mode), standard deviation (2nd mode), skew (3rd mode), 
# kurtosis and excess kurtosis (4th mode) of the Voronoi volume distribution
"""

def voronoi_stats_py(filename):
  # get volume information
  traj = gsd.hoomd.open(filename, 'r')
  Lbox = traj[-1].configuration.box[:3]
  L_X = Lbox[0]
  L_Y = Lbox[1]
  L_Z = Lbox[2]
  real_vol = L_X * L_Y * L_Z
  buffered_vol = (L_X+(2*Lbuffer)) * (L_Y+(2*Lbuffer)) * (L_Z+(2*Lbuffer))

  f = np.genfromtxt(data_outpath+'/voronoi-volumes.txt', skip_header=1)
  f_c1 = np.genfromtxt(data_outpath+'/voronoi-volumes-c1only.txt', skip_header=1)
  f_c2 = np.genfromtxt(data_outpath+'/voronoi-volumes-c2only.txt', skip_header=1)
  nframes = len(np.unique(f[:,0]))
  ncolloids = len(np.unique(f[:,1]))
  ncolloid1 = len(np.unique(f_c1[:,1]))
  ncolloid2 = len(np.unique(f_c2[:,1]))
  if (ncolloids != ncolloid1+ncolloid2):
    print("ERROR: Total number of colloids does not equal the sum of subpopulations in Voronoi stats; stats NOT calculated")
    return

  # reshape data into a 3D array of volumes
  # access individual volumes with: vorvol_array[frame-number][colloidID-number][0]
  vorvol_array = np.reshape(f[:,4],(nframes,ncolloids,1))
  free_vorvol_array = np.reshape(f[:,5],(nframes,ncolloids,1))

  vorvol_array_c1 = np.reshape(f_c1[:,4],(nframes,ncolloid1,1))
  free_vorvol_array_c1 = np.reshape(f_c1[:,5],(nframes,ncolloid1,1))

  vorvol_array_c2 = np.reshape(f_c2[:,4],(nframes,ncolloid2,1))
  free_vorvol_array_c2 = np.reshape(f_c2[:,5],(nframes,ncolloid2,1))

  # 1. Calculate the mean volume for each frame (1st mode/moment)
  # The mean volume is the same throughout (total volume and number of particles do not change)
  mean_volume = np.mean(vorvol_array, axis=(1, 2)) * (real_vol / buffered_vol)
  mean_free_volume = np.mean(free_vorvol_array, axis=(1, 2)) * (real_vol / buffered_vol)
  mean_volume_c1 = np.mean(vorvol_array_c1, axis=(1, 2))
  mean_free_volume_c1 = np.mean(free_vorvol_array_c1, axis=(1, 2))
  mean_volume_c2 = np.mean(vorvol_array_c2, axis=(1, 2))
  mean_free_volume_c2 = np.mean(free_vorvol_array_c2, axis=(1, 2))

  # 2. Calculate the divergence of volumes for each frame (2nd mode/moment)
  sd_volume = np.std(vorvol_array, axis=(1, 2))
  sd_free_volume = np.std(free_vorvol_array, axis=(1, 2))
  sd_volume_c1 = np.std(vorvol_array_c1, axis=(1, 2))
  sd_free_volume_c1 = np.std(free_vorvol_array_c1, axis=(1, 2))
  sd_volume_c2 = np.std(vorvol_array_c2, axis=(1, 2))
  sd_free_volume_c2 = np.std(free_vorvol_array_c2, axis=(1, 2))

  # 3. Calculate the skewness of the volume distribution (3rd mode/moment)
  # How asymmetrical is the distribution?
  # NEGATIVE -- shift right from normal distribution
  # POSITIVE -- shift left from normal distribution
  skew_volume = np.zeros(nframes)
  skew_free_volume = np.zeros(nframes)
  skew_volume_c1 = np.zeros(nframes)
  skew_free_volume_c1 = np.zeros(nframes)
  skew_volume_c2 = np.zeros(nframes)
  skew_free_volume_c2 = np.zeros(nframes)
  for i in range(nframes):
    skew_volume[i] = skew(vorvol_array[i,:,0])
    skew_free_volume[i] = skew(free_vorvol_array[i,:,0])
    skew_volume_c1[i] = skew(vorvol_array_c1[i,:,0])
    skew_free_volume_c1[i] = skew(free_vorvol_array_c1[i,:,0])
    skew_volume_c2[i] = skew(vorvol_array_c2[i,:,0])
    skew_free_volume_c2[i] = skew(free_vorvol_array_c2[i,:,0])

  # 4. Calculate the kurtosis of the volumes (4th mode/moment)
  # Is there an outlier problem?
  # Two kinds:
  #   - Pearson kurtosis: unscaled "bare" 4th moment, better for multivariate distributions
  #     (i.e. comparing relationships between two or more measurements)
  #   - Fisher kurtosis AKA excess kurtosis: kurtosis deviation from a normal distribution
  # NEGATIVE == thin tails; fewer and/or less extreme outliers than the normal distribution
  # POSITIVE == fat tails; more outliers than the normal distribution
  # a. Pearson
  kurtosis_volume = np.zeros(nframes)
  kurtosis_free_volume = np.zeros(nframes)
  kurtosis_volume_c1 = np.zeros(nframes)
  kurtosis_free_volume_c1 = np.zeros(nframes)
  kurtosis_volume_c2 = np.zeros(nframes)
  kurtosis_free_volume_c2 = np.zeros(nframes)
  # b. Fisher
  excesskurtosis_volume = np.zeros(nframes)
  excesskurtosis_free_volume = np.zeros(nframes)
  excesskurtosis_volume_c1 = np.zeros(nframes)
  excesskurtosis_free_volume_c1 = np.zeros(nframes)
  excesskurtosis_volume_c2 = np.zeros(nframes)
  excesskurtosis_free_volume_c2 = np.zeros(nframes)
  for i in range(nframes):
    # a. Pearson
    kurtosis_volume[i] = kurtosis(vorvol_array[i,:,0], fisher=False)
    kurtosis_free_volume[i] = kurtosis(free_vorvol_array[i,:,0], fisher=False)
    kurtosis_volume_c1[i] = kurtosis(vorvol_array_c1[i,:,0], fisher=False)
    kurtosis_free_volume_c1[i] = kurtosis(free_vorvol_array_c1[i,:,0], fisher=False)
    kurtosis_volume_c2[i] = kurtosis(vorvol_array_c2[i,:,0], fisher=False)
    kurtosis_free_volume_c2[i] = kurtosis(free_vorvol_array_c2[i,:,0], fisher=False)
    # b. Fisher
    excesskurtosis_volume[i] = kurtosis(vorvol_array[i,:,0])
    excesskurtosis_free_volume[i] = kurtosis(free_vorvol_array[i,:,0])
    excesskurtosis_volume_c1[i] = kurtosis(vorvol_array_c1[i,:,0])
    excesskurtosis_free_volume_c1[i] = kurtosis(free_vorvol_array_c1[i,:,0])
    excesskurtosis_volume_c2[i] = kurtosis(vorvol_array_c2[i,:,0])
    excesskurtosis_free_volume_c2[i] = kurtosis(free_vorvol_array_c2[i,:,0])

  # 5. convert to DF to save
  stats_df = pd.DataFrame(np.arange(nframes), columns=['frame'])
  stats_df['totalv-mean'] = mean_volume
  stats_df['freev-mean'] = mean_volume
  stats_df['totalv-divergence'] = sd_volume
  stats_df['freev-divergence'] = sd_volume
  stats_df['totalv-skew'] = skew_volume
  stats_df['freev-skew'] = skew_volume
  stats_df['totalv-kurtosis'] = kurtosis_volume
  stats_df['freev-kurtosis'] = kurtosis_volume
  stats_df['totalv-excess-kurtosis'] = excesskurtosis_volume
  stats_df['freev-excess-kurtosis'] = excesskurtosis_volume

  stats_df.to_csv(data_outpath+'/voronoi-stats.csv', index=False)

  print('Mean, divergence, skew, kurtosis, and excess kurtosis caluclated for '+str(nframes)+' frames and '+str(ncolloids)+' monodisperse colloids.')

  # network of only colloid 1
  stats_df_c1 = pd.DataFrame(np.arange(nframes), columns=['frame'])
  stats_df_c1['totalv-mean'] = mean_volume_c1
  stats_df_c1['freev-mean'] = mean_volume_c1
  stats_df_c1['totalv-divergence'] = sd_volume_c1
  stats_df_c1['freev-divergence'] = sd_volume_c1
  stats_df_c1['totalv-skew'] = skew_volume_c1
  stats_df_c1['freev-skew'] = skew_volume_c1
  stats_df_c1['totalv-kurtosis'] = kurtosis_volume_c1
  stats_df_c1['freev-kurtosis'] = kurtosis_volume_c1
  stats_df_c1['totalv-excess-kurtosis'] = excesskurtosis_volume_c1
  stats_df_c1['freev-excess-kurtosis'] = excesskurtosis_volume_c1

  stats_df_c1.to_csv(data_outpath+'/voronoi-stats-c1only.csv', index=False)

  print('Mean, divergence, skew, kurtosis, and excess kurtosis caluclated for '+str(nframes)+' frames and a network of '+str(ncolloid1)+' type 1 colloids.')

  # network of only colloid 2
  stats_df_c2 = pd.DataFrame(np.arange(nframes), columns=['frame'])
  stats_df_c2['totalv-mean'] = mean_volume_c2
  stats_df_c2['freev-mean'] = mean_volume_c2
  stats_df_c2['totalv-divergence'] = sd_volume_c2
  stats_df_c2['freev-divergence'] = sd_volume_c2
  stats_df_c2['totalv-skew'] = skew_volume_c2
  stats_df_c2['freev-skew'] = skew_volume_c2
  stats_df_c2['totalv-kurtosis'] = kurtosis_volume_c2
  stats_df_c2['freev-kurtosis'] = kurtosis_volume_c2
  stats_df_c2['totalv-excess-kurtosis'] = excesskurtosis_volume_c2
  stats_df_c2['freev-excess-kurtosis'] = excesskurtosis_volume_c2

  stats_df_c2.to_csv(data_outpath+'/voronoi-stats-c2only.csv', index=False)

  print('Mean, divergence, skew, kurtosis, and excess kurtosis caluclated for '+str(nframes)+' frames and a network of '+str(ncolloid2)+' type 2 colloids.')

######


######
## Voronoi distribution
"""
# count/bin how many Voronoi volumes exist at each volume size
"""

def voronoi_distribution_py(filename):

  f = np.genfromtxt(data_outpath+'/voronoi-volumes.txt', skip_header=1)
  f_c1 = np.genfromtxt(data_outpath+'/voronoi-volumes-c1only.txt', skip_header=1)
  f_c2 = np.genfromtxt(data_outpath+'/voronoi-volumes-c2only.txt', skip_header=1)

  # create dataframes with the frames as integer values
  voro_df = pd.DataFrame(f[:,0], columns=['frame'], dtype=int)
  voro_c1_df = pd.DataFrame(f_c1[:,0], columns=['frame'], dtype=int)
  voro_c2_df = pd.DataFrame(f_c2[:,0], columns=['frame'], dtype=int)

  # add remaining data to the dataframes
  voro_df["particles"] = f[:,1]
  voro_df["tyepid"] = f[:,2]
  voro_df["radius"] = f[:,3]
  voro_df["voronoi_v"] = f[:,4]
  voro_df["free_v"] = f[:,5]

  particle_volumes = (4/3)*np.pi*np.array(f[:,3])**3
  scaled_free_volumes = np.array(f[:,5])/particle_volumes
  voro_df["free_v_scaled"] = scaled_free_volumes

  voro_c1_df["particles"] = f_c1[:,1]
  voro_c1_df["typeid"] = f_c1[:,2]
  voro_c1_df["radius"] = f_c1[:,3]
  voro_c1_df["voronoi_v"] = f_c1[:,4]
  voro_c1_df["free_v"] = f_c1[:,5]
  voro_c1_df["free_v_scaled"] = f_c1[:,5]

  particle_volumes_c1 = (4/3)*np.pi*np.array(f_c1[:,3])**3
  scaled_free_volumes_c1 = np.array(f_c1[:,5])/particle_volumes_c1
  voro_c1_df["free_v_scaled"] = scaled_free_volumes_c1

  voro_c2_df["particles"] = f_c2[:,1]
  voro_c2_df["typeid"] = f_c2[:,2]
  voro_c2_df["radius"] = f_c2[:,3]
  voro_c2_df["voronoi_v"] = f_c2[:,4]
  voro_c2_df["free_v"] = f_c2[:,5]
  voro_c2_df["free_v_scaled"] = f_c2[:,5]

  particle_volumes_c2 = (4/3)*np.pi*np.array(f_c2[:,3])**3
  scaled_free_volumes_c2 = np.array(f_c2[:,5])/particle_volumes_c2
  voro_c2_df["free_v_scaled"] = scaled_free_volumes_c2

  # get number of frames
  framelist = pd.unique(voro_df['frame'])
  nframes = len(framelist)

  # get the number of colloids
  ncolloids = len(pd.unique(voro_df['particles']))
  ncolloid1 = len(pd.unique(voro_c1_df['particles']))
  ncolloid2 = len(pd.unique(voro_c2_df['particles']))

  # create bins for the volumes
  max_vol_all = max(voro_df['voronoi_v'].to_numpy())
  max_vol_c1 = max(voro_c1_df['voronoi_v'].to_numpy())
  max_vol_c2 = max(voro_c2_df['voronoi_v'].to_numpy())
  max_vol = max(max_vol_all, max_vol_c1, max_vol_c2)

  max_free_vol_all = max(voro_df['free_v'].to_numpy())
  max_free_vol_c1 = max(voro_c1_df['free_v'].to_numpy())
  max_free_vol_c2 = max(voro_c2_df['free_v'].to_numpy())
  max_free_vol = max(max_free_vol_all, max_free_vol_c1, max_free_vol_c2)

  max_free_vol_scaled_all = max(voro_df['free_v_scaled'].to_numpy())
  max_free_vol_scaled_c1 = max(voro_c1_df['free_v_scaled'].to_numpy())
  max_free_vol_scaled_c2 = max(voro_c2_df['free_v_scaled'].to_numpy())
  max_free_vol_scaled = max(max_free_vol_scaled_all, max_free_vol_scaled_c1, max_free_vol_scaled_c2)


  """
  #### MANUAL BINNING
  bins_voro = np.linspace(0,max_vol,int(max_vol+1)).astype(int)
  bins_free = np.linspace(0,max_free_vol,int(max_free_vol+1)).astype(int)
  bins_free_scaled = np.linspace(0,max_free_vol_scaled,int(max_free_vol_scaled+1)).astype(int)

  # make array to count particles per bin per frame
  histbins_voro = np.zeros((nframes,int(max_vol+1),3))
  histbins_free = np.zeros((nframes,int(max_free_vol+1),3))
  histbins_free_scaled = np.zeros((nframes,int(max_free_vol_scaled+1),3))

  # make array for data normalized by total number of particles
  normcounts_voro = np.zeros((nframes,int(max_vol+1),3))
  normcounts_free = np.zeros((nframes,int(max_free_vol+1),3))
  normcounts_free_scaled = np.zeros((nframes,int(max_free_vol_scaled+1),3))

  # loop through the unique frames
  for k in range(nframes):
    frame = framelist[k]
    # create a new data frame for each time
    currframe_all = voro_df.loc[voro_df['frame'] == frame]
    currframe_c1 = voro_c1_df.loc[voro_c1_df['frame'] == frame]
    currframe_c2 = voro_c2_df.loc[voro_c2_df['frame'] == frame]

    # extract the volumes
    voro_volumes_all = list(currframe_all['voronoi_v'])
    free_volumes_all = list(currframe_all['free_v'])
    free_volumes_scaled_all = list(currframe_all['free_v_scaled'])
    voro_volumes_c1 = list(currframe_c1['voronoi_v'])
    free_volumes_c1 = list(currframe_c1['free_v'])
    free_volumes_scaled_c1 = list(currframe_c1['free_v_scaled'])
    voro_volumes_c2 = list(currframe_c2['voronoi_v'])
    free_volumes_c2 = list(currframe_c2['free_v'])
    free_volumes_scaled_c2 = list(currframe_c2['free_v_scaled'])

    # bin the volumes
    for curr_bin in range(0,int(max_vol)):

      for j in range(len(voro_volumes_all)):
        if curr_bin <= voro_volumes_all[j] < curr_bin+1:
          histbins_voro[k,curr_bin,0] += 1
      for j in range(len(free_volumes_all)):
        if curr_bin <= free_volumes_all[j] < curr_bin+1:
          histbins_free[k,curr_bin,0] += 1
      for j in range(len(free_volumes_scaled_all)):
        if curr_bin <= free_volumes_scaled_all[j] < curr_bin+1:
          histbins_free_scaled[k,curr_bin,0] += 1

      for s in range(len(voro_volumes_c1)):
        if curr_bin <= voro_volumes_c1[s] < curr_bin+1:
          histbins_voro[k,curr_bin,1] += 1
      for s in range(len(free_volumes_c1)):
        if curr_bin <= free_volumes_c1[s] < curr_bin+1:
          histbins_free[k,curr_bin,1] += 1
      for s in range(len(free_volumes_scaled_c1)):
        if curr_bin <= free_volumes_scaled_c1[s] < curr_bin+1:
          histbins_free_scaled[k,curr_bin,1] += 1

      for l in range(len(voro_volumes_c2)):
        if curr_bin <= voro_volumes_c2[l] < curr_bin+1:
          histbins_voro[k,curr_bin,2] += 1
      for l in range(len(free_volumes_c2)):
        if curr_bin <= free_volumes_c2[l] < curr_bin+1:
          histbins_free[k,curr_bin,2] += 1
      for l in range(len(free_volumes_scaled_c2)):
        if curr_bin <= free_volumes_scaled_c2[l] < curr_bin+1:
          histbins_free_scaled[k,curr_bin,2] += 1

    # normalize the bin values (0-1)
    normcounts_voro[k,:,0] = histbins_voro[k,:,0]/ncolloids
    normcounts_free[k,:,0] = histbins_free[k,:,0]/ncolloids
    normcounts_free_scaled[k,:,0] = histbins_free_scaled[k,:,0]/ncolloids
    normcounts_voro[k,:,1] = histbins_voro[k,:,1]/ncolloid1
    normcounts_free[k,:,1] = histbins_free[k,:,1]/ncolloid1
    normcounts_free_scaled[k,:,1] = histbins_free_scaled[k,:,1]/ncolloid1
    normcounts_voro[k,:,2] = histbins_voro[k,:,2]/ncolloid2
    normcounts_free[k,:,2] = histbins_free[k,:,2]/ncolloid2
    normcounts_free_scaled[k,:,2] = histbins_free_scaled[k,:,2]/ncolloid2
  """

  ##### AUTO BINNING
  nbins = 100
  histbins_voro = np.zeros((nframes,nbins,3))
  histbins_free = np.zeros((nframes,nbins,3))
  histbins_free_scaled = np.zeros((nframes,nbins,3))
  normcounts_voro = np.zeros((nframes,nbins,3))
  normcounts_free = np.zeros((nframes,nbins,3))
  normcounts_free_scaled = np.zeros((nframes,nbins,3))

  # loop through the unique frames
  for k in range(nframes):
    frame = framelist[k]
    # create a new data frame for each time
    currframe_all = voro_df.loc[voro_df['frame'] == frame]
    currframe_c1 = voro_c1_df.loc[voro_c1_df['frame'] == frame]
    currframe_c2 = voro_c2_df.loc[voro_c2_df['frame'] == frame]

    # extract the volumes
    voro_volumes_all = list(currframe_all['voronoi_v'])
    free_volumes_all = list(currframe_all['free_v'])
    free_volumes_scaled_all = list(currframe_all['free_v_scaled'])
    voro_volumes_c1 = list(currframe_c1['voronoi_v'])
    free_volumes_c1 = list(currframe_c1['free_v'])
    free_volumes_scaled_c1 = list(currframe_c1['free_v_scaled'])
    voro_volumes_c2 = list(currframe_c2['voronoi_v'])
    free_volumes_c2 = list(currframe_c2['free_v'])
    free_volumes_scaled_c2 = list(currframe_c2['free_v_scaled']) 

    histbins_voro_all_frame, bins_voro = np.histogram(voro_volumes_all, bins=nbins, range=(0,max_vol+1))
    histbins_voro_c1_frame, _ = np.histogram(voro_volumes_c1, bins=nbins, range=(0,max_vol+1))
    histbins_voro_c2_frame, _ = np.histogram(voro_volumes_c2, bins=nbins, range=(0,max_vol+1))
    normcounts_voro_all_frame = histbins_voro_all_frame/ncolloids
    normcounts_voro_c1_frame = histbins_voro_c1_frame/ncolloid1
    normcounts_voro_c2_frame = histbins_voro_c2_frame/ncolloid2
    histbins_voro[k,:,0] = histbins_voro_all_frame
    histbins_voro[k,:,1] = histbins_voro_c1_frame
    histbins_voro[k,:,2] = histbins_voro_c2_frame
    normcounts_voro[k,:,0] = normcounts_voro_all_frame
    normcounts_voro[k,:,1] = normcounts_voro_c1_frame
    normcounts_voro[k,:,2] = normcounts_voro_c2_frame

    histbins_free_all_frame, bins_free = np.histogram(free_volumes_all, bins=nbins, range=(0,max_free_vol+1))
    histbins_free_c1_frame, bins_free = np.histogram(free_volumes_c1, bins=nbins, range=(0,max_free_vol+1))
    histbins_free_c2_frame, bins_free = np.histogram(free_volumes_c2, bins=nbins, range=(0,max_free_vol+1))
    normcounts_free_all_frame = histbins_free_all_frame/ncolloids
    normcounts_free_c1_frame = histbins_free_c1_frame/ncolloid1
    normcounts_free_c2_frame = histbins_free_c2_frame/ncolloid2
    histbins_free[k,:,0] = histbins_free_all_frame
    histbins_free[k,:,1] = histbins_free_c1_frame
    histbins_free[k,:,2] = histbins_free_c2_frame
    normcounts_free[k,:,0] = normcounts_free_all_frame
    normcounts_free[k,:,1] = normcounts_free_c1_frame
    normcounts_free[k,:,2] = normcounts_free_c2_frame

    histbins_freescaled_all_frame, bins_free_scaled = np.histogram(free_volumes_scaled_all, bins=nbins, range=(0,max_free_vol_scaled+1))
    histbins_freescaled_c1_frame, bins_free_scaled = np.histogram(free_volumes_scaled_c1, bins=nbins, range=(0,max_free_vol_scaled+1))
    histbins_freescaled_c2_frame, bins_free_scaled = np.histogram(free_volumes_scaled_c2, bins=nbins, range=(0,max_free_vol_scaled+1))
    normcounts_freescaled_all_frame = histbins_freescaled_all_frame/ncolloids
    normcounts_freescaled_c1_frame = histbins_freescaled_c1_frame/ncolloid1
    normcounts_freescaled_c2_frame = histbins_freescaled_c2_frame/ncolloid2
    histbins_free_scaled[k,:,0] = histbins_freescaled_all_frame
    histbins_free_scaled[k,:,1] = histbins_freescaled_c1_frame
    histbins_free_scaled[k,:,2] = histbins_freescaled_c2_frame
    normcounts_free_scaled[k,:,0] = normcounts_freescaled_all_frame
    normcounts_free_scaled[k,:,1] = normcounts_freescaled_c1_frame
    normcounts_free_scaled[k,:,2] = normcounts_freescaled_c2_frame

  bins_voro = bins_voro[0:-1]
  bins_free = bins_free[0:-1]
  bins_free_scaled = bins_free_scaled[0:-1]


  # convert arrays to data frames for easy saving 
  framedf_voro_list = []

  for k in range(nframes):
    frame = framelist[k]

    frame_df = pd.DataFrame(histbins_voro[k,:])
    frame_df = frame_df.rename(columns = {0:'allc-counts'})
    frame_df = frame_df.rename(columns = {1:'c1-counts'})
    frame_df = frame_df.rename(columns = {2:'c2-counts'})
    frame_df.insert(0, "vol-bins", bins_voro)
    frame_df.insert(0, "frame", frame)
    frame_df['allc-norm'] = normcounts_voro[frame,:,0]
    frame_df['c1-norm'] = normcounts_voro[frame,:,1]
    frame_df['c2-norm'] = normcounts_voro[frame,:,2]
    framedf_voro_list.append(frame_df)

  vorvols_binned = pd.concat(framedf_voro_list)

  vorvols_binned.to_csv(data_outpath+'/voronoi_binned.csv', index=False)


  framedf_free_list = []

  for k in range(nframes):
    frame = framelist[k]

    frame_df = pd.DataFrame(histbins_free[k,:])
    frame_df = frame_df.rename(columns = {0:'allc-counts'})
    frame_df = frame_df.rename(columns = {1:'c1-counts'})
    frame_df = frame_df.rename(columns = {2:'c2-counts'})
    frame_df.insert(0, "vol-bins", bins_free)
    frame_df.insert(0, "frame", frame)
    frame_df['allc-norm'] = normcounts_free[frame,:,0]
    frame_df['c1-norm'] = normcounts_free[frame,:,1]
    frame_df['c2-norm'] = normcounts_free[frame,:,2]
    framedf_free_list.append(frame_df)

  freevols_binned = pd.concat(framedf_free_list)

  freevols_binned.to_csv(data_outpath+'/freevol_binned.csv', index=False)


  framedf_free_scaled_list = []

  for k in range(nframes):
    frame = framelist[k]

    frame_df = pd.DataFrame(histbins_free_scaled[k,:])
    frame_df = frame_df.rename(columns = {0:'allc-counts'})
    frame_df = frame_df.rename(columns = {1:'c1-counts'})
    frame_df = frame_df.rename(columns = {2:'c2-counts'})
    frame_df.insert(0, "vol-bins", bins_free_scaled)
    frame_df.insert(0, "frame", frame)
    frame_df['allc-norm'] = normcounts_free_scaled[frame,:,0]
    frame_df['c1-norm'] = normcounts_free_scaled[frame,:,1]
    frame_df['c2-norm'] = normcounts_free_scaled[frame,:,2]
    framedf_free_scaled_list.append(frame_df)

  freevols_scaled_binned = pd.concat(framedf_free_scaled_list)

  freevols_scaled_binned.to_csv(data_outpath+'/freevol_scaled_binned.csv', index=False)


  print("Voronoi volumes and free volumes binned and normalized by number of particles for 2 particle types in "+str(nframes)+" frames")
######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':
  posCSV_calc(filepath)	
  voronoi_volume_py(filepath)
  voronoi_stats_py(filepath)
  voronoi_distribution_py(filepath)

  print("\nAll Voronoi analyses complete.")
