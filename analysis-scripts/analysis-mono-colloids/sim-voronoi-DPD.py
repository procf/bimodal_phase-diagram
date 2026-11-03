## Extract data to CSV and run Voronoi volume analysis on all colloids
## in a DPD simulation
## NOTE: this code assumes 1 solvent type (typeid=0) and 
##		             1 colloid type (typeid=1)
##
"""
## This code performs the following analyses:
##    - extracts colloid position data to from GSD to CSV files for all frames
##    - calculates:
##         * Voronoi volume around each colloid in the simulation for each frame
##         * statistical measures of the Voronoi volume distribution for each frame
##            - mean (1st mode)
##            - standard deviation (2nd mode)
##            - skew (3rd mode)
##            - kurtosis and excess kurtosis (4th mode)
##         * bins the distribution of Voronoi volumes for each frame 
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

  colloid_typeid = colloid1_typeid


elif from_path == False:
  # manual simulation parameters
  filepath = '../Gelation-DPD.gsd'
  data_outpath = 'data'
  colloid_typeid = 1




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
  colloids = np.where(traj[-1].particles.typeid == [1])[0]
  ncolloids = len(colloids)
  typeid = traj[-1].particles.typeid[colloids]
  radii = 0.5*traj[-1].particles.diameter[colloids]
  rcut = 0.1
  lbox = traj[-1].configuration.box[:3]

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

  lastframe_df = allpos_df.loc[allpos_df['frame']==nframes-1]

  # number of particles
  ncolloids = len(lastframe_df.loc[lastframe_df['typeid'] == 1]) # number of colloids

  # index numbers for all particles (extracted from first column of the csv file)
  colloids = lastframe_df["tag"] 

  # positions
  allpos = allpos_df[['frame','tag','x', 'y', 'z', 'radius']] # all colloid positions

  # create files to save the voronoi volumes
  f_all=open(data_outpath+'/voronoi-volumes.txt', 'w')
  f_all.write('frame particle radius volume\n')

  #vorvols = np.zeros((nframes,ncolloids))

  for frame in range(0,nframes):

    frameradii = allpos.loc[allpos['frame']==frame]['radius']
    frameradii_array = frameradii.to_numpy()
 
    # a little number manipulation to get the first AND last frame
    pos_c = allpos.loc[allpos['frame']==frame][['x', 'y', 'z']]
    framepos_c = pos_c.to_numpy()
   
    cells_all = pyvoro.compute_voronoi(
      framepos_c, # point positions
      [[(-L_X/2)-Lbuffer, (L_X/2)+Lbuffer], [(-L_Y/2)-Lbuffer, (L_Y/2)+Lbuffer], [(-L_Z/2)-Lbuffer, (L_Z/2)+Lbuffer]], # limits 
      2.0, # block size
      radii=frameradii_array # matching list of particle radii -- required for polydisperse / radical tessellation 
    )
    
    volumes_all = np.zeros(len(cells_all)) 
    for j in range(len(volumes_all)):
      volumes_all[j] = cells_all[j].get("volume")
    
    for particle in range(ncolloids):
      f_all.write('{0} {1} {2} {3}\n'.format(int(frame), particle, frameradii_array[particle], volumes_all[particle])) 	

    #print(vorvols[1])

  print('Voronoi volume calculation complete for '+str(nframes)+' frames and '+str(ncolloids)+' colloid particles.')
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
  buffer_vol = (L_X+(2*Lbuffer)) * (L_Y+(2*Lbuffer)) * (L_Z+(2*Lbuffer))

  f = np.genfromtxt(data_outpath+'/voronoi-volumes.txt', skip_header=1)
  nframes = len(np.unique(f[:,0]))
  ncolloids = len(np.unique(f[:,1]))

  # reshape data into a 3D array of volumes
  # access individual volumes with: vorvol_array[frame-number][colloidID-number][0]
  vorvol_array = np.reshape(f[:,3],(nframes,ncolloids,1))

  # 1. Calculate the mean volume for each frame (1st mode/moment)
  # The mean volume is the same throughout (total volume and number of particles do not change)
  mean_volume = np.mean(vorvol_array, axis=(1, 2)) * (real_vol / buffer_vol)

  # 2. Calculate the divergence of volumes for each frame (2nd mode/moment)
  sd_volume = np.std(vorvol_array, axis=(1, 2))

  # 3. Calculate the skewness of the volume distribution (3rd mode/moment)
  # How asymmetrical is the distribution?
  # NEGATIVE -- shift right from normal distribution
  # POSITIVE -- shift left from normal distribution
  skew_volume = np.zeros(nframes)
  for i in range(nframes):
    skew_volume[i] = skew(vorvol_array[i,:,0])

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
  # b. Fisher
  excesskurtosis_volume = np.zeros(nframes)
  for i in range(nframes):
    # a. Pearson
    kurtosis_volume[i] = kurtosis(vorvol_array[i,:,0], fisher=False)
    # b. Fisher
    excesskurtosis_volume[i] = kurtosis(vorvol_array[i,:,0])

  # 5. convert to DF to save
  stats_df = pd.DataFrame(np.arange(nframes), columns=['frame'])
  stats_df['mean'] = mean_volume
  stats_df['divergence'] = sd_volume
  stats_df['skew'] = skew_volume
  stats_df['kurtosis'] = kurtosis_volume
  stats_df['excess-kurtosis'] = excesskurtosis_volume

  stats_df.to_csv(data_outpath+'/voronoi-stats.csv', index=False)

  print('Mean, divergence, skew, kurtosis, and excess kurtosis caluclated for '+str(nframes)+' frames and '+str(ncolloids)+' monodisperse colloids.')
######


######
## Voronoi distribution
"""
# count/bin how many Voronoi volumes exist at each volume size
"""

def voronoi_distribution_py(filename):

  f = np.genfromtxt(data_outpath+'/voronoi-volumes.txt', skip_header=1)

  # create dataframes with the frames as integer values
  voro_df = pd.DataFrame(f[:,0], columns=['frame'], dtype=int)

  # add remaining data to the dataframes
  voro_df["particles"] = f[:,1]
  voro_df["radius"] = f[:,2]
  voro_df["voronoi_v"] = f[:,3]

  # get number of frames
  framelist = pd.unique(voro_df['frame'])
  nframes = len(framelist)

  # get the number of colloids
  ncolloids = len(pd.unique(voro_df['particles']))

  # create bins for the volumes
  max_vol = max(voro_df['voronoi_v'].to_numpy())
  bins = np.linspace(0,max_vol,int(max_vol+1)).astype(int)

  # make array to count particles per bin per frame
  histbins = np.zeros((nframes,int(max_vol+1)))

  # make array for data normalized by total number of particles
  normcounts = np.zeros((nframes,int(max_vol+1)))

  # loop through the unique frames
  for k in range(nframes):
    frame = framelist[k]
    # create a new data frame for each time
    currframe = voro_df.loc[voro_df['frame'] == frame]

    # extract the volumes
    volumes = list(currframe['voronoi_v'])

    # bin the volumes
    for curr_bin in range(0,int(max_vol)):
      for j in range(len(volumes)):
        if curr_bin <= volumes[j] < curr_bin+1:
          histbins[k,curr_bin] += 1

    # normalize the bin values (0-1)
    normcounts[k,:] = histbins[k,:]/ncolloids

  # convert arrays to data frames for easy saving 
  framedf_list = []

  for k in range(nframes):
    frame = framelist[k]

    frame_df = pd.DataFrame(histbins[k,:])
    frame_df = frame_df.rename(columns = {0:'mono-counts'})
    frame_df.insert(0, "vol-bins", bins)
    frame_df.insert(0, "frame", frame)
    frame_df['mono-norm'] = normcounts[k,:] 
    framedf_list.append(frame_df)

  vorvols_binned = pd.concat(framedf_list)

  vorvols_binned.to_csv(data_outpath+'/voronoi_binned.csv', index=False)

  print("Voronoi volumes binned and normalized by number of particles for 1 particle type in "+str(nframes)+" frames")
######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':
  posCSV_calc(filepath)	
  voronoi_volume_py(filepath)
  voronoi_stats_py(filepath)
  voronoi_distribution_py(filepath)
