## Extract data to CSV and run primary network analysis on the results 
## of a DPD colloid simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 solvent type (typeid=0) and 
##		             1 colloid type (typeid=1)
##
"""
## This code performs the following analyses:
## - Extracts colloid position data to from GSD to CSV files for all frames
## - Extracts the list of all bonded colloid pairs (i.e network edges) from
##     GSD to CSV files for all frames
## - Calculates primary network analysis metrics for all frames
## 	- number of connected components 
##	- average degree (AKA average coordination number in gels)
##      - largest connected component (LCC)
##      - average clustering coefficient
##      - average square clustering coefficient 
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import gsd.hoomd
import math
import fortranmod as module 
import pandas as pd
import networkx as nx
import os
import re
from statistics import mean
from scipy.spatial.distance import pdist
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
    filepath, data_outpath, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_netx_variables(data_path_name, gsd_sim_path)
  elif len(sys.argv) > 1:
    data_path_name = sys.argv[1]
    filepath, data_outpath, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_netx_variables(data_path_name)
  else:
    filepath, data_outpath, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_netx_variables()

  colloid_typeid = colloid1_typeid

elif from_path == False:
  # manual simulation parameters
  filepath = '../Gelation-DPD.gsd'
  data_outpath = 'data'
  colloid_typeid = 1
  kappa = 30

cut_off = round(3/kappa,2)


###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)

# define a function to recenter the dims on the center of mass, accounting for periodic boundaries in the cluster length
def pbc_centered_positions(positions, box):
    """
    Center the positions around their center of mass (COM),
    accounting for periodic boundaries using the minimum image convention.

    Returns:
    - centered positions (N, 3), shifted such that COM is at the origin (0,0,0)
    """
    # Convert positions to [0, L) for safety
    positions = positions % box

    # Compute COM using angles (Fourier trick for PBC-aware averaging)
    theta = 2 * np.pi * positions / box
    com_theta = np.arctan2(np.mean(np.sin(theta), axis=0),
                           np.mean(np.cos(theta), axis=0))
    com = (com_theta % (2 * np.pi)) * box / (2 * np.pi)

    # Shift all positions relative to COM using minimum image convention
    disp = positions - com
    disp -= box * np.round(disp / box)  # Apply minimum image convention
    return disp

#######
## posCSV
"""
# create a CSV of particle (i.e. node) position information for each frame
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
      print('position data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [1])[0]
  ncolloids = len(colloids)
  typeid = traj[-1].particles.typeid[colloids]
  radii = 0.5*traj[-1].particles.diameter[colloids]
  rcut = cut_off
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
## edgelistCSV
"""
# create a CSV file of all the bonded particle pairs (i.e. edges of the network) for each frame
# i,j position for each bond/edge
"""

def edgelistCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  
  # set path and filename
  edge_dir_path = data_outpath+'/frame-edges'
  edge_output = edge_dir_path+'/edgelist' # + <#>.csv in python loop
  
  # check for existing CSV data 
  if os.path.exists(edge_dir_path) == False:
    os.mkdir(edge_dir_path)
  if os.path.exists(edge_dir_path) == True:
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'edgelist\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(edge_dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1
    
    if nframes == len(traj):
      print('edgelist data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [1])[0]
  ncolloids = len(colloids)
  radii = 0.5*traj[-1].particles.diameter[colloids]
  rcut = cut_off
  lbox = traj[-1].configuration.box[:3]
  inv_lbox = 1/lbox

  # create an array of xyz positon of all colloids in all frames    
  allpos = np.zeros((nframes,ncolloids,3))
  for i in range(0,nframes):
    allpos[i,:,:] = traj[i].particles.position[colloids] 
 
  module.edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)
  print("Edgelist calculation complete for "+str(nframes)+" frames")
#######


#######
## networkx analysis
"""
# use the networkx package to calculate for all frames:
#   - the number of connected components
#   - average degree (AKA average coordination number)
#   - largest connected component (LCC)
#   - average clustering coefficient
#   - average square clustering coefficient
"""

def primary_networkx_calc(filename):

  edge_output = data_outpath+'/frame-edges/edgelist' # + <#>.csv in f90

  # get the number of particles
  traj = gsd.hoomd.open(filename, 'r')
  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [1])[0]
  ncolloids = len(colloids)

  # get radii for backbone classification
  radii = 0.5*traj[-1].particles.diameter[colloids]

  # create empty dataframe for all network data  
  network_frames_df = pd.DataFrame()
  
  # import all data into one dataframe
  frame_dfs = []
  for frame in range(nframes):
    # loop through all frames
    filepath = edge_output+str(frame)+'.csv'
    # import CSV data
    df = pd.read_csv(filepath)
    # rename colums as needed
    df = df.rename(columns={"i": "source", "j": "target"})
    df.insert(loc=0, column='frame', value=frame)
    frame_dfs.append(df)

  alledge_df = pd.concat(frame_dfs, ignore_index=True) 

  # network analysis
  for frame in range(nframes):

    # get particle positions for percolation measurement
    pos = traj[frame].particles.position
    # box diagonal (for percolation)
    Lbox = traj[frame].configuration.box[:3]
    max_diameter = np.sqrt(Lbox[0]**2 + Lbox[1]**2 + Lbox[2]**2)

    df = alledge_df[alledge_df['frame'] == frame][["source", "target"]]

    # create the network from edge list
    g = nx.from_pandas_edgelist(df)

    # if a node is not in the network, add it
    for particle in range(ncolloids):
      if (  not(   g.has_node(particle)   )  ):
        g.add_node(particle)

    # calculate the total number of edges      
    nedges = nx.number_of_edges(g)

    # calculate the average degree of the network, i.e. avg contact number in gels
    avg_degree = 2 * nedges / ncolloids

    # number of connected components
    n_cc = nx.number_connected_components(g)
    components = list(nx.connected_components(g))

    # return the indices of the nodes of the largest connected components
    lcc_nodes = max(components, key=len)
    #lcc_nodes = max(nx.connected_components(g), key=len)

    # calculate the size of the largest connected component
    lcc_size = len(lcc_nodes)

    # --- species composition of network / LCC / backbone (allc only) ---
    if len(pd.unique(radii)) != 1:
        print(f'WARNING: system contains more than one size of particle: {pd.unique(radii)}, backbone classification is only designed to handle monodisperse size 1 and size 2')
    size = pd.unique(radii)[0] 
    if int(size) == 1:
        # LCC composition (fraction of LCC that is type 1)
        lcc_frac_type1 = 1.0
        bond_frac_11 = 1.0
        bond_frac_12 = 1.0
        bond_frac_22 = 0.0

        # 2-core = degree-based backbone proxy (iteratively strip nodes with deg<2)
        # pruned dangling ends = nodes with degree exactly 1
        core2 = nx.k_core(g, k=2)
        core_nodes = list(core2.nodes())
        core2_size = len(core_nodes)
        core2_frac_type1 = 1.0 if core2_size > 0 else np.nan
        degs = dict(g.degree())
        dangling = [n for n, d in degs.items() if d == 1]
        n_dangling = len(dangling)
        dangling_frac_type1 = 1.0 if n_dangling > 0 else np.nan

    elif int(size) == 2:
        # LCC composition (fraction of LCC that is type 1)
        lcc_frac_type1 = 0.0
        bond_frac_11 = 0.0
        bond_frac_12 = 0.0
        bond_frac_22 = 1.0

        # 2-core = degree-based backbone proxy (iteratively strip nodes with deg<2)
        # pruned dangling ends = nodes with degree exactly 1
        core2 = nx.k_core(g, k=2)
        core_nodes = list(core2.nodes())
        core2_size = len(core_nodes)
        core2_frac_type1 = 0.0 if core2_size > 0 else np.nan
        degs = dict(g.degree())
        dangling = [n for n, d in degs.items() if d == 1]
        n_dangling = len(dangling)
        dangling_frac_type1 = 0.0 if n_dangling > 0 else np.nan

    else:
        print(f'WARNING: particle size "{size}" not recognized, no backbone topology recorded')
        lcc_frac_type1 = bond_frac_11 = bond_frac_12 = bond_frac_22 = np.nan
        core2_size = np.nan; core2_frac_type1 = np.nan
        n_dangling = np.nan; dangling_frac_type1 = np.nan

    # --- mean cluster size (susceptibility) chi, LCC excluded ---
    # chi = sum'(s^2) / sum'(s) over all clusters EXCEPT the single largest.
    # This is the particle-weighted mean finite-cluster size; peaks at percolation.
    comp_sizes = np.array([len(c) for c in components])
    finite_sizes = np.delete(comp_sizes, np.argmax(comp_sizes))  # drop ONE largest
    finite_particles = finite_sizes.sum()
    if finite_particles > 0:
        chi = (finite_sizes.astype(float)**2).sum() / finite_particles
        mean_finite_size = finite_sizes.mean()   # = (N - lcc_size)/(n_cc - 1)
    else:
        chi = 0.0            # everything is in the LCC
        mean_finite_size = 0.0

    # physical diameter = maximum distance in the lcc
    #lcc_pos = pos_curr[list(lcc_nodes)]
    lcc_pos = pos[list(lcc_nodes)]
    # use corrected positions centered on the center of mass
    centered_pos = pbc_centered_positions(lcc_pos, Lbox)
    physical_diameter = pdist(centered_pos).max() if len(centered_pos) > 1 else 0.0

    # --- per-axis extent of the LCC (crude anisotropy-aware span) ---
    # Reuses centered_pos so it shares a basis with physical_diameter.
    # extent/L -> 1 means the LCC fills that axis (1D-spanning proxy);
    # a compact blob has all three ~equal and < 1.
    if len(centered_pos) > 1:
        axis_extent = centered_pos.max(axis=0) - centered_pos.min(axis=0)  # (3,)
    else:
        axis_extent = np.zeros(3)
    axis_span_frac = axis_extent / Lbox            # per-axis, each in [0,1]
    max_axis_span_frac = float(axis_span_frac.max())
    n_axes_spanning = int((axis_span_frac > 0.85).sum())  # threshold tunable

    # calculate the avg clustering coefficient of the network
    avg_clustering_coeff = nx.average_clustering(g)

    # calculate the avg square clustering coefficient of the network
    square_clustering_per_node = nx.square_clustering(g)
    avg_sq_clustering_coeff = mean(square_clustering_per_node[jj] for jj in square_clustering_per_node)

    # compile outputs
    data={'network-type'           :['allc'],
          'frame'                  :[frame],
          'cut_off'                :[cut_off],
          'n_components'           :[n_cc], 
          'lcc_size'               :[lcc_size],
          'ncolloids'              :[ncolloids],
          'f_lcc'                  :[lcc_size/ncolloids],
          'lcc_frac_type1'         :[lcc_frac_type1],
          'bond_frac_11'           :[bond_frac_11],
          'bond_frac_12'           :[bond_frac_12],
          'bond_frac_22'           :[bond_frac_22],
          'core2_size'             :[core2_size], # topological backbone proxy, not a rigid/stress backbone
          'core2_frac_type1'       :[core2_frac_type1],
          'n_dangling'             :[n_dangling],
          'dangling_frac_type1'    :[dangling_frac_type1],
          'lcc_span'               :[physical_diameter],
          'box_span'               :[max_diameter],
          'percolation'            :[physical_diameter/max_diameter],
          'chi_meanclustersize'    :[chi],
          'mean_finite_size'       :[mean_finite_size],
          'lcc_extent_x'           :[float(axis_span_frac[0])],
          'lcc_extent_y'           :[float(axis_span_frac[1])],
          'lcc_extent_z'           :[float(axis_span_frac[2])],
          'max_axis_span'          :[max_axis_span_frac],
          'n_axes_spanning'        :[n_axes_spanning],
          'avg_degree'             :[avg_degree],
          'avg_clustering_coeff'   :[avg_clustering_coeff],
          'avg_sq_clustering_coeff':[avg_sq_clustering_coeff]}

    res_df = pd.DataFrame(data)
    network_frames_df = pd.concat([network_frames_df, res_df])

  # write all data to one CSV file
  network_frames_df.to_csv(data_outpath+'/networkx-allframes.csv',index = False)
  print("Primary network analysis complete")
#######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':
  posCSV_calc(filepath)	
  edgelistCSV_calc(filepath)	
  primary_networkx_calc(filepath)	
