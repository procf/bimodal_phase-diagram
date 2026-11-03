# this is a module for configuring sim-specific variables from information stored in the current path

import numpy as np
import glob
import os
import re

"""
how to call from the module:
 import config_from_path
  filepath, data_outpath, period, kT, eta0, Lbox_shortest, R_C1, R_C2, R_min, R_max, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_gel_variables()
  filepath, data_outpath, colloid1_typeid, colloid2_typeid, kappa = config_from_path.call_netx_variables()
  filepath, data_outpath, colloid1_typeid, colloid2_typeid = config_from_path.call_voronoi_variables()
  shear_filepath, gel_filepath, data_outpath, shear_style, shear_rate_const, Amp, omega, kT, eta0, frames_per_strain, dt_Integration, theta, Iter_cycle, period, MDtime, colloid1_typeid, colloid2_typeid = config_from_path.call_shear_variables()
"""

def call_gel_variables(data_path_name=None, gsd_sim_path='..'):

  if gsd_sim_path == '..':
    print(f"[config_from_path] Current working directory: {os.getcwd()}\n")
    param_path = os.getcwd()
  else:
    print(f"[config_from_path] Simulation directory: {gsd_sim_path}\n")
    param_path = gsd_sim_path

  # select the last Gelation sim file one folder up
  file_list = []
  for filepath in glob.glob(f'{gsd_sim_path}/Gelation*.gsd'):
    if 'lastframe' not in filepath:
        file_list.append(filepath)
  if not file_list:
    print("ERROR: no GSD files found")
    exit()
  else:
    file_list.sort(key=os.path.getctime) # sort list by file creation time
    filepath = file_list[-1]

  # filepath to folder where data files will be created
  if data_path_name == None:
    data_outpath = 'data'
  else:
    data_outpath = data_path_name

  ## simulation specific parameters
  # extract parameters from most recent simulation script:
  def find_params(simfile):
    # Check if the file exists
    if not os.path.exists(simfile):
      print("ERROR: no sim files found matching", simfile)
      exit()

    # Open the file and find parameters
    with open(simfile, 'r') as f:
      period = kT = eta0 = None
      param_dict = {"period":period,
                    #"KT":kT,
                    "kT":kT,
                    "eta0":eta0
      }
      lines = f.readlines()
      for line in lines:
        for param, value in param_dict.items():
          if param + " = " in line:
            param_dict[param] = line.split('=')[1].split('#')[0].strip()
        if "L_Y = " in line:
          print(f"Warning: L_Y value is set in {simfile}, double-check Lbox_shortest")
        if "L_Z = " in line:
          print(f"Warning: L_Z value is set in {simfile}, double-check Lbox_shortest")

      # check for missing values
      for param, value in param_dict.items():
        if value == None:
          print(f"ERROR: no {param} value in {simfile}")
      if None in param_dict.values():
        quit()

      return (int(param_dict["period"]),
              #float(param_dict["KT"]),
              float(param_dict["kT"]),
              float(param_dict["eta0"])
      )

  # use find_params to get params
  sim_list = []
  for simpath in glob.glob(f'{gsd_sim_path}/sim-*.py'):
    if 'lastframe' not in simpath and 'convert-to' not in simpath and 'scaleup' not in simpath:
      sim_list.append(simpath)
  if not sim_list:
    print("ERROR: no sim scripts found to pull parameters from")
    exit()
  else:
    sim_list.sort(key=os.path.getctime) # sort list by file creation time
    sim_script = sim_list[-1]
    period, kT, eta0 = find_params(sim_script)


  # extract remaining params from path
  def config_frompath(current_dir):
    required_patterns = {
      "radii" : r".*/r\d+(\.\d+|-\d+(\.\d+)?)?/.*/.*",
      "box size (L)" : r".*/(?:L\d+/.*/.*|\d+large)",
      #"interaction range (kappa)" : r'.*/kappa(\d+).*'
      "interaction range (kappa)" : r'(?:kappa|_delta)(\d+(?:\.\d+)?)' 
    }
    for pattern_name, pattern in required_patterns.items():
      if not re.search(pattern, current_dir):
        print(f"Error: Current path is missing the expected pattern for '{pattern_name}'")
        print(f"       Cannot extract all properties from path.")
        exit()
    print("Filepath seems to contain required parameters, extracting values...")

    # Extract R_C1 and R_C2 from current directory
    radii = current_dir.split('/r')[1].split('/')[0]
    first_r, second_r = map(float, radii.split('-'))
    R_C1 = first_r
    R_C2 = second_r
    if R_C1 == 0 or R_C2 == 0:
      if rank == 0 or comm.size == 1:
        print("ERROR: monomodal system detected: R_C1="+str(R_C1)+", R_C2="+str(R_C2))
        exit()

    # set colloid typeid
    #if 'DPD' in current_dir:
    if 'DPD' in param_path:
      colloid1_typeid = 1
      colloid2_typeid = 2
    #elif 'BD' in current_dir:
    elif 'BD' in param_path:
      colloid1_typeid = 0
      colloid2_typeid = 1
    else:
      print(f"ERROR: 'DPD' and 'BD' not found, cannot set colloid type for {param_path}") 
      exit()

    # Extract L_X from current directory
    #L_X = int(current_dir.split('/L')[1].split('/')[0])
    L_match = re.search(r'/L(\d+)', current_dir)
    L_X = int(L_match.group(1)) if L_match else None
    if L_X == None:
      Nlarge_match = re.search(r'/(\d+)large', current_dir)
      Nlarge = int(Nlarge_match.group(1)) if Nlarge_match else None
      if Nlarge == None:
        print('ERROR: no box size L or particle number provided')
        exit(1)
      else:
        comp_match = re.search(r'/(\d+)-(\d+)', current_dir)
        large_comp = int(comp_match.group(2)) if comp_match else None
        if large_comp == None:
          print('ERROR: no bimodal composition provided, cannot caluclate L from particle number')
          exit(1)
        else:
          phi_match = re.search(r'/phi(\d+)', current_dir)
          total_phi = float(phi_match.group(1)) if phi_match else None
          if total_phi == None:
            print('ERROR: no volume fraction in path, cannot calculate L from particle number')
            exit(1)
          else:
            phi_C2 = (large_comp/100) * (total_phi/100)
            L_X = (Nlarge * ((4/3)*np.pi*max(R_C1,R_C2)**3) / phi_C2)**(1./3.)

    kappa_match = re.search(r'.*/kappa(\d+).*', current_dir)
    kappa = int(kappa_match.group(1)) if kappa_match else None
    if kappa == None:
      delta_match = re.search(r'.*_delta(\d+(?:\.\d+)?)', current_dir)
      delta = float(delta_match.group(1)) if delta_match else None
      if delta != None:
        kappa = int(3/delta)
      elif delta == None:
        kappa = 30
        print("Warning: no Morse attraction range (kappa) set, use default kappa=30")

    return L_X, R_C1, R_C2, colloid1_typeid, colloid2_typeid, kappa

  current_dir = os.getcwd()
  if gsd_sim_path == '..':
    L_X, R_C1, R_C2, colloid1_typeid, colloid2_typeid, kappa = config_frompath(current_dir)
  else:
    L_X, R_C1, R_C2, colloid1_typeid, colloid2_typeid, kappa = config_frompath(gsd_sim_path)

  Lbox_shortest = L_X # the shortest side of the simulation box
  R_min = min(R_C1, R_C2)
  R_max = max(R_C1, R_C2)

  return filepath, data_outpath, period, kT, eta0, Lbox_shortest, R_C1, R_C2, R_min, R_max, colloid1_typeid, colloid2_typeid, kappa


########

def call_netx_variables(data_path_name=None, gsd_sim_path='..'):

  if gsd_sim_path == '..':
    print(f"[config_from_path] Current working directory: {os.getcwd()}\n")
    param_path = os.getcwd()
  else:
    print(f"[config_from_path] Simulation directory: {gsd_sim_path}\n")
    param_path = gsd_sim_path

  # select the last Gelation sim file one folder up
  file_list = []
  for filepath in glob.glob(f'{gsd_sim_path}/Gelation*.gsd'):
    if 'lastframe' not in filepath:
        file_list.append(filepath)
  if not file_list:
    print("ERROR: no GSD files found")
    exit()
  else:
    file_list.sort(key=os.path.getctime) # sort list by file creation time
    filepath = file_list[-1]

  # filepath to folder where data files will be created
  if data_path_name == None:
    data_outpath = 'data'
  else:
    data_outpath = data_path_name

  # extract remaining params from path
  def config_frompath(current_dir):
    required_patterns = {
      #"interaction range (kappa)" : r'.*/kappa(\d+).*'
      "interaction range (kappa)" : r'(?:kappa|_delta)(\d+(?:\.\d+)?)'
    }
    for pattern_name, pattern in required_patterns.items():
      if not re.search(pattern, current_dir):
        print(f"Error: Current path is missing the expected pattern for '{pattern_name}'")
        print(f"       Cannot extract all properties from path.")
        exit()
    #print("Filepath seems to contain required parameters, extracting values...")

    # set colloid typeid
    #if 'DPD' in current_dir:
    if 'DPD' in param_path:
      colloid1_typeid = 1
      colloid2_typeid = 2
    #elif 'BD' in current_dir:
    elif 'BD' in param_path:
      colloid1_typeid = 0
      colloid2_typeid = 1
    else:
      print(f"ERROR: 'DPD' and 'BD' not found, cannot set colloid type for {param_path}") 
      exit()

    kappa_match = re.search(r'.*/kappa(\d+).*', current_dir)
    kappa = int(kappa_match.group(1)) if kappa_match else None
    if kappa == None:
      delta_match = re.search(r'.*_delta(\d+(?:\.\d+)?)', current_dir)
      delta = float(delta_match.group(1)) if delta_match else None
      if delta != None:
        kappa = int(3/delta)
      elif delta == None:
        kappa = 30
        print("Warning: no Morse attraction range (kappa) set, use default kappa=30")

    return colloid1_typeid, colloid2_typeid, kappa

  current_dir = os.getcwd()
  if gsd_sim_path == '..':
    colloid1_typeid, colloid2_typeid, kappa = config_frompath(current_dir)
  else:
    colloid1_typeid, colloid2_typeid, kappa = config_frompath(gsd_sim_path)

  return filepath, data_outpath, colloid1_typeid, colloid2_typeid, kappa


########

def call_voronoi_variables(data_path_name=None, gsd_sim_path='..'):

  if gsd_sim_path == '..':
    print(f"[config_from_path] Current working directory: {os.getcwd()}\n")
    param_path = os.getcwd()
  else:
    print(f"[config_from_path] Simulation directory: {gsd_sim_path}\n")
    param_path = gsd_sim_path

  # select the last Gelation sim file one folder up
  file_list = []
  for filepath in glob.glob(f'{gsd_sim_path}/Gelation*.gsd'):
    if 'lastframe' not in filepath:
        file_list.append(filepath)
  if not file_list:
    print("ERROR: no GSD files found")
    exit()
  else:
    file_list.sort(key=os.path.getctime) # sort list by file creation time
    filepath = file_list[-1]

  # filepath to folder where data files will be created
  if data_path_name == None:
    data_outpath = 'data'
  else:
    data_outpath = data_path_name

  # extract remaining params from path
  def config_frompath(current_dir):
    # set colloid typeid
    #if 'DPD' in current_dir:
    if 'DPD' in param_path:
      colloid1_typeid = 1
      colloid2_typeid = 2
    #elif 'BD' in current_dir:
    elif 'BD' in param_path:
      colloid1_typeid = 0
      colloid2_typeid = 1
    else:
      print(f"ERROR: 'DPD' and 'BD' not found, cannot set colloid type for {param_path}")
      exit()

    return colloid1_typeid, colloid2_typeid

  current_dir = os.getcwd()
  if gsd_sim_path == '..':
    colloid1_typeid, colloid2_typeid = config_frompath(current_dir)
  else:
    colloid1_typeid, colloid2_typeid = config_frompath(gsd_sim_path)

  return filepath, data_outpath, colloid1_typeid, colloid2_typeid


########

def call_shear_variables(data_path_name=None, gsd_sim_path='../..', shear_sim_path='..'):

  if (gsd_sim_path == '../..') or (shear_sim_path == '..'):
    print(f"[config_from_path] Current working directory: {os.getcwd()}\n")
    print(f"[config_from_path] Shear directory: {shear_sim_path}\n")
    print(f"[config_from_path] Gelation directory: {gsd_sim_path}\n")
  else:
    print(f"[config_from_path] Shear directory: {shear_sim_path}\n")
    print(f"[config_from_path] Gelation directory: {gsd_sim_path}\n")

  # select the last Shear sim file one folder up
  shear_file_list = []
  for filepath in glob.glob(f'{shear_sim_path}/Shear*.gsd'):
      #if not filepath.startswith('{shear_sim_path}/Shear_Colloids'):
      shear_file_list.append(filepath)
  if not shear_file_list:
    print("ERROR: no Shear GSD files found")
    exit()
  else:
    shear_file_list.sort(key=os.path.getctime) # sort list by file creation time
    shear_filepath = shear_file_list[-1]

  # select the last Gelation sim file one folder up
  gel_file_list = []
  for gel_path in glob.glob(f'{gsd_sim_path}/Gelation*.gsd'):
    if 'lastframe' not in gel_path:
        gel_file_list.append(gel_path)
  if not gel_file_list:
    print("ERROR: no GSD files found")
    exit()
  else:
    gel_file_list.sort(key=os.path.getctime) # sort list by file creation time
    gel_filepath = gel_file_list[-1]

  # filepath to folder where data files will be created
  if data_path_name == None:
    data_outpath = 'data'
  else:
    data_outpath = data_path_name


  ## simulation specific parameters
  # extract remaining params from path
  def config_frompath(current_dir):
    required_patterns = {
      "radii" : r".*/r\d+(\.\d+|-\d+(\.\d+)?)?/.*/.*",
      "box size (L)" : r".*/L\d+/.*/.*",
      "shear parameters" : r".*shear-.*"
    }
    for pattern_name, pattern in required_patterns.items():
      if not re.search(pattern, current_dir):
        print(f"Error: Current path is missing the expected pattern for '{pattern_name}'")
        print(f"       Cannot extract all properties from path.")
        exit()
    print("Filepath seems to contain required parameters, extracting values...")

    # Extract R_C1 and R_C2 from current directory
    radii = current_dir.split('/r')[1].split('/')[0]
    first_r, second_r = map(float, radii.split('-'))
    R_C1 = first_r
    R_C2 = second_r
    if R_C1 == 0 or R_C2 == 0:
      if rank == 0 or comm.size == 1:
        print("ERROR: monomodal system detected: R_C1="+str(R_C1)+", R_C2="+str(R_C2))
        exit()

    # set colloid typeid
    #if 'DPD' in os.path.basename(shear_filepath):
    if 'DPD' in current_dir:
      colloid1_typeid = 1
      colloid2_typeid = 2
    #elif 'BD' in os.path.basename(shear_filepath):
    elif 'BD' in current_dir:
      colloid1_typeid = 0
      colloid2_typeid = 1
    else:
      print(f"ERROR: 'DPD' and 'BD' not found, cannot set colloid type for {current_dir}")
      exit()

    # Extract L_X from current directory
    #L_X = int(current_dir.split('/L')[1].split('/')[0])
    L_match = re.search(r'/L(\d+)', current_dir)
    L_X = int(L_match.group(1)) if L_match else None
    if L_X == None:
      Nlarge_match = re.search(r'/(\d+)large', current_dir)
      Nlarge = int(Nlarge_match.group(1)) if Nlarge_match else None
      if Nlarge == None:
        print('ERROR: no box size L or particle number provided')
        exit(1)
      else:
        comp_match = re.search(r'/(\d+)-(\d+)', current_dir)
        large_comp = int(comp_match.group(2)) if comp_match else None
        if large_comp == None:
          print('ERROR: no bimodal composition provided, cannot caluclate L from particle number')
          exit(1)
        else:
          phi_match = re.search(r'/phi(\d+)', current_dir)
          total_phi = float(phi_match.group(1)) if phi_match else None
          if total_phi == None:
            print('ERROR: no volume fraction in path, cannot calculate L from particle number')
            exit(1)
          else:
            phi_C2 = (large_comp/100) * (total_phi/100)
            L_X = (Nlarge * ((4/3)*np.pi*max(R_C1,R_C2)**3) / phi_C2)**(1./3.)

    ## Extract shear parameters (style + shear rate or amplitude/frequency)
    # check if any shear information exists
    shear_match = re.search(r'shear-.*', current_dir)
    if not shear_match:
      if rank == 0 or comm.size == 1:
        print("Error: no shear information found in the path")
        exit()

    # extract shear_style
    shear_style_match = re.search(r'shear-(.*?)-', current_dir)
    if not shear_style_match:
      if rank == 0 or comm.size == 1:
        print("Error: shear style could not be extracted from shear infromation")
        exit()

    shear_style = shear_style_match.group(1)

    # extract values based on shear_style
    if shear_style == 'const':
      Amp = 0.0
      omega = 0.0
      const_match = re.search(r'shear-const-([\d\.]+)-', current_dir)
      if const_match:
        shear_rate_const = float(const_match.group(1))
      else:
        if rank == 0 or comm.size == 1:
          print("Error: constant shear dectected, but no shear rate is provided")
          exit()

    elif shear_style == 'cosin':
      shear_rate_const = 0.0
      cosin_match = re.search(r'shear-cosin-A([\d\.]+)-omega([\d\.]+)-', current_dir)
      if cosin_match:
        Amp = float(cosin_match.group(1))
        omega = float(cosin_match.group(2))
      else:
        if rank == 0 or comm.size == 1:
          print("Error: oscillatory (cosin) shear dectected, but amplitude and frequency are not provided")
          exit()

    else:
      if rank == 0 or comm.size == 1:
        print(f"Error: unknown shear style detected, cannot extract shear parameters from '{shear_style}'")
        exit()

    return L_X, R_C1, R_C2, colloid1_typeid, colloid2_typeid, shear_style, shear_rate_const, Amp, omega

  current_dir = os.getcwd()
  if shear_sim_path == '..':
     L_X, R_C1, R_C2, colloid1_typeid, colloid2_typeid, shear_style, shear_rate_const, Amp, omega = config_frompath(current_dir)
  else:
     L_X, R_C1, R_C2, colloid1_typeid, colloid2_typeid, shear_style, shear_rate_const, Amp, omega = config_frompath(shear_sim_path)
  L_Y = L_X

  ## simulation specific parameters
  # extract parameters from most recent simulation script:
  def find_params(simfile):
    # Check if the file exists
    if not os.path.exists(simfile):
      print("ERROR: no sim files found matching", simfile)
      exit()

    # Open the file and find parameters
    with open(simfile, 'r') as f:
      kT = eta0 = frames_per_strain = dt_Integration = theta = None
      param_dict = {"kT":kT,
                    "eta0":eta0,
                    "frames_per_strain":frames_per_strain,
                    "dt_Integration":dt_Integration,
                    "theta":theta
      }
      lines = f.readlines()
      for line in lines:
        for param, value in param_dict.items():
          if (param + " = " in line) and ("#" + param not in line):
            param_dict[param] = line.split('=')[1].split('#')[0].strip()
        if line.startswith("L_Y") and ("L_Y = L_X" not in line):
          print(f"Warning: L_Y value is set in {simfile}, double-check Lbox_shortest")
        if line.startswith("L_Z") and ("L_Z = L_X" not in line):
          print(f"Warning: L_Z value is set in {simfile}, double-check Lbox_shortest")

      # check for missing values
      for param, value in param_dict.items():
        if value == None:
          print(f"ERROR: no {param} value in {simfile}")
      if None in param_dict.values():
        quit()

      return (float(param_dict["kT"]),
              float(param_dict["eta0"]),
              float(param_dict["frames_per_strain"]),
              float(param_dict["dt_Integration"]),
              float(param_dict["theta"])
      )

  # use find_params to get params
  sim_list = []
  for simpath in glob.glob(f'{shear_sim_path}/sim-*.py'):
    if 'lastframe' not in simpath and 'convert-to' not in simpath and 'scaleup' not in simpath:
      sim_list.append(simpath)
  if not sim_list:
    print("ERROR: no sim scripts found to pull parameters from")
    exit()
  else:
    sim_list.sort(key=os.path.getctime) # sort list by file creation time
    sim_script = sim_list[-1]
    kT, eta0, frames_per_strain, dt_Integration, theta = find_params(sim_script)

  # recalcualte additional shear values
  if shear_style == 'const':
    shear_rate = shear_rate_const
    #vinf = shear_rate * L_Y
    t_ramp = theta/shear_rate
    #file_tag = shear_style+'-SR'+str(shear_rate)+'-'+str(n_strains)+'cycles'

  else:
    omega_timesteps = omega*dt_Integration # frequency (1/N_timesteps)
    #vinf = Amp * omega * L_Y
    shear_rate = Amp * omega
    t_ramp = Amp/shear_rate
    #file_tag = shear_style+'-SR'+str(round(shear_rate,2))+"-A"+str(round(Amp,3))+"-w"+str(round(omega,1))+'-'+str(n_strains)+'cycles'

  Iter_cycle = int(t_ramp / dt_Integration)
  period = int(Iter_cycle/frames_per_strain)
  MDtime = period*dt_Integration


  R_min = min(R_C1, R_C2)
  R_max = max(R_C1, R_C2)

  return shear_filepath, gel_filepath, data_outpath, shear_style, shear_rate_const, shear_rate, Amp, omega, kT, eta0, frames_per_strain, dt_Integration, theta, Iter_cycle, period, MDtime, colloid1_typeid, colloid2_typeid
