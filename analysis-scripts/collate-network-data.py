import numpy as np
import pandas as pd
import os


size_ratios = ['1-0', '1-2', '0-2']
volume_fractions = ['05', '10', '15', '20', '25', '30']
small_comps = [100, 75, 50, 25, 0]
compositions = ['100-0', '75-25', '50-50', '25-75', '0-100']
attraction_strengths = ['03', '05', '07', '09', '12']
depletion_scaling = ['bysize']
attraction_range = ['60']
L = 70

net_data_columns = ['network-type','frame','cut_off','n_components',
                    'lcc_size','ncolloids','f_lcc',
                    'lcc_frac_type1','bond_frac_11','bond_frac_12','bond_frac_22',
                    'core2_size','core2_frac_type1','n_dangling','dangling_frac_type1',
                    'lcc_span','box_span','percolation',
                    'chi_meanclustersize','mean_finite_size',
                    'lcc_extent_x','lcc_extent_y','lcc_extent_z','max_axis_span','n_axes_spanning',
                    'avg_degree',
                    #'avg_clustering_coeff','avg_sq_clustering_coeff'
                    ]


net_data_list = []
for size in size_ratios:
  radii = size.split('-')
  R_C1 = float(radii[0])
  R_C2 = float(radii[1])
  for phi in volume_fractions:
      for comp_1 in small_comps:
          comp_2 = 100 - comp_1
          comp = f'{comp_1}-{comp_2}'
          # skip invalid fake mixtures for monodisperse case
          #print(size, comp, R_C1, R_C2, f'0-{int(R_C2)}', (size == f'0-{int(R_C2)}'))
          if size != '1-0' and comp == '100-0':
              #print('...skip bimodal size or mono_large size for mono_small comp')
              continue
          if size == '1-0' and comp != '100-0':
              #print('---skip mono_small size for bimodal comp or mono_large comp')
              continue
          if size != f'0-{int(R_C2)}' and comp == '0-100':
              #print('++++skip bimodal size or mono_small size for mono_large comp')
              continue
          if size == f'0-{int(R_C2)}' and comp != '0-100':
              #print('////skip mono_large size for bimodal comp or mono_small comp')
              continue

          path_size = size
          if size == f'0-{int(R_C2)}' and comp == '0-100':
              path_size = f'{int(R_C2)}-0'
              #print('reformat size for mono_large size and mono_large comp')
              if R_C2 != float(int(R_C2)):
                  print('ERROR: large monomodal data processing script cannot accept float particle size: {R_C2}')
                  continue
          #print(size, comp)
          for D0 in attraction_strengths:
              for scaling in depletion_scaling:
                  for kappa in attraction_range:
                      if (path_size == '1-0') and (comp == '100-0'):
                         data_path = (f'/projects/props/Rob/colloids/bimodal/r{path_size}/DPD/L{L}/mono/seedNone/phi{phi}'
                            f'/potential-morse-bysize/{D0}kT/kappa{kappa}/')
                         data_outpath=f'{data_path}analysis-DPD/data/'
                         #print(data_path)
                      elif (path_size == f'{int(R_C2)}-0') and (comp == '0-100'):
                         data_path = (f'/projects/props/Rob/colloids/bimodal/r{path_size}/DPD/L{L}/phi{phi}'
                            f'/potential-morse-bysize/{D0}kT/kappa{kappa}/')
                         data_outpath=f'{data_path}analysis-DPD/data/'
                         #print(data_outpath)
                      elif (comp != '100-0') and (comp != '0-100'):
                         data_path = (f'/projects/props/Rob/colloids/bimodal/r{path_size}'
                            f'/DPD/L70/poly0.0/seedNone/phi{phi}/{comp}/potential-morse-{scaling}/'
                            f'{D0}kT/kappa{kappa}/')
                         data_outpath=f'{data_path}analysis-bi-colloids/data/'
                      else:
                         continue

                      net_path = f"{data_outpath}networkx-allframes.csv"
                      if os.path.exists(net_path) == False:
                          print(f"No networkx-allframes.csv file: {data_outpath}")
                          continue

                      Zavg_path = f"{data_outpath}Zavg.csv"
                      if os.path.exists(Zavg_path) == False:
                          print(f"No Z-counts.csv file: {data_outpath}")
                          continue

                      msd_path = f"{data_outpath}msd.txt"
                      if os.path.exists(msd_path) == False:
                          print(f"No msd.txt file: {data_outpath}")
                          continue

                      df_net = pd.read_csv(net_path)

                      if 'lcc_frac_type1' not in df_net.keys(): 
                          #'bond_frac_11', 'bond_frac_12', 'bond_frac_22', 'core2_size', 'core2_frac_type1', 'n_dangling', 'dangling_frac_type1'] not in df_net.keys():
                          print(f'NO BAKCBONE DATA: phi={phi}, {comp}, {D0}kT')
                          if '0' in path_size:
                              if str(int(R_C1)) in path_size:
                                  df_net['lcc_frac_type1'] = 1
                                  df_net['bond_frac_11'] = 1
                                  df_net['bond_frac_12'] = 0
                                  df_net['bond_frac_22'] = 0
                                  df_net['core2_size'] = None
                                  df_net['core2_frac_type1'] = 1
                                  df_net['n_dangling'] = None
                                  df_net['dangling_frac_type1'] = 1
                              if str(int(R_C1)) in path_size:
                                  df_net['lcc_frac_type1'] = 0
                                  df_net['bond_frac_11'] = 0
                                  df_net['bond_frac_12'] = 0
                                  df_net['bond_frac_22'] = 1
                                  df_net['core2_size'] = None
                                  df_net['core2_frac_type1'] = 0
                                  df_net['n_dangling'] = None
                                  df_net['dangling_frac_type1'] = 0

                      df = df_net[net_data_columns].copy()

                      if '0' not in path_size:
                          net_c1_path = f"{data_outpath}networkx-allframes-c1only.csv"
                          if os.path.exists(net_c1_path) == False:
                              print(f"No networkx-allframes-c1only.csv file: {data_outpath}")
                              continue

                          net_c2_path = f"{data_outpath}networkx-allframes-c2only.csv"
                          if os.path.exists(net_c1_path) == False:
                              print(f"No networkx-allframes-c2only.csv file: {data_outpath}")
                              continue

                          df_net_c1 = pd.read_csv(net_c1_path)
                          df_c1 = df_net_c1[net_data_columns].copy()

                          df_net_c2 = pd.read_csv(net_c2_path)
                          df_c2 = df_net_c2[net_data_columns].copy()

                          df = pd.concat([df, df_c1, df_c2], ignore_index=True) 

                      meta = {
                          "phi": float(phi)/100,
                          "comp": comp,
                          "D0": int(D0),
                          "kappa": int(kappa),
                      }


                      """Zavg + msd, one row per frame. msd's missing frame 0 -> 0."""
                      zavg = pd.read_csv(Zavg_path).rename(columns={"simframe": "frame"})
                      zavg["frame"] = zavg["frame"].astype(int)
                      zavg = zavg.rename(
                          columns={c: ("Zavg_" + c) for c in zavg.columns if c != "frame"}
                      )

                      if '0' in path_size:
                          if str(int(R_C1)) in path_size:
                              zavg['Zavg_Z_11'] = zavg['Zavg_Z_any']
                              zavg['Zavg_Z_12'] = 0
                              zavg['Zavg_Z_21'] = 0
                              zavg['Zavg_Z_22'] = 0
                              zavg['Zavg_Z_1any'] = zavg['Zavg_Z_any']
                              zavg['Zavg_Z_2any'] = 0
                          elif str(int(R_C2)) in path_size:
                              zavg['Zavg_Z_11'] = 0
                              zavg['Zavg_Z_12'] = 0
                              zavg['Zavg_Z_21'] = 0
                              zavg['Zavg_Z_22'] = zavg['Zavg_Z_any'] 
                              zavg['Zavg_Z_1any'] = 0
                              zavg['Zavg_Z_2any'] = zavg['Zavg_Z_any']
                          else:
                              print(f'ERROR: path_size contains 0 but not {R_C1} or {R_C2}: {path_size}')

                      msd = pd.read_csv(msd_path, sep=r"\s+")
                      msd["frame"] = msd["frame"].astype(int)
                      msd = msd.drop(columns=[c for c in msd.columns if "solv" in c.lower()])

                      if '0' in path_size:
                          if str(int(R_C1)) in path_size:
                              msd['msd_colloid_1'] = msd['msd_colloid']
                              msd['msdc1_smplstddev'] = msd['msd_coll_smplstddev']
                              msd['msd_colloid_2'] = 0
                              msd['msdc2_smplstddev'] = 0
                          elif str(int(R_C2)) in path_size:
                              msd['msd_colloid_1'] = 0 
                              msd['msdc1_smplstddev'] = 0 
                              msd['msd_colloid_2'] = msd['msd_colloid'] 
                              msd['msdc2_smplstddev'] = msd['msd_coll_smplstddev'] 
                          else:
                              print(f'ERROR: path_size contains 0 but not {R_C1} or {R_C2}: {path_size}')

                      msd_value_cols = [c for c in msd.columns if c != "frame"]

                      frames = np.union1d(zavg["frame"].to_numpy(), msd["frame"].to_numpy())
                      out = pd.DataFrame({"frame": np.sort(frames)})
                      out = out.merge(zavg, on="frame", how="left")
                      out = out.merge(msd, on="frame", how="left")

                      # MSD(t=0) = 0: fill any missing msd values (the frame-0 gap) with 0.
                      out[msd_value_cols] = out[msd_value_cols].fillna(0.0)

                      df_combined = df.merge(out, on="frame", how="left")

                      """Insert metadata columns at the front, in order."""
                      for i, (k, v) in enumerate(meta.items()):
                          df_combined.insert(i, k, v)

                      net_data_list.append(df_combined)

          print(f"data processed for phi={phi}, {comp}, {D0}kT")

outfile = f"phase_diagram_net_data.parquet" 
pd.concat(net_data_list, ignore_index=True).to_parquet(outfile, index=False)
print(f"{outfile} created from all available data")
