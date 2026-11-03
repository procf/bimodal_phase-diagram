"""
Collect one run's analysis outputs into 3 tables and save each as its own parquet, sorted by total volume fraction:

  1. per-frame          : Zavg.csv + msd.txt        -> <prefix>_perframe.parquet
                          (msd missing frame 0 -> filled with 0, since MSD(t=0)=0)
  2. per-pore           : poresize(-all).csv         -> <prefix>_perpore.parquet
  3. per-particle-frame : voronoi-volumes.txt + Z-counts.csv (joined on
                          frame + particle, where Z-counts.colloidID == particle)
                                                     -> <prefix>_perparticle.parquet

Metadata (dataset_id, system_type, phi, comp, D0, lbB2_ct) is included on every
row.
Monodisperse vs bimodal is auto-detected from the Z-counts columns.
"""

import os
import numpy as np
import pandas as pd


# ---------------------------------------------------------------- helpers
def _read_ws(path):
    return pd.read_csv(path, sep=r"\s+")


def _exists(path):
    return path is not None and os.path.isfile(path)


def _detect_bimodal(data_dir):
    zc = pd.read_csv(os.path.join(data_dir, "Z-counts.csv"), nrows=1)
    return "Z_11" in zc.columns


def _stamp(df, meta):
    """Insert metadata columns at the front, in order."""
    for i, (k, v) in enumerate(meta.items()):
        df.insert(i, k, v)
    return df


# ---------------------------------------------------------------- loaders
def load_per_frame(data_dir, meta):
    """Zavg + msd, one row per frame. msd's missing frame 0 -> 0."""
    p = lambda n: os.path.join(data_dir, n)

    zavg = pd.read_csv(p("Zavg.csv")).rename(columns={"simframe": "frame"})
    zavg["frame"] = zavg["frame"].astype(int)
    zavg = zavg.rename(
        columns={c: ("Zavg_" + c) for c in zavg.columns if c != "frame"}
    )

    msd = _read_ws(p("msd.txt"))
    msd["frame"] = msd["frame"].astype(int)
    msd = msd.drop(columns=[c for c in msd.columns if "solv" in c.lower()])
    msd_value_cols = [c for c in msd.columns if c != "frame"]

    frames = np.union1d(zavg["frame"].to_numpy(), msd["frame"].to_numpy())
    out = pd.DataFrame({"frame": np.sort(frames)})
    out = out.merge(zavg, on="frame", how="left")
    out = out.merge(msd, on="frame", how="left")

    # MSD(t=0) = 0: fill any missing msd values (the frame-0 gap) with 0.
    out[msd_value_cols] = out[msd_value_cols].fillna(0.0)

    return _stamp(out, meta)


def load_per_pore(data_dir, meta):
    """poresize, one row per detected void. Empty table if the file is absent."""
    p = lambda n: os.path.join(data_dir, n)
    bimodal = meta["system_type"] == "bimodal"
    pore_path = p("poresize-allc.csv") if bimodal else p("poresize.csv")

    if not _exists(pore_path):
        cols = ["frame", #"probe_posx", "probe_posy", "probe_posz",
                "voidcenter_x", "voidcenter_y", "voidcenter_z",
                #"porediameter_T", 
                #"porediameter_G"]
                "void-diameter_G"]
        print(f"  [warn] no pore file at {pore_path} -> empty per-pore table")
        return _stamp(pd.DataFrame(columns=cols), meta)

    pore = pd.read_csv(pore_path)
    pore["frame"] = pore["frame"].astype(int)
    pore = pore.rename(columns={"porediameter_G": "void-diameter_G"})
    pore = pore.drop(columns=["probe_posx", "probe_posy","probe_posz","porediameter_T"])
    pore = pore.reset_index(drop=True)
    return _stamp(pore, meta)



# decomposed-Z schema shared by bimodal and (synthesized) monodisperse tables
_Z_DECOMP = ["Z_11", "Z_12", "Z_21", "Z_22", "Z_1any", "Z_2any"]
 
 
def load_per_particle(data_dir, meta):
    """
    voronoi-volumes + Z-counts joined on (frame, particle).
    Z-counts.colloidID is renamed to `particle` to match voronoi.
    One row per particle per frame.
 
    Coordination schema (uniform across mono and bimodal):
      Z_any                      : total coordination number (was Z_any / Z)
      Z_11,Z_12,Z_21,Z_22        : type-resolved pair coordinations
      Z_1any,Z_2any              : per-species total coordinations
 
    Monodisperse runs carry no type-resolved counts, so they are mapped into the
    shared schema from `meta['comp']`:
      comp == '100-0' -> all particles are species 1: typeid=1,
                         Z_11 = Z_1any = Z_any, others 0
      comp == '0-100' -> all particles are species 2: typeid=2,
                         Z_22 = Z_2any = Z_any, others 0
    """
    p = lambda n: os.path.join(data_dir, n)
    bimodal = meta["system_type"] == "bimodal"
 
    vor = _read_ws(p("voronoi-volumes.txt"))

    # Standardize the volume column label across the three file layouts:
    #   0-100 (mono large): frame particle radius volume
    #   100-0 (mono small): frame particle radius total-volume free-volume
    #   bimodal           : frame particle typeid radius total-volume free-volume
    # In every case the per-particle *total* Voronoi cell volume is what we keep
    # as voro_total (the bare "volume" column in the 0-100 layout is that total).
    if "total-volume" in vor.columns:
        vor = vor.rename(columns={"total-volume": "voro_total"})
    elif "volume" in vor.columns:
        vor = vor.rename(columns={"volume": "voro_total"})
    else:
        raise KeyError(
            f"voronoi-volumes.txt in {data_dir} has neither 'total-volume' or "
            f"'volume'; columns = {list(vor.columns)}"
        )

    vor["frame"] = vor["frame"].astype(int)
    vor["particle"] = vor["particle"].astype(int)

    # free volme is incorrectly calculated for monodisperse and not in large mono
    # ...use original voro only
    vor = vor.drop(columns=[c for c in ["free-volume", "voro_free"] if c in vor.columns])
 
    zc = pd.read_csv(p("Z-counts.csv")).rename(columns={"colloidID": "particle"})
    zc["frame"] = zc["frame"].astype(int)
    zc["particle"] = zc["particle"].astype(int)
    # typeid appears in both; keep voronoi's and drop the duplicate from Z-counts
    zc = zc.drop(columns=[c for c in ["typeid"] if c in zc.columns])
 
    mono_typeid = None  # set below for monodisperse, used to override typeid post-merge
    #if bimodal:
    #    # canonical total-coordination name
    #    zc = zc.rename(columns={"Z_any": "Z_all"})
    #else:
    if not bimodal:
        zc = zc.rename(columns={"Z": "Z_any"})
        # synthesize the decomposed schema from composition
        for c in _Z_DECOMP:
            zc[c] = 0
        comp = str(meta.get("comp", "")).strip()
        if comp == "100-0":
            mono_typeid = 1
            zc["Z_11"] = zc["Z_any"]
            zc["Z_1any"] = zc["Z_any"]
        elif comp == "0-100":
            mono_typeid = 2
            zc["Z_22"] = zc["Z_any"]
            zc["Z_2any"] = zc["Z_any"]
        else:
            print(f"  [warn] monodisperse comp={comp!r} is not '100-0' or "
                  f"'0-100'; decomposed Z left as 0 and typeid not overridden")
 
    merged = vor.merge(
        zc, on=["frame", "particle"], how="outer", indicator=True
    )
    n_unmatched = int((merged["_merge"] != "both").sum())
    if n_unmatched:
        print(f"  [warn] {n_unmatched} particle-frame rows not in both "
              f"voronoi and Z-counts (kept via outer join)")
    merged = merged.drop(columns=["_merge"])
 
    # monodisperse: force typeid to match the single species present
    if mono_typeid is not None:
        merged["typeid"] = mono_typeid
 
    merged = merged.sort_values(["frame", "particle"]).reset_index(drop=True)

    return _stamp(merged, meta)

# ---------------------------------------------------------------- driver
#def load_all_split(data_dir, phi, comp, D0, dataset_id=None):
def load_all_split(data_dir, phi, comp, D0, lnB2_ct=np.nan, dataset_id=None):
    if dataset_id is None:
        dataset_id = os.path.basename(os.path.normpath(data_dir))
    meta = {
        #"dataset_id": dataset_id,
        "system_type": "bimodal" if _detect_bimodal(data_dir) else "monodisperse",
        "phi": float(phi)/100,
        "comp": comp,
        "D0": int(D0),
        "kappa": int(kappa),
        "lnB2_round": lnB2_ct,
    }
    return {
        #"perframe": load_per_frame(data_dir, meta),
        "perpore": load_per_pore(data_dir, meta),
        "perparticle": load_per_particle(data_dir, meta),
    }


def save_all_split(tables, out_prefix, engine="pyarrow"):
    """Write the three tables. out_prefix e.g. 'data/run01' -> run01_perframe.parquet ..."""
    paths = {}
    for name, df in tables.items():
        path = f"{out_prefix}_{name}.parquet"
        df.to_parquet(path, engine=engine, index=False)
        paths[name] = path
    return paths




#############################
## SET PARAMETERS
#############################

size_ratios = ['1-0', '1-2', '0-2']
#volume_fractions = ['05', '10', '15', '20', '25', '30']
#small_comps = [100, 75, 50, 25, 0]
#compositions = ['100-0', '75-25', '50-50', '25-75', '0-100']
#attraction_strengths = ['03', '05', '07', '09', '12']
volume_fractions = ['20']
small_comps = [100, 25]
compositions = ['100-0', '25-75']
attraction_strengths = ['05', '09']
depletion_scaling = ['bysize']
attraction_range = ['60']
L = 70

## Manually sort by volume fraction
dataset = {"10":[{"size":'1-0', "phi":'10', "comp":'100-0', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-2', "phi":'10', "comp":'25-75', "D0":'03', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-0', "phi":'10', "comp":'100-0', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-2', "phi":'10', "comp":'50-50', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-0', "phi":'10', "comp":'100-0', "D0":'12', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'10', "comp":'75-25', "D0":'07', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'10', "comp":'25-75', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        {"size":'0-2', "phi":'10', "comp":'0-100', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        ], '15':[
        {"size":'1-0', "phi":'15', "comp":'100-0', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-2', "phi":'15', "comp":'25-75', "D0":'03', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-0', "phi":'15', "comp":'100-0', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-2', "phi":'15', "comp":'50-50', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-0', "phi":'15', "comp":'100-0', "D0":'12', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'15', "comp":'75-25', "D0":'07', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'15', "comp":'25-75', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        {"size":'0-2', "phi":'15', "comp":'0-100', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        ], "20":[
        {"size":'1-0', "phi":'20', "comp":'100-0', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-2', "phi":'20', "comp":'25-75', "D0":'03', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-0', "phi":'20', "comp":'100-0', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-2', "phi":'20', "comp":'50-50', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-0', "phi":'20', "comp":'100-0', "D0":'12', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'20', "comp":'75-25', "D0":'07', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'20', "comp":'25-75', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        {"size":'0-2', "phi":'20', "comp":'0-100', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        ],'25':[
        {"size":'1-0', "phi":'25', "comp":'100-0', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-2', "phi":'25', "comp":'25-75', "D0":'03', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-0', "phi":'25', "comp":'100-0', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-2', "phi":'25', "comp":'50-50', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-0', "phi":'25', "comp":'100-0', "D0":'12', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'25', "comp":'75-25', "D0":'07', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'25', "comp":'25-75', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        {"size":'0-2', "phi":'25', "comp":'0-100', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        ], "30":[
        {"size":'1-0', "phi":'30', "comp":'100-0', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-2', "phi":'30', "comp":'25-75', "D0":'03', "scaling":'bysize', "kappa":'60', "lnB2_ct": 6.0,},
        {"size":'1-0', "phi":'30', "comp":'100-0', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-2', "phi":'30', "comp":'50-50', "D0":'05', "scaling":'bysize', "kappa":'60', "lnB2_ct": 13.0,},
        {"size":'1-0', "phi":'30', "comp":'100-0', "D0":'12', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'30', "comp":'75-25', "D0":'07', "scaling":'bysize', "kappa":'60', "lnB2_ct": 19.0,},
        {"size":'1-2', "phi":'30', "comp":'25-75', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        {"size":'0-2', "phi":'30', "comp":'0-100', "D0":'09', "scaling":'bysize', "kappa":'60', "lnB2_ct": 30.0,},
        ]}


for p, phi_label in enumerate(['10','15','20','25','30']):
    phi_dataset = dataset[phi_label]    
    void_df_list = []
    particle_df_list = []
    for params in phi_dataset:
                      ##for params in dataset:
                      size = params['size']
                      radii = size.split('-')
                      R_C1 = float(radii[0])
                      R_C2 = float(radii[1])
                      phi = params["phi"]
                      comp = params["comp"]
                      path_size = size
                      D0 = params["D0"]
                      scaling = params["scaling"]
                      kappa = params["kappa"]
                      lnB2 = params["lnB2_ct"]

                      if (path_size == '1-0') and (comp == '100-0'):
                         data_path = (f'/projects/props/Rob/colloids/bimodal/r{path_size}/DPD/L{L}/mono/seedNone/phi{phi}'
                            f'/potential-morse-bysize/{D0}kT/kappa{kappa}/')
                         data_outpath=f'{data_path}analysis-DPD/data/'
                         #print(data_path)

                      elif (path_size == f'0-{int(R_C2)}') and (comp == '0-100'):
                         data_path = (f'/projects/props/Rob/colloids/bimodal/r{int(R_C2)}-0/DPD/L{L}/phi{phi}'
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

                      tabs = load_all_split(data_outpath, phi, comp, D0, lnB2_ct=lnB2)
                      pervoid_df = tabs['perpore']
                      void_df_list.append(pervoid_df.drop(columns='system_type'))

                      perparticle_df = tabs['perparticle']
                      particle_df_list.append(perparticle_df.drop(columns='system_type'))
                      print(f" - loaded ln(B2)={lnB2}: phi={phi}%, {comp}, {D0}kT")

    void_df = pd.concat(void_df_list, ignore_index=True)
    dist_df = pd.concat(particle_df_list, ignore_index=True)

    void_df.to_parquet(f"voidsize_phi{phi}.parquet", engine="pyarrow", index=False)
    dist_df.to_parquet(f"distributions_phi{phi}.parquet", engine="pyarrow", index=False)


    print(f'all data for phi={phi_label} saved to voidsize_phi{phi_label}.parquet and distributions_phi{phi_label}.parquet')
