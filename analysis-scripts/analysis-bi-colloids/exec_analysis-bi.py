#!/usr/bin/env python3
import os
import sys
import subprocess
from datetime import datetime

# file/dir names
analysis_target_dir = 'analysis-bi-colloids'
data_outpath = 'data'

# Path to your centralized analysis scripts and compiled .so
#ANALYSIS_HOME = "/projects/props/Rob/colloids/gel-scripts/hoomd4.2.1-scripts/config-from-path/2-analysis/analysis-bi-colloids"
ANALYSIS_HOME = str(sys.argv[1]) 

# List of analysis scripts to run (must accept --input and --output args)
SCRIPTS = [
    "config_from_path.py",
    "sim-analysis-bi-colloids.py",
    "sim-networkCSV-bi-colloids.py",
    "sim-voronoi-bi-colloids.py"
]


# Add analysis home to PYTHONPATH
if ANALYSIS_HOME not in sys.path:
    sys.path.insert(0, ANALYSIS_HOME)

def write_version_info(output_dir):
    """Write a file with timestamp and script version info."""
    version_log = os.path.join(output_dir, "analysis_version.txt")
    with open(version_log, "w") as f:
        f.write(f"Analysis run on: {datetime.now()}\n")
        f.write(f"Scripts from: {ANALYSIS_HOME}\n")
        try:
            git_hash = subprocess.check_output(
                ["git", "-C", ANALYSIS_HOME, "rev-parse", "HEAD"], text=True
            ).strip()
            f.write(f"Git commit: {git_hash}\n")
        except:
            f.write("Git commit: [unknown or not a git repo]\n")

def main():
    # Start in current directory (data directory)
    data_dir = os.getcwd()

    # Create and move into 'analysis/' subdirectory
    analysis_dir = os.path.join(data_dir, analysis_target_dir)
    os.makedirs(analysis_dir, exist_ok=True)
    os.chdir(analysis_dir)  # change working directory

    # Define data input and output paths relative to parent data directory
    output_dir = analysis_dir  # save everything in 'analysis/'
    data_output_dir = os.path.join(analysis_dir, data_outpath)

    # Log version info
    write_version_info(output_dir)

    # Run each analysis script from ANALYSIS_HOME
    for script in SCRIPTS:
        full_path = os.path.join(ANALYSIS_HOME, script)
        if not os.path.exists(full_path):
            print(f"ERROR, Script not found: {full_path}")
            continue

        if script != 'config_from_path.py':
          print(f"▶ Running: {script}")
          cmd = ["python3", "-u", full_path] #, "--input", data_file, "--output", data_output_dir]
          subprocess.run(cmd, check=True)
          #print("Done:", script)

        print("----")

if __name__ == "__main__":
    main()

