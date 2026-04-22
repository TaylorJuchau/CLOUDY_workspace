import os
import shutil
import numpy as np
import pandas as pd
import subprocess
import re
import json
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
import re
import sys
import os
import glob
from Functions import *
# -------------------------------
# Paths
# -------------------------------
work_dir = "/project/galaxies/tjuchau/projects/CLOUDY_scripts/Galaxy_class_project2/work/"
output_dir = "/project/galaxies/tjuchau/projects/CLOUDY_scripts/Galaxy_class_project2/outputs/"

os.makedirs(work_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# -------------------------------
# Parameter grid
# -------------------------------
log_nH_grid = np.linspace(1, 5, 10)
log_U_grid = np.linspace(-4, -1, 10)
temp_grid = np.linspace(35000, 50000, 10)
# -------------------------------
# Write CLOUDY input file
# -------------------------------
def write_cloudy_input(model_name, log_nH, log_U, temp):
    filename = os.path.join(work_dir, model_name + ".in")

    with open(filename, 'w') as f:
        f.write(f"title Orion model log_nH={log_nH:.2f} log_U={log_U:.2f}\n")
        f.write(f"blackbody {temp} K\n")
        f.write(f"hden {log_nH}\n")
        f.write(f"ionization parameter {log_U}\n")
        f.write("abundances GASS\n")
        f.write("sphere\n")
        f.write("stop temperature 4000 K\n")
        f.write("stop column density 22\n")
        f.write("iterate to convergence\n")
        f.write(f'save overview "{model_name}.txt"\n')

def write_AGN_cloudy(model_name, log_nH, log_U):
    filename = os.path.join(work_dir, model_name + ".in")


    with open(filename, 'w') as f:
        f.write(f"title AGN log_nH={log_nH:.2f} log_U={log_U:.2f}\n")
        f.write(f"hden {log_nH}\n")
        f.write(f'table AGN\n')
        f.write(f'stop efrac 0.01\n')

        f.write(f"ionization parameter {log_U}\n")
        f.write("abundances GASS\n")
        f.write("sphere\n")
        f.write("stop temperature 4000 K\n")
        f.write("stop column density 22\n")
        f.write("iterate to convergence\n")
        f.write(f'save overview "{model_name}.txt"\n')
# -------------------------------
# Run CLOUDY
# -------------------------------
def run_cloudy(model_name, cloudy_path):
    result = subprocess.run(
        [cloudy_path, "-r", model_name],
        cwd=work_dir,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"CLOUDY failed for {model_name}")

# -------------------------------
# Move outputs
# -------------------------------
def move_outputs(model_name):
    for ext in [".out"]:
        src = os.path.join(work_dir, model_name + ext)
        dst = os.path.join(output_dir, model_name + ext)
        
        if os.path.exists(src):
            shutil.move(src, dst)
        else:
            print(f"WARNING: {src} not found")

# -------------------------------
# Helper
# -------------------------------
def try_float(x):
    try:
        return float(x)
    except:
        return None

def parse_target_line(target):
    """
    "O  3 5006.84A" -> ("O  3", 5006.84e-10)
    """
    parts = target.split()
    species = f"{parts[0]}  {parts[1]}"
    wavelength = float(parts[2].replace("A", "")) * 1e-10
    return species, wavelength

# -------------------------------
# Parse CLOUDY output
# -------------------------------
def parse_emission_lines(filename, target_lines, wl_tol=5e-10):
    """
    Extract specific emission lines from CLOUDY .out file.

    Args:
        filename: path to .out file
        target_lines: list like ["H  1 4861.33A", "O  3 5006.84A"]
        wl_tol: wavelength tolerance in meters

    Returns:
        dict mapping line name -> flux
    """

    # Preprocess targets
    targets = []
    for t in target_lines:
        species, wl = parse_target_line(t)
        targets.append((t, species, wl))

    results = {t: None for t in target_lines}

    iteration_pattern = re.compile(r'Iteration\s+(\d+)\s+of\s+(\d+)')
    final_iteration = False
    line_type = None

    with open(filename, 'r') as file:
        for line in file:

            # Wait for final iteration
            if not final_iteration:
                match = iteration_pattern.search(line)
                if match:
                    n, m = match.groups()
                    if n == m:
                        final_iteration = True
                continue

            # Track section
            if 'Intrinsic line intensities' in line:
                line_type = 'intrinsic'
                continue
            if 'Emergent line intensities' in line:
                line_type = 'emergent'
                continue

            # Only use emergent lines
            if line_type != 'emergent':
                continue

            # Try matching each target
            for label, species, target_wl in targets:

                if species in line:
                    start = line.index(species)
                    parts = line[start + len(species):].split()

                    if len(parts) < 3:
                        continue

                    wl_str = parts[0]
                    flux = try_float(parts[1])

                    if flux is None:
                        continue

                    # Convert wavelength
                    if 'A' in wl_str:
                        wl = float(wl_str.replace('A', '')) * 1e-10
                    elif 'm' in wl_str:
                        wl = float(wl_str.replace('m', '')) * 1e-6
                    else:
                        continue

                    # Match wavelength
                    if abs(wl - target_wl) < wl_tol:
                        results[label] = flux

    return results

# -------------------------------
# Compute ratios
# -------------------------------
def compute_ratios(lines):
    return {
        "OIII_Hb": lines["O  3 5006.84A"] / lines["H  1 4861.33A"],
        "NII_Ha": lines["N  2 6583.45A"] / lines["H  1 6562.80A"],
        "SII_Ha": (
            lines["S  2 6716.44A"] + lines["S  2 6730.82A"]
        ) / lines["H  1 6562.80A"],
        "SII_SII": lines["S  2 6716.44A"] / lines["S  2 6730.82A"],
        'OI_Ha': lines["O  1 6300.30A"]/lines["H  1 6562.80A"],
        "HeII_Hb": lines["He 2 4685.68A"]/lines["H  1 4861.33A"]
    }

# -------------------------------
# Main grid loop
# -------------------------------
def run_grid(cloudy_path, model):

    target_lines = [
        "H  1 4861.33A",
        "H  1 6562.80A",
        "O  3 5006.84A",
        "N  2 6583.45A",
        "S  2 6716.44A",
        "S  2 6730.82A",
        "O  1 6300.30A",
        "He 2 4685.68A"
    ]

    files = []
    for log_nH in log_nH_grid:
        for log_U in log_U_grid:
            for temp in temp_grid:
                
                model_name = f"{model}model_nH{log_nH:.2f}_U{log_U:.2f}_temp{temp:.2f}".replace("-", "m")
                if glob.glob(f'/project/galaxies/tjuchau/projects/CLOUDY_scripts/Galaxy_class_project2/outputs/{model_name}.out') != []:
                    files.append(glob.glob(f'/project/galaxies/tjuchau/projects/CLOUDY_scripts/Galaxy_class_project2/outputs/{model_name}.out')[0])
                    print(f'Model {model_name} already exists, skipping...')
                    continue
                try:
                    if model == 'HII':
                        write_cloudy_input(model_name, log_nH, log_U, temp)
                    elif model == 'AGN':
                        write_AGN_cloudy(model_name, log_nH, log_U)

                    run_cloudy(model_name, cloudy_path)
                    move_outputs(model_name)

                    
                except Exception as e:
                    print(f"FAILED: {model_name} -> {e}")
    return files

def analyze_grid(model_type):
    directory = '/project/galaxies/tjuchau/projects/CLOUDY_scripts/Galaxy_class_project2/outputs/'
    results = []
    files = glob.glob(directory+f"{model_type}*")
    for i, file in enumerate(files):
        trial_info = {}
        lines = parse_emission_lines(outfile, target_lines)
        # Skip incomplete models
        if any(v is None for v in lines.values()):
            print(f"Missing lines in {model_name}")
            continue
        ratios = compute_ratios(lines)

        trail_info["log_nH"] = log_nH
        trial_info["log_U"] = log_U
        trial_info['temp'] = temp
        trial_info["lines"] = lines
        trial_info['ratios'] = ratios
        results.append(trial_info)
    results.to_csv(f"cloudy_{model_type}_grid.csv", index=False)
    

    return pd.DataFrame(results)
# -------------------------------
# Run everything
# -------------------------------
cloudy_path = "/project/galaxies/tjuchau/software/Cloudy/c25.00/source/cloudy.exe"
model_type = "AGN"
run_grid(cloudy_path, model_type)
model_type = 'HII'
run_grid(cloudy_path, model_type)

def find_best(model_type):
    best = np.inf
    chi2=[]
    targets = {}
    if model_type == "AGN":
        targets['NII_Ha'] = 0.438
        targets['SII_Ha'] = 0.374
        targets['OI_Ha'] = 0.0868
        targets['HeII_Hb'] = 0.38
    elif model_type == "HII":
        targets['OIII_Hb'] = 2.820
        targets['NII_Ha'] = 0.19
        targets['SII_Ha'] = 0.034
        targets['SII_SII'] = 0.683
    results = analyze_grid(model_type)
    for i, trial in enumerate(results):
        if model_type == "AGN":
            c = ((np.log10(trial['ratios']["NII_Ha"])  - np.log10(targets['NII_Ha']))**2 +
            (np.log10(trial['ratios']["SII_Ha"])  - np.log10(targets['SII_Ha']))**2 +
            (np.log10(trial['ratios']["OI_Ha"])   - np.log10(targets['OI_Ha']))**2 +
            (np.log10(trial['ratios']["HeII_Hb"]) - np.log10(targets['HeII_Hb']))**2)
            chi2.append(c)
            if c < best:
                best = c
                best_trial = trial
        elif model_type == "HII":
            c = ((np.log10(trial['ratios']["OIII_Hb"])  - np.log10(targets['OIII_Hb']))**2 +
            (np.log10(trial['ratios']["NII_Ha"])  - np.log10(targets['NII_Ha']))**2 +
            (np.log10(trial['ratios']["SII_Ha"])   - np.log10(targets['SII_Ha']))**2 +
            (np.log10(trial['ratios']["SII_SII"]) - np.log10(targets['SII_SII']))**2)
            chi2.append(c)
            if c < best:
                best = c
                best_trial = trial
    return best_trial
best_AGN = find_best('AGN')
best_HII = find_best('HII')
