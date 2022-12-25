import os
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

os.system("mkdir -p ../plots")
os.chdir("../runs/")

BBH_separation_range_to_analyze = (0.1, 0.4)
perturber_separation_range_to_analyze = (2.5, 3.5)

def run_finished(directory):
    return os.path.exists(os.path.join(directory, "outcome.json"))

def retrieve_data(run_directory, data_name, data_location, source):
    if source == "data_file":
        with pd.HDFStore(os.path.join(run_directory, "result", "data.h5"), complevel = 5, complib = "zlib") as data_file:
            data = data_file[data_name]
    elif source == "outcome_file":
        with open(os.path.join(run_directory, "outcome.json")) as outcome_file:
            data = json.load(outcome_file)[data_location][data_name]

    return data

def get_final_values(data_name, data_location, source = "data_file"):

    final_values = []

    for BBH_separation in [float(directory.split("_")[-1]) for directory in os.listdir() if "BBH_separation" in directory]:
        if not (BBH_separation_range_to_analyze[0] <= BBH_separation <= BBH_separation_range_to_analyze[1]):
            continue
        for perturber_separation in [float(directory.split("_")[-1]) for directory in os.listdir(f"BBH_separation_{BBH_separation}") if "perturber_separation" in directory]:
            if not (perturber_separation_range_to_analyze[0] <= perturber_separation <= perturber_separation_range_to_analyze[1]):
                continue
            for run_number in os.listdir(os.path.join(f"BBH_separation_{BBH_separation}", f"perturber_separation_{perturber_separation}")):
                directory = os.path.join(f"BBH_separation_{BBH_separation}", f"perturber_separation_{perturber_separation}", run_number)
                if not run_finished(directory):
                    continue

                data = retrieve_data(directory, data_name, data_location, source)
                final_value = data[-1]
                final_values.append(final_value)

    return np.array(final_values)









