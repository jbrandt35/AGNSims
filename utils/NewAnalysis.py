import os
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd
import seaborn as sns
from itertools import cycle
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
            data = data_file[data_location][data_name].tolist()
    elif source == "outcome_file":
        with open(os.path.join(run_directory, "outcome.json")) as outcome_file:
            data = json.load(outcome_file)[data_name]

    return data

def get_ith_values(data_name, i, data_location = None, source = "data_file"):

    ith_values = []

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
                ith_value = data[i]
                ith_values.append(ith_value)

    return np.array(ith_values)

def get_ith_vectors(i, data_name, data_location):
    x_components = get_ith_values(f"{data_name}_x", i, data_location = data_location)
    y_components = get_ith_values(f"{data_name}_y", i, data_location = data_location)
    z_components = get_ith_values(f"{data_name}_z", i, data_location = data_location)
    vectors = np.stack([x_components, y_components, z_components], axis = 1)
    return vectors

def get_inclination(v_1, v_2):
    return np.degrees(np.arccos(np.dot(v_1, v_2)))

def get_ith_inclinations(i, data_name_1, data_location_1, data_name_2, data_location_2):
    vectors_1 = get_ith_vectors(i, data_name_1, data_location_1)
    vectors_2 = get_ith_vectors(i, data_name_2, data_location_2)

    inclinations = []

    for j in range(len(vectors_1)):
        v_1 = np.array(vectors_1[j])
        v_2 = np.array(vectors_2[j])
        inclination = get_inclination(v_1, v_2)
        inclinations.append(inclination)

    return inclinations


def create_histogram(*data, **hist_kwargs):
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = cycle(prop_cycle.by_key()["color"])
    for datum in data:
        color = next(colors)
        plt.hist(datum, density = True, color = color)
        sns.kdeplot(datum, color = color)
    plt.title(hist_kwargs["title"])
    plt.xlabel(hist_kwargs["x_label"])
    plt.ylabel(hist_kwargs["y_label"])
    plt.legend()
    plt.savefig(f"../plots/{hist_kwargs['title']}.jpg", bbox_inches = "tight")
    plt.close()


#create_histogram(get_final_values("x", data_location = "Positions/binary"), get_final_values("y", data_location = "Positions/binary"), title = "Title", x_label = "x", y_label = "y")









