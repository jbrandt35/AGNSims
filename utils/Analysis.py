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

perturber_separation_range_to_analyze = (5.0, 5.0)
BBH_separation_range_to_analyze = (0.3, 0.3)

def run_finished(directory):
    return os.path.exists(os.path.join(directory, "outcome.json"))


def retrieve_data(run_directory, data_name, data_location, source):
    if source == "data_file":
        with pd.HDFStore(os.path.join(run_directory, "result", "data.h5"), complevel = 5, complib = "zlib") as data_file:
            data = data_file[data_location][data_name].tolist()
    elif source == "outcome_file":
        try:
            with open(os.path.join(run_directory, "outcome.json")) as outcome_file:
                data = json.load(outcome_file)[data_name]
        except FileNotFoundError:
            #Return the result of an unfinished run
            if data_name == "Result":
                data = "Run Didn't Finish"
            #If trying to access data that doesn't exist, throw an error
            else:
                raise FileNotFoundError[f"The file {os.path.join(run_directory, 'outcome.json')} couldn't be opened"]
    return data


def run_ended_in_merger(directory):
    return "Collision Encountered" in retrieve_data(directory, "Result", None, "outcome_file")


#Returns whether BBH_1 and BBH_2 shared a bounded orbit at the end of the simulation
def run_ended_with_bound_binary(directory):
    # Preferably, pull it from outcome.json
    try:
        final_eccentricity = retrieve_data(directory, "Final Binary Eccentricity", None, "outcome_file")
    # If not available there, take the last recorded eccentricity
    except KeyError:
        final_eccentricity = retrieve_data(directory, "binary eccentricity", "Misc", "data_file")[-1]
    return bool(final_eccentricity < 1)


def get_data_values(data_name, i, data_location = None, source = "data_file", only_include_finished = True, only_include_merged = False, exclude_merged = False, exclude_unbound_binaries = False):
    values = []

    for BBH_separation in [float(directory.split("_")[-1]) for directory in os.listdir() if "BBH_separation" in directory]:
        if not (BBH_separation_range_to_analyze[0] <= BBH_separation <= BBH_separation_range_to_analyze[1]):
            continue
        for perturber_separation in [float(directory.split("_")[-1]) for directory in os.listdir(f"BBH_separation_{BBH_separation}") if "perturber_separation" in directory]:
            if not (perturber_separation_range_to_analyze[0] <= perturber_separation <= perturber_separation_range_to_analyze[1]):
                continue
            for run_number in os.listdir(os.path.join(f"BBH_separation_{BBH_separation}", f"perturber_separation_{perturber_separation}")):
                directory = os.path.join(f"BBH_separation_{BBH_separation}", f"perturber_separation_{perturber_separation}", run_number)

                run_is_unfinished = not run_finished(directory)

                #Skip this run if only want finished runs and this one isn't finished
                if only_include_finished and run_is_unfinished:
                    continue

                #If the run finished, evaluate the ending conditions
                if not run_is_unfinished:

                    run_ended_merged = run_ended_in_merger(directory)
                    binary_ended_unbound = not run_ended_with_bound_binary(directory)

                    #Skip this run if only want merged runs and this one isn't merged
                    if only_include_merged and not run_ended_merged:
                        continue
                    #Skip this run if only want un-merged and this one ended merged
                    if exclude_merged and run_ended_merged:
                        continue
                    #Skip this run if only want bound binaries and this one ended unbound
                    if exclude_unbound_binaries and binary_ended_unbound:
                        continue

                data = retrieve_data(directory, data_name, data_location, source)

                if i is not None:
                    value = data[i]
                else:
                    value = data

                values.append(value)

    return np.array(values)


def get_ith_vectors(i, data_name, data_location):
    x_components = get_data_values(f"{data_name}_x", i, data_location = data_location, only_include_merged = True)
    y_components = get_data_values(f"{data_name}_y", i, data_location = data_location, only_include_merged = True)
    z_components = get_data_values(f"{data_name}_z", i, data_location = data_location, only_include_merged = True)
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


def create_bar_labels(data, bars):
    labels = []
    for bar in bars:
        bar_height, bar_width = bar.get_height(), bar.get_width()
        if bar_height == 0:
            labels.append("")
        else:
            labels.append(str(int(np.round(bar_height * bar_width * len(data)))))
    return labels


def create_histogram(*data, title, x_label, y_label):
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = cycle(prop_cycle.by_key()["color"])
    for datum in data:
        color = next(colors)
        counts, edges, bars = plt.hist(np.array(datum), color = color, density = True)
        plt.bar_label(bars, labels = create_bar_labels(datum, bars))
        sns.kdeplot(datum, color = color)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.savefig(f"../plots/{title}.jpg", bbox_inches = "tight")
    plt.close()


def generate_outcome_piechart():

    outcomes = get_data_values("Result", None, source = "outcome_file", only_include_finished = False)

    aggregate_outcome_data = {}
    for outcome in outcomes:
        if ":" in outcome:
            outcome = outcome.split(":")[0]
        if outcome not in list(aggregate_outcome_data.keys()):
            aggregate_outcome_data[outcome] = 1
        else:
            aggregate_outcome_data[outcome] += 1

    plt.pie(list(aggregate_outcome_data.values()), labels = list(aggregate_outcome_data.keys()), normalize = True, autopct = "%.2f%%")
    plt.title(f"Distribution of Simulation Outcomes (N = {len(outcomes)})")
    plt.savefig("../plots/outcome_piechart.jpg", bbox_inches = "tight")
    plt.close()

generate_outcome_piechart()

#Create Spin Histograms
create_histogram(get_ith_inclinations(-1, "S", "Positions/binary", "L_Binary", "Positions/binary"), title = "Final Inclination of Spin to Binary-SMBH L", x_label = "Inclination [deg]", y_label = "Distribution")
create_histogram(get_ith_inclinations(-1, "S", "Positions/binary", "L_BBH2", "Positions/binary"), title = "Final Inclination of Spin to $m1_b$-$m1_a$ L", x_label = "Inclination [deg]", y_label = "Distribution")
create_histogram(get_ith_inclinations(-1, "L_Binary", "Positions/binary", "L_BBH2", "Positions/binary"), title = "Final Inclination of Binary-SMBH L to $m1_b$-$m1_a$ L", x_label = "Inclination [deg]", y_label = "Distribution")

initial_t_GW_values = get_data_values("t_GW", 0, data_location = "Positions/binary", exclude_merged = True, exclude_unbound_binaries = True)
print(initial_t_GW_values)







