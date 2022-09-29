import os
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd
import seaborn as sns
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

run_dir = "/storage/home/hhive1/jbrandt35/data/AGN_sims/gr_heatmap/runs"
plot_dir = "/".join(run_dir.split("/")[:-1]) + "/plots"

perturber_separation_range = (2.5, 10, 0.25)
BBH_separation_range = (0.1, 0.5, 0.02)

perturber_separation_range = np.array(list(range(int(100 * perturber_separation_range[0]), int(100 * (perturber_separation_range[1] + perturber_separation_range[2])), int(100 * perturber_separation_range[2])))) / 100
BBH_separation_range = np.array(list(range(int(100 * BBH_separation_range[0]), int(100 * (BBH_separation_range[1] + BBH_separation_range[2])), int(100 * BBH_separation_range[2])))) / 100

os.system("mkdir -p ../plots")
#########################Create Data Frames to Hold Data##########################################

empty_data = np.empty((len(perturber_separation_range), len(BBH_separation_range)))
empty_data[:] = np.nan

data = pd.DataFrame(data = empty_data)
data.columns = BBH_separation_range
data.index = perturber_separation_range

bounded_data = data.copy(deep = True)
t_GW_data = data.copy(deep = True)
relative_t_GW_data = data.copy(deep = True)
min_BBH_distance = data.copy(deep = True)
perturber_distance_data = data.copy(deep = True)
merged_data = data.copy(deep = True)
ejected_data = data.copy(deep = True)

#########################  Get Data  ##########################################

bounded_data_dict = {}
t_GW_data_dict = {}
relative_t_GW_data_dict = {}
BBH_separation_data_dict = {}
perturber_BBH_distance_data_dict = {}
merged_data_dict = {}
ejected_data_dict = {}

pairs = []

total_outcomes = {"Reached $10^5$ Periods": 0, "Binary Unbound ($e>1$)": 0, "Didn't Finish": 0, "$t_{\mathregular{GW}}$ < Period": 0, "Event Horizons Crossed": 0, "$m_2$ merged with binary": 0}
all_min_t_GWs = []


for BBH_separation in BBH_separation_range:
    for perturber_separation in perturber_separation_range:

        pair = (BBH_separation, perturber_separation)

        pairs.append((float(pair[0]), float(pair[1])))

        directory = os.path.join(run_dir, f"BBH_separation_{BBH_separation}", f"perturber_separation_{perturber_separation}")
        run_numbers = np.arange(1, 101, 1)

        bound_outcomes = []
        t_GWs = []
        relative_t_GWs = []
        perturber_distances = []
        BBH_distances = []
        merged_outcomes = []
        ejected_outcomes = []

        for i in run_numbers:

            outcome_file = os.path.join(directory, str(i), "outcome.json")

            if os.path.exists(outcome_file):

                outcome = json.load(open(outcome_file))

                ended_bound = True

            #  Check Bounded-ness  #

                if "Unbound" in outcome["Result"]:
                    bound_outcomes.append(1)
                    merged_outcomes.append(0)
                    ejected_outcomes.append(1)
                    total_outcomes["Binary Unbound ($e>1$)"] += 1
                    ended_bound = False
                elif "The simulation ended with the binary bound and no mergers." in outcome["Result"]:
                    bound_outcomes.append(0)
                    merged_outcomes.append(0)
                    ejected_outcomes.append(0)
                    total_outcomes["Reached $10^5$ Periods"] += 1
                elif "Collision Encountered" in outcome["Result"]:
                    bound_outcomes.append(1)
                    merged_outcomes.append(1)
                    ejected_outcomes.append(0)
                    if "when the sum of event horizon radii was" in outcome["Result"]:
                        print(f"A collision occured! For BBH separation {BBH_separation} with m2 separation {perturber_separation}, run {i}: m1a and m1b crossed event horizons.")
                        total_outcomes["Event Horizons Crossed"] += 1
                    if "t_GW was" in outcome["Result"]:
                        print(
                            f"A collision occured! For BBH separation {BBH_separation} with m2 separation {perturber_separation}, run {i}: t_GW fell below the period.")
                        total_outcomes["$t_{\mathregular{GW}}$ < Period"] += 1
                    if "The distance between m1_a and m2" in outcome["Result"]:
                        print(
                            f"A collision occured! For BBH separation {BBH_separation} with m2 separation {perturber_separation}, run {i}: m1a and m2 merged.")
                        total_outcomes["$m_2$ merged with binary"] += 1
                    if "he distance between m1_b and m2" in outcome["Result"]:
                        print(
                            f"A collision occured! For BBH separation {BBH_separation} with m2 separation {perturber_separation}, run {i}: m1b and m2 merged.")
                        total_outcomes["$m_2$ merged with binary"] += 1


                #  Store t_GW  #
                if ended_bound:
                    min_t_GW = outcome["Minimum t_GW"]
                    t_GWs.append(min_t_GW)
                    all_min_t_GWs.append(min_t_GW)

                    min_relative_t_GW = outcome["Minimum relative t_GW"]
                    relative_t_GWs.append(min_relative_t_GW)

                min_distance_binary_perturber = outcome["Minimum Distance Between Binary COM and Perturber"]
                perturber_distances.append(min_distance_binary_perturber)

                min_binary_dist = outcome["Minimum Distance Between BBHs"]
                BBH_distances.append(min_binary_dist)

            else:
                print(f"For BBH separation {BBH_separation} with m2 separation {perturber_separation}, run {i} didn't finish!")
                total_outcomes["Didn't Finish"] += 1

        bounded_data_dict[pair] = bound_outcomes
        t_GW_data_dict[pair] = t_GWs
        relative_t_GW_data_dict[pair] = relative_t_GWs
        BBH_separation_data_dict[pair] = BBH_distances
        perturber_BBH_distance_data_dict[pair] = perturber_distances
        merged_data_dict[pair] = merged_outcomes
        ejected_data_dict[pair] = ejected_outcomes


#########################   Format Data  ##########################################

def format_average(dict):
    formated_dict = {}
    for key in dict.keys():
        average = np.average(dict[key])

        new_key = (float(key[0]), float(key[1]))
        formated_dict[new_key] = average

    return formated_dict

def format_prop(dict):
    formated_data_dict = {}
    for key in dict.keys():
        lst = dict[key]

        if len(lst) > 0:
            mergers = lst.count(1)
            prop = mergers / len(lst)
        else:
            prop = None

        new_key = (float(key[0]), float(key[1]))
        formated_data_dict[new_key] = prop

    return formated_data_dict

formated_bounded_data_dict = format_prop(bounded_data_dict)
formated_merged_data_dict = format_prop(merged_data_dict)
formated_ejected_data_dict = format_prop(ejected_data_dict)

t_GW_data_dict = format_average(t_GW_data_dict)
relative_t_GW_data_dict = format_average(relative_t_GW_data_dict)
BBH_separation_data_dict = format_average(BBH_separation_data_dict)
perturber_BBH_distance_data_dict = format_average(perturber_BBH_distance_data_dict)

#########################  Fill In Data   ##########################################

for (BBH_separation, perturber_separation) in pairs:
    bounded_data.loc[perturber_separation, BBH_separation] = formated_bounded_data_dict[(BBH_separation, perturber_separation)]
    t_GW_data.loc[perturber_separation, BBH_separation] = t_GW_data_dict[(BBH_separation, perturber_separation)]
    relative_t_GW_data.loc[perturber_separation, BBH_separation] = relative_t_GW_data_dict[(BBH_separation, perturber_separation)]
    min_BBH_distance.loc[perturber_separation, BBH_separation] = BBH_separation_data_dict[(BBH_separation, perturber_separation)]
    perturber_distance_data.loc[perturber_separation, BBH_separation] = perturber_BBH_distance_data_dict[(BBH_separation, perturber_separation)]
    merged_data.loc[perturber_separation, BBH_separation] = formated_merged_data_dict[(BBH_separation, perturber_separation)]
    ejected_data.loc[perturber_separation, BBH_separation] = formated_ejected_data_dict[(BBH_separation, perturber_separation)]

#########################  Create Plots  ##########################################

cmap = sns.dark_palette("#69d", reverse = True, as_cmap = True)

########################## Bounded Heatmaps ############################

fig, (merged_plot, ejected_plot, bounded_plot) = plt.subplots(1, 3, figsize = (18, 5), sharey = True)

fig.suptitle("$10^5$ Periods")

bounded_heatmap = sns.heatmap(bounded_data, cmap = cmap, cbar = True, vmin = 0, vmax = 1, ax = bounded_plot)
bounded_heatmap.invert_yaxis()
bounded_plot.set_xlabel("Binary Separation [$R_{\mathregular{hill, m_1}}$]")
bounded_plot.set_title("Proportion of Unstable Binaries")

merged_heatmap = sns.heatmap(merged_data, cmap = cmap, cbar = False, ax = merged_plot, vmin = 0, vmax = 1)
merged_heatmap.invert_yaxis()
merged_plot.set_xlabel("Binary Separation [$R_{\mathregular{hill, m_1}}$]")
merged_plot.set_ylabel("Perturber Distance from Binary [$R_{\mathregular{hill,} m_1-m_2}}$]")
merged_plot.set_title("Proportion of Merged Binaries")

ejected_heatmap = sns.heatmap(ejected_data, cmap = cmap, cbar = False, vmin = 0, vmax = 1, ax = ejected_plot)
ejected_heatmap.invert_yaxis()
ejected_plot.set_xlabel("Binary Separation [$R_{\mathregular{hill, m_1}}$]")
ejected_plot.set_title("Proportion of Binary Ejections")

plt.savefig(f"{plot_dir}/bounded_heatmaps.jpg", bbox_inches = "tight")
plt.close()

########################## t_GW Plots ############################

from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize = (20, 12))
gs = GridSpec(nrows = 2, ncols = 4)
gs.update(wspace = 0.25)

t_GW_plot = fig.add_subplot(gs[0, :2])
t_GW_heatmap = sns.heatmap(t_GW_data, cmap = cmap, cbar = True, norm = LogNorm(), ax = t_GW_plot)
t_GW_heatmap.collections[0].colorbar.ax.set_title("[yr]")
t_GW_heatmap.invert_yaxis()
t_GW_plot.set_xlabel("Binary Separation [$R_{\mathregular{hill, m_1}}$]")
t_GW_plot.set_ylabel("Perturber Distance from Binary [$R_{\mathregular{hill,} m_1-m_2}}$]")
t_GW_plot.set_title("Minimum $t_{\mathregular{GW}}$")

relative_t_GW_plot = fig.add_subplot(gs[0, 2:])
relative_t_GW_heatmap = sns.heatmap(relative_t_GW_data, cmap = cmap, cbar = True, norm = LogNorm(), ax = relative_t_GW_plot)
relative_t_GW_heatmap.invert_yaxis()
relative_t_GW_plot.set_xlabel("Binary Separation [$R_{\mathregular{hill, m_1}}$]")
relative_t_GW_plot.set_ylabel("Perturber Distance from Binary [$R_{\mathregular{hill,} m_1-m_2}}$]")
relative_t_GW_plot.set_title("Minimum $t_{\mathregular{GW}} / T_{\mathregular{Binary}}$")

t_GW_histogram = fig.add_subplot(gs[1, 1:3])
t_GW_histogram.hist(all_min_t_GWs, bins = 150, align = "mid")
t_GW_histogram.set_xscale("log")
t_GW_histogram.set_title("Distribution of Minimum $t_{\mathregular{GW}}$")
t_GW_histogram.set_xlabel("Minimum $t_{\mathregular{GW}}$")
t_GW_histogram.set_ylabel("Count")

plt.savefig(f"{plot_dir}/t_GW_plots.jpg", bbox_inches = "tight")
plt.close()

########################## Distance Heatmaps ############################

fig, (BBH_distance_plot, perturber_distance_plot) = plt.subplots(1, 2, figsize = (12, 5), sharey = True)

BBH_heatmap = sns.heatmap(min_BBH_distance, cmap = cmap, cbar = True, ax = BBH_distance_plot)
BBH_heatmap.invert_yaxis()
BBH_heatmap.collections[0].colorbar.ax.set_title("[AU]")
BBH_distance_plot.set_ylabel("Perturber Distance from Binary [$R_{\mathregular{hill,} m_1-m_2}}$]")
BBH_distance_plot.set_xlabel("Binary Separation [$R_{\mathregular{hill, m_1}}$]")
BBH_distance_plot.set_title("Minimum Distance Between BBHs")

perturber_heatmap = sns.heatmap(perturber_distance_data, cmap = cmap, cbar = True, ax = perturber_distance_plot)
perturber_heatmap.invert_yaxis()
perturber_heatmap.collections[0].colorbar.ax.set_title("[AU]")
perturber_distance_plot.set_xlabel("Binary Separation [$R_{\mathregular{hill, m_1}}$]")
perturber_distance_plot.set_title("Minimum Distance Between Binary COM & $m_2$")

plt.savefig(f"{plot_dir}/distance_heatmaps.jpg", bbox_inches = "tight")
plt.close()

########################## Outcome Pie Chart  ############################

outcomes = list(total_outcomes.items())
outcomes_for_pie_labels = [outcome for outcome in outcomes if outcome[1] != 0]
outcome_labels_for_pie, counts_for_pie = map(list, zip(*outcomes_for_pie_labels))
outcome_labels_for_legend, counts_for_legend = map(list, zip(*outcomes))
plt.pie(counts_for_pie, labels = outcome_labels_for_pie, normalize = True, autopct = "%.2f%%")
plt.gca().set_prop_cycle(None)
plt.title("Distribution of Simulation Outcomes")
plt.legend(plt.pie(counts_for_legend)[0], outcome_labels_for_legend, loc = "upper left", bbox_to_anchor = (-0.5, 1))
plt.savefig(f"{plot_dir}/outcome_pie.jpg", bbox_inches = "tight")
plt.close()
