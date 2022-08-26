import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np
from numpy.linalg import norm as mag
import os


data_objects = ["BBH_1", "BBH_2", "SMBH", "perturber", "binary"]


def get_position_data(particle):
    data = pd.read_hdf("result/data.h5", key = f"Positions/{particle}")
    return data


def get_misc_data():
    data = pd.read_hdf("result/data.h5", key = "Misc")
    return data


def get_xyz(data_frame):
    return np.array(data_frame[["x", "y", "z"]].to_numpy().tolist())


def get_relative_xyz(data_frame_target, data_frame_reference):
    target_positions = get_xyz(data_frame_target)
    reference_positions = get_xyz(data_frame_reference)
    return target_positions - reference_positions


def extract_coordinate(positions, coordinate):
    if coordinate == "x":
        coordinate = 0
    elif coordinate == "y":
        coordinate = 1
    elif coordinate == "z":
        coordinate = 2
    coordinates = [position[coordinate] for position in positions]
    return np.array(coordinates)


def split_coordinates(data, mode = "x,y"):
    if mode == "x,y":
        return extract_coordinate(get_xyz(data), "x"), extract_coordinate(get_xyz(data), "y")


def get_last_time(df):
    return df["time"].iloc[-1]


def generate_data():
    d = {"Misc": get_misc_data()}
    for name in data_objects:
        df = get_position_data(name)
        try:
            df["time"] = d["Misc"]["time"].values
        except ValueError:
            df["time"] = d["Misc"]["time"].values[:-1]
        d[name] = df
    return d


data_dict = generate_data()

def get_approp_t_end():
    end_times = []
    for df in data_dict.values():
        end_times.append(get_last_time(df))
    return min(end_times)


def time_extract_data(t_start = 0, t_end = get_approp_t_end()):
    for name, df in data_dict.items():
        data_dict[name] = df[(df["time"] >= t_start) & (df["time"] <= t_end)].reset_index(drop = True)

time_extract_data()

SMBH_x, SMBH_y = split_coordinates(data_dict["SMBH"])
BBH_1_x, BBH_1_y = split_coordinates(data_dict["BBH_1"])
BBH_2_x, BBH_2_y = split_coordinates(data_dict["BBH_2"])
perturber_x, perturber_y = split_coordinates(data_dict["perturber"])
BBH_1_barycentric_x, BBH_1_barycentric_y = extract_coordinate(get_relative_xyz(data_dict["BBH_1"], data_dict["binary"]), "x"), extract_coordinate(get_relative_xyz(data_dict["BBH_1"], data_dict["binary"]), "y")
BBH_2_barycentric_x, BBH_2_barycentric_y = extract_coordinate(get_relative_xyz(data_dict["BBH_2"], data_dict["binary"]), "x"), extract_coordinate(get_relative_xyz(data_dict["BBH_2"], data_dict["binary"]), "y")


def distance(p1, p2):
    p1, p2 = get_xyz(p1), get_xyz(p2)
    distances = []
    for i, _ in enumerate(p1):
        distances.append(mag(p1[i] - p2[i]))
    return np.array(distances)


def construct_plots():

    os.system("mkdir -p plots")

    plt.plot(data_dict["Misc"]["time"], data_dict["Misc"]["time-step"])
    plt.xlabel("Time [yr]")
    plt.ylabel("dt [yr]")
    plt.title("Time-step over Time")
    plt.yscale("log", base = 10)
    plt.savefig("plots/time_vs_dt.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(data_dict["BBH_1"]["time"], distance(data_dict["BBH_1"], data_dict["BBH_2"]))
    plt.xlabel("Time [yr]")
    plt.ylabel("Distance [AU]")
    plt.title("Distance Between BHs in Binary over Time")
    plt.savefig("plots/time_vs_binary_distance.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(data_dict["binary"]["time"], data_dict["Misc"]["binary eccentricity"])
    plt.xlabel("Time [yr]")
    plt.ylabel("Eccentricity")
    plt.title("Eccentricity of the Binary in the Barycentric Frame over Time")
    plt.savefig("plots/time_vs_eccentricity.jpg", bbox_inches = "tight")
    plt.close()

    plt.scatter(extract_coordinate(get_relative_xyz(data_dict["BBH_1"], data_dict["binary"]), "x"), extract_coordinate(get_relative_xyz(data_dict["BBH_1"], data_dict["binary"]), "y"), label = "m1_a")
    plt.scatter(extract_coordinate(get_relative_xyz(data_dict["BBH_2"], data_dict["binary"]), "x"), extract_coordinate(get_relative_xyz(data_dict["BBH_2"], data_dict["binary"]), "y"), label = "m1_b")
    plt.title("Trajectory of Binary BBHs in Barycentric Frame")
    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.legend()
    plt.savefig("plots/BinaryTrajectory.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(data_dict["SMBH"]["time"], distance(data_dict["SMBH"], data_dict["binary"]))
    plt.xlabel("Time [yr]")
    plt.ylabel("Distance [AU]")
    plt.title("Distance from SMBH to $m_1$ over Time")
    plt.savefig("plots/BinaryDisttoSMBH.jpg", bbox_inches = "tight")
    plt.close()

    plt.scatter(extract_coordinate(get_xyz(data_dict["BBH_1"]), "x"), extract_coordinate(get_xyz(data_dict["BBH_1"]), "y"), label = "m1_a")
    plt.scatter(extract_coordinate(get_xyz(data_dict["BBH_2"]), "x"), extract_coordinate(get_xyz(data_dict["BBH_2"]), "y"), label = "m1_b")
    plt.scatter(extract_coordinate(get_xyz(data_dict["perturber"]), "x"), extract_coordinate(get_xyz(data_dict["perturber"]), "y"), label = "perturber")
    plt.scatter(extract_coordinate(get_xyz(data_dict["SMBH"]), "x"), extract_coordinate(get_xyz(data_dict["SMBH"]), "y"), label = "SMBH")
    plt.title("System Trajectory")
    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.legend()
    plt.savefig("plots/SystemTrajectory.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(data_dict["Misc"]["time"],data_dict["Misc"]["perturber-binary separation"])
    plt.title("Separation Between Orbits of Perturber and Binary \n $\\tau = 10^5 T_{m_2}$")
    plt.xlabel("Time [yr]")
    plt.ylabel("Separation [AU]")
    plt.xscale("log")
    plt.axhline(y = 4, color = "black", linestyle = "--", alpha = 0.5)
    plt.axhline(y = 2.5, color = "black", linestyle = "--", alpha = 0.5)
    plt.savefig("plots/DistPerturbBinary.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(data_dict["Misc"]["time"], data_dict["Misc"]["a_perturber"])
    plt.title("Semi-major Axis of Perturber")
    plt.xlabel("Time [yr]")
    plt.ylabel("a [AU]")
    plt.savefig("plots/PerturberSemiMajor.jpg", bbox_inches = "tight")
    plt.close()


def animate_system_trajectory(i, ax, SMBH_scatter, perturber_scatter, BBH_1_scatter, BBH_2_scatter, trail_args):
    SMBH_scatter.set_offsets([SMBH_x[i], SMBH_y[i]])
    perturber_scatter.set_offsets([perturber_x[i], perturber_y[i]])
    BBH_1_scatter.set_offsets([BBH_1_x[i], BBH_1_y[i]])
    BBH_2_scatter.set_offsets([BBH_2_x[i], BBH_2_y[i]])

    trail, points = trail_args
    if i > 0:
        if points >= i:
            points = i
        trail.set_offsets(np.array([np.array(i) for i in zip(perturber_x[i-points:i], perturber_y[i-points:i])]))

    ax.set_title(f"Animated System Trajectory \n t = {np.round(data_dict['BBH_1']['time'][i], 3)} yr")


def animate_binary_trajectory(i, ax, BBH_1_scatter, BBH_2_scatter, trail_args):
    BBH_1_scatter.set_offsets([BBH_1_barycentric_x[i], BBH_1_barycentric_y[i]])
    BBH_2_scatter.set_offsets([BBH_2_barycentric_x[i], BBH_2_barycentric_y[i]])

    BBH_1_trail, BBH_2_trail, trail_points = trail_args

    if i > 0:
        if trail_points >= i:
            trail_points = i
        BBH_1_trail.set_offsets(np.array([np.array(i) for i in zip(BBH_1_barycentric_x[i - trail_points:i], BBH_1_barycentric_y[i - trail_points:i])]))
        BBH_2_trail.set_offsets(np.array([np.array(i) for i in zip(BBH_2_barycentric_x[i - trail_points:i], BBH_2_barycentric_y[i - trail_points:i])]))

    ax.set_title(f"Animated Binary Trajectory \n t = {np.round(data_dict['BBH_1']['time'][i], 3)} yr")


def generate_binary_trajectory_animation():

    plt.close()
    lim = 0.04
    fig, ax = plt.subplots(figsize = (5, 5))
    ax.set_aspect("equal")
    ax.set(xlim = (-lim, lim), ylim = (-lim, lim))
    trail_points = 6

    BBH_1_scatter = ax.scatter(BBH_1_barycentric_x[0], BBH_1_barycentric_y[0], s = 10, facecolor = "blue", zorder = 3, label = "BBH_1")
    BBH_2_scatter = ax.scatter(BBH_2_barycentric_x[0], BBH_2_barycentric_y[0], s = 10, facecolor = "green", zorder = 3, label = "BBH_2")
    BBH_1_trail = ax.scatter(BBH_1_barycentric_x[0] * trail_points,  BBH_1_barycentric_y[0] * trail_points, marker = ".", alpha = 0.35, color = "blue")
    BBH_2_trail = ax.scatter(BBH_2_barycentric_x[0] * trail_points, BBH_2_barycentric_y[0] * trail_points, marker = ".", alpha = 0.35, color = "green")

    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.title("Animated Binary Trajectory")
    plt.legend(loc = "upper right")

    animation = FuncAnimation(fig, animate_binary_trajectory, frames = len(BBH_1_barycentric_x), interval = 100, fargs = (ax, BBH_1_scatter, BBH_2_scatter, (BBH_1_trail, BBH_2_trail, trail_points)))

    animation.save("plots/BinaryTrajectory.mp4")


def generate_system_trajectory_animation():

    plt.close()
    lim = 15
    fig, ax = plt.subplots(figsize = (5, 5))
    ax.set_aspect("equal")
    ax.set(xlim = (-lim, lim), ylim = (-lim, lim))
    trail_points = 50

    SMBH_scatter = ax.scatter(SMBH_x[0], SMBH_y[0], s = 35, marker = "*", facecolor = "yellow", zorder = 3, label = "SMBH")
    perturber_scatter = ax.scatter(perturber_x[0], perturber_y[0], s = 10, facecolor = "red", zorder = 3, label = "Perturber")
    perturber_trail = ax.scatter(perturber_x[0] * trail_points,  perturber_y[0] * trail_points, marker = ".", alpha = 0.15, color = "red")
    BBH_1_scatter = ax.scatter(BBH_1_x[0], BBH_1_y[0], s = 10, facecolor = "blue", zorder = 3, label = "BBH_1")
    BBH_2_scatter = ax.scatter(BBH_2_x[0], BBH_2_y[0], s = 10, facecolor = "green", zorder = 3, label = "BBH_2")

    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.title("Animated System Trajectory")
    plt.legend(loc = "upper right")

    animation = FuncAnimation(fig, animate_system_trajectory, frames = len(SMBH_x), interval = 500, fargs = (ax, SMBH_scatter, perturber_scatter, BBH_1_scatter, BBH_2_scatter, (perturber_trail, trail_points)))

    animation.save("plots/SystemTrajectory.mp4")


construct_plots()
generate_binary_trajectory_animation()
generate_system_trajectory_animation()