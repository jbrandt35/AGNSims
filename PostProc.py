import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy.linalg import norm as mag
import os
from matplotlib.animation import FuncAnimation


def read_data(object):
    data = pd.read_csv(f"data/{object}_data.csv")
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


def distance(p1, p2):
    p1, p2 = get_xyz(p1), get_xyz(p2)
    distances = []
    for i, pos in enumerate(p1):
        dist = mag(p1[i] - p2[i])
        distances.append(dist)
    return np.array(distances)

def construct_plots():

    plt.plot(other_data["time"], other_data["Time-Step"])
    plt.xlabel("Time [yr]")
    plt.ylabel("dt [yr]")
    plt.title("Time-step over Time")
    plt.yscale("log", base = 10)
    plt.savefig("plots/time_vs_dt.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(BBH_1_data["time"], distance(BBH_1_data, BBH_2_data))
    plt.xlabel("Time [yr]")
    plt.ylabel("Distance [AU]")
    plt.title("Distance Between BHs in Binary over Time")
    plt.savefig("plots/time_vs_binary_distance.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(binary_data["time"], binary_data["Eccentricity"])
    plt.xlabel("Time [yr]")
    plt.ylabel("Eccentricity")
    plt.title("Eccentricity of the Binary in the Barycentric Frame over Time")
    plt.savefig("plots/time_vs_eccentricity.jpg", bbox_inches = "tight")
    plt.close()

    plt.scatter(extract_coordinate(get_relative_xyz(BBH_1_data, binary_data), "x"), extract_coordinate(get_relative_xyz(BBH_1_data, binary_data), "y"), label = "m1_a")
    plt.scatter(extract_coordinate(get_relative_xyz(BBH_2_data, binary_data), "x"), extract_coordinate(get_relative_xyz(BBH_2_data, binary_data), "y"), label = "m1_b")
    plt.title("Trajectory of Binary BBHs in Barycentric Frame")
    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.legend()
    plt.savefig("plots/BinaryTrajectory.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(SMBH_data["time"], distance(SMBH_data, binary_data))
    plt.xlabel("Time [yr]")
    plt.ylabel("Distance [AU]")
    plt.title("Distance from SMBH to $m_1$ over Time")
    plt.savefig("plots/BinaryDisttoSMBH.jpg", bbox_inches = "tight")
    plt.close()

    plt.scatter(extract_coordinate(get_xyz(BBH_1_data), "x"), extract_coordinate(get_xyz(BBH_1_data), "y"), label = "m1_a")
    plt.scatter(extract_coordinate(get_xyz(BBH_2_data), "x"), extract_coordinate(get_xyz(BBH_2_data), "y"), label = "m1_b")
    plt.scatter(extract_coordinate(get_xyz(perturber_data), "x"), extract_coordinate(get_xyz(perturber_data), "y"), label = "perturber")
    plt.scatter(extract_coordinate(get_xyz(SMBH_data), "x"), extract_coordinate(get_xyz(SMBH_data), "y"), label = "SMBH")
    plt.title("System Trajectory")
    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.legend()
    plt.savefig("plots/SystemTrajectory.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(other_data["time"], other_data["Perturber-Binary Separation"])
    plt.title("Separation Between Orbits of Perturber and Binary \n $\\tau = 10^5 T_{m_2}$")
    plt.xlabel("Time [yr]")
    plt.ylabel("Separation [AU]")
    plt.xscale("log")
    plt.axhline(y = 4, color = "black", linestyle = "--", alpha = 0.5)
    plt.axhline(y = 2.5, color = "black", linestyle = "--", alpha = 0.5)
    plt.savefig("plots/DistPerturbBinary.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(other_data["time"], other_data["a_perturber"])
    plt.title("Semi-major Axis of Perturber")
    plt.xlabel("Time [yr]")
    plt.ylabel("a [AU]")
    plt.savefig("plots/PerturberSemiMajor.jpg", bbox_inches = "tight")
    plt.close()

def split_coordinates(data, mode = "x,y"):
    if mode == "x,y":
        return extract_coordinate(get_xyz(data), "x"), extract_coordinate(get_xyz(data), "y")


BBH_1_data = read_data("BBH_1")
BBH_2_data = read_data("BBH_2")
SMBH_data = read_data("SMBH")
perturber_data = read_data("perturber")
binary_data = read_data("binary")
other_data = read_data("other")

SMBH_x, SMBH_y = split_coordinates(SMBH_data)
BBH_1_x, BBH_1_y = split_coordinates(BBH_1_data)
BBH_2_x, BBH_2_y = split_coordinates(BBH_2_data)
perturber_x, perturber_y = split_coordinates(perturber_data)
BBH_1_barycentric_x, BBH_1_barycentric_y = extract_coordinate(get_relative_xyz(BBH_1_data, binary_data), "x"), extract_coordinate(get_relative_xyz(BBH_1_data, binary_data), "y")
BBH_2_barycentric_x, BBH_2_barycentric_y = extract_coordinate(get_relative_xyz(BBH_2_data, binary_data), "x"), extract_coordinate(get_relative_xyz(BBH_2_data, binary_data), "y")


def animate_system_trajectory(i, SMBH_scatter, perturber_scatter, BBH_1_scatter, BBH_2_scatter, trail_args):
    SMBH_scatter.set_offsets([SMBH_x[i], SMBH_y[i]])
    perturber_scatter.set_offsets([perturber_x[i], perturber_y[i]])
    BBH_1_scatter.set_offsets([BBH_1_x[i], BBH_1_y[i]])
    BBH_2_scatter.set_offsets([BBH_2_x[i], BBH_2_y[i]])

    trail, points = trail_args
    if i > 0:
        if points >= i:
            points = i
        trail.set_offsets(np.array([np.array(i) for i in zip(perturber_x[i-points:i], perturber_y[i-points:i])]))



def animate_binary_trajectory(i, BBH_1_scatter, BBH_2_scatter, trail_args):
    BBH_1_scatter.set_offsets([BBH_1_barycentric_x[i], BBH_1_barycentric_y[i]])
    BBH_2_scatter.set_offsets([BBH_2_barycentric_x[i], BBH_2_barycentric_y[i]])

    BBH_1_trail, BBH_2_trail, trail_points = trail_args

    if i > 0:
        if trail_points >= i:
            trail_points = i
        BBH_1_trail.set_offsets(np.array([np.array(i) for i in zip(BBH_1_barycentric_x[i-trail_points:i],  BBH_1_barycentric_y[i-trail_points:i])]))
        BBH_2_trail.set_offsets(np.array([np.array(i) for i in zip(BBH_2_barycentric_x[i - trail_points:i], BBH_2_barycentric_y[i - trail_points:i])]))



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

    animation = FuncAnimation(fig, animate_binary_trajectory, interval = 100, fargs = (BBH_1_scatter, BBH_2_scatter, (BBH_1_trail, BBH_2_trail, trail_points)))

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

    animation = FuncAnimation(fig, animate_system_trajectory, interval = 30, fargs = (SMBH_scatter, perturber_scatter, BBH_1_scatter, BBH_2_scatter, (perturber_trail, trail_points)))

    animation.save("plots/SystemTrajectory.mp4")




generate_binary_trajectory_animation()
generate_system_trajectory_animation()