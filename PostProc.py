import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
        dist = np.linalg.norm(p1[i] - p2[i])
        distances.append(dist)
    return np.array(distances)

def construct_plots():

    other_data = read_data("other")
    BBH_1_data = read_data("BBH_1")
    BBH_2_data = read_data("BBH_2")
    SMBH_data = read_data("SMBH")
    perturber_data = read_data("perturber")
    binary_data = read_data("binary")

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

construct_plots()


