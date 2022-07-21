import numpy as np
from Tools import *
import matplotlib.pyplot as plt

perturber_a = 4
binary_separation = 0.3

#####################################################################################

sim = create_simulation()

w = populate_simulation(sim, perturber_a = perturber_a, binary_separation = binary_separation)

binary_period, SMBH_period = get_binary_period(sim), get_SMBH_period(sim)
sim.dt = 0.05 * binary_period

#####################################################################################

outcome_record = {"Result": None, "Minimum Distance Between Binary COM and Perturber": np.inf,
                  "Minimum Distance Between BBHs": np.inf, "Minimum t_GW": np.inf, "Minimum relative t_GW": np.inf}

total_time_steps_completed = 0

time_steps, binary_distances, SMBH_distances, times, eccentricities = [], [], [], [], []

m1_a_xs, m1_a_ys, m1_b_xs, m1_b_ys = [], [], [], []

perturber_xs, perturber_ys = [], []
SMBH_xs, SMBH_ys = [], []
m1_a_helio_xs, m1_a_helio_ys = [], []
m1_b_helio_xs, m1_b_helio_ys = [], []

perturber_m1a_distance, perturber_m1b_distance = [], []


#####################################################################################

while sim.t <= 10**5 * SMBH_period:

    binary_barycenter = barycenter(sim.particles["BBH_1"], sim.particles["BBH_2"])
    distance_perturber_to_binary_COM = np.linalg.norm(np.array(sim.particles["perturber"].xyz) - binary_barycenter)

    SMBH_distance = np.linalg.norm(np.array(sim.particles["SMBH"].xyz) - binary_barycenter)

    distance_between_BBHs = sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"]).d

    t_gw = t_GW(sim, sim.particles["BBH_1"], sim.particles["BBH_2"])
    relative_tGW = t_gw / get_binary_period(sim)

    if distance_between_BBHs < outcome_record["Minimum Distance Between BBHs"]:
        outcome_record["Minimum Distance Between BBHs"] = distance_between_BBHs
    if distance_perturber_to_binary_COM < outcome_record["Minimum Distance Between Binary COM and Perturber"]:
        outcome_record["Minimum Distance Between Binary COM and Perturber"] = distance_perturber_to_binary_COM
    if t_gw < outcome_record["Minimum t_GW"]:
        outcome_record["Minimum t_GW"] = t_gw
    if relative_tGW < outcome_record["Minimum relative t_GW"]:
        outcome_record["Minimum relative t_GW"] = relative_tGW

    if total_time_steps_completed % 10000 == 0:
        times.append(sim.t)
        time_steps.append(sim.dt)
        binary_distances.append(distance_between_BBHs)
        m1_a_xs.append(sim.particles["BBH_1"].x - binary_barycenter[0])
        m1_a_ys.append(sim.particles["BBH_1"].y - binary_barycenter[1])
        m1_b_xs.append(sim.particles["BBH_2"].x - binary_barycenter[0])
        m1_b_ys.append(sim.particles["BBH_2"].y - binary_barycenter[1])
        SMBH_distances.append(SMBH_distance)
        eccentricities.append(sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"]).e)

        SMBH_xs.append(sim.particles["SMBH"].x)
        SMBH_ys.append(sim.particles["SMBH"].y)

        perturber_xs.append(sim.particles["perturber"].x)
        perturber_ys.append(sim.particles["perturber"].y)

        m1_a_helio_xs.append(sim.particles["BBH_1"].x)
        m1_a_helio_ys.append(sim.particles["BBH_1"].y)

        m1_b_helio_xs.append(sim.particles["BBH_2"].x)
        m1_b_helio_ys.append(sim.particles["BBH_2"].y)

        perturber_m1a_distance.append(sim.particles["perturber"].calculate_orbit(primary = sim.particles["BBH_1"]).d)
        perturber_m1b_distance.append(sim.particles["perturber"].calculate_orbit(primary=sim.particles["BBH_2"]).d)

        #---------------Plotting------------------------------#

        plt.plot(times, time_steps)
        plt.xlabel("Time [yr]")
        plt.ylabel("dt [yr]")
        plt.title("Time-step over Time")
        plt.yscale("log", base = 10)
        plt.savefig("plots/timeVsdt.jpg", bbox_inches = "tight")
        plt.close()

        plt.plot(times, binary_distances)
        plt.xlabel("Time [yr]")
        plt.ylabel("Distance [AU]")
        plt.title("Distance Between BHs in Binary over Time")
        plt.savefig("plots/timeVSDistance.jpg", bbox_inches = "tight")
        plt.close()

        plt.plot(times, eccentricities)
        plt.xlabel("Time [yr]")
        plt.ylabel("Eccentricity")
        plt.title("Eccentricity of the Binary in the Barycentric Frame over Time")
        plt.savefig("plots/timeVSEccentricity.jpg", bbox_inches = "tight")
        plt.close()

        plt.scatter(m1_a_xs, m1_a_ys, label = "m1_a")
        plt.scatter(m1_b_xs, m1_b_ys, label = "m1_b")
        plt.title("Trajectory of Binary BBHs in Barycentric Frame")
        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        plt.legend()
        plt.savefig("plots/BinaryTrajectory.jpg", bbox_inches = "tight")
        plt.close()

        plt.plot(times, SMBH_distances)
        plt.xlabel("Time [yr]")
        plt.ylabel("Distance [AU]")
        plt.title("Distance from SMBH to m1 over Time")
        plt.savefig("plots/BinaryDisttoSMBH.jpg", bbox_inches = "tight")
        plt.close()

        plt.scatter(np.array(m1_a_helio_xs), np.array(m1_a_helio_ys), label = "m1_a")
        plt.scatter(np.array(m1_b_helio_xs), np.array(m1_b_helio_ys), label = "m1_b")
        plt.scatter(np.array(perturber_xs), np.array(perturber_ys), label = "perturber")
        plt.scatter(0, 0, label = "SMBH")
        plt.title("System Trajectory")
        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        plt.legend()
        plt.savefig("plots/SystemTrajectory.jpg", bbox_inches = "tight")
        plt.close()

        plt.plot(times, perturber_m1a_distance)
        plt.title("Distance Between m1_a and m2 over Time")
        plt.xlabel("Time [yr]")
        plt.ylabel("Distance [AU]")
        plt.savefig("plots/DistPerturbm1a.jpg", bbox_inches = "tight")
        plt.close()

        plt.plot(times, perturber_m1b_distance)
        plt.title("Distance Between m1_b and m2 over Time")
        plt.xlabel("Time [yr]")
        plt.ylabel("Distance [AU]")
        plt.savefig("plots/DistPerturbm1b.jpg", bbox_inches = "tight")
        plt.close()

    # ---------------Check for Boundedness of Binary------------------------------#
    check_binary_bound(sim, outcome_record)
    #---------------Check for Collisions------------------------------------------#
    check_for_collisions(sim, w, outcome_record)

    sim.steps(1)
    total_time_steps_completed += 1

outcome_record["Result"] = "Ended with no mergers and bound binary"
dump_record(outcome_record)


