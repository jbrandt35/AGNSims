import json
import rebound
import numpy as np
from numpy.linalg import norm as mag
import matplotlib.pyplot as plt
from astropy import constants
import pandas as pd

m0, m1_a, m1_b, m2 = 10 ** 6, 10, 10, 20
c = constants.c.to("AU/yr").value


#Returns a rebound simulation object with the appropriate units
def create_simulation():
    sim = rebound.Simulation()
    sim.units = ["Msun", "yr", "AU"]
    return sim


# Returns the velocity of the COM of the binary given its orbital elements. Note: SMBH needs to be in sim.particles
def get_binary_COM_data(m_SMBH, m_binary, a, e = 0, M = 0, inc = 0):

    #Create a simulation just for getting these values
    virtual_sim = create_simulation()

    #Put in the SMBH
    virtual_sim.add(m = m_SMBH, hash = "SMBH")

    #Add a virtual particle to represent the COM of the binary
    virtual_sim.add(primary = virtual_sim.particles["SMBH"], m = m_binary, a = a, e = e, M = M, inc = inc, hash = "binary_COM")

    #Get its position
    position = np.array(virtual_sim.particles["binary_COM"].xyz)

    #Get its velocity
    velocity = np.array(virtual_sim.particles["binary_COM"].vxyz)

    #Get its Hill Radius
    R_hill = virtual_sim.particles["binary_COM"].rhill

    return position, velocity, R_hill


def get_BBH_data(m1, m2, a, e = 0, M = 0, inc = 0):

    # Create a simulation just for getting these values
    virtual_sim = create_simulation()

    #Add m1 as the primary object
    virtual_sim.add(m = m1, hash = "BH_1")

    #Add m2 around it
    virtual_sim.add(primary = virtual_sim.particles["BH_1"], m = m2, a = a, e = e, M = M, inc = inc, hash = "BH_2")

    #Switch to COM frame
    virtual_sim.move_to_com()

    #Get positions
    BH_1_position = np.array(virtual_sim.particles["BH_1"].xyz)
    BH_2_position = np.array(virtual_sim.particles["BH_2"].xyz)

    #Get velocities
    BH_1_velocity = np.array(virtual_sim.particles["BH_1"].vxyz)
    BH_2_velocity = np.array(virtual_sim.particles["BH_2"].vxyz)

    BH_1_data = BH_1_position, BH_1_velocity
    BH_2_data = BH_2_position, BH_2_velocity

    return BH_1_data, BH_2_data


#Takes in sim object and orbital elements, populates the simulation with the objects accordingly
#SMBH_a in units of Gravitational Radii, binary_separation in units of m1 Hill radii, perturber_a in units of m1+m2 Hill radii
#Returns w - for calculating mergers between m1 and m2
#Each particle is hashed with a name useful for other functions - i.e. see period calculation functions
def populate_simulation(sim, binary_separation = 0.1, binary_eccentricity = 0, binary_M = 0, binary_inc = 0, SMBH_a = 1000, SMBH_eccentricity = 0, SMBH_M = 0, SMBH_inc = 0, perturber_a = 20, perturber_e = 0, perturber_M = 0, perturber_inc = 0, randomize_M = False):

    #Calculate Gravitational Radius and convert units of SMBH_a
    Rg = sim.G * m0 / c ** 2
    SMBH_a *= Rg

    #Calculate m1+m2 Hill radius and find the semi-major axis of the perturber
    mass_factor = np.power((m1_a+m1_b+m2)/(3*m0), 1/3)
    perturber_factor = (1 + (perturber_a/2) * mass_factor) / (1 - (perturber_a/2) * mass_factor)
    perturber_a = SMBH_a * perturber_factor

    #If randomized mean anamolies are wanted, generate three numbers 0-2pi and set mean anamolies
    if randomize_M:
        SMBH_M, binary_M, perturber_M = 2 * np.pi * np.random.rand(3)

    #Add the SMBH to the sim
    sim.add(m = m0, hash = "SMBH")

    #Get position, velocity, and Hill radius of binary's COM
    binary_COM_position, binary_COM_velocity, m1_R_hill = get_binary_COM_data(m0, m1_a + m1_b, SMBH_a, e = SMBH_eccentricity, M = SMBH_M, inc = SMBH_inc)

    #Change units of binary_separation
    binary_separation *= m1_R_hill

    #Get data of BHs in binary (in binary COM F.O.R.)
    BH_1_data, BH_2_data = get_BBH_data(m1_a, m1_b, binary_separation, e = binary_eccentricity, M = binary_M, inc = binary_inc)

    #Unpack BBH data
    BH_1_position, BH_1_velocity = BH_1_data
    BH_2_position, BH_2_velocity = BH_2_data

    #Create particles for BBHs
    BH_1 = rebound.Particle(m = m1_a, hash = "BBH_1")
    BH_1.xyz = BH_1_position + binary_COM_position
    BH_1.vxyz = BH_1_velocity + binary_COM_velocity

    BH_2 = rebound.Particle(m = m1_b, hash = "BBH_2")
    BH_2.xyz = BH_2_position + binary_COM_position
    BH_2.vxyz = BH_2_velocity + binary_COM_velocity

    #Add BBHs into sim
    sim.add(BH_1)
    sim.add(BH_2)

    #Add perturber into sim
    sim.add(primary = sim.particles["SMBH"], m = m2, a = perturber_a, e = perturber_e, M = perturber_M, inc = perturber_inc, hash = "perturber")

    #Move to COM frame of simulation
    sim.move_to_com()

    return np.sqrt(sim.G * m0 / SMBH_a) - np.sqrt(sim.G * m0 / perturber_a) + np.sqrt(
        sim.G * (m1_a + m1_b) / (binary_separation ** 3)) * binary_separation * (m1_b/(m1_a+m1_b))


def get_binary_period(sim):
    orbit = sim.particles["BBH_1"].calculate_orbit(primary = sim.particles["BBH_2"])
    return orbit.P


def get_binary_SMBH_period(sim):
    binary_barycenter = sim.calculate_com(first = 1, last = 3)
    virtual_sim = create_simulation()
    virtual_sim.add(m = m0)
    virtual_sim.add(binary_barycenter)
    orbit = virtual_sim.particles[1].calculate_orbit(primary = virtual_sim.particles[0])
    return orbit.P


def get_perturber_period(sim):
    orbit = sim.particles["perturber"].calculate_orbit(primary = sim.particles["SMBH"])
    return orbit.P


def is_bound(primary, secondary):
    orbit = secondary.calculate_orbit(primary = primary)
    return bool(orbit.e < 1)


def t_GW(sim, particle1, particle2):

    orbit = particle1.calculate_orbit(primary = particle2)

    m1 = particle1.m
    m2 = particle2.m
    a = orbit.a
    e = orbit.e

    first_term = -64/5
    second_term = (sim.G ** 3) * (m1 + m2) * m1 * m2 / ((c**5) * (a**3) * ((1-e**2)**(7/2)))
    third_term = 1 + ((73/24)*(e**2)) + ((37/94)*(e**4))

    da_dt = first_term * second_term * third_term

    t_gw = np.abs(a/da_dt)
    return t_gw


def check_for_collisions(sim, w, record):

    ##############Check Binary Orbits###############

    BH_a = sim.particles["BBH_1"]
    BH_b = sim.particles["BBH_2"]

    sum_of_event_horizons = 2 * sim.G * (BH_a.m + BH_b.m) / (c ** 2)
    t_gw = t_GW(sim, BH_a, BH_b)

    orbit = BH_a.calculate_orbit(primary = BH_b)
    distance, period = orbit.d, orbit.P

    if distance < sum_of_event_horizons:
        record["Result"] = f"Collision Encountered: The distance between m1_a and m1_b was {distance} AU when the sum of event horizon radii was {sum_of_event_horizons} AU."
        dump_record(record)
        raise CollisionException(record["Result"])
    elif t_gw < period:
        record["Result"] = f"Collision Encountered: t_GW was {t_gw} when the period was {period}."
        dump_record(record)
        raise CollisionException(record["Result"])

    ##############Check Perturber Orbits###############

    perturber = sim.particles["perturber"]

    n = (BH_a.m * perturber.m) / ((BH_a.m + perturber.m)**2)
    M_tot = BH_a.m + perturber.m

    first_term = (85 * np.pi / (6 * np.sqrt(2)))**(2/7)
    second_term = M_tot * (n**(2/7)) / (w**(4/7))
    third_term = 1 - (1/4) * ((85 * np.pi / 3)**(2/7)) * ((4*n)**(2/7)) * (w**(10/7))

    r_p = first_term * second_term * third_term * (sim.G/(c**2))

    m1_a_distance = BH_a.calculate_orbit(primary = perturber).d
    m1_b_distance = BH_b.calculate_orbit(primary = perturber).d

    if m1_a_distance < r_p:
        record["Result"] = f"Collision Encountered: The distance between m1_a and m2 was {m1_a_distance} AU when r_p was {r_p} AU."
        dump_record(record)
        raise CollisionException(record["Result"])

    if m1_b_distance < r_p:
        record["Result"] = f"Collision Encountered: The distance between m1_b and m2 was {m1_b_distance} AU when r_p was {r_p} AU."
        dump_record(record)
        raise CollisionException(record["Result"])


def dump_record(record):
    with open("outcome.json", "w") as file:
        json.dump(record, file)


def check_binary_bound(sim, record):
    if not is_bound(sim.particles["BBH_1"], sim.particles["BBH_2"]):
        record["Result"] = f"Unbound: Binary eccentricity reached {sim.particles['BBH_2'].calculate_orbit(primary = sim.particles['BBH_1']).e}"
        dump_record(record)
        raise UnboundException(record["Result"])


def dist_between(particle_1, particle_2):
    particle_1_position, particle_2_position = np.array(particle_1.xyz), np.array(particle_2.xyz)
    return mag(particle_1_position - particle_2_position)

def get_perturber_binary_separation(binary_COM, perturber):

    virtual_sim = create_simulation()
    virtual_sim.add(m = m0)
    virtual_sim.add(binary_COM)

    a_binary_COM = virtual_sim.particles[1].calculate_orbit(primary = virtual_sim.particles[0]).a
    a_perturber = perturber.calculate_orbit().a

    del virtual_sim

    mass_factor = np.power(3*m0 / (m1_a + m1_b + m2), 1/3)
    a_factor = (a_perturber - a_binary_COM) / (a_perturber + a_binary_COM)

    n = 2 * mass_factor * a_factor

    return n


def check_and_assign_minimums(sim, record):
    # Get binary COM as a rebound.Particle
    binary_barycenter = sim.calculate_com(first = 1, last = 3)

    distance_perturber_to_binary_COM = dist_between(sim.particles["perturber"], binary_barycenter)
    distance_between_BBHs = sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"]).d

    t_gw = t_GW(sim, sim.particles["BBH_1"], sim.particles["BBH_2"])
    relative_tGW = t_gw / get_binary_period(sim)

    if distance_between_BBHs < record["Minimum Distance Between BBHs"]:
        record["Minimum Distance Between BBHs"] = distance_between_BBHs
    if distance_perturber_to_binary_COM < record["Minimum Distance Between Binary COM and Perturber"]:
        record["Minimum Distance Between Binary COM and Perturber"] = distance_perturber_to_binary_COM
    if t_gw < record["Minimum t_GW"]:
        record["Minimum t_GW"] = t_gw
    if relative_tGW < record["Minimum relative t_GW"]:
        record["Minimum relative t_GW"] = relative_tGW

def initialize_data_collection():
    for BH in ["SMBH", "BBH_1", "BBH_2", "perturber"]:
        with open(f"data/{BH}_data.csv", "w") as file:
            file.write(",".join(["time", "x", "y", "z"]) + "\n")


def establish_dataframes():
    BH_df, mixed_df = pd.DataFrame(columns = ["time", "x", "y", "z"]), pd.DataFrame()

    SMBH_df = BH_df.copy(deep = True)
    BBH_1_df = BH_df.copy(deep = True)
    BBH_2_df = BH_df.copy(deep = True)
    perturber_df = BH_df.copy(deep = True)

    return [("SMBH", SMBH_df), ("BBH_1", BBH_1_df), ("BBH_2", BBH_2_df), ("perturber", perturber_df)]

def save_in_frame(frame, data):
    frame.loc[len(frame.index)] = data


def save_data(sim):
    for hash, frame in establish_dataframes():
        data = [sim.t] + sim.particles[hash].xyz
        save_in_frame(frame, data)
        frame.to_csv(f"data/{hash}_data.csv", mode = "a", header = False)


def save_plotting_data(sim):
    binary_barycenter = sim.calculate_com(first = 1, last = 3)

    SMBH_distance = dist_between(sim.particles["SMBH"], binary_barycenter)
    distance_between_BBHs = sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"]).d

    times.append(sim.t)
    time_steps.append(sim.dt)
    binary_distances.append(distance_between_BBHs)
    m1_a_xs.append(sim.particles["BBH_1"].x - binary_barycenter.x)
    m1_a_ys.append(sim.particles["BBH_1"].y - binary_barycenter.y)
    m1_b_xs.append(sim.particles["BBH_2"].x - binary_barycenter.x)
    m1_b_ys.append(sim.particles["BBH_2"].y - binary_barycenter.y)
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
    perturber_m1b_distance.append(sim.particles["perturber"].calculate_orbit(primary = sim.particles["BBH_2"]).d)
    perturber_binary_separation.append(get_perturber_binary_separation(binary_barycenter, sim.particles["perturber"]))
    perturber_a.append(sim.particles["perturber"].calculate_orbit(primary = sim.particles["SMBH"]).a)


def construct_plots():

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

    plt.plot(times, perturber_binary_separation)
    plt.title("Separation Between Orbits of Perturber and Binary \n $\\tau = 10^5 T_{m_2}$")
    plt.xlabel("Time [yr]")
    plt.ylabel("Separation [AU]")
    plt.xscale("log")
    plt.axhline(y = 4, color = "black", linestyle = "--", alpha = 0.5)
    plt.axhline(y = 2.5, color = "black", linestyle = "--", alpha = 0.5)
    plt.savefig("plots/DistPerturbBinary.jpg", bbox_inches = "tight")
    plt.close()

    plt.plot(times, perturber_a)
    plt.title("Semi-major Axis of Perturber")
    plt.xlabel("Time [yr]")
    plt.ylabel("a [AU]")
    plt.savefig("plots/PerturberSemiMajor.jpg", bbox_inches = "tight")
    plt.close()


class UnboundException(Exception):
    pass

class CollisionException(Exception):
    pass


#####################################################################################

outcome_record = {"Result": None, "Minimum Distance Between Binary COM and Perturber": np.inf,
                  "Minimum Distance Between BBHs": np.inf, "Minimum t_GW": np.inf, "Minimum relative t_GW": np.inf}

total_time_steps_completed = 0

time_steps, binary_distances, SMBH_distances, times, eccentricities = [], [], [], [], []

perturber_binary_separation = []
perturber_a = []

m1_a_xs, m1_a_ys, m1_b_xs, m1_b_ys = [], [], [], []

perturber_xs, perturber_ys = [], []
SMBH_xs, SMBH_ys = [], []
m1_a_helio_xs, m1_a_helio_ys = [], []
m1_b_helio_xs, m1_b_helio_ys = [], []

perturber_m1a_distance, perturber_m1b_distance = [], []

#####################################################################################

