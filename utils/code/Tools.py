import json
import rebound
import numpy as np
from numpy.linalg import norm as mag
from astropy import constants
import pandas as pd
import os
import config
import copy


data_objects = ["BBH_1", "BBH_2", "SMBH", "perturber"]

m0, m1_a, m1_b, m2 = 10 ** 6, 10, 10, 20
c = constants.c.to("AU/yr").value


#Returns a rebound simulation object with the appropriate units
def create_simulation():
    sim = rebound.Simulation()
    sim.units = ["Msun", "yr", "AU"]
    return sim

def calc_drag_force(tau, vxyz, G, M, r_0, r_pert):
    vx, vy, vz = vxyz
    r = [r_pert[0] - r_0[0], r_pert[1] - r_0[1], r_pert[2] - r_0[2]]
    x, y, z = r
    d = np.sqrt(x**2 + y**2 + z**2)
    beta = np.sqrt(G*M/(d**5))
    ax = (-1/tau)*(vx + y*beta)
    ay = (-1/tau)*(vy - x*beta)
    az = (-1/tau)*(vz)
    return ax, ay, az

def calc_trap_force(tau, G, M, a_bin, r_pert, r_SMBH):
    omega = np.sqrt(G*M/(a_bin**3))
    r = [r_pert[0] - r_SMBH[0], r_pert[1] - r_SMBH[1], r_pert[2] - r_SMBH[2]]
    x, y, z = r
    d = np.sqrt(x**2 + y**2 + z**2)
    F_mag = -omega*(d - a_bin)/tau
    ax = F_mag*(-y/d)
    ay = F_mag*(x/d)
    az = 0
    return ax, ay, az


#Secunda Distributions for Orbital Elements
def sample_inclination_distribution():
    inc = np.random.normal(loc = 0, scale = np.radians(3))
    return np.abs(inc)

def sample_eccentricity_distribution():
    while True:
        e = np.random.normal(loc = 0.05, scale = 0.02)
        if e > 0:
            return e

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

def unit_mutual_orbital_angular_momentum(particle_1, particle_2):
    r = np.array(particle_1.xyz) - np.array(particle_2.xyz)
    v = np.array(particle_1.vxyz) - np.array(particle_2.vxyz)
    L_hat = normalize(np.cross(r, particle_1.m * v))
    return L_hat


#Takes in sim object and orbital elements, populates the simulation with the objects accordingly
#SMBH_a in units of Gravitational Radii, binary_separation in units of m1 Hill radii, perturber_a in units of m1+m2 Hill radii
#Returns w - for calculating mergers between m1 and m2
#Each particle is hashed with a name useful for other functions - i.e. see period calculation functions
def populate_simulation(sim, binary_separation = 0.1, binary_eccentricity = 0, binary_M = 0, binary_inc = 0, SMBH_a = 1000, SMBH_eccentricity = 0, SMBH_M = 0, SMBH_inc = 0, perturber_a = 20, perturber_e = 0, perturber_M = 0, perturber_inc = 0, randomize_M = False, ignore_perturber = False, randomize_eccentricities = False, randomize_binary_inc = False):

    #Calculate Gravitational Radius and convert units of SMBH_a
    Rg = sim.G * m0 / c ** 2
    SMBH_a *= Rg

    #Calculate m1+m2 Hill radius and find the semi-major axis of the perturber
    mass_factor = np.power((m1_a+m1_b+m2)/(3*m0), 1/3)
    perturber_factor = (1 + (perturber_a/2) * mass_factor) / (1 - (perturber_a/2) * mass_factor)
    perturber_a = SMBH_a * perturber_factor

    if randomize_M:
        SMBH_M, binary_M, perturber_M = 2 * np.pi * np.random.rand(3)
    if randomize_eccentricities:
        SMBH_eccentricity = sample_eccentricity_distribution()
        binary_eccentricity = sample_eccentricity_distribution()
        perturber_e = sample_eccentricity_distribution()
    if randomize_binary_inc:
        binary_inc = sample_inclination_distribution()

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

    if not ignore_perturber:
        #Add perturber into sim
        sim.add(primary = sim.particles["SMBH"], m = m2, a = perturber_a, e = perturber_e, M = perturber_M, inc = perturber_inc, hash = "perturber")
    else:
        #Add perturber as test particle
        sim.add(primary=sim.particles["SMBH"], a = perturber_a, e = perturber_e, M = perturber_M, inc = perturber_inc, hash = "perturber")
        sim.N_active = 3

    #Move to COM frame of simulation
    sim.move_to_com()

    if config.mode == "initial_spin_aligned_with_L_of_BBH2":
        config.spin = unit_mutual_orbital_angular_momentum(sim.particles["BBH_2"], sim.particles["BBH_1"])
    elif config.mode == "initial_spin_aligned_with_L_of_Binary":
        config.spin = np.array([0, 0, 1])

    config.outcome_record["Initial Inclination of BBH_2 around BBH_1"] = np.degrees(binary_inc)

    return np.sqrt(sim.G * m0 / SMBH_a) - np.sqrt(sim.G * m0 / perturber_a) + np.sqrt(
        sim.G * (m1_a + m1_b) / (binary_separation ** 3)) * binary_separation * (m1_b/(m1_a+m1_b))


def get_binary_period(sim):
    orbit = sim.particles["BBH_1"].calculate_orbit(primary = sim.particles["BBH_2"])
    return orbit.P

def normalize(vector):
    return np.array(vector)/mag(vector)

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


def t_GW(sim):

    orbit = sim.particles["BBH_1"].calculate_orbit(primary = sim.particles["BBH_2"])

    m1 = sim.particles["BBH_1"].m
    m2 = sim.particles["BBH_2"].m
    a = orbit.a
    e = orbit.e

    first_term = -64/5
    second_term = (sim.G ** 3) * (m1 + m2) * m1 * m2 / ((c**5) * (a**3) * ((1-e**2)**(7/2)))
    third_term = 1 + ((73/24)*(e**2)) + ((37/94)*(e**4))

    da_dt = first_term * second_term * third_term

    t_gw = np.abs(a/da_dt)

    if e < 1:
        return t_gw
    else:
        return -1

def system_ejected(primary, secondary, limit):
    d = dist_between(primary, secondary)
    return bool(d > limit)

def save_initial_data(sim):
    config.outcome_record["Initial t_GW"] = t_GW(sim)

def save_final_data(sim):

    binary_orbit = sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"])
    final_binary_eccentricity = binary_orbit.e

    final_t_GW = t_GW(sim)

    config.outcome_record["Final t_GW"] = final_t_GW
    config.outcome_record["Final Binary Eccentricity"] = final_binary_eccentricity
    config.outcome_record["Final Binary Semimaj"] = binary_orbit.a

def check_for_collisions(sim, w):

    ##############Check Binary Orbit for Event Horizon Crossings###############

    BH_a = sim.particles["BBH_1"]
    BH_b = sim.particles["BBH_2"]

    sum_of_event_horizons = 2 * sim.G * (BH_a.m + BH_b.m) / (c ** 2)

    orbit = BH_a.calculate_orbit(primary = BH_b)
    distance, period = orbit.d, orbit.P

    if distance < sum_of_event_horizons:
        config.outcome_record["Result"] = f"Collision Encountered: The distance between m1_a and m1_b was {distance} AU when the sum of event horizon radii was {sum_of_event_horizons} AU."
        save_final_data(sim)
        dump_record()
        raise CollisionException(config.outcome_record["Result"])
    elif orbit.e < 1:
        t_gw = t_GW(sim)
        if t_gw < period:
            config.outcome_record["Result"] = f"GW Collision: At t = {sim.t}, t_GW was {t_gw} when the period was {period}."
            save_final_data(sim)
            dump_record()
            raise CollisionException(config.outcome_record["Result"])


    ##############Check Perturber-Binary BBH Orbits for Event Horizon Crossings###############

    perturber = sim.particles["perturber"]

    for BBH in [BH_a, BH_b]:

        sum_of_event_horizons = 2 * sim.G * (perturber.m + BBH.m) / (c ** 2)
        distance = perturber.calculate_orbit(primary = BBH).d

        if distance < sum_of_event_horizons:
            config.outcome_record["Result"] = f"Collision Encountered: The distance between the perturber and a BBH was {distance} AU when the sum of event horizon radii was {sum_of_event_horizons} AU."
            save_final_data(sim)
            dump_record()
            raise CollisionException(config.outcome_record["Result"])


    ##############Check Perturber Orbit Flags###############

    n = (BH_a.m * perturber.m) / ((BH_a.m + perturber.m)**2)
    M_tot = BH_a.m + perturber.m

    first_term = (85 * np.pi / (6 * np.sqrt(2)))**(2/7)
    second_term = M_tot * (n**(2/7)) / (w**(4/7))
    third_term = 1 - (1/4) * ((85 * np.pi / 3)**(2/7)) * ((4*n)**(2/7)) * (w**(10/7))

    r_p = first_term * second_term * third_term * (sim.G/(c**2))

    m1_a_distance = perturber.calculate_orbit(primary = BH_a).d
    m1_b_distance = perturber.calculate_orbit(primary = BH_b).d

    if m1_a_distance < r_p:
        config.outcome_record["Events"].append(f"At t = {sim.t}, the distance between m1_a and m2 was {m1_a_distance} AU when r_p was {r_p} AU.")

    if m1_b_distance < r_p:
        config.outcome_record["Events"].append(f"At t = {sim.t}, the distance between m1_b and m2 was {m1_b_distance} AU when r_p was {r_p} AU.")

def check_bound(sim):

    Rg = sim.G * m0 / c ** 2
    SMBH_a = 1000 * Rg
    limit = 10 * SMBH_a

    SMBH, BBH_1, BBH_2, perturber = sim.particles

    record = config.outcome_record

    if system_ejected(SMBH, BBH_1, limit):
        if is_bound(BBH_2, perturber):
            record["Result"] = f"Ejected: BBH_1 left the system and was replaced by perturber in binary, BBH_2 + pertuber ecc: {BBH_2.calculate_orbit(primary = perturber).e}"
            save_final_data(sim)
            dump_record()
        else:
            record["Result"] = f"Ejected: BBH_1 left the system without replacement, BBH_2 + pertuber ecc: {BBH_2.calculate_orbit(primary = perturber).e}"
            save_final_data(sim)
            dump_record()
        raise UnboundException(record["Result"])

    if system_ejected(SMBH, BBH_2, limit):
        if is_bound(BBH_1, perturber):
            record["Result"] = f"Ejected: BBH_2 left the system and was replaced by perturber in binary, BBH_1 + pertuber ecc: {BBH_1.calculate_orbit(primary = perturber).e}"
            save_final_data(sim)
            dump_record()
        else:
            record["Result"] = f"Ejected: BBH_2 left the system without replacement, BBH_1 + pertuber ecc: {BBH_1.calculate_orbit(primary = perturber).e}"
            save_final_data(sim)
            dump_record()
        raise UnboundException(record["Result"])

    if system_ejected(SMBH, perturber, limit):
        if is_bound(BBH_1, BBH_2):
            record["Result"] = f"Ejected: perturber left the system and binary remained intact, binary ecc: {BBH_2.calculate_orbit(primary = BBH_1).e}"
            save_final_data(sim)
            dump_record()
        else:
            record["Result"] = f"Ejected: perturber left the system and binary was separated, binary ecc: {BBH_2.calculate_orbit(primary = BBH_1).e}"
            save_final_data(sim)
            dump_record()
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


def check_and_assign_minimums(sim):
    # Get binary COM as a rebound.Particle
    binary_barycenter = sim.calculate_com(first = 1, last = 3)

    record = config.outcome_record

    distance_perturber_to_binary_COM = dist_between(sim.particles["perturber"], binary_barycenter)
    distance_between_BBHs = sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"]).d

    t_gw = t_GW(sim)
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
    os.system("rm -r -f result")
    os.system("mkdir result")


def save_to_frame(df, data):
    df.loc[len(df.index)] = data


def save_data_to_buffer(sim):
    binary_COM = sim.calculate_com(first = 1, last = 3)
    binary_COM_velocity = np.array(binary_COM.vxyz)
    binary_COM_position = np.array(binary_COM.xyz)
    binary_COM_L = np.cross(binary_COM_position, binary_COM.m * binary_COM_velocity)
    binary_COM_L_hat = binary_COM_L / mag(binary_COM_L)

    BBH_2_L_hat = unit_mutual_orbital_angular_momentum(sim.particles["BBH_2"], sim.particles["BBH_1"])

    save_to_frame(buffers["Misc"], [sim.t, sim.dt, sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"]).e, get_perturber_binary_separation(binary_COM, sim.particles["perturber"]), sim.particles["perturber"].calculate_orbit(primary = sim.particles["SMBH"]).a])
    save_to_frame(buffers["binary"], [binary_COM.x, binary_COM.y, binary_COM.z, t_GW(sim), sim.particles["BBH_1"].calculate_orbit(primary = sim.particles["BBH_2"]).a, binary_COM.vx, binary_COM.vy, binary_COM.vz, config.spin[0], config.spin[1], config.spin[2], BBH_2_L_hat[0], BBH_2_L_hat[1], BBH_2_L_hat[2], binary_COM_L_hat[0], binary_COM_L_hat[1], binary_COM_L_hat[2]])

    for particle in data_objects:
        save_to_frame(buffers[particle], [sim.particles[particle].x, sim.particles[particle].y, sim.particles[particle].z, sim.particles[particle].vx, sim.particles[particle].vy, sim.particles[particle].vz])


def save_data_to_disk():
    with pd.HDFStore("result/data.h5") as data_file:
        data_file.append("Misc", buffers["Misc"], complib = "zlib", complevel = 5)
        data_file.append(f"Positions/binary", buffers["binary"], complib = "zlib", complevel = 5)
        for particle in data_objects:
            data_file.append(f"/Positions/{particle}", buffers[particle], complib = "zlib", complevel = 5)


def clear_buffer():
    buffers["Misc"] = pd.DataFrame(columns = ["time", "time-step", "binary eccentricity", "perturber-binary separation", "a_perturber"])
    buffers["binary"] = pd.DataFrame(columns = ["x", "y", "z", "t_GW", "a", "vx", "vy", "vz", "S_x", "S_y", "S_z", "L_BBH2_x", "L_BBH2_y", "L_BBH2_z", "L_Binary_x", "L_Binary_y", "L_Binary_z"])
    for particle in data_objects:
        buffers[particle] = pd.DataFrame(columns = ["x", "y", "z", "vx", "vy", "vz"])


def save_data(sim):
    save_data_to_buffer(sim)
    if buffers["Misc"].shape[0] >= 10000:
        save_data_to_disk()
        clear_buffer()

def dump_record():
    with open("outcome.json", "w") as file:
        json.dump(config.outcome_record, file)
    save_data_to_disk()

def check_for_binary_swaps(sim):
    updated_bound_pairs = []
    interacting_BHs = ["BBH_1", "BBH_2", "perturber"]
    for object_a, object_b in zip(interacting_BHs, interacting_BHs[1:]):
        if is_bound(sim.particles[object_a], sim.particles[object_b]):
            pair = (object_a, object_b)
            updated_bound_pairs.append(pair)
    if not set(updated_bound_pairs) == set(config.bound_pairs):
        new_event = f"Swap Event: At t = {sim.t}, the bound pairs went from {config.bound_pairs} to {updated_bound_pairs}."
        config.bound_pairs = copy.deepcopy(updated_bound_pairs)
        for (object_a, object_b) in updated_bound_pairs:
            orbit = sim.particles[object_a].calculate_orbit(primary = sim.particles[object_b])
            new_event += f" The orbit ({object_a}, {object_b}) has a = {orbit.a} AU and e = {orbit.e}."
        config.outcome_record["Events"].append(new_event)

class UnboundException(Exception):
    pass


class CollisionException(Exception):
    pass

class SwapException(Exception):
    pass


#####################################################################################

buffers = {
    "Misc": pd.DataFrame(columns = ["time", "time-step", "binary eccentricity", "perturber-binary separation", "a_perturber"]),
    "binary": pd.DataFrame(columns = ["x", "y", "z", "t_GW", "a", "vx", "vy", "vz", "S_x", "S_y", "S_z", "L_BBH2_x", "L_BBH2_y", "L_BBH2_z", "L_Binary_x", "L_Binary_y", "L_Binary_z"])
}
for particle in data_objects:
    buffers[particle] = pd.DataFrame(columns = ["x", "y", "z", "vx", "vy", "vz"])
#####################################################################################

