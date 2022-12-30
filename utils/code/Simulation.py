from Tools import *
import numpy as np

perturber_a = PERTSEPARATION
binary_separation = BINSEPARATION

mode = "initial_spin_aligned_with_L_of_Binary"

#####################################################################################

global inclination_of_binary

sim = create_simulation()
w = populate_simulation(sim, mode, perturber_a = perturber_a, binary_separation = binary_separation, ignore_perturber = False, randomize_M = True, binary_inc = inclination_of_binary)

binary_period, SMBH_period, perturber_period = get_binary_period(sim), get_binary_SMBH_period(sim), get_perturber_period(sim)
sim.dt = 0.05 * binary_period

rebx = reboundx.Extras(sim)
gr_radiation = rebx.load_force("gr_radiation")
rebx.add_force(gr_radiation)
gr_radiation.params["c"] = c
gr_radiation.params["gr_rad_part1"] = 1
gr_radiation.params["gr_rad_part2"] = 2

spin_ode = sim.create_ode(length = 3, needs_nbody = True)

global BBH_2_L_hat
if mode == "initial_spin_aligned_with_L_of_BBH2":
    spin_ode.y[0], spin_ode.y[1], spin_ode.y[2] = BBH_2_L_hat[0], BBH_2_L_hat[1], BBH_2_L_hat[2]
elif mode == "initial_spin_aligned_with_L_of_Binary":
    spin_ode.y[0], spin_ode.y[1], spin_ode.y[2] = 0, 0, 1

sim.integrator = "BS"

def spin_ode_update(ode, ds_dt, s_hat, t):

    s_hat = np.array([s_hat[0], s_hat[1], s_hat[2]])

    global spin
    spin[0] = s_hat[0]
    spin[1] = s_hat[1]
    spin[2] = s_hat[2]

    BBH_1, BBH_2 = sim.particles["BBH_1"], sim.particles["BBH_2"]
    r = np.array(BBH_2.xyz) - np.array(BBH_1.xyz)
    v = np.array(BBH_2.vxyz) - np.array(BBH_1.vxyz)
    L_hat = normalize(np.cross(r, BBH_2.m * v))

    global BBH_2_L_hat
    BBH_2_L_hat[0] = L_hat[0]
    BBH_2_L_hat[1] = L_hat[1]
    BBH_2_L_hat[2] = L_hat[2]

    binary_orbit = BBH_2.calculate_orbit(primary = BBH_1)

    n = binary_orbit.n
    G = sim.G
    a2 = binary_orbit.a
    mu = BBH_1.m * BBH_2.m / (BBH_1.m + BBH_2.m)
    m2 = BBH_2.m

    #========================

    # n2 = 2*np.pi/get_binary_SMBH_period(sim)
    # m3 = sim.particles["SMBH"].m
    # m1 = BBH_1.m
    # a3 = np.cbrt((G * m3 * get_binary_SMBH_period(sim)**2/(4 * np.pi**2)))
    #
    # de_sitter_magnitude = ( 3 * G * n *(m2 + mu/3) ) / (2 * c**2 * a2)
    #
    # l_precession_magnitude = ( 3 * n2 * m3 * a2**3) / ( 4*(m1 + m2) * a3**3)
    #
    # print(de_sitter_magnitude / l_precession_magnitude)

    #==========================


    ds_dt_vector = ((3 * G * n * (m2 + mu/3)) / (2 * c**2 * a2)) * np.cross(L_hat, s_hat)

    ds_dt[0] = ds_dt_vector[0]
    ds_dt[1] = ds_dt_vector[1]
    ds_dt[2] = ds_dt_vector[2]


spin_ode.derivatives = spin_ode_update


initialize_data_collection()
#####################################################################################

while sim.t <= 10**5 * SMBH_period:

    check_and_assign_minimums(sim, outcome_record)

    global total_time_steps_completed
    if total_time_steps_completed % 100 == 0:
        save_data(sim)

    check_binary_bound(sim, outcome_record)
    check_for_collisions(sim, w, outcome_record)

    total_time_steps_completed += 1
    sim.step()

#####################################################################################

outcome_record["Result"] = "Ended with no mergers and bound binary"
dump_record(outcome_record)


