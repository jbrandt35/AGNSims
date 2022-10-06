from Tools import *

perturber_a = 2.0
binary_separation = 0.1

#####################################################################################

sim = create_simulation()
w = populate_simulation(sim, perturber_a = perturber_a, binary_separation = binary_separation, randomize_M = True)

binary_period, SMBH_period, perturber_period = get_binary_period(sim), get_binary_SMBH_period(sim), get_perturber_period(sim)
sim.dt = 0.05 * binary_period

rebx = reboundx.Extras(sim)
gr_radiation = rebx.load_force("gr_radiation")
rebx.add_force(gr_radiation)
gr_radiation.params["c"] = c
gr_radiation.params["gr_rad_part1"] = 1
gr_radiation.params["gr_rad_part2"] = 2

spin_ode = sim.create_ode(length = 3, needs_nbody = True)

spin_ode.y[0], spin_ode.y[1], spin_ode.y[2] = 0, 0, 1

def spin_ode_update(ode, ds_dt, s_hat, t):

    s_hat = np.array([s_hat[0], s_hat[1], s_hat[2]])

    BBH_1, BBH_2 = sim.particles["BBH_1"], sim.particles["BBH_2"]

    binary_orbit = BBH_2.calculate_orbit(primary = BBH_1)

    r = np.array(BBH_2.xyz) - np.array(BBH_1.xyz)
    v = np.array(BBH_2.vxyz) - np.array(BBH_1.vxyz)

    L_hat = normalize(np.cross(r, BBH_2.m * v))

    n = binary_orbit.n
    G = sim.G
    a2 = binary_orbit.a
    mu = BBH_1.m * BBH_2.m / (BBH_1.m + BBH_2.m)
    m2 = BBH_2.m

    ds_dt_vector = ((3 * G * n * (m2 + mu/3)) / (2 * c**2 * a2)) * np.cross(L_hat, s_hat)

    ds_dt[0], ds_dt[1], ds_dt[2] = ds_dt_vector


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


