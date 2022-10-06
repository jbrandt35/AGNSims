from Tools import *
from math import pi

perturber_a = 2.0
binary_separation = 0.1

#####################################################################################

sim = create_simulation()
w = populate_simulation(sim, perturber_a = perturber_a, binary_separation = binary_separation, binary_inc = pi, randomize_M = True)

binary_period, SMBH_period, perturber_period = get_binary_period(sim), get_binary_SMBH_period(sim), get_perturber_period(sim)
sim.dt = 0.05 * binary_period

rebx = reboundx.Extras(sim)
gr_radiation = rebx.load_force("gr_radiation")
rebx.add_force(gr_radiation)
gr_radiation.params["c"] = c
gr_radiation.params["gr_rad_part1"] = 1
gr_radiation.params["gr_rad_part2"] = 2

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


