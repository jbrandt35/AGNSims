from Tools import *

perturber_a = 4
binary_separation = 0.3

#####################################################################################

def migration_force(reb_sim):
    sim.particles["perturber"].ax -= sim.particles["perturber"].vx/tau
    sim.particles["perturber"].ay -= sim.particles["perturber"].vy/tau
    sim.particles["perturber"].az -= sim.particles["perturber"].vz/tau


#####################################################################################

sim = create_simulation()
w = populate_simulation(sim, perturber_a = perturber_a, binary_separation = binary_separation)

binary_period, SMBH_period, perturber_period = get_binary_period(sim), get_binary_SMBH_period(sim), get_perturber_period(sim)
sim.dt = 0.05 * binary_period

tau = 10**5 * perturber_period
sim.additional_forces = migration_force
sim.force_is_velocity_dependent = 1

initialize_data_collection()
#####################################################################################

while sim.t <= 10**3 * SMBH_period:

    check_and_assign_minimums(sim, outcome_record)

    global total_time_steps_completed
    if total_time_steps_completed % 1 == 0:
        save_data(sim)

    check_binary_bound(sim, outcome_record)
    check_for_collisions(sim, w, outcome_record)

    total_time_steps_completed += 1
    sim.step()

#####################################################################################

outcome_record["Result"] = "Ended with no mergers and bound binary"
dump_record(outcome_record)


