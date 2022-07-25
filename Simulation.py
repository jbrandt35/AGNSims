import numpy as np
from Tools import *

perturber_a = 4
binary_separation = 0.5

#####################################################################################

def migration_force(reb_sim):
    sim.particles["perturber"].ax -= sim.particles["perturber"].vx/tau
    sim.particles["perturber"].ay -= sim.particles["perturber"].vy/tau
    sim.particles["perturber"].az -= sim.particles["perturber"].vz/tau

def heartbeat(reb_sim):
    check_and_assign_minimums(sim, outcome_record)

    if total_time_steps_completed % 10000 == 0:
        save_plotting_data(sim)
        construct_plots()

    check_binary_bound(sim, outcome_record)
    check_for_collisions(sim, w, outcome_record)


#####################################################################################

sim = create_simulation()
w = populate_simulation(sim, perturber_a = perturber_a, binary_separation = binary_separation)
sim.heartbeat = heartbeat

binary_period, SMBH_period = get_binary_period(sim), get_SMBH_period(sim)
sim.dt = 0.05 * binary_period

tau = 1000.
sim.additional_forces = migration_force
sim.force_is_velocity_dependent = 1

sim.integrate(10**5 * SMBH_period)

#####################################################################################

outcome_record["Result"] = "Ended with no mergers and bound binary"
dump_record(outcome_record)


