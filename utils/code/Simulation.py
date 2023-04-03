from Tools import *
import numpy as np
import reboundx
import config

############################   Settings    ###########################################

perturber_a = PERTSEPARATION
binary_separation = BINSEPARATION

config.mode = "initial_spin_aligned_with_L_of_Binary"

include_drag_force = True
tau_drag = 5e6

include_trap_force = False
tau_trap = 5e5

include_rebound_migration = False
tau_mig = 1e5

##############################    Initializing   #####################################

sim = create_simulation()
w = populate_simulation(sim, perturber_a = perturber_a, binary_separation = binary_separation, randomize_M = True, randomize_binary_inc = True, randomize_eccentricities = True)

binary_period, SMBH_period, perturber_period = get_binary_period(sim), get_binary_SMBH_period(sim), get_perturber_period(sim)
sim.dt = 0.05 * binary_period

##############################    Post-Newtonian Effects   ###########################

rebx = reboundx.Extras(sim)

#2.5 PN
gr_radiation = rebx.load_force("gr_radiation")
rebx.add_force(gr_radiation)
gr_radiation.params["c"] = c
gr_radiation.params["gr_rad_part1"] = 1
gr_radiation.params["gr_rad_part2"] = 2

#1 PN
# gr_full = rebx.load_force("gr_full")
# rebx.add_force(gr_full)
# gr_full.params["c"] = c

##############################    Tracking Spin   #####################################

sim.integrator = "BS"

spin_ode = sim.create_ode(length = 3, needs_nbody = True)
spin_ode.y[0], spin_ode.y[1], spin_ode.y[2] = config.spin.tolist()

def spin_ode_update(ode, ds_dt, s_hat, t):

    unit_spin_vector = np.array([s_hat[0], s_hat[1], s_hat[2]])
    config.spin = np.copy(unit_spin_vector)

    BBH_1, BBH_2 = sim.particles["BBH_1"], sim.particles["BBH_2"]
    L_hat = unit_mutual_orbital_angular_momentum(BBH_2, BBH_1)

    binary_orbit = BBH_2.calculate_orbit(primary = BBH_1)

    n = binary_orbit.n
    G = sim.G
    a2 = binary_orbit.a
    mu = BBH_1.m * BBH_2.m / (BBH_1.m + BBH_2.m)
    m2 = BBH_2.m

    ds_dt_vector = ((3 * G * n * (m2 + mu/3)) / (2 * c**2 * a2)) * np.cross(L_hat, unit_spin_vector)

    ds_dt[0], ds_dt[1], ds_dt[2] = ds_dt_vector.tolist()

spin_ode.derivatives = spin_ode_update


##############################   Drag Effects   #####################################

def dragForce(reb_sim, rebx_force, particles, N):
    sim = reb_sim.contents
    dragforce = rebx_force.contents
    vxyz = particles[3].vxyz
    r_0 = particles[0].xyz
    r_pert = particles[3].xyz
    tau_drag = dragforce.params["c"]
    ax, ay, az = calc_drag_force(tau_drag*SMBH_period, vxyz, sim.G, particles[0].m, r_0, r_pert)
    particles[3].ax += ax
    particles[3].ay += ay
    particles[3].az += az

def trapForce(reb_sim, rebx_force, particles, N):
    sim = reb_sim.contents
    trapforce = rebx_force.contents
    r_SMBH = particles[0].xyz
    r_pert = particles[3].xyz
    tau_trap = trapforce.params["c"]
    ax, ay, az = calc_trap_force(tau_trap*SMBH_period, sim.G, particles[0].m, 1000*(sim.G * m0 / c ** 2), r_pert, r_SMBH)
    particles[3].ax += ax
    particles[3].ay += ay
    particles[3].az += az

if include_drag_force:
    myforce = rebx.create_force("drag")
    myforce.force_type = "vel"
    myforce.update_accelerations = dragForce
    rebx.add_force(myforce)
    myforce.params["c"] = tau_drag

if include_trap_force:
    myforce = rebx.create_force("trap")
    myforce.force_type = "vel"
    myforce.update_accelerations = trapForce
    rebx.add_force(myforce)
    myforce.params["c"] = tau_trap

if include_rebound_migration:
    mig = rebx.load_force("modify_orbits_forces")
    rebx.add_force(mig)
    sim.particles['perturber'].params['tau_a'] = -tau_mig

#####################################################################################

initialize_data_collection()
save_initial_data(sim)

while sim.t <= 10**5 * SMBH_period:

    check_and_assign_minimums(sim)

    if config.total_time_steps_completed % 100 == 0:
        save_data(sim)

    check_bound(sim)
    check_for_collisions(sim, w)

    config.total_time_steps_completed += 1
    sim.step()

#####################################################################################

if is_bound(sim.particles["BBH_1"], sim.particles["BBH_2"]):
    config.outcome_record["Result"] = "Ended with no mergers, bound binary, and all particles within limit to SMBH"
else:
    config.outcome_record["Result"] = "Ended with no mergers, unbound binary, and all particles within limit to SMBH"

save_final_data(sim)
dump_record()


