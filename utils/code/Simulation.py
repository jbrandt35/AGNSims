from Tools import *
import numpy as np

############################   Settings    ###########################################

perturber_a = PERTSEPARATION
binary_separation = BINSEPARATION

mode = "initial_spin_aligned_with_L_of_Binary"

include_drag_force = True
tau_drag = 5e6

include_trap_force = False
tau_trap = 5e5

include_rebound_migration = False
tau_mig = 1e5

##############################    Initializing   #####################################

global inclination_of_binary

sim = create_simulation()
w = populate_simulation(sim, mode, perturber_a = perturber_a, binary_separation = binary_separation, ignore_perturber = False, randomize_M = True, binary_inc = inclination_of_binary)

binary_period, SMBH_period, perturber_period = get_binary_period(sim), get_binary_SMBH_period(sim), get_perturber_period(sim)
sim.dt = 0.05 * binary_period

rebx = reboundx.Extras(sim)

#2.5 PN
gr_radiation = rebx.load_force("gr_radiation")
rebx.add_force(gr_radiation)
gr_radiation.params["c"] = c
gr_radiation.params["gr_rad_part1"] = 1
gr_radiation.params["gr_rad_part2"] = 2

#1 PN
gr_full = rebx.load_force("gr_full")
rebx.add_force(gr_full)
gr_full.params["c"] = c

##############################    Tracking Spin   #####################################

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

    ds_dt_vector = ((3 * G * n * (m2 + mu/3)) / (2 * c**2 * a2)) * np.cross(L_hat, s_hat)

    ds_dt[0] = ds_dt_vector[0]
    ds_dt[1] = ds_dt_vector[1]
    ds_dt[2] = ds_dt_vector[2]


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

while sim.t <= 10**5 * SMBH_period:

    check_and_assign_minimums(sim, outcome_record)

    global total_time_steps_completed
    if total_time_steps_completed % 100 == 0:
        save_data(sim)

    check_bound(sim, outcome_record)
    check_for_collisions(sim, w, outcome_record)

    total_time_steps_completed += 1
    sim.step()

#####################################################################################

if is_bound(sim.particles[1], sim.particles[2]):
    outcome_record["Result"] = "Ended with no mergers, bound binary, and all particles within limit to SMBH"
else:
    outcome_record["Result"] = "Ended with no mergers, unbound binary, and all particles within limit to SMBH"

save_final_data(outcome_record, sim)
dump_record(outcome_record)


