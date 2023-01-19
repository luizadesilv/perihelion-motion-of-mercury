from vpython import *


# Define the coordinate and velocity update
# euler-cromer
def evolve_mercury(vec_rM_old, vec_vM_old, alpha, beta):
    # Compute the strength of the acceleration
    temp = 1 + alpha * rS / vec_rM_old.mag + beta * rL2 / vec_rM_old.mag ** 2
    aMS = c_a * temp / vec_rM_old.mag ** 2
    # Multiply by the direction
    vec_aMS = - aMS * (vec_rM_old / vec_rM_old.mag)
    # Update velocity vector
    vec_vM_new = vec_vM_old + vec_aMS * dt
    # Update position vector
    vec_rM_new = vec_rM_old + vec_vM_new * dt

    # error = (1 * vec_rM_new ** 2) * 0.5 + (1 * vec_vM_new ** 2) * 0.5
    # print(error)

    return vec_rM_new, vec_vM_new


# euler-explicito
def evolve_euler_explicito(vec_rM_old, vec_vM_old, alpha, beta):
    temp = 1 + alpha * rS / vec_rM_old.mag + beta * rL2 / vec_rM_old.mag ** 2
    aMS = c_a * temp / vec_rM_old.mag ** 2
    vec_aMS = - aMS * (vec_rM_old / vec_rM_old.mag)
    vec_rM_new = vec_rM_old + vec_vM_old * dt
    vec_vM_new = vec_vM_old + vec_aMS * dt

    return vec_rM_new, vec_vM_new


# velocity verlet
def evolve_velocity_verlet(vec_rM_old, vec_vM_old, alpha, beta):
    temp = 1 + alpha * rS / vec_rM_old.mag + beta * rL2 / vec_rM_old.mag ** 2
    aMS = c_a * temp / vec_rM_old.mag ** 2
    vec_aMS = - aMS * (vec_rM_old / vec_rM_old.mag)
    vec_rM_new = vec_rM_old + (vec_vM_old * dt) + (0.5 * vec_aMS * dt ** 2)

    temp1 = 1 + alpha * rS / vec_rM_new.mag + beta * rL2 / vec_rM_new.mag ** 2
    aMS1 = c_a * temp1 / vec_rM_new.mag ** 2
    vec_aMS1 = - aMS1 * (vec_rM_new / vec_rM_new.mag)
    vec_vM_new = vec_vM_old + (0.5 * (vec_aMS + vec_aMS1) * dt)

    return vec_rM_new, vec_vM_new


# Runge-Kutta 2 - Ponto Central
def evolve_runge_kutta2(vec_rM_old, vec_vM_old, alpha, beta):
    temp = 1 + alpha * rS / vec_rM_old.mag + beta * rL2 / vec_rM_old.mag ** 2
    aMS = c_a * temp / vec_rM_old.mag ** 2
    vec_aMS = - aMS * (vec_rM_old / vec_rM_old.mag)

    k1r = vec_vM_old
    k1v = vec_aMS
    k2r = vec_vM_old + k1v * 0.5 * dt

    vec_rM_old2 = vec_rM_old + k1r * 0.5 * dt
    temp2 = 1 + alpha * rS / vec_rM_old2.mag + beta * rL2 / vec_rM_old2.mag ** 2
    aMS2 = c_a * temp2 / (vec_rM_old + k1r * 0.5 * dt).mag ** 2
    k2v = -aMS2 * vec_rM_old2 / vec_rM_old2.mag

    vec_vM_new = vec_vM_old + k2v * dt
    vec_rM_new = vec_rM_old + k2r * dt

    return vec_rM_new, vec_vM_new


if __name__ == '__main__':
    # Definition of physical parameters
    rM0 = 4.60  # Initial radius of Mercury orbit,in units of R0
    vM0 = 5.10e-1  # Initial orbital speed of Mercury, in units of R0/T0
    c_a = 9.90e-1  # Base acceleration of Mercury, in units of R0**3/T0**2
    TM = 8.80e+1  # Orbit period of Mercury
    rS = 2.95e-7  # Schwarzschild radius of Sun, in units of R0
    rL2 = 8.19e-7  # Specific angular momentum, in units of R0**2

    vec_rM0 = vector(0, rM0, 0)
    vec_vM0 = vector(vM0, 0, 0)

    # Definition of the time step
    dt = 2 * vM0 / c_a / 20

    # Define the initial coordinates; M = mercury, S = Sun
    M = sphere(pos=vector(0, rM0, 0), radius=0.5, color=color.red)
    S = sphere(pos=vector(0, 0, 0), radius=1.5, color=color.yellow)

    # And the initial velocities
    M.velocity = vector(vM0, 0, 0)
    S.velocity = vector(0, 0, 0)

    # Add a visible trajectory to mercury
    M.trajectory = curve(color=color.white)

    ##perihelion_motion()
    t = 0.0
    alpha = 0.0
    beta = 0.0

    # Execute the loop as long as t < 2*TM
    while t < 2 * TM:
        # Set the frame rate (you can choose a higher rate to accelerate the program)
        rate(500)
        # Update the drawn trajectory with the current position
        M.trajectory.append(pos=M.pos)

        # Update velocity and position
        # M.pos, M.velocity = evolve_mercury(M.pos, M.velocity, alpha, beta)
        # M.pos, M.velocity = evolve_euler_explicito(M.pos, M.velocity, alpha, beta)
        # M.pos, M.velocity = evolve_velocity_verlet(M.pos, M.velocity, alpha, beta)
        M.pos, M.velocity = evolve_runge_kutta2(M.pos, M.velocity, alpha, beta)

        # Advance time by one step
        t = t + dt

