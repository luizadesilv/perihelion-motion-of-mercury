import numpy as np
from matplotlib import pyplot as plt
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

    # Calculate energy error
    error = ((0.5 * vec_vM_new.mag ** 2) - (GM / vec_rM_new.mag) - (0.5 * (vM0 ** 2)) + (GM / rM0)) / (
                0.5 * (vM0 ** 2) - (GM / rM0))

    return vec_rM_new, vec_vM_new, abs(error)


# euler-explicito
def evolve_euler_explicito(vec_rM_old, vec_vM_old, alpha, beta):
    temp = 1 + alpha * rS / vec_rM_old.mag + beta * rL2 / vec_rM_old.mag ** 2
    aMS = c_a * temp / vec_rM_old.mag ** 2
    vec_aMS = - aMS * (vec_rM_old / vec_rM_old.mag)
    vec_rM_new = vec_rM_old + vec_vM_old * dt
    vec_vM_new = vec_vM_old + vec_aMS * dt

    # Calculate energy error
    error = ((0.5 * vec_vM_new.mag ** 2) - (GM / vec_rM_new.mag) - (0.5 * (vM0 ** 2)) + (GM / rM0)) / (
                0.5 * (vM0 ** 2) - (GM / rM0))

    return vec_rM_new, vec_vM_new, abs(error)


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

    # Calculate energy error
    error = ((0.5 * vec_vM_new.mag ** 2) - (GM / vec_rM_new.mag) - (0.5 * (vM0 ** 2)) + (GM / rM0)) / (
                0.5 * (vM0 ** 2) - (GM / rM0))

    return vec_rM_new, vec_vM_new, abs(error)


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

    # Calculate energy error
    error = ((0.5 * vec_vM_new.mag ** 2) - (GM / vec_rM_new.mag) - (0.5 * (vM0 ** 2)) + (GM / rM0)) / (0.5 * (vM0 ** 2) - (GM / rM0))

    return vec_rM_new, vec_vM_new, abs(error)


if __name__ == '__main__':
    # Definition of physical parameters
    rM0 = 4.60  # Initial radius of Mercury orbit,in units of R0
    vM0 = 5.10e-1  # Initial orbital speed of Mercury, in units of R0/T0
    c_a = 9.90e-1  # Base acceleration of Mercury, in units of R0**3/T0**2
    TM = 8.80e+1  # Orbit period of Mercury
    rS = 2.95e-7  # Schwarzschild radius of Sun, in units of R0
    rL2 = 8.19e-7  # Specific angular momentum, in units of R0**2
    GM = 0.9906  # The same as c_a
    print(GM)

    vec_rM0 = vector(0, rM0, 0)
    vec_vM0 = vector(vM0, 0, 0)

    # Definition of the time step
    dt = 2 * vM0 / c_a / 20

    # Define the initial coordinates; M = mercury, S = Sun
    # Euler-cromer
    M1 = sphere(pos=vector(0, rM0, 0), radius=0.5, color=color.red)
    S1 = sphere(pos=vector(0, 0, 0), radius=1.5, color=color.yellow)
    # Euler-explicito
    M2 = sphere(pos=vector(0, rM0, 0), radius=0.5, color=color.green)
    S2 = sphere(pos=vector(0, 0, 0), radius=1.5, color=color.yellow)
    # Velocity-verlet
    M3 = sphere(pos=vector(0, rM0, 0), radius=0.5, color=color.blue)
    S3 = sphere(pos=vector(0, 0, 0), radius=1.5, color=color.yellow)
    # Runge-kutta2
    M4 = sphere(pos=vector(0, rM0, 0), radius=0.5, color=color.magenta)
    S4 = sphere(pos=vector(0, 0, 0), radius=1.5, color=color.yellow)

    # And the initial velocities
    # Euler-cromer
    M1.velocity = vector(vM0, 0, 0)
    S1.velocity = vector(0, 0, 0)
    # Euler-explicito
    M2.velocity = vector(vM0, 0, 0)
    S2.velocity = vector(0, 0, 0)
    # Velocity-verlet
    M3.velocity = vector(vM0, 0, 0)
    S3.velocity = vector(0, 0, 0)
    # Runge-kutta2
    M4.velocity = vector(vM0, 0, 0)
    S4.velocity = vector(0, 0, 0)

    # Add a visible trajectory to mercury
    # Euler-cromer
    M1.trajectory = curve(color=color.white)
    # Euler-explicito
    M2.trajectory = curve(color=color.white)
    # Velocity-verlet
    M3.trajectory = curve(color=color.white)
    # Runge-kutta2
    M4.trajectory = curve(color=color.white)

    t = 0.0
    alpha = 0.0
    beta = 0.0

    energy_error = [[], [], [], []]

    # Execute the loop as long as t < 2*TM
    while t < 2 * TM:
        # Set the frame rate (you can choose a higher rate to accelerate the program)
        rate(500)

        # Update the drawn trajectory with the current position
        # Euler-cromer
        M1.trajectory.append(pos=M1.pos)
        # Euler-explicito
        M2.trajectory.append(pos=M2.pos)
        # Velocity-verlet
        M3.trajectory.append(pos=M3.pos)
        # Runge-kutta2
        M4.trajectory.append(pos=M4.pos)

        # Update velocity and position
        # Euler-cromer
        M1.pos, M1.velocity, error1 = evolve_mercury(M1.pos, M1.velocity, alpha, beta)
        energy_error[0].append(error1)
        # Euler-explicito
        M2.pos, M2.velocity, error2 = evolve_euler_explicito(M2.pos, M2.velocity, alpha, beta)
        energy_error[1].append(error2)
        # Velocity-verlet
        M3.pos, M3.velocity, error3 = evolve_velocity_verlet(M3.pos, M3.velocity, alpha, beta)
        energy_error[2].append(error3)
        # Runge-kutta2
        M4.pos, M4.velocity, error4 = evolve_runge_kutta2(M4.pos, M4.velocity, alpha, beta)
        energy_error[3].append(error4)

        # Advance time by one step
        t = t + dt

    # Verify the energy error values
    print("1", energy_error[0])
    print("2", energy_error[1])
    print("3", energy_error[2])
    print("4", energy_error[3])

    # Plot image energy error x time relates to two turns
    x = np.arange(0, len(energy_error[0])*dt, dt)
    plt.plot(x, energy_error[0], color="b", label="Euler Cromer")
    plt.plot(x, energy_error[1], color="r", label="Euler Explícito")
    plt.plot(x, energy_error[2], color="g", label="Velocity Verlet")
    plt.plot(x, energy_error[3], color="y", label="Runge Kutta 2")
    plt.yscale("log")
    plt.xlabel("Time")
    plt.ylabel("Energy Error")
    plt.title("Energy error")
    plt.legend()
    plt.show()
    plt.savefig(r"energy_error.png")
