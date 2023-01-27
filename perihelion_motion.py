import numpy as np
import pandas as pd
from vpython import *
import seaborn as sns
import matplotlib.pyplot as plt


# Euler-Cromer
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

    return vec_rM_new, vec_vM_new


# Calculate the angle between the perihelion for each pair of successive turns
def angle_between(v1, v2):
    return acos(dot(v1, v2) / (v1.mag * v2.mag)) * 180. / pi


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
    dt = 2 * vM0 / c_a / 200

    # Define the initial coordinates; M = mercury, S = Sun
    M = sphere(pos=vec_rM0, radius=0.5, color=color.red)
    S = sphere(pos=vector(0, 0, 0), radius=1.5, color=color.yellow)

    # And the initial velocities
    M.velocity = vec_vM0
    S.velocity = vector(0, 0, 0)

    # Add a visible trajectory to Mercury
    M.trajectory = curve(color=color.white)

    # Extract multiple position sof the perihelion
    beta = [0, 1.e4, 2.e4, 3.e4, 4.e4, 5.e4, 6.e4, 7.e4, 8.e4, 9.e4, 1.e5]
    alpha = 0
    # alpha = [0, 1.e4, 2.e4, 3.e4, 4.e4, 5.e4, 6.e4, 7.e4, 8.e4, 9.e4, 1.e5]
    # beta = 0
    sum_angle_list = []

    # For each beta, run 10 turns
    for j in range(len(beta)):
        vec_rM_last = vec_rM0
        max_i = 10
        sum_angle = 0
        list_perih = []
        i = 0
        while i < max_i:
            vec_rM_before_last = vec_rM_last
            vec_rM_last = vector(M.pos)
            rate(5000)
            M.trajectory.append(pos=M.pos)
            M.pos, M.velocity = evolve_mercury(M.pos, M.velocity, alpha,  beta[j])
            if (vec_rM_last.mag < M.pos.mag) & (vec_rM_last.mag < vec_rM_before_last.mag):
                list_perih.append(vec_rM_last)
                i = i + 1
                if i > 1:
                    sphere(color=color.green, radius=0.2, pos=vec_rM_last)
                    print(f'Turn: {i}, perihelion growth: {angle_between(list_perih[i-2], list_perih[i-1])}')
                    sum_angle = sum_angle + angle_between(list_perih[i-2], list_perih[i-1])

        # Display the average Perihelion Growth for each beta/alpha
        print(f'Perihelion growth: {sum_angle / (len(list_perih)-1)}')
        sum_angle_list.append(sum_angle/9)

    pg_minutes = (sum_angle_list[-1] / 1.e5 * 3) * 60 * 60
    pg_100years = pg_minutes * 415
    print(f'Perihelion growth: {pg_minutes} minutes')
    print(f'Perihelion motion per 100 earth years: {pg_100years}')

    # Plot image δΘ(α = 0, β) x β and δΘ(α, β = 0) x α
    data = pd.DataFrame({"beta": beta, "perihelion growth": sum_angle_list})
    plot = sns.lmplot(x="beta", y="perihelion growth", data=data)
    plt.legend([f'δΘ(α = 0, β) = {sum_angle_list[-1]} x β'])
    plot.savefig(r'C:\Users\Luiza D. da Silva\Documents\GitHub\Prova1MetCompB\beta_linear_reg.jpg')

    print('End')
