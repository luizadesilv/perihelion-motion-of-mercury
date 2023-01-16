import pandas as pd
from vpython import *
import seaborn as sns


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


# Calculate the angle between the perihelions for each pair of successive turns
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
    dt = 2 * vM0 / c_a / 20

    # Define the initial coordinates; M = mercury, S = Sun
    M = sphere(pos=vec_rM0, radius=0.5, color=color.red)
    S = sphere(pos=vector(0, 0, 0), radius=1.5, color=color.yellow)

    # And the initial velocities
    M.velocity = vec_vM0
    S.velocity = vector(0, 0, 0)

    # Add a visible trajectory to Mercury
    M.trajectory = curve(color=color.white)

    # Extract multiple position sof the perihelion

    alpha = [0, 1.e4, 2.e4, 3.e4, 4.e4, 5.e4, 6.e4, 7.e4, 8.e4, 9.e4]
    beta = 0

    # alpha = 0
    # beta = [0, 1.e4, 2.e4, 3.e4, 4.e4, 5.e4, 6.e4, 7.e4, 8.e4, 9.e4]

    vec_rM_last = vec_rM0
    max_i = 10
    sum_angle = 0
    list_perih = []
    sum_angle_list = [0]
    i = 0
    while i < max_i:
        vec_rM_before_last = vec_rM_last
        vec_rM_last = vector(M.pos)
        rate(500)
        M.trajectory.append(pos=M.pos)
        M.pos, M.velocity = evolve_mercury(M.pos, M.velocity, alpha[i],  beta)
        if (vec_rM_last.mag < M.pos.mag) & (vec_rM_last.mag < vec_rM_before_last.mag):
            list_perih.append(vec_rM_last)
            i = i + 1
            if i > 1:
                print(list_perih)
                sphere(color=color.green, radius=0.2, pos=vec_rM_last)
                print(f'Turn: {i}, perihelion growth: {angle_between(list_perih[-2], list_perih[-1])}')
                sum_angle = sum_angle + angle_between(list_perih[-2], list_perih[-1])
                sum_angle_list.append(sum_angle / (max_i - 1))
                print(sum_angle_list)
    # Display the average
    print(f'Average perihelion growth: {sum_angle / (len(list_perih)-1 ) * 3. / 3600 * 4.15 * 100}')

    # plot image

    print(len(alpha), len(sum_angle_list))
    data = pd.DataFrame({"alpha": alpha, "perihelion growth": sum_angle_list})
    plot = sns.lmplot(x="alpha", y="perihelion growth", data=data)
    plot.savefig('alpha.jpg')

    # print(len(beta), len(sum_angle_list))
    # data = pd.DataFrame({"beta": beta, "perihelion growth": sum_angle_list})
    # plot = sns.lmplot(x="beta", y="perihelion growth", data=data)
    # plot.savefig('beta.jpg')
