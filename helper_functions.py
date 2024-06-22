from vpython import *
from math import *

mu_0 = 4 * pi * (10 ** -7)


# Returns the unit vector between from_vector and to_vector
def unit_vect(from_vector, to_vector):
    direction_vector = to_vector - from_vector
    return direction_vector / mag(direction_vector)


# Calculates the center of mass of the system of masses
def calc_CM(massList):
    pos_sum = 0
    tot_mass = 0

    for whatMass in massList:
        tot_mass += whatMass.mass
        pos_sum += whatMass.pos.y * whatMass.mass

    return pos_sum / tot_mass


# a is a scalar and b is a vector representing an interval. Returns 1 if a is contained in interval, 0 otherwise
def kronecker(a, b):
    if b[0] <= a <= b[1]:
        return 1
    else:
        return 0


# equilibrium position of mass_num-mass in the presence of gravity and all masses are the same - DOES NOT WORK
def get_eq_grav(tot_len, tot_mass, tot_k, N, mass_num):
    coil_k = N * tot_k
    coil_mass = tot_mass / N
    coil_len = tot_len / (N - 1)
    g = 9.81

    sub = 0
    for i in range(mass_num):
        sub += i

    eq_pos = vector(0, - (mass_num + 1) * coil_len - mag(g) * coil_mass / coil_k * ((mass_num + 1) * N - sub), 0)

    return eq_pos


# calculates the lorentz force on a mass given the whole spring. mass_num is the index of the specific mass
def calc_lorentz(mass_list, mass_num, current=0, r=0, linear=False):

    alpha = mu_0 * (current ** 2) * r

    if linear:
        return 0

    else:

        # Calculate the sum of the forces
        sigma = 0
        for i in range(len(mass_list)):
            if i != mass_num:
                sigma += 1 / (mass_list[i].pos.y - mass_list[mass_num].pos.y)

        # Total lorentz force acting on the coil
        F_lorentz = alpha * sigma

    return F_lorentz


"""
def calc_lorentz(mass_list, mass_num, current=0, r=1.0):

    alpha = mu_0 * (current ** 2) * r

    # Calculate the sum of the forces
    sigma = 0
    for i in [mass_num-2, mass_num-1, mass_num+1, mass_num+2]:
        if 0 < i < len(mass_list):
            sigma += 1 / (mass_list[i].pos.y - mass_list[mass_num].pos.y)

    # Total lorentz force acting on the coil
    F_lorentz = alpha * sigma

    return F_lorentz
"""