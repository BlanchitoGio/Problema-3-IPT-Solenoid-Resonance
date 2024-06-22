"""
MISURA DELLA RISPOSTA IN FREQUENZA (ampiezza fissa, varia la frequenza)
"""

from vpython import *
from random import uniform
from math import sin, cos
from numpy import sign
from helper_functions import *

#outfile=open('spettro.dat','w')
scene = canvas()

# Create two instances of a graph object
frequency_graph = graph(xtitle="frequency", ytitle="amplitude")
position_graph = graph(title='Y-Position vs. Time', xtitle='Time', ytitle='Y-Position', fast=False)

position_curve = gdots(color=color.blue, graph=position_graph)
frequency_curve = gcurve(color=color.red, graph=frequency_graph)

# Turns on (1) and off (0) gravity and current
GRAV = 1

# system parameters
N = 27                  # total number of coils of spring
k_tot = 0.32            # total elastic const. of spring
m_tot = 0.0072          # total mass of spring
rest_length = 0         # rest length of whole spring
spring_rad = 0.055      # radius of spring
damping_const = 0.001    # damping constant (> 0)

g = vector(0, -9.81*GRAV, 0)

# derived system parameters - DO NOT TOUCH THESE
coil_k = N * k_tot
coil_mass = m_tot / N
coil_len = rest_length / (N-1)
mass_radius = spring_rad * 0.5

# Arrays containing masses and springs
masses = []
springs = []

initial_positions = [0.0, -0.012167881944444264, -0.024032986111110752, -0.03559531249999946, -0.046854861111110396, -0.057811631944443546, -0.06846562499999892, -0.07881684027777651, -0.08886527777777634, -0.09861093749999837, -0.10805381944444263, -0.11719392361110911, -0.12603124999999782, -0.13456579861110876, -0.14279756944444194, -0.15072656249999736, -0.158352777777775, -0.1656762152777749, -0.17269687499999703, -0.17941475694444137, -0.18582986111110794, -0.19194218749999675, -0.1977517361111078, -0.2032585069444411, -0.2084624999999966, -0.21336371527777437, -0.21796215277777437]

# Generate masses and springs
for i in range(N):
    newMass = sphere(radius=mass_radius, color=color.red)
    newMass.pos = vector(0, initial_positions[i], 0)  # EQUILIBRIUM WITH GRAVITY
    # newMass.pos = vector(0, -i*coil_len, 0)   # EQUILIBRIUM WITHOUT GRAVITY
    newMass.mass = coil_mass
    newMass.force = vector(0, 0, 0)
    newMass.vel = vector(0, 0, 0)
    newMass.acc = vector(0, 0, 0)

    masses.append(newMass)

    if i != 0:
        newSpring = helix(pos=masses[i-1].pos, axis=masses[i].pos-masses[i-1].pos, radius=spring_rad, coils=1)
        newSpring.k = coil_k
        newSpring.restLen = coil_len
        springs.append(newSpring)

masses[N-1].color = color.blue
masses[N-1].mass = 0.00405


#pulses = [0.2, 0.3, 0.4, 0.45, 0.5, 0.9, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3]        # current frequencies (Hz) converted to (rad/s)
pulses = [0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.53,0.535,0.54,0.545,0.55,0.56,0.565,0.57,0.575,0.58,0.59,0.6,0.65,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.775,1.8,1.825,1.85,1.875,1.9,1.925,1.95,1.975,2.0,2.025,2.05,2.075,2.1,2.125,2.15,2.175,2.2,2.3,2.4,2.5]
for i in range(len(pulses)):
    pulses[i] = pulses[i]*2*pi

curr_ampl = 5                           # current amplitude (A)

mechanical_oscillations = []
temp = []                       # temporary array to contain positions of last mass

for curr_pulse in pulses:

    # RESET PARAMETERS
    for i in range(N):
        masses[i].pos.y = initial_positions[i]
        masses[i].vel.y = 0
    temp = []

    timespan = 5 + 3 * (2*pi)/curr_pulse

    # WHILE LOOP
    dt = 0.005
    time = 0

    while time <= timespan:
        rate(200)

        # CURRENT
        current = curr_ampl * sin(curr_pulse*time)

        # Calculate the force experienced by each mass
        for i in range(1, N):

            # damping force
            damping_force = - masses[i].vel * damping_const

            # lorentz force
            lorentz_force = vector(0, calc_lorentz(masses, i, current, spring_rad), 0)   # MOST ACCURATE MODEL
            # lorentz_force = vector(0, calc_lorentz(masses, i, current, spring_rad, linear="True"), 0)   # MOST ACCURATE MODEL

            # LAST MASS
            if i == N-1:

                # distance from the upper mass
                upper_dist = masses[N-2].pos.y - masses[N-1].pos.y
                upper_spring_force = springs[N-2].k * (upper_dist - springs[N-2].restLen) * vector(0, 1, 0)

                masses[i].force = upper_spring_force + damping_force + masses[i].mass*g + lorentz_force

            # MIDDLE MASSES
            else:
                upper_dist = masses[i-1].pos.y - masses[i].pos.y
                lower_dist = masses[i].pos.y - masses[i+1].pos.y

                lower_spring_force = -springs[i].k * (lower_dist - springs[i].restLen) * vector(0, 1, 0)
                upper_spring_force = springs[i-1].k * (upper_dist - springs[i-1].restLen) * vector(0, 1, 0)

                masses[i].force = upper_spring_force + lower_spring_force + damping_force + masses[i].mass*g + lorentz_force

        # Evaluate displacement experienced by each mass
        for i in range(N):

            masses[i].acc = masses[i].force / masses[i].mass
            masses[i].vel += masses[i].acc * dt
            masses[i].pos += masses[i].vel * dt

            if i != N-1:
                springs[i].pos = masses[i].pos
                springs[i].axis = masses[i+1].pos - masses[i].pos

        if time > timespan/2:
            temp.append(masses[N-1].pos.y)

        # Plots
        position_curve.plot(time, masses[N-1].pos.y)

        time += dt

    osc = max(temp) - min(temp)
    mechanical_oscillations.append(max(temp) - min(temp))
    
    with open('frequenze', 'a') as file:
        file.write(('%f\n')%((curr_pulse)/(2*pi)))
    with open('ampiezze', 'a') as file:
        file.write(('%f\n')%(osc))
    
    #outfile.write(('%f %f\n')%(curr_pulse/(2*pi),osc))

    frequency_curve.plot(curr_pulse/(2*pi), osc)

exit()