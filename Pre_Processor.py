
'''
A program to Solve the Electromagnetic Wave Equation in a closed domain.
'''


import numpy as np
import matplotlib.pyplot as plt
from Wave_Eqn_Class_Defs import *


            

# Domain parameters
Domain_length = np.array([0.5, 0.5])
# An n dimensional array of the domain lengths in Meters.
# Later used to determine number of dimensions.
Simulation_Duration = 1* 10** -8   		# Seconds

dom = Domain(Domain_length, Simulation_Duration)



# Input signal parameters
frequency = 2*10**9     # in Hz
input_duration = 13     # Measured in timesteps

in_f = InputFrequency(frequency, input_duration)



# Mesh and timestep sizing
Elements_per_Wavelength = 14        # Based on input wavelength
TimeSteps_per_Wavelength = 25

grid = GridSettings(Elements_per_Wavelength, TimeSteps_per_Wavelength)


sim = Simulation(dom, in_f, grid)


file = "Input"

sim.save(file)



