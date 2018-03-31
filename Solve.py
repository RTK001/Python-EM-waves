
'''
A program to Solve the Electromagnetic Wave Equation in a closed domain.
'''


import numpy as np
import matplotlib as mp
from matplotlib.pyplot import plot
from matplotlib import animation
from Wave_Eqn_Class_Defs import *



# Extract input file parameters

# Load the sim file
input_file_name = "Input.npy"
sim = np.load(input_file_name).item(0)



# Create input function
# Creates a sine wave of the input function
input_function = np.array([np.sin(2*np.pi* sim.InputFrequency.frequency * sim.timestep_size*t) for t in range(sim.InputFrequency.input_duration)], np.double)
input_function.resize(sim.no_of_timesteps)                  # Resize the array to match the simulation size, adding zeros where appropriate to indicate no signal.


# Test

# Create Array
A = Element_Array_builder(sim)

# Setup Boundary conditions - input pulse
A[0, 1 ].SetE(0 ,input_function)


# Solve
for t in range(3,sim.no_of_timesteps):
    for index, elem in np.ndenumerate(A):
        A[index].Calculate(t)


E_field = np.array([elem.E[:] for index, elem in np.ndenumerate(A)])
E_field = E_field.reshape(*sim.no_of_elements, sim.no_of_timesteps, sim.no_of_dimensions)
B_field = np.array([elem.B[:] for index, elem in np.ndenumerate(A)])
B_field = B_field.reshape(*sim.no_of_elements, sim.no_of_timesteps, sim.no_of_dimensions)


np.save("Res", np.array([E_field, B_field]))

