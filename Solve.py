
'''
A program to Solve the Electromagnetic Wave Equation in a closed domain.
'''




import numpy as np
import matplotlib as mp
from matplotlib.pyplot import plot
from matplotlib import animation
from Wave_Eqn_Class_Defs import *

from Pre_Processor import Simulation


# Used for debugging of numerical errors
import pdb
np.seterr(all='raise')



# Extract input file parameters

# Load the sim file
input_file_name = "Input.npy"
sim = np.load(input_file_name).item(0)




# Create input function
input_function = np.array([np.sin(2*np.pi* sim.InputFrequency.frequency * sim.timestep_size*t) for t in range(sim.InputFrequency.input_duration)], np.double)
# Creates a sine wave of the input function
input_function.resize(sim.no_of_timesteps)                  # Resize the array to match the simulation size, adding zeros where appropriate to indicate no signal.



# pre-compute constant variables to reduce solve time
speed_of_light_squared_x_timestep_sq = (sim.speed_of_light**2) * (sim.timestep_size**2)
# pre-calculated alpha woithout the element length term, which is added later in the individual element. 




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
B_field = np.array([elem.B[:] for index, elem in np.ndenumerate(A)])

np.save("Res", np.array([E_field, B_field]))

