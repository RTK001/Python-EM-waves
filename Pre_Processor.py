
'''
A program to Solve the Electromagnetic Wave Equation in a closed domain.
'''


import numpy as np


# Used for debugging of numerical errors
import pdb
np.seterr(all='raise')


class Domain:
    '''
    Takes and stores the Domain Length (m) and the simulation duration (s).
    The Domain Length should be a numpy array; the length of the Domain in each dimension(x,y,z,etc..)
    '''
    def __init__(self, Length, Duration):
        self.Domain_length = Length
        self.Simulation_Duration = Duration


class InputFrequency:
    '''
    A class to create an input frequency for the simulation. It requires the frequency, and the Duration in timesteps.
    '''
    def __init__(self, frequency, Duration):
        self.frequency = frequency
        self.input_duration = Duration

class GridSettings:
    ''' A Class to contain the settings for the Grid to be used. As these will mostly be numeric parameters, they have default values.'''
    def __init__(self, Elements_Per_WaveLength = 14, Elements_Per_Timestep = 25):
        self.Elements_per_Wavelength = Elements_Per_WaveLength
        self.TimeSteps_per_Wavelength = Elements_Per_Timestep


class Simulation():
    ''' A container class to contain the simulation subclasses '''
    def __init__(self, domain, input_freq, settings):
        self.Domain = domain
        self.InputFrequency = input_freq
        self.GridSettings = settings

    def save(self, filename):
        np.save(filename, self)



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



file = "Sim"
sim = Simulation(dom, in_f, grid)
sim.save(file)



