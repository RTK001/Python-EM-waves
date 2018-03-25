
'''
A program to Solve the Electromagnetic Wave Equation in a closed domain.
'''


import numpy as np
import matplotlib.pyplot as plt


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

        # Calculate number of dimensions from Domain_length
        self.no_of_dimensions = np.size(self.Domain.Domain_length)

        #Define Physical Constant
        self.speed_of_light = 3*10**8  # Meters / Second

        # Calculate and output solver spatial resolution
        # Calcualtes the wavelength of the input frequency, sizing the Element to match
        self.Element_Length = np.ones(self.no_of_dimensions) * self.speed_of_light / [(self.InputFrequency.frequency * self.GridSettings.Elements_per_Wavelength )]
                                
        # Rounds the number of elements to an integer to approximate the domain size
        self.no_of_elements = ( self.Domain.Domain_length / self.Element_Length ).astype(int)

        # Calculate and output solver time resolution
        self.timestep_size = 1/(self.InputFrequency.frequency * self.GridSettings.TimeSteps_per_Wavelength )              # Sizes the timesteps as desired
        self.no_of_timesteps = int(self.Domain.Simulation_Duration / self.timestep_size)            # rounds the number of timesteps
        
        # Calculate and output solver stability checks
        self.alpha = (self.speed_of_light * self.timestep_size) **2 / sum(self.Element_Length)**2
        # Alpha is the coefficient linking the previous timestep with the current during solve.
        # Shouldn't be higher that 0.5 [Von-Neumann]

        # create an empty (0,0,..) vector
        self.cartesian_vector = np.zeros(self.no_of_dimensions)

        # Initial Values
        self.Initial_E = self.cartesian_vector        # used to initialise fields in elements
        self.Initial_B = self.cartesian_vector


        
                
    def save(self, filename):
        ''' Saves the input file in the given filename in the current directory '''
        np.save(filename, self)

    def print_solver_params(self):
        ''' prints key solver parameters to the screen '''
        
        print("Element Length: {}".format(self.Element_Length))
        print("Number of Elements: {}".format(self.no_of_elements))
        print("CFL no.: {}".format(sum(self.speed_of_light * self.timestep_size / self.Element_Length)))
        print("alpha:{}".format(self.alpha))

    def plot_input_function(self):
        ''' plots the input function as a matplotlib plot '''
        
        # Creates a sine wave of the input function
        input_function = np.array([np.sin(2*np.pi* self.InputFrequency.frequency * self.timestep_size*t)
                                   for t in range(self.InputFrequency.input_duration)], np.double)
        
        input_function.resize(self.no_of_timesteps)
        for dim in range(self.no_of_dimensions):         
            plt.plot(np.arange(sim.InputFrequency.input_duration) * self.timestep_size, 
                 input_function[:sim.InputFrequency.input_duration])


            

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



