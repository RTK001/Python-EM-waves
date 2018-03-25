import numpy as np
import matplotlib.pyplot as plt


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





class BaseElement():
    '''
    Base class for all Elemenets to inherit from
    '''
    Last_ID = 0             # ID to identify each element
    def __init__(self, sim_obj, coords, length):
        
        # Element IDs
        BaseElement.Last_ID += 1        
        self.ID = BaseElement.Last_ID
        
        # Coordinates
        self.coordinates = coords
        # Length is typically Element_Legnth variable, stored in sim, though this equation formulation can support variable element lengths.
        self.length = length
        # Alpha is the coefficient used to relate the spatial terms of the wave equation to the time terms.
        self.alpha = (sim_obj.speed_of_light * sim_obj.timestep_size) **2 / (sim_obj.no_of_dimensions * self.length) **2
        
        # Electromagnetic Properties
        self.E = np.array([sim_obj.Initial_E for i in range(sim_obj.no_of_timesteps)])  # empty solution array, size: (cartesian vector by timesteps)
        self.B = np.array([sim_obj.Initial_B for i in range(sim_obj.no_of_timesteps)])
        
        self.positive_neighbours = np.array([None for i in range(sim_obj.no_of_dimensions)], dtype = object)    # empty array to hold neighbouring elements (+ve direction along each axis)
        self.negative_neighbours = self.positive_neighbours.copy()              # Similar for -ve direction in each axis.
        
        
    def link(self, elem, direction):
        '''
        Links two elements for calculations.
        Directions are to be specified as: [x+, y+, z+]
        '''
        self.positive_neighbours[direction] = elem  
        elem.negative_neighbours[direction] = self
    



class BoundaryElement(BaseElement):
    '''
    Element with pre-determined values through time
    '''
    def __init__(self, sim_obj, coords, length):
        BaseElement.__init__(self, sim_obj, coords, length)

    def SetE(self, direction, E_Array):
        ''' Takes an array of Eelectric field in a given direction and sets the element value to that throughout the simualtion.
        '''
        self.E[:,direction] = E_Array
    
    def SetB(self, direction, B_Array):
        ''' Takes an array of MAgnetic field in a given direction and sets the element value to that throughout the simualtion.
        '''
        self.B[:][direction] = B_Array 
        
    def Calculate(self, timestep_to_calculate):
        ''' Does nothing for a boundary element - values should be already calculated.
        '''
        pass    



class Element(BaseElement):
    '''
    Element to be solved in the simulation
    '''
    def __init__(self, sim_obj, coords, length):
        BaseElement.__init__(self, sim_obj, coords, length)
        self.LastSolvedTimestep = 1
        
        
        
    def Calculate(self, timestep_to_calculate):
        
        # re-define for shorter syntax
        t = timestep_to_calculate
        prev_t =  t-1

        # initialises gradient matrices
        dEdt2 = np.zeros(np.shape(self.E[0]))
        dBdt2 = np.zeros(np.shape(self.B[0]))
        
        for dim in range(len(self.E[0])):

            # Calculate Spatial Derivatives
            # These these are evaluated at the element boundary - ie/ the mid-point between the two elements.
            forward_derivatives_E = (self.positive_neighbours[dim].E[prev_t] - self.E[prev_t])
            backward_derivatives_E = ( - self.negative_neighbours[dim].E[prev_t] + self.E[prev_t])
            forward_derivatives_B = (self.positive_neighbours[dim].B[prev_t] - self.B[prev_t])
            backward_derivatives_B = ( - self.negative_neighbours[dim].B[prev_t] + self.B[prev_t])
            
            # Calculate Second time gradient
            # Evaluated across the current element - ie/ between the two boundaries.
            dEdt2 += self.alpha * (forward_derivatives_E - backward_derivatives_E)
            dBdt2 += self.alpha * (forward_derivatives_B - backward_derivatives_B)
            
               
        # Calculate New E and B values for the timestep
        # Derived from Taylor expansion
        self.E[t] = dEdt2 + 2*self.E[prev_t] - self.E[t-2]
        self.B[t] = dBdt2 + 2*self.B[prev_t] - self.B[t-2]

        self.LastSolvedTImestep = t
                
       


def Element_Array_builder(sim_obj):
    '''
    Builds an Element array of elements ready to calculate.
    Creates and populates an array of coordiantes.
    Creates a mask array to select boundary elements at the edges, using the highest and lowest domain values.
        (only valid for linear, square, cube, etc. domains.)
    Creates an array of elements, then set relevant elements as boundary elements. 
    '''
    # Create Coordinate Array as an array of X-coordiante, Y-Coordinate and Z-Coordinate grids of element x element
    Coordinate_Array = np.array([np.ones(sim_obj.no_of_elements) for dim in range(sim_obj.no_of_dimensions)])

    # Populate the Coordinate matrix, assuming constant Element Length
    for dim in range (sim_obj.no_of_dimensions):
        axis = Coordinate_Array[dim]
        for indices, point in np.ndenumerate(axis):
                axis[indices] = sim_obj.Element_Length[dim] * indices[dim]

    # Assume maximum and minimum values - used to determine boundary elements
    Min_vals = np.zeros(sim_obj.no_of_dimensions)
    Max_vals = sim_obj.Element_Length * (sim_obj.no_of_elements-1)
    
    
    # Create a Boundary Element mask array
    Boundary_Array = np.zeros(sim_obj.no_of_elements, dtype = bool)
    for dim in range(sim_obj.no_of_dimensions):
        Boundary_Array += (Coordinate_Array[dim] == Min_vals[dim]) + (Coordinate_Array[dim] == Max_vals[dim])

    # Initialise the Element Array
    Element_Array = np.empty(sim_obj.no_of_elements , dtype = object)
    
    # Populate the Element Array with either Boundary Elements or regular Elements, using the Boundary_Array mask
    for indices, coords in np.ndenumerate(Boundary_Array):
        
        # Gather the coordinates for the Element
        coordinates = np.zeros(sim_obj.no_of_dimensions)
        for dim in range(sim_obj.no_of_dimensions):
            coordinates[dim] = Coordinate_Array[dim][indices]
        if Boundary_Array[indices]:
            Element_Array[indices] = BoundaryElement( sim_obj, coordinates, sim_obj.Element_Length)
        else:
            Element_Array[indices] = Element( sim_obj, coordinates, sim_obj.Element_Length)


    
    # Link elements
    # Iterate through all elements and link each one with the one ahead.
    # The link fucntion is symmetric, so the element ahead's negative neighbour will be set to the one behind.
    for index, elem in np.ndenumerate(Element_Array):
        for dim in range(sim_obj.no_of_dimensions):
            if index[dim] == sim_obj.no_of_elements[dim]-1:
                continue
            neighbour_index = list(index)   #Get the index of the current element
            neighbour_index[dim] += 1       # increment it to get the next element in current axis
            neighbour_index = tuple(neighbour_index)    # make it a tuple to use it as an index
            elem.link(Element_Array[neighbour_index], dim)      # Link the current element with the one at the neighbour index.

        
    return Element_Array


