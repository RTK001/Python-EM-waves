
'''
A program to Solve the Electromagnetic Wave Equation in a closed domain.
'''


import numpy as np
import matplotlib as mp
from matplotlib.pyplot import plot
from matplotlib import animation

# Used for debugging of numerical errors
import pdb
np.seterr(all='raise')



# Extract input file parameters

# Load the sim file
input_file_name = "Sim.npy"
sim = np.load(input_file_name).item(0)



# Set Physical Constant
speed_of_light = 3*10**8        # Meters / Second


# Calculate number of dimensions from Domain_length
no_of_dimensions = np.size( sim.Domain.Domain_length )


# Calculate and output solver spatial resolution
# Calcualtes the wavelength of the input frequency, sizing the Element to match
Element_Length = np.ones(no_of_dimensions) * speed_of_light / [(sim.InputFrequency.frequency * sim.GridSettings.Elements_per_Wavelength )]       
no_of_elements = ( sim.Domain.Domain_length / Element_Length ).astype(int)             # Rounds the number of elements to an integer to approximate the domain size
print("Element Length: {}".format(Element_Length))
print("Number of Elements: {}".format(no_of_elements))

# Calculate and output solver time resolution
timestep_size = 1/(sim.InputFrequency.frequency * sim.GridSettings.TimeSteps_per_Wavelength )              # Sizes the timesteps as desired
no_of_timesteps = int(sim.Domain.Simulation_Duration / timestep_size)            # rounds the number of timesteps

# Calculate and output solver stability checks
print("CFL no.: {}".format(sum(speed_of_light * timestep_size / Element_Length)))   # CFL number is useful a stability metric for parabolic diff. eqns.
alpha = (speed_of_light * timestep_size / Element_Length[0])**2
# Alpha is the coefficient linking the previous timestep with the current during solve.
# Shouldn't be higher that 0.5 [Von-Neumann]
print("alpha:{}".format(alpha))
print("alpha:{} dB".format(10*np.log10((speed_of_light * timestep_size / Element_Length[0])**2)))


# Create input function
input_function = np.array([np.sin(2*np.pi* sim.InputFrequency.frequency * timestep_size*t) for t in range(sim.InputFrequency.input_duration)], np.double)          # Creates a sine wave of the input function
input_function.resize(no_of_timesteps)                  # Resize the array to match the simulation size, adding zeros where appropriate to indicate no signal.

# Plot created input function
for dim in range(no_of_dimensions):         
    plot(range(sim.InputFrequency.input_duration),input_function[:sim.InputFrequency.input_duration])
    


# pre-compute constant variables to reduce solve time
speed_of_light_squared_x_timestep_sq = (speed_of_light**2) * (timestep_size**2)
# pre-calculated alpha woithout the element length term, which is added later in the individual element. 
cartesian_vector = np.zeros(no_of_dimensions)               # create an empty (0,0,..) vector


# Initial Values
Initial_E = cartesian_vector        # used to initialise fields in elements
Initial_B = cartesian_vector




class BaseElement():
    '''
    Base class for all Elemenets to inherit from
    '''
    Last_ID = 0             # ID to identify each element
    def __init__(self, coords, length):
        
        # Element IDs
        BaseElement.Last_ID += 1        
        self.ID = BaseElement.Last_ID
        
        # Coordinates
        self.coordinates = coords
        self.length = length
        self.reciprocal_length = 1/self.length      # pre-calculated to speed up calculations
        
        # Electromagnetic Properties
        self.E = np.array([Initial_E for i in range(no_of_timesteps)])  # empty solution array, size: (cartesian vector by timesteps)
        self.B = np.array([Initial_B for i in range(no_of_timesteps)])
        
        self.positive_neighbours = np.array([None for i in range(no_of_dimensions)], dtype = object)    # empty array to hold neighbouring elements (+ve direction along each axis)
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
    def __init__(self, coords, length):
        BaseElement.__init__(self, coords, length)

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
    def __init__(self, coords, length):
        BaseElement.__init__(self, coords, length)
        self.LastSolvedTimestep = 1
        
        
        
    def Calculate(self, timestep_to_calculate):
        
        # re-define for shorter syntax
        t = timestep_to_calculate
        prev_t =  t-1

        # initialises gradient matrices
        dEdt2 = np.zeros(np.shape(Initial_E))
        dBdt2 = np.zeros(np.shape(Initial_B))
        
        for dim in range(no_of_dimensions):

            # Calculate Spatial Derivatives
            # These these are evaluated at the element boundary - ie/ the mid-point between the two elements.
            forward_derivatives_E = (self.positive_neighbours[dim].E[prev_t] - self.E[prev_t]) * self.reciprocal_length
            backward_derivatives_E = ( - self.negative_neighbours[dim].E[prev_t] + self.E[prev_t]) * self.reciprocal_length
            forward_derivatives_B = (self.positive_neighbours[dim].B[prev_t] - self.B[prev_t]) * self.reciprocal_length
            backward_derivatives_B = ( - self.negative_neighbours[dim].B[prev_t] + self.B[prev_t]) * self.reciprocal_length
            
            # Calculate Second time gradient
            # Evaluated across the current element - ie/ between the two boundaries.
            dEdt2 += speed_of_light_squared_x_timestep_sq * (forward_derivatives_E - backward_derivatives_E) * self.reciprocal_length
            dBdt2 += speed_of_light_squared_x_timestep_sq * (forward_derivatives_B - backward_derivatives_B) * self.reciprocal_length
            
               
        # Calculate New E and B values for the timestep
        # Derived from Taylor expansion
        self.E[t] = dEdt2 + 2*self.E[prev_t] - self.E[t-2]
        self.B[t] = dBdt2 + 2*self.B[prev_t] - self.B[t-2]
                
       


def Element_Array_builder():
    '''
    Builds an Element array of elements ready to calculate.
    Creates and populates an array of coordiantes.
    Creates a mask array to select boundary elements at the edges, using the highest and lowest domain values.
        (only valid for linear, square, cube, etc. domains.)
    Creates an array of elements, then set relevant elements as boundary elements. 
    '''
    # Create Coordinate Array as an array of X-coordiante, Y-Coordinate and Z-Coordinate grids of element x element
    Coordinate_Array = np.array([np.ones(no_of_elements) for dim in range(no_of_dimensions)])

    # Populate the Coordinate matrix, assuming constant Element Length
    for dim in range (no_of_dimensions):
        axis = Coordinate_Array[dim]
        for indices, point in np.ndenumerate(axis):
                axis[indices] = Element_Length[dim] * indices[dim]

    # Assume maximum and minimum values - used to determine boundary elements
    Min_vals = np.zeros(no_of_dimensions)
    Max_vals = Element_Length*(no_of_elements-1)
    
    
    # Create a Boundary Element mask array
    Boundary_Array = np.zeros(no_of_elements, dtype = bool)
    for dim in range(no_of_dimensions):
        Boundary_Array += (Coordinate_Array[dim] == Min_vals[dim]) + (Coordinate_Array[dim] == Max_vals[dim])

    # Initialise the Element Array
    Element_Array = np.empty(no_of_elements , dtype = object)
    
    # Populate the Element Array with either Boundary Elements or regular Elements, using the Boundary_Array mask
    for indices, coords in np.ndenumerate(Boundary_Array):
        
        # Gather the coordinates for the Element
        coordinates = np.zeros(no_of_dimensions)
        for dim in range(no_of_dimensions):
            coordinates[dim] = Coordinate_Array[dim][indices]
        if Boundary_Array[indices]:
            Element_Array[indices] = BoundaryElement(coordinates, Element_Length)
        else:
            Element_Array[indices] = Element(coordinates, Element_Length)


    
    # Link elements
    # Iterate through all elements and link each one with the one ahead.
    # The link fucntion is symmetric, so the element ahead's negative neighbour will be set to the one behind.
    for index, elem in np.ndenumerate(Element_Array):
        for dim in range(no_of_dimensions):
            if index[dim] == no_of_elements[dim]-1:
                continue
            neighbour_index = list(index)   #Get the index of the current element
            neighbour_index[dim] += 1       # increment it to get the next element in current axis
            neighbour_index = tuple(neighbour_index)    # make it a tuple to use it as an index
            elem.link(Element_Array[neighbour_index], dim)      # Link the current element with the one at the neighbour index.

        
    return Element_Array



# Test

# Create Array
A = Element_Array_builder()

# Setup Boundary conditions - input pulse
A[0, 1 ].SetE(0 ,input_function)


# Solve
for t in range(3,no_of_timesteps):
    for index, elem in np.ndenumerate(A):
        A[index].Calculate(t)


E_field = np.array([elem.E[:] for index, elem in np.ndenumerate(A)])
B_field = np.array([elem.B[:] for index, elem in np.ndenumerate(A)])

np.save("Res", np.array([E_field, B_field]))

