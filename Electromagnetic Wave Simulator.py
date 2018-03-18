
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



# Domain parameters
Domain_length = np.array([0.5, 0.5])
# An n dimensional array of the domain lengths in Meters.
# Later used to determine number of dimensions.
Simulation_Duration = 1* 10** -8   		# Seconds


# Set Physical Constant
speed_of_light = 3*10**8        # Meters / Second


# Input signal parameters
frequency = 2*10**9     # in Hz
input_duration = 13     # Measured in timesteps


# Calculate number of dimensions from Domain_length
no_of_dimensions = np.size(Domain_length)


# Mesh and timestep sizing
Elements_per_Wavelength = 14        # Based on input wavelength
TimeSteps_per_Wavelength = 25

# Calculate and output solver spatial resolution
Element_Length = np.ones(no_of_dimensions) * speed_of_light / [(frequency * Elements_per_Wavelength)]       # Calcualtes the wavelength of the input frequency, sizing the Element to match
no_of_elements = (Domain_length/Element_Length).astype(int)             # Rounds the number of elements to an integer to approximate the domain size
print("Element Length: {}".format(Element_Length))
print("Number of Elements: {}".format(no_of_elements))

# Calculate and output solver time resolution
timestep_size = 1/(frequency*TimeSteps_per_Wavelength)              # Sizes the timesteps as desired
no_of_timesteps = int(Simulation_Duration/timestep_size)            # rounds the number of timesteps

# Calculate and output solver stability checks
print("CFL no.: {}".format(sum(speed_of_light * timestep_size / Element_Length)))   # CFL number is useful a stability metric for parabolic diff. eqns.
alpha = (speed_of_light * timestep_size / Element_Length[0])**2
# Alpha is the coefficient linking the previous timestep with the current during solve.
# Shouldn't be higher that 0.5 [Von-Neumann]
print("alpha:{}".format(alpha))
print("alpha:{} dB".format(10*np.log10((speed_of_light * timestep_size / Element_Length[0])**2)))


# Create input function
input_function = np.array([np.sin(2*np.pi*frequency * timestep_size*t) for t in range(input_duration)], np.double)          # Creates a sine wave of the input function
input_function.resize(no_of_timesteps)                  # Resize the array to match the simulation size, adding zeros where appropriate to indicate no signal.

# Plot created input function
for dim in range(no_of_dimensions):         
    plot(range(input_duration),input_function[:input_duration])
    


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



# Calcualte Gain

# This is to check stability of calculation.
# Ideally, should approximate peaks and troughs at each wave in time
gain = lambda t,x,y: (A[x,y].E[t][0] / A[x,y].E[t-1][0])
print_gain = lambda t,x,y:print("{} \n".format(gain(t,x,y)) )

rng = [10,no_of_timesteps-10]
x = 2
y = 1

res = [gain(i,x,y) for i in range(rng[0], rng[1])]

mp.pyplot.figure(figsize = (20,10))
mp.pyplot.axis([0, rng[1], 0, res[0] ])
mp.pyplot.plot(range(rng[0], rng[1]), res)

test1 = np.array([i+1 for i in range(rng[0],rng[1])])
mp.pyplot.plot(test1, np.e**(-test1*alpha) +3.5)

print_gain(rng[0],x,y)
print_gain(rng[1],x,y)

print_gain(30,x,y)



# Plot the E field in time for each element for 1D Domain
if no_of_dimensions == 1:
    t_steps = 5         # number of timesteps between plotted curves
    for t in range(0,int(no_of_timesteps/40), t_steps):
        mp.pyplot.figure(figsize = (20,10))
        mp.pyplot.axis([0,xelements, -2,2])
        Es = [elem.E[t] for elem in A]
        plot(range(0,len(A)), Es)


# Plot the E Field for rows of elements at column intervals
if no_of_dimensions > 1:

    xelements = no_of_elements[0]       # number of elements in row to plot
    yelements = no_of_elements[1]       # number of rows in column to plot
    plots = 5                           # number of plots
    timestep_min = 0                
    timestep_max = no_of_timesteps

    for plt in range(plots):

        y_axis_to_plot = int(plt*yelements/plots)
        mp.pyplot.figure(figsize = (20,10))
        mp.pyplot.axis([0,xelements, -2,2])

        

        for indices, element in np.ndenumerate( A[0:xelements, y_axis_to_plot] ):
            mp.pyplot.plot(range(timestep_min,timestep_max ), element.E[timestep_min:timestep_max], label = indices)
            mp.pyplot.title("Y_Elements = {0}".format(y_axis_to_plot))
        mp.pyplot.legend()


# Calculate magnitude of E_field

E_field = np.array([np.zeros(np.shape(A)) for t in range(no_of_timesteps)])     # Initialise array

for t in range(no_of_timesteps):    
    for index, elem in np.ndenumerate(E_field[t,:,:]):
        E_field[t][index] = sum(np.round_((A[index].E[t]),decimals = 10))       # Rounded to avoid overflow/underflow errors when summing E-field components

E_field.astype(np.float32)          # Reduce the size of the matrix in memory, to speed up plotting and animation
        
print(E_field[20])      # Print value, for checks




def time_varying_contour(x):
    ''' Plot a contour at the given timestep. '''
    mp.pyplot.figure(figsize = (20,10))
    the_contourf = mp.pyplot.contour(E_field[x])
    mp.pyplot.colorbar()
    




from mpl_toolkits.mplot3d import Axes3D

coordinates_x = np.array([elem.coordinates[0] for index,elem in np.ndenumerate(A)]).reshape(no_of_elements)
coordinates_y = np.array([elem.coordinates[1] for index, elem in np.ndenumerate(A)]).reshape(no_of_elements)

    
def time_varying_3d_contour(t):
    ''' Creates a countour plot at time t, of the E_field, and coordiantes grid. '''
    fig = mp.pyplot.figure(figsize = (20,10))
    ax = mp.pyplot.gcf().add_subplot(1,1,1, projection = '3d')
    ax.view_init(30,-35)
    the_contourf = ax.plot_surface( coordinates_x, coordinates_y, E_field[t])
    

    

def manual(t = 0, elevation = 30, azimuth = 30):
    ''' Plots the contour at a specific elevationa nd azimuth '''
    fig = mp.pyplot.figure(figsize = (20,10))
    ax = fig.add_subplot(1,1,1,projection = '3d')
    ax.view_init(elevation,azimuth)
    the_contourf = ax.plot_surface( coordinates_x, coordinates_y, E_field[t])
    


# Create Animated E-field plot.

fig = mp.pyplot.figure(figsize = (10,5))
ax = mp.pyplot.gcf().add_subplot(1,1,1,projection = '3d')
ax.view_init(60,-60)
ax.set_zlim(-0.125, 0.125)
surface = ax.plot_surface( coordinates_x, coordinates_y, E_field[0], cmap = mp.cm.plasma)


def animate_t(t):
    ax.cla()
    surface = ax.plot_surface( coordinates_x, coordinates_y, E_field[t], cmap = mp.cm.plasma)
    ax.set_zlim(-0.25, 0.25)
    
        
ani = mp.animation.FuncAnimation(fig, animate_t, np.arange(10,no_of_timesteps,1),  interval = 1)

