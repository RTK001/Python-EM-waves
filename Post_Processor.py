

'''
A program to Solve the Electromagnetic Wave Equation in a closed domain.
'''

import numpy as np
import matplotlib as mp
from matplotlib.pyplot import plot
from matplotlib import animation

from Wave_Eqn_Class_Defs import *





Test = PostProcessor()
Test.load_sim("Input.npy", "Res.npy")


Test.element_gain( (10,10))

Test.E_field_time_plot((10,10))

Test.E_field_contour_at_t(50)

Test.E_field_surface_at_t(50)

ani = Test.animate_t( Test.E_field_contour_at_t)

plt.show()
