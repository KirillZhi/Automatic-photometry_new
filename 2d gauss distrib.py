import pandas as pd
from numpy import zeros, empty, append, arange, max, absolute, array, delete, count_nonzero, argmax, reshape, unique
from numpy import min as arr_min
from numpy import max as arr_max

from pandas import DataFrame
from itertools import combinations

from math import sqrt
from frame_functions import *
from coord_functions import *

from astropy.io.fits import getdata, open
from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename

from os import listdir, remove
from os.path import isfile, join

from re import search
from matplotlib import pyplot as plt
from time import time
import dateutil
import datetime
import openpyxl

def gauss_rad_distr(radius):
    #return radius * math.exp(-radius ** 2 / 2)
    return 2 / math.sqrt(2 * 3.14159) * math.exp(-radius ** 2 / 2)

# Making a table of values for a 2d gauss distribution in a radius form
r_values = np.linspace(0, 10, 2001)           # Values of equalized radius
distrib_density_values = np.zeros_like(r_values)    # Values of density distribution
distrib_function_values = np.zeros_like(r_values)
for i in range(1, len(r_values)):
    distrib_density_values[i] = gauss_rad_distr(r_values[i])
    distrib_function_values[i] = np.trapz(distrib_density_values[:i + 1], r_values[:i + 1])
data_for_df = {"Equalized radius": r_values, "Distribution function": distrib_function_values}
twod_dist_df = pd.DataFrame(data_for_df)
print('Choose where to save table of values')
filepath = askdirectory(initialdir="D:/Photometry/")
twod_dist_df.to_excel("{}".format(filepath) + "/one_d_gauss_radius_distribution.xlsx")
#plt.plot(r_values, distrib_function_values)
#plt.show()
#plt.close()


