from numpy import zeros, empty, append, arange, max, absolute, array, delete, count_nonzero, argmax, reshape, unique
from numpy import min as arr_min
from numpy import max as arr_max

from pandas import DataFrame
from itertools import combinations

from math import sqrt
from frame_functions import *
#from coord_functions import *
from not_finished_coord_functions import *

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

print("Choose directory with frames")
print('')
filepath = askdirectory(initialdir="D:/Photometry/") + '/'
photo_frame_path = np.array([])
res = []
for path in listdir(filepath):
    if isfile(join(filepath, path)):
        res.append(path)
for i in res:
    x = search(".*fit$", i)
    y = search(".*fits$", i)
    if x or y:
        photo_frame_path = np.append(photo_frame_path, i)
sum_frame = np.zeros_like(getdata("{}".format(filepath + photo_frame_path[0])))
for frame_num in range(len(photo_frame_path)):
    check_frame = getdata("{}".format(filepath + photo_frame_path[frame_num]))
    sum_frame += check_frame
sum_frame = np.rot90(sum_frame)
im = plt.imshow(sum_frame)
plt.colorbar(im)
mean = astropy.stats.sigma_clipped_stats(sum_frame, sigma=1.5)[0]
c_min = mean
plt.clim([c_min, c_min + 8 * sqrt(c_min)])
plt.gca().invert_yaxis()
#plt.gca().invert_xaxis()
plt.show()
plt.close()