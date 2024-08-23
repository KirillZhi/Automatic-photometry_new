import pandas as pd
from frame_functions import *
from coord_functions import *
from astropy.io import fits
from Tkinter import Tk
from tkFileDialog import askdirectory

import os
import re
from matplotlib import pyplot as plt

# Filepath for fits frames

readout_noise = 9
root = Tk()
root.withdraw()
filepath = askdirectory(initialdir="D:/Photometry")
photometry_path = askdirectory(initialdir=filepath)
filepath = filepath + "/"
photometry_path = photometry_path + "/"
root.update()
root.destroy()

print(filepath)
print(photometry_path)

# Search of fits frames in a folder
photo_frame_path = np.array([])
res = []
for path in os.listdir(filepath):
    if os.path.isfile(os.path.join(filepath, path)):
        res.append(path)
for i in res:
    # print(i)
    x = re.search(".*fit$", i)
    y = re.search(".*fits$", i)
    if x or y:
        photo_frame_path = np.append(photo_frame_path, i)
print(photo_frame_path)

# Reading of a test fits frame
filename = photo_frame_path[0]
check_path = filepath + filename
check_frame = fits.getdata(check_path)

check_shape = check_frame.shape
check_mid_x = check_shape[0] / 2
check_mid_y = check_shape[1] / 2
# Minimum value of the frame (positive)
test_ground = np.where(check_frame > 0, check_frame, np.max(check_frame))
c_min = np.percentile(test_ground, 5)
print(c_min)
#Looking for sources
c_min = 2300

centered_frame = np.copy(check_frame)
centered_frame = np.rot90(centered_frame, 1)
#circle_check(check_x, check_y, 5, 3, 12, centered_frame)
plt.figure(2)

im = plt.imshow(centered_frame - c_min)
plt.colorbar(im)
#plt.clim([0, math.sqrt(c_min + readout_noise**2) * 5])
plt.clim([140,300])
plt.show()
while (1):
    check_for_adding = int(
        raw_input("Do you want to stop adding checks? 1 to stop adding, 0 to keep adding: ", ))
    if check_for_adding == 1 or check_for_adding == 0:
        break
    else:
        check_for_adding = raw_input(
            "Unregistered command, repeat. Do you want to stop adding checks? 1 to stop adding, 0 to keep adding: ", )
plt.close(2)
