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
filepath = askdirectory(initialdir="C:/Users/kirio/Desktop/astronomy/GRB/")
photometry_path = askdirectory(initialdir="C:/Users/kirio/Desktop/astronomy/GRB/")
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

# Search of reference and check stars
plt.ion()
plt.figure(1)
im = plt.imshow(check_frame)
plt.colorbar(im)
# c_min=20
plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
#plt.clim([9500, 11000])
# plt.xlim([155, 355])
# plt.ylim([155, 355])
plt.show()
# plt.close(fig)

# Inputting coordinates and magnitude of an object
while (1):
    ell_circle_check = 1
    while (1):
        if ell_circle_check == 1:
            annul_check = int(input("An ellipse annulus or a circle one? [1 for ellipse, 0 for circle]: "))
            ell_circle_check = 0
        else:
            annul_check = int(
                input("An unregistered command. An ellipse annulus or a circle one? [1 for ellipse, 0 for circle]: "))
        if annul_check == 1 or annul_check == 0:
            break
    source_x = int(raw_input("Y coordinate of an object: ", ))
    source_y = int(raw_input("X coordinate of an object: ", ))
    source_color = float(raw_input("GAIA BP-RP color of an object (even an estimation): ", ))
    ap_radius = int(raw_input("Radius of an aperture (for source, reference, and checks): ", ))
    source_gap = int(raw_input("Length of a gap between aperture and background annulus (source): ", ))
    if annul_check == 0:
        source_annul = int(raw_input("Outer radius of a background annulus (source): ", ))
        test_frame = np.copy(check_frame)
        circle_ref(source_x, source_y, ap_radius, source_gap, source_annul,
                   test_frame)  # Drawing aperture and annulus on a frame
    else:
        source_annul1 = int(raw_input("1st outer radius of a background annulus (source): ", ))
        source_annul2 = int(raw_input("2nd outer radius of a background annulus (source): ", ))
        source_angle = int(raw_input("Angle of a background annulus ellipse(source): ", ))
        test_frame = np.copy(check_frame)
        ellipse_source(source_x, source_y, ap_radius, source_gap, source_annul1, source_annul2, source_angle,
                       test_frame)  # Drawing aperture and annulus on a frame
    plt.figure(2)
    im = plt.imshow(test_frame)
    plt.colorbar(im)
    plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
    plt.show()
    cont_check = 1
    while (1):
        if cont_check == 1:
            cont_check = int(raw_input("Do you want to change anything? [1 for changing, 0 to not change]: ", ))
        else:
            cont_check = raw_input(
                "Unregistered command, repeat. Do you want to change anything? [1 for changing, 0 to not change]: ", )
        if cont_check == 1 or cont_check == 0:
            break
    if cont_check == 0:
        plt.close(2)
        break
    else:
        plt.close(2)

final_frame = np.copy(test_frame)

#Inputting coordinates of checks
check_x = []
check_y = []
check_ra = []
check_dec = []
check_true_mag = []
check_color = []
check_for_gaps = 1
secondary_frame = np.copy(final_frame)
check_for_adding = 0
while (1):
    test_x = input("Y coordinate of a check {}: ".format(np.shape(check_x)[0] + 1), )  # X coordinate for a check
    test_y = input("X coordinate of a check {}: ".format(np.shape(check_x)[0] + 1), )  # Y coordinate for a check
    test_ra = input("RA coordinate of a check {}: ".format(np.shape(check_x)[0] + 1), )  # RA coordinate for a check
    test_dec = input("DEC coordinate of a check {}: ".format(np.shape(check_x)[0] + 1), )  # DEC coordinate for a check
    test_mag = input(
        "GAIA G magnitude of a check {}: ".format(np.shape(check_x)[0] + 1), )  # GAIA G magnitude of a star
    test_color = input(
        "GAIA BP-RP color of a check {}: ".format(np.shape(check_x)[0] + 1), )  # GAIA BP-RP color of a star
    if check_for_gaps == 1:
        check_gap = int(
            raw_input("Length of a gap between aperture and background annulus (common for all checks): ", ))
        check_annul = int(raw_input("Outer radius of a background annulus (common for all checks): ", ))
    check_for_gaps += 1

    test_frame = np.copy(secondary_frame)
    circle_ref(test_x, test_y, ap_radius, check_gap, check_annul, test_frame)  # Drawing aperture and annulus on a frame
    plt.figure(2)
    im = plt.imshow(test_frame)
    plt.colorbar(im)
    plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
    plt.show()
    while (1):
        cont_check = int(raw_input("Do you want to change anything? [1 for changing, 0 to not change]: ", ))
        if cont_check == 1 or cont_check == 0:
            break
        else:
            cont_check = raw_input(
                "Unregistered command, repeat. Do you want to change anything? [1 for changing, 0 to not change]: ", )
    if cont_check == 0:
        plt.close(2)
        check_x = np.append(check_x, test_x)
        check_y = np.append(check_y, test_y)
        check_ra = np.append(check_ra, test_ra)
        check_dec = np.append(check_dec, test_dec)
        check_true_mag = np.append(check_true_mag, test_mag)
        check_color = np.append(check_color, test_color)
        secondary_frame = np.copy(test_frame)
        while (1):
            check_for_adding = int(
                raw_input("Do you want to stop adding checks? 1 to stop adding, 0 to keep adding: ", ))
            if check_for_adding == 1 or check_for_adding == 0:
                break
            else:
                check_for_adding = raw_input(
                    "Unregistered command, repeat. Do you want to stop adding checks? 1 to stop adding, 0 to keep adding: ", )
    else:
        plt.close(2)
        check_for_gaps -= 1
    if check_for_adding == 1:
        break
print("Source and checks are chosen")
check_x_centre = np.copy(check_x)
check_y_centre = np.copy(check_y)
centered_frame = fits.getdata(check_path)
for i in range(check_x.shape[0]):
    #print(check_x[i])
    #print(check_y[i])
    check_x_centre[i], check_y_centre[i] = obj_centroid(check_x[i], check_y[i], ap_radius, check_gap, check_annul, centered_frame)
print(check_x_centre, check_y_centre)

#Using checks for making a coordinate net
tester_ra, tester_dec = 12.52242886860, -56.06811904258

x_vector_param, y_vector_param = coord_net_fit(check_x_centre, check_y_centre, check_path, check_ra, check_dec)
print(x_vector_param, y_vector_param)
#Finding x, y of a tester
tester_x, tester_y = coord_net_use(x_vector_param, y_vector_param, check_path, tester_ra, tester_dec)

circle_ref(tester_x, tester_y, ap_radius, check_gap, check_annul, centered_frame)
plt.figure(2)
im = plt.imshow(centered_frame)
plt.colorbar(im)
plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
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
