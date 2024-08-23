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

def axis_fit_distortion(input_x, input_y, supposed_x_0, supposed_y_0, input_parameters):
    # input_x, input_y - distorted coordinates of anchors on the frame
    # supposed_x, supposed_y - undistorted (supposedly) coordinates of anchors on the frame
    # input_parameters - rotation, shift and scale parameters
    # Initial center
    frame_centre_init_x = input_parameters[0] * cos(input_parameters[2]) + input_parameters[1] * sin(
        input_parameters[2])
    frame_centre_init_y = - input_parameters[0] * sin(input_parameters[2]) + input_parameters[1] * cos(
        input_parameters[2])

    supposed_x = supposed_x_0 * cos(input_parameters[11]) + supposed_y_0 * sin(input_parameters[11]) + input_parameters[
        9]
    supposed_y = -supposed_x_0 * sin(input_parameters[11]) + supposed_y_0 * cos(input_parameters[11]) + \
                 input_parameters[10]
    deviation_x = supposed_x - input_x  # Deviation between distorted x and non-distorted x
    deviation_y = supposed_y - input_y  # Deviation between distorted y and non-distorted y
    frame_centre_x = frame_centre_init_x
    frame_centre_y = frame_centre_init_y
    centered_x = input_x - frame_centre_x  # x as measured from the center of the full frame
    centered_y = input_y - frame_centre_y  # y as measured from the center of the full frame

    centered_radius = arr_sqrt(centered_x ** 2 + centered_y ** 2)  # Distance from the center of the frame (distorted)
    # Additional variables
    a_k = centered_radius ** 2 + 2 * centered_x ** 2
    b_k = 2 * centered_x * centered_y
    c_k = centered_radius ** 2 + 2 * centered_y ** 2

    radius_six_sum = sum(centered_radius ** 6)
    radius_eight_sum = sum(centered_radius ** 8)
    radius_ten_sum = sum(centered_radius ** 10)

    # Sums used in the calculation (refer to the calculations or Kirill for the explanations)
    sum_one = sum(centered_radius ** 2 * centered_x * a_k + centered_radius ** 2 * centered_y * b_k)
    sum_two = sum(centered_radius ** 2 * centered_x * b_k + centered_radius ** 2 * centered_y * c_k)
    sum_three = sum(centered_radius ** 2 * centered_x * deviation_x + centered_radius ** 2 * centered_y * deviation_y)
    sum_four = sum(centered_radius ** 4 * centered_x * a_k + centered_radius ** 4 * centered_y * b_k)
    sum_five = sum(centered_radius ** 4 * centered_x * b_k + centered_radius ** 4 * centered_y * c_k)
    sum_six = sum(centered_radius ** 4 * centered_x * deviation_x + centered_radius ** 4 * centered_y * deviation_y)
    sum_seven = sum(a_k ** 2 + b_k ** 2)
    sum_eight = sum(deviation_x * a_k + deviation_y * b_k)
    sum_nine = sum(b_k ** 2 + c_k ** 2)
    sum_ten = sum(deviation_x * b_k + deviation_y * c_k)
    sum_eleven = sum(4 * centered_radius ** 2 * b_k)

    # Setting up arrays for the system of equations
    left_array = array(
        [[radius_six_sum, radius_eight_sum, sum_one, sum_two], [radius_eight_sum, radius_ten_sum, sum_four, sum_five],
         [sum_one, sum_four, sum_seven, sum_eleven], [sum_two, sum_five, sum_eleven, sum_nine]])
    right_array = array([[sum_three], [sum_six], [sum_eight], [sum_ten]])
    gauss_left_array, gauss_right_array = gauss_elimination(left_array, right_array)

    # Estimating distortion parameters from the gaussed arrays
    # print('All 4 parameters')
    for parameter_ind in range(gauss_left_array.shape[0] - 1, -1, -1):
        substract_value = 0
        if parameter_ind != gauss_left_array.shape[0] - 1:
            for substract_ind in range(parameter_ind + 1, gauss_left_array.shape[0]):
                substract_value += input_parameters[5 + substract_ind] * gauss_left_array[parameter_ind][substract_ind]
        right_value = gauss_right_array[parameter_ind] - substract_value
        input_parameters[5 + parameter_ind] = right_value / gauss_left_array[parameter_ind][parameter_ind]
    # Estimating sum of squared deviations after correction
    min_value = 0
    for min_ind in range(len(input_x)):
        f_fun = deviation_x[min_ind] - centered_x[min_ind] * (
                input_parameters[5] * centered_radius[min_ind] ** 2 + input_parameters[6] * centered_radius[
            min_ind] ** 4)

        f_fun -= input_parameters[7] * a_k[min_ind] + input_parameters[8] * b_k[min_ind]
        g_fun = deviation_y[min_ind] - centered_y[min_ind] * (
                input_parameters[5] * centered_radius[min_ind] ** 2 + input_parameters[6] * centered_radius[
            min_ind] ** 4)
        g_fun -= input_parameters[7] * b_k[min_ind] + input_parameters[8] * c_k[min_ind]
        min_value += f_fun ** 2 + g_fun ** 2

    return min_value, input_parameters
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
im = plt.imshow(check_frame - c_min)
plt.colorbar(im)
# c_min=20
plt.clim([0, math.sqrt(c_min + readout_noise**2) * 8])
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


print("Source and checks are chosen")

centered_frame = fits.getdata(check_path)
source_x_centre, source_y_centre = obj_centroid(source_x, source_y, ap_radius, source_gap, source_annul, centered_frame)


# Photometry
photo_len = len(photo_frame_path)
ref_photo = np.zeros(photo_len)
frame_datetime = np.zeros_like(ref_photo)
frame_exp = np.zeros_like(ref_photo)
frame_tube = np.array([])
frame_filter = np.array([])

source_est_mag = np.zeros(photo_len)
source_photo = np.zeros(photo_len)
source_err = np.zeros(photo_len)
noise_lim_intensity = np.zeros(photo_len)
noise_lim_mag = np.zeros(photo_len)
noise_lim_pixel = np.zeros(photo_len)
# Conducting the photometry
for frame_num in range(len(photo_frame_path)):
    test_frame = fits.open("{}".format(filepath + photo_frame_path[frame_num]))
    test_frame_header = test_frame[0].header

    photo_frame = fits.getdata("{}".format(filepath + photo_frame_path[frame_num]))
    noise_lim_mag[frame_num] = test_frame_header['MLIMIT']
    if annul_check == 0:
        source_photo[frame_num], source_err[frame_num], noise_lim_intensity[frame_num], noise_lim_pixel[frame_num] = photometry_source(source_x, source_y, ap_radius, source_gap,
                                                                           source_annul, photo_frame)
    else:
        source_photo[frame_num], source_err[frame_num] = photometry_source_ell(source_x, source_y, ap_radius,
                                                                               source_gap, source_annul1, source_annul2,
                                                                               source_angle, photo_frame)
    source_est_mag[frame_num] = magnitude_obj(noise_lim_mag[frame_num], noise_lim_intensity[frame_num], source_photo[frame_num])
print(source_est_mag)
print(noise_lim_intensity, noise_lim_pixel)
print(source_photo / source_err)
input()
