import pandas as pd
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
import openpyxl
# Parameters of the initial fit (parameters of anchors)

num_of_default_anchors = 7  # Number for anchors of the initial net
offset_from_peak_flux = 0  # Number for the offset of initial anchors from the peak flux stars in both
# catalog and frame
# Parameters for a preferred magnitude deviation fit
color_param_num = 0
max_coord_poly_degr = 0         # 0 for both parameters if you just want to correct reference magnitude errors
max_param_num = (max_coord_poly_degr ** 2 + 3 * max_coord_poly_degr) / 2 + color_param_num + 1
# Parameters of the studied source
source_ra, source_dec = 197.645625, -21.718055555556 # Celestial (equatorial) coordinates of the studied source
#source_ra, source_dec = 159.40894640347, -28.37316991257
source_color = 0.31  # Color BP-RP of the studied source
supposed_magnitude = 0  # 0 if you don't know minimum possible magnitude of the studied source (peak brightness).
                        # If you do know (or guessed), it will not use low limit frames from the sample
source_mag_interval = np.array([12., 17.], dtype=float)
saturation_mag = 12     # Magnitude at which stars begin to be saturated (overflowing)
source_centroid_radius = 3  # Radius of the aperture that is used to calculate centroid (multiply by 3 to get correct radius)
bright_centre_radius = 3
centroid_radius = 6.0         # Radius for centroid aperture
source_aperture_radius = 0  # Change from 0 to any positive integer less than 7, it will serve as an aperture radius. If 0, adaptive aperture will be used
mask_check = 0              # If there is a relatively bright star in the vicinity (1-5 pixels) of the source,
                            # put 1. It will decrease its influence on the photometry
discard_level = 0.5         # Proportion of zero pixels in the frame upon which the frame is dropped.
                            # If number of zero pixels is more than total number of pixels * discard_level, the frame is dropped
                            # Generally, it should be 0.3. But the frame is important, leave at 0.5
                            # Parameters of catalog checks, ref and anchors
search_radius = 16.         # Radius in arcminutes where we look for checks and anchors in the catalog
check_mag_up = 11        # Upper limit on the flux of checks
mag_lower = 17.              # Lower limit on the flux of checks and anchors
anchor_mag_up = 11      # Upper limit on the flux of anchors

num_of_initial_anchors = copy(num_of_default_anchors)
# Arrays for parameter numbers
coord_num_arr = np.zeros(0)
color_num_arr = np.zeros(0)

mid_ra = source_ra
mid_dec = source_dec

# Filepath for fits frames
readout_noise = 9  # Readout noise of camera
root = Tk()
root.withdraw()
print("Choose directory with frames")
print('')
filepath = askdirectory(initialdir="D:/Photometry/")
print("Choose directory for exporting data")
photometry_path = askdirectory(initialdir=filepath)
print("Choose star catalog")
catalog_path = askopenfilename(initialdir=photometry_path)
filepath = filepath + "/"
photometry_path = photometry_path + "/"
root.update()
root.destroy()

# Determining limits of the frames
print('Pre-analysis of frames')
# print('Catalog analyzed, beginning operations with the frames')

# Search of fits frames in a folder
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

# Determining 3 sigma limits of frames
list_for_deletion = np.zeros(0)
limit_arr = zeros(0)
corrected_photo_frame_path = zeros(0)
for frame_num in range(len(photo_frame_path)):
    # Reading of a test fits frame

    filename = photo_frame_path[frame_num]
    # Parameters of a frame
    print(filename)
    check_frame = getdata("{}".format(filepath + photo_frame_path[frame_num]))
    if frame_num == 0:
        frame_shape = check_frame.shape
    zero_pixel_num = count_nonzero(check_frame == 0)
    if zero_pixel_num >= check_frame.shape[0] * check_frame.shape[1] * discard_level:
        if zero_pixel_num >= check_frame.shape[0] * check_frame.shape[1] * 0.5:
            list_for_deletion = append(list_for_deletion, filename)
        continue
    else:
        test_frame = open("{}".format(filepath + photo_frame_path[frame_num]))  # Opening the frame
        test_frame_header = test_frame[0].header  # Reading the header
        frame_limit = float(test_frame_header['MLIMIT'])
        print(frame_limit, supposed_magnitude)
        # Discarding frames with low limit
        if frame_limit + 0.8 >= supposed_magnitude:
            limit_arr = append(limit_arr, float(test_frame_header['MLIMIT']))
            corrected_photo_frame_path = append(corrected_photo_frame_path, filename)
        if frame_limit + 0.8 <= supposed_magnitude:
            list_for_deletion = append(list_for_deletion, filename)
        test_frame.close()

#check_frame = getdata("{}".format(filepath + corrected_photo_frame_path[0]))
check_frame = 0
# Deleting empty frames
print('Deleting ' + "{}".format(len(list_for_deletion)) + ' frames either for being empty or low limit')
for del_ind in range(len(list_for_deletion)):
    #print(list_for_deletion[del_ind])
    remove("{}".format(filepath + list_for_deletion[del_ind]))

if len(corrected_photo_frame_path) == 0:
    print('No frames to be analysed, quitting the program')
    quit()
# Deleting duplicate frames (with the same datetime of observation
pre_final_photo_frame_path = array([corrected_photo_frame_path[0]])
list_for_deletion = np.zeros(0)
duplicate_skip_check = int(raw_input('Do you want to skip checking for duplicates? Advised if you already checked: 0 to skip, 1 to check: '))
if duplicate_skip_check == 1:
    for frame_num in range(len(corrected_photo_frame_path)):
        discard_check = 0
        test_frame = open("{}".format(filepath + corrected_photo_frame_path[frame_num]))  # Opening the frame

        test_frame_header = test_frame[0].header  # Reading the header
        frame_inter_datetime = test_frame_header['DATE-OBS']
        for check_frame_num in range(len(pre_final_photo_frame_path)):
            check_frame = open("{}".format(filepath + pre_final_photo_frame_path[check_frame_num]))  # Opening the frame
            check_frame_header = check_frame[0].header  # Reading the header
            check_frame_inter_datetime = check_frame_header['DATE-OBS']
            if frame_inter_datetime == check_frame_inter_datetime and pre_final_photo_frame_path[check_frame_num] != corrected_photo_frame_path[frame_num]:
                print('discarding', frame_inter_datetime, check_frame_inter_datetime)
                discard_check = 1
                break
            check_frame.close()
        if discard_check == 0:
            pre_final_photo_frame_path = append(pre_final_photo_frame_path, corrected_photo_frame_path[frame_num])
        else:
            list_for_deletion = append(list_for_deletion, corrected_photo_frame_path[frame_num])
        test_frame.close()

    print('Deleting ' + "{}".format(len(list_for_deletion)) + ' frames for being duplicates')

if duplicate_skip_check == 0:
    pre_final_photo_frame_path = copy(corrected_photo_frame_path)

for del_ind in range(len(list_for_deletion)):
    print('delete')
    remove("{}".format(filepath + list_for_deletion[del_ind]))

photo_frame_path = copy(pre_final_photo_frame_path)
photo_frame_path = unique(photo_frame_path)

corrected_photo_frame_path = zeros(0)
print('Picked correct frames, analysing catalog')
print(arr_min(limit_arr))
mag_lower = arr_min(limit_arr)
#mag_lower = 15.2
print('Limit for checks: ' + "{}".format(mag_lower))
# Making quads from the catalogue
bright_stars, interm_bright_stars, bright_anchor_stars, catalog_stars, super_bright_stars, poss_check_stars = gaia_bright(
    catalog_path,
    num_of_initial_anchors,
    offset_from_peak_flux,
    search_radius,
    source_ra,
    source_dec, check_mag_up, anchor_mag_up, mag_lower, source_mag_interval, saturation_mag)

prev_num_of_bright_anchors = len(bright_anchor_stars['Gmag'])
print(prev_num_of_bright_anchors)

while len(interm_bright_stars['Gmag']) <= 20:
    print('aaa')
    mag_lower += 1
    bright_stars, interm_bright_stars, bright_anchor_stars, catalog_stars, super_bright_stars, poss_check_stars = gaia_bright(
        catalog_path,
        num_of_initial_anchors,
        offset_from_peak_flux,
        search_radius,
        source_ra,
        source_dec, check_mag_up, anchor_mag_up, mag_lower, source_mag_interval, saturation_mag)
    if prev_num_of_bright_anchors == len(interm_bright_stars['Gmag']) and len(interm_bright_stars['Gmag']) != 0:
        supposed_magnitude = mag_lower + 1.0
        # Again deleting frames with low limit
        list_for_deletion = np.zeros(0)
        limit_arr = zeros(0)
        corrected_photo_frame_path = zeros(0)
        for frame_num in range(len(photo_frame_path)):
            # Reading of a test fits frame
            filename = photo_frame_path[frame_num]
            # Parameters of a frame
            check_frame = getdata("{}".format(filepath + photo_frame_path[frame_num]))

            test_frame = open("{}".format(filepath + photo_frame_path[frame_num]))  # Opening the frame
            test_frame_header = test_frame[0].header  # Reading the header
            frame_limit = float(test_frame_header['MLIMIT'])
            # Discarding frames with low limit
            if frame_limit >= supposed_magnitude:
                limit_arr = append(limit_arr, float(test_frame_header['MLIMIT']))
                corrected_photo_frame_path = append(corrected_photo_frame_path, filename)
            if frame_limit <= supposed_magnitude:
                list_for_deletion = append(list_for_deletion, filename)
            del check_frame, test_frame_header
            test_frame.close()

        #for del_ind in range(len(list_for_deletion)):
        #    print('delete')
        #    remove("{}".format(filepath + list_for_deletion[del_ind]))

        photo_frame_path = copy(corrected_photo_frame_path)
        photo_frame_path = unique(photo_frame_path)
        corrected_photo_frame_path = zeros(0)
        break
    print(len(interm_bright_stars['Gmag']))

    prev_num_of_bright_anchors = len(interm_bright_stars['Gmag'])
    #print(len(bright_anchor_stars['Gmag']))

print(photo_frame_path)
print('Limit for checks: ' + "{}".format(mag_lower))
print(len(interm_bright_stars['Gmag']))
print("Number of anchors from the catalog")
print(len(bright_anchor_stars['Gmag']))



for x_ind in range(len(corrected_photo_frame_path)):
    for j_ind in range(len(list_for_deletion)):
        if corrected_photo_frame_path[x_ind] == list_for_deletion[j_ind]:
            print('aaaaa', list_for_deletion[j_ind])

# print(bright_stars)
bright_quad_comb, bright_stars_RA, bright_stars_DEC = quad_catalog(bright_stars)
print(bright_stars['Gmag'])
#print(bright_stars)
# Picking reference and checks

check_ref_stars, check_ref_color_array = check_ref_search(interm_bright_stars, source_ra, source_dec)
#print(check_stars['Gmag'])
num_of_check_refs = len(check_ref_stars['Gmag'])
if num_of_check_refs > 120:  # Removing an excess of check stars
    while num_of_check_refs > 120:
        color_max_num = argmax(check_ref_color_array)
        check_ref_color_array[color_max_num] -= 1
        num_of_check_refs -= 1
    min_color = arr_min(check_ref_stars['BP-RP'])
    max_color = arr_max(check_ref_stars['BP-RP'])
    color_number_array = zeros(int(round((max_color - min_color) / 0.1)) + 1)
    dict_stars = []
    dict_ref = []
    array_columns = check_ref_stars.columns.values.tolist()
    for check_ind in range(len(check_ref_stars['BP-RP'])):
        star_color = check_ref_stars['BP-RP'].values[check_ind]
        color_ind = int((star_color - min_color) / 0.1)
        if color_number_array[color_ind] < check_ref_color_array[color_ind]:
            dict1 = dict((col_name, check_ref_stars[col_name].values[check_ind]) for col_name in array_columns)
            color_number_array[color_ind] += 1
            dict_stars.append(dict1)
    check_ref_stars = pd.DataFrame(dict_stars, columns=array_columns)
#print("Number of checks")
#print(len(check_stars['Gmag']))
original_check_ref_stars = pd.DataFrame.copy(check_ref_stars)


print('Frames are picked, building the coordinate net')
print('Approximate time till completion: ' + "{}".format(len(photo_frame_path) * 20. / 60 / 2.4) + ' minutes')
start_time = time()
# Photometry parameters
ref_gap = 5
ref_annul = 20
check_gap = 5
check_annul = 20
photo_len = len(photo_frame_path)



check_ref_len = len(check_ref_stars['Gmag'])
check_ref_x = zeros(check_ref_len)
check_ref_y = zeros(check_ref_len)

every_frame_param = empty((0, 7), float)
initial_parameters = zeros(7)
initial_parameters[2] = 1000
initial_parameters[5] = 5

list_for_deletion = np.zeros(0)         # List for deleting frames where we can't make quads for some reason (likely fault with an hourclock machine)

dict_ref = []
frame_ind = 0
mag_dev_fit_param = np.zeros(0)
mag_dev_fit_param_err = np.zeros(0)
ref_photo = np.zeros(0)
source_photo_arr = np.zeros(0)
source_photo_err_arr = np.zeros(0)
source_x_arr = np.zeros(0)
source_y_arr = np.zeros(0)
photo_limit_arr = np.zeros(0)
mag_ref = np.zeros(0)
ap_radius_arr = np.zeros(0)
for frame_num in range(len(photo_frame_path)):
    num_of_initial_anchors = num_of_default_anchors
    discard_frame = 0
    # Reading of a test fits frame
    filename = photo_frame_path[frame_num]
    print('Analysing image ' + "{}".format(filename))
    check_path = filepath + filename
    check_frame = getdata(check_path)

    # Shape of a frame
    check_shape = check_frame.shape
    check_mid_x = check_shape[0] / 2
    check_mid_y = check_shape[1] / 2
    test_frame = open(check_path)
    frame_header = test_frame[0].header
    frame_inter_datetime = test_frame_header['DATE-OBS']
    frame_obs_lon = float(test_frame_header['LONGITUD'])
    frame_obs_lat = float(test_frame_header['LATITUDE'])
    #frame_gain = float(test_frame_header['GAIN'])
    #print(frame_obs_lat, frame_obs_lon)
    #mid_ra = float(frame_header["CRVAL1"])
    #mid_dec = float(frame_header["CRVAL2"])

    mid_ra = source_ra
    mid_dec = source_dec

    # Finding minimal positive value of the frame for pyplot
    check_frame = np.where(check_frame < 100000, check_frame, 100000)
    test_ground = np.where(check_frame > 0, check_frame, np.max(check_frame))
    c_min = np.percentile(test_ground, 5)
    frame_max = np.max(check_frame)
    #frame_max = 50000
    print(frame_max)
    background_param = astropy.stats.sigma_clipped_stats(check_frame, sigma=3)
    mean, stddev = background_param[0], background_param[2]
    frame_equiv_gain = stddev ** 2 / mean
    print(frame_equiv_gain, stddev)
    print('Background level: ' + "{}".format(mean))
    #if mean <0 or mean >= frame_max * 0.9:
    #    continue
    if mean >= frame_max * 0.9:
        continue
    if mean <=0:
        check_frame -= 2 * mean
    test_frame = copy(check_frame)
    # Finding 30 brightest stars on the frame (hot pixel may exist in the sample, stars are chosen by their photolevel and SNR)
    bright_stars_y, bright_stars_x, corrected_frame = BrightStarFinder(test_frame, 30)
    centered_frame = copy(corrected_frame)
    #circle_check(bright_stars_x, bright_stars_y, 5, 10, 20, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min * frame_equiv_gain + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    print('Bright stars on frame are picked, analysing them')

    centered_frame = copy(check_frame)
    mean = astropy.stats.sigma_clipped_stats(centered_frame, sigma=1.5)[0]

    # Array of the 6 brightest stars on the frame (with an offset)
    coincidence_check = 0
    while num_of_initial_anchors <= 15 and coincidence_check !=1:
        bright_source_x = []
        bright_source_y = []
        print(num_of_initial_anchors)
        for i in range(offset_from_peak_flux, len(bright_stars_x)):
            if len(bright_source_x) < num_of_initial_anchors:
                discard_check = 0
                if len(bright_source_x) != 0:
                    x_interim, y_interim = bright_stars_x[i], bright_stars_y[i]
                    bright_source_x = np.append(bright_source_x, x_interim)
                    bright_source_y = np.append(bright_source_y, y_interim)
                else:
                    x_interim, y_interim = bright_stars_x[i], bright_stars_y[i]
                    bright_source_x = np.append(bright_source_x, x_interim)
                    bright_source_y = np.append(bright_source_y, y_interim)
            else:
                break
        #centered_frame = copy(corrected_frame)
        #circle_check(bright_source_x, bright_source_y, 5, 10, 20, centered_frame)
        #plt.figure(2)
        #im = plt.imshow(centered_frame)
        #plt.colorbar(im)
        #c_min = mean
        #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
        #plt.show()
        #plt.close()
        # Creating combinations of stars on a frame
        print('Making quads out of frame stars')
        index_arr = np.arange(0, len(bright_source_x))
        index_comb = list(combinations(index_arr, 4))
        comb_num = len(index_comb)
        comb_x = np.zeros((comb_num, 5))
        comb_y = np.zeros((comb_num, 5))

        for i in range(comb_num):
            comb_x[i][0] = i
            comb_y[i][0] = i
            discard_check = 0
            for j in range(4):
                for close_ind in range(4):
                    x_sep = (bright_source_x[index_comb[i][j]] - bright_source_x[index_comb[i][close_ind]]) ** 2
                    y_sep = (bright_source_y[index_comb[i][j]] - bright_source_y[index_comb[i][close_ind]]) ** 2
                    total_sep = sqrt(x_sep + y_sep)
                    if total_sep < 100 and total_sep != 0:
                        discard_check = 1
                        break
                    x_sep = (bright_source_x[index_comb[i][j]] - check_mid_x) ** 2
                    y_sep = (bright_source_y[index_comb[i][j]] - check_mid_y) ** 2
                    total_sep = sqrt(x_sep + y_sep)
                    if total_sep <= 10:
                        discard_check = 1
                        break
                    if abs(bright_source_x[index_comb[i][j]] - check_mid_x) > check_mid_x - 10 or abs(
                            bright_source_y[index_comb[i][j]] - check_mid_y) > check_mid_y - 10:
                        discard_check = 1
                        break
            if discard_check == 0:
                for j in range(4):
                    comb_x[i][j + 1] = bright_source_x[index_comb[i][j]]
                    comb_y[i][j + 1] = bright_source_y[index_comb[i][j]]
            else:
                for j in range(4):
                    comb_x[i][j + 1] = -100
                    comb_y[i][j + 1] = -100
        a_ind = 0
        skip_frame_check = 0
        if comb_x[a_ind][1] == -100:
            while comb_x[a_ind][1] == -100:
                a_ind += 1
                if a_ind >= comb_x.shape[0]:
                    print('Unanalysable frame, probable fault of the hourclock machine, eligible for deletion')
                    list_for_deletion = append(list_for_deletion, filename)
                    skip_frame_check = 1
                    break
        if skip_frame_check == 1:
            continue
        # Finding a coincidence between catalogs
        coincidence_check = 0  # Check if we were able to find a proper coordinate net
        sum_for_min_param = 10 ** 6
        for ind in range(comb_x.shape[0]):
            continuation_check = 0
            if comb_x[ind][1] == -100:
                continue
            bright_frame_quad = quad_frame(comb_x[ind], comb_y[ind])
            for j in range(bright_quad_comb.shape[0]):
                #print('check quad 1', bright_quad_comb[j][5] * rad_conversion * 1900, bright_frame_quad[9])
                #if bright_frame_quad[9] / (bright_quad_comb[j][5] * rad_conversion * 1900) < 0.7 or bright_frame_quad[
                #        9] / (bright_quad_comb[j][5] * rad_conversion * 1900) > 1.5:
                #    print('skipping a quad')
                #    continue
                sum_check = 0  # Check of a number of coincidences
                check_type = 0  # Check of a type of coincidence
                check_arr = zeros(4)  # Array for stocking indexes of already checked stars in the catalog
                pair_check = zeros(2)  # Necessary for discerning the coincident pair
                for i in range(1, 8, 2):
                    test_value = bright_frame_quad[i]  # Value of the frame quad that we test
                    necessary_precis = bright_frame_quad[i + 1]
                    if sum_check == 2:
                        pair_check[0] = 1
                        pair_check[1] = 1
                    for k in range(1, 5):
                        # Skipping the check if we already checked the value or if it's not in the corresponding pair
                        if check_arr[k - 1] == 1:
                            continue
                        if k < 3:
                            pair_index = 0
                        else:
                            pair_index = 1
                        if check_type != 0 and pair_check[pair_index] == 0:
                            continue
                        # Value for checking
                        check_value = bright_quad_comb[j][k]
                        symm_value = 1 - check_value
                        if check_type == 0:
                            if abs(test_value - check_value) <= necessary_precis:
                                check_type = 1
                                sum_check += 1
                                check_arr[k - 1] = 1
                                if k < 3:
                                    pair_check[0] = 1
                                else:
                                    pair_check[1] = 1
                                #print('break e')
                                break
                            if abs(test_value - symm_value) <= necessary_precis:
                                check_type = 2
                                sum_check += 1
                                check_arr[k - 1] = 1
                                if k < 3:
                                    pair_check[0] = 1
                                else:
                                    pair_check[1] = 1
                                #print('break f')
                                break
                        if check_type == 1:
                            if abs(test_value - check_value) <= necessary_precis:
                                sum_check += 1
                                check_arr[k - 1] = 1
                                #print('break g')
                                break
                        if check_type == 2:
                            if abs(test_value - symm_value) <= necessary_precis:
                                sum_check += 1
                                check_arr[k - 1] = 1
                                #print('break h')
                                break
                #print(sum_check, 'aaaa')
                if sum_check >= 4:
                    right_catalog_quad = bright_quad_comb[j]

                    quad_ind = int(right_catalog_quad[0])

                    # Identify stars with the quad between catalog and frame
                    index_correlation_catal_frame = np.zeros(
                        4)  # Index of a cell corresponds to an index of a star in the frame combination - 1,
                    # number in a cell corresponds to an index of a star in the catalog combination

                    check_type = 0
                    correlation_index = 0
                    symmetry_index = 0  # Checks whether farthest stars are mirrored, 0 if not, 1 if yes
                    reversed_index = 0  # Checks whether inner stars are reversed in order, 0 if not, 1 if yes.
                    # Meaning if the first inner star in the catalog is the first in the frame
                    # index is 0. 1 is not
                    for i in range(1):
                        test_value = bright_quad_comb[quad_ind][i + 1]
                        additive = 0
                        pair_ind = 0
                        while abs(test_value - 0.5) <= 0.1:
                            additive += 1
                            test_value = bright_quad_comb[quad_ind][i + 1 + additive]
                        if additive > 1:
                            pair_ind = 1
                        for j in range(1, 8, 2):
                            check_value = bright_frame_quad[j]
                            symm_value = 1 - bright_frame_quad[j]
                            necessary_precis = bright_frame_quad[j + 1]
                            if check_type == 0:
                                if abs(test_value - check_value) <= necessary_precis:
                                    check_type = 1
                                    if pair_ind == 0:
                                        if j <= 3:
                                            reversed_index = 0
                                        else:
                                            reversed_index = 1
                                        #print('break d')
                                        break
                                    else:
                                        if pair_ind == 1:
                                            if j >= 5:
                                                reversed_index = 0
                                            else:
                                                reversed_index = 1
                                if abs(test_value - symm_value) <= necessary_precis:
                                    check_type = 2
                                    symmetry_index = 1
                                    if pair_ind == 0:
                                        if j <= 3:
                                            reversed_index = 0
                                        else:
                                            reversed_index = 1
                                        #print('break c')
                                        break
                                    else:
                                        if pair_ind == 1:
                                            if j >= 5:
                                                reversed_index = 0
                                            else:
                                                reversed_index = 1
                            if check_type == 1:
                                if abs(test_value - check_value) <= necessary_precis:
                                    #print('break a')
                                    break
                            if check_type == 2:
                                if abs(test_value - symm_value) <= necessary_precis:
                                    #print('break b')
                                    break
                    # Placing indexes in the array (correct placements for stars in the frame and the catalog
                    if symmetry_index == 1:
                        index_correlation_catal_frame[int(bright_frame_quad[10]) - 1] = int(bright_quad_comb[quad_ind][7])
                        index_correlation_catal_frame[int(bright_frame_quad[11]) - 1] = int(bright_quad_comb[quad_ind][6])
                    else:
                        index_correlation_catal_frame[int(bright_frame_quad[10]) - 1] = int(bright_quad_comb[quad_ind][6])
                        index_correlation_catal_frame[int(bright_frame_quad[11]) - 1] = int(bright_quad_comb[quad_ind][7])

                    if reversed_index == 0:
                        index_correlation_catal_frame[int(bright_frame_quad[12]) - 1] = int(bright_quad_comb[quad_ind][8])
                        index_correlation_catal_frame[int(bright_frame_quad[13]) - 1] = int(bright_quad_comb[quad_ind][9])
                    else:
                        index_correlation_catal_frame[int(bright_frame_quad[12]) - 1] = int(bright_quad_comb[quad_ind][9])
                        index_correlation_catal_frame[int(bright_frame_quad[13]) - 1] = int(bright_quad_comb[quad_ind][8])
                    index_correlation_catal_frame = index_correlation_catal_frame.astype(int)

                    # Making arrays of corresponding coordinates
                    coord_corresp_x = np.zeros(4)
                    coord_corresp_y = np.zeros(4)
                    coord_corresp_RA = np.zeros(4)
                    coord_corresp_DEC = np.zeros(4)
                    for i in range(1, 5):
                        coord_corresp_x[i - 1] = comb_x[ind][i]
                        coord_corresp_y[i - 1] = comb_y[ind][i]
                        coord_corresp_RA[i - 1] = bright_stars_RA[quad_ind][index_correlation_catal_frame[i - 1]]
                        coord_corresp_DEC[i - 1] = bright_stars_DEC[quad_ind][index_correlation_catal_frame[i - 1]]
                    print('Quad has been picked, determining net parameters')
                    frame_param, initial_min_sum = axis_fit_new(coord_corresp_x, coord_corresp_y, mid_ra, mid_dec,
                                                                coord_corresp_RA,
                                                                coord_corresp_DEC, initial_parameters)
                    #print(coord_corresp_RA, coord_corresp_DEC, coord_corresp_x, coord_corresp_y)
                    #print(right_catalog_quad, bright_frame_quad)
                    if initial_min_sum < sum_for_min_param:
                        min_param = copy(frame_param)
                        sum_for_min_param = initial_min_sum
                    if initial_min_sum / 4 / 2 <= 1:
                        coincidence_check = 1
                        break
            if coincidence_check == 1:
                #print('break')
                break
        if sum_for_min_param / 4 / 2 >= 5:
            print('Skipping the frame because no proper quads', sum_for_min_param)
            num_of_initial_anchors += 1
            continue
        num_of_initial_anchors += 1
    # Creating a rough coordinate net (first version on the frame

    #circle_check(coord_corresp_x, coord_corresp_y, 5, 10, 20, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    frame_param = min_param
    print("Sum of squared deviation for the minimal parameters: " + "{}".format(sum_for_min_param))
    print(frame_param)

    # Adjusting a coordinate net using anchor stars
    # Picking anchors close to the center
    print('Picking anchor stars')
    print(len(bright_anchor_stars['RA_ICRS']))
    inter_anchor_data = {'RA': [], 'DEC': [], 'Gmag': [], 'x': [], 'y': []}  # Series of data for anchor stars
    center_anchor_stars = DataFrame(data=inter_anchor_data)
    for i in range(len(bright_anchor_stars['Gmag'])):
        x_interim, y_interim = axis_fit_use(bright_anchor_stars['RA_ICRS'].values[i],
                                            bright_anchor_stars['DE_ICRS'].values[i], mid_ra, mid_dec, frame_param)
        #circle_ref(x_interim, y_interim, 5, 10, 20, centered_frame)
        dev_x_inter = abs(x_interim - check_mid_x)
        dev_y_inter = abs(y_interim - check_mid_y)

        bright_closeness = 0  # Check whether a star is close to the brightest on the frame
        twin_check = 0  # Check whether a star has a twin already
        if dev_x_inter <= check_mid_x * 0.9 and dev_y_inter <= check_mid_y * 0.9:
            for j in range(len(super_bright_stars['RA_ICRS'])):
                if check_frame[int(x_interim)][int(y_interim)] != 0 and check_frame[int(x_interim) - 20][int(y_interim)] != 0 and check_frame[int(x_interim) + 20][int(y_interim)] != 0 and check_frame[int(x_interim)][int(y_interim) - 20] != 0 and check_frame[int(x_interim)][int(y_interim) + 20] != 0:
                    #print(photometry_snr(x_interim, y_interim, 5, 4, 20, corrected_frame, mean, frame_equiv_gain))
                    if j ==0:
                        if photometry_snr(x_interim, y_interim, 5, 4, 20, centered_frame, mean, frame_equiv_gain) < 6:
                            print('snr too low', photometry_snr(x_interim, y_interim, 5, 4, 20, centered_frame, mean, frame_equiv_gain))
                            bright_closeness = 1
                else:
                    bright_closeness = 1
                    #print('test 2')

            if bright_closeness == 0:
                #print('test 1')
                if len(center_anchor_stars['Gmag']) == 0:
                    inter_series = {'RA': bright_anchor_stars['RA_ICRS'].values[i],
                                    'DEC': bright_anchor_stars['DE_ICRS'].values[i],
                                    'Gmag': bright_anchor_stars['Gmag'].values[i], 'x': x_interim, 'y': y_interim}
                    center_anchor_stars = center_anchor_stars.append(inter_series, ignore_index=True)
                else:
                    for k in range(len(center_anchor_stars['Gmag'])):
                        if abs(x_interim - center_anchor_stars['x'].values[k]) < 2 and abs(
                                y_interim - center_anchor_stars['y'].values[k]) < 2:
                            twin_check = 1
                            break
                    if twin_check != 1:
                        inter_series = {'RA': bright_anchor_stars['RA_ICRS'].values[i],
                                        'DEC': bright_anchor_stars['DE_ICRS'].values[i],
                                        'Gmag': bright_anchor_stars['Gmag'].values[i], 'x': x_interim, 'y': y_interim}
                        center_anchor_stars = center_anchor_stars.append(inter_series, ignore_index=True)
                        if len(center_anchor_stars['RA']) >= 50:
                            break
    #print(len(center_anchor_stars['RA']), 'aaa')
    #circle_check(center_anchor_stars['x'], center_anchor_stars['y'], 5, 10, 20, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    #print(len(center_anchor_stars['RA']), 'aaa')
    dict_stars = []
    array_columns = center_anchor_stars.columns.values.tolist()
    # Adjusting coordinates of anchor stars
    for i in range(len(center_anchor_stars['Gmag'])):
        anchor_snr = photometry_snr(center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i],
                                    centroid_radius, 2, centroid_radius * 2, centered_frame, 0, frame_equiv_gain)
        #print(anchor_snr)
        if anchor_snr <= 4:
            continue
        center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i] = obj_centroid(
            center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i], centroid_radius, 2, centroid_radius * 2, centered_frame, 0)
        anchor_snr = photometry_snr(center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i], centroid_radius, 2, centroid_radius * 2, centered_frame, 0, frame_equiv_gain)
        if anchor_snr >=5:
            dict1 = dict((col_name, center_anchor_stars[col_name].values[i]) for col_name in array_columns)
            dict_stars.append(dict1)
    center_anchor_stars = pd.DataFrame(data = dict_stars, columns=array_columns)
    #circle_check(center_anchor_stars['x'], center_anchor_stars['y'], 5, 10, 20, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    # Estimating parameters of a coordinate net once more (more precise, hopefully)
    print('Number of anchor stars: ' + "{}".format(len(center_anchor_stars['Gmag'])))
    if len(center_anchor_stars['Gmag']) <= 6:
        print('Skipping the frame')
        continue
    inter_frame_param, param_sum = adv_axis_fit_complete(center_anchor_stars['x'], center_anchor_stars['y'], mid_ra,
                                                     mid_dec,
                                                     center_anchor_stars['RA'], center_anchor_stars['DEC'], frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
    if param_sum / len(center_anchor_stars['Gmag']) / 2 >=4:
        inter_frame_param, param_sum = axis_fit_new(center_anchor_stars['x'], center_anchor_stars['y'], mid_ra,
                                                     mid_dec,
                                                     center_anchor_stars['RA'], center_anchor_stars['DEC'], frame_param)
    x_cel = zeros_like(center_anchor_stars['x'])
    y_cel = zeros_like(center_anchor_stars['x'])

    # Making a celestial projection array
    for cel_ind in range(len(x_cel)):
        x_cel[cel_ind] = rad_conversion * (
                cos(center_anchor_stars['DEC'].values[cel_ind] / rad_conversion) * sin(mid_dec / rad_conversion) * cos(
            (center_anchor_stars['RA'].values[cel_ind] - mid_ra) / rad_conversion) - sin(center_anchor_stars['DEC'].values[cel_ind] / rad_conversion) * cos(
            mid_dec / rad_conversion))
        y_cel[cel_ind] = rad_conversion * (
                sin((center_anchor_stars['RA'].values[cel_ind] - mid_ra) / rad_conversion) * cos(center_anchor_stars['DEC'].values[cel_ind] / rad_conversion))
    min_value = min_fun(center_anchor_stars['x'], center_anchor_stars['y'], x_cel, y_cel, min_param)
    #print('test 3', min_value)
    print(inter_frame_param)
    if param_sum / len(center_anchor_stars['Gmag']) / 2 >=6:
        continue
    #for i in range(len(center_anchor_stars['Gmag'])):
    #    test_anchor_stars_x, test_anchor_stars_y = adv_axis_fit_use(center_anchor_stars['RA'].values[i], center_anchor_stars['DEC'].values[i], mid_ra, mid_dec,
    #                                            inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
    #    circle_ref(test_anchor_stars_x, test_anchor_stars_y, 5, 10, 20, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    if param_sum / len(center_anchor_stars['Gmag']) / 2 >=10:
        print('Skipping the frame due to bad coordinate net')
        continue

    # Discarding frames where source is outside the valuable space (where signal is counted)
    source_x_test, source_y_test = adv_axis_fit_use(source_ra, source_dec, mid_ra, mid_dec,
                                                inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
    #print(source_x_test)
    if check_frame[int(source_x_test)][int(source_y_test)] <= 0:
        continue



    # Limiting catalog of bright stars (initial anchors) to those that are inside the frame
    dict_stars = []
    array_columns = bright_stars.columns.values.tolist()
    if frame_num == -1:
        for bright_ind in range(len(bright_stars['Gmag'])):
            bright_catalog_stars_x, bright_catalog_stars_y = adv_axis_fit_use(
                bright_stars['RA_ICRS'].values[bright_ind],
                bright_stars['DE_ICRS'].values[bright_ind], mid_ra, mid_dec,
                inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
            if abs(bright_catalog_stars_x - check_mid_x) <= check_mid_x - 10 and abs(
                    bright_catalog_stars_y - check_mid_y) <= check_mid_y - 10:
                dict1 = dict((col_name, bright_stars[col_name].values[bright_ind]) for col_name in array_columns)
                dict_stars.append(dict1)
            #circle_ref(bright_catalog_stars_x, bright_catalog_stars_y, 10, 5, 25, centered_frame)
        bright_stars = pd.DataFrame(dict_stars, columns=array_columns)

        bright_quad_comb, bright_stars_RA, bright_stars_DEC = quad_catalog(bright_stars)

    # Choosing check refs that are inside valuable space
    dict_check = []
    check_columns = original_check_ref_stars.columns.values.tolist()
    check_len = len(original_check_ref_stars['Gmag'])
    print('Number of check refs: ' + "{}".format(check_len))
    for i in range(check_len):
        original_check_ref_x, original_check_ref_y = adv_axis_fit_use(original_check_ref_stars['RA_ICRS'].values[i],
                                                          original_check_ref_stars['DE_ICRS'].values[i], mid_ra,
                                                          mid_dec,
                                                          inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
       #circle_ref(original_check_ref_x, original_check_ref_y, 5, 10, 20, corrected_frame)
        if abs(original_check_ref_x - check_mid_x) >= check_mid_x * 0.98 or abs(
                original_check_ref_y - check_mid_y) >= check_mid_y * 0.98:
            #print('awooga', original_check_ref_x, original_check_ref_y)
            continue
        if check_frame[original_check_ref_x][original_check_ref_y] != 0 and check_frame[original_check_ref_x][
            original_check_ref_y] != 0:
            check_snr = photometry_snr(original_check_ref_x, original_check_ref_y, centroid_radius, 2, centroid_radius * 2, centered_frame, 0, frame_equiv_gain)
            #print(check_snr)
            #print(check_snr, 'aaa', original_check_ref_x, original_check_ref_y)
            if check_snr >= 5:
                original_check_ref_x, original_check_ref_y = obj_bright_centre(original_check_ref_x, original_check_ref_y, bright_centre_radius, 2, bright_centre_radius * 2, centered_frame,
                                                              0)
                #circle_ref(original_check_ref_x, original_check_ref_y, 3, 4, 10, centered_frame)
                #print(original_check_ref_x, original_check_ref_y)
                check_snr = photometry_snr(original_check_ref_x, original_check_ref_y, centroid_radius, 2, centroid_radius * 2, centered_frame, 0, frame_equiv_gain)
                #print(check_snr)
                max_brightness = pixel_brightness(original_check_ref_x, original_check_ref_y, 5, centered_frame)
                if max_brightness > frame_max * 0.9:
                    print('skipping the check ref', frame_max, max_brightness)
                    continue
                # print(check_x[i], check_y[i], check_snr)
                # print(check_snr)
                if check_snr >= 10:
                    frame_sample = centered_frame[original_check_ref_x - 5:original_check_ref_x + 5][original_check_ref_y - 5:original_check_ref_y + 5]
                    if np.any(frame_sample == 0):
                        continue
                    if abs(original_check_ref_x - check_mid_x) <= check_mid_x * 0.9 and abs(
                            original_check_ref_y - check_mid_y) <= check_mid_y * 0.9:
                        dict1 = dict((col_name, original_check_ref_stars[col_name].values[i]) for col_name in check_columns)
                        dict_check.append(dict1)
                        #circle_ref(original_check_ref_x, original_check_ref_y, 5, 10, 20, centered_frame)
    check_ref_stars = pd.DataFrame(dict_check, columns=check_columns)
    #print(len(dict_check), 'Check refs')
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min * frame_equiv_gain + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    check_len = len(check_ref_stars['Gmag'])
    print('Number of check refs on the frame: ' + "{}".format(check_len))
    if check_len <= 5:
        continue

    every_frame_param = append(every_frame_param, inter_frame_param)

    corrected_photo_frame_path = append(corrected_photo_frame_path, filename)

    dict2 = dict((col_name, check_ref_stars[col_name].values[0]) for col_name in check_columns)
    check_ref_stars = check_ref_stars.drop([0])
    dict_ref.append(dict2)
    ref_star = pd.DataFrame(dict_ref, columns=check_columns)
    check_ref_stars.reset_index(drop=True)
    check_len = len(check_ref_stars['Gmag'])
    #print(check_len, 'aaaaaaaaaaaaaaaaaaa')
    coord_poly_num = int(
        (-3 + sqrt(1 + 8 * (check_len - color_param_num))) / 2 - 1)  # Maximum polynom number
    #print(coord_poly_num, sqrt(1 + 8 * (check_len - color_param_num)), coord_poly_num, (-3 + sqrt(1 + 8 * (check_len - color_param_num))) / 2 - 1)
    if coord_poly_num > max_coord_poly_degr:
        coord_poly_num = max_coord_poly_degr
    elif coord_poly_num < 0:
        coord_poly_num = 0
    # coord_poly_num = 0
    print('Degree of a coordinate polynom: ' + "{}".format(coord_poly_num))
    print('Polynom degree of color function: ' + "{}".format(color_param_num))
    coord_num_arr = np.append(coord_num_arr, coord_poly_num)
    color_num_arr = np.append(color_num_arr, color_param_num)

    # Conducting photometry of fit ref stars
    source_x, source_y = adv_axis_fit_use(source_ra, source_dec, mid_ra, mid_dec,
                                                               inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
    source_x_inter, source_y_inter = obj_centroid(source_x, source_y,
                                                       source_centroid_radius, 3, 10,
                                                       check_frame, 0)
    source_x_inter, source_y_inter = obj_bright_centre(source_x, source_y,
                                                  bright_centre_radius, 2, bright_centre_radius * 2,
                                                  check_frame, 0)
    #centered_frame = copy()
    #circle_ref(source_x, source_y, 4, 3, 10, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    if abs(source_x_inter - source_x) <= 4 and abs(source_y_inter - source_y) <= 4:
        source_x = source_x_inter
        source_y = source_y_inter
    if source_aperture_radius == 0:
        source_photo, source_err, ap_radius, photo_limit = adaptive_photometry_source(
            source_x,
            source_y,
            check_frame, frame_equiv_gain)
        if ap_radius > 6:
            ap_radius = 6
    else:
        ap_radius = source_aperture_radius
        source_photo, source_err, photo_limit = photometry_source(source_x, source_y, ap_radius, 3, 10, check_frame, frame_equiv_gain)

    source_x_arr = np.append(source_x_arr, source_x)
    source_y_arr = np.append(source_y_arr, source_y)
    ap_radius_arr = np.append(ap_radius_arr, ap_radius)
    source_photo_arr = np.append(source_photo_arr, source_photo)
    source_photo_err_arr = np.append(source_photo_err_arr, source_err)
    photo_limit_arr = np.append(photo_limit_arr, photo_limit)
    check_annul = ap_radius + check_gap + 5
    check_ref_x = np.zeros_like(check_ref_stars['RA_ICRS'])
    check_ref_y = np.zeros_like(check_ref_stars['RA_ICRS'])
    check_ref_photo = np.zeros_like(check_ref_stars['RA_ICRS'])
    check_ref_photo_err = np.zeros_like(check_ref_stars['RA_ICRS'])

    for i in range(check_len):
        check_ref_x[i], check_ref_y[i] = adv_axis_fit_use(check_ref_stars['RA_ICRS'].values[i],
                                                                    check_ref_stars['DE_ICRS'].values[i], mid_ra, mid_dec,
                                                                    inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
        check_ref_x[i], check_ref_y[i] = obj_centroid(check_ref_x[i], check_ref_y[i], centroid_radius, 2,
                                                           centroid_radius * 2,
                                                           check_frame, 0)
        check_ref_x[i], check_ref_y[i] = obj_bright_centre(check_ref_x[i], check_ref_y[i], bright_centre_radius, 2, bright_centre_radius * 2,
                                                                    check_frame, 0)

    check_ref_photo, check_ref_photo_err = photometry_check(check_ref_x, check_ref_y, ap_radius, check_gap, check_annul,
                                                            check_frame, frame_equiv_gain)
    #print(check_ref_x[23], check_ref_y[23])
    #print(check_frame.shape[0])
    #print(check_ref_photo, 'ccccccccccccccccccc')
    #test_photo, test_photo_err = photometry_test(check_ref_x[23], check_ref_y[23], ap_radius, check_gap, check_annul,
    #                                      check_frame)
    ref_x, ref_y = adv_axis_fit_use(ref_star['RA_ICRS'].values[frame_ind], ref_star['DE_ICRS'].values[frame_ind],
                                                          mid_ra, mid_dec, inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
    #print(check_ref_photo)
    #circle_check(check_ref_x, check_ref_y, ap_radius, check_gap, check_annul, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
    #plt.clim([40000, 65000])
    #plt.show()
    #plt.close()
    ref_x, ref_y = obj_centroid(ref_x, ref_y, bright_centre_radius, 2, bright_centre_radius * 2,
                                     check_frame, 0)
    ref_x, ref_y = obj_bright_centre(ref_x, ref_y, centroid_radius, 2, centroid_radius * 2,
                                                      check_frame, 0)

    ref_photo_frame = photometry_ref(ref_x, ref_y, ap_radius, ref_gap, ref_annul,
                                          check_frame)
    ref_photo = np.append(ref_photo, ref_photo_frame)
    #print(ref_photo)
    #print(ref_star['Gmag'])
    if ap_radius <= 2:
        mag_ref = np.append(mag_ref, ref_star['Gmag_2'].values[frame_ind])
        check_ref_true_mag = np.copy(check_ref_stars['Gmag_2'])
    elif ap_radius > 2 and ap_radius <= 4:
        mag_ref = np.append(mag_ref, ref_star['Gmag_4'].values[frame_ind])
        check_ref_true_mag = np.copy(check_ref_stars['Gmag_4'])
    elif ap_radius > 4 and ap_radius <= 6:
        mag_ref = np.append(mag_ref, ref_star['Gmag_6'].values[frame_ind])
        check_ref_true_mag = np.copy(check_ref_stars['Gmag_6'])
    elif ap_radius > 6 and ap_radius <= 8:
        mag_ref = np.append(mag_ref, ref_star['Gmag_8'].values[frame_ind])
        check_ref_true_mag = np.copy(check_ref_stars['Gmag_8'])


    check_ref_color = np.copy(check_ref_stars['BP-RP'])
    print(check_ref_photo)

    check_ref_mag, check_ref_mag_err = magnitude_check_ref(mag_ref[frame_ind], ref_photo_frame, check_ref_photo, check_ref_photo_err)
    x_par_sing, x_par_err_sing = sing_frame_full_var_correct(check_ref_mag, check_ref_mag_err,
                                        check_ref_true_mag, check_ref_color, check_ref_x, check_ref_y, frame_shape,
                                        coord_poly_num, color_param_num, max_param_num)
    mag_dev_fit_param = np.append(mag_dev_fit_param, x_par_sing)
    mag_dev_fit_param_err = np.append(mag_dev_fit_param_err, x_par_err_sing)

    dict_check = []
    check_columns = check_ref_stars.columns.values.tolist()
    check_len = len(poss_check_stars['Gmag'])
    for i in range(check_len):
        if poss_check_stars['Gmag'].values[i]==ref_star['Gmag'].values[frame_ind]:
            continue
        check_x, check_y = adv_axis_fit_use(poss_check_stars['RA_ICRS'].values[i],
                                                          poss_check_stars['DE_ICRS'].values[i], mid_ra,
                                                          mid_dec,
                                                          inter_frame_param, frame_inter_datetime, frame_obs_lat, frame_obs_lon)
        # circle_ref(check_x[i], check_y[i], 5, 10, 20, centered_frame)
        if abs(check_x - check_mid_x) >= check_mid_x * 0.9 or abs(
                check_y - check_mid_y) >= check_mid_y * 0.9:
            continue
        max_brightness = pixel_brightness(check_x, check_y, 5, centered_frame)
        #if max_brightness > frame_max * 0.9:
        #    continue
        frame_sample = centered_frame[check_x - 10:check_x + 10][check_y - 10:check_y + 10]
        if np.any(frame_sample == 0):
            continue
        dict1 = dict((col_name, poss_check_stars[col_name].values[i]) for col_name in check_columns)
        dict_check.append(dict1)
    poss_check_stars = pd.DataFrame(dict_check, columns=check_columns)
    print("Number of remaining check stars: " + "{}".format(len(poss_check_stars['Gmag'])))
    frame_ind += 1
    #print(check_stars)

print('Cant use ' + "{}".format(len(list_for_deletion)) + ' frames for being unanalyzable, likely fault with an hourclock machine')
#for del_ind in range(len(list_for_deletion)):
    #remove("{}".format(filepath + list_for_deletion[del_ind]))
# Picking the brightest ref




min_color = arr_min(original_check_ref_stars['BP-RP'])
max_color = arr_max(original_check_ref_stars['BP-RP'])
#print(check_stars)


# Performing photometry
photo_frame_path = copy(corrected_photo_frame_path)
photo_len = len(photo_frame_path)

frame_exp = zeros_like(ref_photo)
frame_tube = []
frame_filter = []
frame_datetime = []


# Poss check stars
poss_check_stars = poss_check_stars.loc[(poss_check_stars['BP-RP'] <= max_color) & (poss_check_stars['BP-RP'] >= min_color)]
check_len = len(poss_check_stars['Gmag'])
poss_check_stars = poss_check_stars.reset_index()
check_photo = zeros((photo_len, check_len))
check_photo_err = zeros((photo_len, check_len))
check_x = zeros((photo_len, check_len))
check_y = zeros((photo_len, check_len))
print(ref_photo)
print(ref_star['Gmag'])

# Readjusting check star arrays
check_ref_len = len(check_ref_stars['Gmag'])
print('Number of check stars: ' + "{}".format(check_len))

check_color = poss_check_stars['BP-RP']
#print(check_stars['Gmag'])
every_frame_param = reshape(every_frame_param, (photo_len, 7))
x_par = reshape(mag_dev_fit_param, (photo_len, max_param_num))
x_par_err = reshape(mag_dev_fit_param_err, (photo_len, max_param_num))
for frame_num in range(photo_len):
    # Parameters of a frame
    test_frame = open("{}".format(filepath + photo_frame_path[frame_num]))  # Opening the frame
    check_frame = getdata("{}".format(filepath + photo_frame_path[frame_num]))  # Getting pixel data
    test_frame_header = test_frame[0].header  # Reading the header
    frame_exp[frame_num] = test_frame_header['EXPTIME']  # Reading exposure time
    frame_inter_datetime = test_frame_header['DATE-OBS']
    inter_datetime = dateutil.parser.parse(frame_inter_datetime)
    inter_JD_date = inter_datetime.toordinal()
    inter_JD_date_datetime = datetime.datetime.fromordinal(inter_JD_date)
    inter_JD_date += 1721424.5
    inter_JD_date_full = inter_JD_date + (inter_datetime - inter_JD_date_datetime).total_seconds() / 86400
    frame_datetime = append(frame_datetime, inter_JD_date_full + frame_exp[frame_num] / 2 / 60 / 60 / 24)
    frame_obs_lon = float(test_frame_header['LONGITUD'])
    frame_obs_lat = float(test_frame_header['LATITUDE'])
    #frame_tube = append(frame_tube, 'Vostryakovo')  # Reading tube
    #frame_filter = append(frame_filter, 'Vostryakovo')  # Reading filter
    frame_tube = append(frame_tube, test_frame_header['AXIS'])  # Reading tube
    frame_filter = append(frame_filter, test_frame_header['FILTER'])  # Reading filter


    # Coordinates of source, ref, checks

    # Making a mask array and erasing the close bright star
    if mask_check != 0:
        check_masks_arr = zeros((check_len, 20, 20))
        check_mask = zeros((20, 20))
        for i in range(check_len):
            check_ref_x[frame_num][i], check_ref_y[frame_num][i], check_masks_arr[i] = check_mask_center(check_ref_x[frame_num][i],
                                                                                                 check_ref_y[frame_num][i],
                                                                                                 3,
                                                                                                 5,
                                                                                                 15,
                                                                                                 check_frame)
            check_mask += check_masks_arr[i]
        check_mask = check_mask / max(check_mask)
        check_peak_x, check_peak_y = obj_peak(10, 10, 10, check_mask)

        # Erasing the bright star near GRB
        near_star_x, near_star_y = obj_peak(512, 512, 15, check_frame)
        near_star_peak = check_frame[near_star_x][near_star_y] - back_level(512, 512, 10, 20, check_frame)
        d_x, d_y = near_star_x - check_peak_x, near_star_y - check_peak_y
        for x_ind in range(d_x, d_x + 20):
            for y_ind in range(d_y, d_y + 20):
                check_frame[x_ind][y_ind] -= near_star_peak * check_mask[x_ind - d_x][y_ind - d_y]

    # The photometry


    check_annul = ap_radius_arr[frame_num] + check_gap + 5


    # Checks
    for i in range(check_len):
        check_x[frame_num][i], check_y[frame_num][i] = adv_axis_fit_use(poss_check_stars['RA_ICRS'].values[i],
                                                                                poss_check_stars['DE_ICRS'].values[i],
                                                                                mid_ra, mid_dec,
                                                                                every_frame_param[frame_num], frame_inter_datetime, frame_obs_lat, frame_obs_lon)
        check_snr = photometry_snr(check_x[frame_num][i], check_y[frame_num][i], 6, 2, 15, check_frame, 0, frame_equiv_gain)
        if check_snr >=3:
            check_x[frame_num][i], check_y[frame_num][i] = obj_centroid(check_x[frame_num][i],
                                                                             check_y[frame_num][i], centroid_radius, 2,
                                                                             centroid_radius * 2,
                                                                             check_frame, 0)
            check_x[frame_num][i], check_y[frame_num][i] = obj_bright_centre(check_x[frame_num][i],
                                                                                check_y[frame_num][i], centroid_radius, 2, centroid_radius * 2,
                                                                                check_frame, 0)
    #mean = astropy.stats.sigma_clipped_stats(check_frame, sigma=1.5)[0]
    #centered_frame = copy(check_frame)
    # circle_check(check_ref_x[frame_num], check_ref_y[frame_num], 5, 10, 20, centered_frame)
    #circle_check(check_x[frame_num], check_y[frame_num], 5, 10, 20, centered_frame)
    # circle_ref(ref_x[frame_num], ref_y[frame_num], 5, 10, 20, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
    #plt.show()
    #plt.close()
    check_photo[frame_num], check_photo_err[frame_num] = photometry_check(check_x[frame_num], check_y[frame_num],
                                                                          ap_radius_arr[frame_num], check_gap,
                                                                          check_annul, check_frame, frame_equiv_gain)
    if photo_frame_path[frame_num] == 'suck':
        #print(check_ref_x[frame_num])
        #print(check_ref_y[frame_num])
        mean = astropy.stats.sigma_clipped_stats(check_frame, sigma=1.5)[0]
        centered_frame = copy(check_frame)
        #circle_check(check_ref_x[frame_num], check_ref_y[frame_num], 5, 10, 20, centered_frame)
        circle_check(check_x[frame_num], check_y[frame_num], 5, 10, 20, centered_frame)
        #circle_ref(ref_x[frame_num], ref_y[frame_num], 5, 10, 20, centered_frame)
        plt.figure(2)
        im = plt.imshow(centered_frame)
        plt.colorbar(im)
        c_min = mean
        plt.clim([c_min, c_min + 8 * sqrt(c_min + readout_noise ** 2)])
        plt.show()
        plt.close()
    print(photo_frame_path[frame_num])
    print('Radius of adaptive aperture: ' + "{}".format(ap_radius))




# Determining magnitude
source_mag, source_mag_signal_err, check_mag = magnitude(mag_ref, ref_photo, source_photo_arr, source_photo_err_arr, check_photo, check_photo_err)

limit_mag = magnitude_obj(mag_ref, ref_photo, photo_limit_arr)

# Determining magnitude errors for source and checks
source_mag_check_var_err, check_mag_var_err, check_mag_signal_err, source_check_var_num, check_average_mag, source_bright_check_var_err = adv_magnitude_err(source_mag,
                                                                                                   check_mag,
                                                                                                   check_photo,
                                                                                                   check_photo_err)

# Finding fit parameters for all the frames

source_mag_corr = zeros_like(source_mag)


source_mag_corr = source_mag_correct(source_mag, source_color, source_x_arr, source_y_arr, x_par, frame_shape, coord_num_arr, color_num_arr)

limit_mag_corr = source_mag_correct(limit_mag, source_color, source_x_arr, source_y_arr, x_par, frame_shape, coord_num_arr, color_num_arr)

check_mag_corr = check_mag_correct(check_mag, check_color, check_x, check_y, x_par, frame_shape, coord_num_arr, color_num_arr)

source_mag_check_corr_var_err, check_mag_corr_var_err, check_mag_signal_err, source_mag_check_corr_var_num, check_corr_average_mag, source_mag_bright_check_corr_var_err = adv_magnitude_err(source_mag_corr,
                                                                                                   check_mag_corr,
                                                                                                   check_photo,
                                                                                                   check_photo_err)


# Estimating source mag error from check signal noise
#source_mag_check_signal_err = np.zeros_like(source_mag)
#signal_err_const = 0
#for k in range(photo_len):
#    check_vis_num = 0
#    for i in range(check_len):
#        if np.isnan(check_mag_signal_err[k][i]) == 0:
#            signal_err_const += check_mag_signal_err[k][i] ** 2
#            check_vis_num += 1
#    source_mag_check_signal_err[k] = math.sqrt(signal_err_const / (check_vis_num - 1))


# Determining deviation from true mag for every ref
dev_data = {"Frame": photo_frame_path}






# Fit data in a

fit = {"Frame": photo_frame_path}
fit_data = pd.DataFrame(data=fit, index=frame_datetime)
coord_par_num = (x_par.shape[1] - 1 - color_param_num) / 2
for m in range(x_par.shape[1]):
    export_x_par = np.zeros(x_par.shape[0])
    export_x_par_err = np.zeros(x_par.shape[0])
    for k in range(x_par.shape[0]):
        export_x_par[k] = x_par[k][m]
        export_x_par_err[k] = x_par_err[k][m]
    if m <= color_param_num:
        export_column = "Color parameter " + "{}".format(m)
        export_err_column = "Color parameter " + "{}".format(m) + "error"
    if m > color_param_num and m <= coord_par_num + color_param_num:
        export_column = "X coord parameter " + "{}".format(m - color_param_num)
        export_err_column = "X coord parameter " + "{}".format(m - color_param_num) + "error"
    elif m > coord_par_num + color_param_num:
        export_column = "Y coord parameter " + "{}".format(m - coord_par_num - color_param_num)
        export_err_column = "Y coord parameter " + "{}".format(m - coord_par_num - color_param_num) + "error"
    df = {export_column: export_x_par,
          export_err_column: export_x_par_err}
    df_export = pd.DataFrame(data=df, index=frame_datetime)
    fit_data = fit_data.join(df_export)
export_column = "Coord polynom max degree for a frame"
export_err_column = "Color num for a frame"
df = {export_column: coord_num_arr,
      export_err_column: color_num_arr}
df_export = pd.DataFrame(data=df, index=frame_datetime)
fit_data = fit_data.join(df_export)

fit_data.index.name = 'Julian Days, middle of exposure'

# Exporting results of the photometry into excel in the same folder
data = {"Frame": photo_frame_path, 'Tube': frame_tube, 'Filter': frame_filter,
        'Frame exposure time, seconds': frame_exp, "Source mag (corrected)": source_mag_corr,
        '3 sigma limit, mag': limit_mag,
        '3 sigma corrected limit, mag': limit_mag_corr,
        'Source mag error (variability of checks, not corrected)': source_mag_check_var_err,
        'Source mag error (variability of bright checks, not corrected)': source_bright_check_var_err,
        'Source mag error (variability of checks, corrected)': source_mag_check_corr_var_err,
        'Source mag error (variability of bright checks, corrected)': source_mag_bright_check_corr_var_err,
        "Source mag error (signal noise)": source_mag_signal_err,
#        "Source mag error (signal noise of checks)": source_mag_check_signal_err,
        "Source x coord": source_x_arr, "Source y coord": source_y_arr}

full_data = {"Frame": photo_frame_path, 'Tube': frame_tube, 'Filter': frame_filter,
             'Frame exposure time, seconds': frame_exp, 'Source mag': source_mag,
             '3 sigma limit, mag': limit_mag,
             '3 sigma corrected limit, mag': limit_mag_corr,
             "Source mag (corrected)": source_mag_corr,
             'Source mag error (variability of checks, not corrected)': source_mag_check_var_err,
             'Source mag error (variability of bright checks, not corrected)': source_bright_check_var_err,
             'Source mag error (variability of checks, corrected)': source_mag_check_corr_var_err,
             'Source mag error (variability of bright checks, corrected)': source_mag_bright_check_corr_var_err,
             "Source mag error (signal noise)": source_mag_signal_err,
             #"Source mag error (signal noise of checks)": source_mag_check_signal_err,
             "Source x coord": source_x, "Source y coord": source_y}

useful_data = {"Frame": photo_frame_path, 'Tube': frame_tube, 'Filter': frame_filter,
             'Frame exposure time, seconds': frame_exp, 'Source mag': source_mag,
             '3 sigma limit, mag': limit_mag,
             '3 sigma corrected limit, mag': limit_mag_corr,
             "Source mag (corrected)": source_mag_corr,
             'Source mag error (variability of checks, not corrected)': source_mag_check_var_err,
             'Source mag error (variability of bright checks, not corrected)': source_bright_check_var_err,
             'Source mag error (variability of checks, corrected)': source_mag_check_corr_var_err,
             'Source mag error (variability of bright checks, corrected)': source_mag_bright_check_corr_var_err,
             "Source mag error (signal noise)": source_mag_signal_err,
             #"Source mag error (signal noise of checks)": source_mag_check_signal_err,
             "Source x coord": source_x, "Source y coord": source_y}

photo_data = pd.DataFrame(data=data, index=frame_datetime)
full_photo_data = pd.DataFrame(data=full_data, index=frame_datetime)
useful_photo_data = pd.DataFrame(data=full_data, index=frame_datetime)
photo_data.index.name = 'Julian Days, middle of exposure'
full_photo_data.index.name = 'Julian Days, middle of exposure'
useful_photo_data.index.name = 'Julian Days, middle of exposure'

# Adding checks
for k in range(check_mag.shape[1]):
    check_mag_data = np.zeros_like(source_mag)
    check_mag_corr_data = np.zeros_like(source_mag)
    check_mag_var_err_data = np.zeros_like(source_mag)
    check_mag_corr_var_err_data = np.zeros_like(source_mag)
    check_mag_signal_err_data = np.zeros_like(source_mag)
    for j in range(check_mag.shape[0]):
        check_mag_data[j] = check_mag[j][k]
        check_mag_corr_data[j] = check_mag_corr[j][k]
        check_mag_var_err_data[j] = check_mag_var_err[j][k]
        check_mag_corr_var_err_data[j] = check_mag_corr_var_err[j][k]
        check_mag_signal_err_data[j] = check_mag_signal_err[j][k]
    df = {'Check star ' + "{}".format(k + 1) + ' mag': check_mag_data,
          'Check star ' + "{}".format(k + 1) + 'corrected mag': check_mag_corr_data,
          'Check star ' + "{}".format(k + 1) + ' mag error (variability, not corrected)': check_mag_var_err_data,
          'Check star ' + "{}".format(k + 1) + ' mag error (variability, corrected)': check_mag_corr_var_err_data,
          'Check star ' + "{}".format(k + 1) + ' mag error (signal noise)': check_mag_signal_err_data}
    df2 = pd.DataFrame(data=df, index=frame_datetime)
    full_photo_data = full_photo_data.join(df2)
    useful_photo_data = useful_photo_data.join(df2)

# Getting rid of bad frames in useful data
file_col = useful_photo_data.columns.values.tolist()
photo_dict = []
frame_datetime = np.empty(0)
for frame_num in range(len(useful_photo_data['Tube'])):
    if useful_photo_data['Source mag error (variability of checks, not corrected)'].values[frame_num] < 10:
        dict1 = dict((col_name, useful_photo_data[col_name].values[frame_num]) for col_name in file_col)
        frame_datetime = np.append(frame_datetime, useful_photo_data.index[frame_num])
        photo_dict.append(dict1)
useful_photo_data = pd.DataFrame(data=photo_dict, columns=file_col, index=frame_datetime)
useful_photo_data.index.name = 'Julian Days, middle of exposure'
#print(useful_photo_data)
# Resetting mag errors
useful_photo_data = export_magnitude_correct(useful_photo_data)


test_frame = open("{}".format(filepath + filename))
hdr = test_frame[0].header
#hdr['AXIS'] = 'Vostryakovo'
#hdr['FILTER'] = 'Vostryakovo'
full_photo_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(
    hdr['FILTER']) + ' full photometry data (source and checks).xlsx'
useful_photo_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(
    hdr['FILTER']) + ' useful photometry data (source and checks).xlsx'
photo_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' photometry data (source).xlsx'
fit_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' fit data.xlsx'
first_frame_dev_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' first frame deviation data.xlsx'
deviation_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' deviation data.xlsx'
check_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' check data.xlsx'
photo_data.to_excel("{}".format(photometry_path + photo_file))
full_photo_data.to_excel("{}".format(photometry_path + full_photo_file))
fit_data.to_excel("{}".format(photometry_path + fit_file))
poss_check_stars.to_excel("{}".format(photometry_path + check_file))
useful_photo_data.to_excel("{}".format(photometry_path + useful_photo_file))

print("Exporting over")
print('Program took ' + "{}".format((time() - start_time) / 60) + ' minutes')
