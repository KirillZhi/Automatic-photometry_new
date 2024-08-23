import pandas as pd
from numpy import zeros, empty, append, arange, max, absolute, array, delete

from pandas import DataFrame
#import photutils
from itertools import combinations

from math import sqrt
from frame_functions import *
from coord_functions import *
import dateutil
import datetime
from astropy.io.fits import getdata, open
from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename

from os import listdir
from os.path import isfile, join

from re import search
from matplotlib import pyplot as plt
from time import time

# Parameters of the initial fit (parameters of anchors)
# print(combined_magnitude(16.5, 13))
num_of_initial_anchors = 7          # Answers for anchors of the initial net
offset_from_peak_flux = 0           # Answers for the offset of initial anchors from the peak flux stars in both catalog and frame
# 168.57846 -36.74690
# source_ra, source_dec = 311.0217083, -47.663472
source_ra, source_dec = 86.566791666667, 23.244722222222
 # Celestial (equatorial) coordinates of the studied source
source_centroid_radius = 1          # Radius of the aperture that is used to calculate centroid (multiply by 3 to get correct radius)
source_aperture_radius = 0          # Change from 0 to any positive integer less than 7, it will serve as an aperture radius
# source_ra, source_dec = 319.71831, 33.85044

# source_ra, source_dec = 243.91836, 14.39909
# 168.57846 -36.74690

source_color = 0.96                 # Color BP-RP of the studied source
search_radius = 14.                 # Radius in arcminutes where we look for checks and anchors in the catalog
check_mag_up = 10                   # Upper limit on the flux of checks
mag_lower = 15                      # Lower limit on the flux of checks and anchors
anchor_mag_up = 10                  # Upper limit on the flux of anchors
desired_ref = check_mag_up + 2      # Lower limit on the flux of reference
# Filepath for fits frames
readout_noise = 9                   # Readout noise of camera
root = Tk()
root.withdraw()
print("Choose directory with frames")
filepath = askdirectory(initialdir="D:/Photometry/")
# photometry_path = askdirectory(initialdir="F:/Astronomy data/")
# catalog_path = askopenfilename(initialdir="F:/Astronomy data/")
print("Choose directory for exporting data")
photometry_path = askdirectory(initialdir=filepath)
catalog_path = askopenfilename(initialdir=photometry_path)
# print(angular_separation(source_ra, source_dec, 168.22493680432, -36.78760242461) * 60)
# print(angular_separation(1., -1., 0., 1.))
filepath = filepath + "/"
photometry_path = photometry_path + "/"
root.update()
root.destroy()

# Making quads from the catalogue
bright_stars, interm_bright_stars, bright_anchor_stars, catalog_stars, super_bright_stars= gaia_bright(
    catalog_path,
    num_of_initial_anchors,
    offset_from_peak_flux,
    search_radius,
    source_ra,
    source_dec, check_mag_up, anchor_mag_up, mag_lower)
print("Number of anchors")
print(len(bright_anchor_stars['Gmag']))
print(len(bright_stars['Gmag']))
print(bright_stars)
bright_quad_comb, bright_stars_RA, bright_stars_DEC = quad_catalog(bright_stars)
# print(bright_anchor_stars)
#print(bright_stars)
# print(bright_stars['Gmag'])

# Picking reference and checks

check_stars, ref_star = check_ref_search(interm_bright_stars, desired_ref, source_ra, source_dec)
#print(ref_star)
#print(check_stars['Gmag'])
print("Number of checks")
print(len(check_stars['Gmag']))
# print(check_stars['Gmag_2'] - check_stars['Gmag_4'])
check_true_mag = check_stars['Gmag']
# print(check_stars)
check_color = check_stars['BP-RP']
# print(check_color)

# mag_ref = ref_star['Gmag'].values[0]
print('Catalog analyzed, beginning operations with the frames')

# Coordinates of a test star
tester_ra = 12.21290956884
tester_dec = -56.04868649664
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

print('Frames are picked, building the coordinate net')
print('Approximate time till completion: ' + "{}".format(len(photo_frame_path) * 10. / 60) + ' minutes')
start_time = time()
# Photometry parameters
ref_gap = 8
ref_annul = 20
check_gap = 8
check_annul = 20
photo_len = len(photo_frame_path)
ref_photo = zeros(photo_len)
frame_datetime = zeros_like(ref_photo)
frame_exp = zeros_like(ref_photo)
frame_tube = []
frame_filter = []

source_x = zeros(photo_len)
source_y = zeros(photo_len)

ref_x = zeros(photo_len)
ref_y = zeros(photo_len)
check_len = len(check_stars['Gmag'])
mag_ref = zeros(photo_len)
check_photo = zeros((photo_len, check_len))
check_photo_err = zeros((photo_len, check_len))
check_x = zeros((photo_len, check_len))
check_y = zeros((photo_len, check_len))
source_photo = zeros(photo_len)
source_err = zeros(photo_len)

every_frame_param = zeros((photo_len, 12))
initial_parameters = zeros(12)
initial_parameters[2] = 1000
initial_parameters[4] = 5
for frame_num in range(len(photo_frame_path)):
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
    mid_ra = float(frame_header["CRVAL1"])
    mid_dec = float(frame_header["CRVAL2"])
    # mid_ra = mid_ra.replace(' ', ':')
    # mid_dec = mid_dec.replace(' ', ':')
    # mid_cel_coord = mid_ra + ' ' + mid_dec
    # mid_cel_coord = SkyCoord(mid_cel_coord, frame='icrs', unit=(u.hourangle, u.deg))
    # mid_ra = mid_cel_coord.ra.deg
    # mid_dec = mid_cel_coord.dec.deg
    # Minimum value of the frame (positive)
    test_ground = np.where(check_frame > 0, check_frame, np.max(check_frame))
    c_min = np.percentile(test_ground, 5)
    frame_max = np.max(check_frame)
    # a, b = np.where(check_frame >= frame_max * 0.5)
    # print(a)
    # Looking for sources in the frame
    # print('eeee')
    mean = astropy.stats.sigma_clipped_stats(check_frame, sigma=1.5)[0]
    # print('ggggg')
    test_frame = copy(check_frame)

    bright_stars_y, bright_stars_x, corrected_frame = BrightStarFinder(test_frame, 30)
    print('Bright stars on frame are picked, analysing them')
    # print(bright_stars_x, bright_stars_y)
    # print(len(bright_stars_y))
    centered_frame = copy(corrected_frame)
    mean = astropy.stats.sigma_clipped_stats(centered_frame, sigma=1.5)[0]
    #circle_check(bright_stars_x, bright_stars_y, 15, 5, 25, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    # daofind = photutils.detection.DAOStarFinder(fwhm=3.0, threshold=15. * std, sigma_radius=3, exclude_border=True)
    # sources = daofind(check_frame - mean)

    # Picking 6 brightest sources (a large portion of code)
    # sorted_sources = sources.group_by('peak')
    # Precise photometry of the sources

    # discard_check = 0
    # print(photometry_ref(359, 319, 5, 5, 20, check_frame))
    # Correcting fluxes and centers
    # outer_end = len(sorted_sources['flux']) - 40 - offset_from_peak_flux
    # if len(sorted_sources['flux']) - 40 - offset_from_peak_flux < 1:
    #    outer_end = 1
    # for i in range(len(sorted_sources['flux']) - 1, outer_end, -1):
    #    sorted_sources['ycentroid'][i], sorted_sources['xcentroid'][i] = obj_centroid(sorted_sources['ycentroid'][i],
    #                                                                                  sorted_sources['xcentroid'][i], 5,
    #                                                                                  5,
    #                                                                                  20, check_frame)

    # if abs(sorted_sources['ycentroid'][i] - 359) < 8 and abs(sorted_sources['xcentroid'][i] - 319) < 8:
    #    print(sorted_sources['flux'][i])
    #    sorted_sources['flux'][i] = photometry_ref(sorted_sources['ycentroid'][i], sorted_sources['xcentroid'][i], 5,
    #                                               5,
    #                                               20, check_frame)
    # if abs(sorted_sources['ycentroid'][i] - 359) < 8 and abs(sorted_sources['xcentroid'][i] - 319) < 8:
    # print(sorted_sources['flux'][i])
    # sorted_sources = sorted_sources.group_by('flux')

    # Discarding twins of stars
    # outer_end = len(sorted_sources['flux']) - 40 - offset_from_peak_flux
    # if len(sorted_sources['flux']) - 40 - offset_from_peak_flux < 1:
    #    outer_end = 1
    # for i in range(len(sorted_sources['flux']) - 1, outer_end, -1):
    #    discard_check = 0
    #    # first_star = 0
    #    for j in range(len(sorted_sources['flux']) - 1, outer_end, -1):
    #        x_sep = (sorted_sources['xcentroid'][j] - sorted_sources['xcentroid'][i]) ** 2
    #        y_sep = (sorted_sources['ycentroid'][j] - sorted_sources['ycentroid'][i]) ** 2
    #        total_sep = sqrt(x_sep + y_sep)
    #        if total_sep < 4 and i < j:
    #            discard_check = 1
    #            break
    #        # if total_sep < 4 and first_star == 0 and i != j:
    #        #    first_star = 1
    #    if discard_check == 1:
    #        sorted_sources['flux'][i] = -100
    # if abs(sorted_sources['ycentroid'][i] - 359) < 8 and abs(sorted_sources['xcentroid'][i] - 319) < 8:
    # print(sorted_sources['flux'][i])
    # print(i)

    # sorted_sources = sorted_sources.group_by('flux')
    # centered_frame = copy(check_frame)
    # centered_frame = copy(check_frame)
    # circle_check(sorted_sources['ycentroid'][len(sorted_sources['flux']) - 1:len(sorted_sources['flux']) - 20:-1], sorted_sources['xcentroid'][len(sorted_sources['flux']) - 1:len(sorted_sources['flux']) - 20:-1], 15, 5, 25, centered_frame)
    # print(sorted_sources)
    # plt.figure(2)
    # im = plt.imshow(centered_frame)
    # plt.colorbar(im)
    # c_min = mean
    # plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    # plt.show()
    # plt.close
    # Array of the 6 brightest stars on the frame (with an offset)
    # print(bright_stars_x[1], bright_stars_y[1])
    bright_source_x = []
    bright_source_y = []
    # outer_end = len(sorted_sources['flux']) - 20 - offset_from_peak_flux
    # if len(bright_stars_x) - 20 - offset_from_peak_flux < 1:
    #    outer_end = 1
    for i in range(offset_from_peak_flux, len(bright_stars_x)):
        if len(bright_source_x) < num_of_initial_anchors:
            discard_check = 0
            if len(bright_source_x) != 0:
                x_interim, y_interim = bright_stars_x[i], bright_stars_y[i]
                bright_source_x = np.append(bright_source_x, x_interim)
                bright_source_y = np.append(bright_source_y, y_interim)
            else:
                # x_interim, y_interim = obj_centroid(sorted_sources['ycentroid'][i], sorted_sources['xcentroid'][i], 5, 10,
                #                                    20, check_frame)
                x_interim, y_interim = bright_stars_x[i], bright_stars_y[i]
                bright_source_x = np.append(bright_source_x, x_interim)
                bright_source_y = np.append(bright_source_y, y_interim)
        else:
            break
    # centered_frame = copy(check_frame)
    #centered_frame = copy(corrected_frame)
    #circle_check(bright_source_x, bright_source_y, 15, 5, 25, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    #Creating combinations of stars on a frame
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
    centered_frame = np.copy(check_frame)
    a_ind = 0
    if comb_x[a_ind][1] == -100:
        while comb_x[a_ind][1] == -100:
            a_ind += 1
    # print(comb_x[a_ind], comb_y[a_ind])
    # print(quad_frame(comb_x[a_ind], comb_y[a_ind]))
    # print(first_cel_comb)
    # print(first_cel_ra)
    # print(first_cel_dec)
    # circle_ref(bright_source_x[index_comb[a_ind][0]], bright_source_y[index_comb[a_ind][0]], 15, 5, 25, centered_frame)
    # circle_ref(bright_source_x[index_comb[a_ind][1]], bright_source_y[index_comb[a_ind][1]], 15, 5, 25, centered_frame)
    # circle_ref(bright_source_x[index_comb[a_ind][2]], bright_source_y[index_comb[a_ind][2]], 15, 5, 25, centered_frame)
    # circle_ref(bright_source_x[index_comb[a_ind][3]], bright_source_y[index_comb[a_ind][3]], 15, 5, 25, centered_frame)
    #circle_check(comb_x[a_ind][1:5], comb_y[a_ind][1:5], 15, 5, 25, centered_frame)
    #circle_check(sorted_sources['ycentroid'][len(sorted_sources['flux']) - 1:len(sorted_sources['flux']) - 1 - 6:-1], sorted_sources['xcentroid'][len(sorted_sources['flux']) - 1:len(sorted_sources['flux']) - 1 - 6:-1], 15, 5, 25, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close
    # Finding a coincidence between catalogs
    # for ind in range(comb_x.shape[0]):
    #    continuation_check = 0
    #    if comb_x[ind][1] == -100:
    #        continue
    #    bright_frame_quad = quad_frame(comb_x[ind], comb_y[ind])
    #    for j in range(bright_quad_comb.shape[0]):
    #        #print(bright_frame_quad[9], bright_quad_comb[j][5] * rad_conversion * 1930)
    #        if bright_frame_quad[9] / (bright_quad_comb[j][5] * rad_conversion * 1900) < 0.9 or bright_frame_quad[9] / (bright_quad_comb[j][5] * rad_conversion * 1900) > 1.025:
    #            continue
    #        #else:
    #        #    print('yeeehaaaaw')
    #        sum_check = 0  # Check of a number of coincidences
    #        check_type = 0  # Check of a type of coincidence
    #        check_arr = zeros(4)  # Array for stocking indexes of already checked stars in the catalog
    #        for i in range(1, 8, 2):
    #            test_value = bright_frame_quad[i]  # Value of the frame quad that we test
    #            necessary_precis = bright_frame_quad[i + 1]
    #            for k in range(1, 5):
    #                if check_arr[k - 1] == 1:
    #                    continue
    #                # Value for checking
    #                check_value = bright_quad_comb[j][k]
    #                symm_value = 1 - check_value
    #                if check_type == 0:
    #                    if abs(test_value - check_value) <= necessary_precis:
    #                        check_type = 1
    #                        sum_check += 1
    #                        check_arr[k - 1] = 1
    #                        break
    #                    if abs(test_value - symm_value) <= necessary_precis:
    #                        check_type = 2
    #                        sum_check += 1
    #                        check_arr[k - 1] = 1
    #                        break
    #                if check_type == 1:
    #                    if abs(test_value - check_value) <= necessary_precis:
    #                        sum_check += 1
    #                        check_arr[k - 1] = 1
    #                        break
    #                if check_type == 2:
    #                    if abs(test_value - symm_value) <= necessary_precis:
    #                        sum_check += 1
    #                        check_arr[k - 1] = 1
    #                        break
    #        if sum_check == 4:
    #            right_catalog_quad = bright_quad_comb[j]
    #            #if
    #            break
    #    if sum_check == 4:
    #        break

    # Finding a coincidence between catalogs
    coincidence_check = 0       # Check if we were able to find a proper coordinate net
    sum_for_min_param = 10 ** 6
    for ind in range(comb_x.shape[0]):
        continuation_check = 0
        if comb_x[ind][1] == -100:
            continue
        bright_frame_quad = quad_frame(comb_x[ind], comb_y[ind])
        for j in range(bright_quad_comb.shape[0]):
            # print(bright_frame_quad[9], bright_quad_comb[j][5] * rad_conversion * 1930)
            if bright_frame_quad[9] / (bright_quad_comb[j][5] * rad_conversion * 1900) < 0.8 or bright_frame_quad[
                9] / (bright_quad_comb[j][5] * rad_conversion * 1900) > 1.2:
                continue
            # else:
            #    print('yeeehaaaaw')
            sum_check = 0  # Check of a number of coincidences
            check_type = 0  # Check of a type of coincidence
            check_arr = zeros(4)  # Array for stocking indexes of already checked stars in the catalog
            pair_check = zeros(2)  # Necessary for discerning the coincident pair
            for i in range(1, 8, 2):
                test_value = bright_frame_quad[i]  # Value of the frame quad that we test
                necessary_precis = bright_frame_quad[i + 1]
                if sum_check == 2:
                    # pair_check = pair_check.fill(pair_check, 1)
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
                            break
                        if abs(test_value - symm_value) <= necessary_precis:
                            check_type = 2
                            sum_check += 1
                            check_arr[k - 1] = 1
                            if k < 3:
                                pair_check[0] = 1
                            else:
                                pair_check[1] = 1
                            break
                    if check_type == 1:
                        if abs(test_value - check_value) <= necessary_precis:
                            sum_check += 1
                            check_arr[k - 1] = 1
                            break
                    if check_type == 2:
                        if abs(test_value - symm_value) <= necessary_precis:
                            sum_check += 1
                            check_arr[k - 1] = 1
                            break
            #if sum_check >=4 and ind <=6:
            #    print(sum_check)
            if sum_check >= 4:
                right_catalog_quad = bright_quad_comb[j]
                # if
                #break
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
                                # print(j)
                                if pair_ind == 0:
                                    if j <= 3:
                                        reversed_index = 0
                                    else:
                                        reversed_index = 1
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
                                # print(j)
                                if pair_ind == 0:
                                    if j <= 3:
                                        reversed_index = 0
                                    else:
                                        reversed_index = 1
                                    break
                                else:
                                    if pair_ind == 1:
                                        if j >= 5:
                                            reversed_index = 0
                                        else:
                                            reversed_index = 1
                        if check_type == 1:
                            # print(j)
                            if abs(test_value - check_value) <= necessary_precis:
                                break
                        if check_type == 2:
                            # print(j)
                            if abs(test_value - symm_value) <= necessary_precis:
                                break
                # Placing indexes in the array
                # print(symmetry_index, reversed_index)
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
                frame_param, initial_min_sum = axis_fit_new(coord_corresp_x, coord_corresp_y, check_path, coord_corresp_RA,
                                           coord_corresp_DEC, initial_parameters)
                if initial_min_sum < sum_for_min_param:
                    min_param = copy(frame_param)
                    sum_for_min_param = initial_min_sum
                if initial_min_sum / 4 / 2 <= 1:
                    coincidence_check = 1
                    break
                #print(coord_corresp_x, coord_corresp_y)
                #print(coord_corresp_RA, coord_corresp_DEC)
                # print(bright_stars_RA[quad_ind], bright_stars_DEC[quad_ind])
                #print(right_catalog_quad)
                #print(bright_frame_quad)
                #centered_frame = np.copy(check_frame)
                #circle_check(coord_corresp_x, coord_corresp_y, 10, 5, 20, centered_frame)
                #plt.figure(2)
                #im = plt.imshow(centered_frame)
                #plt.colorbar(im)
                #c_min = mean
                #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
                #plt.show()
                #plt.close()
        if coincidence_check == 1:
            break


    # Creating a rough coordinate net (first version on the frame



    # frame_param = axis_fit_complete(coord_corresp_x, coord_corresp_y, mid_ra, mid_dec, coord_corresp_RA,
    #                           coord_corresp_DEC, initial_parameters)
    frame_param = min_param
    print("Sum of squared deviation for the minimal parameters: " + "{}".format(sum_for_min_param))
    print(frame_param)
    # ra_line_arr = arange(12, 13.2, 0.2)
    # ra_line_arr = arange(16, 25, 2)
    # print(axis_fit_use(18.93973, 85.094416, check_path, frame_param))
    # print(axis_fit_use(12.5, -55.9, check_path, frame_param))
    # ra_ref = 20.36050356461
    # dec_ref = 84.96567586615
    # centered_frame = np.copy(check_frame)
    # source_x, source_y = axis_fit_use(319.734086, 33.681089, check_path, frame_param)
    # circle_ref(source_x, source_y, 10, 5, 20, centered_frame)
    # plt.figure(2)
    # im = plt.imshow(centered_frame)
    # plt.colorbar(im)
    # c_min = mean
    # plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    # plt.show()
    # plt.close()

    # Adjusting a coordinate net using anchor stars
    # Picking anchors close to the center
    dev_x = max(absolute(coord_corresp_x - check_mid_x))
    dev_y = max(absolute(coord_corresp_y - check_mid_y))
    inter_anchor_data = {'RA': [], 'DEC': [], 'Gmag': [], 'x': [], 'y': []}  # Series of data for anchor stars
    center_anchor_stars = DataFrame(data=inter_anchor_data)
    for i in range(len(bright_anchor_stars['Gmag'])):
        x_interim, y_interim = axis_fit_use(bright_anchor_stars['RA_ICRS'].values[i],
                                            bright_anchor_stars['DE_ICRS'].values[i], mid_ra, mid_dec, frame_param)
        # circle_ref(x_interim, y_interim, 10, 5, 20, centered_frame)
        dev_x_inter = abs(x_interim - check_mid_x)
        dev_y_inter = abs(y_interim - check_mid_y)

        bright_closeness = 0  # Check whether a star is close to the brightest on the frame
        twin_check = 0  # Check whether a star has a twin already
        if dev_x_inter <= check_mid_x - 100 and dev_y_inter <= check_mid_y - 100:
            for j in range(len(super_bright_stars['RA_ICRS'])):
                #inter_bright_ang_separation = angular_separation(bright_anchor_stars['RA_ICRS'].values[i] / 57.3,
                #                                                 bright_anchor_stars['DE_ICRS'].values[i] / 57.3,
                #                                                 super_bright_stars['RA_ICRS'].values[j] / 57.3,
                #                                                 super_bright_stars['DE_ICRS'].values[j] / 57.3)
                #inter_bright_ang_separation = inter_bright_ang_separation * 57.3 * 60
                if check_frame[int(x_interim)][int(y_interim)] != 0:
                    if photometry_snr(x_interim, y_interim, 10, 4, 20, check_frame) < 3:
                        bright_closeness = 1
                        # x_test, y_test = x_interim, y_interim
                        # print(bright_anchor_stars['RA_ICRS'].values[i], bright_anchor_stars['DE_ICRS'].values[i])
                else:
                    bright_closeness = 1
                #if inter_bright_ang_separation < 1:
                #    bright_closeness = 1
                #    break
            if bright_closeness == 0:
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
                        if len(center_anchor_stars['RA']) >=30:
                            break
    #circle_check(center_anchor_stars['x'], center_anchor_stars['y'], 10, 5, 25, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    # Adjusting coordinates of anchor stars
    for i in range(len(center_anchor_stars['Gmag'])):
        center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i] = obj_centroid(
            center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i], 5, 5, 20, check_frame)

    # Estimating parameters of a coordinate net once more (more precise, hopefully)
    print('Number of anchor stars: ' + "{}".format(len(center_anchor_stars['Gmag'])))

    # every_frame_param[frame_num] = axis_fit_complete(center_anchor_stars['x'], center_anchor_stars['y'], mid_ra, mid_dec,
    #                           center_anchor_stars['RA'], center_anchor_stars['DEC'], frame_param)
    # frame_center_x = frame_header['CRPIX2'] - frame_header['YSTART']
    # frame_center_y = frame_header['CRPIX1'] - frame_header['XSTART']
    frame_center_x = frame_header['CRPIX2']
    frame_center_y = frame_header['CRPIX1']
    # every_frame_param[frame_num] = axis_fit_all_and_more(center_anchor_stars['x'], center_anchor_stars['y'],
    #                                                 center_anchor_stars['RA'], center_anchor_stars['DEC'], mid_ra,
    #                                                 mid_dec, frame_center_x, frame_center_y)

    every_frame_param[frame_num] = axis_fit_complete(center_anchor_stars['x'], center_anchor_stars['y'], mid_ra,
                                                     mid_dec,
                                                     center_anchor_stars['RA'], center_anchor_stars['DEC'], frame_param)
    # test_ra, test_dec = 319.90876155880, 33.73959957568
    # test_x, test_y = axis_fit_use_all(mid_ra, mid_dec, test_ra, test_dec, every_frame_param[frame_num])
    # test_x_1, test_y_1 = axis_fit_use(test_ra, test_dec, mid_ra, mid_dec, every_frame_param[frame_num])
    # x_array = np.arange(-10, 10, 0.1)
    # y_array = np.arange(-10, 10, 0.1)
    # z_array = np.zeros((len(x_array), len(y_array)))
    # for i in range(len(x_array)):
    #   for j in range(len(y_array)):
    #       inter_undistorted_x = undistorted_x + x_array[i]
    #       inter_undistorted_y = undistorted_y + y_array[j]
    #       z_array[i][j], inter_parameters = axis_fit_distortion(distorted_x, distorted_y, inter_undistorted_x, inter_undistorted_y, input_parameters)
    # = distortion_min_sum(distorted_x, distorted_y, inter_undistorted_x, inter_undistorted_y, inter_parameters)
    # every_frame_param[frame_num] = axis_fit_new(center_anchor_stars['x'], center_anchor_stars['y'], check_path,
    #                                                 center_anchor_stars['RA'], center_anchor_stars['DEC'], frame_param)

    print(every_frame_param[frame_num])
    # Limiting catalog of bright stars (initial anchors) to those that are inside the frame
    dict_stars = []
    array_columns = bright_stars.columns.values.tolist()
    if frame_num == 0:
        for bright_ind in range(len(bright_stars['Gmag'])):
            bright_catalog_stars_x, bright_catalog_stars_y = axis_fit_use_all(bright_stars['RA_ICRS'].values[bright_ind],
                                                                        bright_stars['DE_ICRS'].values[bright_ind], mid_ra, mid_dec,
                                                                        every_frame_param[frame_num])
            if abs(bright_catalog_stars_x - check_mid_x) <= check_mid_x - 10 and abs(bright_catalog_stars_y - check_mid_y) <= check_mid_y - 10:
                dict1 = dict((col_name, bright_stars[col_name].values[bright_ind]) for col_name in array_columns)
                dict_stars.append(dict1)
            circle_ref(bright_catalog_stars_x, bright_catalog_stars_y, 10, 5, 25, centered_frame)
        bright_stars = pd.DataFrame(dict_stars, columns=array_columns)
        #print(bright_stars)
        bright_quad_comb, bright_stars_RA, bright_stars_DEC = quad_catalog(bright_stars)
        #print(bright_quad_comb)
        #print(bright_quad_comb.shape)
    #else:
        #for bright_ind in range(len(bright_stars['Gmag'])):
        #    bright_catalog_stars_x, bright_catalog_stars_y = axis_fit_use_all(bright_stars['RA_ICRS'].values[bright_ind],
        #                                                                bright_stars['DE_ICRS'].values[bright_ind], mid_ra, mid_dec,
        #                                                                every_frame_param[frame_num])
        #    circle_ref(bright_catalog_stars_x, bright_catalog_stars_y, 10, 5, 25, centered_frame)
    # print('aaa')
    #circle_check(center_anchor_stars['x'], center_anchor_stars['y'], 10, 5, 25, centered_frame)
    # circle_ref(test_x, test_y, 10, 5, 25, centered_frame)
    # circle_ref(test_x_1, test_y_1, 10, 5, 25, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    #print(center_anchor_stars['x'])
    #print(center_anchor_stars['y'])
    #for i in range(len(center_anchor_stars['x'])):
    #    # print(i)
    #    center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i] = axis_fit_use_all(center_anchor_stars['RA'].values[i],
    #                                                                center_anchor_stars['DEC'].values[i], mid_ra, mid_dec,
    #                                                                every_frame_param[frame_num])
    #print(center_anchor_stars['x'])
    #print(center_anchor_stars['y'])
    #i = 0
    #center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i] = axis_fit_use_all(
    #    center_anchor_stars['RA'].values[i],
    #    center_anchor_stars['DEC'].values[i], mid_ra, mid_dec,
    #    every_frame_param[frame_num])
    #print(center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i])
    #print(center_anchor_stars['x'])
    #print(center_anchor_stars['y'])
    #centered_frame = copy(check_frame)
    # x_test, y_test = axis_fit_use(tester_ra, tester_dec, check_path, frame_param)
    # circle_ref(x_test, y_test, 10, 5, 25, centered_frame)
    # print(center_anchor_stars['RA'], center_anchor_stars['DEC'])
    #circle_check(center_anchor_stars['x'], center_anchor_stars['y'], 10, 5, 25, centered_frame)
    # circle_ref(test_x, test_y, 10, 5, 25, centered_frame)
    # circle_ref(test_x_1, test_y_1, 10, 5, 25, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    # Correcting array of check stars
    dict_check = []
    check_columns = check_stars.columns.values.tolist()
    check_len = len(check_stars['Gmag'])
    for i in range(check_len):
        # print(i)
        check_x[frame_num][i], check_y[frame_num][i] = axis_fit_use_all(check_stars['RA_ICRS'].values[i],
                                                                    check_stars['DE_ICRS'].values[i], mid_ra, mid_dec,
                                                                    every_frame_param[frame_num])
        #print(check_x[frame_num][i], check_y[frame_num][i])
        if check_frame[int(check_x[frame_num][i])][int(check_y[frame_num][i])] != 0 and abs(
                check_x[frame_num][i] - check_mid_x) <= check_mid_x - 50 and abs(
                check_y[frame_num][i] - check_mid_y) <= check_mid_y - 50 and check_frame[int(check_x[frame_num][i])][
            int(check_y[frame_num][i])] != 0:
            dict1 = dict((col_name, check_stars[col_name].values[i]) for col_name in check_columns)
            dict_check.append(dict1)
    check_stars = pd.DataFrame(dict_check, columns=check_columns)
    #print(check_x[frame_num])
    centered_frame = copy(check_frame)
    #circle_check(check_x[frame_num], check_y[frame_num], 2, check_gap, check_annul, centered_frame)
    # circle_ref(ref_x[frame_num], ref_y[frame_num], ap_radius, ref_gap, ref_annul, centered_frame)
    # circle_ref(source_x[frame_num], source_y[frame_num], ap_radius * 2, ref_gap, ref_annul, centered_frame)
    #plt.figure()
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    dict_ref = []
    for i in range(len(ref_star)):
        ref_x[frame_num], ref_y[frame_num] = axis_fit_use(ref_star['RA_ICRS'].values[i], ref_star['DE_ICRS'].values[i],
                                                          mid_ra, mid_dec, every_frame_param[frame_num])
        if check_frame[int(ref_x[frame_num])][int(ref_y[frame_num])] != 0:
            dict2 = dict((col_name, ref_star[col_name].values[i]) for col_name in check_columns)
            dict_ref.append(dict2)
    ref_star = pd.DataFrame(dict_ref, columns=check_columns)
    #mean = astropy.stats.sigma_clipped_stats(check_frame, sigma=1.5)[0]
    #centered_frame = copy(check_frame)
    #circle_check(check_x[frame_num], check_y[frame_num], 3, 5, 25, centered_frame)
    #plt.figure(2)
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
# Picking the brightest ref
# dict_check = []
# if len(ref_star['Gmag']) != 0:
#    for i in range(len(ref_star)):
#        dict1 = dict((col_name, ref_star[col_name].values[i]) for col_name in check_columns)
#        dict_check.append(dict1)
# add_check_stars = pd.DataFrame(dict_check)
# check_stars = pd.concat([check_stars, add_check_stars], axis=0)
dict_ref = []
# print(ref_star)
dict2 = dict((col_name, ref_star[col_name].values[0]) for col_name in check_columns)
ref_star = ref_star.drop([0])
# print(ref_star)
check_stars = pd.concat([check_stars, ref_star])
#print(check_stars)
check_stars.reset_index(drop = True)
if len(check_stars['Gmag']) >= 30:
    check_stars = check_stars.iloc[:30]
#print(check_stars)
dict_ref.append(dict2)
ref_star = pd.DataFrame(dict_ref, columns=check_columns)
print("Coordinates of a reference star: RA - " + "{}".format(
    ref_star['RA_ICRS'].values[0]) + ' deg, DE - ' + "{}".format(ref_star['DE_ICRS'].values[0]) + ' deg.')
# Performing photometry
# Readjusting check star arrays
check_len = len(check_stars['Gmag'])
check_photo = zeros((photo_len, check_len))
check_photo_err = zeros((photo_len, check_len))
check_x = zeros((photo_len, check_len))
check_y = zeros((photo_len, check_len))
check_true_mag = zeros((photo_len, check_len))
check_color = check_stars['BP-RP']
# mag_ref = zeros
# print(check_stars)
for frame_num in range(photo_len):
    # Parameters of a frame
    test_frame = open("{}".format(filepath + photo_frame_path[frame_num]))
    check_frame = getdata("{}".format(filepath + photo_frame_path[frame_num]))
    test_frame_header = test_frame[0].header
    frame_exp[frame_num] = test_frame_header['EXPTIME']
    #frame_datetime[frame_num] = test_frame_header['JD'] + frame_exp[frame_num] / 2 / 60 / 60 / 24
    frame_inter_datetime = test_frame_header['DATE-OBS']
    inter_datetime = dateutil.parser.parse(frame_inter_datetime)
    inter_JD_date = inter_datetime.toordinal()
    inter_JD_date_datetime = datetime.datetime.fromordinal(inter_JD_date)
    inter_JD_date += 1721424.5
    inter_JD_date_full = inter_JD_date + (inter_datetime - inter_JD_date_datetime).total_seconds() / 86400
    frame_datetime[frame_num] = inter_JD_date_full + frame_exp[frame_num] / 2 / 60 / 60 / 24
    frame_tube = np.append(frame_tube, test_frame_header['AXIS'])
    frame_filter = np.append(frame_filter, test_frame_header['FILTER'])
    photo_frame = getdata("{}".format(filepath + photo_frame_path[frame_num]))
    mid_ra = float(test_frame_header["CRVAL1"])
    mid_dec = float(test_frame_header["CRVAL2"])
    # Coordinates of source, ref, checks
    # print('bbb')

    # Making a mask array
    check_masks_arr = zeros((check_len, 20, 20))
    check_mask = zeros((20, 20))
    # circle_ref(source_x[frame_num], source_y[frame_num], 10, 5, 25, centered_frame)
    # plt.figure(2)
    # im = plt.imshow(centered_frame)
    # plt.colorbar(im)
    # c_min = mean
    # plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    # plt.show()
    # plt.close()
    ref_x[frame_num], ref_y[frame_num] = axis_fit_use(ref_star['RA_ICRS'].values[0], ref_star['DE_ICRS'].values[0],
                                                      mid_ra, mid_dec, every_frame_param[frame_num])
    ref_x[frame_num], ref_y[frame_num] = obj_centroid(ref_x[frame_num], ref_y[frame_num], 3, 5, 15,
                                                      check_frame)
    for i in range(check_len):
        check_x[frame_num][i], check_y[frame_num][i] = axis_fit_use(check_stars['RA_ICRS'].values[i],
                                                                    check_stars['DE_ICRS'].values[i], mid_ra, mid_dec,
                                                                    every_frame_param[frame_num])
        check_x[frame_num][i], check_y[frame_num][i], check_masks_arr[i] = check_mask_center(check_x[frame_num][i],
                                                                                             check_y[frame_num][i], 3,
                                                                                             5,
                                                                                             15,
                                                                                             check_frame)
        check_mask += check_masks_arr[i]
    check_mask = check_mask / max(check_mask)
    check_peak_x, check_peak_y = obj_peak(10, 10, 10, check_mask)
    # plt.figure(2)
    # im = plt.imshow(check_mask)
    # plt.colorbar(im)
    # c_min = mean
    # plt.clim([0, 1])
    # plt.show()
    # plt.close()
    # Erasing the bright star near GRB
    near_star_x, near_star_y = obj_peak(512, 512, 15, check_frame)
    # print(near_star_x, near_star_y)
    near_star_peak = check_frame[near_star_x][near_star_y] - back_level(512, 512, 10, 20, check_frame)
    d_x, d_y = near_star_x - check_peak_x, near_star_y - check_peak_y
    #for x_ind in range(d_x, d_x + 20):
    #    for y_ind in range(d_y, d_y + 20):
    #        check_frame[x_ind][y_ind] -= near_star_peak * check_mask[x_ind - d_x][y_ind - d_y]
    # plt.figure(2)
    # im = plt.imshow(check_frame)
    # plt.colorbar(im)
    # mean = astropy.stats.sigma_clipped_stats(check_frame, sigma=1.5)[0]
    # c_min = mean
    # plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    # plt.show()
    # plt.close()
    # The photometry
    # print(check_x[frame_num])
    # print(check_y[frame_num])
    # print('ddd')
    source_x[frame_num], source_y[frame_num] = axis_fit_use(source_ra, source_dec, mid_ra, mid_dec,
                                                            every_frame_param[frame_num])
    #source_x[frame_num], source_y[frame_num] = 512, 512
    source_x[frame_num], source_y[frame_num] = obj_centroid(source_x[frame_num], source_y[frame_num], source_centroid_radius, 5, 15,
                                                            check_frame)
    # source_x[frame_num], source_y[frame_num] = obj_peak(source_x[frame_num], source_y[frame_num], 2, check_frame)
    if source_aperture_radius == 0:
        source_photo[frame_num], source_err[frame_num], ap_radius = adaptive_photometry_source(source_x[frame_num],
                                                                                               source_y[frame_num],
                                                                                               check_frame)
        if ap_radius > 6:
            ap_radius = 6
    else:
        ap_radius = source_aperture_radius
        source_photo[frame_num], source_err[frame_num] = photometry_source(source_x[frame_num], source_y[frame_num],
                                                                           ap_radius, 4, 10, photo_frame)

    ref_photo[frame_num] = photometry_ref(ref_x[frame_num], ref_y[frame_num], ap_radius, ref_gap, ref_annul,
                                          photo_frame)
    check_photo[frame_num], check_photo_err[frame_num] = photometry_check(check_x[frame_num], check_y[frame_num],
                                                                          ap_radius, check_gap,
                                                                          check_annul, photo_frame)
    print('Radius of adaptive aperture: ' + "{}".format(ap_radius))
    # if frame_num == 0:
    #    print(ap_radius)

    if ap_radius <= 2:
        mag_ref[frame_num] = ref_star['Gmag_2'].values[0]
        check_true_mag[frame_num] = check_stars['Gmag_2']
    elif ap_radius > 2 and ap_radius <= 4:
        mag_ref[frame_num] = ref_star['Gmag_4'].values[0]
        check_true_mag[frame_num] = check_stars['Gmag_4']
    elif ap_radius > 4 and ap_radius <= 6:
        mag_ref[frame_num] = ref_star['Gmag_6'].values[0]
        check_true_mag[frame_num] = check_stars['Gmag_6']
    elif ap_radius > 6 and ap_radius <= 8:
        mag_ref[frame_num] = ref_star['Gmag_8'].values[0]
        check_true_mag[frame_num] = check_stars['Gmag_8']

    #centered_frame = np.copy(check_frame)
    #circle_check(check_x[frame_num], check_y[frame_num], ap_radius, check_gap, check_annul, centered_frame)
    # circle_ref(ref_x[frame_num], ref_y[frame_num], ap_radius, ref_gap, ref_annul, centered_frame)
    # circle_ref(source_x[frame_num], source_y[frame_num], ap_radius * 2, ref_gap, ref_annul, centered_frame)
    #plt.figure()
    #im = plt.imshow(centered_frame)
    #plt.colorbar(im)
    #c_min = mean
    #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
# print(check_true_mag)
# print(mag_ref)
# Determining magnitude
source_mag, check_mag = magnitude(mag_ref[frame_num], ref_photo, source_photo, check_photo)
# print(check_true_mag)
# Deleting apparently overly variable stars
new_check_photo = copy(check_photo)
new_check_x = copy(check_x)
new_check_y = copy(check_y)
new_check_mag = copy(check_mag)
new_check_true_mag = copy(check_true_mag)
new_check_color = copy(check_color)
new_check_photo_err = copy(check_photo_err)
# Checking the first frame
del_ind = 0
for j_ind in range(len(check_photo[0])):
    # print(check_mag[0][j_ind])
    # print(check_true_mag[j_ind])
    if abs(check_mag[0][j_ind] - check_true_mag[0][j_ind]) >= 0.5:
        # new_check_photo = append(new_check_photo, check_photo[:, j_ind], axis = 1)
        # new_check_x = append(new_check_x, check_x[:, j_ind], axis = 1)
        # new_check_y = append(new_check_y, check_y[:, j_ind], axis = 1)
        # new_check_mag = append(new_check_mag, check_mag[:, j_ind], axis = 1)
        # new_check_true_mag = append(new_check_true_mag, check_true_mag[:, j_ind], axis = 1)
        # new_check_color = append(new_check_color, check_color[j_ind])
        # new_check_photo_err = append(new_check_photo_err, check_photo_err[:, j_ind], axis = 1)
        new_check_photo = delete(new_check_photo, j_ind - del_ind, axis=1)
        new_check_x = delete(new_check_x, j_ind - del_ind, axis=1)
        new_check_y = delete(new_check_y, j_ind - del_ind, axis=1)
        new_check_mag = delete(new_check_mag, j_ind - del_ind, axis=1)
        new_check_true_mag = delete(new_check_true_mag, j_ind - del_ind, axis=1)
        new_check_color = delete(new_check_color, j_ind - del_ind)
        new_check_photo_err = delete(new_check_photo_err, j_ind - del_ind, axis=1)
        del_ind += 1
check_photo = copy(new_check_photo)
check_photo_err = copy(new_check_photo_err)
check_x = copy(new_check_x)
check_y = copy(new_check_y)
check_mag = copy(new_check_mag)
check_true_mag = copy(new_check_true_mag)
check_color = copy(new_check_color)
# print(check_color)
# print(check_mag - check_true_mag)
# print(del_ind)
check_len = len(check_color)
# print(check_len)
frame_shape = check_frame.shape
# print(check_shape)
# Determining magnitude errors for source and checks
source_mag_var_err, source_mag_signal_err, check_mag_var_err, check_mag_signal_err = magnitude_err(source_photo,
                                                                                                   source_err,
                                                                                                   check_mag,
                                                                                                   check_photo,
                                                                                                   check_photo_err)

# Finding fit parameters for all the frames
# x_par = zeros((photo_len, 18))
# x_par_err = zeros((photo_len, 18))
source_mag_corr = zeros_like(source_mag)
check_mag_corr = zeros_like(check_mag)
# Finding parameter for correcting for frame inadequacies and atmospheric effects(Forbes effect)
# Maximum polynom number
# print(len(check_photo[0]))
print('Number of checks: ' + "{}".format(len(check_color)))
coord_poly_num = int((-3 + sqrt(1 + 8 * (len(check_mag[0]) - color_param_num))) / 2 - 2)
if coord_poly_num >3:
    coord_poly_num = 3
# print(coord_poly_num)
# print(check_mag[0] - check_true_mag[0])
x_par, x_par_err = full_var_correct(check_mag, check_mag_signal_err,
                                    check_true_mag, check_color, check_x, check_y, frame_shape, coord_poly_num)
# print(x_par)
# Correcting checks and source
source_mag_corr = source_mag_correct(source_mag, source_color, source_x, source_y, x_par, frame_shape, coord_poly_num)
check_mag_corr = check_mag_correct(check_mag, check_color, check_x, check_y, x_par, frame_shape, coord_poly_num)

source_mag_corr_var_err, source_mag_signal_err, check_mag_corr_var_err, check_mag_signal_err = magnitude_err(
    source_photo, source_err, check_mag_corr, check_photo, check_photo_err)

# Estimating source mag error from check signal noise
source_mag_check_signal_err = np.zeros_like(source_mag)
signal_err_const = 0
for k in range(photo_len):
    for i in range(check_len):
        signal_err_const += check_mag_signal_err[k][i] ** 2
    source_mag_check_signal_err[k] = math.sqrt(signal_err_const / (check_len - 1))

# Determining deviation from true mag for every check
dev_data = {"Frame": photo_frame_path}
# print(dev_data)
deviation_dataframe = pd.DataFrame(data=dev_data, index=frame_datetime)
deviation_dataframe.index.name = 'Julian Days, middle of exposure'
# print(deviation_dataframe)
# print(deviation_dataframe)
# print(check_mag.shape[1])
for k in range(check_mag.shape[1]):
    check_mag_deviation_data = np.zeros_like(source_mag)
    check_mag_signal_err_data = np.zeros_like(source_mag)
    check_mag_color = np.zeros_like(source_mag)
    check_mag_true_mag = np.zeros_like(source_mag)
    for j in range(check_mag.shape[0]):
        check_mag_deviation_data[j] = check_mag[j][k] - check_true_mag[j][k]
        check_mag_signal_err_data[j] = check_mag_signal_err[j][k]
        check_mag_color[j] = check_color[k]
        check_mag_true_mag[j] = check_true_mag[j][k]

    df = {'Check star ' + "{}".format(k + 1) + ' deviation, mag': check_mag_deviation_data,
          'Check star ' + "{}".format(k + 1) + 'mag error': check_mag_signal_err_data,
          'Check star ' + "{}".format(k + 1) + ' GAIA G mag': check_mag_true_mag,
          'Check star ' + "{}".format(k + 1) + ' BP-RP color': check_mag_color}
    # print(df)
    df2 = pd.DataFrame(data=df, index=frame_datetime)
    # print(df2)
    deviation_dataframe = deviation_dataframe.join(df2)
    # print(deviation_dataframe)
    # print(deviation_dataframe)
    # deviation_dataframe.drop_duplicates(subset='Frame', keep='first')

# Determining deviation from true mag for every check
check_dev_first_frame = np.zeros_like(check_color)
check_signal_err_first_frame = np.zeros_like(check_color)

for k in range(check_mag.shape[1]):
    j = 0
    check_dev_first_frame[k] = check_mag[j][k] - check_true_mag[j][k]
    check_signal_err_first_frame[k] = check_mag_signal_err[j][k]

check_color_first_frame = np.copy(check_color)
# check_true_mag_first_frame = np.copy(check_col)
df1 = {'Check stars BP-RP color': check_color_first_frame,
       'Check stars GAIA G mag': check_true_mag[0],
       'Check stars deviation on first frame, mag': check_dev_first_frame,
       'Check star signal error, mag': check_signal_err_first_frame,
       'Check star x coord - x mid': check_x[0] - check_mid_x,
       'Check star y coord - y mid': check_y[0] - check_mid_y,
       'Check star x coord': check_x[0],
       'Check star y coord': check_y[0]}
first_frame_dev_dataframe = pd.DataFrame(data=df1)

# Fit data in a

fit = {"Frame": photo_frame_path}
fit_data = pd.DataFrame(data=fit, index=frame_datetime)
coord_par_num = (x_par.shape[1] - 1) / 2
for m in range(x_par.shape[1]):
    export_x_par = np.zeros(x_par.shape[0])
    export_x_par_err = np.zeros(x_par.shape[0])
    for k in range(x_par.shape[0]):
        export_x_par[k] = x_par[k][m]
        export_x_par_err[k] = x_par_err[k][m]
    if m <= coord_par_num:
        export_column = "X coord parameter " + "{}".format(m)
        export_err_column = "X coord parameter " + "{}".format(m) + "error"
    else:
        export_column = "Y coord parameter " + "{}".format(m - coord_par_num)
        export_err_column = "Y coord parameter " + "{}".format(m - coord_par_num) + "error"
    df = {export_column: export_x_par,
          export_err_column: export_x_par_err}
    df_export = pd.DataFrame(data=df, index=frame_datetime)
    fit_data = fit_data.join(df_export)

# Determining deviation of theoretical prediction from the data and comparing with empirical (squared sum of errors)
diff_empiric_dev = np.zeros(check_mag.shape[0])  # Difference between deviation of data and theoretical deviation
diff_check_err = np.zeros(check_mag.shape[0])
diff_sigma = np.zeros(check_mag.shape[0])
for k in range(check_mag.shape[0]):
    for m in range(check_mag.shape[1]):
        diff_check_dev = check_mag[k][m] - check_true_mag[k][m]
        diff_check_theor_dev = x_par[k][0]
        power_check = 1
        for i in range(1, color_param_num + 1):
            diff_check_theor_dev += check_color[m] ** power_check * x_par[k][i]
            power_check += 1
        power_check = 1
        # for i in range(color_param_num + 1, x_par.shape[1], 2):
        #    diff_check_theor_dev += x_par[k][i] * (check_x[k][m] - check_mid_x) ** power_check
        #    diff_check_theor_dev += x_par[k][i + 1] * (check_y[k][m] - check_mid_y) ** power_check
        #    power_check += 1
        # for i in range(color_param_num + 1, x_par.shape[1]):
        #    diff_check_theor_dev += x_par[k][i] * (check_x[k][m] - check_mid_x + check_y[k][m] - check_mid_y) ** power_check
        #    power_check += 1
        place_index = color_param_num + 1
        for poly_num in range(1, coord_poly_num + 1):
            for coord_ind in range(poly_num + 1):
                diff_check_theor_dev += x_par[k][place_index] * (check_x[k][i] - check_mid_x) ** (
                            poly_num - coord_ind) * (
                                                check_y[k][i] - check_mid_y) ** (coord_ind)
                place_index += 1
        diff_check_err[k] += check_mag_signal_err[k][m] ** 2
        diff_empiric_dev[k] += (diff_check_dev - diff_check_theor_dev) ** 2
    diff_check_err[k] = math.sqrt(diff_check_err[k])
    diff_empiric_dev[k] = math.sqrt((diff_empiric_dev[k]))
    diff_sigma[k] = diff_empiric_dev[k] / diff_check_err[k]
diff_predf = {"Root of a squared sum of check signal errors": diff_check_err,
              "Root of a squared sum of differences between deviation of check mag from the GAIA mag and theoretical deviation": diff_empiric_dev,
              "Ratio of a squared sum of check signal errors and differences": diff_sigma}
diff_df = pd.DataFrame(data=diff_predf, index=frame_datetime)

fit_data = fit_data.join(diff_df)
fit_data.index.name = 'Julian Days, middle of exposure'

# Exporting results of the photometry into excel in the same folder
data = {"Frame": photo_frame_path, 'Tube': frame_tube, 'Filter': frame_filter,
        'Frame exposure time, seconds': frame_exp, "Source mag (corrected)": source_mag_corr,
        'Source mag error (variability, corrected)': source_mag_corr_var_err,
        "Source mag error (signal noise)": source_mag_signal_err,
        "Source mag error (signal noise of checks)": source_mag_check_signal_err,
        "Source x coord": source_x, "Source y coord": source_y}

full_data = {"Frame": photo_frame_path, 'Tube': frame_tube, 'Filter': frame_filter,
             'Frame exposure time, seconds': frame_exp, 'Source mag': source_mag,
             "Source mag (corrected)": source_mag_corr,
             'Source mag error (variability, not corrected)': source_mag_var_err,
             'Source mag error (variability, corrected)': source_mag_corr_var_err,
             "Source mag error (signal noise)": source_mag_signal_err,
             "Source mag error (signal noise of checks)": source_mag_check_signal_err,
             "Source x coord": source_x, "Source y coord": source_y}
photo_data = pd.DataFrame(data=data, index=frame_datetime)
full_photo_data = pd.DataFrame(data=full_data, index=frame_datetime)
photo_data.index.name = 'Julian Days, middle of exposure'
full_photo_data.index.name = 'Julian Days, middle of exposure'

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

# Adding checks

test_frame = open("{}".format(filepath + filename))
hdr = test_frame[0].header
full_photo_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(
    hdr['FILTER']) + ' full photometry data (source and checks).xlsx'
photo_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' photometry data (source).xlsx'
fit_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' fit data.xlsx'
first_frame_dev_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' first frame deviation data.xlsx'
deviation_file = "{}".format(hdr['AXIS']) + ' ' + "{}".format(hdr['FILTER']) + ' deviation data.xlsx'
photo_data.to_excel("{}".format(photometry_path + photo_file))
full_photo_data.to_excel("{}".format(photometry_path + full_photo_file))
fit_data.to_excel("{}".format(photometry_path + fit_file))
deviation_dataframe.to_excel("{}".format(photometry_path + deviation_file))
first_frame_dev_dataframe.to_excel("{}".format(photometry_path + first_frame_dev_file))

print("Exporting over")
print('Program took ' + "{}".format((time() - start_time) / 60) + ' minutes')
