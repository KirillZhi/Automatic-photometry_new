# import numpy as np
from numpy import zeros, empty, append, arange, max, absolute
# import pandas as pd
from pandas import DataFrame
import photutils
from itertools import combinations
# photutils
from frame_functions import *
from coord_functions import *
# from astropy.io import fits
from astropy.io.fits import getdata, open
from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename
from time import sleep
# import os
from os import listdir
from os.path import isfile, join
# import re
from re import search
from matplotlib import pyplot as plt
from math import sqrt

# Filepath for fits frames
readout_noise = 9
root = Tk()
root.withdraw()
filepath = askdirectory(initialdir="C:/Users/kirio/Desktop/astronomy/GRB/")
photometry_path = askdirectory(initialdir="C:/Users/kirio/Desktop/astronomy/GRB/")
catalog_path = askopenfilename(initialdir="C:/Users/kirio/Desktop/astronomy/GRB/")

filepath = filepath + "/"
photometry_path = photometry_path + "/"
root.update()
root.destroy()

# Making quads from the catalogue
bright_stars, interm_bright_stars, bright_anchor_stars, catalog_stars = gaia_bright(catalog_path)
bright_quad_comb, bright_stars_RA, bright_stars_DEC = quad_catalog(bright_stars)
print('Catalog analyzed, beginning operations with the frames')
# Coordinates of a test star
tester_ra = 12.82938663368
tester_dec = -55.73110941053
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
#print(photo_frame_path)

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

    # Minimum value of the frame (positive)
    test_ground = np.where(check_frame > 0, check_frame, np.max(check_frame))
    c_min = np.percentile(test_ground, 5)

    # Looking for sources in the frame
    mean = astropy.stats.sigma_clipped_stats(check_frame, sigma=1.5)[0]
    std = math.sqrt(mean + readout_noise ** 2)
    daofind = photutils.detection.DAOStarFinder(fwhm=3.0, threshold=15. * std, sigma_radius=2, exclude_border=True)
    sources = daofind(check_frame - mean)

    # Picking 6 brightest sources
    sorted_sources = sources.group_by('peak')

    # Precise photometry of the sources
    for i in range(len(sorted_sources['flux']) - 1, len(sorted_sources['flux']) - 30, -1):
        sorted_sources['flux'][i] = photometry_ref(sorted_sources['ycentroid'][i], sorted_sources['xcentroid'][i], 5, 10,
                                                   20, check_frame)

    sorted_sources = sorted_sources.group_by('flux')
    # Array of the 6 brightest stars on the frame
    bright_source_x = []
    bright_source_y = []
    for i in range(len(sorted_sources['flux']) - 1, len(sorted_sources['flux']) - 20, -1):
        if len(bright_source_x) < 6:
            discard_check = 0
            if len(bright_source_x) != 0:
                for j in range(len(bright_source_x)):
                    if abs(sorted_sources['ycentroid'][i] - bright_source_x[j]) < 10 and abs(
                            sorted_sources['xcentroid'][i] - bright_source_y[j]) < 10:
                        discard_check = 1
                        break
                if discard_check == 1:
                    continue
                x_interim, y_interim = obj_centroid(sorted_sources['ycentroid'][i], sorted_sources['xcentroid'][i], 1, 10,
                                                    20, check_frame)
                bright_source_x = np.append(bright_source_x, x_interim)
                bright_source_y = np.append(bright_source_y, y_interim)
            else:
                x_interim, y_interim = obj_centroid(sorted_sources['ycentroid'][i], sorted_sources['xcentroid'][i], 5, 10,
                                                    20, check_frame)
                bright_source_x = np.append(bright_source_x, x_interim)
                bright_source_y = np.append(bright_source_y, y_interim)
        else:
            break
    # Creating combinations of stars on a frame
    index_arr = np.arange(0, len(bright_source_x))
    index_comb = list(combinations(index_arr, 4))
    comb_num = len(index_comb)
    comb_x = np.zeros((comb_num, 5))
    comb_y = np.zeros((comb_num, 5))
    for i in range(comb_num):
        comb_x[i][0] = i
        comb_y[i][0] = i
        for j in range(4):
            comb_x[i][j + 1] = bright_source_x[index_comb[i][j]]
            comb_y[i][j + 1] = bright_source_y[index_comb[i][j]]

    for ind in range(comb_x.shape[0]):
        bright_frame_quad = quad_frame(comb_x[ind], comb_y[ind])
        for j in range(bright_quad_comb.shape[0]):
            sum_check = 0                       # Check of a number of coincidences
            check_type = 0                      # Check of a type of coincidence
            for i in range(1, 8, 2):
                test_value = bright_frame_quad[i]                  # Value of the frame quad that we test
                necessary_precis = bright_frame_quad[i + 1] * 2
                for k in range(1, 5):
                    # Value for checking
                    check_value = bright_quad_comb[j][k]
                    symm_value = 1 - check_value
                    if check_type == 0:
                        if abs(test_value - check_value) <= necessary_precis:
                            check_type = 1
                            sum_check += 1
                            break
                        if abs(test_value - symm_value) <= necessary_precis:
                            check_type = 2
                            sum_check += 1
                            break
                    if check_type == 1:
                        if abs(test_value - check_value) <= necessary_precis:
                            sum_check += 1
                            break
                    if check_type == 2:
                        if abs(test_value - symm_value) <= necessary_precis:
                            sum_check += 1
                            break
            if sum_check == 4:
                right_catalog_quad = bright_quad_comb[j]
                break
        if sum_check == 4:
            break

    quad_ind = int(right_catalog_quad[0])

    # Identify stars with the quad between catalog and frame
    index_correlation_catal_frame = np.zeros(4)  # Index of a cell corresponds to an index of a star in the frame combination - 1,
                                                 # number in a cell corresponds to an index of a star in the catalog combination

    check_type = 0
    correlation_index = 0
    symmetry_index = 0  # Checks whether farthest stars are mirrored, 0 if not, 1 if yes
    reversed_index = 0  # Checks whether inner stars are reversed in order, 0 if not, 1 if yes.
                        # Meaning if the first inner star in the catalog is the first in the frame
                        # index is 0. 1 is not
    for i in range(1):
        test_value = bright_quad_comb[quad_ind][i + 1]
        for j in range(1, 8, 2):
            check_value = bright_frame_quad[j]
            symm_value = 1 - bright_frame_quad[j]
            necessary_precis = bright_frame_quad[j + 1]
            if check_type == 0:
                if abs(test_value - check_value) <= necessary_precis:
                    check_type = 1
                    if j <= 3:
                        reversed_index = 0
                    else:
                        reversed_index = 1
                    break
                if abs(test_value - symm_value) <= necessary_precis:
                    check_type = 2
                    symmetry_index = 1
                    if j <= 3:
                        reversed_index = 0
                    else:
                        reversed_index = 1
                    break
            if check_type == 1:
                if abs(test_value - check_value) <= necessary_precis:
                    break
            if check_type == 2:
                if abs(test_value - symm_value) <= necessary_precis:
                    break
    # Placing indexes in the array
    if symmetry_index == 1:
        index_correlation_catal_frame[int(bright_frame_quad[9]) - 1] = int(bright_quad_comb[quad_ind][6])
        index_correlation_catal_frame[int(bright_frame_quad[10]) - 1] = int(bright_quad_comb[quad_ind][5])
    else:
        index_correlation_catal_frame[int(bright_frame_quad[9]) - 1] = int(bright_quad_comb[quad_ind][5])
        index_correlation_catal_frame[int(bright_frame_quad[10]) - 1] = int(bright_quad_comb[quad_ind][6])

    if reversed_index == 0:
        index_correlation_catal_frame[int(bright_frame_quad[11]) - 1] = int(bright_quad_comb[quad_ind][7])
        index_correlation_catal_frame[int(bright_frame_quad[12]) - 1] = int(bright_quad_comb[quad_ind][8])
    else:
        index_correlation_catal_frame[int(bright_frame_quad[11]) - 1] = int(bright_quad_comb[quad_ind][8])
        index_correlation_catal_frame[int(bright_frame_quad[12]) - 1] = int(bright_quad_comb[quad_ind][7])
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
    # Creating a rough coordinate net (first version on the frame
    x_vector_param, y_vector_param = coord_net_fit(coord_corresp_x, coord_corresp_y, check_path, coord_corresp_RA,
                                                   coord_corresp_DEC, 1)
    line_ind = 0
    centered_frame = np.copy(check_frame)
    #ra_line_arr = arange(12, 13.2, 0.2)
    ra_line_arr = arange(16, 25, 2)
    print(ra_line_arr)
    #if line_ind % 3 == 0:
        #for ra_line in ra_line_arr:
            #dec_line(x_vector_param, y_vector_param, -55.9, ra_line, centered_frame, check_path)
            #dec_line(x_vector_param, y_vector_param, 85, ra_line, centered_frame, check_path)
        # plt.ioff()
        #plt.figure()
        #im = plt.imshow(centered_frame)
        #plt.colorbar(im)
        #c_min = mean
        #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
        # plt.clim([50000, 75000])
        # plt.show(block=False)
        #plt.show()
    # Adjusting a coordinate net using anchor stars
    # Picking anchors close to the center
    dev_x = max(absolute(coord_corresp_x - check_mid_x))
    dev_y = max(absolute(coord_corresp_y - check_mid_y))
    line_ind = 0
    while dev_x < check_mid_x - 30 and dev_y < check_mid_y - 30:
        print(dev_x, dev_y)
        dev_max_x = 0
        dev_max_y = 0
        inter_anchor_data = {'RA': [], 'DEC': [], 'Gmag': [], 'x': [], 'y': []}         #Series of data for anchor stars
        center_anchor_stars = DataFrame(data=inter_anchor_data)
        for i in range(len(bright_anchor_stars['Gmag'])):
            #print(i)
            x_interim, y_interim = coord_net_use(x_vector_param, y_vector_param, check_path,
                                                 bright_anchor_stars['RA_ICRS'].values[i],
                                                 bright_anchor_stars['DE_ICRS'].values[i])
            dev_x_inter = abs(x_interim - check_mid_x)
            dev_y_inter = abs(y_interim - check_mid_y)

            bright_closeness = 0                #Check whether a star is close to the brightest on the frame
            twin_check = 0                      #Check whehter a star has a twin already
            if dev_x_inter <= dev_x and dev_y_inter <= dev_y and 1:
                for j in range(len(bright_stars['RA_ICRS'])):
                    #print(j)
                    inter_bright_ang_separation = angular_separation(bright_anchor_stars['RA_ICRS'].values[i] / 57.3,
                                                                     bright_anchor_stars['DE_ICRS'].values[i] / 57.3,
                                                                     bright_stars['RA_ICRS'].values[j] / 57.3,
                                                                     bright_stars['DE_ICRS'].values[j] / 57.3)
                    inter_bright_ang_separation = inter_bright_ang_separation * 57.3 * 60
                    if photometry_snr(x_interim, y_interim, 4, 4, 12, check_frame) < 15:
                        bright_closeness = 1
                    if inter_bright_ang_separation < 1.5:
                        bright_closeness = 1
                        break
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
                            if dev_x_inter > dev_max_x:
                                dev_max_x = dev_x_inter
                            if dev_y_inter > dev_max_y:
                                dev_max_y = dev_y_inter

        # Adjusting coordinates of anchor stars
        for i in range(len(center_anchor_stars['Gmag'])):
            center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i] = obj_centroid(
                center_anchor_stars['x'].values[i], center_anchor_stars['y'].values[i], 5, 5, 20, check_frame)
        # Estimating parameters of a coordinate net once more (more precise, hopefully)
        # Finding number of polynoms
        anchor_poly_num = int((len(center_anchor_stars['Gmag']) - 1) / 2) - 2
        print(center_anchor_stars)
        #anchor_poly_num = int((-3 + sqrt(1 + 8 *len(center_anchor_stars['Gmag'])))/2)
        print('Number of anchor stars: ' + "{}".format(len(center_anchor_stars['Gmag'])))
        print('Highest order of the polynom: ' + "{}".format(anchor_poly_num))
        print('Max x dev: ' + "{}".format(dev_max_x))
        print('Max y dev: ' + "{}".format(dev_max_y))
        #if anchor_poly_num <= 2:
        #   anchor_poly_num = 3
        #if anchor_poly_num > 4:
           #anchor_poly_num = 4
        x_vector_param, y_vector_param = coord_net_fit(center_anchor_stars['x'], center_anchor_stars['y'], check_path, center_anchor_stars['RA'], center_anchor_stars['DEC'], anchor_poly_num)
        centered_frame = np.copy(check_frame)
        print(center_anchor_stars)
        circle_check(center_anchor_stars['x'], center_anchor_stars['y'], 5, 10, 25, centered_frame)
        plt.figure()
        im = plt.imshow(centered_frame)
        plt.colorbar(im)
        c_min = mean
        plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
        # plt.clim([50000, 75000])
        # plt.show(block=False)
        plt.show()
        # plt.pause(15)
        plt.close()
        #ra_line_arr = arange(11, 13, 0.2)
        ra_line_arr = arange(16, 25, 2)
        print(ra_line_arr)
        #if line_ind % 3 == 0:
            #for ra_line in ra_line_arr:
                #dec_line(x_vector_param, y_vector_param, -55.9, ra_line, centered_frame, check_path)
                #dec_line(x_vector_param, y_vector_param, 85, ra_line, centered_frame, check_path)
            #plt.ioff()
            #plt.figure()
            #im = plt.imshow(centered_frame)
            #plt.colorbar(im)
            #c_min = mean
            #plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
            # plt.clim([50000, 75000])
            #plt.show(block=False)
            #plt.show()
        #plt.pause(5)
            #plt.close()
        dev_x += 40
        dev_y += 40
        line_ind += 1
    tester_ra = 12.98801957447
    tester_dec = -56.04946635794
    #tester_ra = 16.39263137587
    #tester_dec = 85.12594280207
    #tester_x, tester_y = coord_net_use(x_vector_param, y_vector_param, check_path, tester_ra, tester_dec)
    centered_frame = np.copy(check_frame)
    #circle_check(center_anchor_stars['x'], center_anchor_stars['y'], 5, 10, 25, centered_frame)
    #print(tester_x, tester_y)
    #circle_ref(tester_x, tester_y, 5, 10, 25, centered_frame)
    #plt.ion()
    plt.ioff()
    plt.figure()
    im = plt.imshow(centered_frame)
    plt.colorbar(im)
    c_min = mean
    plt.clim([c_min, c_min + math.sqrt(mean + readout_noise ** 2) * 8])
    # plt.clim([50000, 75000])
    #plt.show(block=False)
    plt.show()
    #plt.pause(15)
    plt.close()
