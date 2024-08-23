# Functions for automatic photometry
# Circles and ellipses over stars
import math
import numpy as np
import astropy.stats
import pandas as pd
from coord_functions import *
from matplotlib import pyplot as plt

readout_noise = 9
coord_param_num = 5  # number of polynom elements except the first which is just a number
color_param_num = 0  # number of polynom elements except the first which is just a number


# full_param_num = coord_param_num * 2 + color_param_num + 1
# full_param_num = coord_param_num + color_param_num + 1


def circle_ref(x0, y0, radius, gap, annul, ground):
    x0 = int(x0)
    y0 = int(y0)
    radius = int(radius)
    gap = int(gap)
    annul = int(annul)
    x_top = x0 + annul + 1
    x_bottom = x0 - annul
    y_top = y0 + annul + 1
    y_bottom = y0 - annul
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if (rad_val >= radius - 0.5 and rad_val <= radius + 0.5) or (
                    rad_val >= radius + gap - 0.5 and rad_val <= radius + gap + 0.5) or (
                    rad_val >= annul - 0.5 and rad_val <= annul + 0.5):
                ground[i][j] += 10 ** 4


def circle_check(x0_arr, y0_arr, radius, gap, annul, ground):
    x0_arr = np.int64(x0_arr)
    y0_arr = np.int64(y0_arr)
    radius = int(radius)
    gap = int(gap)
    annul = int(annul)
    for k in range(len(x0_arr)):
        x_top = x0_arr[k] + annul + 1
        x_bottom = x0_arr[k] - annul
        y_top = y0_arr[k] + annul + 1
        y_bottom = y0_arr[k] - annul
        if x_top > ground.shape[0]:
            x_top = ground.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > ground.shape[1]:
            y_top = ground.shape[1]
        if y_bottom < 0:
            y_bottom = 0
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0_arr[k] - i) ** 2 + (y0_arr[k] - j) ** 2))
                if (rad_val >= radius - 0.5 and rad_val <= radius + 0.5) or (
                        rad_val >= radius + gap - 0.5 and rad_val <= radius + gap + 0.5) or (
                        rad_val >= annul - 0.5 and rad_val <= annul + 0.5):
                    ground[i][j] += 10 ** 4


def ellipse_source(x0, y0, radius, gap, annul1, annul2, angle, ground):  # Ellipse-like annulus
    x0 = int(x0)
    y0 = int(y0)
    radius = int(radius)
    gap = int(gap)
    angle = float(angle)
    annul1 = int(annul1)
    annul2 = int(annul2)
    ext_radius = max(annul1, annul2)
    x_top = x0 + int(ext_radius) + 1
    x_bottom = x0 - int(ext_radius)
    y_top = y0 + int(ext_radius) + 1
    y_bottom = y0 - int(ext_radius)
    angle = angle / 180 * 3.14159
    print(angle)
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            axis1 = ((x0 - i) ** 2 * math.cos(angle) ** 2 - 2 * (x0 - i) * (y0 - j) * math.cos(angle) * math.sin(
                angle) + (y0 - j) ** 2 * math.sin(angle) ** 2) / annul1 ** 2
            axis2 = ((x0 - i) ** 2 * math.sin(angle) ** 2 + 2 * (x0 - i) * (y0 - j) * math.cos(angle) * math.sin(
                angle) + (y0 - j) ** 2 * math.cos(angle) ** 2) / annul2 ** 2
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            ell_val = axis1 + axis2
            if (ell_val >= 0.95 and ell_val <= 1.05) or (rad_val >= radius - 0.5 and rad_val <= radius + 0.5) or (
                    rad_val >= radius + gap - 0.5 and rad_val <= radius + gap + 0.5):
                # if (rad_val >= radius - 0.5 and rad_val <= radius + 0.5) or (rad_val >= radius + gap - 0.5 and rad_val <= radius + gap + 0.5):
                ground[i][j] += 10 ** 4


# Photometry functions
def photometry_ref(x0, y0, inn_radius, intermediate_increment, ext_radius, ground):
    source_pixel = 0
    total_signal = 0
    background_arr = []
    x0 = int(round(x0))
    y0 = int(round(y0))
    inn_radius = int(inn_radius)
    intermediate_increment = int(intermediate_increment)
    ext_radius = int(ext_radius)
    x_top = x0 + ext_radius + 1
    x_bottom = x0 - ext_radius
    y_top = y0 + ext_radius + 1
    y_bottom = y0 - ext_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if ground[i][j] > 0:
                if rad_val <= inn_radius:
                    total_signal += ground[i][j]
                    source_pixel += 1
                if rad_val > inn_radius + intermediate_increment and rad_val <= ext_radius + 1:
                    background_arr = np.append(background_arr, ground[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]
    source_signal = total_signal - background_per_pixel * source_pixel
    return source_signal


def adaptive_photometry_source(x0, y0, ground):  # Photometry with an automatic choice of aperture parameters(max SNR)
    max_check = 0
    ext_radius = 20
    gap = 5
    inn_radius = 1
    prev_snr = 0
    initial_check = 0
    output_inn_radius = 0
    x0 = int(round(x0))
    y0 = int(round(y0))
    x_top = x0 + ext_radius + 1
    x_bottom = x0 - ext_radius
    y_top = y0 + ext_radius + 1
    y_bottom = y0 - ext_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    while max_check == 0:
        source_pixel = 0
        background_arr = []
        source_arr = []
        total_signal = 0
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if rad_val <= inn_radius:
                        total_signal += ground[i][j]
                        source_arr = np.append(source_arr, ground[i][j])
                        source_pixel += 1
                    if rad_val > inn_radius + gap and rad_val <= ext_radius + 1:
                        background_arr = np.append(background_arr, ground[i][j])
        background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]
        # Part for discerning star-like shape of a studied source
        source_err_arr = np.sqrt(source_arr + readout_noise ** 2)
        SNR_arr = (source_arr - background_per_pixel) / source_err_arr
        one_point_five_sigma_pixels_arr = np.count_nonzero(SNR_arr >= 1.5)
        if one_point_five_sigma_pixels_arr / (source_arr.shape[0]) <= 0.5:
            output_inn_radius = inn_radius
            break
        source_signal = total_signal - background_per_pixel * source_pixel
        SNR = source_signal / math.sqrt(total_signal + readout_noise ** 2 * source_pixel)
        if inn_radius >3 and SNR <=2:
            output_inn_radius = 3
            break
        if initial_check == 0:
            prev_snr = SNR
            initial_check = 1
            inn_radius += 1
        else:
            # print(prev_snr, SNR)
            if prev_snr > SNR:
                max_check = 1
                output_inn_radius = inn_radius - 1
            else:
                prev_snr = SNR
                inn_radius += 1
    # Estimating a max SNR photometry
    inn_radius = output_inn_radius
    source_pixel = 0
    background_arr = []
    total_signal = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if ground[i][j] > 0:
                if rad_val <= inn_radius:
                    total_signal += ground[i][j]
                    source_pixel += 1
                if rad_val > inn_radius + gap and rad_val <= ext_radius + 1:
                    background_arr = np.append(background_arr, ground[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]
    source_signal = total_signal - background_per_pixel * source_pixel
    source_err = math.sqrt(total_signal + readout_noise ** 2 * source_pixel)

    # Calculating a 3 sigma limit
    snr_limit = 3
    #source_pixel = 5
    signal_limit = (snr_limit ** 2 + sqrt(snr_limit ** 4 + 4 * source_pixel * snr_limit ** 2 * (background_per_pixel + readout_noise ** 2))) / 2
    return source_signal, source_err, output_inn_radius, signal_limit


def photometry_source(x0, y0, inn_radius, intermediate_increment, ext_radius, ground):
    source_pixel = 0
    background_arr = []
    total_signal = 0
    x0 = int(round(x0))
    y0 = int(round(y0))
    inn_radius = int(inn_radius)
    intermediate_increment = int(intermediate_increment)
    ext_radius = int(ext_radius)
    x_top = x0 + ext_radius + 1
    x_bottom = x0 - ext_radius
    y_top = y0 + ext_radius + 1
    y_bottom = y0 - ext_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if ground[i][j] > 0:
                if rad_val <= inn_radius:
                    total_signal += ground[i][j]
                    source_pixel += 1
                if rad_val > inn_radius + intermediate_increment and rad_val <= ext_radius + 1:
                    background_arr = np.append(background_arr, ground[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]
    source_signal = total_signal - background_per_pixel * source_pixel
    SNR = source_signal / math.sqrt(total_signal + readout_noise ** 2 * source_pixel)
    source_err = source_signal / SNR

    # Calculating a 3 sigma limit
    snr_limit = 3
    source_signal = 5
    signal_limit = (snr_limit ** 2 + sqrt(
        snr_limit ** 4 + 4 * source_pixel * snr_limit ** 2 * (background_per_pixel + readout_noise ** 2))) / 2
    return source_signal, source_err, signal_limit


def photometry_snr(x0, y0, inn_radius, intermediate_increment, ext_radius, ground, back_level):
    source_pixel = 0
    background_arr = []
    total_signal = 0
    x0 = int(round(x0))
    y0 = int(round(y0))
    inn_radius = int(inn_radius)
    intermediate_increment = int(intermediate_increment)
    ext_radius = int(ext_radius)
    x_top = x0 + ext_radius + 1
    x_bottom = x0 - ext_radius
    y_top = y0 + ext_radius + 1
    y_bottom = y0 - ext_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    if back_level != 0:
        background_per_pixel = back_level
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if rad_val <= inn_radius:
                        total_signal += ground[i][j]
                        source_pixel += 1
    else:
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if rad_val <= inn_radius:
                        total_signal += ground[i][j]
                        source_pixel += 1
                    if rad_val > inn_radius + intermediate_increment and rad_val <= ext_radius + 1:
                        background_arr = np.append(background_arr, ground[i][j])
        background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]

    source_signal = total_signal - background_per_pixel * source_pixel
    #print(background_per_pixel, source_signal)
    SNR = source_signal / math.sqrt(total_signal + readout_noise ** 2 * source_pixel)
    return SNR

def prec_photometry_snr(x0, y0, inn_radius, intermediate_increment, ext_radius, ground):
    source_pixel = 0
    background_arr = []
    total_signal = 0
    x0 = int(round(x0))
    y0 = int(round(y0))
    inn_radius = int(inn_radius)
    intermediate_increment = int(intermediate_increment)
    ext_radius = int(ext_radius)
    x_top = x0 + ext_radius + 1
    x_bottom = x0 - ext_radius
    y_top = y0 + ext_radius + 1
    y_bottom = y0 - ext_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    if back_level != 0:
        background_per_pixel = back_level
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if rad_val <= inn_radius:
                        total_signal += ground[i][j]
                        source_pixel += 1
    else:
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if rad_val <= inn_radius:
                        total_signal += ground[i][j]
                        source_pixel += 1
                    if rad_val > inn_radius + intermediate_increment and rad_val <= ext_radius + 1:
                        background_arr = np.append(background_arr, ground[i][j])
        background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]

    source_signal = total_signal - background_per_pixel * source_pixel
    #print(background_per_pixel, source_signal)
    SNR = source_signal / math.sqrt(total_signal + readout_noise ** 2 * source_pixel)
    return SNR

def pixel_brightness(x0, y0, inn_radius, ground):
    source_pixel = 0
    background_arr = []
    total_signal = 0
    x0 = int(round(x0))
    y0 = int(round(y0))
    inn_radius = int(inn_radius)
    x_top = x0 + inn_radius + 1
    x_bottom = x0 - inn_radius
    y_top = y0 + inn_radius + 1
    y_bottom = y0 - inn_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    peak_brightness = 0
    if back_level != 0:
        background_per_pixel = back_level
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                    if ground[i][j] >= peak_brightness:
                        peak_brightness = ground
    return peak_brightness

def photometry_source_ell(x0, y0, inn_radius, gap, annul1, annul2, angle, ground):
    source_signal = 0
    source_pixel = 0
    background_per_pixel = 0
    background_arr = []
    total_signal = 0
    x0 = int(x0)
    y0 = int(y0)
    inn_radius = int(inn_radius)
    gap = int(gap)
    annul1 = int(annul1)
    annul2 = int(annul2)
    ext_radius = max(annul1, annul2)
    angle = float(angle) / 180 * 3.14159
    x_top = x0 + ext_radius + 1
    x_bottom = x0 - ext_radius
    y_top = y0 + ext_radius + 1
    y_bottom = y0 - ext_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            axis1 = ((x0 - i) ** 2 * math.cos(angle) ** 2 - 2 * (x0 - i) * (y0 - j) * math.cos(angle) * math.sin(
                angle) + (y0 - j) ** 2 * math.sin(angle) ** 2) / annul1 ** 2
            axis2 = ((x0 - i) ** 2 * math.sin(angle) ** 2 + 2 * (x0 - i) * (y0 - j) * math.cos(angle) * math.sin(
                angle) + (y0 - j) ** 2 * math.cos(angle) ** 2) / annul2 ** 2
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            ell_val = axis1 + axis2
            if ground[i][j] > 0:
                if rad_val <= inn_radius:
                    total_signal += ground[i][j]
                    source_pixel += 1
                if rad_val > inn_radius + gap and ell_val <= 1:
                    background_arr = np.append(background_arr, ground[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]
    source_signal = total_signal - background_per_pixel * source_pixel
    # if source_signal<0:
    # source_signal=
    SNR = source_signal / math.sqrt(total_signal)
    source_err = source_signal / SNR
    # print(source_pixel)
    return source_signal, source_err


def photometry_ell_source(x0, y0, inn_radius2, inn_radius1, intermediate_increment, ext_radius, ground):
    source_signal = 0
    source_pixel = 0
    background_per_pixel = 0
    background_arr = []
    total_signal = 0
    x0 = int(x0)
    y0 = int(y0)
    inn_radius1 = int(inn_radius1)
    inn_radius2 = int(inn_radius2)
    inn_radius = max(inn_radius1, inn_radius2)
    intermediate_increment = int(intermediate_increment)
    ext_radius = int(ext_radius)
    x_top = x0 + ext_radius + 1
    x_bottom = x0 - ext_radius
    y_top = y0 + ext_radius + 1
    y_bottom = y0 - ext_radius
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            axis1 = (x0 - i) ** 2 / inn_radius1 ** 2
            axis2 = (y0 - j) ** 2 / inn_radius2 ** 2
            rad_val = axis1 + axis2
            centre_radius = math.sqrt((x0 - i) ** 2 + (y0 - j) ** 2)
            if ground[i][j] > 0:
                if rad_val <= 1:
                    total_signal += ground[i][j]
                    source_pixel += 1
                if centre_radius > inn_radius + intermediate_increment and centre_radius <= ext_radius + 1:
                    background_arr = np.append(background_arr, ground[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]
    source_signal = total_signal - background_per_pixel * source_pixel
    # if source_signal<0:
    # source_signal=
    SNR = source_signal / math.sqrt(total_signal)
    source_err = source_signal / SNR
    return source_signal, source_err


def photometry_check(x0_arr, y0_arr, inn_radius, intermediate_increment, ext_radius, ground):
    x0_arr = np.int64(np.round(x0_arr))
    y0_arr = np.int64(np.round(y0_arr))
    # print(x0_arr)
    source_signal = np.zeros_like(x0_arr)
    source_signal_err = np.zeros_like(x0_arr)
    inn_radius = int(inn_radius)
    intermediate_increment = int(intermediate_increment)
    ext_radius = int(ext_radius)
    for k in range(len(x0_arr)):
        source_pixel = 0
        background_per_pixel = 0
        background_arr = []
        total_signal = 0
        x_top = x0_arr[k] + ext_radius + 1
        x_bottom = x0_arr[k] - ext_radius
        y_top = y0_arr[k] + ext_radius + 1
        y_bottom = y0_arr[k] - ext_radius
        if x_top > ground.shape[0]:
            x_top = ground.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > ground.shape[1]:
            y_top = ground.shape[1]
        if y_bottom < 0:
            y_bottom = 0
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0_arr[k] - i) ** 2 + (y0_arr[k] - j) ** 2))
                if ground[i][j] > 0:
                    if rad_val <= inn_radius:
                        total_signal += ground[i][j]
                        source_pixel += 1
                    if rad_val > inn_radius + intermediate_increment and rad_val <= ext_radius + 1:
                        background_arr = np.append(background_arr, ground[i][j])
        background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=1)[0]
        # print(total_signal, background_per_pixel, source_pixel, x0_arr[k], y0_arr[k])
        source_signal[k] = total_signal - background_per_pixel * source_pixel
        #if source_signal[k] <=0:
            #print(background_arr)
            #print(total_signal)
        #print(source_signal[k])
        source_signal_err[k] = math.sqrt(total_signal + readout_noise ** 2 * source_pixel)
    return source_signal, source_signal_err


# Calculation of magnitude and its variation
def magnitude(mag_ref, photo_ref, photo_source, photo_check):
    mag_source = np.zeros_like(photo_ref)
    mag_check = np.zeros_like(photo_check)
    mag_source = mag_ref - 2.5 * np.log10(photo_source / photo_ref)
    #print(np.log10(photo_source / photo_ref))
    for k in range(photo_check.shape[1]):
        for j in range(photo_check.shape[0]):
            if photo_check[j][k] > 0:
                mag_check[j][k] = mag_ref - 2.5 * math.log10(photo_check[j][k] / photo_ref[j])

    return mag_source, mag_check


def magnitude_obj(mag_ref, photo_ref, photo_source):
    mag_source = mag_ref - 2.5 * math.log10(photo_source / photo_ref)
    return mag_source


def magnitude_err(photo_source, photo_err, mag_check, check_photo, check_err):
    SNR = photo_source / photo_err
    SNR_check = check_photo / check_err
    check_average = np.zeros(mag_check.shape[1])
    mag_check_var_err = np.zeros_like(mag_check)
    mag_check_signal_err = np.zeros_like(mag_check)
    mag_source_var_err = np.zeros_like(photo_source)
    mag_source_signal_err = np.zeros_like(photo_source)
    for k in range(mag_check.shape[1]):
        check_total_mag = 0
        for j in range(mag_check.shape[0]):
            check_total_mag += mag_check[j][k]
        check_average[k] = check_total_mag / mag_check.shape[0]
        for j in range(mag_check.shape[0]):
            mag_check_var_err[j][k] = abs(check_average[k] - mag_check[j][k])
            mag_check_signal_err[j][k] = 1 / SNR_check[j][k]
    # print(check_average)
    for j in range(mag_check.shape[0]):
        for k in range(mag_check.shape[1]):
            mag_source_var_err[j] += (mag_check[j][k] - check_average[k]) ** 2
        mag_source_var_err[j] = math.sqrt(mag_source_var_err[j] / (mag_check.shape[1] - 1))
    mag_source_signal_err = 1 / SNR
    return mag_source_var_err, mag_source_signal_err, mag_check_var_err, mag_check_signal_err


# Taking into account atmospheric variations
def var_correct(check_mag, check_mag_err, check_true_mag,
                check_color):  # Fitting by a 2 order polynom (a0+a1*x+a2*x**2) taking into account weights
    param_num = 4
    # Defining arrays
    check_len = check_color.shape[0]
    photo_len = check_mag.shape[0]
    x_arr = np.zeros((check_len, param_num))
    var_mag = np.zeros_like(check_color)
    weight_arr = np.zeros((check_len, check_len))
    # Fit parameters
    alpha = np.zeros(photo_len)
    beta = np.zeros(photo_len)
    gamma = np.zeros(photo_len)
    delta = np.zeros(photo_len)
    alpha_err = np.zeros(photo_len)
    beta_err = np.zeros(photo_len)
    gamma_err = np.zeros(photo_len)
    delta_err = np.zeros(photo_len)
    # Main body
    for k in range(photo_len):
        for i in range(check_len):
            if k == 0:
                x_arr[i][0] = 1
                x_arr[i][1] = check_color[i]
                x_arr[i][2] = check_color[i] ** 2
                x_arr[i][3] = check_color[i] ** 3
            var_mag[i] = check_mag[k][i] - check_true_mag[i]
            weight_arr[i][i] = 1 / check_mag_err[k][i] ** 2
        x_arr_transpose = np.transpose(x_arr)
        x_mod_arr = np.dot(x_arr_transpose, weight_arr)
        y_mod_arr = np.dot(x_mod_arr, var_mag)
        x_mod_arr = np.dot(x_mod_arr, x_arr)
        x_inv_arr = np.linalg.inv(x_mod_arr)
        vector_param = np.dot(x_inv_arr, y_mod_arr)
        alpha[k] = vector_param[0]
        beta[k] = vector_param[1]
        gamma[k] = vector_param[2]
        delta[k] = vector_param[3]
        # Error estimation
        sum_fun = 0
        for i in range(check_len):
            sum_fun += weight_arr[i][i] * (
                    var_mag[i] - alpha[k] - beta[k] * x_arr[i][1] - gamma[k] * x_arr[i][2] - delta[k] * x_arr[i][
                3]) ** 2
        parameter_variance = sum_fun / (check_len - param_num) * x_inv_arr
        alpha_err[k] = math.sqrt(parameter_variance[0][0])
        beta_err[k] = math.sqrt(parameter_variance[1][1])
        gamma_err[k] = math.sqrt(parameter_variance[2][2])
        delta_err[k] = math.sqrt(parameter_variance[3][3])

    return alpha, alpha_err, beta, beta_err, gamma, gamma_err, delta, delta_err


# Taking into account CCD sensitivity variations (imperfect flats or something like that)
def coord_var_correct(check_mag, check_mag_err, check_true_mag, check_x, check_y,
                      frame_shape):  # Fitting by a 4 order polynom on each coordinate (
    # a0+a1*x+a2*x**2+a3*x**3+a4*x**4) taking into account weights
    param_num = 7
    check_mid_x = frame_shape[0] / 2
    check_mid_y = frame_shape[1] / 2
    # Defining arrays
    check_len = check_true_mag.shape[0]
    photo_len = check_mag.shape[0]
    x_arr = np.zeros((check_len, param_num * 2 - 1))
    coord_var_mag = np.zeros_like(check_true_mag)
    weight_arr = np.zeros((check_len, check_len))
    # Fit parameters
    coord_par = np.zeros((photo_len, param_num * 2 - 1))
    coord_par_err = np.zeros((photo_len, param_num * 2 - 1))
    # Main body
    for k in range(photo_len):
        for i in range(check_len):
            if k == 0:
                x_arr[i][0] = 1
                power_check = 1  # for a correct power index in creating the X array
                for m in range(1, param_num * 2 - 1, 2):
                    x_arr[i][m] = (check_x[i] - check_mid_x) ** power_check
                    x_arr[i][m + 1] = (check_y[i] - check_mid_y) ** power_check
                    power_check += 1
            coord_var_mag[i] = check_mag[k][i] - check_true_mag[i]
            weight_arr[i][i] = 1 / check_mag_err[k][i] ** 2
        # print(x_arr)
        x_arr_transpose = np.transpose(x_arr)
        x_mod_arr = np.dot(x_arr_transpose, weight_arr)
        y_mod_arr = np.dot(x_mod_arr, coord_var_mag)
        x_mod_arr = np.dot(x_mod_arr, x_arr)
        x_inv_arr = np.linalg.inv(x_mod_arr)
        vector_param = np.dot(x_inv_arr, y_mod_arr)
        for m in range(0, param_num * 2 - 1):
            coord_par[k][m] = vector_param[m]
        # Error estimation
        sum_arr = coord_var_mag - np.dot(x_arr, vector_param)
        sum_arr_trans = np.transpose(sum_arr)
        minim_sum = np.dot(np.dot(sum_arr_trans, weight_arr), sum_arr)
        parameter_variance = minim_sum / (check_len - param_num * 2 + 1) * x_inv_arr
        for m in range(0, param_num * 2 - 1):
            if parameter_variance[m][m] >= 0:
                coord_par_err[k][m] = math.sqrt(parameter_variance[m][m])
            else:
                coord_par_err[k][m] = -math.sqrt(-parameter_variance[m][m])

    return coord_par, coord_par_err


# Taking into account CCD sensitivity variations as well as atmospheric effects
def full_var_correct(check_mag, check_mag_err, check_true_mag, check_color, check_x, check_y,
                     frame_shape, coord_poly_num, color_param_num):  # Fitting by a 4 order polynom on each coordinate (
    # a0+a1*x+a2*x**2+a3*x**3+a4*x**4) taking into account weights
    print(color_param_num)
    full_param_num = (coord_poly_num ** 2 + 3 * coord_poly_num) / 2 + color_param_num + 1
    print("Number of all parameters: " + "{}".format(full_param_num))
    check_mid_x = frame_shape[0] / 2
    check_mid_y = frame_shape[1] / 2
    # Defining arrays
    check_len = len(check_color)
    photo_len = check_mag.shape[0]
    alpha = 10 ** (-3)



    # Fit parameters
    coord_par = np.zeros((photo_len, full_param_num))
    coord_par_err = np.zeros((photo_len, full_param_num))
    # Main body
    for k in range(photo_len):
        x_input_arr = np.zeros((check_len, full_param_num))
        coord_var_mag = np.zeros_like(check_color)
        weight_arr = np.zeros((check_len, check_len))
        for i in range(check_len):
            x_input_arr[i][0] = 1
            power_check = 1  # for a correct power index in creating the X array
            for m in range(1, color_param_num + 1):
                x_input_arr[i][m] = check_color[i] ** power_check
                power_check += 1
            # print(x_arr[0])
            power_check = 1
            # for m in range(color_param_num + 1, full_param_num, 2):
            #    x_arr[i][m] = (check_x[k][i] - check_mid_x) ** power_check
            #    x_arr[i][m + 1] = (check_y[k][i] - check_mid_y) ** power_check
            #   power_check += 1
            # for m in range(color_param_num + 1, full_param_num):
            #    x_arr[i][m] = (check_x[k][i] - check_mid_x + check_y[k][i] - check_mid_y) ** power_check
            #    power_check += 1
            place_index = color_param_num + 1
            for poly_num in range(1, coord_poly_num + 1):
                for coord_ind in range(poly_num + 1):
                    x_input_arr[i][place_index] = (check_x[k][i] - check_mid_x) ** (poly_num - coord_ind) * (
                                check_y[k][i] - check_mid_y) ** (coord_ind)
                    place_index += 1
            coord_var_mag[i] = check_mag[k][i] - check_true_mag[k][i]
            weight_arr[i][i] = 1 / check_mag_err[k][i] ** 2
        if k == 0:
            print(coord_var_mag)
        # print(x_arr[0])
        # print(x_arr)
        x_arr_transpose = np.transpose(x_input_arr)
        x_mod_arr = np.dot(x_arr_transpose, weight_arr)
        y_mod_arr = np.dot(x_mod_arr, coord_var_mag)
        x_mod_arr = np.dot(x_mod_arr, x_input_arr)
        x_inv_arr = np.linalg.inv(x_mod_arr + alpha * np.identity(full_param_num))
        #vector_param = np.dot(x_inv_arr, y_mod_arr)
        # Gauss elimination
        x_converted_arr, y_converted_arr = gauss_elimination(x_mod_arr, y_mod_arr)
        for parameter_ind in range(x_converted_arr.shape[0] - 1, -1, -1):
            substract_value = 0
            if parameter_ind != x_converted_arr.shape[0] - 1:
                for substract_ind in range(parameter_ind + 1, x_converted_arr.shape[1]):
                    substract_value += coord_par[k][substract_ind] * x_converted_arr[parameter_ind][
                        substract_ind]
            right_value = y_converted_arr[parameter_ind] - substract_value
            coord_par[k][parameter_ind] = right_value / x_converted_arr[parameter_ind][parameter_ind]
        #for m in range(0, full_param_num):
        #    coord_par[k][m] = vector_param[m]
        # Error estimation
        sum_arr = coord_var_mag - np.dot(x_input_arr, coord_par[k])
        sum_arr_trans = np.transpose(sum_arr)
        minim_sum = np.dot(np.dot(sum_arr_trans, weight_arr), sum_arr)
        #x_error_converted,
        parameter_variance = minim_sum / (check_len - full_param_num) * x_inv_arr
        for m in range(0, full_param_num):
            if parameter_variance[m][m] >= 0:
                coord_par_err[k][m] = math.sqrt(parameter_variance[m][m])
            else:
                coord_par_err[k][m] = -math.sqrt(-parameter_variance[m][m])
    print('Sum of deviation: ' + "{}".format(minim_sum))
    return coord_par, coord_par_err


# Taking into account CCD sensitivity variations as well as atmospheric effects
def all_frame_full_var_correct(check_mag, check_mag_err, check_true_mag, check_color, check_x, check_y,
                               frame_shape):  # Fitting by a 4 order polynom on each coordinate (a0+a1*x+a2*x**2+a3*x**3+a4*x**4) taking into account weights
    full_param_num = 19
    check_mid_x = frame_shape[0] / 2
    check_mid_y = frame_shape[1] / 2
    # Defining arrays
    check_len = check_true_mag.shape[0]
    photo_len = check_mag.shape[0]
    alpha = 10 ** -3
    x_arr = np.zeros((check_len, full_param_num))
    coord_var_mag = np.zeros_like(check_true_mag)
    weight_arr = np.zeros((check_len, check_len))
    # Fit parameters
    coord_par = np.zeros((photo_len, full_param_num))
    coord_par_err = np.zeros((photo_len, full_param_num))
    # Main body
    for k in range(photo_len):
        for i in range(check_len):
            if k == 0:
                x_arr[i][0] = 1
                power_check = 1  # for a correct power index in creating the X array
                for m in range(1, color_param_num + 1):
                    x_arr[i][m] = check_color[i] ** power_check
                    power_check += 1
                power_check = 1
                for m in range(color_param_num + 1, full_param_num, 2):
                    x_arr[i][m] = (check_x[i] - check_mid_x) ** power_check
                    x_arr[i][m + 1] = (check_y[i] - check_mid_y) ** power_check
                    power_check += 1
            coord_var_mag[i] = check_mag[k][i] - check_true_mag[i]
            weight_arr[i][i] = 1 / check_mag_err[k][i] ** 2
        # print(x_arr)
        x_arr_transpose = np.transpose(x_arr)
        x_mod_arr = np.dot(x_arr_transpose, weight_arr)
        y_mod_arr = np.dot(x_mod_arr, coord_var_mag)
        x_mod_arr = np.dot(x_mod_arr, x_arr)
        x_inv_arr = np.linalg.inv(x_mod_arr + alpha * np.identity(full_param_num))
        vector_param = np.dot(x_inv_arr, y_mod_arr)
        for m in range(0, full_param_num):
            coord_par[k][m] = vector_param[m]
        # Error estimation
        sum_arr = coord_var_mag - np.dot(x_arr, vector_param)
        sum_arr_trans = np.transpose(sum_arr)
        minim_sum = np.dot(np.dot(sum_arr_trans, weight_arr), sum_arr)
        parameter_variance = minim_sum / (check_len - full_param_num) * x_inv_arr
        for m in range(0, full_param_num):
            if parameter_variance[m][m] >= 0:
                coord_par_err[k][m] = math.sqrt(parameter_variance[m][m])
            else:
                coord_par_err[k][m] = -math.sqrt(-parameter_variance[m][m])
    print('Sum of squared deviations of magnitude (last frame): ' + "{}".format(minim_sum))
    return coord_par, coord_par_err


# Correcting magnitudes using obtained fit
def check_mag_correct(check_mag, check_color, check_x, check_y, corr_param, frame_shape, coord_poly_num, color_param_num):
    full_param_num = (coord_poly_num ** 2 + 3 * coord_poly_num) / 2 + color_param_num + 1
    check_len = check_mag.shape[1]
    photo_len = check_mag.shape[0]
    check_mid_x = frame_shape[0] / 2
    check_mid_y = frame_shape[1] / 2
    check_mag_corr = np.copy(check_mag)
    for k in range(photo_len):
        for i in range(check_len):
            check_mag_corr[k][i] -= corr_param[k][0]
            color_power_index = 1
            coord_power_index = 1
            for m in range(1, color_param_num + 1):
                check_mag_corr[k][i] -= corr_param[k][m] * check_color[i] ** color_power_index
                color_power_index += 1
            # for m in range(color_param_num + 1, full_param_num, 2):
            #    check_mag_corr[k][i] -= corr_param[k][m] * (check_x[k][i] - check_mid_x) ** coord_power_index
            #    check_mag_corr[k][i] -= corr_param[k][m + 1] * (check_y[k][i] - check_mid_y) ** coord_power_index
            #    coord_power_index += 1
            # for m in range(color_param_num + 1, full_param_num):
            #    check_mag_corr[k][i] -= corr_param[k][m] * (check_x[k][i] - check_mid_x + check_y[k][i] - check_mid_y) ** coord_power_index
            #    coord_power_index += 1
            place_index = color_param_num + 1

            for poly_num in range(1, coord_poly_num + 1):
                for coord_ind in range(poly_num + 1):
                    check_mag_corr[k][i] -= corr_param[k][place_index] * (check_x[k][i] - check_mid_x) ** (
                                poly_num - coord_ind) * (
                                                    check_y[k][i] - check_mid_y) ** (coord_ind)
                    place_index += 1
            if k == 0:
                if i == 0:
                    print(place_index)
    return check_mag_corr


def source_mag_correct(source_mag, source_color, source_x, source_y, corr_param, frame_shape, coord_poly_num, color_param_num):
    photo_len = source_mag.shape[0]
    check_mid_x = frame_shape[0] / 2
    check_mid_y = frame_shape[1] / 2
    source_mag_corr = np.copy(source_mag)
    for k in range(photo_len):
        source_mag_corr[k] -= corr_param[k][0]
        color_power_index = 1
        coord_power_index = 1
        for m in range(1, color_param_num + 1):
            source_mag_corr[k] -= corr_param[k][m] * source_color ** color_power_index
            color_power_index += 1
        # for m in range(color_param_num + 1, full_param_num, 2): source_mag_corr[k] -= corr_param[k][m] * (source_x[
        # k] - check_mid_x) ** coord_power_index source_mag_corr[k] -= corr_param[k][m + 1] * (source_y[k] -
        # check_mid_y) ** coord_power_index coord_power_index += 1 for m in range(color_param_num + 1,
        # full_param_num): source_mag_corr[k] -= corr_param[k][m] * (source_x[k] - check_mid_x + source_y[k] -
        # check_mid_y) ** coord_power_index coord_power_index += 1
        place_index = color_param_num + 1
        for poly_num in range(1, coord_poly_num + 1):
            for coord_ind in range(poly_num + 1):
                source_mag_corr[k] -= corr_param[k][place_index] * (source_x[k] - check_mid_x) ** (
                        poly_num - coord_ind) * (
                                              source_y[k] - check_mid_y) ** (coord_ind)
                place_index += 1
    return source_mag_corr


# Finding centroid of an object
def obj_centroid(check_x, check_y, radius, gap, annul, ground, back_level):
    x0 = int(round(check_x))
    y0 = int(round(check_y))
    radius = int(radius)
    gap = int(gap)
    annul = int(annul)
    x_top = x0 + annul + 1
    x_bottom = x0 - annul
    y_top = y0 + annul + 1
    y_bottom = y0 - annul
    background_arr = []
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    # Calculating background
    if back_level !=0:
        background_per_pixel = back_level
    else:
        for i in range(x_bottom, x_top):
           for j in range(y_bottom, y_top):
               rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
               if ground[i][j] > 0:
                   if rad_val > radius + gap and rad_val <= annul + 1:
                       background_arr = np.append(background_arr, ground[i][j])
        background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=2)[0]
    for iter in range(2):
        # if x0 == 539 and (y0 == 168 or y0 == 167):
        #    print(background_arr)
        # print(background_per_pixel)
        x_top = x0 + radius + gap + 1
        x_bottom = x0 - radius - gap
        y_top = y0 + radius + gap + 1
        y_bottom = y0 - radius - gap
        if x_top > ground.shape[0]:
            x_top = ground.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > ground.shape[1]:
            y_top = ground.shape[1]
        if y_bottom < 0:
            y_bottom = 0
        # Calculating centroid
        source_signal = 0
        x_signal = 0
        y_signal = 0
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if (rad_val <= radius and rad_val <= radius + gap):
                        source_signal += ground[i][j] - background_per_pixel
                        x_signal += (ground[i][j] - background_per_pixel) * i
                        y_signal += (ground[i][j] - background_per_pixel) * j
        if source_signal <= 0.1:
            x_centre = x0
            y_centre = y0
            break
        else:
            x_centre = x_signal / source_signal
            y_centre = y_signal / source_signal
            # x0 = x_centre
            # y0 = y_centre
        # print(x_centre - check_x)
        # print(y_centre - check_y)
        check_x = x_centre
        check_y = y_centre
        # print(x_centre, y_centre)
    return x_centre, y_centre


def test_centroid(check_x, check_y, radius, gap, annul, ground):
    x0 = int(round(check_x))
    y0 = int(round(check_y))
    radius = int(radius)
    gap = int(gap)
    annul = int(annul)
    x_top = x0 + annul + 1
    x_bottom = x0 - annul
    y_top = y0 + annul + 1
    y_bottom = y0 - annul
    background_arr = []
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    # Calculating background
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if ground[i][j] > 0:
                if rad_val > radius + gap and rad_val <= annul + 1:
                    background_arr = np.append(background_arr, ground[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=2)[0]
    #print(background_per_pixel)
    for iter in range(3):
        x_top = int(x0 + radius + gap + 1)
        x_bottom = int(x0 - radius - gap)
        y_top = int(y0 + radius + gap + 1)
        y_bottom = int(y0 - radius - gap)
        if x_top > ground.shape[0]:
            x_top = ground.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > ground.shape[1]:
            y_top = ground.shape[1]
        if y_bottom < 0:
            y_bottom = 0
        # Calculating centroid
        source_signal = 0
        x_signal = 0
        y_signal = 0
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if (rad_val <= radius * 3 and rad_val <= radius + gap):
                        pixel_signal = ground[i][j] - background_per_pixel
                        source_signal += pixel_signal
                        x_signal += pixel_signal * i
                        y_signal += pixel_signal * j
        #print(source_signal)
        #print(x_signal)
        #print(y_signal)
        if source_signal <= 0.1:
            x_centre = x0
            y_centre = y0
            break
        else:
            x_centre = x_signal / source_signal
            y_centre = y_signal / source_signal
            # x0 = x_centre
            # y0 = y_centre
        # print(x_centre - check_x)
        # print(y_centre - check_y)
        x0 = x_centre
        y0 = y_centre
        # print(x_centre, y_centre)
    return x_centre, y_centre


# Functions for finding the brightest sources on a frame

def fast_obj_centroid(input_x, input_y, back_level, frame_data):
    x0 = int(round(input_x))
    y0 = int(round(input_y))
    radius = 4
    gap = 5

    for iter in range(3):
        x_top = int(x0 + radius + gap + 1)
        x_bottom = int(x0 - radius - gap)
        y_top = int(y0 + radius + gap + 1)
        y_bottom = int(y0 - radius - gap)
        if x_top > frame_data.shape[0]:
            x_top = frame_data.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > frame_data.shape[1]:
            y_top = frame_data.shape[1]
        if y_bottom < 0:
            y_bottom = 0
        # Calculating centroid
        source_signal = 0
        x_signal = 0
        y_signal = 0
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if frame_data[i][j] > 0:
                    if (rad_val <= radius * 3 and rad_val <= radius + gap):
                        pixel_signal = frame_data[i][j] - back_level
                        source_signal += pixel_signal
                        x_signal += pixel_signal * i
                        y_signal += pixel_signal * j
        if source_signal <= 1:
            x_centre = x0
            y_centre = y0
            break
        else:
            x_centre = x_signal / source_signal
            y_centre = y_signal / source_signal
        x0 = x_centre
        y0 = y_centre
    return x_centre, y_centre


def fill_area(input_x, input_y, frame_data, back_level):
    x0 = int(round(input_x))
    y0 = int(round(input_y))
    annul = 25

    x_top = int(x0 + annul + 1)
    x_bottom = int(x0 - annul)
    y_top = int(y0 + annul + 1)
    y_bottom = int(y0 - annul)
    if x_top > frame_data.shape[0]:
        x_top = frame_data.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > frame_data.shape[1]:
        y_top = frame_data.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            frame_data[i][j] = back_level


def BrightStarFinder(frame_data, num_of_bright):
    # print('')
    initial_back_level = np.mean(frame_data)
    back_level_arr = np.zeros((10, 10))
    # Creating a background signal array
    # Check if the background is flat
    #background_1 = np.mean(frame_data[0:100, 0:100])
    #background_2 = np.mean(frame_data[900:1000, frame_data.shape[1] - 100:frame_data.shape[1]])
    frame_size_x = frame_data.shape[0]
    frame_size_y = frame_data.shape[1]
    background_1 = astropy.stats.sigma_clipped_stats(frame_data[0:int(frame_size_x / 10), 0: int(frame_size_y / 10)], sigma=1)[0]
    background_2 = astropy.stats.sigma_clipped_stats(frame_data[int(frame_size_x * 0.9): int(frame_size_x), int(frame_size_y * 0.9): frame_size_x], sigma=1)[0]
    frame_mid_x = frame_data.shape[0] / 2
    frame_mid_y = frame_data.shape[1] / 2
    median_level = astropy.stats.sigma_clipped_stats(frame_data, sigma=1.5)[0]
    frame_data = np.where(frame_data > median_level * 0.1, frame_data, median_level)
    back_level = astropy.stats.sigma_clipped_stats(frame_data, sigma=1.5)[0]
    #plt.figure(2)
    #im = plt.imshow(frame_data)
    #plt.colorbar(im)
    #c_min = back_level
    #plt.clim([c_min, c_min + math.sqrt(back_level + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    if abs(background_2 - background_1) / math.sqrt(background_2 + readout_noise ** 2) >= 3:
        print('Correcting background level into uniform one')
        for i_ind in range(10):
            for j_ind in range(10):
                x_lower = i_ind * int(frame_size_x / 10)
                x_top = (i_ind + 1) * int(frame_size_x / 10)
                y_lower = j_ind * int(frame_size_y / 10)
                y_top = (j_ind + 1) * int(frame_size_y / 10)
                if y_top >= frame_size_y:
                    y_top = frame_size_y - 1
                if x_top >= frame_size_x:
                    x_top = frame_size_x - 1
                # print(x_lower, x_top, y_lower, y_top)
                back_level_arr[i_ind][j_ind] = np.mean(frame_data[x_lower:x_top, y_lower:y_top])
        # Finding a back_level for every pixel on the frame using least squares
        # print(back_level_arr)
        # print('aaaa')
        back_level_par = np.zeros(6)
        x_arr = np.zeros((100, 6))
        y_arr = np.zeros(100)
        for i_ind in range(10):
            for j_ind in range(10):
                x_arr[i_ind * 10 + j_ind][0] = 1
                x_arr[i_ind * 10 + j_ind][1] = i_ind * int(frame_size_x / 10) + int(frame_size_x / 10 / 2)
                x_arr[i_ind * 10 + j_ind][2] = j_ind * int(frame_size_y / 10) + int(frame_size_y / 10 / 2)
                x_arr[i_ind * 10 + j_ind][3] = (i_ind * int(frame_size_x / 10) + int(frame_size_x / 10 / 2)) ** 2
                x_arr[i_ind * 10 + j_ind][4] = (j_ind * int(frame_size_y / 10) + int(frame_size_y / 10 / 2)) * (i_ind * int(frame_size_x / 10) + int(frame_size_x / 10 / 2))
                x_arr[i_ind * 10 + j_ind][5] = (j_ind * int(frame_size_y / 10) + int(frame_size_y / 10 / 2)) ** 2
                y_arr[i_ind * 10 + j_ind] = back_level_arr[i_ind][j_ind]
        # print(x_arr)
        # print('bbbb')
        x_arr_transpose = np.transpose(x_arr)
        y_mod_arr = np.dot(x_arr_transpose, y_arr)
        x_mod_arr = np.dot(x_arr_transpose, x_arr)
        x_inv_arr = np.linalg.inv(x_mod_arr)
        vector_param = np.dot(x_inv_arr, y_mod_arr)
        # print('cccc')
        for m in range(0, len(vector_param)):
            back_level_par[m] = vector_param[m]
        # print(back_level_par)
        # Correct the initial frame from the changing background
        initial_frame_data = np.copy(frame_data)


        for i_ind in range(frame_mid_x * 2):
            for j_ind in range(frame_mid_y * 2):
                if frame_data[i_ind][j_ind] != 0:
                    frame_data[i_ind][j_ind] -= back_level_par[1] * i_ind + back_level_par[
                        2] * j_ind + back_level_par[3] * i_ind ** 2 + back_level_par[4] * i_ind * j_ind + \
                                                back_level_par[5] * j_ind ** 2
    back_level = astropy.stats.sigma_clipped_stats(frame_data, sigma=1.5)[0]
    #plt.figure(2)
    #im = plt.imshow(frame_data)
    #plt.colorbar(im)
    #c_min = back_level
    #plt.clim([c_min, c_min + math.sqrt(initial_back_level + readout_noise ** 2) * 8])
    #plt.show()
    #plt.close()
    frame_max = np.max(frame_data)
    #print(frame_max)
    if back_level <=500:
        frame_max = 1000
    # print(frame_mid_y, frame_mid_x)
    found_stars = 0
    search_frame = np.copy(frame_data)

    # Filling the border with background
    for x_ind in range(10):
        for y_ind in range(frame_mid_y * 2):
            search_frame[x_ind][y_ind] = back_level
    for y_ind in range(10):
        for x_ind in range(frame_mid_y * 2):
            search_frame[x_ind][y_ind] = back_level
    for x_ind in range(frame_mid_x * 2 - 10, frame_mid_x * 2):
        for y_ind in range(frame_mid_y * 2):
            search_frame[x_ind][y_ind] = back_level
    for y_ind in range(frame_mid_x * 2 - 10, frame_mid_x * 2):
        for x_ind in range(frame_mid_y * 2):
            search_frame[x_ind][y_ind] = back_level

    # Array of bright star coordinates
    bright_star_coord_x = []
    bright_star_coord_y = []
    # bright_star_peak = []
    bright_star_photo = []
    iter_num = 0
    print('Picking bright stars')
    while found_stars < num_of_bright:
        search_level = (frame_max - back_level) * (100 - iter_num - 3) / 100 + back_level
        #print(search_level)
        if search_level <=0:
            break
        if search_level <= back_level * 1:
            break
        poss_star_coord_x, poss_star_coord_y = np.where(search_frame > search_level)

        def_star_coord_x, def_star_coord_y = [10], [10]
        if len(poss_star_coord_x) == 1:
            def_star_coord_x, def_star_coord_y = [poss_star_coord_x[0]], [poss_star_coord_y[0]]
            def_star_coord_x[0], def_star_coord_y[0] = fast_obj_centroid(poss_star_coord_x[0], poss_star_coord_y[0],
                                                                         back_level, frame_data)
            def_star_photo = fast_photo(def_star_coord_x[0], def_star_coord_y[0],
                                        back_level, frame_data)
            # def_star_peak = obj_peak(def_star_coord_x[0], def_star_coord_y[0], search_frame)
            bright_star_coord_x = np.append(bright_star_coord_x, def_star_coord_x[0])
            bright_star_coord_y = np.append(bright_star_coord_y, def_star_coord_y[0])
            # bright_star_peak = np.append(bright_star_peak, def_star_peak)
            bright_star_photo = np.append(bright_star_photo, def_star_photo)
            # print(def_star_coord_x, def_star_coord_y)
            found_stars += 1
        elif len(poss_star_coord_x) > 1:
            def_star_coord_x, def_star_coord_y = [0], [0]

            def_star_coord_x[0], def_star_coord_y[0] = fast_obj_centroid(poss_star_coord_x[0], poss_star_coord_y[0],
                                                                         back_level, frame_data)
            def_star_photo = fast_photo(def_star_coord_x[0], def_star_coord_y[0],
                                        back_level, frame_data)
            # def_star_peak = obj_peak(def_star_coord_x[0], def_star_coord_y[0], search_frame)
            bright_star_coord_x = np.append(bright_star_coord_x, def_star_coord_x[0])
            bright_star_coord_y = np.append(bright_star_coord_y, def_star_coord_y[0])
            # bright_star_peak = np.append(bright_star_peak, def_star_peak)
            bright_star_photo = np.append(bright_star_photo, def_star_photo)
            found_stars += 1
            if found_stars >= num_of_bright:
                break
            for poss_ind in range(1, len(poss_star_coord_y)):
                distance_check = 0  # From how many definite stars our object is far away
                for def_ind in range(len(def_star_coord_x)):
                    if abs(poss_star_coord_x[poss_ind] - def_star_coord_x[def_ind]) > 20 or abs(
                            poss_star_coord_y[poss_ind] - def_star_coord_y[def_ind]) > 20:
                        distance_check += 1
                if distance_check == len(def_star_coord_x):
                    inter_star_coord_x, inter_star_coord_y = fast_obj_centroid(poss_star_coord_x[poss_ind],
                                                                               poss_star_coord_y[poss_ind],
                                                                               back_level, frame_data)
                    def_star_photo = fast_photo(inter_star_coord_x, inter_star_coord_y,
                                                back_level, frame_data)
                    # inter_star_peak = obj_peak(inter_star_coord_x, inter_star_coord_y, search_frame)
                    def_star_coord_x = np.append(def_star_coord_x, inter_star_coord_x)
                    def_star_coord_y = np.append(def_star_coord_y, inter_star_coord_y)
                    bright_star_coord_x = np.append(bright_star_coord_x, inter_star_coord_x)
                    bright_star_coord_y = np.append(bright_star_coord_y, inter_star_coord_y)
                    # bright_star_peak = np.append(bright_star_peak, inter_star_peak)
                    bright_star_photo = np.append(bright_star_photo, def_star_photo)
                    found_stars += 1
                    # if found_stars >= num_of_bright:
                    #    break
        if found_stars >= num_of_bright:
            #print('Too many stars')
            break
        iter_num += 1
        #print('aaaa')
        if def_star_coord_x[0] == 10:
            continue

        for def_ind in range(len(def_star_coord_x)):
            fill_area(def_star_coord_x[def_ind], def_star_coord_y[def_ind], search_frame, back_level)
    bright_star_series = {'x': bright_star_coord_x, 'y': bright_star_coord_y, 'photo': bright_star_photo}
    bright_star_df = pd.DataFrame(bright_star_series)
    bright_star_df = bright_star_df.sort_values(by='photo', ascending=False)
    bright_star_df = bright_star_df.reset_index(drop=True)
    # print(bright_star_df)
    return bright_star_df['y'], bright_star_df['x'], frame_data


def obj_peak(input_x, input_y, radius, frame_data):
    x0 = int(round(input_x))
    y0 = int(round(input_y))
    for iter in range(3):
        x_top = int(x0 + radius + 1)
        x_bottom = int(x0 - radius)
        y_top = int(y0 + radius + 1)
        y_bottom = int(y0 - radius)
        if x_top > frame_data.shape[0]:
            x_top = frame_data.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > frame_data.shape[1]:
            y_top = frame_data.shape[1]
        if y_bottom < 0:
            y_bottom = 0
        # Calculating centroid
        max = 0
        max_x, max_y = 0, 0
        #background_arr = []
        # max = np.max(frame_data[x_bottom:x_top][y_bottom:y_top])
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if frame_data[i][j] > 0:
                    if (rad_val <= radius) and frame_data[i][j] > max:
                        max = frame_data[i][j]
                        max_x = i
                        max_y = j
                    #if rad_va
                    #background_arr = np.append(background_arr, ground[i][j])
        #background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=2)[0]
    return max_x, max_y


def back_level(input_x, input_y, radius, annul, frame_data):           # Find background signal near the position
    x0 = int(round(input_x))
    y0 = int(round(input_y))
    x_top = int(x0 + radius + annul + 1)
    x_bottom = int(x0 - radius - annul)
    y_top = int(y0 + radius + annul + 1)
    y_bottom = int(y0 - radius - annul)
    if x_top > frame_data.shape[0]:
        x_top = frame_data.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > frame_data.shape[1]:
        y_top = frame_data.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    # Calculating centroid
    background_arr = []
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if frame_data[i][j] > 0:
                if (rad_val > radius):
                    background_arr = np.append(background_arr, frame_data[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=2)[0]
    return background_per_pixel


def fast_photo(input_x, input_y, back_level, frame_data):
    x0 = int(round(input_x))
    y0 = int(round(input_y))
    radius = 4

    x_top = int(x0 + radius + 1)
    x_bottom = int(x0 - radius)
    y_top = int(y0 + radius + 1)
    y_bottom = int(y0 - radius)
    if x_top > frame_data.shape[0]:
        x_top = frame_data.shape[0]
    if x_bottom < 0:
        x_bottom = 0
    if y_top > frame_data.shape[1]:
        y_top = frame_data.shape[1]
    if y_bottom < 0:
        y_bottom = 0
    # Calculating centroid
    source_signal = 0
    # max = np.max(frame_data[x_bottom:x_top][y_bottom:y_top])
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if frame_data[i][j] > 0:
                if (rad_val <= radius):
                    source_signal += frame_data[i][j] - back_level

    return source_signal

# Creating a mask of all stars (a star template)
def check_mask_center(check_x, check_y, radius, gap, annul, ground):
    mask_discard = 0
    x0 = int(round(check_x))
    y0 = int(round(check_y))
    radius = int(radius)
    gap = int(gap)
    annul = int(annul)
    x_top = x0 + annul + 1
    x_bottom = x0 - annul
    y_top = y0 + annul + 1
    y_bottom = y0 - annul
    background_arr = []
    if x_top > ground.shape[0]:
        x_top = ground.shape[0]
        mask_discard = 1
    if x_bottom < 0:
        x_bottom = 0
        mask_discard = 1
    if y_top > ground.shape[1]:
        y_top = ground.shape[1]
        mask_discard = 1
    if y_bottom < 0:
        y_bottom = 0
        mask_discard = 1
    # Calculating background
    for i in range(x_bottom, x_top):
        for j in range(y_bottom, y_top):
            rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
            if ground[i][j] > 0:
                if rad_val > radius + gap and rad_val <= annul + 1:
                    background_arr = np.append(background_arr, ground[i][j])
    background_per_pixel = astropy.stats.sigma_clipped_stats(background_arr, sigma=2)[0]
    # Finding centre
    for iter in range(3):
        x0 = int(round(check_x))
        y0 = int(round(check_y))
        radius = int(radius)
        gap = int(gap)
        annul = int(annul)

        x_top = x0 + (radius + gap) + 1
        x_bottom = x0 - (radius + gap)
        y_top = y0 + (radius + gap) + 1
        y_bottom = y0 - (radius + gap)
        if x_top > ground.shape[0]:
            x_top = ground.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > ground.shape[1]:
            y_top = ground.shape[1]
        if y_bottom < 0:
            y_bottom = 0
        # Calculating centroid
        source_signal = 0
        x_signal = 0
        y_signal = 0
        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                rad_val = math.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))
                if ground[i][j] > 0:
                    if (rad_val <= radius * 3 and rad_val <= radius + gap):
                        source_signal += ground[i][j] - background_per_pixel
                        x_signal += (ground[i][j] - background_per_pixel) * i
                        y_signal += (ground[i][j] - background_per_pixel) * j
        if source_signal <= 0.1:
            x_centre = x0
            y_centre = y0
            break
        else:
            x_centre = x_signal / source_signal
            y_centre = y_signal / source_signal

        check_x = x_centre
        check_y = y_centre
    mask_arr = zeros((20, 20))
    if mask_discard == 0:
        # Making a mask
        x_top = x0 + 10
        x_bottom = x0 - 10
        y_top = y0 + 10
        y_bottom = y0 - 10
        if x_top > ground.shape[0]:
            x_top = ground.shape[0]
        if x_bottom < 0:
            x_bottom = 0
        if y_top > ground.shape[1]:
            y_top = ground.shape[1]
        if y_bottom < 0:
            y_bottom = 0

        for i in range(x_bottom, x_top):
            for j in range(y_bottom, y_top):
                mask_arr[i-x_bottom][j-y_bottom] = ground[i][j] - background_per_pixel
    check_max = np.max(ground[x_bottom:x_top, y_bottom:y_top] - background_per_pixel)
    mask_arr = mask_arr / check_max
    return x_centre, y_centre, mask_arr

# Erasing a star with a given peak and a given form(through mask)

#def star_eraser()
