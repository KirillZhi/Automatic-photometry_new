from math import sin, cos, asin, acos, atan2, tan, sqrt, isnan, log10

import numpy as np
import pandas as pd
from numpy import zeros, transpose, dot, arange, append, copy, zeros_like, sum, min, array, meshgrid, linspace, absolute
from numpy import sqrt as arr_sqrt
from numpy import max as arr_max
from numpy.linalg import inv as invert_array
import math
from numpy.linalg import inv
import matplotlib.pyplot as plt
from pandas import read_csv
from astropy.io.fits import open, getdata
from astropy.coordinates import SkyCoord, Angle

from itertools import combinations
from astropy.coordinates.angle_utilities import angular_separation, position_angle
import astropy.units as u

rad_conversion = 57.295828
scale = 1939


# scale = 1901
# lon1 = 0
# lat1 = 0
# lon2 = 2
# lat2 = 2
# print(angular_separation(lon1, lat1, lon2, lat2))

def min_fun(exp_x, exp_y, theoretical_x, theoretical_y, parameters):  # Function we minimize
    init_fun = 0
    f_fun = 0
    g_fun = 0
    for fun_ind in range(len(exp_y)):
        f_fun = exp_x[fun_ind] - parameters[3] * theoretical_x[fun_ind] * cos(parameters[2]) - parameters[3] * \
                theoretical_y[fun_ind] * sin(
            parameters[2]) - parameters[0] * cos(parameters[2]) - parameters[1] * sin(parameters[2])
        g_fun = exp_y[fun_ind] + parameters[3] * theoretical_x[fun_ind] * sin(parameters[2]) - parameters[3] * \
                theoretical_y[fun_ind] * cos(
            parameters[2]) + parameters[0] * sin(parameters[2]) - parameters[1] * cos(parameters[2])
        init_fun += f_fun ** 2 + g_fun ** 2
    del f_fun, g_fun
    return init_fun

def adv_min_fun(exp_x, exp_y, theoretical_x, theoretical_y, parameters):  # Function we minimize
    init_fun = 0
    f_fun = 0
    g_fun = 0
    for fun_ind in range(len(exp_y)):
        f_fun = exp_x[fun_ind] - parameters[3] * theoretical_x[fun_ind] * cos(parameters[2]) - parameters[3] * \
                theoretical_y[fun_ind] * sin(
            parameters[2]) - parameters[0] * cos(parameters[2]) - parameters[1] * sin(parameters[2]) - parameters[4]
        g_fun = exp_y[fun_ind] + parameters[3] * theoretical_x[fun_ind] * sin(parameters[2]) - parameters[3] * \
                theoretical_y[fun_ind] * cos(
            parameters[2]) + parameters[0] * sin(parameters[2]) - parameters[1] * cos(parameters[2]) - parameters[5]
        init_fun += f_fun ** 2 + g_fun ** 2
    del f_fun, g_fun
    return init_fun

def adv_fun_min_search(input_x, input_y, cel_x, cel_y, init_teta):  # Search for a minimum in coordinate deviation
    min = 100000000
    teta_min = init_teta * rad_conversion
    teta_step = [60, 30, 15, 7.5, 3.75, 1.875, 0.938, 0.469, 0.234, 0.117, 0.0585, 0.02925, 0.014625, 0.00365625]
    # teta_step = [30, 15, 7.5, 3.75, 1.875, 0.938, 0.469, 0.234, 0.117, 0.0585, 0.02925, 0.014625, 0.00365625]
    #teta_step = [0.25]
    frame_inter_parameters = zeros(7)

    export_parameters = copy(frame_inter_parameters)
    if init_teta == 1000:

        for teta_step_ind in range(len(teta_step)):

            if teta_step_ind == 0:

                angle_step = teta_step[teta_step_ind]
                teta_lower = 0
                teta_upper = 360 + angle_step
                teta_steps = int((teta_upper - teta_lower) / angle_step)

            else:
                angle_step = teta_step[teta_step_ind - 1]
                teta_lower = teta_min - angle_step
                teta_upper = teta_min + angle_step * 1.2
                teta_steps = int((teta_upper - teta_lower) / teta_step[teta_step_ind])
                angle_step = teta_step[teta_step_ind]

            for step_ind in range(0, teta_steps):
                teta_value = (teta_lower + angle_step * step_ind) / rad_conversion
                # Determining scale
                x_cel_sum = sum(cel_x)
                y_cel_sum = sum(cel_y)
                x_transform = input_x * cos(teta_value) - input_y * sin(teta_value)
                y_transform = input_x * sin(teta_value) + input_y * cos(teta_value)
                x_transform_sum = sum(x_transform)
                y_transform_sum = sum(y_transform)
                anchor_num = len(cel_x)
                A_arr = np.array([
                    [anchor_num, 0, x_cel_sum],
                    [0, anchor_num, y_cel_sum],
                    [x_cel_sum, y_cel_sum, sum(cel_x ** 2 + cel_y ** 2)]], dtype=float)
                #print(A_arr, np.linalg.det(A_arr))
                #print(A_arr.shape)
                Y_arr = np.array([[x_transform_sum],
                                  [y_transform_sum],
                                  [sum((input_x * cel_x + input_y * cel_y) * cos(
                                      teta_value) + (input_x * cel_y - input_y * cel_x) * sin(teta_value))]])
                A_inv = np.linalg.inv(A_arr)
                X_arr = np.dot(A_inv, Y_arr)

                frame_inter_parameters[0] = X_arr[0]
                frame_inter_parameters[1] = X_arr[1]
                frame_inter_parameters[2] = teta_value
                frame_inter_parameters[3] = X_arr[2]



                # Determining other parameters
                function_value = adv_min_fun(input_x, input_y, cel_x, cel_y, frame_inter_parameters)

                if function_value < min:
                    min = function_value
                    teta_min = teta_value * rad_conversion
                    export_parameters = copy(frame_inter_parameters)
    else:
        teta_step = [3.75, 1.875, 0.938, 0.469, 0.234, 0.117, 0.0585, 0.02925, 0.014625, 0.00365625]

        for teta_step_ind in range(len(teta_step)):
            angle_step = teta_step[teta_step_ind]
            teta_lower = teta_min - angle_step * 1.5
            teta_upper = teta_min + angle_step * 1.5
            teta_steps = int((teta_upper - teta_lower) / teta_step[teta_step_ind])

            for step_ind in range(0, teta_steps):
                teta_value = (teta_lower + angle_step * step_ind) / rad_conversion

                x_cel_sum = sum(cel_x)
                y_cel_sum = sum(cel_y)
                x_transform = input_x * cos(teta_value) - input_y * sin(teta_value)
                y_transform = input_x * sin(teta_value) + input_y * cos(teta_value)
                x_transform_sum = sum(x_transform)
                y_transform_sum = sum(y_transform)
                anchor_num = len(cel_x)
                A_arr = np.array([
                    [anchor_num, 0,  x_cel_sum],
                    [0, anchor_num,  y_cel_sum],
                    [x_cel_sum, y_cel_sum, sum(cel_x ** 2 + cel_y ** 2)]], dtype=float)
                #print(A_arr, np.linalg.det(A_arr))
                #print(A_arr.shape)
                Y_arr = np.array([[x_transform_sum],
                         [y_transform_sum],
                         [sum((input_x * cel_x + input_y * cel_y) * cos(
                             teta_value) + (input_x * cel_y - input_y * cel_x) * sin(teta_value))]])
                A_inv = np.linalg.inv(A_arr)
                X_arr = np.dot(A_inv, Y_arr)

                frame_inter_parameters[0] = X_arr[0]
                frame_inter_parameters[1] = X_arr[1]
                frame_inter_parameters[2] = teta_value
                frame_inter_parameters[3] = X_arr[2]

                # Determining other parameters
                function_value = adv_min_fun(input_x, input_y, cel_x, cel_y, frame_inter_parameters)
                if function_value < min:
                    min = function_value
                    teta_min = teta_value * rad_conversion
                    export_parameters = copy(frame_inter_parameters)

    del function_value, teta_value, angle_step, teta_upper, teta_lower, teta_steps
    del teta_step_ind, A_arr, A_inv, Y_arr, X_arr, x_transform, y_transform
    return min, export_parameters


def fun_min_search(input_x, input_y, cel_x, cel_y, init_teta):  # Search for a minimum in coordinate deviation
    min = 100000000
    teta_min = init_teta * rad_conversion
    teta_step = [60, 30, 15, 7.5, 3.75, 1.875, 0.938, 0.469, 0.234, 0.117, 0.0585, 0.02925, 0.014625, 0.00365625]
    # teta_step = [30, 15, 7.5, 3.75, 1.875, 0.938, 0.469, 0.234, 0.117, 0.0585, 0.02925, 0.014625, 0.00365625]
    #teta_step = [0.25]
    frame_inter_parameters = zeros(12)

    export_parameters = copy(frame_inter_parameters)
    if init_teta == 1000:

        for teta_step_ind in range(len(teta_step)):

            if teta_step_ind == 0:

                angle_step = teta_step[teta_step_ind]
                teta_lower = 0
                teta_upper = 360 + angle_step
                teta_steps = int((teta_upper - teta_lower) / angle_step)

            else:
                angle_step = teta_step[teta_step_ind - 1]
                teta_lower = teta_min - angle_step
                teta_upper = teta_min + angle_step * 1.2
                teta_steps = int((teta_upper - teta_lower) / teta_step[teta_step_ind])
                angle_step = teta_step[teta_step_ind]

            for step_ind in range(0, teta_steps):
                teta_value = (teta_lower + angle_step * step_ind) / rad_conversion
                # Determining scale
                x_cel_sum = sum(cel_x)
                y_cel_sum = sum(cel_y)
                x_transform = input_x * cos(teta_value) - input_y * sin(teta_value)
                y_transform = input_x * sin(teta_value) + input_y * cos(teta_value)
                x_transform_sum = sum(x_transform)
                y_transform_sum = sum(y_transform)
                scale_left = sum(
                    cel_x ** 2 + cel_y ** 2 - cel_x / len(input_x) * x_cel_sum - cel_y / len(input_x) * y_cel_sum)
                scale_right = sum(
                    (input_x * cel_x + input_y * cel_y) * cos(teta_value) + (input_x * cel_y - cel_x * input_y) * sin(
                        teta_value) - cel_x * x_transform_sum / len(input_x) - cel_y * y_transform_sum / len(input_x))
                frame_inter_parameters[3] = scale_right / scale_left
                # Determining other parameters
                frame_inter_parameters[0] = (x_transform_sum - frame_inter_parameters[3] * x_cel_sum) / len(input_x)
                frame_inter_parameters[1] = (y_transform_sum - frame_inter_parameters[3] * y_cel_sum) / len(input_x)
                frame_inter_parameters[2] = teta_value
                function_value = min_fun(input_x, input_y, cel_x, cel_y, frame_inter_parameters)

                if function_value < min:
                    min = function_value
                    teta_min = teta_value * rad_conversion
                    export_parameters = copy(frame_inter_parameters)
    else:
        teta_step = [3.75, 1.875, 0.938, 0.469, 0.234, 0.117, 0.0585, 0.02925, 0.014625, 0.00365625]

        for teta_step_ind in range(len(teta_step)):
            angle_step = teta_step[teta_step_ind]
            teta_lower = teta_min - angle_step * 1.5
            teta_upper = teta_min + angle_step * 1.5
            teta_steps = int((teta_upper - teta_lower) / teta_step[teta_step_ind])

            for step_ind in range(0, teta_steps):
                teta_value = (teta_lower + angle_step * step_ind) / rad_conversion

                # Determining scale
                x_cel_sum = sum(cel_x)
                y_cel_sum = sum(cel_y)
                x_transform = input_x * cos(teta_value) - input_y * sin(teta_value)
                y_transform = input_x * sin(teta_value) + input_y * cos(teta_value)
                x_transform_sum = sum(x_transform)
                y_transform_sum = sum(y_transform)
                scale_left = sum(
                    cel_x ** 2 + cel_y ** 2 - cel_x / len(input_x) * x_cel_sum - cel_y / len(input_x) * y_cel_sum)
                scale_right = sum(
                    (input_x * cel_x + input_y * cel_y) * cos(teta_value) + (input_x * cel_y - cel_x * input_y) * sin(
                        teta_value) - cel_x * x_transform_sum / len(input_x) - cel_y * y_transform_sum / len(input_x))
                frame_inter_parameters[3] = scale_right / scale_left
                # Determining other parameters
                frame_inter_parameters[0] = (x_transform_sum - frame_inter_parameters[3] * x_cel_sum) / len(input_x)
                frame_inter_parameters[1] = (y_transform_sum - frame_inter_parameters[3] * y_cel_sum) / len(input_x)
                frame_inter_parameters[2] = teta_value
                function_value = min_fun(input_x, input_y, cel_x, cel_y, frame_inter_parameters)
                if function_value < min:
                    min = function_value
                    teta_min = teta_value * rad_conversion
                    export_parameters = copy(frame_inter_parameters)

    del function_value, teta_value, angle_step, teta_upper, teta_lower, teta_steps
    del teta_step_ind, scale_left, scale_right, x_transform, y_transform
    return min, export_parameters


def axis_fit_new(input_x, input_y, mid_ra, mid_dec, input_ra, input_dec, input_parameters):
    # input_x, input_y - x and y coordinate on the frame (distorted)
    # frame_path - path to a frame
    # input_ra, input_dec - corresponding celestial coordinates
    # axis_direction - check if a frame is mirrored in some way
    axis_direction = input_parameters[4]
    x_cel = zeros_like(input_x)
    y_cel = zeros_like(input_y)

    # Making a celestial projection array
    for cel_ind in range(len(x_cel)):
        x_cel[cel_ind] = rad_conversion * (
                cos(input_dec[cel_ind] / rad_conversion) * sin(mid_dec / rad_conversion) * cos(
            (input_ra[cel_ind] - mid_ra) / rad_conversion) - sin(input_dec[cel_ind] / rad_conversion) * cos(
            mid_dec / rad_conversion))
        y_cel[cel_ind] = rad_conversion * (
                sin((input_ra[cel_ind] - mid_ra) / rad_conversion) * cos(input_dec[cel_ind] / rad_conversion))

    min_boundary = 2 * len(input_x) * 15 ** 2

    min_function = 100000000

    min_axis = 6
    if axis_direction == 5:
        for axis_direction in range(2):

            x_cel_inter = copy(x_cel)
            y_cel_inter = copy(y_cel)
            if axis_direction == 1:
                x_cel_inter = x_cel * (-1)
            min_value, inter_parameters = fun_min_search(input_x, input_y, x_cel_inter, y_cel_inter, 1000)

            check_value = min_fun(input_x, input_y, x_cel_inter, y_cel_inter, inter_parameters)

            if min_value <= min_function:
                min_function = min_value

                min_axis = axis_direction
                x_cel_final = copy(x_cel_inter)
                y_cel_final = copy(y_cel_inter)
                output_parameters = inter_parameters
                x_cel_final = copy(x_cel_inter)
                y_cel_final = copy(y_cel_inter)

            if min_value <= min_boundary:
                min_function = min_value
                min_axis = axis_direction
                output_parameters = inter_parameters
                x_cel_final = copy(x_cel_inter)
                y_cel_final = copy(y_cel_inter)

                break
    else:
        x_cel_inter = copy(x_cel)
        y_cel_inter = copy(y_cel)
        if axis_direction == 1:
            x_cel_inter = x_cel * (-1)
        min_value, output_parameters = fun_min_search(input_x, input_y, x_cel_inter, y_cel_inter, input_parameters[2])
        x_cel_final = copy(x_cel_inter)
        y_cel_final = copy(y_cel_inter)
        min_axis = axis_direction
    # Estimating final parameters

    output_parameters[4] = min_axis
    x_cel_sum = sum(x_cel_final)
    y_cel_sum = sum(y_cel_final)
    x_transform = input_x * cos(output_parameters[2]) - input_y * sin(output_parameters[2])
    y_transform = input_x * sin(output_parameters[2]) + input_y * cos(output_parameters[2])
    x_transform_sum = sum(x_transform)
    y_transform_sum = sum(y_transform)
    scale_left = sum(x_cel_final ** 2 + y_cel_final ** 2 - x_cel_final / len(input_x) * x_cel_sum - y_cel_final / len(
        input_x) * y_cel_sum)
    scale_right = sum((input_x * x_cel_final + input_y * y_cel_final) * cos(output_parameters[2]) + (
                input_x * y_cel_final - x_cel_final * input_y) * sin(
        output_parameters[2]) - x_cel_final * x_transform_sum / len(input_x) - y_cel_final * y_transform_sum / len(
        input_x))
    output_parameters[3] = scale_right / scale_left

    # Determining other parameters
    output_parameters[0] = (x_transform_sum - output_parameters[3] * x_cel_sum) / len(input_x)
    output_parameters[1] = (y_transform_sum - output_parameters[3] * y_cel_sum) / len(input_x)

    min_value = min_fun(input_x, input_y, x_cel_final, y_cel_final, output_parameters)
    print("Sum of squared deviations: " + "{}".format(min_value))
    del x_cel, y_cel
    del x_cel_final, x_cel_inter, y_cel_inter, y_cel_final, min_axis, min_function

    return output_parameters, min_value


def axis_fit_complete(input_x, input_y, mid_ra, mid_dec, input_ra, input_dec,
                      input_parameters):  # Fitting distortion and rotation
    # input_x, input_y - x and y coordinate on the frame (distorted)
    # frame_path - path to a frame
    # input_ra, input_dec - corresponding celestial coordinates
    # axis_direction - check if a frame is mirrored in some way
    x_cel = zeros_like(input_x)
    y_cel = zeros_like(input_y)

    # Making a celestial projection array
    for cel_ind in range(len(x_cel)):
        x_cel[cel_ind] = rad_conversion * (
                cos(input_dec[cel_ind] / rad_conversion) * sin(mid_dec / rad_conversion) * cos(
            (input_ra[cel_ind] - mid_ra) / rad_conversion) - sin(input_dec[cel_ind] / rad_conversion) * cos(
            mid_dec / rad_conversion))
        y_cel[cel_ind] = rad_conversion * (
                sin((input_ra[cel_ind] - mid_ra) / rad_conversion) * cos(input_dec[cel_ind] / rad_conversion))
    x_cel_use = copy(x_cel)
    y_cel_use = copy(y_cel)
    non_distorted_parameters, non_distorted_sum = axis_fit_nondistorted(input_x, input_y, x_cel_use, y_cel_use, input_parameters)
    x_cel_use = copy(x_cel)
    y_cel_use = copy(y_cel)
    inter_parameters = zeros(12)
    for inter_ind in range(5):
        inter_parameters[inter_ind] = non_distorted_parameters[inter_ind]
    if non_distorted_sum / len(input_x) / 2 > 0.25:

        non_distorted_x, non_distorted_y = fit_use_nondistorted(x_cel_use, y_cel_use, non_distorted_parameters)
        x_cel_use = copy(x_cel)
        y_cel_use = copy(y_cel)

        # Estimating parameters of the distortion

        check_value, inter_parameters = axis_fit_distortion(input_x, input_y, non_distorted_x, non_distorted_y,
                                                            inter_parameters)
        print("Sum of squared deviations: " + "{}".format(check_value))
        if check_value <= non_distorted_sum:
            if check_value / len(input_x) / 2 > 0.25:
                min_sum, output_parameters = frame_centre_conjugate_grad(input_x, input_y, non_distorted_x, non_distorted_y,
                                                                         inter_parameters)
                print("Sum of squared deviations: " + "{}".format(min_sum))
            else:
                min_sum, output_parameters = check_value, inter_parameters

        else:
            inter_parameters = zeros(12)
            for inter_ind in range(5):
                inter_parameters[inter_ind] = non_distorted_parameters[inter_ind]
            min_sum, output_parameters = non_distorted_sum, inter_parameters
    else:
        min_sum, output_parameters = non_distorted_sum, inter_parameters
    return output_parameters, min_sum

def adv_axis_fit_complete(input_x, input_y, mid_ra, mid_dec, input_ra, input_dec,
                      input_parameters):  # Fitting distortion and rotation
    # input_x, input_y - x and y coordinate on the frame (distorted)
    # frame_path - path to a frame
    # input_ra, input_dec - corresponding celestial coordinates
    # axis_direction - check if a frame is mirrored in some way
    x_cel = zeros_like(input_x)
    y_cel = zeros_like(input_y)

    # Making a celestial projection array
    for cel_ind in range(len(x_cel)):
        x_cel[cel_ind] = rad_conversion * (
                cos(input_dec[cel_ind] / rad_conversion) * sin(mid_dec / rad_conversion) * cos(
            (input_ra[cel_ind] - mid_ra) / rad_conversion) - sin(input_dec[cel_ind] / rad_conversion) * cos(
            mid_dec / rad_conversion))
        y_cel[cel_ind] = rad_conversion * (
                sin((input_ra[cel_ind] - mid_ra) / rad_conversion) * cos(input_dec[cel_ind] / rad_conversion))
    x_cel_use = copy(x_cel)
    y_cel_use = copy(y_cel)
    non_distorted_parameters, non_distorted_sum = adv_axis_fit_nondistorted(input_x, input_y, x_cel_use, y_cel_use, input_parameters)
    print("Sum of squared deviations: " + "{}".format(non_distorted_sum))
    x_cel_use = copy(x_cel)
    y_cel_use = copy(y_cel)
    inter_parameters = zeros(7)
    for inter_ind in range(7):
        inter_parameters[inter_ind] = non_distorted_parameters[inter_ind]
    min_sum, output_parameters = non_distorted_sum, inter_parameters

    return output_parameters, min_sum

def adv_axis_fit_use(input_ra, input_dec, mid_ra, mid_dec, frame_parameters):
    # Making a celestial projection array
    x_cel = 0
    y_cel = 0
    x_cel = frame_parameters[3] * rad_conversion * (
                cos(input_dec / rad_conversion) * sin(mid_dec / rad_conversion) * cos(
            (input_ra - mid_ra) / rad_conversion) - sin(input_dec / rad_conversion) * cos(
            mid_dec / rad_conversion))
    y_cel = frame_parameters[3] * rad_conversion * (
                sin((input_ra - mid_ra) / rad_conversion) * cos(input_dec / rad_conversion))

    # Finding coordinates
    if frame_parameters[4] == 1:
        x_cel = x_cel * (-1)
    elif frame_parameters[4] == 2:
        y_cel = y_cel * (-1)
    elif frame_parameters[4] == 3:
        y_cel = y_cel * (-1)
        x_cel = x_cel * (-1)

    x_theor = x_cel * cos(frame_parameters[2]) + y_cel * sin(frame_parameters[2])
    y_theor = -x_cel * sin(frame_parameters[2]) + y_cel * cos(frame_parameters[2])
    x_output = x_theor + frame_parameters[0] * cos(frame_parameters[2]) + frame_parameters[1] * sin(frame_parameters[2]) + frame_parameters[5]
    y_output = y_theor - frame_parameters[0] * sin(frame_parameters[2]) + frame_parameters[1] * cos(frame_parameters[2]) + frame_parameters[6]

    del x_theor, y_theor, x_cel, y_cel, mid_ra, mid_dec
    return int(round(x_output)), int(round(y_output))

def adv_axis_fit_nondistorted(input_x, input_y, cel_x, cel_y, input_parameters):  # Fitting only rotation, shift and scale
    min_boundary = 2 * len(input_x) * 15 ** 2
    axis_direction = input_parameters[4]
    # input_parameters[2] = input_parameters[2] * rad_conversion
    # print(input_parameters[2])
    # print(axis_direction)
    if axis_direction == 5:
        for axis_direction in range(4):
            x_cel_inter = copy(cel_x)
            y_cel_inter = copy(cel_y)
            if axis_direction == 1:
                x_cel_inter = cel_x * (-1)
            elif axis_direction == 2:
                y_cel_inter = cel_y * (-1)
            elif axis_direction == 3:
                y_cel_inter = cel_y * (-1)
                x_cel_inter = cel_x * (-1)
            min_value, inter_parameters = adv_fun_min_search(input_x, input_y, x_cel_inter, y_cel_inter, 1000)

            if min_value <= min_boundary:
                min_function = min_value
                min_axis = axis_direction
                x_cel_final = copy(x_cel_inter)
                y_cel_final = copy(y_cel_inter)
                inter_parameters[4] = min_axis
                break
    else:
        x_cel_inter = copy(cel_x)
        y_cel_inter = copy(cel_y)
        if axis_direction == 1:
            x_cel_inter = cel_x * (-1)
        elif axis_direction == 2:
            y_cel_inter = cel_y * (-1)
        elif axis_direction == 3:
            y_cel_inter = cel_y * (-1)
            x_cel_inter = cel_x * (-1)
        # print(input_parameters[2])9
        min_value, inter_parameters = adv_fun_min_search(input_x, input_y, x_cel_inter, y_cel_inter, input_parameters[2])
        x_cel_final = copy(x_cel_inter)
        y_cel_final = copy(y_cel_inter)
        min_axis = axis_direction
        inter_parameters[4] = min_axis
    print("Sum of squared deviations: " + "{}".format(min_value))
    # print(inter_parameters)
    return inter_parameters, min_value

def fit_use_nondistorted(input_cel_x, input_cel_y, input_parameters):  # Use fit parameters (nondistorted)
    cel_x = copy(input_cel_x)
    cel_y = copy(input_cel_y)
    if input_parameters[4] == 1:
        cel_x = input_cel_x * (-1)
    # elif input_parameters[4] == 2:
    #    cel_y = input_cel_y * (-1)

    x_theor = cel_x * cos(input_parameters[2]) + cel_y * sin(input_parameters[2])
    y_theor = -cel_x * sin(input_parameters[2]) + cel_y * cos(input_parameters[2])

    x_theor = x_theor * input_parameters[3]
    y_theor = y_theor * input_parameters[3]

    x_output = array(
        x_theor + input_parameters[0] * cos(input_parameters[2]) + input_parameters[1] * sin(input_parameters[2]))
    y_output = array(
        y_theor - input_parameters[0] * sin(input_parameters[2]) + input_parameters[1] * cos(input_parameters[2]))

    return x_output, y_output


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


# Finding centre of cadre through conjugate gradients

def frame_centre_conjugate_grad(input_x, input_y, supposed_x, supposed_y, input_parameters):
    check_value = distortion_min_sum(input_x, input_y, supposed_x, supposed_y, input_parameters)
    pix_dev = 0.2
    prev_value = check_value
    init_ind = 0
    # prev_descent
    output_parameters = copy(input_parameters)
    while check_value >= pix_dev ** 2 * 2 * len(input_x):
        # Calculating gradient of the minimal sum
        # We calculate gradients through additives to the undistorted axes
        # Derivative through x axis
        min_dev_x0, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x - 0.1, supposed_y,
                                                                 output_parameters)
        min_dev_x1, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y,
                                                                 output_parameters)
        min_dev_x2, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x + 0.1, supposed_y,
                                                                 output_parameters)
        sum_x_der = (4 * min_dev_x1 - 3 * min_dev_x0 - min_dev_x2) / (2 * 0.1)
        sum_x_sec_der = (min_dev_x2 - 2 * min_dev_x1 + min_dev_x0) / (0.1 ** 2)
        # Derivative through y axis
        min_dev_y0, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y - 0.1,
                                                                 output_parameters)
        min_dev_y1, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y,
                                                                 output_parameters)
        min_dev_y2, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y + 0.1,
                                                                 output_parameters)
        sum_y_der = (4 * min_dev_y1 - 3 * min_dev_y0 - min_dev_y2) / (2 * 0.1)
        sum_y_sec_der = (min_dev_y2 - 2 * min_dev_y1 + min_dev_y0) / (0.1 ** 2)

        # Derivative through teta (angle)
        inter_parameters = copy(output_parameters)
        inter_parameters[11] -= 0.01 / rad_conversion
        min_dev_teta0, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y,
                                                                    inter_parameters)

        inter_parameters = copy(output_parameters)
        min_dev_teta1, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y,
                                                                    inter_parameters)

        inter_parameters = copy(output_parameters)
        inter_parameters[11] += 0.01 / rad_conversion
        min_dev_teta2, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y,
                                                                    inter_parameters)

        sum_teta_der = (4 * min_dev_teta1 - 3 * min_dev_teta0 - min_dev_teta2) / (
                    2 * 0.01 / rad_conversion * 3438)  # In arcminutes

        sum_teta_sec_der = (min_dev_teta2 - 2 * min_dev_teta1 + min_dev_teta0) / ((0.01 / rad_conversion * 3438) ** 2)
        # Calculating direction of descent and size of a step
        sum_grad = array([sum_x_der, sum_y_der, sum_teta_der])
        sum_descent = - sum_grad / sqrt(sum(sum_grad ** 2))
        if init_ind != 0:
            w_coef = sum(sum_grad ** 2) / sum(prev_grad ** 2)
            sum_descent += w_coef * prev_descent
        sum_descent = sum_descent / sqrt(sum(sum_descent ** 2))
        prev_descent = copy(sum_descent)
        prev_grad = copy(sum_grad)

        step_multiplier = sqrt(sum(sum_grad ** 2)) / sqrt(
            sum_x_sec_der ** 2 + sum_y_sec_der ** 2 + sum_teta_sec_der ** 2)

        teta_step = 10 / 3438 * sum_descent[2] / abs(sum_descent[2])
        x_step = sum_descent[0] / abs(sum_descent[2]) * 10 / 3438
        y_step = sum_descent[1] / abs(sum_descent[2]) * 10 / 3438
        inter_parameters = copy(output_parameters)
        inter_parameters[9] += x_step
        inter_parameters[10] += y_step
        inter_parameters[11] += teta_step
        descent_check_value, unnecessary_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y,
                                                                          output_parameters)
        if descent_check_value >= prev_value:
            sum_descent = - sum_grad / sqrt(sum(sum_grad ** 2))
        check_value, output_parameters = descent_minimum(input_x, input_y, supposed_x, supposed_y, output_parameters,
                                                         sum_descent, step_multiplier)
        if abs(check_value - prev_value) <= 0.01:
            break
        else:
            prev_value = check_value
        if init_ind == 0:
            init_ind += 1
    return check_value, output_parameters


def descent_minimum(input_x, input_y, supposed_x, supposed_y, input_parameters, descent_direction, step_multiplier):
    min_check = 0

    step_multiplier = 5
    x_step = descent_direction[0] * step_multiplier
    y_step = descent_direction[1] * step_multiplier
    teta_step = descent_direction[2] * step_multiplier / 3438

    inter_parameters = copy(input_parameters)
    prev_min_value = distortion_min_sum(input_x, input_y, supposed_x, supposed_y, inter_parameters)

    while min_check == 0:
        if abs(x_step) < 0.1 or abs(y_step) < 0.1:
            min_check = 1
            break

        inter_parameters[9] += x_step
        inter_parameters[10] += y_step
        inter_parameters[11] += teta_step
        min_value_1, inter_parameters = axis_fit_distortion(input_x, input_y, supposed_x, supposed_y,
                                                            inter_parameters)
        if min_value_1 > prev_min_value:
            inter_parameters[9] -= x_step
            inter_parameters[10] -= y_step
            inter_parameters[11] -= teta_step
            x_step = x_step / 2
            y_step = y_step / 2
            teta_step = teta_step / 2
        else:
            prev_min_value = min_value_1

    # print(inter_parameters)
    return prev_min_value, inter_parameters

def axis_fit_nondistorted(input_x, input_y, cel_x, cel_y, input_parameters):  # Fitting only rotation, shift and scale
    min_boundary = 2 * len(input_x) * 15 ** 2
    axis_direction = input_parameters[4]
    # input_parameters[2] = input_parameters[2] * rad_conversion
    # print(input_parameters[2])
    # print(axis_direction)
    if axis_direction == 5:
        for axis_direction in range(4):
            x_cel_inter = copy(cel_x)
            y_cel_inter = copy(cel_y)
            if axis_direction == 1:
                x_cel_inter = cel_x * (-1)
            elif axis_direction == 2:
                y_cel_inter = cel_y * (-1)
            elif axis_direction == 3:
                y_cel_inter = cel_y * (-1)
                x_cel_inter = cel_x * (-1)
            min_value, inter_parameters = fun_min_search(input_x, input_y, x_cel_inter, y_cel_inter, 1000)

            if min_value <= min_boundary:
                min_function = min_value
                min_axis = axis_direction
                x_cel_final = copy(x_cel_inter)
                y_cel_final = copy(y_cel_inter)
                inter_parameters[4] = min_axis
                break
    else:
        x_cel_inter = copy(cel_x)
        y_cel_inter = copy(cel_y)
        if axis_direction == 1:
            x_cel_inter = cel_x * (-1)
        elif axis_direction == 2:
            y_cel_inter = cel_y * (-1)
        elif axis_direction == 3:
            y_cel_inter = cel_y * (-1)
            x_cel_inter = cel_x * (-1)
        # print(input_parameters[2])
        min_value, inter_parameters = fun_min_search(input_x, input_y, x_cel_inter, y_cel_inter, input_parameters[2])
        x_cel_final = copy(x_cel_inter)
        y_cel_final = copy(y_cel_inter)
        min_axis = axis_direction
        inter_parameters[4] = min_axis
    print("Sum of squared deviations: " + "{}".format(min_value))
    # print(inter_parameters)
    return inter_parameters, min_value



def axis_fit_use_distorted(supposed_x, supposed_y, input_parameters):
    # Convert from celestial projection coordinates into frame coordinates
    # input_x, input_y - distorted coordinates of anchors on the frame
    # Initial center
    frame_centre_init_x = input_parameters[0] * cos(input_parameters[2]) + input_parameters[1] * sin(
        input_parameters[2])
    frame_centre_init_y = - input_parameters[0] * sin(input_parameters[2]) + input_parameters[1] * cos(
        input_parameters[2])
    # Converting celestial coordinates to corrected celestial ones
    supposed_x_1 = supposed_x * cos(input_parameters[11]) + supposed_y * sin(input_parameters[11]) + input_parameters[
        9]
    supposed_y_1 = -supposed_x * sin(input_parameters[11]) + supposed_y * cos(input_parameters[11]) + \
                 input_parameters[10]
    #supposed_x_1 = (supposed_x - input_parameters[9]) * cos(- input_parameters[11]) + (
    #            supposed_y - input_parameters[10]) * sin(- input_parameters[11])
    #supposed_y_1 = -(supposed_x - input_parameters[9]) * sin(- input_parameters[11]) + (
    #            supposed_y - input_parameters[10]) * cos(- input_parameters[11])
    frame_centre_x = frame_centre_init_x
    frame_centre_y = frame_centre_init_y
    init_frame_x = supposed_x_1
    init_frame_y = supposed_y_1

    prev_coord_arr = array([[init_frame_x], [init_frame_y]])
    coord_diff = 50
    while coord_diff >= 0.01:
        centered_x = prev_coord_arr[0][0] - frame_centre_x  # x as measured from the center of the full frame
        centered_y = prev_coord_arr[1][0] - frame_centre_y  # y as measured from the center of the full frame
        deviation_x = supposed_x_1 - prev_coord_arr[0][0]
        deviation_y = supposed_y_1 - prev_coord_arr[1][0]
        centered_radius = sqrt(centered_x ** 2 + centered_y ** 2)  # Distance from the center of the frame (distorted)
        # Additional variables
        a_k = centered_radius ** 2 + 2 * centered_x ** 2
        b_k = 2 * centered_x * centered_y
        c_k = centered_radius ** 2 + 2 * centered_y ** 2
        d_k = centered_radius ** 4 + 4 * centered_radius ** 2 * centered_x ** 2
        e_k = centered_radius ** 4 + 4 * centered_radius ** 2 * centered_y ** 2
        #f_fun = -init_frame_x + prev_coord_arr[0][0] + centered_x * (
        #            input_parameters[5] * centered_radius ** 2 + input_parameters[6] * centered_radius ** 4) + (
        #                    input_parameters[7] * (centered_radius ** 2 + 2 * centered_x ** 2) + input_parameters[
        #                8] * 2 * centered_x * centered_y)
        #g_fun = - init_frame_y + prev_coord_arr[1][0] + centered_y * (
        #        input_parameters[5] * centered_radius ** 2 + input_parameters[6] * centered_radius ** 4) + (
        #                input_parameters[8] * (centered_radius ** 2 + 2 * centered_y ** 2) + input_parameters[
        #            7] * 2 * centered_x * centered_y)

        f_fun = deviation_x - centered_x * (
                input_parameters[5] * centered_radius ** 2 + input_parameters[6] * centered_radius ** 4)

        f_fun -= input_parameters[7] * a_k + input_parameters[8] * b_k
        g_fun = deviation_y - centered_y * (
                input_parameters[5] * centered_radius ** 2 + input_parameters[6] * centered_radius ** 4)
        g_fun -= input_parameters[7] * b_k + input_parameters[8] * c_k
        f_arr = array([[f_fun], [g_fun]]) * (-1)
        f_fun_x_der = 1 + input_parameters[5] * a_k + input_parameters[6] * d_k + 6 * input_parameters[7] * centered_x + 2 * input_parameters[8] * centered_y
        g_fun_y_der = 1 + input_parameters[5] * c_k + input_parameters[6] * e_k + 2 * input_parameters[7] * centered_x + 6 * \
                      input_parameters[8] * centered_y
        f_g_fun_common_der = input_parameters[5] * b_k + 2 * input_parameters[6] * centered_radius ** 2 * b_k + 2 * input_parameters[7] * centered_y + 2 * input_parameters[8] * centered_x

        der_arr = array([[f_fun_x_der, f_g_fun_common_der], [f_g_fun_common_der, g_fun_y_der]])
        #print(der_arr)
        #print(f_arr)
        der_arr_inv = invert_array(der_arr)
        right_side_arr = dot(der_arr_inv, f_arr)
        coord_arr = prev_coord_arr - right_side_arr
        coord_arr_diff = absolute(coord_arr - prev_coord_arr)
        coord_diff = max(coord_arr_diff)
        prev_coord_arr = copy(coord_arr)
    return coord_arr[0], coord_arr[1]

def axis_fit_use(input_ra, input_dec, mid_ra, mid_dec, frame_parameters):
    # Making a celestial projection array
    x_cel = 0
    y_cel = 0
    x_cel = frame_parameters[3] * rad_conversion * (
                cos(input_dec / rad_conversion) * sin(mid_dec / rad_conversion) * cos(
            (input_ra - mid_ra) / rad_conversion) - sin(input_dec / rad_conversion) * cos(
            mid_dec / rad_conversion))
    y_cel = frame_parameters[3] * rad_conversion * (
                sin((input_ra - mid_ra) / rad_conversion) * cos(input_dec / rad_conversion))

    # Finding coordinates
    if frame_parameters[4] == 1:
        x_cel = x_cel * (-1)
    elif frame_parameters[4] == 2:
        y_cel = y_cel * (-1)
    elif frame_parameters[4] == 3:
        y_cel = y_cel * (-1)
        x_cel = x_cel * (-1)

    x_theor = x_cel * cos(frame_parameters[2]) + y_cel * sin(frame_parameters[2])
    y_theor = -x_cel * sin(frame_parameters[2]) + y_cel * cos(frame_parameters[2])
    x_output = x_theor + frame_parameters[0] * cos(frame_parameters[2]) + frame_parameters[1] * sin(frame_parameters[2])
    y_output = y_theor - frame_parameters[0] * sin(frame_parameters[2]) + frame_parameters[1] * cos(frame_parameters[2])

    del x_theor, y_theor, x_cel, y_cel, mid_ra, mid_dec
    return x_output, y_output

def axis_fit_use_all(input_ra, input_dec, mid_ra, mid_dec, input_parameters):       # For a singular object
    # input_ra, input_dec - corresponding celestial coordinates
    # axis_direction - check if a frame is mirrored in some way
    x_cel = rad_conversion * (cos(input_dec / rad_conversion) * sin(mid_dec / rad_conversion) * cos(
        (input_ra - mid_ra) / rad_conversion) - sin(input_dec / rad_conversion) * cos(
        mid_dec / rad_conversion))
    y_cel = rad_conversion * (sin((input_ra - mid_ra) / rad_conversion) * cos(input_dec / rad_conversion))
    # Making a celestial projection array
    supposed_x, supposed_y = axis_fit_use(input_ra, input_dec, mid_ra, mid_dec,  input_parameters)
    final_frame_x, final_frame_y = axis_fit_use_distorted(supposed_x, supposed_y, input_parameters)
    return final_frame_x, final_frame_y


# Functions for defining quads
def quad_coord(length_diagonal, coord_c, coord_a, diag_angle):
    # Positional angles for inner stars
    length_1 = sqrt((coord_c[0] - coord_a[0]) ** 2 + (coord_c[1] - coord_a[1]) ** 2)
    angle1 = Angle(atan2(-(coord_c[1] - coord_a[1]), -(coord_c[0] - coord_a[0])),
                   u.rad).wrap_at('360d').rad
    # angle1 = math.asin((coord_c[1] - coord_a[1]) / length_1)
    anglex = diag_angle - 45. / 180. * 3.14159
    angle_1_x = angle1 - anglex
    corr_x_1 = length_1 * cos(angle_1_x)
    corr_y_1 = length_1 * sin(angle_1_x)
    side = length_diagonal * sqrt(2) / 2
    quad_x_1 = corr_x_1 / side
    quad_y_1 = corr_y_1 / side
    return quad_x_1, quad_y_1


def quad_coord_err(length_diagonal, coord_c, coord_a, diag_angle, coord_err):
    quad_1_x, quad_1_y = quad_coord(length_diagonal, coord_c, coord_a, diag_angle)

    coord_c_devia = zeros(2)
    coord_c_devia[0] = coord_c[0] + coord_err
    coord_c_devia[1] = coord_c[1] + coord_err
    coord_a_devia = zeros(2)
    coord_a_devia[0] = coord_a[0] + coord_err
    coord_a_devia[1] = coord_a[1] + coord_err
    quad_2_c_x, quad_2_c_y = quad_coord(length_diagonal, coord_c_devia, coord_a, diag_angle)
    quad_2_a_x, quad_2_a_y = quad_coord(length_diagonal, coord_c, coord_a_devia, diag_angle)
    error_x = sqrt((quad_2_a_x - quad_1_x) ** 2 + (quad_2_c_x - quad_1_x) ** 2)
    error_y = sqrt((quad_2_a_y - quad_1_y) ** 2 + (quad_2_c_y - quad_1_y) ** 2)
    return error_x, error_y


# A quad
def quad_frame(coord_x, coord_y):
    index_arr = arange(0, 4)
    index_arr += 1
    index_comb = list(combinations(index_arr, 2))
    max_len = 0

    for a in index_comb:
        x1 = coord_x[a[0]]
        x2 = coord_x[a[1]]
        y1 = coord_y[a[0]]
        y2 = coord_y[a[1]]
        length = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        if length > max_len:
            max_comb = a
            max_len = length

    quad_index = zeros(4)
    # Picking the farthest stars
    if coord_x[max_comb[0]] >= coord_x[max_comb[1]]:
        max_x1 = coord_x[max_comb[0]]
        max_x2 = coord_x[max_comb[1]]
        max_y1 = coord_y[max_comb[0]]
        max_y2 = coord_y[max_comb[1]]
        quad_index[0] = max_comb[0]
        quad_index[1] = max_comb[1]
    else:
        max_x1 = coord_x[max_comb[1]]
        max_x2 = coord_x[max_comb[0]]
        max_y1 = coord_y[max_comb[1]]
        max_y2 = coord_y[max_comb[0]]
        quad_index[0] = max_comb[1]
        quad_index[1] = max_comb[0]
    diag_pos_angle = Angle(atan2(-(max_y2 - max_y1), -(max_x2 - max_x1)), u.rad).wrap_at(
        '360d').rad

    # Picking the inner stars
    inner_index = []
    for i in range(1, 5):
        if i != int(max_comb[0]) and i != int(max_comb[1]):
            inner_index = append(inner_index, i)
    inner_index = inner_index.astype(int)
    quad_index[2] = inner_index[0]
    quad_index[3] = inner_index[1]

    coord_a = [max_x1, max_y1]
    coord_c = [coord_x[inner_index[0]], coord_y[inner_index[0]]]
    coord_d = [coord_x[inner_index[1]], coord_y[inner_index[1]]]
    # Positional angles for inner stars
    quad_x_1, quad_y_1 = quad_coord(max_len, coord_c, coord_a, diag_pos_angle)
    quad_x_2, quad_y_2 = quad_coord(max_len, coord_d, coord_a, diag_pos_angle)
    quad_x_1_err, quad_y_1_err = quad_coord_err(max_len, coord_c, coord_a, diag_pos_angle, 3)
    quad_x_2_err, quad_y_2_err = quad_coord_err(max_len, coord_d, coord_a, diag_pos_angle, 3)
    quad_arr = [coord_x[0], quad_x_1, quad_x_1_err, quad_y_1, quad_y_1_err, quad_x_2, quad_x_2_err, quad_y_2,
                quad_y_2_err, max_len]
    quad_arr = append(quad_arr, quad_index)
    return quad_arr


# Quads from the GAIA catalogue
def gaia_bright(catalog_path, anchor_num, offset, radius, obj_ra, obj_dec, check_up, anchor_up, mag_low, source_mag_int):
    center_ra = obj_ra
    center_dec = obj_dec

    data = read_csv(catalog_path)
    data = data.sort_values('Gmag')
    # print(radius)
    # Pick the brightest 10
    if offset >=2:
        data_bright = data.nsmallest(30, 'Gmag')[offset - 2:offset + 3 * anchor_num]  # The brightest with an offset
    else:
        data_bright = data.nsmallest(30, 'Gmag')[:3 * anchor_num]
    # Discarding stars from the data_bright that are too close to each other (less than 1 minute)
    dict_stars = []
    array_columns = data_bright.columns.values.tolist()

    for bright_first_ind in range(len(data_bright['Gmag'])):
        add_check = 1
        for bright_second_ind in range(len(data_bright['Gmag'])):
            ang_sep = angular_separation(data_bright['RA_ICRS'].values[bright_first_ind] / rad_conversion,
                                         data_bright['DE_ICRS'].values[bright_first_ind] / rad_conversion,
                                         data_bright['RA_ICRS'].values[bright_second_ind] / rad_conversion,
                                         data_bright['DE_ICRS'].values[bright_second_ind] / rad_conversion)
            #print(ang_sep * rad_conversion * 60)
            if ang_sep != 0. and ang_sep <= 0.5 / 60 / rad_conversion:
                add_check = 0
        if add_check == 1:
            dict1 = dict((col_name, data_bright[col_name].values[bright_first_ind]) for col_name in array_columns)
            dict_stars.append(dict1)
    data_bright = pd.DataFrame(dict_stars, columns=array_columns)       # Stars for initial anchors

    super_bright = data.loc[(data['Gmag'] <= 10.5)]  # The brightest stars
    # Stars with GMAG between 12 and 16
    inter_bright = data.loc[(data['Gmag'] <= mag_low) & (data['Gmag'] >= check_up)]
    if len(inter_bright['Gmag']) >=200:                                 # Limiting the number of check stars
        inter_bright = inter_bright.nsmallest(200, 'Gmag')
    if mag_low < 15:
        anchor_stars = inter_bright.loc[(inter_bright['Gmag'] <= mag_low) & (inter_bright['Gmag'] >= anchor_up)]
    else:
        anchor_stars = inter_bright.loc[(inter_bright['Gmag'] <= 15) & (inter_bright['Gmag'] >= anchor_up)]
    if len(anchor_stars['Gmag']) >=200:
        anchor_stars = anchor_stars.nsmallest(200, 'Gmag')
    print("Number of initial anchors: " + "{}".format(len(anchor_stars)))
    stars_for_checking = data.loc[(data['Gmag'] <= mag_low + 2) & (data['Gmag'] >= 11)]
    stars_for_checking.sort_values('Gmag')
    array_columns = anchor_stars.columns.values.tolist()
    # Picking possible check stars
    poss_check_stars = data.loc[(data['Gmag'] <= source_mag_int[1]) & (data['Gmag'] >= source_mag_int[0])]
    # Discarding check stars that are close to superbright ones
    check_dict = []
    for check_ind in range(len(poss_check_stars['Gmag'])):
        adding_check = 0
        if isnan(poss_check_stars['BP-RP'].values[check_ind]) == 0 and poss_check_stars['BP-RP'].values[check_ind]<=4:
            for super_bright_ind in range(len(super_bright['Gmag'])):
                ang_sep = angular_separation(poss_check_stars['RA_ICRS'].values[check_ind] / rad_conversion,
                                             poss_check_stars['DE_ICRS'].values[check_ind] / rad_conversion,
                                             super_bright['RA_ICRS'].values[super_bright_ind] / rad_conversion,
                                             super_bright['DE_ICRS'].values[super_bright_ind] / rad_conversion)
                if ang_sep <= 1. / 60 / rad_conversion and ang_sep != 0:
                    adding_check = 1
                    break
            if adding_check == 0:
                ang_sep = angular_separation(poss_check_stars['RA_ICRS'].values[check_ind] / rad_conversion,
                                             poss_check_stars['DE_ICRS'].values[check_ind] / rad_conversion,
                                             center_ra / rad_conversion, center_dec / rad_conversion)
                if ang_sep > radius / 60 / rad_conversion:
                    adding_check = 1
            if adding_check == 0:
                ang_sep = angular_separation(poss_check_stars['RA_ICRS'].values[check_ind] / rad_conversion,
                                             poss_check_stars['DE_ICRS'].values[check_ind] / rad_conversion,
                                             obj_ra / rad_conversion, obj_dec / rad_conversion)
                if ang_sep < 0.5 / 60 / rad_conversion:
                    adding_check = 1
            if adding_check == 0:
                dict1 = dict((col_name, poss_check_stars[col_name].values[check_ind]) for col_name in array_columns)
                check_dict.append(dict1)
    check_stars = pd.DataFrame(check_dict, columns=array_columns)
    # Discarding check stars that are close to stars with comparable brightness
    check_dict = []
    for check_ind in range(len(check_stars['Gmag'])):
        adding_check = 0
        for checking_ind in range(len(stars_for_checking['Gmag'])):
            mag_diff = check_stars['Gmag'].values[check_ind] - stars_for_checking['Gmag'].values[checking_ind]
            if mag_diff >= -3:
                ang_sep = angular_separation(check_stars['RA_ICRS'].values[check_ind] / rad_conversion,
                                             check_stars['DE_ICRS'].values[check_ind] / rad_conversion,
                                             stars_for_checking['RA_ICRS'].values[checking_ind] / rad_conversion,
                                             stars_for_checking['DE_ICRS'].values[checking_ind] / rad_conversion)
                if ang_sep <= 0.3 / 60 / rad_conversion and ang_sep != 0:
                    adding_check = 1
                    break
        if adding_check == 0:
            for check_two_ind in range(len(check_stars['Gmag'])):
                mag_diff = check_stars['Gmag'].values[check_ind] - check_stars['Gmag'].values[check_two_ind]
                if mag_diff >= -3:
                    ang_sep = angular_separation(check_stars['RA_ICRS'].values[check_ind] / rad_conversion,
                                                 check_stars['DE_ICRS'].values[check_ind] / rad_conversion,
                                                 check_stars['RA_ICRS'].values[check_two_ind] / rad_conversion,
                                                 check_stars['DE_ICRS'].values[check_two_ind] / rad_conversion)
                    if ang_sep <= 0.3 / 60 / rad_conversion and ang_sep != 0:
                        adding_check = 1
                        break
        if adding_check == 0:
            dict1 = dict((col_name, check_stars[col_name].values[check_ind]) for col_name in array_columns)
            check_dict.append(dict1)
    check_stars = pd.DataFrame(data=check_dict, columns=array_columns)
    # Discarding checks in magnitude intervals (of 0.5 mag)
    mag_limits = np.arange(source_mag_int[0], source_mag_int[1], 0.25)
    check_dict = []
    for mag_ind in range(len(mag_limits)):
        mag_small = float(mag_limits[mag_ind])
        mag_big = mag_small + 0.25
        #print(mag_small)
        sample_check_stars = check_stars.loc[(check_stars['Gmag'] >= mag_small) & (check_stars['Gmag'] <= mag_big)]
        sample_check_stars = sample_check_stars.nsmallest(20, 'Gmag')
        for star_ind in range(len(sample_check_stars['Gmag'])):
            dict1 = dict((col_name, sample_check_stars[col_name].values[star_ind]) for col_name in array_columns)
            check_dict.append(dict1)
    check_stars = pd.DataFrame(data=check_dict, columns=array_columns)
    # Discarding check stars that are too close to super bright stars and to the object, as well as those that are beyond the frame
    dict_stars = []
    array_columns = anchor_stars.columns.values.tolist()

    for intermediate_ind in range(len(inter_bright['Gmag'])):
        adding_check = 0
        if isnan(inter_bright['BP-RP'].values[intermediate_ind]) == 0 and inter_bright['BP-RP'].values[intermediate_ind]<=4:
            for super_bright_ind in range(len(super_bright['Gmag'])):
                ang_sep = angular_separation(inter_bright['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                             inter_bright['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                             super_bright['RA_ICRS'].values[super_bright_ind] / rad_conversion,
                                             super_bright['DE_ICRS'].values[super_bright_ind] / rad_conversion)
                if ang_sep <= 1. / 60 / rad_conversion and ang_sep != 0:
                    adding_check = 1
                    break
            for check_ind in range(len(check_stars['Gmag'])):
                ng_sep = angular_separation(inter_bright['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                            inter_bright['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                            check_stars['RA_ICRS'].values[check_ind] / rad_conversion,
                                            check_stars['DE_ICRS'].values[check_ind] / rad_conversion)
                if ang_sep <= 0.2 / 60 / rad_conversion:
                    adding_check = 1
                    break
            if adding_check == 0:
                ang_sep = angular_separation(inter_bright['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                             inter_bright['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                             center_ra / rad_conversion, center_dec / rad_conversion)
                if ang_sep > radius / 60 / rad_conversion:
                    adding_check = 1
            if adding_check == 0:
                ang_sep = angular_separation(inter_bright['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                             inter_bright['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                             obj_ra / rad_conversion, obj_dec / rad_conversion)
                if ang_sep < 0.5 / 60 / rad_conversion:
                    adding_check = 1
            if adding_check == 0:
                dict1 = dict((col_name, inter_bright[col_name].values[intermediate_ind]) for col_name in array_columns)
                dict_stars.append(dict1)

    inter_bright = pd.DataFrame(dict_stars, columns=array_columns)
    # Discarding anchors that are too close to bright intermediary stars, super bright ones and beyond the center
    dict_stars = []
    final_anchor_num = 0
    for anchor_ind in range(len(anchor_stars['Gmag'])):
        adding_check = 0
        if final_anchor_num >=40:
            break
        for intermediate_ind in range(len(stars_for_checking['Gmag'])):
            ang_sep = angular_separation(stars_for_checking['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                         stars_for_checking['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                         anchor_stars['RA_ICRS'].values[anchor_ind] / rad_conversion,
                                         anchor_stars['DE_ICRS'].values[anchor_ind] / rad_conversion)
            mag_diff = stars_for_checking['Gmag'].values[intermediate_ind] - anchor_stars['Gmag'].values[anchor_ind]
            if ang_sep <= 0.2 / 60 / rad_conversion and ang_sep != 0 and mag_diff <= 1:
                adding_check = 1
                #print('cccc')
                break
        if adding_check == 0:
            for super_bright_ind in range(len(super_bright)):
                ang_sep = angular_separation(anchor_stars['RA_ICRS'].values[anchor_ind] / rad_conversion,
                                             anchor_stars['DE_ICRS'].values[anchor_ind] / rad_conversion,
                                             super_bright['RA_ICRS'].values[super_bright_ind] / rad_conversion,
                                             super_bright['DE_ICRS'].values[super_bright_ind] / rad_conversion)
                # print(ang_sep)
                if ang_sep <= 2. / 60 / rad_conversion and ang_sep != 0:
                    #print('bbbb')
                    adding_check = 1
                    break
            if adding_check == 0:
                ang_sep = angular_separation(anchor_stars['RA_ICRS'].values[anchor_ind] / rad_conversion,
                                             anchor_stars['DE_ICRS'].values[anchor_ind] / rad_conversion,
                                             center_ra / rad_conversion, center_dec / rad_conversion)

                if ang_sep > radius / 60 / rad_conversion:
                    #print('ddddd')
                    adding_check = 1
            if adding_check == 0:
                dict1 = dict((col_name, anchor_stars[col_name].values[anchor_ind]) for col_name in array_columns)
                dict_stars.append(dict1)
                final_anchor_num +=1
    inter_anchor_stars = pd.DataFrame(dict_stars, columns=array_columns)
    #print(inter_anchor_stars['RA_ICRS'])
    # Discarding anchors that are too close to each other
    dict_stars = []
    for x_ind in range(len(inter_anchor_stars['Gmag'])):
        adding_check = 1
        for y_ind in range(len(inter_anchor_stars['Gmag'])):
            ang_sep = angular_separation(inter_anchor_stars['RA_ICRS'].values[x_ind] / rad_conversion,
                                         inter_anchor_stars['DE_ICRS'].values[x_ind] / rad_conversion,
                                         inter_anchor_stars['RA_ICRS'].values[y_ind] / rad_conversion,
                                         inter_anchor_stars['DE_ICRS'].values[y_ind] / rad_conversion)
            #print(ang_sep)
            if ang_sep <= 0.5 / 60. / rad_conversion and ang_sep !=0.0:
                adding_check = 0
                break
        if adding_check == 1:
            dict1 = dict((col_name, inter_anchor_stars[col_name].values[x_ind]) for col_name in array_columns)
            dict_stars.append(dict1)
    final_anchor_stars = pd.DataFrame(dict_stars, columns=array_columns)


    # Discarding ref stars that are too close to each other
    dict_stars = []
    for intermediate_ind in range(len(inter_bright['Gmag'])):
        adding_check = 0
        for check_ind in range(len(stars_for_checking['Gmag'])):
            mag_diff = inter_bright['Gmag'].values[intermediate_ind] - stars_for_checking['Gmag'].values[check_ind]
            if mag_diff >= -3:
                ang_sep = angular_separation(inter_bright['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                             inter_bright['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                             stars_for_checking['RA_ICRS'].values[check_ind] / rad_conversion,
                                             stars_for_checking['DE_ICRS'].values[check_ind] / rad_conversion)
                if ang_sep <= 0.5 / 60 / rad_conversion and ang_sep != 0:
                    adding_check = 1
                    break
        if adding_check == 0:
            for inter_ind in range(len(inter_bright['Gmag'])):
                ang_sep = angular_separation(inter_bright['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                             inter_bright['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                             inter_bright['RA_ICRS'].values[inter_ind] / rad_conversion,
                                             inter_bright['DE_ICRS'].values[inter_ind] / rad_conversion)
                if ang_sep <= 0.5 / 60 / rad_conversion and ang_sep != 0:
                    adding_check = 1
                    break
        if adding_check == 0:
            dict1 = dict((col_name, inter_bright[col_name].values[intermediate_ind]) for col_name in array_columns)
            dict_stars.append(dict1)
    inter_bright = pd.DataFrame(dict_stars, columns=array_columns)

    # Correcting magnitudes of inter_bright for neighbouring stars (for aperture radiuses 2, 4, 6 pixels)
    inter_bright_gmag_2_pix = copy(inter_bright['Gmag'])
    inter_bright['Gmag_2'] = copy(inter_bright_gmag_2_pix)
    inter_bright['Gmag_4'] = copy(inter_bright_gmag_2_pix)
    inter_bright['Gmag_6'] = copy(inter_bright_gmag_2_pix)
    inter_bright['Gmag_8'] = copy(inter_bright_gmag_2_pix)
    for intermediate_ind in range(len(inter_bright['Gmag'])):
        for star_ind in range(len(data['Gmag'])):
            ang_sep = angular_separation(inter_bright['RA_ICRS'].values[intermediate_ind] / rad_conversion,
                                         inter_bright['DE_ICRS'].values[intermediate_ind] / rad_conversion,
                                         data['RA_ICRS'].values[star_ind] / rad_conversion,
                                         data['DE_ICRS'].values[star_ind] / rad_conversion)
            if ang_sep !=0:
                if ang_sep <= 4./206265:
                    inter_bright['Gmag_2'].values[intermediate_ind] = combined_magnitude(inter_bright['Gmag'].values[intermediate_ind], data['Gmag'].values[star_ind])

                if ang_sep >= 4./206265 and ang_sep <= 8./206265:
                    inter_bright['Gmag_4'].values[intermediate_ind] = combined_magnitude(inter_bright['Gmag'].values[intermediate_ind], data['Gmag'].values[star_ind])

                if ang_sep >= 8./206265 and ang_sep <= 12./206265:
                    inter_bright['Gmag_6'].values[intermediate_ind] = combined_magnitude(inter_bright['Gmag'].values[intermediate_ind], data['Gmag'].values[star_ind])

                if ang_sep >= 12. / 206265 and ang_sep <= 16. / 206265:
                    inter_bright['Gmag_8'].values[intermediate_ind] = combined_magnitude(inter_bright['Gmag'].values[intermediate_ind],
                                                                data['Gmag'].values[star_ind])

    check_stars_gmag_2_pix = copy(check_stars['Gmag'])
    check_stars['Gmag_2'] = copy(check_stars_gmag_2_pix)
    check_stars['Gmag_4'] = copy(check_stars_gmag_2_pix)
    check_stars['Gmag_6'] = copy(check_stars_gmag_2_pix)
    check_stars['Gmag_8'] = copy(check_stars_gmag_2_pix)
    for check_ind in range(len(check_stars['Gmag'])):
        for star_ind in range(len(data['Gmag'])):
            ang_sep = angular_separation(check_stars['RA_ICRS'].values[check_ind] / rad_conversion,
                                         check_stars['DE_ICRS'].values[check_ind] / rad_conversion,
                                         data['RA_ICRS'].values[star_ind] / rad_conversion,
                                         data['DE_ICRS'].values[star_ind] / rad_conversion)
            if ang_sep != 0:
                if ang_sep <= 4. / 206265:
                    check_stars['Gmag_2'].values[check_ind] = combined_magnitude(
                        check_stars['Gmag'].values[check_ind], data['Gmag'].values[star_ind])

                if ang_sep >= 4. / 206265 and ang_sep <= 8. / 206265:
                    check_stars['Gmag_4'].values[check_ind] = combined_magnitude(
                        check_stars['Gmag'].values[check_ind], data['Gmag'].values[star_ind])

                if ang_sep >= 8. / 206265 and ang_sep <= 12. / 206265:
                    check_stars['Gmag_6'].values[check_ind] = combined_magnitude(
                        check_stars['Gmag'].values[check_ind], data['Gmag'].values[star_ind])

                if ang_sep >= 12. / 206265 and ang_sep <= 16. / 206265:
                    check_stars['Gmag_8'].values[check_ind] = combined_magnitude(
                        check_stars['Gmag'].values[check_ind],
                        data['Gmag'].values[star_ind])

    color_arr = inter_bright.nsmallest(1, 'BP-RP')
    #print(color_arr)
    min_color = color_arr['BP-RP'].values[0]
    color_arr = inter_bright.nlargest(1, 'BP-RP')
    max_color = color_arr['BP-RP'].values[0]

    check_stars = check_stars.loc[(check_stars['BP-RP'] <= max_color) & (check_stars['BP-RP'] >= min_color)]
    return data_bright, inter_bright, final_anchor_stars, data, super_bright, check_stars

# Calculating a combined magnitude of 2 stars through Pogson

def combined_magnitude(mag1, mag2):
    if isnan(mag2) == 0:
        I1 = 1
        I2 = 10 ** (0.4 * (mag1 - mag2))
        I_comb = I1 + I2
        m_comb = mag1 - 2.5 * log10(I_comb)
    else:
        m_comb = mag1
    return m_comb


def check_ref_search(star_catalog, center_ra, center_dec):
    # print(star_catalog)
    #print(star_catalog['BP-RP'])
    min_color = min(star_catalog['BP-RP'])
    max_color = arr_max(star_catalog['BP-RP'])
    print('Maximum BP-RP among check ref stars: ' + "{}".format(max_color))
    color_number_array = zeros(int(round((max_color - min_color) / 0.1)) + 1)
    # print(len(color_number_array))
    dict_stars = []
    dict_ref = []
    array_columns = star_catalog.columns.values.tolist()
    for check_ind in range(len(star_catalog['BP-RP'])):
        star_color = star_catalog['BP-RP'].values[check_ind]
        if isnan(star_color):
            continue
        else:
            # print(star_color, star_catalog['Gmag'].values[check_ind])
            color_ind = int((star_color - min_color) / 0.1)
            #ang_sep = angular_separation(star_catalog['RA_ICRS'].values[check_ind] / rad_conversion,
            #                             star_catalog['DE_ICRS'].values[check_ind] / rad_conversion,
            #                             center_ra / rad_conversion, center_dec / rad_conversion)
            if color_number_array[color_ind] < 10:
                dict1 = dict((col_name, star_catalog[col_name].values[check_ind]) for col_name in array_columns)
                color_number_array[color_ind] += 1
                dict_stars.append(dict1)
    check_catalog = pd.DataFrame(dict_stars, columns=array_columns)
    return check_catalog, color_number_array


def quad_cel(input_RA, input_DEC):
    # Picking the farthest stars

    coord_RA = copy(input_RA)
    coord_DEC = copy(input_DEC)
    index_arr = arange(1, 5)
    index_comb = list(combinations(index_arr, 2))
    max_len = 0
    rad_conversion = 180. / 3.14159

    # Converting to radian
    for i in range(1, 5):
        coord_RA[i] = coord_RA[i] / rad_conversion
        coord_DEC[i] = coord_DEC[i] / rad_conversion
    for a in index_comb:
        RA_1 = coord_RA[a[0]]
        RA_2 = coord_RA[a[1]]
        DEC_1 = coord_DEC[a[0]]
        DEC_2 = coord_DEC[a[1]]
        length = angular_separation(RA_1, DEC_1, RA_2, DEC_2)
        if length > max_len:
            max_comb = a
            max_len = length

    # Picking the farthest stars
    quad_index = zeros(4)
    if coord_DEC[max_comb[0]] <= coord_DEC[max_comb[1]]:
        max_RA_1 = coord_RA[max_comb[0]]
        max_RA_2 = coord_RA[max_comb[1]]
        max_DEC_1 = coord_DEC[max_comb[0]]
        max_DEC_2 = coord_DEC[max_comb[1]]
        quad_index[0] = max_comb[0]
        quad_index[1] = max_comb[1]
    else:
        max_RA_1 = coord_RA[max_comb[1]]
        max_RA_2 = coord_RA[max_comb[0]]
        max_DEC_1 = coord_DEC[max_comb[1]]
        max_DEC_2 = coord_DEC[max_comb[0]]
        quad_index[0] = max_comb[1]
        quad_index[1] = max_comb[0]
    diag_pos_angle = position_angle(max_RA_1, max_DEC_1, max_RA_2, max_DEC_2).rad

    # Picking the inner stars
    inner_index = []
    for i in range(1, 5):
        if i != int(max_comb[0]) and i != int(max_comb[1]):
            inner_index = append(inner_index, i)
    inner_index = inner_index.astype(int)
    quad_index[2] = inner_index[0]
    quad_index[3] = inner_index[1]

    # Positional angles for inner stars
    len1 = angular_separation(coord_RA[inner_index[0]], coord_DEC[inner_index[0]], max_RA_1, max_DEC_1)
    len2 = angular_separation(coord_RA[inner_index[1]], coord_DEC[inner_index[1]], max_RA_1, max_DEC_1)
    angle1 = position_angle(max_RA_1, max_DEC_1, coord_RA[inner_index[0]], coord_DEC[inner_index[0]]).rad
    angle2 = position_angle(max_RA_1, max_DEC_1, coord_RA[inner_index[1]], coord_DEC[inner_index[1]]).rad
    anglex = diag_pos_angle - 45. / 180. * 3.14159

    angle_1_x = angle1 - anglex
    angle_2_x = angle2 - anglex
    corr_x_1 = len1 * cos(angle_1_x)
    corr_y_1 = len1 * sin(angle_1_x)
    corr_x_2 = len2 * cos(angle_2_x)
    corr_y_2 = len2 * sin(angle_2_x)
    side = max_len * sqrt(2) / 2
    quad_x_1 = corr_x_1 / side
    quad_y_1 = corr_y_1 / side
    quad_x_2 = corr_x_2 / side
    quad_y_2 = corr_y_2 / side
    quad_arr = [coord_RA[0], quad_x_1, quad_y_1, quad_x_2, quad_y_2, max_len]
    quad_arr = append(quad_arr, quad_index)
    return quad_arr


def quad_catalog(input_data):
    index_arr = arange(len(input_data['RA_ICRS']))
    index_comb = list(combinations(index_arr, 4))
    comb_num = len(index_comb)
    comb_RA = zeros((comb_num, 5))
    comb_DEC = zeros((comb_num, 5))
    for i in range(comb_num):
        comb_RA[i][0] = i
        comb_DEC[i][0] = i
        for j in range(4):
            comb_RA[i][j + 1] = input_data['RA_ICRS'].values[index_comb[i][j]]
            comb_DEC[i][j + 1] = input_data['DE_ICRS'].values[index_comb[i][j]]
            # comb_RA[i][j + 1] = input_data['RAJ2000'].values[index_comb[i][j]]
            # comb_DEC[i][j + 1] = input_data['DEJ2000'].values[index_comb[i][j]]
    quad_all_arr = zeros((comb_num, 10))
    for i in range(comb_num):
        intermediate_arr = quad_cel(comb_RA[i], comb_DEC[i])
        for j in range(intermediate_arr.shape[0]):
            quad_all_arr[i][j] = intermediate_arr[j]

    return quad_all_arr, comb_RA, comb_DEC


def dec_line(frame_parameters, dec, ra, ground, frame_path):
    net_param = copy(frame_parameters)
    dec_lower = dec - 0.3
    dec_upper = dec + 0.3
    # x_painted = []
    # y_painted = []
    line_mid_x = ground.shape[0] / 2
    line_mid_y = ground.shape[1] / 2
    dec_arr = arange(dec_lower, dec_upper, 0.001)
    for i in dec_arr:
        # print(i)
        x_inter, y_inter = axis_fit_use(ra, i, frame_path, net_param)
        print(x_inter, y_inter)
        if abs(x_inter - line_mid_x) < line_mid_x - 10 and abs(y_inter - line_mid_y) < line_mid_y - 10:
            ground[int(x_inter)][int(y_inter)] += 10 ** 4
    return 0


# Gauss elimination


def gauss_elimination(a_array, b_array):
    # div_array = zeros(a_array.shape[0])     # Array for storing divisors
    a_array_copy = copy(a_array)
    b_array_copy = copy(b_array)
    for a_ind in range(a_array.shape[0] - 1):

        for inner_ind in range(a_ind, a_array.shape[0]):
            div_value = a_array_copy[inner_ind][a_ind]
            if div_value == 0:
                break
            # print(div_value)
            for recurs_ind in range(a_ind, a_array.shape[0]):
                a_array_copy[inner_ind][recurs_ind] = a_array_copy[inner_ind][recurs_ind] / div_value
            b_array_copy[inner_ind] = b_array_copy[inner_ind] / div_value
        if div_value == 0:
            continue
        substract_row = a_array_copy[a_ind]
        substract_value = b_array_copy[a_ind]

        for another_ind in range(a_ind + 1, a_array.shape[0]):
            a_array_copy[another_ind] -= substract_row
            b_array_copy[another_ind] -= substract_value
    return a_array_copy, b_array_copy


def distortion_min_sum(input_x, input_y, supposed_x_0, supposed_y_0, input_parameters):
    # input_x, input_y - distorted coordinates of anchors on the frame
    # supposed_x, supposed_y - undistorted (supposedly) coordinates of anchors on the frame
    # input_parameters - rotation, shift and scale parameters
    # print(input_x, input_y)
    # print(supposed_x, supposed_y)
    # print(supposed_x - input_x, supposed_y - input_y)
    # Initial center
    frame_centre_init_x = input_parameters[0] * cos(input_parameters[2]) + input_parameters[1] * sin(
        input_parameters[2])
    frame_centre_init_y = - input_parameters[0] * sin(input_parameters[2]) + input_parameters[1] * cos(
        input_parameters[2])
    # print(input_parameters)

    # supposed_x = (supposed_x_0) * cos(input_parameters[11]) + (supposed_y_0 + input_parameters[10]) * sin(input_parameters[11]) - input_parameters[9]
    # supposed_y = - (supposed_x_0 + input_parameters[9]) * sin(input_parameters[11]) + (supposed_y_0 + input_parameters[10]) * cos(input_parameters[11])
    supposed_x = supposed_x_0 * cos(input_parameters[11]) + supposed_y_0 * sin(input_parameters[11]) + input_parameters[
        9]
    supposed_y = -supposed_x_0 * sin(input_parameters[11]) + supposed_y_0 * cos(input_parameters[11]) + \
                 input_parameters[10]
    # print(supposed_x, supposed_y)
    deviation_x = supposed_x - input_x  # Deviation between distorted x and non-distorted x
    deviation_y = supposed_y - input_y  # Deviation between distorted y and non-distorted y
    frame_centre_x = frame_centre_init_x
    frame_centre_y = frame_centre_init_y
    # print(frame_centre_x, frame_centre_y)
    centered_x = input_x - frame_centre_x  # x as measured from the center of the full frame
    centered_y = input_y - frame_centre_y  # y as measured from the center of the full frame
    # print(centered_x, centered_y)
    centered_radius = arr_sqrt(centered_x ** 2 + centered_y ** 2)  # Distance from the center of the frame (distorted)
    # Additional variables
    a_k = centered_radius ** 2 + 2 * centered_x ** 2
    b_k = 2 * centered_x * centered_y
    c_k = centered_radius ** 2 + 2 * centered_y ** 2
    # print(a_k, b_k, c_k)
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
        # print(deviation_x[min_ind], deviation_y[min_ind], f_radial, f_tangential, g_radial, g_tangential, f_fun, g_fun)
    # print('Sum of squared deviations from distortion: ' + "{}".format(min_value))
    return min_value

