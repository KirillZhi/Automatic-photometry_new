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
from mpl_toolkits.mplot3d import Axes3D
from time import time
import dateutil
import datetime
import openpyxl
# Parameters of the initial fit (parameters of anchors)




# Filepath for fits frames

root = Tk()
root.withdraw()
print("Choose file with photometry")
print('')
photopath = askopenfilename(initialdir="D:/Photometry/")
print("Choose file with check star data")
check_path = askopenfilename(initialdir=photopath)

check_photo = pd.read_csv(photopath)
check_arr = pd.read_excel(check_path)
#print(check_photo)
#print(check_photo.columns)
check_photo_arr = np.zeros(0)
check_photo_err_arr = np.zeros_like(check_photo_arr)
check_b = np.zeros(0)
check_r = np.zeros(0)
#column_name = 'Check star ' + j + ' mag'
for i in range(len(check_photo['Frame'])):
    for j in range(1, len(check_arr['Gmag']) + 1):
        #print(j)
        column_name = 'Check star ' + "{}".format(j) + 'corrected mag'
        #print(column_name)
        if abs(check_photo[column_name].values[i] - check_arr['B1mag'].values[j - 1]) <= 2 and abs(check_photo[column_name].values[i] - check_arr['R1mag'].values[j - 1]) <= 2:
            check_photo_arr = np.append(check_photo_arr, check_photo[column_name].values[i])
            column_name = 'Check star ' + "{}".format(j) + ' mag error (signal noise)'
            check_photo_err_arr = np.append(check_photo_err_arr, check_photo[column_name].values[i])
            check_b = np.append(check_b, check_arr['B1mag'].values[j - 1])
            check_r = np.append(check_r, check_arr['R1mag'].values[j - 1])
print(check_photo_arr)

#check_b = check_arr['B1mag'] - check_arr['B1mag'].values[0]
#check_r = check_arr['R1mag'] - check_arr['R1mag'].values[0]
#check_b = check_arr['BPmag']
print(len(check_photo_arr))
check_len = len(check_b)
photo_len = len(check_photo['Frame'])
check_photo_arr = np.reshape(check_photo_arr, (photo_len, check_len))
check_photo_err_arr = np.reshape(check_photo_err_arr, (photo_len, check_len))
#check_photo_arr -= check_r
#check_b = check_b - check_r
full_param_num = 2
g_check = 1
# Defining arrays
check_len = len(check_b)
photo_len = len(check_photo['Frame'])
alpha = 10 ** -3
x_arr = np.zeros((check_len * photo_len, full_param_num))
coord_var_mag = np.zeros(check_len * photo_len)
weight_arr = np.zeros((check_len * photo_len, check_len * photo_len))
# Fit parameters
coord_par = np.zeros(full_param_num)
coord_par_err = np.zeros(full_param_num)
# Main body
for k in range(photo_len):
    for i in range(check_len):
        if full_param_num == 3:
            x_arr[i + check_len * k][0] = 1
            x_arr[i + check_len * k][1] = check_b[i]
            x_arr[i + check_len * k][2] = check_r[i]
        if full_param_num == 2 and g_check != 1:
            x_arr[i + check_len * k][0] = check_b[i]
            x_arr[i + check_len * k][1] = check_r[i]
        if full_param_num == 2 and g_check == 1:
            x_arr[i + check_len * k][0] = 1
            x_arr[i + check_len * k][1] = check_b[i]
        coord_var_mag[i + check_len * k] = check_photo_arr[k][i]
        #print(i + check_len * k, check_len * photo_len)
        weight_arr[i + check_len * k][i + check_len * k] = 1 / check_photo_err_arr[k][i] ** 2
    # print(x_arr)
x_arr_transpose = np.transpose(x_arr)
x_mod_arr = np.dot(x_arr_transpose, weight_arr)
y_mod_arr = np.dot(x_mod_arr, coord_var_mag)
x_mod_arr = np.dot(x_mod_arr, x_arr)
x_inv_arr = np.linalg.inv(x_mod_arr + alpha * np.identity(full_param_num))
vector_param = np.dot(x_inv_arr, y_mod_arr)
for m in range(0, full_param_num):
    coord_par[m] = vector_param[m]
print(coord_par)
# Error estimation
sum_arr = coord_var_mag - np.dot(x_arr, vector_param)
sum_arr_trans = np.transpose(sum_arr)
minim_sum = np.dot(np.dot(sum_arr_trans, weight_arr), sum_arr)
parameter_variance = minim_sum / (check_len - full_param_num) * x_inv_arr
for m in range(0, full_param_num):
    if parameter_variance[m][m] >= 0:
        coord_par_err[m] = math.sqrt(parameter_variance[m][m])
    else:
        coord_par_err[m] = -math.sqrt(-parameter_variance[m][m])
print(coord_par_err)
print('Sum of squared deviations of magnitude (last frame): ' + "{}".format(minim_sum))
for k in range(photo_len):
    for i in range(check_len):
        weight_arr[i + check_len * k][i + check_len * k] = 1

x_arr_transpose = np.transpose(x_arr)
x_mod_arr = np.dot(x_arr_transpose, weight_arr)
y_mod_arr = np.dot(x_mod_arr, coord_var_mag)
x_mod_arr = np.dot(x_mod_arr, x_arr)
x_inv_arr = np.linalg.inv(x_mod_arr + alpha * np.identity(full_param_num))
vector_param = np.dot(x_inv_arr, y_mod_arr)
for m in range(0, full_param_num):
    coord_par[m] = vector_param[m]
print(coord_par)
# Error estimation
sum_arr = coord_var_mag - np.dot(x_arr, vector_param)
sum_arr_trans = np.transpose(sum_arr)
minim_sum = np.dot(np.dot(sum_arr_trans, weight_arr), sum_arr)
parameter_variance = minim_sum / (check_len - full_param_num) * x_inv_arr
for m in range(0, full_param_num):
    if parameter_variance[m][m] >= 0:
        coord_par_err[m] = math.sqrt(parameter_variance[m][m])
    else:
        coord_par_err[m] = -math.sqrt(-parameter_variance[m][m])
print(coord_par_err)
#coord_par[1] = 0.11
print('Sum of squared deviations of magnitude (last frame): ' + "{}".format(minim_sum))
if full_param_num == 2 and g_check == 1:
    b_min = np.min(check_b)
    b_max = np.max(check_b)

    #coord_par[0] = 0.2
    #coord_par[1] = 0.8
    b_arr = np.linspace(b_min, b_max, 30)

    c_arr = np.zeros(len(b_arr))

    for i in range(len(b_arr)):
        c_arr[i] = coord_par[1] * b_arr[i] + coord_par[0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(b_arr, c_arr, alpha=0.5)

    ax.scatter(check_b, check_photo_arr)

    ax.set_xlabel('B1mag')
    ax.set_ylabel('MASTER mag')

    plt.show()
if (full_param_num == 2 and g_check !=1) or full_param_num == 3:
    b_min = np.min(check_b)
    b_max = np.max(check_b)
    r_min = np.min(check_r)
    r_max = np.max(check_r)
    #coord_par[0] = 0.2
    #coord_par[1] = 0.8
    b_arr = np.linspace(b_min, b_max, 30)
    r_arr = np.linspace(r_min, r_max, 30)

    c_arr = np.zeros((len(b_arr), len(r_arr)))

    for i in range(len(b_arr)):
        for j in range(len(r_arr)):
            if full_param_num == 2:
                c_arr[i][j] = coord_par[0] * b_arr[i] + coord_par[1] * r_arr[j]
            if full_param_num == 3:
                c_arr[i][j] = coord_par[1] * b_arr[i] + coord_par[2] * r_arr[j] + coord_par[0]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X,Y=np.meshgrid(b_arr,r_arr)
    ax.scatter(X, Y, c_arr, alpha=0.5)

    ax.scatter(check_b, check_r, check_photo_arr)

    ax.set_xlabel('B1mag')
    ax.set_ylabel('R1mag')
    ax.set_zlabel('MASTER mag')

    plt.show()
