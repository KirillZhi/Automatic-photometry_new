import pandas as pd
from numpy import zeros, empty, append, arange, max, absolute, array, delete, count_nonzero, argmax, reshape, unique
from numpy import min as arr_min
from numpy import max as arr_max

from pandas import DataFrame
from itertools import combinations

from math import sqrt
from frame_functions import *
from coord_functions import *

from astropy.io.fits import getdata, open
from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename

from os import listdir, remove
from os.path import isfile, join

from re import search
from matplotlib import pyplot as plt
import time
import dateutil
import datetime
import openpyxl


# Reading a datafile
print("Choose data with original excel")
data_path = askopenfilename(initialdir="D:/")
data = pd.read_excel(data_path)
data_columns = data.columns.values.tolist()
print(data[data_columns[0]])


print("Choose file with corrected data")
complete_data_path = askopenfilename(initialdir="D:/")
complete_data = pd.read_excel(complete_data_path)

# Estimating syst error thru least squares
syst_err_arr = np.arange(0., 5., 0.1)
best_syst_err_arr = 0
min_sq_sum = 10 ** 6
for syst_ind in range(len(syst_err_arr)):
    new_source_dist = np.zeros_like(source_dist)
    new_one_sigma_dist = np.zeros_like(source_dist)
    for i in range(len(source_dist)):
        new_one_sigma_dist[i] = math.sqrt(one_sigma_dist[i] ** 2 + syst_err_arr[syst_ind] ** 2)
    for i in range(len(source_dist)):
        new_source_dist[i] = source_dist[i] / new_one_sigma_dist[i]
    new_source_dist = np.sort(new_source_dist)
    if syst_ind == 0:
        print(new_source_dist)
    # Making empirical cumulative distribution function
    true_distrib_fun = np.zeros_like(distribution_function)
    cum_distr_fun = np.zeros_like(distribution_function)
    distribution_function_ninety_perc_conf = np.zeros_like(distribution_function)
    for i in range(len(distribution_function)):
        value_num = 0.0
        for j in range(len(new_source_dist)):
            if r_values.loc[i] >= new_source_dist[j]:
                value_num += 1
        cum_distr_fun[i] = (value_num) / len(one_sigma_dist)
        distribution_function_ninety_perc_conf[i] = math.sqrt(value_num) / len(one_sigma_dist) * 1.65
        #print(distribution_function_ninety_perc_conf[i], i + syst_err_arr[syst_ind])
    #print(true_distrib_fun)
    #line1 = plt.plot(r_values, distribution_function, c='r')
    #print(new_source_dist, cum_distr_fun)
    if syst_ind == 100:
        print(new_source_dist)
        #line1 = plt.plot(new_source_dist, true_distrib_fun, c='r')
        #line2 = plt.plot(new_source_dist, cum_distr_fun, c='b')
        line2 = plt.plot(r_values, cum_distr_fun, c='b')
        line3 = plt.plot(r_values, distribution_function, c='green')
        plt.fill_between(r_values, cum_distr_fun - distribution_function_ninety_perc_conf,
                         cum_distr_fun + distribution_function_ninety_perc_conf, color='b', alpha=.1)
        plt.figure(1)
        plt.show()
        plt.close()
    #time.sleep(2)
    # Least squares
    #print(true_distrib_fun)
    #LS_sum = np.sum(((cum_distr_fun - true_distrib_fun) / distribution_function_ninety_perc_conf * 1.65) ** 2)
    LS_sum = np.sum(((cum_distr_fun - distribution_function) / 1) ** 2)
    print(LS_sum)
    if LS_sum <= min_sq_sum:
        best_syst_err_arr = syst_err_arr[syst_ind]
        min_sq_sum = LS_sum
        print("Best guess of syst error: " + "{}".format(best_syst_err_arr) + ", sum of dev: " + "{}".format(min_sq_sum))
    else:
        break
# Estimating error of syst error
syst_err = 1.5
new_source_dist = np.zeros_like(source_dist)
test_source_dist = np.zeros_like(source_dist)
test_one_sigma_dist = np.zeros_like(source_dist)
for i in range(len(source_dist)):
    new_one_sigma_dist[i] = math.sqrt(one_sigma_dist[i] ** 2 + syst_err ** 2)
    test_one_sigma_dist[i] = math.sqrt(one_sigma_dist[i] ** 2 + (syst_err + 0.01) ** 2)
for i in range(len(source_dist)):
    new_source_dist[i] = source_dist[i] / new_one_sigma_dist[i]
    test_source_dist[i] = source_dist[i] / test_one_sigma_dist[i]


cum_distr_fun = np.zeros_like(distribution_function)
test_cum_distr_fun = np.zeros_like(distribution_function)
distribution_function_ninety_perc_conf = np.zeros_like(distribution_function)
df_fun = np.zeros_like(distribution_function)
for i in range(len(distribution_function)):
    value_num = 0.0
    test_value_num = 0.0
    for j in range(len(new_source_dist)):
        if r_values.loc[i] >= new_source_dist[j]:
            value_num += 1
        if r_values.loc[i] >= test_source_dist[j]:
            test_value_num += 1
    cum_distr_fun[i] = (value_num) / len(one_sigma_dist)
    test_cum_distr_fun[i] = (test_value_num) / len(one_sigma_dist)
    df_fun[i] = (test_cum_distr_fun[i] - cum_distr_fun[i]) / 0.01
#print(df_fun)
print(np.dot(np.transpose(df_fun), df_fun))
syst_err_variance = min_sq_sum / (np.dot(np.transpose(df_fun), df_fun))
print("Systematic error variance: " + "{}".format(syst_err_variance))
print(math.sqrt(syst_err_variance))

print("End of calculations")
print("Best guess of syst error: " + "{}".format(best_syst_err_arr) + ", sum of dev: " + "{}".format(min_sq_sum))
line1 = plt.plot(r_values, distribution_function, c='r')
line2 = plt.plot(r_values, cum_distr_fun, c='b')
plt.fill_between(r_values, cum_distr_fun - distribution_function_ninety_perc_conf, cum_distr_fun + distribution_function_ninety_perc_conf, color='b', alpha = .1)
plt.show()
plt.close()
#print(one_sigma_dist)