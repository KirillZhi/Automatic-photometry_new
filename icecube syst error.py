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
print("Choose data with source distances")
data_path = askopenfilename(initialdir="D:/")
data = pd.read_excel(data_path)
print(data)
# 1.645 0.68, 2.4
# 1.66, 2.8
data_cols = data.columns.values.tolist()
#data = data.
one_sigma_dist = data[data_cols[6]] / 2.145
source_dist = data[data_cols[7]]
source_dist = np.append(source_dist, 0)
one_sigma_dist = np.append(one_sigma_dist, 99)
print(source_dist)
print("Choose file with gauss distribution")
catalog_path = askopenfilename(initialdir="D:/")
gauss_distr_data = pd.read_excel(catalog_path)
distribution_function = gauss_distr_data['Distribution function']
r_values = gauss_distr_data['Equalized radius']

# Estimating syst error thru least squares
syst_err_arr = np.arange(0, 1.6, 0.05)
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
    # Making empirical cumulative distribution function
    true_distrib_fun = np.zeros_like(one_sigma_dist)
    cum_distr_fun = np.zeros_like(one_sigma_dist)
    distribution_function_ninety_perc_conf = np.zeros_like(one_sigma_dist)
    for i in range(len(one_sigma_dist)):
        cum_distr_fun[i] = float(i) / (len(one_sigma_dist) - 1)
        distribution_function_ninety_perc_conf[i] = math.sqrt(i) / (len(one_sigma_dist) - 1)
        value_num = 0
        for j in range(len(distribution_function)):
            if i != len(one_sigma_dist) - 1 and i != 0:
                if r_values.loc[j] < (new_source_dist[i] + new_source_dist[i + 1]) / 2 and r_values.loc[j] > (new_source_dist[i] + new_source_dist[i - 1]) / 2:
                    true_distrib_fun[i] += distribution_function[j]
                    value_num += 1
                if r_values.loc[j] > (new_source_dist[i] + new_source_dist[i + 1]) / 2:
                    break
            if i == 0:
                if r_values.loc[j] < (new_source_dist[i] + new_source_dist[i + 1]) / 2 and r_values.loc[j] >= new_source_dist[i] - (new_source_dist[i + 1] - new_source_dist[i]) / 2:
                    true_distrib_fun[i] += distribution_function[j]
                    value_num += 1
                if r_values.loc[j] > (new_source_dist[i] + new_source_dist[i + 1]) / 2:
                    break
            if i == len(one_sigma_dist) - 1:
                if r_values.loc[j] < new_source_dist[i] + (new_source_dist[i] - new_source_dist[i - 1]) / 2 and r_values.loc[j] >= (new_source_dist[i - 1] + new_source_dist[i]) / 2:
                    true_distrib_fun[i] += distribution_function[j]
                    value_num += 1
                if r_values.loc[j] > new_source_dist[i] + (new_source_dist[i] - new_source_dist[i - 1]) / 2:
                    break
        true_distrib_fun[i] = true_distrib_fun[i] / value_num

        #print(distribution_function_ninety_perc_conf[i], i + syst_err_arr[syst_ind])
    #line1 = plt.plot(r_values, distribution_function, c='r')
    #print(new_source_dist, cum_distr_fun)
    if syst_ind == 0:
        line1 = plt.plot(r_values, distribution_function, c='r', label='Model distribution')
        line2 = plt.plot(new_source_dist, cum_distr_fun, c='b', label='Empirical distribution')
        line4 = plt.fill_between(new_source_dist, cum_distr_fun - distribution_function_ninety_perc_conf,
                                 cum_distr_fun + distribution_function_ninety_perc_conf, color='b', alpha=.1,
                                 label='1 sigma error')
        plt.xlabel('r normalised')
        plt.ylabel('p')
        # plt.legend([line1, line2, line3])
        plt.legend()
        plt.ylim([0, 1.05])
        plt.show()
        plt.close()
    #time.sleep(2)
    # Least squares
    #print(true_distrib_fun)
    #LS_sum = np.sum(((cum_distr_fun - true_distrib_fun) / distribution_function_ninety_perc_conf * 1.65) ** 2)
    LS_sum = np.sum(((cum_distr_fun - true_distrib_fun) / 1) ** 2)
    print(LS_sum)
    if LS_sum <= min_sq_sum:
        best_syst_err_arr = syst_err_arr[syst_ind]
        min_sq_sum = LS_sum
        print("Best guess of syst error: " + "{}".format(best_syst_err_arr) + ", sum of dev: " + "{}".format(min_sq_sum))
# Estimating error of syst error
syst_err = 1.4
new_source_dist = np.zeros_like(source_dist)
new_one_sigma_dist = np.zeros_like(source_dist)
for i in range(len(source_dist)):
    new_one_sigma_dist[i] = math.sqrt(one_sigma_dist[i] ** 2 + syst_err ** 2)
for i in range(len(source_dist)):
    new_source_dist[i] = source_dist[i] / new_one_sigma_dist[i]
new_source_dist = np.sort(new_source_dist)
# Making empirical cumulative distribution function
true_distrib_fun = np.zeros_like(one_sigma_dist)
cum_distr_fun = np.zeros_like(one_sigma_dist)
distribution_function_ninety_perc_conf = np.zeros_like(one_sigma_dist)
for i in range(len(one_sigma_dist)):
    cum_distr_fun[i] = float(i) / (len(one_sigma_dist) - 1)
    distribution_function_ninety_perc_conf[i] = math.sqrt(i) / (len(one_sigma_dist) - 1)
    value_num = 0
    for j in range(len(distribution_function)):
        if i != len(one_sigma_dist) - 1 and i != 0:
            if r_values.loc[j] < (new_source_dist[i] + new_source_dist[i + 1]) / 2 and r_values.loc[j] > (new_source_dist[i] + new_source_dist[i - 1]) / 2:
                true_distrib_fun[i] += distribution_function[j]
                value_num += 1
            if r_values.loc[j] > (new_source_dist[i] + new_source_dist[i + 1]) / 2:
                break
        if i == 0:
            if r_values.loc[j] < (new_source_dist[i] + new_source_dist[i + 1]) / 2 and r_values.loc[j] >= new_source_dist[i] - (new_source_dist[i + 1] - new_source_dist[i]) / 2:
                true_distrib_fun[i] += distribution_function[j]
                value_num += 1
            if r_values.loc[j] > (new_source_dist[i] + new_source_dist[i + 1]) / 2:
                break
        if i == len(one_sigma_dist) - 1:
            if r_values.loc[j] < new_source_dist[i] + (new_source_dist[i] - new_source_dist[i - 1]) / 2 and r_values.loc[j] >= (new_source_dist[i - 1] + new_source_dist[i]) / 2:
                true_distrib_fun[i] += distribution_function[j]
                value_num += 1
            if r_values.loc[j] > new_source_dist[i] + (new_source_dist[i] - new_source_dist[i - 1]) / 2:
                break
    true_distrib_fun[i] = true_distrib_fun[i] / value_num
# Second variant
syst_err = 1.5
test_source_dist = np.zeros_like(source_dist)
test_one_sigma_dist = np.zeros_like(source_dist)
for i in range(len(source_dist)):
    test_one_sigma_dist[i] = math.sqrt(one_sigma_dist[i] ** 2 + syst_err ** 2)
for i in range(len(source_dist)):
    test_source_dist[i] = source_dist[i] / test_one_sigma_dist[i]
test_source_dist = np.sort(test_source_dist)
# Making empirical cumulative distribution function
test_true_distrib_fun = np.zeros_like(one_sigma_dist)
test_cum_distr_fun = np.zeros_like(one_sigma_dist)
test_distribution_function_ninety_perc_conf = np.zeros_like(one_sigma_dist)
for i in range(len(one_sigma_dist)):
    test_cum_distr_fun[i] = float(i) / (len(one_sigma_dist) - 1)
    test_distribution_function_ninety_perc_conf[i] = math.sqrt(i) / (len(one_sigma_dist) - 1)
    value_num = 0
    for j in range(len(distribution_function)):
        if i != len(one_sigma_dist) - 1 and i != 0:
            if r_values.loc[j] < (test_source_dist[i] + test_source_dist[i + 1]) / 2 and r_values.loc[j] > (test_source_dist[i] + test_source_dist[i - 1]) / 2:
                test_true_distrib_fun[i] += distribution_function[j]
                value_num += 1
            if r_values.loc[j] > (test_source_dist[i] + test_source_dist[i + 1]) / 2:
                break
        if i == 0:
            if r_values.loc[j] < (test_source_dist[i] + test_source_dist[i + 1]) / 2 and r_values.loc[j] >= test_source_dist[i] - (test_source_dist[i + 1] - test_source_dist[i]) / 2:
                test_true_distrib_fun[i] += distribution_function[j]
                value_num += 1
            if r_values.loc[j] > (test_source_dist[i] + test_source_dist[i + 1]) / 2:
                break
        if i == len(one_sigma_dist) - 1:
            if r_values.loc[j] < test_source_dist[i] + (test_source_dist[i] - test_source_dist[i - 1]) / 2 and r_values.loc[j] >= (test_source_dist[i - 1] + test_source_dist[i]) / 2:
                test_true_distrib_fun[i] += distribution_function[j]
                value_num += 1
            if r_values.loc[j] > test_source_dist[i] + (test_source_dist[i] - test_source_dist[i - 1]) / 2:
                break
    test_true_distrib_fun[i] = test_true_distrib_fun[i] / value_num
df_fun = np.zeros_like(test_true_distrib_fun)
for i in range(len(test_true_distrib_fun)):
    df_fun[i] = (test_true_distrib_fun[i] - true_distrib_fun[i]) / 0.1
    #df_fun[i] = (test_cum_distr_fun[i] - cum_distr_fun[i]) / 0.1
print(df_fun)
syst_err_variance = min_sq_sum / (np.dot(np.transpose(df_fun), df_fun))
print("Systematic error variance: " + "{}".format(syst_err_variance))
print(math.sqrt(syst_err_variance))
print("End of calculations")
print("Best guess of syst error: " + "{}".format(best_syst_err_arr) + ", sum of dev: " + "{}".format(min_sq_sum))
line1 = plt.plot(r_values, distribution_function, c='r', label='Model distribution')
line2 = plt.plot(new_source_dist, cum_distr_fun, c='b', label='Empirical distribution')
line3 = plt.plot(new_source_dist, test_true_distrib_fun, c='g', label='Binned model distribution')
line4 = plt.fill_between(new_source_dist, cum_distr_fun - test_distribution_function_ninety_perc_conf, cum_distr_fun + test_distribution_function_ninety_perc_conf, color='b', alpha = .1, label='1 sigma error')
plt.xlabel('r normalised')
plt.ylabel('p')
#plt.legend([line1, line2, line3])
plt.legend()
plt.ylim([0, 1.05])
plt.show()
plt.close()
#print(one_sigma_dist)