from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename
import numpy as np
import math
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import dateutil
import datetime


start_point_str = '2020-01-07 19:26:22'         # t0 for plots such as GRBs
# Reading a datafile
print('Picking data for the first plot')
end_check=0
line_ind = 0
label_arr = []
color_arr = ['red', 'blue', 'green', 'orange', 'grey', 'saddlebrown', 'blueviolet', 'cyan']
font_small = 14
font_medium = 20
font_large = 24
#plt.rcParams.update({'font.size': 22})
#plt.rc('font', size = font_small)
plt.rc('axes', titlesize = font_medium)
plt.rc('axes', labelsize = font_medium)
plt.rc('xtick', labelsize = font_small)
plt.rc('ytick', labelsize = font_small)
plt.rc('legend', fontsize = font_small)
plt.rc('figure', titlesize = font_large)

comb_mag_data = np.zeros(0)
comb_date_data = np.zeros(0)
comb_mag_err_data = np.zeros(0)
while end_check == 0:
    if line_ind == 0:
        data_path = askopenfilename(initialdir="F:/Astronomy data/")
    else:
        data_path = askopenfilename(initialdir=data_path)

    data = pd.read_excel(data_path)
    data = data.sort_values(by='Julian Days, middle of exposure')
    # Discarding points with a large error (signal is less than 3 sigma
    dict_data = []
    array_columns = data.columns.values.tolist()
    for j_ind in range(len(data['Source mag (corrected)'])):
        if abs(float(data['Source mag error (signal noise)'].values[j_ind])) < 0.114 and abs(float(data['Source mag error (variability, corrected)'].values[j_ind])) < 0.114:
            dict1 = dict((col_name, data[col_name].values[j_ind]) for col_name in array_columns)
            dict_data.append(dict1)
    data = pd.DataFrame(dict_data, columns=array_columns)
    #print(data)
    #data = data.drop(np.where(data['Source mag error (signal noise)'] >= 0.31))
    data = data.reset_index(drop=True)
    mag_data = data['Source mag (corrected)']
    date_data = data['Julian Days, middle of exposure']
    exp_data = data['Frame exposure time, seconds']
    mag_err = np.sqrt(data['Source mag error (variability, corrected)'] ** 2 + data['Source mag error (signal noise)'] ** 2)

    comb_mag_data = np.append(comb_mag_data, mag_data)
    comb_date_data = np.append(comb_date_data, date_data)
    comb_mag_err_data = np.append(comb_mag_err_data, mag_err)
    end_check = int(raw_input('Do you want to end picking data? 0 to continue, everything else to end: '))

# Making correlated pairs
pair_time_gap_data = np.zeros(0)
ucdf_pair_data = np.zeros(0)

mag_average = np.average(comb_mag_data)
mag_stdev = np.std(comb_mag_data)
print(mag_stdev)
for i_ind in range(len(comb_mag_data)):
    if i_ind % 100 == 0:
        print(i_ind)
    #if i_ind == 500:
    #    break
    for j_ind in range(len(comb_mag_data)):
        if i_ind == j_ind:
            continue
        ucdf_value = (comb_mag_data[i_ind] - mag_average) * (comb_mag_data[j_ind] - mag_average) / np.sqrt((mag_stdev ** 2 - comb_mag_err_data[i_ind] ** 2) * (mag_stdev ** 2 - comb_mag_err_data[j_ind] ** 2))
        ucdf_pair_data = np.append(ucdf_pair_data, ucdf_value)
        time_gap = abs(comb_date_data[i_ind] - comb_date_data[j_ind]) * 24
        pair_time_gap_data = np.append(pair_time_gap_data, time_gap)
df = {"UCDF": ucdf_pair_data, "Time gap": pair_time_gap_data}
ucdf_pairs = pd.DataFrame(data=df)
print(ucdf_pairs)
#ucdf_pairs = ucdf_pairs.sort_values(by='Time gap')
#print(ucdf_pairs)
# Averaging UCDF's
bin = 1.0         # in hours
time_lag = 1    # in hours
max_time_lag = int(np.max(pair_time_gap_data / 4)) # in hours
print(max_time_lag)
tau_arr = np.zeros(0)
dcf_arr = np.zeros(0)
dcf_err_arr = np.zeros(0)
pair_bin_num = np.zeros(0)
for tau_ind in range(0, max_time_lag, time_lag):
    bin_num = 0
    dcf_value = 0
    dcf_err = 0
    ucdf_values = np.zeros(0)
    for pair_ind in range(len(pair_time_gap_data)):
        if pair_ind % 10000 == 0:
            print(tau_ind, pair_ind, pair_time_gap_data[pair_ind], tau_ind - bin / 2, tau_ind + bin / 2)
        if ucdf_pairs['Time gap'][pair_ind] >= tau_ind - bin/2 and ucdf_pairs['Time gap'][pair_ind] < tau_ind + bin/2:
            #print('ooga booga')
            bin_num += 1
            dcf_value += ucdf_pairs['UCDF'][pair_ind]
            ucdf_values = np.append(ucdf_values, ucdf_pairs['UCDF'][pair_ind])
    if bin_num == 0 or bin_num == 1:
        continue
    dcf_value = dcf_value / bin_num
    for err_ind in range(len(ucdf_values)):
        dcf_err += (ucdf_values[err_ind] - dcf_value) ** 2
    dcf_err = math.sqrt(dcf_err) / (bin_num - 1)
    dcf_arr = np.append(dcf_arr, dcf_value)
    dcf_err_arr = np.append(dcf_err_arr, dcf_err)
    tau_arr = np.append(tau_arr, float(tau_ind) / 24)
# Making a plot
plt.figure(2)
continuation_check = 1
dcf_df = {"ACF": dcf_arr, "ACF_err": dcf_err_arr, "Tau (d)": tau_arr}
dcf_dataframe = pd.DataFrame(data = dcf_df)
print('Choose where to save autocorrelation data')
dataframe_path = askdirectory(initialdir=data_path)
dcf_dataframe.to_excel("{}".format(dataframe_path) + "auto_fun_data.xlsx")
print(dcf_arr)
#label_string = label_string[:len(label_string) - 1]
line_1 = plt.errorbar(tau_arr, dcf_arr, yerr=dcf_err_arr, marker='o', c=color_arr[line_ind], ms=10, ls='none', elinewidth=3, capsize=6)
# Handling titles, legend, axes of the plot
plt.legend()
min_x = -1
max_x = 1

plt.ylabel('Magnitude')
title_string = raw_input('Title of the plot: ')
#title_string = title_string[:len(title_string) - 1]
plt.title(title_string)
#xlabel_string = raw_input('Name of the GRB gamma-detector: ')
#xlabel_string = xlabel_string[:len(xlabel_string) - 1]

plt.xlabel("Tau (d)")
#plt.xscale('log')
#plt.gca().invert_yaxis()

plt.show()
plt.close()