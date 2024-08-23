from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import dateutil
import datetime
import xlrd
import openpyxl
from os import listdir, remove, mkdir
from os.path import isfile, join, isdir
from re import search
import time
line_point_str = '2021-12-08 20:02:51.1'         # t0 for plots such as GRBs
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
data_path = askdirectory(initialdir="D:/Astronomy data/")
for file_path in listdir(data_path):
    print(data_path + file_path)
    x = search("job", file_path)
    if x:
        data = pd.read_csv(data_path + '/' + file_path, sep='\s+')
    else:
        continue

    # Discarding points with a large error (signal is less than 3 sigma)
    dict_data = []
    dict_limit = []
    array_columns = data.columns.values.tolist()
    for j_ind in range(len(data['m'])):
        if data['uJy'].values[j_ind] > 0.0 and data['uJy'].values[j_ind]/data['duJy'].values[j_ind] >= 5.0 and data['F'].values[j_ind] == 'o' and data['Sky'].values[j_ind] >= 16.0:
            dict1 = dict((col_name, data[col_name].values[j_ind]) for col_name in array_columns)
            dict_data.append(dict1)

    data = pd.DataFrame(dict_data, columns=array_columns)
    data = data.reset_index(drop=True)
    mag_data = data['m']
    date_data = data['###MJD']
    mag_err = data['dm']

    # Correcting dates
    # Making a plot
    plt.figure(2)
    #plt.scatter(date_limit_data, limit_data, marker='v', s=100, c='black')

    line_1 = plt.errorbar(date_data, mag_data, yerr=mag_err, marker='o', c=color_arr[line_ind], ms=5, ls='none', elinewidth=1, capsize=5)

    if line_ind == 0:
        min_x = np.min(date_data)
        max_x = np.max(date_data)
    else:
        if min_x > np.min(date_data):
            min_x = np.min(date_data)
        if max_x < np.max(date_data):
            max_x = np.max(date_data)

    plt.legend()
    min_x = min_x * 0.8
    max_x = max_x * 1.2

    plt.ylabel('Magnitude')
    title_string = str(data['RA'].values[0]) + '_' + str(data['Dec'].values[0])
    plt.title(title_string)

    plt.xlabel("MJD")
    plt.gca().invert_yaxis()
    #plt.ion()
    plt.savefig(data_path + '/' + title_string + '.png')
    data.to_excel(data_path + '/' + title_string + '.xlsx')
    plt.clf()

# Handling titles, legend, axes of the plot
