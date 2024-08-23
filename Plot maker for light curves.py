from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import dateutil
import datetime
start_point_str = '2023-02-04 21:47:51'         # t0 for plots such as GRBs
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
        if abs(float(data['Source mag error (signal noise)'].values[j_ind])) < 0.314 and abs(float(data['Source mag error (variability of bright checks, corrected)'].values[j_ind])) < 0.314:
            dict1 = dict((col_name, data[col_name].values[j_ind]) for col_name in array_columns)
            dict_data.append(dict1)
    data = pd.DataFrame(dict_data, columns=array_columns)
    #print(data)
    #data = data.drop(np.where(data['Source mag error (signal noise)'] >= 0.31))
    data = data.reset_index(drop=True)
    mag_data = data['Source mag (corrected)']
    date_data = data['Julian Days, middle of exposure']
    exp_data = data['Frame exposure time, seconds']
    mag_err = np.sqrt(data['Source mag error (variability of bright checks, corrected)'] ** 2 + data['Source mag error (signal noise)'] ** 2)

    # Part for GRB lc's

    start_datetime = dateutil.parser.parse(start_point_str)
    start_JD_date = start_datetime.toordinal()
    start_JD_date_datetime = datetime.datetime.fromordinal(start_JD_date)
    print(start_datetime)
    print(start_JD_date_datetime)
    # Correcting dates
    corrected_date_data = np.zeros_like(date_data)
    for j in range(len(corrected_date_data)):
        corrected_date_data[j] = (date_data[j] - start_JD_date - 1721424.5) * 86400 - (start_datetime - start_JD_date_datetime).total_seconds()
    # Making a plot
    plt.figure(2)
    continuation_check = 1
    while (1):
        label_string = raw_input('Name for the plot: ')
        continuation_check = int(input('Do you want to change name of the plot? 1 to change, anything else to not: '))
        if continuation_check != 1:
            break
    label_string = label_string[:len(label_string)]
    line_1 = plt.errorbar(corrected_date_data, mag_data, yerr=mag_err, xerr=exp_data / 2, marker='o', c=color_arr[line_ind], ms=3, ls='none', label=label_string)
    if line_ind == 0:
        min_x = np.min(corrected_date_data)
        max_x = np.max(corrected_date_data)
    else:
        if min_x > np.min(corrected_date_data):
            min_x = np.min(corrected_date_data)
        if max_x < np.max(corrected_date_data):
            max_x = np.max(corrected_date_data)
    line_ind += 1
    end_check = int(raw_input('Do you want to end picking data? 0 to continue, everything else to continue: '))

# Handling titles, legend, axes of the plot
plt.legend()
min_x = min_x * 0.8
max_x = max_x * 1.2

plt.ylabel('Magnitude')
title_string = raw_input('Title of the plot: ')
title_string = title_string[:len(title_string)]
plt.title(title_string)
xlabel_string = raw_input('Name of the GRB gamma-detector: ')
xlabel_string = xlabel_string[:len(xlabel_string)]

plt.xlabel("Time since " + "{}".format(xlabel_string) + " trigger, seconds")
plt.xscale('log')
plt.gca().invert_yaxis()

plt.show()
plt.close()