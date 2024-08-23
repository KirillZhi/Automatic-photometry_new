import pandas as pd

from frame_functions import *
from astropy.io import fits
from Tkinter import Tk
from tkFileDialog import askdirectory

import os
import re
from matplotlib import pyplot as plt

# Filepath for fits frames


root = Tk()
root.withdraw()
filepath = askdirectory(initialdir="C:/Users/kirio/Desktop/astronomy/GRB/")
photometry_path = askdirectory(initialdir="C:/Users/kirio/Desktop/astronomy/GRB/")
filepath = filepath + "/"
photometry_path = photometry_path + "/"
root.update()
root.destroy()

print(filepath)
print(photometry_path)

# Search of fits frames in a folder
photo_frame_path = np.array([])
res = []
for path in os.listdir(filepath):
    if os.path.isfile(os.path.join(filepath, path)):
        res.append(path)
for i in res:
    # print(i)
    x = re.search(".*fit$", i)
    y = re.search(".*fits$", i)
    if x or y:
        photo_frame_path = np.append(photo_frame_path, i)
print(photo_frame_path)

# Reading of a test fits frame
filename = photo_frame_path[0]
check_path = filepath + filename
check_frame = fits.getdata(check_path)

check_shape = check_frame.shape
check_mid_x = check_shape[0] / 2
check_mid_y = check_shape[1] / 2
# Minimum value of the frame (positive)
test_ground = np.where(check_frame > 0, check_frame, np.max(check_frame))
c_min = np.percentile(test_ground, 5)
print(c_min)

# Search of reference and check stars
plt.ion()
plt.figure(1)
im = plt.imshow(check_frame)
plt.colorbar(im)
# c_min=20
plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
#plt.clim([9500, 11000])
# plt.xlim([155, 355])
# plt.ylim([155, 355])
plt.show()
# plt.close(fig)

# Inputting coordinates and magnitude of an object
while (1):
    ell_circle_check = 1
    while (1):
        if ell_circle_check == 1:
            annul_check = int(input("An ellipse annulus or a circle one? [1 for ellipse, 0 for circle]: "))
            ell_circle_check = 0
        else:
            annul_check = int(
                input("An unregistered command. An ellipse annulus or a circle one? [1 for ellipse, 0 for circle]: "))
        if annul_check == 1 or annul_check == 0:
            break
    source_x = int(raw_input("Y coordinate of an object: ", ))
    source_y = int(raw_input("X coordinate of an object: ", ))
    source_color = float(raw_input("GAIA BP-RP color of an object (even an estimation): ", ))
    ap_radius = int(raw_input("Radius of an aperture (for source, reference, and checks): ", ))
    source_gap = int(raw_input("Length of a gap between aperture and background annulus (source): ", ))
    if annul_check == 0:
        source_annul = int(raw_input("Outer radius of a background annulus (source): ", ))
        test_frame = np.copy(check_frame)
        circle_ref(source_x, source_y, ap_radius, source_gap, source_annul,
                   test_frame)  # Drawing aperture and annulus on a frame
    else:
        source_annul1 = int(raw_input("1st outer radius of a background annulus (source): ", ))
        source_annul2 = int(raw_input("2nd outer radius of a background annulus (source): ", ))
        source_angle = int(raw_input("Angle of a background annulus ellipse(source): ", ))
        test_frame = np.copy(check_frame)
        ellipse_source(source_x, source_y, ap_radius, source_gap, source_annul1, source_annul2, source_angle,
                       test_frame)  # Drawing aperture and annulus on a frame
    plt.figure(2)
    im = plt.imshow(test_frame)
    plt.colorbar(im)
    plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
    plt.show()
    cont_check = 1
    while (1):
        if cont_check == 1:
            cont_check = int(raw_input("Do you want to change anything? [1 for changing, 0 to not change]: ", ))
        else:
            cont_check = raw_input(
                "Unregistered command, repeat. Do you want to change anything? [1 for changing, 0 to not change]: ", )
        if cont_check == 1 or cont_check == 0:
            break
    if cont_check == 0:
        plt.close(2)
        break
    else:
        plt.close(2)

final_frame = np.copy(test_frame)
# Inputting coordinates and magnitude of a reference
while (1):
    ref_x = int(raw_input("Y coordinate of a reference: ", ))
    ref_y = int(raw_input("X coordinate of a reference: ", ))
    mag_ref = float(raw_input("Magnitude of a reference star: ", ))
    ref_gap = int(raw_input("Length of a gap between aperture and background annulus (reference): ", ))
    ref_annul = int(raw_input("Outer radius of a background annulus (reference): ", ))
    test_frame = np.copy(final_frame)
    circle_ref(ref_x, ref_y, ap_radius, ref_gap, ref_annul, test_frame)  # Drawing aperture and annulus on a frame
    plt.figure(2)
    im = plt.imshow(test_frame)
    plt.colorbar(im)
    plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
    plt.show()
    while (1):
        cont_check = int(raw_input("Do you want to change anything? [1 for changing, 0 to not change]: ", ))
        if cont_check == 1 or cont_check == 0:
            break
        else:
            cont_check = raw_input(
                "Unregistered command, repeat. Do you want to change anything? [1 for changing, 0 to not change]: ", )
    if cont_check == 0:
        plt.close(2)
        break
    else:
        plt.close(2)
final_frame = np.copy(test_frame)

#Inputting coordinates of checks
check_x = []
check_y = []
check_true_mag = []
check_color = []
check_for_gaps = 1
secondary_frame = np.copy(final_frame)
check_for_adding = 0
while (1):
    test_x = input("Y coordinate of a check {}: ".format(np.shape(check_x)[0] + 1), )  # X coordinate for a check
    test_y = input("X coordinate of a check {}: ".format(np.shape(check_x)[0] + 1), )  # Y coordinate for a check
    test_mag = input(
        "GAIA G magnitude of a check {}: ".format(np.shape(check_x)[0] + 1), )  # GAIA G magnitude of a star
    test_color = input(
        "GAIA BP-RP color of a check {}: ".format(np.shape(check_x)[0] + 1), )  # GAIA BP-RP color of a star
    if check_for_gaps == 1:
        check_gap = int(
            raw_input("Length of a gap between aperture and background annulus (common for all checks): ", ))
        check_annul = int(raw_input("Outer radius of a background annulus (common for all checks): ", ))
    check_for_gaps += 1

    test_frame = np.copy(secondary_frame)
    circle_ref(test_x, test_y, ap_radius, check_gap, check_annul, test_frame)  # Drawing aperture and annulus on a frame
    plt.figure(2)
    im = plt.imshow(test_frame)
    plt.colorbar(im)
    plt.clim([c_min, c_min + math.sqrt(c_min + readout_noise**2) * 8])
    plt.show()
    while (1):
        cont_check = int(raw_input("Do you want to change anything? [1 for changing, 0 to not change]: ", ))
        if cont_check == 1 or cont_check == 0:
            break
        else:
            cont_check = raw_input(
                "Unregistered command, repeat. Do you want to change anything? [1 for changing, 0 to not change]: ", )
    if cont_check == 0:
        plt.close(2)
        check_x = np.append(check_x, test_x)
        check_y = np.append(check_y, test_y)
        check_true_mag = np.append(check_true_mag, test_mag)
        check_color = np.append(check_color, test_color)
        secondary_frame = np.copy(test_frame)
        while (1):
            check_for_adding = int(
                raw_input("Do you want to stop adding checks? 1 to stop adding, 0 to keep adding: ", ))
            if check_for_adding == 1 or check_for_adding == 0:
                break
            else:
                check_for_adding = raw_input(
                    "Unregistered command, repeat. Do you want to stop adding checks? 1 to stop adding, 0 to keep adding: ", )
    else:
        plt.close(2)
        check_for_gaps -= 1
    if check_for_adding == 1:
        break
print("Source, reference and checks are chosen")

# Photometry
photo_len = len(photo_frame_path)
ref_photo = np.zeros(photo_len)
frame_datetime = np.zeros_like(ref_photo)
frame_exp = np.zeros_like(ref_photo)
frame_tube = np.array([])
frame_filter = np.array([])

check_len = len(check_x)
check_photo = np.zeros((photo_len, check_len))
check_photo_err = np.zeros((photo_len, check_len))
source_photo = np.zeros(photo_len)
source_err = np.zeros(photo_len)
# Conducting the photometry
for frame_num in range(len(photo_frame_path)):
    test_frame = fits.open("{}".format(filepath + photo_frame_path[frame_num]))
    test_frame_header = test_frame[0].header
    frame_exp[frame_num] = test_frame_header['EXPTIME']
    frame_datetime[frame_num] = test_frame_header['JD'] + frame_exp[frame_num] / 2 / 60 / 60 / 24
    frame_tube = np.append(frame_tube, test_frame_header['AXIS'])
    frame_filter = np.append(frame_filter, test_frame_header['FILTER'])
    photo_frame = fits.getdata("{}".format(filepath + photo_frame_path[frame_num]))
    ref_photo[frame_num] = photometry_ref(ref_x, ref_y, ap_radius, ref_gap, ref_annul, photo_frame)
    if annul_check == 0:
        source_photo[frame_num], source_err[frame_num] = photometry_source(source_x, source_y, ap_radius, source_gap,
                                                                           source_annul, photo_frame)
    else:
        source_photo[frame_num], source_err[frame_num] = photometry_source_ell(source_x, source_y, ap_radius,
                                                                               source_gap, source_annul1, source_annul2,
                                                                               source_angle, photo_frame)
    check_photo[frame_num], check_photo_err[frame_num] = photometry_check(check_x, check_y, ap_radius, check_gap,
                                                                          check_annul, photo_frame)

print(source_photo)
print(ref_photo)
# Determining magnitude
source_mag, check_mag = magnitude(mag_ref, ref_photo, source_photo, check_photo)
# Determining magnitude errors for source and checks
source_mag_var_err, source_mag_signal_err, check_mag_var_err, check_mag_signal_err = magnitude_err(source_photo,
                                                                                                   source_err,
                                                                                                   check_mag,
                                                                                                   check_photo,
                                                                                                   check_photo_err)

# Finding parameter for correcting for frame inadequacies and atmospheric effects(Forbes effect)
x_par, x_par_err = full_var_correct(check_mag, check_mag_signal_err,
                                                       check_true_mag, check_color, check_x, check_y, check_shape)

#Correcting checks and source
source_mag_corr = source_mag_correct(source_mag, source_color, source_x, source_y, x_par, check_shape)
check_mag_corr = check_mag_correct(check_mag, check_color, check_x, check_y, x_par, check_shape)

source_mag_corr_var_err, source_mag_signal_err, check_mag_corr_var_err, check_mag_signal_err = magnitude_err(
    source_photo, source_err, check_mag_corr, check_photo, check_photo_err)

# Estimating source mag error from check signal noise
source_mag_check_signal_err = np.zeros_like(source_mag)
signal_err_const = 0
for k in range(photo_len):
    for i in range(check_len):
        signal_err_const += check_mag_signal_err[k][i] ** 2
    source_mag_check_signal_err[k] = math.sqrt(signal_err_const / (check_len - 1))


#Determining deviation from true mag for every check
dev_data = {"Frame": photo_frame_path}
deviation_dataframe = pd.DataFrame(data=dev_data, index=frame_datetime)
deviation_dataframe.index.name = 'Julian Days, middle of exposure'


for k in range(check_mag.shape[1]):
    check_mag_deviation_data = np.zeros_like(source_mag)
    check_mag_signal_err_data = np.zeros_like(source_mag)
    check_mag_color = np.zeros_like(source_mag)
    check_mag_true_mag = np.zeros_like(source_mag)
    for j in range(check_mag.shape[0]):
        check_mag_deviation_data[j] = check_mag[j][k] - check_true_mag[k]
        check_mag_signal_err_data[j] = check_mag_signal_err[j][k]
        check_mag_color[j] = check_color[k]
        check_mag_true_mag[j] = check_true_mag[k]

    df = {'Check star ' + "{}".format(k + 1) + ' deviation, mag': check_mag_deviation_data,
          'Check star ' + "{}".format(k + 1) + 'mag error': check_mag_signal_err_data,
          'Check star ' + "{}".format(k + 1) + ' GAIA G mag': check_mag_true_mag,
          'Check star ' + "{}".format(k + 1) + ' BP-RP color': check_mag_color}

    df2 = pd.DataFrame(data=df, index=frame_datetime)
    deviation_dataframe = deviation_dataframe.join(df2)

#Determining deviation from true mag for every check
check_dev_first_frame = np.zeros_like(check_true_mag)
check_signal_err_first_frame = np.zeros_like(check_true_mag)

for k in range(check_mag.shape[1]):
        j=0
        check_dev_first_frame[k] = check_mag[j][k] - check_true_mag[k]
        check_signal_err_first_frame[k] = check_mag_signal_err[j][k]

check_color_first_frame = np.copy(check_color)
check_true_mag_first_frame = np.copy(check_true_mag)
df1 = {'Check stars BP-RP color': check_color_first_frame,
       'Check stars GAIA G mag': check_true_mag,
       'Check stars deviation on first frame, mag': check_dev_first_frame,
       'Check star signal error, mag': check_signal_err_first_frame,
       'Check star x coord - x mid': check_x - check_mid_x,
       'Check star y coord - y mid': check_y - check_mid_y,
       'Check star x coord': check_x,
       'Check star y coord': check_y}
first_frame_dev_dataframe = pd.DataFrame(data=df1)


# Fit data in a

fit = {"Frame": photo_frame_path}
fit_data = pd.DataFrame(data=fit, index=frame_datetime)
coord_par_num = (x_par.shape[1]-1)/2
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
        export_column = "Y coord parameter " + "{}".format(m-coord_par_num)
        export_err_column = "Y coord parameter " + "{}".format(m-coord_par_num) + "error"
    df = {export_column: export_x_par,
          export_err_column: export_x_par_err}
    df_export = pd.DataFrame(data=df, index=frame_datetime)
    fit_data = fit_data.join(df_export)

#Determining deviation of theoretical prediction from the data and comparing with empirical (squared sum of errors)
diff_empiric_dev = np.zeros(check_mag.shape[0])         #Difference between deviation of data and theoretical deviation
diff_check_err = np.zeros(check_mag.shape[0])
diff_sigma = np.zeros(check_mag.shape[0])
for k in range(check_mag.shape[0]):
    for m in range(check_mag.shape[1]):
        diff_check_dev = check_mag[k][m]-check_true_mag[m]
        diff_check_theor_dev = x_par[k][0]
        power_check = 1
        for i in range(1, color_param_num + 1):
            diff_check_theor_dev += check_color[m] ** power_check * x_par[k][i]
            power_check += 1
        power_check = 1
        for i in range(color_param_num + 1, x_par.shape[1], 2):
            diff_check_theor_dev += x_par[k][i] * (check_x[m] - check_mid_x) ** power_check
            diff_check_theor_dev += x_par[k][i+1] * (check_y[m] - check_mid_y) ** power_check
            power_check += 1
        diff_check_err[k] += check_mag_signal_err[k][m]**2
        diff_empiric_dev[k] += (diff_check_dev - diff_check_theor_dev)**2
    diff_check_err[k] = math.sqrt(diff_check_err[k])
    diff_empiric_dev[k] = math.sqrt((diff_empiric_dev[k]))
    diff_sigma[k] = diff_empiric_dev[k] / diff_check_err[k]
diff_predf = {"Root of a squared sum of check signal errors": diff_check_err,
              "Root of a squared sum of differences between deviation of check mag from the GAIA mag and theoretical deviation": diff_empiric_dev,
              "Ratio of a squared sum of check signal errors and differences": diff_sigma}
diff_df = pd.DataFrame(data = diff_predf, index = frame_datetime)

fit_data = fit_data.join(diff_df)
fit_data.index.name = 'Julian Days, middle of exposure'

# Exporting results of the photometry into excel in the same folder
data = {"Frame": photo_frame_path, 'Tube': frame_tube, 'Filter': frame_filter,
        'Frame exposure time, seconds': frame_exp, "Source mag (corrected)": source_mag_corr,
        'Source mag error (variability, corrected)': source_mag_corr_var_err,
        "Source mag error (signal noise)": source_mag_signal_err,
        "Source mag error (signal noise of checks)": source_mag_check_signal_err}

full_data = {"Frame": photo_frame_path, 'Tube': frame_tube, 'Filter': frame_filter,
             'Frame exposure time, seconds': frame_exp, 'Source mag': source_mag,
             "Source mag (corrected)": source_mag_corr,
             'Source mag error (variability, not corrected)': source_mag_var_err,
             'Source mag error (variability, corrected)': source_mag_corr_var_err,
             "Source mag error (signal noise)": source_mag_signal_err,
             "Source mag error (signal noise of checks)": source_mag_check_signal_err}
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

test_frame = fits.open("{}".format(filepath + filename))
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
