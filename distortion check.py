import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits import mplot3d
from frame_functions import *
from coord_functions import *
rad_conversion = 57.295828
distortion_parameters = np.array([-5.6 * 10 ** (-8), 2.7 * 10 ** (-14), 9. * 10 ** (-6), 4. * 10 ** (-6)], dtype=np.longdouble)
#distortion_parameters = np.array([-5.6 * 10 ** (-8), 0, 0, 0], dtype=np.longdouble)
distorted_x = np.array([300., 200., 100., 0., -100., -200., -300.]) * 2
distorted_y = np.array([300., 200., 100., 0., -100., -200., -300.]) * 2


undistorted_x = np.zeros_like(distorted_x)
undistorted_y = np.zeros_like(distorted_x)
centre_dev_x = 40
centre_dev_y = 40
radius = np.sqrt((distorted_x - centre_dev_x) ** 2 + (distorted_y - centre_dev_y) ** 2)
#print(distorted_x)
a_k = radius ** 2 + 2 * (distorted_x - centre_dev_x) ** 2
b_k = 2 * (distorted_x - centre_dev_x) * (distorted_y - centre_dev_y)
c_k = radius ** 2 + 2 * (distorted_y - centre_dev_y) ** 2
#print(a_k, b_k, c_k)
#print(distorted_x, distorted_y)
for i in range(distorted_x.shape[0]):
    undistorted_x[i] = distorted_x[i] + (distorted_x[i] - centre_dev_x) * (distortion_parameters[0] * radius[i] ** 2 + distortion_parameters[1] * radius[i] ** 4)
    undistorted_x[i] += distortion_parameters[2] * a_k[i] + distortion_parameters[3] * b_k[i]
    undistorted_y[i] = distorted_y[i] + (distorted_y[i] - centre_dev_y) * (
                distortion_parameters[0] * radius[i] ** 2 + distortion_parameters[1] * radius[i] ** 4)
    undistorted_y[i] += distortion_parameters[2] * b_k[i] + distortion_parameters[3] * c_k[i]
input_parameters = zeros(12)
#print(distorted_x, distorted_y)
#print(axis_fit_distortion(distorted_x, distorted_y, undistorted_x, undistorted_y, input_parameters))
#distortion_angle = 0.01 / 57.3
#dscale = 0.001
distortion_angle = 1. / rad_conversion
dscale = 0.000
scale = 1. + dscale
d_x = 50.
d_y = 0
inter_undistorted_x = copy(undistorted_x)
inter_undistorted_y = copy(undistorted_y)
#print(inter_undistorted_x, inter_undistorted_y)
undistorted_x = scale * inter_undistorted_x * math.cos(distortion_angle) + scale * inter_undistorted_y * math.sin(distortion_angle) + d_x * math.cos(distortion_angle) + d_y * math.sin(distortion_angle)
undistorted_y = - scale * inter_undistorted_x * math.sin(distortion_angle) + scale * inter_undistorted_y * math.cos(distortion_angle) - d_x * math.sin(distortion_angle) + d_y * math.cos(distortion_angle)
print('aaa')
#print(distorted_x, distorted_y)
#undistorted_x = scale * inter_undistorted_x * math.cos(distortion_angle) + scale * inter_undistorted_y * math.sin(distortion_angle)
#undistorted_y = - scale * inter_undistorted_x * math.sin(distortion_angle) + scale * inter_undistorted_y * math.cos(distortion_angle)
input_parameters[0] = centre_dev_x
input_parameters[1] = centre_dev_y
for i in range(len(distortion_parameters)):
    input_parameters[5 + i] = distortion_parameters[i]
#input_parameters[5] = distortion_parameters[0]
input_parameters[9] = -d_x
input_parameters[10] = -d_y
input_parameters[11] = -distortion_angle
#print(input_parameters)
#print(distorted_x, distorted_y)
print(distortion_min_sum(distorted_x, distorted_y, undistorted_x, undistorted_y, input_parameters))
#print(undistorted_x, undistorted_y)
supposed_x_0 = copy(undistorted_x)
supposed_y_0 = copy(undistorted_y)

#print(math.cos(distortion_angle), math.cos(input_parameters[11]))
supposed_x = supposed_x_0 * math.cos(input_parameters[11]) + supposed_y_0 * math.sin(input_parameters[11]) + input_parameters[9]
supposed_y = -supposed_x_0 * math.sin(input_parameters[11]) + supposed_y_0 * math.cos(input_parameters[11]) + input_parameters[10]
#print(supposed_x, supposed_y)
input_parameters = zeros(12)
input_parameters[0] = centre_dev_x
input_parameters[1] = centre_dev_y
#input_parameters[9] = -d_x
#input_parameters[10] = -d_y
#input_parameters[11] = - distortion_angle
print('Starting shifted center')
#print(input_parameters)
min_sum, inter_parameters = axis_fit_distortion(distorted_x, distorted_y, undistorted_x, undistorted_y, input_parameters)
print(min_sum, input_parameters)
print(frame_centre_conjugate_grad(distorted_x, distorted_y, undistorted_x, undistorted_y, inter_parameters))
#print('The end')
# Making a plot of min sums depending on the additive
x_array = np.arange(- d_x - 5, - d_x + 5, 0.1)
y_array = np.arange(- distortion_angle - 0.01, - distortion_angle + 0.01 - 0.0001, 0.0002)
#print(y_array)
z_array = np.zeros((len(x_array), len(y_array)))
inter_undistorted_x = copy(undistorted_x)
inter_undistorted_y = copy(undistorted_y)
for i in range(len(x_array)):
    for j in range(len(y_array)):
        #inter_undistorted_x = undistorted_x + x_array[i]
        #inter_undistorted_y = undistorted_y + y_array[j]
        inter_parameters = copy(input_parameters)
        inter_parameters[9] = x_array[i]
        inter_parameters[11] = y_array[j]
        z_array[i][j], inter_parameters = axis_fit_distortion(distorted_x, distorted_y, inter_undistorted_x, inter_undistorted_y, inter_parameters)
        #=
print(np.min(z_array))
X, Y = np.meshgrid(x_array, y_array)
plt.pcolormesh(X, Y, z_array)
#ax = plt.axes(projection = '3d')
#ax.plot_surface(X, Y, z_array)
plt.colorbar()
plt.clim([0, 1])
plt.show()

plt.close()