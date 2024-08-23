import pandas as pd
import numpy as np
import cmath
from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename
# Choosing a file
print("Choose the file to read")
data_path = askopenfilename()

data_df = pd.read_csv(data_path, sep='\s+')
#print(data_df)
#print(data_df['rlen'])
max_likelihood_distance = np.zeros(0)
# Calculating the max likelihood distance
for x_ind in range(len(data_df['rlen'])):
    parallax_val = data_df['Plx'][x_ind] + 0.029
    parallax_val_err = data_df['e_Plx'][x_ind]
    rlen_val = data_df['rlen'][x_ind]
    # Cordano
    a_val = - 2 * rlen_val
    b_val = parallax_val * rlen_val / (parallax_val_err) ** 2

    rho_val = - a_val ** 2 / 3 + b_val
    q_val = 2 * (a_val / 3) ** 3 - 1 / 3 * a_val * b_val

    q_big_val = cmath.sqrt((rho_val / 3) ** 3 + (q_val / 2) ** 2)

    alpha_val = cmath.exp( cmath.log(- q_val / 2 + q_big_val) * 1 / 3)
    beta_val = cmath.exp( cmath.log(- q_val / 2 - q_big_val) * 1 / 3)
    dist_max_poss = alpha_val + beta_val - a_val / 3

    max_likelihood_distance = np.append(max_likelihood_distance, dist_max_poss.real)
data_df['r_rest'] = max_likelihood_distance
print("Choose where to save data")
export_path = askdirectory()
data_df.to_csv(export_path + "/cepheid_distance.csv")