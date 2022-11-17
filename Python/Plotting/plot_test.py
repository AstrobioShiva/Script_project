import ROOT
from iminuit import Minuit

# we also need a cost function to fit and import the LeastSquares function
from iminuit.cost import LeastSquares

# display iminuit version
import iminuit

import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from numpy import inf
import pandas as pd
import csv
from scipy.optimize import curve_fit
df= pd.read_csv('freq_z_rot.csv', header = 0)
# headers = ['angle', 'frequency']
# df= pd.read_csv('sum_frequncy_z_rot.csv', names= headers)
# df = pd.read_csv('sum_frequncy_z_rot.csv', header=0)
# df_d = pd.read_csv('diff_freq_z_rot.csv', header=0)
# input from file
x = df['angle']
y_s = df['fre_sum']
y_d = df['freq_diff']
y_1D = df['freq_1D'] #freq of |1> -> |0>
y_2D = df['freq_2D'] #freq of |0> -> |-1>
#y = y/10E3
# plt.plot(x,y, '--')
# plt.title(r'$Plot\ of\ |1>\leftrightarrow|0> + |0>\leftrightarrow|-1>$')
# plt.xlabel("Angle")
# plt.ylabel("frequency (Hz)")
# plt.tight_layout()
# plt.savefig("sum_freq_z_rot")

plt.plot(x,y_d, '--')
plt.title(r'$Plot\ of\ |1>\leftrightarrow|0> - |0>\leftrightarrow|-1>$')
plt.xlabel("Angle")
plt.ylabel("frequency (Hz)")
plt.tight_layout()
plt.savefig("diff_freq_z_rot")

# plt.show()

def Hcsa1(theta,A1, B1, C1):
    
    freq_csa1 = A1 + B1*np.cos(2*theta*np.pi/180) + C1*np.sin(2*theta*np.pi/180)
    return freq_csa1

def HQ2 (theta,A2, B2, C2, D2, E2):
    freq_Q2 = A2 + B2*np.cos(2*theta*np.pi/180) + C2*np.sin(2*theta*np.pi/180) + D2*np.cos(4*theta*np.pi/180) + E2*np.sin(4*theta*np.pi/180) 
    return freq_Q2

#Curve Fitting using Iminiuit

# least_squares = LeastSquares(x, y,y/10, Hcsa1)
# #least_squares = LeastSquares(x, y,y/10, HQ2)
# m = Minuit(least_squares, A1 = 1, B1 = 1, C1= 1)
# #m = Minuit(least_squares, A2 = 1, B2 = 1, C2= 1, D2 = 1, E2 = 1)
# # m.fixed["C1"] = True

# m.migrad()  # finds minimum of least_squares function
# m.hesse()   # accurately computes uncertainties
# fit_info = [
#     f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.5f} / {len(x) - m.nfit}",
# ]
# print(fit_info[0])


# for p, v, e in zip(m.parameters, m.values, m.errors):
#     fit_info.append(f"{p} = ${v:.3f} \\pm {e:.3f}$")


# # draw data and fitted line
# plt.errorbar(x, y, y/10., fmt="o", label="data")
# plt.plot(x, Hcsa1(x, *m.values), label="fit")
# #plt.plot(x, HQ2(x, *m.values), label="fit")
# plt.legend(title="\n".join(fit_info))
# plt.savefig("miniuit_fit_Hcsa1")

# print(Hcsa1(0, *m.values))
# #print(HQ2(0, *m.values))
# print(y[0])