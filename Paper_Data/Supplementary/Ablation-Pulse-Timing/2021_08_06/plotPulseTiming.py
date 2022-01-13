# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 14:43:10 2021

@author: brend
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#print(plt.style.available)
plt.style.use('seaborn-colorblind')
# plt.rcParams['figure.facecolor'] = '1'
from scipy.optimize import curve_fit

import glob
import pandas as pd
import re

def getData(data_file, NumSamples = 1):
    raw_data = np.asarray(pd.read_csv(data_file, delimiter = "[\],\[\t]", engine='python',header=None))
    data_counts = raw_data[:, 2:12]
    datas = data_counts
    return datas

def fitLorentz(time, A, peak_time, FWHM):
    return (0.5*FWHM)/(math.pi*((time - peak_time)**2 + (0.5*FWHM)**2))

def fitGaussian(time, A, peak_time, std):
    return A*np.exp(-(time-peak_time)**2/(2*std**2))

data_filepath = r"Ablation-TOF_JustLaser_002*"
data_files = glob.glob(data_filepath)
datas = getData(data_files[0])
data_counts_ave = np.mean(datas, 1)

time_start = 139
time_stop = 146
time_step = 0.5
times_data = np.arange(time_start, time_stop+time_step/2, time_step)

A_init = 1
A_lower = 0.1
A_upper = 10
peak_time_init = 142
peak_time_lower = 141
peak_time_upper = 144
FWHM_init = 1
FWHM_lower = 0.1
FWHM_upper = 3
params_init = np.array([A_init, peak_time_init, FWHM_init])
lower_bounds = np.array([A_lower, peak_time_lower, FWHM_lower])
upper_bounds = np.array([A_upper, peak_time_upper, FWHM_upper])
lu_bounds = np.array([lower_bounds, upper_bounds])
fit_params_1, pcov = curve_fit(fitLorentz, times_data, data_counts_ave, p0=params_init, bounds=lu_bounds)
times_fit_1 = np.arange(time_start, time_stop, 0.01)
counts_fit_1 = fitLorentz(times_fit_1, fit_params_1[0], fit_params_1[1], fit_params_1[2])

A_init = 1
A_lower = 0.1
A_upper = 10
peak_time_init = 142
peak_time_lower = 141
peak_time_upper = 144
std_init = 1
std_lower = 0.1
std_upper = 3
params_init = np.array([A_init, peak_time_init, std_init])
lower_bounds = np.array([A_lower, peak_time_lower, std_lower])
upper_bounds = np.array([A_upper, peak_time_upper, std_upper])
lu_bounds = np.array([lower_bounds, upper_bounds])
fit_params_2, pcov = curve_fit(fitGaussian, times_data, data_counts_ave, p0=params_init, bounds=lu_bounds)
times_fit_2 = np.arange(time_start, time_stop, 0.01)
counts_fit_2 = fitGaussian(times_fit_2, fit_params_2[0], fit_params_2[1], fit_params_2[2])

fig1 = plt.figure(figsize=(20, 10))
ax1 = fig1.add_subplot(111)

data_counts_plot = data_counts_ave
times = times_data
ax1.plot(times, data_counts_plot, label='data')
ax1.set_xlabel(r'Time / $\mu s$',fontsize=30, labelpad=20)
ax1.set_ylabel('Counts / a.u.',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)

data_counts_plot = counts_fit_1
times = times_fit_1
ax1.plot(times, data_counts_plot, label=f'Lorentzian fit: \nA={fit_params_1[0]:0.1f}, t_0={fit_params_1[1]:0.2f}, FWHM={fit_params_1[2]:0.2f}')
ax1.set_xlabel(r'Time / $\mu s$',fontsize=30, labelpad=20)
ax1.set_ylabel('Counts / a.u.',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)

data_counts_plot = counts_fit_2
times = times_fit_2
ax1.plot(times, data_counts_plot, label=f'Gaussian fit: \nA={fit_params_2[0]:0.1f}, t_0={fit_params_2[1]:0.2f}, $\sigma$={fit_params_2[2]:0.2f}')
ax1.set_xlabel(r'Time / $\mu s$',fontsize=30, labelpad=20)
ax1.set_ylabel('Counts / a.u.',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)

ax1.legend(fontsize=20)