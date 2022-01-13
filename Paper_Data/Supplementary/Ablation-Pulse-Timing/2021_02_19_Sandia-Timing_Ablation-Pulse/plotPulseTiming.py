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

import glob
import pandas as pd
import re

def getData(data_file, NumSamples = 1):
    raw_data = np.asarray(pd.read_csv(data_file, delimiter = "[\],\[\t]", engine='python',header=None))
    data_counts = raw_data[:, 2:12]
    datas = data_counts
    return datas

data_filepath = r"Ablation-TOF_7p28x6p11-JustLaser_012*"
data_files = glob.glob(data_filepath)
datas = getData(data_files[0])
data_counts_ave = np.mean(datas, 1)

time_start = 139
time_stop = 155
time_step = 0.05
times_overall = np.arange(time_start, time_stop+time_step/2, time_step)
time_integration = 0.4
time_steps_integrate = int(time_integration/time_step)

y_n = np.zeros(data_counts_ave.shape)

#Attempting to tease out information about each 50 ns time-step, given 400 ns integration... Not working
for i in range(data_counts_ave.size-time_steps_integrate):
    j = data_counts_ave.size-time_steps_integrate-i
    print(j)
    y_n[j] = abs(data_counts_ave[j] - np.sum(data_counts_ave[j:j+time_steps_integrate-1]))/time_steps_integrate

fig1 = plt.figure(figsize=(20, 10))
ax1 = fig1.add_subplot(111)

data_counts_plot = y_n[:y_n.size-time_steps_integrate]
times = times_overall[:y_n.size-time_steps_integrate]
ax1.plot(times, data_counts_plot)
ax1.set_xlabel(r'Time / $\mu s$',fontsize=30, labelpad=20)
ax1.set_ylabel('Counts / a.u.',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.plot(times_overall, data_counts_ave)