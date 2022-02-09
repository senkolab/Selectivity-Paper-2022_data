# -*- coding: utf-8 -*-
"""
Created on Thu May  6 09:11:15 2021

@author: brend
"""

import math
import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
import pandas as pd
plt.style.use('seaborn-colorblind')

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
golden_ratio = (5**.5 - 1) / 2
fig_width_in = fullwidth * inches_per_pt
height_ratio = golden_ratio
fig_height_in = fig_width_in * height_ratio

def getDataCounts(data_path):
    data_file = np.asarray(glob.glob(data_path))[0]
    data = np.asarray(pd.read_csv(data_file, delimiter = "[\]\[\t]", engine='python', skiprows=[0], header=None))
    data_header = []
    for j, data_row in enumerate(data):
        if data_row[0][0] == 'h':
            data_header = data_row
            break
    data = np.delete(data, j, axis=0)
    if data_header.size == 0:
        print("No data_header - assuming freq. in column 1, counts in column 2")
        counts_row = 2
    else:
        counts_row = np.where(data_header == 'IonCounts138')[0]
    counts_data = data[:,counts_row].astype(float).flatten()
    return counts_data

def getAverageCountsSingle(counts, trapping):
    summ = 0
    total_trapped = 0
    for i, count in enumerate(counts):
        if count != 0:
            summ += count
            total_trapped += trapping[i]
    return summ/total_trapped
    

data_path = r'2021_05_05\Load-Effic-Counts_Ba138_140uJ_110uW_1.*'

counts_data = getDataCounts(data_path)

BG = 600

counts_data -= BG

freq1 = -120
freq2 = -80
freq3 = -40
freq4 = 0

optimal_freq = 607.42604

freq1_counts = counts_data[0::4]
freq2_counts = counts_data[1::4]
freq3_counts = counts_data[2::4]
freq4_counts = counts_data[3::4]


freq4_trapping = [0, 1, 1, 3, 1, 2, 0, 1, 1, 0]
freq3_trapping = [3, 1, 1, 2, 2, 3, 0, 3, 2, 1]
freq2_trapping = [1, 3, 1, 1, 2, 0, 3, 2, 3, 1]
freq1_trapping = [0, 3, 2, 1, 1, 2, 3, 2, 2, 1]

freq1_singlecounts = getAverageCountsSingle(freq1_counts, freq1_trapping)
freq2_singlecounts = getAverageCountsSingle(freq2_counts, freq2_trapping)
freq3_singlecounts = getAverageCountsSingle(freq3_counts, freq3_trapping)
freq4_singlecounts = getAverageCountsSingle(freq4_counts, freq4_trapping)

data_ploteff = np.array([np.mean(freq1_trapping), np.mean(freq2_trapping), np.mean(freq3_trapping), np.mean(freq4_trapping)])
data_err_ploteff = [np.std(freq1_trapping)/math.sqrt(10), np.std(freq2_trapping)/math.sqrt(10), np.std(freq3_trapping)/math.sqrt(10), np.std(freq4_trapping)/math.sqrt(10)]
x_plot = [freq1, freq2, freq3, freq4]
data_plotcounts = np.array([freq1_singlecounts, freq2_singlecounts, freq3_singlecounts, freq4_singlecounts])*10

fig = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

ax1.errorbar(x_plot, data_ploteff, yerr=data_err_ploteff, label='Loading rate', fmt='o', markersize=3, mew=0.5, elinewidth=1)
ax2.plot(x_plot, data_plotcounts, color='r', label='Counts', marker='o', linestyle='None', markersize=3)

ax1.set_xlabel('Frequency (MHz + %8.6f THz)'%optimal_freq,fontsize=10)
ax1.set_ylabel('Loading rate (ions/pulse)',fontsize=10)
ax1.tick_params(axis='both', which='major', labelsize=7)
ax1.legend(fontsize=7, loc=2)
ax2.set_ylim(0, 24000)
ax2.tick_params(axis='both', which='major', labelsize=6)
ax2.set_ylabel('Ion brightness (counts/s)', fontsize=10)
ax2.legend(fontsize=7, loc=1)
plt.show()

fig.savefig('LoadingRate-vs-Cooling-Freq_v2.pdf', dpi=300, bbox_inches='tight', format='pdf')