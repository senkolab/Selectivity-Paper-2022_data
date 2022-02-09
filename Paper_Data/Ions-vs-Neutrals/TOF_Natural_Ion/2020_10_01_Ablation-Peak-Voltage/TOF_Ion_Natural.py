# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:41:19 2021

@author: brend
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use('seaborn-colorblind')
# plt.rcParams["font.family"] = "serif"
mpl.rcParams['axes.linewidth'] = 0.5

import glob
import pandas as pd
import re

plot_alpha = 0.5
markersize = 3
markeredgewidth = 0.5
linewidth = 1
elinewidth = 1
legend_fontsize = 8
axis_fontsize = 10
ticks_fontsize = 6
graph_edge_width = 0.5
graph_tick_width = graph_edge_width

plotcolor_ion = 'tab:purple'
plotcolor_bg = 'tab:blue'

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
golden_ratio = (5**.5 - 1) / 2
fig_width_in = fullwidth * inches_per_pt
height_ratio = golden_ratio
fig_height_in = fig_width_in * height_ratio

def get_data(file, samples=1):
    raw_data = np.asarray(pd.read_csv(file, delimiter = "\s\[|\]|\[|,\s", engine='python',header=None))
    counts = raw_data[:,2:2+samples]
    counts = counts.astype(np.float64)
    if samples == 1:
        counts_ave = counts
        counts_std = 0
    else:
        counts_ave = np.mean(counts, 1)
        counts_std = np.std(counts, 1)
    return counts, counts_ave, counts_std


samples = 120
time_ablation = 142.6e-6
time_neutral = 7.5e-6

time_step = 250e-9
time_start = 138e-6
time_stop = 155e-6
times = np.arange(time_start, time_stop, time_step)
times -= time_ablation

data_file_bg = r"TOF_nat_125uJ_no493_no650_001*"
data_file_bg = glob.glob(data_file_bg)[0]
data_raw_no650, data_ave_no650, data_std_no650 = get_data(data_file_bg, samples=samples)

data_file = r"TOF_nat_125uJ_with493_with650_001*"
data_file = glob.glob(data_file)[0]
data_raw, data_ave, data_std = get_data(data_file, samples=samples)



fig = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300)
ax1 = fig.add_subplot(111)

ax1.errorbar(times*1e6, data_ave, yerr=data_std/math.sqrt(samples), fmt='o', color=plotcolor_ion, \
             label='Ion fluorescence', alpha=plot_alpha, lw=linewidth, elinewidth=elinewidth, \
                 ms=markersize, mew=markeredgewidth)
ax1.errorbar(times*1e6, data_ave_no650, yerr=data_std_no650/np.sqrt(samples), fmt='^', color=plotcolor_bg, \
             label='Background (no repump)', alpha=plot_alpha, lw=linewidth, elinewidth=elinewidth, \
                 ms=markersize, mew=markeredgewidth)
line = ax1.vlines(time_neutral*1e6, np.min(data_ave)-5, np.max(data_ave)+5, alpha=0.5, \
                              color='k', label=r'Peak neutral atoms', linestyle='dotted', lw=linewidth)
    
ax1.set_xlabel(r'Time ($\mu s$)',fontsize=axis_fontsize)
ax1.set_ylabel('Fluorescence (counts)',fontsize=axis_fontsize)
ax1.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax1.set(xlim=(-2, 13), ylim=(0, 14))
ax1.legend(fontsize=legend_fontsize, framealpha=0.5, loc=1)
    
fig.savefig('Ion_TOF_v2.pdf', dpi=300, bbox_inches='tight', format='pdf')