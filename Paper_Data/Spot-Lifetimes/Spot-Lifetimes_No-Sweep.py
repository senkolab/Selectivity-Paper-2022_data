# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 17:40:06 2021

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

FocusWaist = 100*1e-6

spot_alpha = 0.5
markersize = 3
markeredgewidth = 0.5
linewidth = 1
elinewidth = 1
legend_fontsize = 5
axis_fontsize = 10
ticks_fontsize = 6
graph_edge_width = 0.5
graph_tick_width = graph_edge_width

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
golden_ratio = (5**.5 - 1) / 2
fig_width_in = columnwidth * inches_per_pt
height_ratio = golden_ratio
fig_height_in = fig_width_in * height_ratio

def getData(data_file, NumSamples = 1, long=False):
    raw_data = np.asarray(pd.read_csv(data_file, delimiter = "[\],\[\t]", engine='python',header=None))
    if NumSamples == 1:
        # print(raw_data)
        counts = raw_data[:,2]
    else:
        if not long:
            raw_counts = raw_data[:,2:2+NumSamples]
            TotalPulses = raw_counts.size
            counts = np.mean(raw_counts, 1)
        else:
            raw_counts = raw_data[:,2:2+NumSamples]
            counts = raw_counts.flatten()
    if not long:
        times = raw_data[:,NumSamples + 17]
        times -= min(times)
        times = FixTimeOffset(times)
        datas = np.vstack((times, counts)).T
    else:
        times = raw_data[:,NumSamples + 17]
        times -= min(times)
        times = np.append(times, times[-1] + 25)
        times = np.array([np.linspace(times[i-1], time, NumSamples) for i, time in enumerate(times) if time != 0])
        times = times.flatten()
        datas = np.vstack((times, counts)).T
    return datas

def FixTimeOffset(data):#Fix timing when somehow it's less than 0.1 s between successive data points
    data = data.astype(float)
    dataDiff = np.diff(data)
    dataDiff2 = np.concatenate((dataDiff, np.array([0])))
    TimeDiff = np.mean(dataDiff)
    # print(f'Time Difference is {TimeDiff}')
    for i, x in enumerate(data):
        if i > 0:
            if data[i] - data[i-1] < 1e-1:
                # print(f"time: {data[i-1]}, time2: {data[i]}")
                data[i:] += TimeDiff
    return data
    
def getSweepTimings(data_file, periods=3):
    times = []
    with open(data_file) as file:
        for line in file:
            linesplit = line.split(' ')
            time = linesplit[1]
            split_time = re.split(':|\.', time)
            time_hour, time_minute, time_second, time_msecond = [float(time) for time in split_time]
            total_seconds = 60*60*time_hour + 60*time_minute + time_second + 1e-3*time_msecond
            times.append(total_seconds)
    times = np.array(times)
    times -= min(times)
    times_diff = np.insert(np.diff(times), 0, 0)
    times = times[times_diff > 0.2]
    steps_per = times.size
    times_expanded = np.array([times+i*max(times) for i in range(periods)]).flatten()
    times_expanded -= min(times_expanded)
    return times_expanded, steps_per

def Down_Sample(data, Num):
    Rem = data.size%Num
    if Rem != 0:
        data = data[:-Rem]
    else:
        data = data
    data = data.reshape(-1, Num).astype(float)
    data_DownSampled = np.mean(data, 1)
    data_DownSampled_std = np.std(data, 1)
    return data_DownSampled, data_DownSampled_std
    
data_filepath = r"2020_11_04_Pitting-Conditioning\Condition-Sacrif1\Neutral_Conditioning_Sacrif1*"
data_files = glob.glob(data_filepath)
data_file_BG = [s for s in data_files if "No-Ablation" in s][0]
data_files_BeforeCond = [s for s in data_files if "Before" in s]
data_files_BeforeCond.pop(0)#Only need 1 example for figure
data_files_AfterCond = [s for s in data_files if "After" in s]
data_file_LongTerm = data_files_AfterCond[-1]
data_files_AfterCond.pop(-1)#Last file was long-term
Condition_Fluences = np.array([130, 190])*1e-6/(math.pi*FocusWaist**2)#In J/m^2


data_BG = getData(data_file_BG)
BG_Ave = np.mean(data_BG[:,1])


fig = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300)
ax2 = fig.add_subplot(111)

PulsesPer = 120
DownSample = 3
data = getData(data_file_LongTerm, PulsesPer, long=True)
data_y = data[:,1]
data_y = data_y[:-1]
data_yD, data_yDsdv = Down_Sample(data_y, 3)
data_x = np.arange(data_yD.size)*DownSample


# ax2.plot(data_x, data_yD, label="Sweep", alpha=0.3, color='grey')


data_file = r"Neutral-Fluorescence_Natural_Pitting_48uJ_7p26x6p13_After-Condition-1_*"
data_file = glob.glob(data_file)[0]
Pit1Num = 5
Pit1 = getData(data_file, Pit1Num)[:,1]

Pit2Num = 10
data_file = r"Neutral-Fluorescence_Natural_Pitting_48uJ_7p265x6p13_After-Condition-*"
data_file = glob.glob(data_file)[0]
Pit2 = getData(data_file, Pit2Num)[:,1]

data_file = r"Neutral_Pitting_7p28x6p125_After-Condition-1_*"
data_file = glob.glob(data_file)[0]
Pit3 = getData(data_file, Pit2Num)[:,1]

data_file = r"Neutral_Pitting_7p265x6p15_After-Condition-2_*"
data_file = glob.glob(data_file)[0]
Pit4 = getData(data_file, Pit2Num)[:,1]

Down = 20
Pit1A, Pit1SD = Down_Sample(Pit1, Down)
Pit2A, Pit2SD = Down_Sample(Pit2, Down)
Pit3A, Pit3SD = Down_Sample(Pit3, Down)
Pit4A, Pit4SD = Down_Sample(Pit4, Down)

Scale1 = 0.7
Scale2 = 0.4
Scale3 = 0.8
Scale4 = 0.6
Scale1 = Scale2 = Scal3 = Scale4 = 1

ax2.errorbar(np.arange(0, Pit1A.size*Down*Pit1Num, Down*Pit1Num), Pit1A*Scale1, yerr=Pit1SD, fmt='^', \
             label='Spot 1', alpha=spot_alpha, lw=linewidth, elinewidth=elinewidth, \
                 ms=markersize, mew=markeredgewidth)
ax2.errorbar(np.arange(0, Pit2A.size*Down*Pit2Num, Down*Pit2Num), Pit2A*Scale2, yerr=Pit2SD, fmt='s', \
             label='Spot 2', alpha=spot_alpha, lw=linewidth, elinewidth=elinewidth, \
                 ms=markersize, mew=markeredgewidth)
ax2.errorbar(np.arange(0, Pit3A.size*Down*Pit2Num, Down*Pit2Num), Pit3A*Scale3, yerr=Pit3SD, fmt='d', \
             label='Spot 3', alpha=spot_alpha, lw=linewidth, elinewidth=elinewidth, \
                 ms=markersize, mew=markeredgewidth)
ax2.errorbar(np.arange(0, Pit4A.size*Down*Pit2Num, Down*Pit2Num), Pit4A*Scale4, yerr=Pit4SD, fmt='o', \
             label='Spot 4', alpha=spot_alpha, lw=linewidth, elinewidth=elinewidth, \
                 ms=markersize, mew=markeredgewidth, color='tab:purple')
average_lifetime = np.mean([Pit1A.size*Down*Pit1Num, Pit2A.size*Down*Pit2Num, Pit3A.size*Down*Pit2Num, Pit4A.size*Down*Pit2Num])
ax2.vlines(average_lifetime, 0, 200, label='Mean lifetime', color='k', linestyle='dotted', lw=linewidth)

## Inset target map from inkscape
im = plt.imread('target_inset.png') # insert local path of the image, old version: 
newax = fig.add_axes([0.45,0.6,0.2,0.2], anchor='W') 
newax.imshow(im)
newax.axis('off')

ax2.set_xlabel('Pulse number',fontsize=axis_fontsize)
ax2.set_ylabel('Counts (arb. units)',fontsize=axis_fontsize)
ax2.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
ax2.xaxis.set_major_locator(MultipleLocator(5000))
ax2.xaxis.set_minor_locator(MultipleLocator(1000))
ax2.yaxis.set_major_locator(MultipleLocator(50))
ax2.yaxis.set_minor_locator(MultipleLocator(10))
ax2.set_xlim(0, 13000)
ax2.set_ylim(0, 150)

# ax2.text(5200, -47, '(b)', ha='center', fontsize=8, font='Times New Roman')

ax2.legend(fontsize=legend_fontsize, framealpha=1)



fig.savefig('Spot_Lifetimes_Final_v3.pdf', dpi=300, bbox_inches='tight', format='pdf')