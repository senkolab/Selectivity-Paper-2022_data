# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 13:51:04 2021

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

FocusWaist = 100*1e-6

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
def plotBackground(ax, Ave):
    [BG_x0, BG_x1] = ax.get_xlim()
    BG_x = np.array([BG_x0, BG_x1])
    BG_y = np.ones(BG_x.shape)*Ave
    ax.plot(BG_x, BG_y, label='Average background', color='k')
    
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

def generateHistogramBins(data_y, count_bins, BG_cutoff=0):
    hist_bins = []
    for i, bin_upper in enumerate(count_bins):
        if i == 0:
            pass
        else:
            bin_lower = count_bins[i-1]
            bin_size = data_y[np.logical_and(data_y > bin_lower,data_y < bin_upper, data_y > BG_cutoff)].size
            hist_bins.append(bin_size)
    return hist_bins
           

x_max = 120
    
    
data_filepath = r"Neutral_Conditioning_Sacrif1*"
data_files = glob.glob(data_filepath)
data_file_BG = [s for s in data_files if "No-Ablation" in s][0]
data_files_BeforeCond = [s for s in data_files if "Before" in s]
data_files_BeforeCond.pop(0)#Only need 1 example for figure
data_files_AfterCond = [s for s in data_files if "After" in s]
data_file_LongTerm = data_files_AfterCond[-1]
data_files_AfterCond.pop(-1)#Last file was long-term
Condition_Fluences = np.array([130, 190])*1e-6/(math.pi*FocusWaist**2)#In J/m^2

data_filepath = r"Ablation-Sweep_Nat-1-Sacrif1 2020-11-04 16-58-30*"
data_file = glob.glob(data_filepath)
data_file = data_file[0]
step_timings, steps_per = getSweepTimings(data_file, periods=5)
period_time = step_timings[steps_per]

bin_step = 5
bin_start = 0
bin_end = 250
count_bins = np.arange(bin_start, bin_end+1, bin_step)

data_BG = getData(data_file_BG)
BG_Ave = np.mean(data_BG[:,1])
BG_var = np.std(data_BG[:,1])/np.sqrt(data_BG.shape[0])


PulsesPer = 120
data_longsweep = getData(data_file_LongTerm, NumSamples = PulsesPer, long=True)

fig = plt.figure(figsize=(20, 10))
ax2 = fig.add_subplot(111)
ax2_in = ax2.inset_axes([35, 55, 160, 140], transform=ax2.transData)



log_histogram = False
transp_histogram = 1
reduction_transp_histogram = 0.3
data_bconds = []
for i, data_file in enumerate(data_files_BeforeCond):
    data_bcond = getData(data_file)
    data_x_bcond = data_bcond[:,0]
    data_y_bcond = data_bcond[:,1]
    label='Neutral fluorescence before conditioning'
    data_bconds.append(data_bcond)
    hist_bins = generateHistogramBins(data_y_bcond, count_bins)
    ax2.hist(data_y_bcond, bins=count_bins.size, range=(bin_start, bin_end), alpha=transp_histogram, log=log_histogram, edgecolor='k', label=label)
    # ax2_in.plot(data_x_bcond/period_time, data_y_bcond, label=label)
    
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()
Last_Offset = 15
data_aconds = []
for i, data_file in enumerate(data_files_AfterCond):
    data_acond = getData(data_file)
    if i == 1:
        data_acond[:,0] -= Last_Offset
        data_acond = data_acond[data_acond[:,0] > 0]
        label = f"After x4 conditioning fluence"
    else:
        label = f"After x3 conditioning fluence"
    Condition_Fluence = Condition_Fluences[i]
    data_x_acond = data_acond[:,0]
    data_y_acond = data_acond[:,1]
    data_aconds.append(data_acond)
    hist_bins = generateHistogramBins(data_y_acond, count_bins)
    ax2.hist(data_y_acond, bins=count_bins.size, range=(bin_start, bin_end), alpha=transp_histogram - (i+1)*reduction_transp_histogram, log=log_histogram, edgecolor='k', label=label)
        
new_start = 60
new_end = 550
data_y_acond_add = data_longsweep[new_start:new_end, 1]
data_x_acond_add = data_longsweep[new_start:new_end, 0]
# data_x_acond_add -= min(data_x_acond_add)
data_x_acond_add += max(data_x_acond) - min(data_x_acond_add) + 0.2

data_x_acond = np.hstack((data_x_acond, data_x_acond_add))
data_y_acond = np.hstack((data_y_acond, data_y_acond_add))
ax2_in.plot(data_x_acond/period_time, data_y_acond, color=colors['color'][2], label=f'After x4 conditioning fluence, cont.')


first_dashed = True
step_color = "tab:purple"
for i, time in enumerate(step_timings):
    label = None
    if i == 0:
        pass
    if i%steps_per == 0:
        vline_style='solid'
        if time == 0:
            label = 'Sweep start $s_0$'
    else:
        vline_style ='dotted'
        if first_dashed:
            label = 'Sweep steps $s_1$, $s_2$, ...'
            first_dashed = False
    ax2_in.vlines(time/period_time, -10, 300, linestyle=vline_style, label=label, color=step_color, alpha=0.5)
    
ax2_in.text(0.97, 258, f'$s_0$', fontsize='14', color=step_color)
ax2_in.text(1.05, 248, f'$s_1$', fontsize='14', color=step_color)
ax2_in.text(1.07, 238, f'$s_2$', fontsize='14', color=step_color)
ax2_in.hlines(BG_Ave, 0, 5, label='Background average', color='k')

with mpl.cbook.get_sample_data(r'C:\Users\brend\Desktop\Repositories\Selectivity-Paper_2021\Paper_Data\Pitting-Conditioning_SacrificialArea\sweep_inset.png') as file:
    im = plt.imread(file, format='png')
ax2_in2 = ax2.inset_axes([150,45,150,150],transform=ax2.transData)
ax2_in2.imshow(im)
ax2_in2.axis('off')

ax2.legend(fontsize=15)
ax2.set_xlabel('Counts (photons)',fontsize=30, labelpad=20)
ax2.set_ylabel('Frequency',fontsize=30, labelpad=23)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set_xlim(0, 250)
# ax2.set_ylim(-2, 255)
ax2_in.legend(fontsize=15)
ax2_in.set_xlabel('Sweep periods',fontsize=30, labelpad=20)
ax2_in.set_ylabel('Counts (a.u.)',fontsize=30, labelpad=23)
ax2_in.tick_params(axis='both', which='major', labelsize=20)
ax2_in.set_xlim(0, 3.5)
ax2_in.set_ylim(-2, 255)
ax2_in.yaxis.labelpad = 4
ax2_in.xaxis.labelpad = 4




plt.show()