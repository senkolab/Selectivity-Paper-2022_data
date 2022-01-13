# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 15:04:51 2021

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

def getData(data_file, NumSamples = 1):
    raw_data = np.asarray(pd.read_csv(data_file, delimiter = "[\],\[\t]", engine='python',header=None))
    if NumSamples == 1:
        # print(raw_data)
        counts = raw_data[:,2]
    else:
        raw_counts = raw_data[:,2:2+NumSamples]
        TotalPulses = raw_counts.size
        counts = np.mean(raw_counts, 1)
    
    times = raw_data[:,NumSamples + 17]
    times -= min(times)
    times = FixTimeOffset(times)
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

data_file = r"Ablation-Sweep_Nat-1-Sacrif1 2020-11-04 16-58-30*"
data_file = glob.glob(data_file)
data_file = data_file[0]
step_timings, steps_per = getSweepTimings(data_file)
period_time = step_timings[steps_per]

data_BG = getData(data_file_BG)
BG_Ave = np.mean(data_BG[:,1])

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(111)

for i, data_file in enumerate(data_files_BeforeCond):
    data = getData(data_file)
    data_x = data[:,0]
    data_y = data[:,1]
    ax1.plot(data_x/period_time, data_y, label='Before conditioning')
    
Last_Offset = 15
for i, data_file in enumerate(data_files_AfterCond):
    data = getData(data_file)
    if i == 1:
        data[:,0] -= Last_Offset
        data = data[data[:,0] > 0]
        label = f"After x3 conditioning fluence"
    else:
        label = f"After x2 conditioning fluence"
    Condition_Fluence = Condition_Fluences[i]
    data_x = data[:,0]
    data_y = data[:,1]
    ax1.plot(data_x/period_time, data_y, label=label)

plotBackground(ax1, BG_Ave)


first_solid = True
first_dashed = True
step_color = "tab:purple"
for i, time in enumerate(step_timings):
    label = None
    if i == 0:
        pass
    if i%steps_per == 0:
        vline_style='solid'
        if first_solid:
            label = 'Sweep end'
            first_solid = False
        elif i < 20:
            ax1.text((time + 0.1)/period_time, 200, 'Start', fontsize='10', color=step_color)
    else:
        vline_style ='dotted'
        if first_dashed:
            label = 'Sweep steps'
            first_dashed = False
        if i < 3:
            ax1.text((time + 0.1)/period_time, 200, f'$s_{i}$', fontsize='10', color=step_color)
    ax1.vlines(time/period_time, -10, 300, linestyle=vline_style, label=label, color=step_color)
    
with mpl.cbook.get_sample_data(r'C:\Users\brend\Desktop\Repositories\Selectivity-Paper_2021\Paper_Data\Pitting-Conditioning_SacrificialArea\sweep_inset.png') as file:
    im = plt.imread(file, format='png')
ax1_in = ax1.inset_axes([1,50,1,130],transform=ax1.transData)
ax1_in.imshow(im)
ax1_in.axis('off')

ax1.legend(fontsize=15)
ax1.set_xlabel('Sweep periods',fontsize=30, labelpad=20)
ax1.set_ylabel('Counts (a.u.)',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlim(0, x_max/period_time)
ax1.set_ylim(-2, 255)



plt.show()





