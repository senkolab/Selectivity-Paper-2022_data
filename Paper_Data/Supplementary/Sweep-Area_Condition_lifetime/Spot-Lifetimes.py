# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 17:40:06 2021

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
ax2 = fig.add_subplot(111)

PulsesPer = 120
DownSample = 3
data = getData(data_file_LongTerm, PulsesPer, long=True)
data_y = data[:,1]
data_y = data_y[:-1]
data_yD, data_yDsdv = Down_Sample(data_y, 3)
data_x = np.arange(data_yD.size)*DownSample


ax2.plot(data_x, data_yD, label="Sweep", alpha=0.3, color='grey')


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

Down = 10
Pit1A, Pit1SD = Down_Sample(Pit1, Down)
Pit2A, Pit2SD = Down_Sample(Pit2, Down)
Pit3A, Pit3SD = Down_Sample(Pit3, Down)
Pit4A, Pit4SD = Down_Sample(Pit4, Down)

Scale1 = 0.7
Scale2 = 0.4
Scale3 = 0.8
Scale4 = 0.6
Scale1 = Scale2 = Scal3 = Scale4 = 1

ax2.plot(np.arange(0, Pit1A.size*Down*Pit1Num, Down*Pit1Num), Pit1A*Scale1, label='Spot 1', alpha=0.85)
ax2.plot(np.arange(0, Pit2A.size*Down*Pit2Num, Down*Pit2Num), Pit2A*Scale2, label='Spot 2', alpha=0.85)
ax2.plot(np.arange(0, Pit3A.size*Down*Pit2Num, Down*Pit2Num), Pit3A*Scale3, label='Spot 3', alpha=0.85)
ax2.plot(np.arange(0, Pit4A.size*Down*Pit2Num, Down*Pit2Num), Pit4A*Scale4, label='Spot 4', alpha=0.85)

average_lifetime = np.mean([Pit1A.size*Down*Pit1Num, Pit2A.size*Down*Pit2Num, Pit3A.size*Down*Pit2Num, Pit4A.size*Down*Pit2Num])

ax2.vlines(average_lifetime, 0, 200, label='Average spot lifetime', color='k', linestyle='dotted')


#Need to get new BG
plotBackground(ax2, 0.55)
ax2.legend(fontsize=20)
ax2.set_xlabel('Pulse number',fontsize=30, labelpad=20)
ax2.set_ylabel('Counts (a.u.)',fontsize=30, labelpad=23)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set_xlim(0, 60000)
ax2.set_ylim(0, 150)