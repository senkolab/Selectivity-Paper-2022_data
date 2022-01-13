# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 15:04:51 2021

@author: brend
"""


import math
import numpy as np
import matplotlib.pyplot as plt
#print(plt.style.available)
plt.style.use('seaborn-colorblind')
# plt.rcParams['figure.facecolor'] = '1'

import glob
import pandas as pd

FocusWaist = 100*1e-6

def getData(data_file, NumSamples = 1):
    raw_data = np.asarray(pd.read_csv(data_file, delimiter = "[\],\[\t]", engine='python',header=None))
    if NumSamples == 1:
        counts = raw_data[:,2]
    else:
        raw_counts = raw_data[:,2:2+NumSamples]
        TotalPulses = raw_counts.size
        print(TotalPulses)
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
    ax.plot(BG_x, BG_y, label='Average background')
    
data_filepath = r"Neutral_Conditioning_Sacrif1*"
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

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(111)

for i, data_file in enumerate(data_files_BeforeCond):
    data = getData(data_file)
    data_x = data[:,0]/60
    data_y = data[:,1]
    ax1.plot(data_x, data_y)
    
Last_Offset = 0.18
for i, data_file in enumerate(data_files_AfterCond):
    data = getData(data_file)
    if i == 1:
        data[:,0] -= Last_Offset
        data = data[data[:,0] > 0]
    Condition_Fluence = Condition_Fluences[i]
    data_x = data[:,0]/60
    data_y = data[:,1]
    ax1.plot(data_x, data_y, label=f"{Condition_Fluence/(100**2):0.3f} $\mu J/cm^2$ Conditioning pulse energy")

plotBackground(ax1, BG_Ave)
ax1.legend(fontsize=20)
ax1.set_xlabel('Time / minute',fontsize=30, labelpad=20)
ax1.set_ylabel('Counts (a.u.)',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlim(0, 2)

plt.show()

fig = plt.figure(figsize=(20, 10))
ax2 = fig.add_subplot(111)

PulsesPer = 120
data = getData(data_file_LongTerm, PulsesPer)
data_y = data[:,1]
data_x = np.arange(data_y.size)*PulsesPer
ax2.plot(data_x, data_y, label="Average fluorescence each sweep")
#Need to get new BG
plotBackground(ax2, BG_Ave)
ax2.legend(fontsize=20)
ax2.set_xlabel('Pulse number',fontsize=30, labelpad=20)
ax2.set_ylabel('Counts (a.u.)',fontsize=30, labelpad=23)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set_xlim(0, max(data_x))