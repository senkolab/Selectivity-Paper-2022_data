# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 13:23:54 2021

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
import re
import json

Ba137Abundance = 0.113

def Sanitize_SelectivityData(data_files, new_file):
    with open(new_file + '.txt', 'w+') as new:
        for i, data_file in enumerate(data_files):
            with open(data_file, 'r') as file:
                for line in file:
                    Words = line.split("\t")
                    if Words[0].isdigit():
                        new.write(line)
    return glob.glob(new_file + '.txt')[0]

def getData(data_file):
    raw_data = np.asarray(pd.read_csv(data_file, delimiter = "[\],\[\t]", engine='python',header=None))
    raw_data[:,2] = np.round(raw_data[:,2].astype(np.float), 6)#Round laser freq. data
    return raw_data

def getMetaData(data_file):
    Metadata = ''
    Header = ''
    with open(data_file, 'r') as file:
        for line in file:
            Words = line.split(",")
            # print(Words)
            if Words[0] == 'Metadata':
                Metadata = Words[1:]
                break
        for line in file:
            Words = line.split("\t")
            if Words[0][0] == 'h':
                Header = Words[1:]
                break
    return Metadata, Header

def countIons(Counts, IonCounts1, IonCountsStd):
    IonsTrapped = np.zeros(Counts.shape)
    for i in range(6):
        Condition = np.array([(IonCounts1 - IonCountsStd)*i < Counts, (IonCounts1 + IonCountsStd)*i > Counts])
        Array = np.all(Condition, axis=0)*i
        IonsTrapped += Array
        # print(np.vstack((Counts, Condition, IonsTrapped)).T)
    return IonsTrapped

def findSingleIonAverage(Counts, BG, BGStd, PlotIonCounts = False):
    CountsNonBG = Counts[Counts > BG + 5*BGStd]
    # print(f'BG:{BG}, BGStd:{BGStd}')
    SmallestNonBGAve = np.mean(CountsNonBG[np.argpartition(CountsNonBG, 3)[:4]])
    CountsSingle = CountsNonBG[CountsNonBG-BG < 2*(SmallestNonBGAve-BG)]
    print(f'SingleCutoff:{2*(SmallestNonBGAve-BG) + BG}')
    AveSingle = np.mean(CountsSingle)
    SteSingle = np.std(CountsSingle)
    if PlotIonCounts:
        fig = plt.figure(figsize=(20,10))
        ax1 = fig.add_subplot(111)
        ax1.plot(CountsNonBG)
        for i in range(1, 5):
            IonsCountBase = i*(AveSingle-BG) + BG
            ax1.plot([0, CountsNonBG.size], [IonsCountBase, IonsCountBase])
    return AveSingle, SteSingle


data_files = "Load-Eff-Spect_Ba137_140uJ_80uW.*"
new_file = "Load-Eff-Spect_Ba137_140uJ_80uW_Sanitized"
data_files = glob.glob(data_files)
# new_file = Sanitize_SelectivityData(data_files, new_file)
new_file = "Load-Eff-Spect_Ba137_140uJ_80uW_Sanitized.txt"

meta_data, Header = getMetaData(data_files[0])
Header = [i.replace('\n', '') for i in Header]
Sweepdata = [float(i) for i in re.split(' |:',meta_data[0])[1:5]]
SweepFreq0 = Sweepdata[-1]
SweepRange = np.arange(Sweepdata[0], Sweepdata[1]+1e-7, Sweepdata[2])
SweepFreqs = SweepFreq0 + np.arange(Sweepdata[0], Sweepdata[1]+1e-7, Sweepdata[2])
CoolingFreq, RepumpFreq, WindowStart, WindowWidth = [float(i.split(':')[1]) for i in meta_data[1:]]
Isotope, PulseEnergy, LaserPower, = re.split("_", re.sub("[^0-9|^_]", "", new_file)[1:-1])

raw_data = getData(new_file)
BG = np.mean(raw_data[:,-2])
BGstd = np.mean(raw_data[:,-1])
IonCounts1, IonCountsStd = findSingleIonAverage(raw_data[:,5], BG, BGstd)
print(f'IonCounts1:{IonCounts1}, IonCountsStd:{IonCountsStd}')
IonCounts1 = 200
IonCountsStd = 100


fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(111)

FreqDisplayOffset = 541.433
plot_x = (SweepFreqs - FreqDisplayOffset)*1e6
plot_y = np.array([])
plot_y_std = np.array([])
plot_y2 = np.array([])
plot_y2_std = np.array([])

for i, Freq in enumerate(SweepRange):
    data = raw_data[i::SweepRange.size, :]
    Counts = data[:,5]
    BG = np.mean(data[:,-2])
    Counts -= BG
    SingleCounts = IonCounts1 - BG
    IonsTrapped = countIons(Counts, SingleCounts, IonCountsStd)
    # print(np.hstack((data, IonsTrapped[:,None])))
    LoadingEfficiency = np.mean(IonsTrapped)
    LoadingEfficiencySte = np.std(IonsTrapped)/np.sqrt(IonsTrapped.size)
    plot_y = np.append(plot_y, LoadingEfficiency)
    plot_y_std = np.append(plot_y_std, LoadingEfficiencySte)
    plot_y2 = np.append(plot_y2, np.mean(data[:, 4]))
    plot_y2_std = np.append(plot_y2_std, np.std(data[:, 4])/np.sqrt(data.shape[0]))

ax1.errorbar(plot_x, plot_y, yerr=plot_y_std, label=f'Ba{Isotope} Loading Efficiency')
ax1.axvline(39+275, color='k', label='Spectrum Peak')
ax2 = ax1.twinx()
ax2.errorbar(plot_x, plot_y2, yerr=plot_y2_std, color='r')
WriteString = f"""
PulseEnergy: {PulseEnergy}
553nmPower: {LaserPower}
CoolingFreq: {CoolingFreq}
RepumpFreq: {RepumpFreq}
WindowStart:{WindowStart}
WindowWidth:{WindowWidth}
DataPoints:{data.shape[0]}"""
ax2.text(320, 22.5, WriteString, fontsize=14)

ax1.set_xlabel(f'Frequency / MHz + {FreqDisplayOffset}',fontsize=30)
ax1.set_ylabel('Loading efficiency', fontsize=30)
ax2.set_ylabel('Neutral Fluorescence (a.u.)',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax1.legend(loc='best', fontsize=14)

plt.show()