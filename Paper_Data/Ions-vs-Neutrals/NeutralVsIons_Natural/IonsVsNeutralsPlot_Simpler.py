# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 20:06:16 2020

@author: brend
"""

import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import math
import re
plt.style.use('seaborn-colorblind')

#Function to get all data from data files
def GetData(Path, Labels=[], Plot=False, OverrideLength=500):
    data_files = glob.glob(Path)
    data_files = np.asarray(data_files)
    if Plot:
        figure(figsize=(30,12))
    data_Avgs = np.array([])
    data_StDevs = np.array([])
    data_Array = np.array([])
    data_Raw = np.array([])
    data_Times_Raw = np.array([])
    data_Energies = np.array([])
    for i in range(0, data_files.size):
        data = pd.read_csv(data_files[i], delimiter = "[\],\[]", engine='python',header=None)
        data = np.asarray(data)
        filename = data_files[i]
        filename_split = re.split(r'_', filename)
        energy = filename_split[1]
        energy = float(re.split(r'uJ', energy)[0])
        
        data_Energies = np.append(data_Energies, energy)
        dataCounts = data[:min(data.shape[0], OverrideLength),2]
        dataTimes= data[:min(data.shape[0], OverrideLength), 18]
        if i == 0:
            data_Times_Raw = dataTimes
            data_Raw = dataCounts
        else:
            data_Times_Raw = np.vstack((data_Times_Raw, dataTimes))
            data_Raw = np.vstack((data_Raw, dataCounts))
        Ave = np.mean(dataCounts)
        Var = np.std(dataCounts)
        data_Avgs = np.append(data_Avgs, Ave)
        data_StDevs = np.append(data_StDevs, Var)
        if Plot:
            if len(Labels) != 0:
                plt.plot(np.arange(0,dataCounts.size), dataCounts, label='%f' %Labels[i], alpha=0.5)
            else:
                plt.plot(np.arange(0,dataCounts.size), dataCounts, alpha=0.5)
    if len(Labels) != 0:
        plt.legend()
    
    return (data_Avgs, data_StDevs, data_Energies, data_Times_Raw, data_Raw)

def GetPressureData(Path, Labels=[], Plot=False, OverrideLength=150):
    data_files = glob.glob(Path)
    data_files = np.asarray(data_files)
    if Plot:
        figure(figsize=(30,12))
    data_Avgs = np.array([])
    data_StDevs = np.array([])
    data_Array = np.array([])
    data_Raw = np.array([])
    data_Energies = np.array([])
    for i in range(0, data_files.size):
        data = pd.read_csv(data_files[i], delimiter = "\t", engine='python',header=None)
        data = np.asarray(data)
        filename = data_files[i]
        filename_split = re.split(r'_', filename)
        energy = filename_split[2]
        energy = float(re.split(r'uJ', energy)[0])
        
        data_Energies = np.append(data_Energies, energy)
        dataTimes= data[:min(data.shape[0], OverrideLength),0]
        data = data[:min(data.shape[0], OverrideLength),1].astype(float)
        if i == 0:
            data_Times_Raw = dataTimes
            data_Raw = data
        else:
            data_Times_Raw = np.vstack((data_Times_Raw, dataTimes))
            data_Raw = np.vstack((data_Raw, data))
        Ave = np.mean(data)
        Var = np.std(data)
        data_Avgs = np.append(data_Avgs, Ave)
        data_StDevs = np.append(data_StDevs, Var)
        if Plot:
            if len(Labels) != 0:
                plt.plot(np.arange(0,data.size), data, label='%f' %Labels[i], alpha=0.5)
            else:
                plt.plot(np.arange(0,data.size), data, alpha=0.5)
    if len(Labels) != 0:
        plt.legend()
    
    return (data_Avgs, data_StDevs, data_Energies, data_Times_Raw, data_Raw)

Neutral_data_path = r"NeutralFluorescence_*"
Neutral_Aves, Neutral_StDevs, Neutral_data_energies, Neutral_times_raw, Neutral_data_raw = GetData(Neutral_data_path, Plot=True)
Neutral_data_fluences = Neutral_data_energies*1e-6/(math.pi*(98*1e-4)**2)
Neutral_Maxes = np.amax(Neutral_data_raw, 1)
Neutral_Vars = Neutral_StDevs/np.sqrt(2) #2 full sweeps

Ion_data_path = r"IonFluorescence_*"
Ion_Aves, Ion_StDevs, Ion_data_energies, Ion_times_raw, Ion_data_raw = GetData(Ion_data_path)
Ion_data_fluences = Ion_data_energies*1e-6/(math.pi*(98*1e-4)**2)
Ion_Maxes = np.amax(Ion_data_raw, 1)
Ion_Vars = Ion_StDevs/np.sqrt(2) #2 full sweeps

BG_data_path = r"BG_*"
BG_Aves, BG_StDevs, BG_data_energies, BG_times_raw, BG_data_raw = GetData(BG_data_path, OverrideLength = 50)

Neutral_Aves = Neutral_Aves - BG_Aves
Neutral_Maxes = Neutral_Maxes - BG_Aves
Ion_Aves = Ion_Aves - BG_Aves
Ion_Maxes = Ion_Maxes - BG_Aves

Pressure_data_path = r"Pressure_*"
Pressure_Aves, Pressure_StDevs, Pressure_data_energies, Pressure_times_raw, Pressure_data_raw = GetPressureData(Pressure_data_path)
Pressure_data_fluences = Pressure_data_energies*1e-6/(math.pi*(98*1e-4)**2)
Pressure_Maxes = np.amax(Pressure_data_raw, 1)
Pressure_Vars = Pressure_StDevs/np.sqrt(Pressure_data_raw.shape[1]) #2 full sweeps


fig1, ax1 = plt.subplots(1, figsize=(15,15))
plott = ax1.errorbar(Neutral_data_fluences, Neutral_Aves, yerr=Neutral_Vars, label='Neutral fluor. ave.', fmt='o', markersize=13)
plotcolor = plott.get_children()[0].get_color()
# ax1.plot(Neutral_data_fluences, Neutral_Maxes, label='Neutral fluor. max', linestyle='', marker='s', color=plotcolor, markersize=10)
plott = ax1.errorbar(Ion_data_fluences, Ion_Aves, yerr=Ion_Vars, label='Ions ave.', fmt='o', markersize=13)
plotcolor = plott.get_children()[0].get_color()
# ax1.plot(Ion_data_fluences, Ion_Maxes, label='Ions Max', linestyle='', marker='s', color=plotcolor, markersize=10)
plt.legend(fontsize=25, loc=2)

# ax2 = ax1.twinx()
# plott = ax2.errorbar(Pressure_data_fluences, Pressure_Aves, yerr=Pressure_Vars, label='Pressure Average', fmt='o', color='#D55E00', markersize=10)
# plotcolor = plott.get_children()[0].get_color()
# ax2.plot(Pressure_data_fluences, Pressure_Maxes, label='Pressure Max', linestyle='', marker='s', color=plotcolor, markersize=10)
# plt.legend(fontsize=25, loc=1)

ax1.tick_params(axis='both', which='major', labelsize=30)
ax1.set_ylim(-10, 500)
ax1.set_xlabel(r'Fluence / $J\cdot cm^{-2}$',fontsize=35)
ax1.set_ylabel('Fluorescence Counts',fontsize=35)

# from matplotlib import ticker
# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True) 
# formatter.set_powerlimits((-1,1)) 
# ax2.yaxis.set_major_formatter(formatter) 
# ax2.yaxis.offsetText.set_fontsize(25)
# ax2.set_ylabel('Pressure / mbar', fontsize=35)
# ax2.tick_params(axis='both', which='major', labelsize=30)
# ax2.set_ylim(0.9e-10, 1e-9)
fig1.tight_layout()

