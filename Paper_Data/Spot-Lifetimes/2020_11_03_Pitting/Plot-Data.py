# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 08:50:26 2020

@author: brend
"""

import math
import glob
import pandas as pd
import numpy as np
from sympy import *
import matplotlib.pyplot as plt



#Function to get all data from data files
def GetPittData(Path, Plot=False):
    data_files = glob.glob(Path)
    data_files = np.asarray(data_files)
    if Plot:
        plt.figure(figsize=(20,12))
    data_Avgs = np.array([])
    TotalCols = 22
    for i in range(0, data_files.size):
        data = pd.read_csv(data_files[i], delimiter = "[\],\[]", engine='python',header=None)#
        data = np.asarray(data)
        dataCounts = data[:data.shape[0],2:data.shape[1]-TotalCols+3]
        dataCounts = dataCounts.flatten()
        dataCounts = dataCounts.astype(int)
        dataCountsAvg = dataCounts
        data_Avgs = np.append(data_Avgs, dataCountsAvg)
        if Plot:
            plt.plot(np.arange(0,dataCountsAvg.size), dataCountsAvg, label='%i' %i, alpha=0.5)
    if Plot:
        plt.legend()
    return (data_Avgs)
#Function to get all data from data files
def GetPittDataOld(Path, Plot=False):
    data_files = glob.glob(Path)
    data_files = np.asarray(data_files)
    if Plot:
        plt.figure(figsize=(20,12))
    data_Avgs = np.array([])
    for i in range(0, data_files.size):
        data = pd.read_csv(data_files[i], delimiter = "\t", engine='python',header=1)#
        data = np.asarray(data)
        dataCounts = data[:data.shape[0],1]
        dataCounts = dataCounts.flatten()
        dataCounts = dataCounts.astype(int)
        dataCountsAvg = dataCounts
        data_Avgs = np.append(data_Avgs, dataCountsAvg)
        if Plot:
            plt.plot(np.arange(0,dataCountsAvg.size), dataCountsAvg, label='%i' %i, alpha=0.5)
    if Plot:
        plt.legend()
    return (data_Avgs)



data_file = r"Neutral-Fluorescence_Natural_Pitting_BG_No-Ablation_*"
BG = GetPittData(data_file)
BG = np.mean(BG)

data_file = r"..\2020_02_12\RawPittingData*"
PitOld = GetPittDataOld(data_file)

data_file = r"Pit1\Neutral-Fluorescence_Natural_Pitting_48uJ_7p26x6p13_After-Condition-1_*"
Pit1 = GetPittData(data_file)

data_file = r"Pit2\Neutral-Fluorescence_Natural_Pitting_48uJ_7p265x6p13_After-Condition-*"
Pit2 = GetPittData(data_file)

data_file = r"..\2020_11_04_Pitting-Conditioning\Pit3\Neutral_Pitting_7p28x6p125_After-Condition-1_*"
Pit3 = GetPittData(data_file)

data_file = r"..\2020_11_04_Pitting-Conditioning\Pit4\Neutral_Pitting_7p265x6p15_After-Condition-2_*"
Pit4 = GetPittData(data_file)

PitOld = PitOld - min(PitOld) - 10
PitOld = PitOld*140/max(PitOld) + BG

DownSample = 10

Rem = PitOld.size%(int(DownSample/2))
if Rem != 0:
    PitOldD = PitOld[:-Rem]
else:
    PitOldD = PitOld
Rem = Pit1.size%DownSample
if Rem != 0:
    Pit1D = Pit1[:-Rem]
else:
    Pit1D = Pit1
Rem = Pit2.size%DownSample
if Rem != 0:
    Pit2D = Pit2[:-Rem]
else:
    Pit2D = Pit2
Rem = Pit3.size%DownSample
if Rem != 0:
    Pit3D = Pit3[:-Rem]
else:
    Pit3D = Pit3
Rem = Pit4.size%DownSample
if Rem != 0:
    Pit4D = Pit4[:-Rem]
else:
    Pit4D = Pit4

PitOldA = np.mean(PitOldD.reshape(-1, int(DownSample/2)), 1)
Pit1A = np.mean(Pit1D.reshape(-1, DownSample), 1)
Pit2A = np.mean(Pit2D.reshape(-1, DownSample), 1)
Pit3A = np.mean(Pit3D.reshape(-1, DownSample), 1)
Pit4A = np.mean(Pit4D.reshape(-1, DownSample), 1)

PitOldSD = np.std(PitOldD.reshape(-1, int(DownSample/2)), 1)
Pit1SD = np.std(Pit1D.reshape(-1, DownSample), 1)
Pit2SD = np.std(Pit2D.reshape(-1, DownSample), 1)
Pit3SD = np.std(Pit3D.reshape(-1, DownSample), 1)
Pit4SD = np.std(Pit4D.reshape(-1, DownSample), 1)

fig , ax1 = plt.subplots(1, figsize=(20, 10))
ax1.plot([0, 13000], [BG, BG], label='BG')

ax1.errorbar(np.arange(0, PitOldA.size*DownSample*2, DownSample*2), PitOldA, PitOldSD/np.sqrt(DownSample), label='PitOld', alpha=0.5)
ax1.errorbar(np.arange(0, Pit1A.size*DownSample, DownSample), Pit1A, Pit1SD/np.sqrt(DownSample), label='Pit1', alpha=0.5)
ax1.errorbar(np.arange(0, Pit2A.size*DownSample, DownSample), Pit2A, Pit2SD/np.sqrt(DownSample), label='Pit2', alpha=0.5)
ax1.errorbar(np.arange(0, Pit3A.size*DownSample, DownSample), Pit3A, Pit3SD/np.sqrt(DownSample), label='Pit3', alpha=0.5)
ax1.errorbar(np.arange(0, Pit4A.size*DownSample, DownSample), Pit4A, Pit4SD/np.sqrt(DownSample), label='Pit4', alpha=0.5)


ax1.set_xlim(6000, 7500)
ax1.set_ylim(10, 55)
ax1.set_xlabel('Pulse Number',fontsize=30, labelpad=20)
ax1.set_ylabel('Counts (a.u.)',fontsize=30, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)

plt.legend(fontsize=20)
plt.show()

PercentageError = (Pit2SD/np.sqrt(DownSample))/Pit2A
PercentageErrorAve = np.mean(PercentageError)
print('DownSample: %i, Error Percentage: %f%%'%(DownSample, PercentageErrorAve*100))

def DownSampleErrors(PitData):
    for DownSample in range(3, 40):
        Rem = PitData.size%(int(DownSample))
        if Rem != 0:
            PitD = PitData[:-Rem]
        else:
            PitD = PitData
        
        PitA = np.mean(PitD.reshape(-1, int(DownSample)), 1)
        
        PitSD = np.std(PitD.reshape(-1, int(DownSample)), 1)

        if DownSample == 3:
            print(PitA)
            print(PitSD)
        PercentageError = (PitSD/np.sqrt(DownSample))/PitA
        PercentageErrorAve = np.mean(PercentageError)
        print('DownSample: %i, Error Percentage: %f%%'%(DownSample, PercentageErrorAve*100))
        
DownSampleErrors(Pit1)
