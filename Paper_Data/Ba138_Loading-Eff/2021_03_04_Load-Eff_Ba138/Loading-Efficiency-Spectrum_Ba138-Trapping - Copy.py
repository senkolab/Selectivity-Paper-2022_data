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

import sys
sys.path.append('C:\\Users\\brend\\Desktop\\Git-Repos\\Ablation-Paper_2020\\Data\\')
from DataAnalysis_Functions import *


Verbose = True

#IonCounts, IonCountsStd, BG, BGStd, NeutralCounts
data_pos = [5, 6, 7, 8, 4]

                                                
data_files = "Load-Eff-Spect_Ba138_100uJ_8uW_?.*"

new_file, meta_data, Header, LinesEachFile = prepareData(data_files)
raw_data_full = getSanitizedData(new_file, LinesEachFile=LinesEachFile)
WriteString = parametersWriteString(meta_data, Verbose=Verbose)
IonsTrapped, LoadingEff, LoadingEffs, NeutralFluorescence, WriteString = analyzeLoading(raw_data_full, data_pos, WriteString, Verbose=Verbose)
# WriteString += f'\nLoading Efficiencies:{LoadingEffs[0][0]*100:0.1f}%, {LoadingEffs[1][0]*100:0.1f}%, {LoadingEffs[2][0]*100:0.1f}%'
neutral_pos = data_pos[-1]
fig, ax1, ax2 = plotIons_Neutrals_VsPulse(raw_data_full, neutral_pos, IonsTrapped, WriteString, [0, 1])

label = 'Ionization Freqs.'
units_factor = 1e-6
scan_array = constructScan(meta_data, label, units_factor)
offset = 541.433
scan_array -= offset
scan_array *= 1e6
xlabel = f'Frequency - {offset} THz / MHz'
fig, ax1, ax2 = plotEfficiency_Neutrals_VsScan(raw_data_full, neutral_pos, scan_array, IonsTrapped, NeutralFluorescence, xlabel)