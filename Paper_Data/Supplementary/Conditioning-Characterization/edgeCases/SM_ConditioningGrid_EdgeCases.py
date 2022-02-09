from cProfile import label
from matplotlib.transforms import Bbox
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.arrayprint import format_float_scientific
from scipy.optimize import curve_fit
import glob
from numpy.core.defchararray import index
import pandas as pd
import statistics
from matplotlib import rc,rcParams
import matplotlib

# matplotlib.rcParams.update(matplotlib.rcParamsDefault)

import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
# graph_edge_width = 0.5
mpl.rcParams.update(mpl.rcParamsDefault)
# mpl.rcParams['axes.linewidth'] = graph_edge_width
# mpl.rcParams.update({'errorbar.capsize': 2})
plt.style.use('seaborn-colorblind')

low_color = 'tab:blue'
high_color = 'tab:purple'
mod_color = 'tab:orange'

## Data Processing ## 

def getDataCounts_General(data_path, start): 
    ## Import the data 
    data_file = np.asarray(glob.glob(data_path))[0]
    data = np.asarray(pd.read_csv(data_file, delimiter = "[\]\[\t]", engine='python', header=None))
    
    counts = data[:,2] ## third row contains the actual PMT counts at each pulse, 
                       ##collected in groups of 10x pulses
    rows = counts.size ## number rows = no. runs

    floats = []
    for i in range(0,rows):
        floats.append([float(x) for x in counts[i].split()]) ## convert from strings to floats

    ## Mitigating for doubles of low/high counts by averaging rows
    tmp = np.zeros((len(floats), 10))
    for i in range(len(floats)):
        for k in range(10):
            tmp[i][k] = floats[i][k]


    floats_low = []
    floats_high = []
    low = 0
    high = 0 
    xlow = []
    xhigh =[]

    rows = len(floats)
    for i in range(0,rows-1,2):
        if start == 'low':
            floats_high.append(floats[i+1])
            floats_low.append(floats[i])
            low += 10 
            high += 10 
        else: 
            floats_high.append(floats[i])
            floats_low.append(floats[i+1]) 
            low += 10 
            high += 10 
        xlow.append(low)
        xhigh.append(high)

    if (start == 'low' and rows%2 != 0): #check if odd
        floats_low.append(floats[rows-1])
        low += 10 
        xlow.append(low)
    if (start == 'high' and rows%2 != 0): #check if odd
        floats_high.append(floats[rows-1])
        high += 10 
        xhigh.append(high)

    floats_low_array = np.array(floats_low)
    floats_high_array = np.array(floats_high)

    low_counts = floats_low_array.astype(np.int)
    high_counts = floats_high_array.astype(np.int)

    return tmp, low_counts, high_counts, np.array(xlow), xhigh

def getAverageCountsLow(low_counts):
    '''
    Averages across the 10 pulses of each row 
    Calculates standard deviation for each row 
    Returns averages, standard deviations 
    '''
    pulses = 10 
    rows = low_counts.shape[0]

    summ = np.zeros(rows)
    avgs = np.zeros(rows)
    sdev = np.zeros(rows)

    for i in range(rows):
        for j in range(pulses):
            summ[i] += low_counts[i][j]
        avgs[i] = summ[i] / pulses
        sdev[i] = statistics.stdev(low_counts[i])

    return avgs, sdev

## Fits : Defining the exponential function form ##

def func(x, a, b, c):
    '''
    Returns fit function 
    '''
    return -a * np.exp(-b * x) + c

## Spots ##
#data paths to all the raw data of the spots data imported 
data_paths = [
    r'Neut-Fluor-Spect_63uJ_145uJ_7o32x_6o21y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o32x_6o14y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o30x_6o08y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o29x_6o20y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o29x_6o18y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o28x_6o18y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o26x_6o02y.txt',## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o25x_6o18y.txt',## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o22x_5o96y.txt',## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o17x_5o90y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o14x_6o04y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o13x_6o02y.txt',## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o11x_6o00y.txt',## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o11x_5o88y.txt',## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o10x_5o90y.txt',## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o09x_5o86y.txt',
    r'Neut-Fluor-Spect_63uJ_145uJ_7o08x_5o92y.txt',
    r'Neut-Fluor-Spect_63uJ_145uJ_7o08x_5o84y.txt',
    r'Neut-Fluor-Spect_63uJ_145uJ_7o07x_5o90y.txt'
] 
##### Spot 1 #########################
data_path_1 = data_paths[0]
test_1, low_1, high_1, xlow_1, xhigh_1 = getDataCounts_General(data_path=data_path_1, start='low')
avgsLow_1, stdLow_1 = getAverageCountsLow(low_1)
avgsHigh_1, stdHigh_1 = getAverageCountsLow(high_1)
popt_1, pcov_1 = curve_fit(func, xlow_1, avgsLow_1, p0=[120, 0.05, 120], maxfev=10000)
## Mitigate: NONE, LEFT, RIM 
## 
##
## 
######################################
##### Spot 2 #########################
data_path_2 = data_paths[1]
test_2, low_2, high_2, xlow_2, xhigh_2 = getDataCounts_General(data_path=data_path_2, start='high')
avgsLow_2, stdLow_2 = getAverageCountsLow(low_2)
avgsHigh_2, stdHigh_2 = getAverageCountsLow(high_2)
popt_2, pcov_2 = curve_fit(func, xlow_2, avgsLow_2, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate: BOTTOM
xlow_2[29] += 10
xlow_2[100] += 10
for i in range(30,len(xlow_2)):
    xlow_2[i] += 10
for i in range(101, len(xlow_2)):
    xlow_2[i] += 10
# ##
##
##
## 
######################################
##### Spot 3 #########################
data_path_3 = data_paths[2]
test_3, low_3, high_3, xlow_3, xhigh_3 = getDataCounts_General(data_path=data_path_3, start='low')
avgsLow_3, stdLow_3 = getAverageCountsLow(low_3)
avgsHigh_3, stdHigh_3 = getAverageCountsLow(high_3)
popt_3, pcov_3 = curve_fit(func, xlow_3, avgsLow_3, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate : bottom
xlow_3[85] += 10
for i in range(86,len(xlow_3)):
    xlow_3[i] += 10
##
##
##
## 
######################################
##### Spot 4 #########################
data_path_4 = data_paths[3]
test_4, low_4, high_4, xlow_4, xhigh_4 = getDataCounts_General(data_path=data_path_4, start='low')
avgsLow_4, stdLow_4 = getAverageCountsLow(low_4)
avgsHigh_4, stdHigh_4 = getAverageCountsLow(high_4)
popt_4, pcov_4 = curve_fit(func, xlow_4, avgsLow_4, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate: LEFT, NONE
##
##
##
## 
######################################
##### Spot 5 #########################
data_path_5 = data_paths[4]
test_5, low_5, high_5, xlow_5, xhigh_5 = getDataCounts_General(data_path=data_path_5, start='low')
avgsLow_5, stdLow_5 = getAverageCountsLow(low_5)
avgsHigh_5, stdHigh_5 = getAverageCountsLow(high_5)
popt_5, pcov_5 = curve_fit(func, xlow_5, avgsLow_5, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate: LEFT
xlow_5[5] += 10
xlow_5[21] += 10
xlow_5[22] += 10
xlow_5[67] += 10
for i in range(6,len(xlow_5)):
    xlow_5[i] += 10
for i in range(22,len(xlow_5)):
    xlow_5[i] += 10
for i in range(23,len(xlow_5)):
    xlow_5[i] += 10
for i in range(68,len(xlow_5)):
    xlow_5[i] += 10
##
##
##
## 
######################################
##### Spot 6 #########################
data_path_6 = data_paths[5]
test_6, low_6, high_6, xlow_6, xhigh_6 = getDataCounts_General(data_path=data_path_6, start='low')
avgsLow_6, stdLow_6 = getAverageCountsLow(low_6)
avgsHigh_6, stdHigh_6 = getAverageCountsLow(high_6)
popt_6, pcov_6 = curve_fit(func, xlow_6, avgsLow_6, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate:  LEFT, NONE

##
##
##
## 
######################################
##### Spot 7 #########################
data_path_7 = data_paths[6]
test_7, low_7, high_7, xlow_7, xhigh_7 = getDataCounts_General(data_path=data_path_7, start='low')
avgsLow_7, stdLow_7 = getAverageCountsLow(low_7)
avgsHigh_7, stdHigh_7 = getAverageCountsLow(high_7)
popt_7, pcov_7 = curve_fit(func, xlow_7, avgsLow_7, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate : NONE, BOTTOM 
##
##
##
## 
######################################
##### Spot 8 #########################
data_path_8 = data_paths[7]
test_8, low_8, high_8, xlow_8, xhigh_8 = getDataCounts_General(data_path=data_path_8, start='low')
avgsLow_8, stdLow_8 = getAverageCountsLow(low_8)
avgsHigh_8, stdHigh_8 = getAverageCountsLow(high_8)
popt_8, pcov_8 = curve_fit(func, xlow_8, avgsLow_8, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate : TOP
xlow_8[95] += 10
xlow_8[107] += 10
for i in range(96,len(xlow_8)):
    xlow_8[i] += 10
for i in range(108,len(xlow_8)):
    xlow_8[i] += 10
##
##
##
## 
######################################
##### Spot 9 #########################
data_path_9 = data_paths[8]
test_9, low_9, high_9, xlow_9, xhigh_9 = getDataCounts_General(data_path=data_path_9, start='low')
avgsLow_9, stdLow_9 = getAverageCountsLow(low_9)
avgsHigh_9, stdHigh_9 = getAverageCountsLow(high_9)
popt_9, pcov_9 = curve_fit(func, xlow_9, avgsLow_9, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate : NONE, BOTTOM (rim)  
##
##
##
## 
######################################
##### Spot 10 #########################
data_path_10 = data_paths[9]
test_10, low_10, high_10, xlow_10, xhigh_10 = getDataCounts_General(data_path=data_path_10, start='low')
avgsLow_10, stdLow_10 = getAverageCountsLow(low_10)
avgsHigh_10, stdHigh_10 = getAverageCountsLow(high_10)
popt_10, pcov_10 = curve_fit(func, xlow_10, avgsLow_10, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate : bottom, rim  

##
##
## 
######################################
##### Spot 11 #########################
data_path_11 = data_paths[10]
test_11, low_11, high_11, xlow_11, xhigh_11 = getDataCounts_General(data_path=data_path_11, start='low')
avgsLow_11, stdLow_11 = getAverageCountsLow(low_11)
avgsHigh_11, stdHigh_11 = getAverageCountsLow(high_11)
popt_11, pcov_11 = curve_fit(func, xlow_11, avgsLow_11, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate: TOP 
xlow_11[100] += 10
xhigh_11[90] += 10
for i in range(101,len(xlow_11)):
    xlow_11[i] += 10
for i in range(91,len(xhigh_11)):
    xhigh_11[i] += 10
##
##
##
## 
######################################
##### Spot 12 #########################
data_path_12 = data_paths[11]
test_12, low_12, high_12, xlow_12, xhigh_12 = getDataCounts_General(data_path=data_path_12, start='low')
avgsLow_12, stdLow_12 = getAverageCountsLow(low_12)
avgsHigh_12, stdHigh_12 = getAverageCountsLow(high_12)
popt_12, pcov_12 = curve_fit(func, xlow_12, avgsLow_12, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate: TOP , NONE
##
##
##
## 
######################################
##### Spot 13 #########################
data_path_13 = data_paths[12]
test_13, low_13, high_13, xlow_13, xhigh_13 = getDataCounts_General(data_path=data_path_13, start='low')
avgsLow_13, stdLow_13 = getAverageCountsLow(low_13)
avgsHigh_13, stdHigh_13 = getAverageCountsLow(high_13)
popt_13, pcov_13 = curve_fit(func, xlow_13, avgsLow_13, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate: NONE 
##
##
##
## 
######################################
##### Spot 14 #########################
data_path_14 = data_paths[13]
test_14, low_14, high_14, xlow_14, xhigh_14 = getDataCounts_General(data_path=data_path_14, start='low')
avgsLow_14, stdLow_14 = getAverageCountsLow(low_14)
avgsHigh_14, stdHigh_14 = getAverageCountsLow(high_14)
popt_14, pcov_14 = curve_fit(func, xlow_14, avgsLow_14, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate:  NONE, RIGHT 
##
##
##
## 
######################################
##### Spot 15 #########################
data_path_15 = data_paths[14]
test_15, low_15, high_15, xlow_15, xhigh_15 = getDataCounts_General(data_path=data_path_15, start='low')
avgsLow_15, stdLow_15 = getAverageCountsLow(low_15)
avgsHigh_15, stdHigh_15 = getAverageCountsLow(high_15)
popt_15, pcov_15 = curve_fit(func, xlow_15, avgsLow_15, p0=[100, 0.01, 100], maxfev=10000)
## Mitigate: TOP, NONE 
##
##
##
## 
######################################
##### Spot 16 #########################
data_path_16 = data_paths[15]
test_16, low_16, high_16, xlow_16, xhigh_16 = getDataCounts_General(data_path=data_path_16, start='low')
avgsLow_16, stdLow_16 = getAverageCountsLow(low_16)
avgsHigh_16, stdHigh_16 = getAverageCountsLow(high_16)
popt_16, pcov_16 = curve_fit(func, xlow_16, avgsLow_16, p0=[100, 0.01, 100], maxfev=10000)
## Mitigate: RIGHT, NONE 
##
##
##
## 
######################################
##### Spot 17 #########################
data_path_17 = data_paths[16]
test_17, low_17, high_17, xlow_17, xhigh_17 = getDataCounts_General(data_path=data_path_17, start='low')
avgsLow_17, stdLow_17 = getAverageCountsLow(low_17)
avgsHigh_17, stdHigh_17 = getAverageCountsLow(high_17)
popt_17, pcov_17 = curve_fit(func, xlow_17, avgsLow_17, p0=[100, 0.01, 100], maxfev=10000)
## Mitigate: TOP, NONE 
##
##
##
## 
######################################
##### Spot 18 #########################
data_path_18 = data_paths[17]
test_18, low_18, high_18, xlow_18, xhigh_18 = getDataCounts_General(data_path=data_path_18, start='low')
avgsLow_18, stdLow_18 = getAverageCountsLow(low_18)
avgsHigh_18, stdHigh_18 = getAverageCountsLow(high_18)
popt_18, pcov_18 = curve_fit(func, xlow_18, avgsLow_18, p0=[100, 0.01, 100], maxfev=10000)
## Mitigate: RIGHT, NONE 
##
##
##
## 
######################################
##### Spot 19 #########################
data_path_19 = data_paths[18]
test_19, low_19, high_19, xlow_19, xhigh_19 = getDataCounts_General(data_path=data_path_19, start='low')
avgsLow_19, stdLow_19 = getAverageCountsLow(low_19)
avgsHigh_19, stdHigh_19 = getAverageCountsLow(high_19)
popt_19, pcov_19 = curve_fit(func, xlow_19, avgsLow_19, p0=[100, 0.01, 100], maxfev=10000)
## Mitigate: TOP RIM , NONE 
##
##
##
## 
######################################

## PLOT ##
# making subplots
fig, ax = plt.subplots(4, 5, constrained_layout = True, figsize=(20, 20), sharex=True, sharey=True) ## , sharex=True

## Plot low counts
lns1 = ax[0, 0].errorbar(xlow_1, avgsLow_1, yerr=stdLow_1, zorder=0, alpha=0.7, capsize=2, fmt='.', label='Neutral Fluor.', color=low_color)
ax[0, 1].errorbar(xlow_2, avgsLow_2, yerr=stdLow_2, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[0, 2].errorbar(xlow_3, avgsLow_3, yerr=stdLow_3, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[0, 3].errorbar(xlow_4, avgsLow_4, yerr=stdLow_4, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[0, 4].errorbar(xlow_5, avgsLow_5, yerr=stdLow_5, zorder=0, alpha=0.7, capsize=2, fmt='.', label='Neutral Fluor.', color=low_color)
ax[1, 0].errorbar(xlow_6, avgsLow_6, yerr=stdLow_6, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[1, 1].errorbar(xlow_7, avgsLow_7, yerr=stdLow_7, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[1, 2].errorbar(xlow_8, avgsLow_8, yerr=stdLow_8, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[1, 3].errorbar(xlow_9, avgsLow_9, yerr=stdLow_9, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[1, 4].errorbar(xlow_10, avgsLow_10, yerr=stdLow_10, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[2, 0].errorbar(xlow_11, avgsLow_11, yerr=stdLow_11, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[2, 1].errorbar(xlow_12, avgsLow_12, yerr=stdLow_12, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[2, 2].errorbar(xlow_13, avgsLow_13, yerr=stdLow_13, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[2, 3].errorbar(xlow_14, avgsLow_14, yerr=stdLow_14, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
# ax[2, 4].errorbar(xlow_15, avgsLow_15, yerr=stdLow_15, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)

ax[2, 4].errorbar(xlow_16, avgsLow_16, yerr=stdLow_16, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[3, 0].errorbar(xlow_17, avgsLow_17, yerr=stdLow_17, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[3, 1].errorbar(xlow_18, avgsLow_18, yerr=stdLow_18, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
ax[3, 2].errorbar(xlow_19, avgsLow_19, yerr=stdLow_19, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)
# ax[3, 4].errorbar(xlow_15, avgsLow_15, yerr=stdLow_15, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)


## Plot high counts 
ax1 = ax[0,0].twinx()
ax1.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax1.tick_params(axis='both', which='major', labelsize=20)
lns2 = ax1.errorbar(xhigh_1, avgsHigh_1, yerr=stdHigh_1, fmt='.', zorder=0, alpha=0.7, capsize=2, label='Conditioning Fluor.', color=high_color)

ax2 = ax[0,1].twinx()
ax2.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.errorbar(xhigh_2, avgsHigh_2, yerr=stdHigh_2, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax3 = ax[0,2].twinx()
ax3.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax3.tick_params(axis='both', which='major', labelsize=20)
ax3.errorbar(xhigh_3, avgsHigh_3, yerr=stdHigh_3, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax4 = ax[0,3].twinx()
ax4.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax4.tick_params(axis='both', which='major', labelsize=20)
ax4.errorbar(xhigh_4, avgsHigh_4, yerr=stdHigh_4, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax5 = ax[0,4].twinx()
ax5.tick_params(axis='y', which='major', labelsize=15)
ax5.errorbar(xhigh_5, avgsHigh_5, yerr=stdHigh_5, fmt='.', zorder=0, alpha=0.7, capsize=2, label='Conditioning Fluor.', color=high_color)

ax6 = ax[1,0].twinx()
ax6.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax6.tick_params(axis='both', which='major', labelsize=20)
ax6.errorbar(xhigh_6, avgsHigh_6, yerr=stdHigh_6, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax7 = ax[1,1].twinx()
ax7.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax7.tick_params(axis='both', which='major', labelsize=20)
ax7.errorbar(xhigh_7, avgsHigh_7, yerr=stdHigh_7, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax8 = ax[1,2].twinx()
ax8.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax8.tick_params(axis='both', which='major', labelsize=20)
ax8.errorbar(xhigh_8, avgsHigh_8, yerr=stdHigh_8, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax9 = ax[1,3].twinx()
ax9.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax9.tick_params(axis='both', which='major', labelsize=20)
ax9.errorbar(xhigh_9, avgsHigh_9, yerr=stdHigh_9, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax10 = ax[1,4].twinx()
ax10.tick_params(axis='y', which='major', labelsize=15)
ax10.errorbar(xhigh_10, avgsHigh_10, yerr=stdHigh_10, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax11 = ax[2,0].twinx()
ax11.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax11.tick_params(axis='both', which='major', labelsize=20)
ax11.errorbar(xhigh_11, avgsHigh_11, yerr=stdHigh_11, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax12 = ax[2,1].twinx()
ax12.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax12.tick_params(axis='both', which='major', labelsize=20)
ax12.errorbar(xhigh_12, avgsHigh_12, yerr=stdHigh_12, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax13 = ax[2,2].twinx()
ax13.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
# ax13.tick_params(axis='both', which='major', labelsize=20)
ax13.errorbar(xhigh_13, avgsHigh_13, yerr=stdHigh_13, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax14 = ax[2,3].twinx()
ax14.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, bottom=False, top=False, left=False, right=False) #bottom=False, top=False, left=False, right=False,
ax14.errorbar(xhigh_14, avgsHigh_14, yerr=stdHigh_14, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

# ax15 = ax[2,4].twinx()
# ax15.tick_params(axis='y', which='major', labelsize=15)
# ax15.errorbar(xhigh_15, avgsHigh_15, yerr=stdHigh_15, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax16 = ax[2,4].twinx()
ax16.tick_params(axis='y', which='major', labelsize=15)
# ax16.tick_params(axis='y', which='major',  
# labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
ax16.errorbar(xhigh_16, avgsHigh_16, yerr=stdHigh_16, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax17 = ax[3,0].twinx()
ax17.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
ax17.errorbar(xhigh_17, avgsHigh_17, yerr=stdHigh_17, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax18 = ax[3,1].twinx()
ax18.tick_params(axis='y', which='major',  
labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
ax18.errorbar(xhigh_18, avgsHigh_18, yerr=stdHigh_18, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

ax19 = ax[3,2].twinx()
# ax19.tick_params(axis='y', which='major',  
# labelbottom=False, labeltop=False, labelleft=False, labelright=False, right=False)
ax19.tick_params(axis='y', which='major', labelsize=15)
ax19.errorbar(xhigh_19, avgsHigh_19, yerr=stdHigh_19, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

# ax20 = ax[3,4].twinx()
# ax20.tick_params(axis='y', which='major', labelsize=15)

## Plot fits 
lns3 = ax[0, 0].plot(xlow_1, func(xlow_1, *popt_1), 'k-', linewidth=2.5, zorder=2, label='Fit')
ax[0, 1].plot(xlow_2, func(xlow_2, *popt_2), 'k-', linewidth=2.5, zorder=2)
ax[0, 2].plot(xlow_3, func(xlow_3, *popt_3), 'k-', linewidth=2.5, zorder=2)
ax[0, 3].plot(xlow_4, func(xlow_4, *popt_4), 'k-', linewidth=2.5, zorder=2)
ax[0, 4].plot(xlow_5, func(xlow_5, *popt_5), 'k-', linewidth=2.5, zorder=2, label='Fit')
ax[1, 0].plot(xlow_6, func(xlow_6, *popt_6), 'k-', linewidth=2.5, zorder=2)
ax[1, 1].plot(xlow_7, func(xlow_7, *popt_7), 'k-', linewidth=2.5, zorder=2)
ax[1, 2].plot(xlow_8, func(xlow_8, *popt_8), 'k-', linewidth=2.5, zorder=2)
ax[1, 3].plot(xlow_9, func(xlow_9, *popt_9), 'k-', linewidth=2.5, zorder=2)
ax[1, 4].plot(xlow_10, func(xlow_10, *popt_10), 'k-', linewidth=2.5, zorder=2)
ax[2, 0].plot(xlow_11, func(xlow_11, *popt_11), 'k-', linewidth=2.5, zorder=2)
ax[2, 1].plot(xlow_12, func(xlow_12, *popt_12), 'k-', linewidth=2.5, zorder=2)
ax[2, 2].plot(xlow_13, func(xlow_13, *popt_13), 'k-', linewidth=2.5, zorder=2)
ax[2, 3].plot(xlow_14, func(xlow_14, *popt_14), 'k-', linewidth=2.5, zorder=2)
# ax[2, 4].plot(xlow_15, func(xlow_15, *popt_15), 'k-', linewidth=2.5, zorder=2)

ax[2, 4].plot(xlow_16, func(xlow_16, *popt_16), 'k-', linewidth=2.5, zorder=2)
ax[3, 0].plot(xlow_17, func(xlow_17, *popt_17), 'k-', linewidth=2.5, zorder=2)
ax[3, 1].plot(xlow_18, func(xlow_18, *popt_18), 'k-', linewidth=2.5, zorder=2)
ax[3, 2].plot(xlow_19 , func(xlow_19, *popt_19), 'k-', linewidth=2.5, zorder=2)


# ax[0, 0].set_ylim(0, 300)
# ax[0, 1].set_ylim(0, 350)
# ax[0, 2].set_ylim(0, 60)
# ax[0, 3].set_ylim(0, 250)
# ax[0, 4].set_ylim(0, 55)
# ax[1, 0].set_ylim(0, 220)
# ax[1, 1].set_ylim(0, 210)
# ax[1, 2].set_ylim(0, 190)
# ax[1, 3].set_ylim(0, 13)
# ax[1, 4].set_ylim(0, 2.1)
# ax[2, 0].set_ylim(0, 160)
# ax[2, 1].set_ylim(0, 65)
# ax[2, 2].set_ylim(0, 210)
# ax[2, 3].set_ylim(0, 10)
# ax[2, 4].set_ylim(0, 250)

ax[0, 0].set_ylim(0, 300)
ax[0, 1].set_ylim(0, 300)
ax[0, 2].set_ylim(0, 300)
ax[0, 3].set_ylim(0, 300)
ax[0, 4].set_ylim(0, 300)
ax[1, 0].set_ylim(0, 300)
ax[1, 1].set_ylim(0, 300)
ax[1, 2].set_ylim(0, 300)
ax[1, 3].set_ylim(0, 300)
ax[1, 4].set_ylim(0, 300)
ax[2, 0].set_ylim(0, 300)
ax[2, 1].set_ylim(0, 300)
ax[2, 2].set_ylim(0, 300)
ax[2, 3].set_ylim(0, 300)
ax[2, 4].set_ylim(0, 300)

ax[3, 0].set_ylim(0, 300)
ax[3, 1].set_ylim(0, 300)
ax[3, 2].set_ylim(0, 300)
# ax[3, 3].set_ylim(0, 300)


# ax1.set_ylim(0, 650)
# ax2.set_ylim(0, 650)
# ax3.set_ylim(0, 650)
# ax4.set_ylim(0, 650)
# ax5.set_ylim(0, 650)
# ax6.set_ylim(0, 650)
# ax7.set_ylim(0, 650)
# ax8.set_ylim(0, 650)
# ax9.set_ylim(0, 650)
# ax10.set_ylim(0, 650)
# ax11.set_ylim(0, 650)
# ax12.set_ylim(0, 650)
# ax13.set_ylim(0, 650)
# ax14.set_ylim(0, 15)
# ax15.set_ylim(0, 650)

ax1.set_ylim(0, 500)
ax2.set_ylim(0, 500)
ax3.set_ylim(0, 500)
ax4.set_ylim(0, 500)
ax5.set_ylim(0, 500)
ax6.set_ylim(0, 500)
ax7.set_ylim(0, 500)
ax8.set_ylim(0, 500)
ax9.set_ylim(0, 500)
ax10.set_ylim(0, 500)
ax11.set_ylim(0, 500)
ax12.set_ylim(0, 500)
ax13.set_ylim(0, 500)
ax14.set_ylim(0, 500)
# ax15.set_ylim(0, 500)

ax16.set_ylim(0, 500)
ax17.set_ylim(0, 500)
ax18.set_ylim(0, 500)
ax19.set_ylim(0, 500)


# ax[0, 0].set_xlim(0, 510)
# ax[0, 1].set_xlim(0, 720)
# ax[0, 2].set_xlim(0, 700)
# ax[0, 3].set_xlim(0, 1030)
# ax[0, 4].set_xlim(0, 1800)
# ax[1, 0].set_xlim(0, 1200)
# ax[1, 1].set_xlim(0, 1200)
# ax[1, 2].set_xlim(0, 1200)
# ax[1, 3].set_xlim(0, 1200)
# ax[1, 4].set_xlim(0, 1200)
# ax[2, 0].set_xlim(0, 1200)
# ax[2, 1].set_xlim(0, 1200)
# ax[2, 2].set_xlim(0, 1200)
# ax[2, 3].set_xlim(0, 1200)
# ax[2, 4].set_xlim(0, 550)

# ax[0, 0].set_xlim(0, 1800)
# ax[0, 1].set_xlim(0, 1800)
# ax[0, 2].set_xlim(0, 1800)
# ax[0, 3].set_xlim(0, 1800)
# ax[0, 4].set_xlim(0, 1800)
# ax[1, 0].set_xlim(0, 1800)
# ax[1, 1].set_xlim(0, 1800)
# ax[1, 2].set_xlim(0, 1800)
# ax[1, 3].set_xlim(0, 1800)
# ax[1, 4].set_xlim(0, 1800)
# ax[2, 0].set_xlim(0, 1800)
# ax[2, 1].set_xlim(0, 1800)
# ax[2, 2].set_xlim(0, 1800)
# ax[2, 3].set_xlim(0, 1800)
# ax[2, 4].set_xlim(0, 1800)

# ax[0, 0].set_xlim(0, 1200)
# ax[0, 1].set_xlim(0, 1200)
# ax[0, 2].set_xlim(0, 1200)
# ax[0, 3].set_xlim(0, 1200)
# ax[0, 4].set_xlim(0, 1200)
# ax[1, 0].set_xlim(0, 1200)
# ax[1, 1].set_xlim(0, 1200)
# ax[1, 2].set_xlim(0, 1200)
# ax[1, 3].set_xlim(0, 1200)
# ax[1, 4].set_xlim(0, 1200)
# ax[2, 0].set_xlim(0, 1200)
# ax[2, 1].set_xlim(0, 1200)
# ax[2, 2].set_xlim(0, 1200)
# ax[2, 3].set_xlim(0, 1200)
# ax[2, 4].set_xlim(0, 1200)

ax[0, 0].set_xlim(0, 1000)
ax[0, 1].set_xlim(0, 1000)
ax[0, 2].set_xlim(0, 1000)
ax[0, 3].set_xlim(0, 1000)
ax[0, 4].set_xlim(0, 1000)
ax[1, 0].set_xlim(0, 1000)
ax[1, 1].set_xlim(0, 1000)
ax[1, 2].set_xlim(0, 1000)
ax[1, 3].set_xlim(0, 1000)
ax[1, 4].set_xlim(0, 1000)
ax[2, 0].set_xlim(0, 1000)
ax[2, 1].set_xlim(0, 1000)
ax[2, 2].set_xlim(0, 1000)
ax[2, 3].set_xlim(0, 1000)
ax[2, 4].set_xlim(0, 1000)

ax[3, 0].set_xlim(0, 1000)
ax[3, 1].set_xlim(0, 1000)
ax[3, 2].set_xlim(0, 1000)
# ax[3, 3].set_xlim(0, 1000)


ax[0,0].text(30, 250, 'A', size=15) #previously: x=1500
ax[0,1].text(30, 250, 'B', size=15) 
ax[0,2].text(30, 250, 'C', size=15)
ax[0,3].text(30, 250, 'D', size=15)
ax[0,4].text(30, 250, 'E', size=15)
ax[1,0].text(30, 250, 'F', size=15)
ax[1,1].text(30, 250, 'G', size=15)
ax[1,2].text(30, 250, 'H', size=15)
ax[1,3].text(30, 250, 'I', size=15)
ax[1,4].text(30, 250, 'J', size=15)
ax[2,0].text(30, 250, 'K', size=15)
ax[2,1].text(30, 250, 'L', size=15)
ax[2,2].text(30, 250, 'M', size=15)
ax[2,3].text(30, 250, 'N', size=15)
ax[2,4].text(30, 250, 'O', size=15)

ax[3,0].text(30, 250, 'P', size=15)
ax[3,1].text(30, 250, 'Q', size=15)
ax[3,2].text(30, 250, 'R', size=15)
# ax[3,3].text(30, 250, 'S', size=15)



# for i in range(3):
#     for j in range(5):
#         ax[i,j].tick_params(axis='both', which='major', labelsize=20)

ax[0,0].tick_params(axis='y', which='major', labelsize=15)
ax[1,0].tick_params(axis='y', which='major', labelsize=15)
ax[2,0].tick_params(axis='y', which='major', labelsize=15)
ax[3,0].tick_params(axis='y', which='major', labelsize=15)

ax[0,1].tick_params(axis='y', left=False)
ax[0,2].tick_params(axis='y', left=False)
ax[0,3].tick_params(axis='y', left=False)
ax[0,4].tick_params(axis='y', left=False)
ax[1,1].tick_params(axis='y', left=False)
ax[1,2].tick_params(axis='y', left=False)
ax[1,3].tick_params(axis='y', left=False)
ax[1,4].tick_params(axis='y', left=False)
ax[2,1].tick_params(axis='y', left=False)
ax[2,2].tick_params(axis='y', left=False)
ax[2,3].tick_params(axis='y', left=False)
ax[2,4].tick_params(axis='y', left=False)

ax[3,0].tick_params(axis='y', left=False)
ax[3,1].tick_params(axis='y', left=False)
ax[3,2].tick_params(axis='y', left=False)
# ax[3,3].tick_params(axis='y', left=False)
# ax[3,4].tick_params(axis='y', left=False)

ax[3,0].tick_params(axis='x', which='major', labelsize=15)
ax[3,1].tick_params(axis='x', which='major', labelsize=15)
ax[3,2].tick_params(axis='x', which='major', labelsize=15)
ax[2,3].tick_params(axis='x', which='major', labelsize=15, labelbottom=True)
ax[2,4].tick_params(axis='x', which='major', labelsize=15, labelbottom=True)

# ax[3,3].tick_params(axis='x', which='major', labelsize=15)
# ax[3,4].tick_params(axis='x', which='major', labelsize=15)

fig.text(0.5, 0.01, 'Pulse Number', ha='center', fontsize=23)
fig.text(0.005, 0.5, 'Neutral fluorescence / arb. units', va='center', rotation='vertical', fontsize=20)
fig.text(0.97, 0.5, 'Conditioning neutral fluorescence / arb. units', va='center', rotation=-90, fontsize=20)

lines, labels = ax[0,0].get_legend_handles_labels()
lines2, labels2 = ax5.get_legend_handles_labels()
ax[0,0].legend(lines + lines2, labels + labels2, loc='upper center', ncol=3, bbox_to_anchor=(2.9, 1.5), fancybox=True, fontsize=15) #previously: bbox_to_anchor=(4.5, 1.5)

# ax[0, 0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.5), fancybox=True, fontsize=15)

plt.subplots_adjust(left=0.08, bottom=0.1, right=0.905, top=0.87, wspace=0.2, hspace=0.2) 
#previously v1: wspace=1.0, hspace=0.3, left=0.07, bottom=0.1, right=0.905, top=0.895,

## bottom=0.15, 
                    # right=0.95, 
                    # top=0.9, 
                    # wspace=0.6, 
                    # hspace=0.09

fig.delaxes(ax[3,4])
fig.delaxes(ax[3,3])

## Inset target map from inkscape
im = plt.imread('grid_targetMap_edgecase.png') # insert local path of the image, old version: 
newax = fig.add_axes([0.62,0.09,0.155,0.155], anchor='W') 
newax.imshow(im)
newax.axis('off')

plt.show()

# fig.savefig('Conditioning-vs-Pulse_Grid_v1.pdf', dpi=300, bbox_inches='tight', format='pdf')

################### TESTING w/ a test file/dataset #######################
############ BASIC DATA IMPORT TEST #########################################
# data_path_test = r'test.txt'
# test, low, high, xlow, xhigh = getDataCounts_General(data_path=data_path_test, start='low')
# avgs_test, std_test = getAverageCountsLow(test)
#############################################################################.errorbar(xlow_1, avgsLow_1, yerr=stdLow_1, zorder=0, alpha=0.7, capsize=2)