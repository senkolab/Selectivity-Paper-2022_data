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
    r'Neut-Fluor-Spect_63uJ_145uJ__1.txt', ## DONE
    r'spot2_Neut-Fluor-Spect_63uJ_145uJ__1.txt', ## DONE
    r'spot3_Neut-Fluor-Spect_63uJ_145uJ__1.txt', ## DONE
    r'spot4_Neut-Fluor-Spect_63uJ_145uJ__1.txt', ## DONE
    r'spot5_Neut-Fluor-Spect_63uJ_145uJ_7o19x_6o05y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o19x_6o02y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o20x_6o04y.txt', ## DONE 
    r'Neut-Fluor-Spect_63uJ_145uJ_7o17x_6o00y.txt', ## DONE 
    r'Neut-Fluor-Spect_63uJ_145uJ_7o27x_6o10y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o15x_5o90y.txt', ## DONE 
    r'Neut-Fluor-Spect_63uJ_145uJ_7o15x_5o96y.txt', ## DONE 
    r'Neut-Fluor-Spect_63uJ_145uJ_7o24x_6o08y.txt', ## DONE
    r'Neut-Fluor-Spect_63uJ_145uJ_7o22x_6o02y.txt', ## DONE
    r'redDownloaded_Neut-Fluor-Spect_63uJ_145uJ_7o13x_5o90y.txt', ## DONE ***accurate?? 
    r'Neut-Fluor-Spect_63uJ_145uJ_7o26x_6o10y.txt' ## DONE, **mirror 
] 
##### Spot 1 #########################
data_path_1 = data_paths[0]
test_1, low_1, high_1, xlow_1, xhigh_1 = getDataCounts_General(data_path=data_path_1, start='low')
avgsLow_1, stdLow_1 = getAverageCountsLow(low_1)
avgsHigh_1, stdHigh_1 = getAverageCountsLow(high_1)
popt_1, pcov_1 = curve_fit(func, xlow_1, avgsLow_1, p0=[120, 0.05, 120], maxfev=10000)
## Mitigate 
xlow_1[1] = 30 
xlow_1[2] = 60
xhigh_1[0] = 30 
xhigh_1[1] = 50 
for i in range(3,len(xlow_1)):
    xlow_1[i] += 30
for i in range(2, len(xhigh_1)):
    xhigh_1[i] += 30 
## 
##
## 
######################################
##### Spot 2 #########################
data_path_2 = data_paths[1]
test_2, low_2, high_2, xlow_2, xhigh_2 = getDataCounts_General(data_path=data_path_2, start='low')
avgsLow_2, stdLow_2 = getAverageCountsLow(low_2)
avgsHigh_2, stdHigh_2 = getAverageCountsLow(high_2)
popt_2, pcov_2 = curve_fit(func, xlow_2, avgsLow_2, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate 
xlow_2[2] = 40
xlow_2[20] += 10
xlow_2[67] += 10
xhigh_2[14] += 10 
for i in range(3,len(xlow_2)):
    xlow_2[i] += 10
for i in range(21, len(xlow_2)):
    xlow_2[i] += 10
for i in range(68, len(xlow_2)):
    xlow_2[i] += 10
for i in range(15, len(xhigh_2)):
    xhigh_2[i] += 10
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
## Mitigate : NONE
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
## Mitigate 
xlow_4[30] += 10
xhigh_4[52] += 10
xhigh_4[72] += 10
for i in range(31,len(xlow_4)):
    xlow_4[i] += 10
for i in range(53, len(xhigh_4)):
    xhigh_4[i] += 10
for i in range(73, len(xhigh_4)):
    xhigh_4[i] += 10
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
## Mitigate : NONE
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
## Mitigate 
xlow_6[93] += 10
xlow_6[96] += 10
for i in range(94,len(xlow_6)):
    xlow_6[i] += 10
for i in range(97,len(xlow_6)):
    xlow_6[i] += 10
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
## Mitigate : NONE 
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
## Mitigate : NONE 
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
## Mitigate : NONE 
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
## Mitigate 
xlow_10[103] += 10
xlow_10[106] += 10
for i in range(104,len(xlow_10)):
    xlow_10[i] += 10
for i in range(107,len(xlow_10)):
    xlow_10[i] += 10
##
##
## 
######################################
##### Spot 11 #########################
data_path_11 = data_paths[10]
test_11, low_11, high_11, xlow_11, xhigh_11 = getDataCounts_General(data_path=data_path_11, start='high')
avgsLow_11, stdLow_11 = getAverageCountsLow(low_11)
avgsHigh_11, stdHigh_11 = getAverageCountsLow(high_11)
popt_11, pcov_11 = curve_fit(func, xlow_11, avgsLow_11, p0=[20, 0.05, 20], maxfev=10000)
## Mitigate: NONE 
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
## Mitigate: NONE 
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
## Mitigate
xlow_13[104] += 10
for i in range(105,len(xlow_13)):
    xlow_13[i] += 10
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
## Mitigate: NONE 
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
## Mitigate: 
for i in range(len(xlow_15)):
    xlow_15[i] += 58
for i in range(len(xhigh_15)):
    xhigh_15[i] += 58
##
##
##
## 
######################################

## PLOT ##
# making subplots
fig, ax = plt.subplots(3, 5, constrained_layout = True, figsize=(20, 20), sharex=True, sharey=True) ## , sharex=True

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
ax[2, 4].errorbar(xlow_15, avgsLow_15, yerr=stdLow_15, zorder=0, alpha=0.7, capsize=2, fmt='.', color=low_color)

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

ax15 = ax[2,4].twinx()
ax15.tick_params(axis='y', which='major', labelsize=15)
ax15.errorbar(xhigh_15, avgsHigh_15, yerr=stdHigh_15, fmt='.', zorder=0, alpha=0.7, capsize=2, color=high_color)

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
ax[2, 4].plot(xlow_15, func(xlow_15, *popt_15), 'k-', linewidth=2.5, zorder=2)

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

ax[0, 0].set_ylim(0, 350)
ax[0, 1].set_ylim(0, 350)
ax[0, 2].set_ylim(0, 350)
ax[0, 3].set_ylim(0, 350)
ax[0, 4].set_ylim(0, 350)
ax[1, 0].set_ylim(0, 350)
ax[1, 1].set_ylim(0, 350)
ax[1, 2].set_ylim(0, 350)
ax[1, 3].set_ylim(0, 350)
ax[1, 4].set_ylim(0, 350)
ax[2, 0].set_ylim(0, 350)
ax[2, 1].set_ylim(0, 350)
ax[2, 2].set_ylim(0, 350)
ax[2, 3].set_ylim(0, 350)
ax[2, 4].set_ylim(0, 350)

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

ax1.set_ylim(0, 650)
ax2.set_ylim(0, 650)
ax3.set_ylim(0, 650)
ax4.set_ylim(0, 650)
ax5.set_ylim(0, 650)
ax6.set_ylim(0, 650)
ax7.set_ylim(0, 650)
ax8.set_ylim(0, 650)
ax9.set_ylim(0, 650)
ax10.set_ylim(0, 650)
ax11.set_ylim(0, 650)
ax12.set_ylim(0, 650)
ax13.set_ylim(0, 650)
ax14.set_ylim(0, 650)
ax15.set_ylim(0, 650)



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

ax[0,0].text(850, 300, 'A', size=15) #previously: x=1500
ax[0,1].text(850, 300, 'B', size=15) 
ax[0,2].text(850, 300, 'C', size=15)
ax[0,3].text(850, 300, 'D', size=15)
ax[0,4].text(850, 300, 'E', size=15)
ax[1,0].text(850, 300, 'F', size=15)
ax[1,1].text(850, 300, 'G', size=15)
ax[1,2].text(850, 300, 'H', size=15)
ax[1,3].text(850, 300, 'I', size=15)
ax[1,4].text(850, 300, 'J', size=15)
ax[2,0].text(850, 300, 'K', size=15)
ax[2,1].text(850, 300, 'L', size=15)
ax[2,2].text(850, 300, 'M', size=15)
ax[2,3].text(850, 300, 'N', size=15)
ax[2,4].text(850, 300, 'O', size=15)


# for i in range(3):
#     for j in range(5):
#         ax[i,j].tick_params(axis='both', which='major', labelsize=20)

ax[0,0].tick_params(axis='y', which='major', labelsize=15)
ax[1,0].tick_params(axis='y', which='major', labelsize=15)
ax[2,0].tick_params(axis='y', which='major', labelsize=15)

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


ax[2,0].tick_params(axis='x', which='major', labelsize=15)
ax[2,1].tick_params(axis='x', which='major', labelsize=15)
ax[2,2].tick_params(axis='x', which='major', labelsize=15)
ax[2,3].tick_params(axis='x', which='major', labelsize=15)
ax[2,4].tick_params(axis='x', which='major', labelsize=15)

fig.text(0.5, 0.01, 'Pulse number', ha='center', fontsize=23)
fig.text(0.005, 0.5, 'Neutral fluorescence (arb. units)', va='center', rotation='vertical', fontsize=20)
fig.text(0.97, 0.5, 'Conditioning neutral fluorescence (arb. units)', va='center', rotation=-90, fontsize=20)

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

plt.show()

# fig.savefig('Conditioning-vs-Pulse_Grid_v1.pdf', dpi=300, bbox_inches='tight', format='pdf')
# fig.savefig('SM_conditioningCentreCases_Feb9Revision.pdf', dpi=300, bbox_inches='tight', format='pdf')

################### TESTING w/ a test file/dataset #######################
############ BASIC DATA IMPORT TEST #########################################
# data_path_test = r'test.txt'
# test, low, high, xlow, xhigh = getDataCounts_General(data_path=data_path_test, start='low')
# avgs_test, std_test = getAverageCountsLow(test)
#############################################################################.errorbar(xlow_1, avgsLow_1, yerr=stdLow_1, zorder=0, alpha=0.7, capsize=2)