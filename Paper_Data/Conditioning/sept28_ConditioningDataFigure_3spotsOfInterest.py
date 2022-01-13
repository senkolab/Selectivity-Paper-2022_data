import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import pandas as pd
import statistics
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
graph_edge_width = 0.5
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['axes.linewidth'] = graph_edge_width
mpl.rcParams.update({'errorbar.capsize': 2})
plt.style.use('seaborn-colorblind')

low_color = 'tab:blue'
high_color = 'tab:purple'
mod_color = 'tab:orange'
markersize = 2
markeredgewidth = 0.5
linewidth = 0.5
elinewidth = 1
legend_fontsize = 5
axis_fontsize = 10
ticks_fontsize = 6
alpha_marker = 0.7
graph_tick_width = graph_edge_width
fluence_label_pos = 0.75

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
golden_ratio = (5**.5 - 1) / 2
fig_width_in = columnwidth * inches_per_pt
height_ratio = golden_ratio
fig_height_in = fig_width_in * height_ratio

## Data Processing ## 

def getDataCounts_General(data_path, low_th, high_th, start): #th = threshold counts
    data_file = np.asarray(glob.glob(data_path))[0]
    data = np.asarray(pd.read_csv(data_file, delimiter = "[\]\[\t]", engine='python', header=None))
    
    counts = data[:,2]
    rows = counts.size

    floats_high = []
    floats_low = []
    floats = []

    x_avg_low = []
    x_avg_high = []

    for i in range(0,rows):
        floats.append([float(x) for x in counts[i].split()])

    pulseNum_low = 0
    pulseNum_high = 0 
    tmp = []
    j = 0

    if start == 'low': 
        marker = 'low'
    else: 
        marker = 'high'

    while (j < rows):
        low_repeats = []
        high_repeats = []
        if (j>50):
            if (floats[j][0] <= low_th): 
                for r in range(j+1, rows):
                    if (floats[r][0] > high_th):
                        marker = 'high'
                        break
                    else:
                        low_repeats.append(r)
                if (len(low_repeats) > 0):
                    for n in range(0, len(low_repeats)):
                        for k in range(10):
                            floats[j][k] = floats[j][k] + floats[low_repeats[n]][k]
                    for k in range(10):
                        floats[j][k] = floats[j][k] / (len(low_repeats)+1)
                    pulseNum_low += (len(low_repeats) + 1) * 10
                    x_avg_low.append(pulseNum_low)
                    tmp.append(floats[j])
                    # print(j)
                    # print(floats[j])
                    j += (len(low_repeats) + 1)
                    marker = 'high'
                else:
                    pulseNum_low += 10
                    tmp.append(floats[j])
                    x_avg_low.append(pulseNum_low) 
                    j += 1
                    marker = 'high'
            elif (floats[j][0] > high_th):
                for r in range(j+1, rows):
                    if (floats[r][0] <= low_th):
                        marker = 'low'
                        break
                    else:
                        high_repeats.append(r)
                if (len(high_repeats) > 0):
                    for n in range(0, len(high_repeats)):
                        for k in range(10):
                            floats[j][k] = floats[j][k] + floats[high_repeats[n]][k]
                    for k in range(10):
                        floats[j][k] = floats[j][k] / (len(high_repeats)+1)
                    pulseNum_high += (len(high_repeats) + 1) * 10
                    x_avg_high.append(pulseNum_high)
                    tmp.append(floats[j])
                    j += (len(high_repeats) + 1)
                    marker = 'low'
                else:
                    pulseNum_high += 10
                    tmp.append(floats[j]) 
                    x_avg_high.append(pulseNum_high) 
                    j += 1  
                    marker = 'low'
            else:
                tmp.append(floats[j])
                if marker == 'low':
                    pulseNum_low += 10
                    x_avg_low.append(pulseNum_low)
                    marker = 'high'
                else: 
                    pulseNum_high += 10
                    x_avg_high.append(pulseNum_high) 
                    marker = 'low'
                j += 1
        else:
            tmp.append(floats[j])
            if marker == 'low':
                pulseNum_low += 10
                x_avg_low.append(pulseNum_low)
                marker = 'high'
            else: 
                pulseNum_high += 10
                x_avg_high.append(pulseNum_high) 
                marker = 'low'
            j += 1

    # print(tmp)
    rows = len(tmp)
    for i in range(0,rows-1,2):
        if start == 'low':
            floats_high.append(tmp[i+1])
            floats_low.append(tmp[i])
        else: 
            # print(tmp[i])   
            floats_high.append(tmp[i])
            floats_low.append(tmp[i+1])      
    if (start == 'low' and rows%2 != 0): #check if odd
        floats_low.append(tmp[rows-1])
    if (start == 'high' and rows%2 != 0): #check if odd
        floats_high.append(tmp[rows-1])

    floats_low_array = np.array(floats_low)
    floats_high_array = np.array(floats_high)

    counts_low_data = floats_low_array.astype(int)
    counts_high_data = floats_high_array.astype(int)
    
    return counts_low_data, counts_high_data, np.array(x_avg_low), np.array(x_avg_high)

def getAverageCountsLow(low_counts):
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

    m = avgs.max()
    y_scale = avgs/m
    return avgs, sdev, y_scale

## Fits : Defining the exponential function form ##

def func(x, a, b, c):
    return -a * np.exp(-b * x) + c

## Spots ##
#data paths to all the raw data of the spots data imported 
numSpots = 3

data_path_0 = r'spot5_Neut-Fluor-Spect_63uJ_145uJ_7o19x_6o05y.txt'
data_path_1 = r'Neut-Fluor-Spect_63uJ_145uJ_7o27x_6o10y.txt'
data_path_2 = r'Neut-Fluor-Spect_63uJ_145uJ_7o15x_5o96y.txt'
data_path_2_b = r'Neut-Fluor-Spect_63uJ_145uJ_7o15x_5o96y_1st4Lines.txt'

low_counts_0, high_counts_0, x_low_counts_0, x_high_counts_0 = getDataCounts_General(data_path_0, low_th=60, high_th=60, start='low') 
low_counts_1, high_counts_1, x_low_counts_1, x_high_counts_1 = getDataCounts_General(data_path_1, low_th=12, high_th=12, start='low') 
low_counts_2, high_counts_2, x_low_counts_2, x_high_counts_2 = getDataCounts_General(data_path_2, low_th=160, high_th=160, start='high') 
low_counts_2_b, high_counts_2_b, x_low_counts_2_b, x_high_counts_2_b = getDataCounts_General(data_path_2_b, low_th=20, high_th=20, start='high') 

avgs_low_0, sdev_low_0, y_scale_0 = getAverageCountsLow(low_counts_0)
avgs_low_1, sdev_low_1, y_scale_1 = getAverageCountsLow(low_counts_1)
avgs_low_2, sdev_low_2, y_scale_2 = getAverageCountsLow(low_counts_2)
avgs_low_2_b, sdev_low_2_b, y_scale_2_b = getAverageCountsLow(low_counts_2_b)

avgs_high_0, sdev_high_0, y_scale_0 = getAverageCountsLow(high_counts_0)
avgs_high_1, sdev_high_1, y_scale_1 = getAverageCountsLow(high_counts_1)
avgs_high_2, sdev_high_2, y_scale_2 = getAverageCountsLow(high_counts_2)
avgs_high_2_b, sdev_high_2_b, y_scale_2_b = getAverageCountsLow(high_counts_2_b)

x_low_counts_2_combo = np.array(list(x_low_counts_2_b) + list(x_low_counts_2))    
avgs_low_2_combo = np.array(list(avgs_low_2_b) + list(avgs_low_2))
sdev_low_2_combo = np.array(list(sdev_low_2_b) + list(sdev_low_2))

x_high_counts_2_combo = np.array(list(x_high_counts_2_b) + list(x_high_counts_2))    
avgs_high_2_combo = np.array(list(avgs_high_2_b) + list(avgs_high_2))
sdev_high_2_combo = np.array(list(sdev_high_2_b) + list(sdev_high_2))

for i in range(2, len(x_high_counts_2_combo)):
    x_high_counts_2_combo[i] = x_high_counts_2_combo[i] + 20 
for i in range(2, len(x_low_counts_2_combo)):
    x_low_counts_2_combo[i] = x_low_counts_2_combo[i] + 20 

popt_0, pcov_0 = curve_fit(func, x_low_counts_0, avgs_low_0, p0=[20, 0.05, 20], maxfev=10000)
popt_1, pcov_1 = curve_fit(func, x_low_counts_1, avgs_low_1, p0=[2, 0.05, 2], maxfev=10000)
popt_2, pcov_2 = curve_fit(func, x_low_counts_2_combo, avgs_low_2_combo, p0=[90, 0.05, 90], maxfev=10000)

## PLOT ##
# making subplots
height_ratio = 0.95
fig_height_in = fig_width_in * height_ratio
fig, ax = plt.subplots(2, 1, constrained_layout = True, figsize=(fig_width_in, fig_height_in), dpi=300, sharex=True) #constrained_layout = True,

# fig.tight_layout()

# fig.text(0.02, 0.5, 'Neutral barium fluorescence / arb. units', va='center', rotation='vertical', fontsize=23)

## Plot low counts
ax[0].errorbar(x_low_counts_2_combo, avgs_low_2_combo, yerr=sdev_low_2_combo, \
               fmt='^', color=high_color, label="High yield", zorder=0, alpha=alpha_marker, ms=markersize, \
                   mew=markeredgewidth, elinewidth=elinewidth)
ax[0].errorbar(x_low_counts_0, avgs_low_0, yerr=sdev_low_0, fmt='s', color=mod_color, label="Moderate yield", \
               zorder=0, alpha=alpha_marker, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
ax[0].errorbar(x_low_counts_1, avgs_low_1, yerr=sdev_low_1, fmt='o', color=low_color, label="Low yield", \
               zorder=0, alpha=alpha_marker, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
## Plot fits 
ax[0].plot(x_low_counts_1, func(x_low_counts_1, *popt_1), 'k-', lw=linewidth, zorder=1)
ax[0].plot(x_low_counts_0, func(x_low_counts_0, *popt_0), 'k-', lw=linewidth, label='Fit', zorder=1)
ax[0].plot(x_low_counts_2_combo, func(x_low_counts_2_combo, *popt_2), 'k-', lw=linewidth, zorder=1)


ax[0].tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax[0].legend(fontsize=legend_fontsize, loc=2) ## loc=1, , loc=(1.04,0)
ax[0].xaxis.set_major_locator(MultipleLocator(200))
ax[0].xaxis.set_minor_locator(MultipleLocator(50))
ax[0].yaxis.set_major_locator(MultipleLocator(50))
ax[0].yaxis.set_minor_locator(MultipleLocator(25))
ax[0].set_ylim(0, 155)
# Cut off at ~1000 
ax[0].set_xlim(0, 1000)
ax[0].annotate('0.20 J/$\mathregular{cm^2}$', xy=(fluence_label_pos*1000,0.9*155), fontsize=legend_fontsize)


## Plot high E counts 
ax[1].errorbar(x_high_counts_2_combo, avgs_high_2_combo, yerr=sdev_high_2_combo, fmt='^', \
               color=high_color, label="High yield", \
               alpha=alpha_marker, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
ax[1].errorbar(x_high_counts_0, avgs_high_0, yerr=sdev_high_0, fmt='s', color=mod_color, label="Moderate yield", \
               alpha=alpha_marker, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
ax[1].errorbar(x_high_counts_1, avgs_high_1, yerr=sdev_high_1, fmt='o', color=low_color, label="Low yield", \
               alpha=alpha_marker, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth) #, mfc='None'

## Inset target map from inkscape
im = plt.imread('target_inset.png') # insert local path of the image, old version: 
newax = fig.add_axes([0.13,0.35,0.135,0.135], anchor='W') 
newax.imshow(im)
newax.axis('off')
    

ax[1].set_xlabel('Pulse number',fontsize=axis_fontsize)
ax[1].set_ylabel('Neutral fluorescence / arb. units',fontsize=axis_fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax[1].xaxis.set_major_locator(MultipleLocator(200))
ax[1].xaxis.set_minor_locator(MultipleLocator(50))
ax[1].yaxis.set_major_locator(MultipleLocator(200))
ax[1].yaxis.set_minor_locator(MultipleLocator(50))
ax[1].yaxis.set_label_coords(-.12, 1.0)
ax[1].set_ylim(0, 500)
# Cut off at ~1000 
ax[1].set_xlim(0, 1000)
ax[1].annotate('0.48 J/$\mathregular{cm^2}$', xy=(fluence_label_pos*1000,0.9*500), fontsize=legend_fontsize)
ax[1].text(425, -210, '(a)', ha='center', fontsize=8, font='Times New Roman')


plt.subplots_adjust(wspace=0.6, hspace=0.09)

plt.show()

############################## TESTING #################################################

## TESTING w/ a test file/dataset ##
# data_path_test = r'test.txt'
# low_test, high_test, x_low_test, x_high_test = getDataCounts_General(data_path_test, low_th=10, high_th=10, start='low') 
# avgs_low_test, std_low_test, y_low_scale_test = getAverageCountsLow(low_test)
# avgs_high_test, std_high_test, y_high_scale_test = getAverageCountsLow(high_test)
# popt_test, pcov_test = curve_fit(func, x_low_test, avgs_low_test, p0=[4, 0.05, 4], maxfev=10000)

# ## PLOT ##
# bbox = dict(boxstyle ="round", fc ="0.9", alpha=0.5)
# window_size = 10

# figure, ax1 = plt.subplots()

# ax1.plot(x_low_test, avgs_low_test, 'ro', label="LOW test")
# plt.plot(x_low_test, func(x_low_test, *popt_test), 'r-', linewidth=2.5)
# ax1.set_xlabel('Number of Pulses', fontsize = 7)
# ax1.set_ylabel('Average Counts', fontsize = 7)

# ax2 = ax1.twinx()
# ax2.plot(x_high_test, avgs_high_test, 'ro', alpha=0.2, label='HIGH test')
# ax2.set_ylabel('Conditioning Counts', fontsize = 7, rotation=-90)

# # lines_1, labels_1 = ax1.get_legend_handles_labels()
# # lines_2, labels_2 = ax2.get_legend_handles_labels()

# # lines = lines_1 + lines_2
# # labels = labels_1 + labels_2

# # ax1.legend(lines, labels, loc=0)

# ax1.legend(loc=(0.015,0.8))
# ax2.legend()

# plt.show()
fig.savefig('Conditioning-vs-Pulse_v1.pdf', dpi=300, bbox_inches='tight', format='pdf')