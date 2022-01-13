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
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
graph_edge_width = 0.5
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['axes.linewidth'] = graph_edge_width
plt.style.use('seaborn-colorblind')

neutralfluor_color = 'tab:blue'
ionfluor_color = 'tab:orange'
pressure_color = 'tab:purple'
markersize = 2
markeredgewidth = 0.5
linewidth = 1
elinewidth = 1
legend_fontsize = 5
axis_fontsize = 10
ticks_fontsize = 6
graph_tick_width = graph_edge_width
regular_fluence_hatch = '-----'
condition_fluence_hatch = '/////'
alpha_marker = 0.7
hash_alpha = 0.3

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
fig_width_in = columnwidth * inches_per_pt
golden_ratio = (5**.5 - 1) / 2

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

PlotAves = False


BG_data_path = r"BG_*"
BG_Aves, BG_StDevs, BG_data_energies, BG_times_raw, BG_data_raw = GetData(BG_data_path, OverrideLength = 50)

Neutral_data_path = r"NeutralFluorescence_*"
Neutral_Aves, Neutral_StDevs, Neutral_data_energies, Neutral_times_raw, Neutral_data_raw = GetData(Neutral_data_path)
Neutral_data_fluences = Neutral_data_energies*1e-6/(math.pi*(98*1e-4)**2)
Neutral_Maxes = np.amax(Neutral_data_raw, 1)
Neutral_Vars = Neutral_StDevs/np.sqrt(2) #2 full sweeps
Neutral_Aves -= BG_Aves
Neutral_Maxes -= BG_Aves
Neutral_Max = max(Neutral_Maxes)
Neutral_Aves /= Neutral_Max
Neutral_Maxes /= Neutral_Max
Neutral_Vars /= Neutral_Max

Ion_data_path = r"IonFluorescence_*"
Ion_Aves, Ion_StDevs, Ion_data_energies, Ion_times_raw, Ion_data_raw = GetData(Ion_data_path)
Ion_data_fluences = Ion_data_energies*1e-6/(math.pi*(98*1e-4)**2)
Ion_Maxes = np.amax(Ion_data_raw, 1)
Ion_Vars = Ion_StDevs/np.sqrt(2) #2 full sweeps
Ion_Aves -= BG_Aves
Ion_Maxes -= BG_Aves
Ion_Max = max(Ion_Maxes)
Ion_Aves /= Ion_Max
Ion_Maxes /= Ion_Max
Ion_Vars /= Ion_Max




Pressure_data_path = r"Pressure_*"
Pressure_Aves, Pressure_StDevs, Pressure_data_energies, Pressure_times_raw, Pressure_data_raw = GetPressureData(Pressure_data_path)
Pressure_data_fluences = Pressure_data_energies*1e-6/(math.pi*(98*1e-4)**2)
Pressure_Maxes = np.amax(Pressure_data_raw, 1)
Pressure_Vars = Pressure_StDevs/np.sqrt(Pressure_data_raw.shape[1]) #2 full sweeps





# #Plot with all on one figure - just maxes
# height_ratio = golden_ratio
# fig_height_in = fig_width_in * height_ratio

# fig1, ax1 = plt.subplots(1, figsize=(fig_width_in, fig_height_in), dpi=300)
# if PlotAves:
#     plot1 = ax1.errorbar(Neutral_data_fluences, Neutral_Aves, yerr=Neutral_Vars, \
#                          label='Neutral fluor. ave.', fmt='o', color='g', alpha=0.6, \
#                              ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
#     plot2 = ax1.plot(Neutral_data_fluences, Neutral_Maxes, label='Neutral fluor. max', \
#                      linestyle='', marker='s', color='g', alpha=0.6, ms=markersize, mew=markeredgewidth)
# else:
#     plot2 = ax1.plot(Neutral_data_fluences, Neutral_Maxes, label='Neutral fluor. max', \
#                      linestyle='', marker='s', color='g', alpha=0.6, ms=markersize, mew=markeredgewidth)
# if PlotAves: 
#     plot3 = ax1.errorbar(Ion_data_fluences, Ion_Aves, yerr=Ion_Vars, label='Ions Average', \
#                          fmt='o', color='b', alpha=0.6, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
#     plot4 = ax1.plot(Ion_data_fluences, Ion_Maxes, label='Ions Max', linestyle='', \
#                      marker='s', color='b', alpha=0.6, ms=markersize, mew=markeredgewidth)
# else:
#     plot4 = ax1.plot(Ion_data_fluences, Ion_Maxes, label='Ions Max', linestyle='', \
#                      marker='s', color='b', alpha=0.6, ms=markersize, mew=markeredgewidth)

# rect = patches.Rectangle((0.1, 0), 0.1, 1.2, linewidth=0, edgecolor='grey', \
#                          facecolor='none', label='Regular fluence', hatch=regular_fluence_hatch, alpha=0.8)
# ax1.add_patch(rect)
# rect = patches.Rectangle((0.3, 0), 0.4, 1.2, linewidth=0, edgecolor='grey', \
#                          facecolor='none', label='Condition fluence', hatch=condition_fluence_hatch, alpha=0.8)
# ax1.add_patch(rect)

# handles, labels = plt.gca().get_legend_handles_labels()
# if PlotAves:
#     order = [0,2,1,3,4,5]
# else:
#     order = [0,1,2,3]
# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=legend_fontsize, loc=2)
# # plt.legend(fontsize=25, loc=2)

# ax2 = ax1.twinx()
# if PlotAves:
#     plott = ax2.errorbar(Pressure_data_fluences, Pressure_Aves, yerr=Pressure_Vars, \
#                          label='Pressure Average', fmt='o', color='#D55E00', alpha=alpha_marker, \
#                              markersize=13, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
#     ax2.plot(Pressure_data_fluences, Pressure_Maxes, label='Pressure Max', linestyle='', \
#              marker='s', color='#D55E00', alpha=alpha_marker, ms=markersize, mew=markeredgewidth)
# else:
#     ax2.plot(Pressure_data_fluences, Pressure_Maxes, label='Pressure Max', linestyle='', 
#              marker='s', color='#D55E00', alpha=alpha_marker, ms=markersize, mew=markeredgewidth)
# plt.legend(fontsize=legend_fontsize, loc=1)

# ax1.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
# ax1.set_ylim(0, 1.2)
# ax1.set_xlim(0.08, 0.42)
# ax1.set_xlabel(r'Fluence / $J\cdot cm^{-2}$',fontsize=axis_fontsize)
# ax1.set_ylabel('Fluor. Counts (arb. units)',fontsize=axis_fontsize)

# from matplotlib import ticker
# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True) 
# formatter.set_powerlimits((-1,1)) 
# ax2.yaxis.set_major_formatter(formatter) 
# ax2.yaxis.offsetText.set_fontsize(ticks_fontsize)
# ax2.set_ylabel('Pressure / mbar', fontsize=axis_fontsize)
# ax2.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
# ax2.set_ylim(0.9e-10, 1e-9)
# fig1.tight_layout()







# #Separate plots for ions, neutrals
# height_ratio = 0.95
# fig_height_in = fig_width_in * height_ratio

# fig1, ax1 = plt.subplots(1, figsize=(fig_width_in, fig_height_in/2), dpi=300)
# plott = ax1.errorbar(Neutral_data_fluences, Neutral_Aves, yerr=Neutral_Vars, \
#                      label='Sweep ave.', fmt='o', color=neutralfluor_color, alpha=alpha_marker, \
#                          ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
# plotcolor = plott.get_children()[0].get_color()
# ax1.plot(Neutral_data_fluences, Neutral_Maxes, label='Sweep max.', linestyle='', \
#          marker='s', color=plotcolor, alpha=alpha_marker, ms=markersize, mew=markeredgewidth)
# plt.legend(fontsize=legend_fontsize, loc=2)
# ax1.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
# ax1.set_ylim(0, 1.1)
# ax1.set_xlim(0.08, 0.42)
# ax1.set_xlabel(r'Fluence / $J\cdot cm^{-2}$',fontsize=axis_fontsize)
# ax1.set_ylabel('Neutral fluor. (arb. units)',fontsize=axis_fontsize)
# rect = patches.Rectangle((0.1, 0), 0.1, 1.2, linewidth=0, edgecolor='grey', \
#                          facecolor='none', label='Regular fluence', hatch=regular_fluence_hatch, alpha=0.8)
# ax1.add_patch(rect)
# rect = patches.Rectangle((0.3, 0), 0.4, 1.2, linewidth=0, edgecolor='grey', \
#                          facecolor='none', label='Condition fluence', hatch=condition_fluence_hatch, alpha=0.8)
# ax1.add_patch(rect)

# fig2, ax2 = plt.subplots(1, figsize=(fig_width_in, fig_height_in/2), dpi=300)
# ln1 = ax2.errorbar(Ion_data_fluences, Ion_Aves, yerr=Ion_Vars, label='Sweep ave.', \
#                    fmt='o', color=ionfluor_color, alpha=alpha_marker, \
#                        ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
# plotcolor = ln1.get_children()[0].get_color()
# ax2.plot(Ion_data_fluences, Ion_Maxes, label='Sweep max.', linestyle='', marker='s', \
#          color=plotcolor, alpha=alpha_marker, ms=markersize, mew=markeredgewidth)
# plt.legend(fontsize=legend_fontsize, loc=2)
# ax3 = ax2.twinx()
# ln3 = ax3.errorbar(Pressure_data_fluences, Pressure_Aves, yerr=Pressure_Vars, label='Pressure ave.', \
#                    fmt='o', color=pressure_color, alpha=alpha_marker, \
#                        ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
# plotcolor = ln3.get_children()[0].get_color()
# ax3.plot(Pressure_data_fluences, Pressure_Maxes, label='Pressure max.', linestyle='', \
#          marker='s', color=plotcolor, alpha=alpha_marker, ms=markersize, mew=markeredgewidth)
# rect = patches.Rectangle((0.1, 0), 0.1, 1.2, linewidth=0, edgecolor='grey', facecolor='none', \
#                          label=None, hatch=regular_fluence_hatch, alpha=0.8)
# ax2.add_patch(rect)
# rect = patches.Rectangle((0.3, 0), 0.4, 1.2, linewidth=0, edgecolor='grey', facecolor='none', \
#                          label=None, hatch=condition_fluence_hatch, alpha=0.8)
# ax2.add_patch(rect)

# # ask matplotlib for the plotted objects and their labels
# lines, labels = ax2.get_legend_handles_labels()
# lines2, labels2 = ax3.get_legend_handles_labels()
# ax2.legend(lines + lines2, labels + labels2, loc=2, fontsize=legend_fontsize)
# # plt.legend(fontsize=legend_fontsize, loc=2)
# ax2.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
# ax2.set_ylim(0, 1.1)
# ax2.set_xlim(0.08, 0.42)
# ax2.set_xlabel(r'Fluence / $J\cdot cm^{-2}$',fontsize=axis_fontsize)
# ax2.set_ylabel('Ion fluor. (arb. units)',fontsize=axis_fontsize)
# from matplotlib import ticker
# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True) 
# formatter.set_powerlimits((-1,1)) 
# ax3.yaxis.set_major_formatter(formatter) 
# ax3.yaxis.offsetText.set_fontsize(ticks_fontsize)
# ax3.set_ylabel('Pressure / mbar', fontsize=axis_fontsize)
# ax3.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
# ax3.set_ylim(0.9e-10, 1e-9)
# fig2.tight_layout()









#Separate plots for ions, neutrals - same figure
height_ratio = 0.95
fig_height_in = fig_width_in * height_ratio

fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(fig_width_in, fig_height_in), dpi=300, sharex=True)
plott = ax1.errorbar(Neutral_data_fluences, Neutral_Aves, yerr=Neutral_Vars, label='Neutral ave.', \
                     fmt='o', color=neutralfluor_color, alpha=alpha_marker, \
                         ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
plotcolor = plott.get_children()[0].get_color()
ax1.plot(Neutral_data_fluences, Neutral_Maxes, label='Neutral max.', linestyle='', \
         marker='s', color=plotcolor, alpha=alpha_marker, ms=markersize, mew=markeredgewidth)
    
rect = patches.Rectangle((0.3, 0), 0.4, 1.2, linewidth=0, edgecolor='grey', facecolor='none', \
                         label='Condition fluence', hatch=condition_fluence_hatch, alpha=hash_alpha)
ax1.add_patch(rect)

handles, labels = ax1.get_legend_handles_labels()
handles2 = [handles[0], handles[2], handles[1]]
labels2 = [labels[0], labels[2], labels[1]]

ax1.legend(fontsize=legend_fontsize, loc=4, framealpha=1, handles=handles2, labels=labels2)
ax1.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.02))
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.set_ylim(0, 1.1)
ax1.set_xlim(0.08, 0.42)
ax1.xaxis.set_label_position('top')
ax1.xaxis.set_ticks_position('top')
# ax1.set_xlabel(r'Fluence / $J\cdot cm^{-2}$',fontsize=axis_fontsize)
# ax1.set_ylabel('Neutral fluor. (arb. units)',fontsize=axis_fontsize)
# rect = patches.Rectangle((0.1, 0), 0.1, 1.2, linewidth=0, edgecolor='grey', facecolor='none', \
#                          label='Regular fluence', hatch=regular_fluence_hatch, alpha=0.8)
# ax1.add_patch(rect)

# fig2, ax2 = plt.subplots(1, figsize=(15,7))

ln1 = ax2.errorbar(Ion_data_fluences, Ion_Aves, yerr=Ion_Vars, label='Ion ave.', \
                   fmt='o', color=ionfluor_color, alpha=alpha_marker, \
                       ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
plotcolor = ln1.get_children()[0].get_color()
ax2.plot(Ion_data_fluences, Ion_Maxes, label='Ion max.', linestyle='', marker='s', \
         color=plotcolor, alpha=alpha_marker, ms=markersize, mew=markeredgewidth)
ax2.legend(fontsize=legend_fontsize, loc=2)

ax3 = ax2.twinx()
ln3 = ax3.errorbar(Pressure_data_fluences, Pressure_Aves, yerr=Pressure_Vars, label='Pressure ave.', \
                   fmt='o', color=pressure_color, alpha=0.9, \
                       ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
plotcolor = ln3.get_children()[0].get_color()
ax3.plot(Pressure_data_fluences, Pressure_Maxes, label='Pressure max.', linestyle='', \
         marker='s', color=plotcolor, alpha=0.9, ms=markersize, mew=markeredgewidth)
# rect = patches.Rectangle((0.1, 0), 0.1, 1.2, linewidth=0, edgecolor='grey', facecolor='none', \
#                          label=None, hatch=regular_fluence_hatch, alpha=0.8)
# ax2.add_patch(rect)
rect = patches.Rectangle((0.3, 0), 0.4, 1.2, linewidth=0, edgecolor='grey', facecolor='none', \
                         label=None, hatch=condition_fluence_hatch, alpha=hash_alpha)
ax2.add_patch(rect)

# ask matplotlib for the plotted objects and their labels
lines, labels = ax2.get_legend_handles_labels()
lines2, labels2 = ax3.get_legend_handles_labels()
ax2.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax2.xaxis.set_major_locator(MultipleLocator(0.1))
ax2.xaxis.set_minor_locator(MultipleLocator(0.02))
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
ax2.set_ylim(0, 1.1)
ax2.set_xlim(0.08, 0.42)
ax2.set_xlabel(r'Fluence / $J\cdot cm^{-2}$',fontsize=axis_fontsize)
ax2.set_ylabel('Fluorescence / arb. units',fontsize=axis_fontsize)
ax2.yaxis.set_label_coords(-.12, 1.0)

ax2.legend(lines + lines2, labels + labels2, loc=2, fontsize=legend_fontsize)

from matplotlib import ticker
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 
ax3.yaxis.set_major_formatter(formatter) 
ax3.yaxis.offsetText.set_fontsize(ticks_fontsize)
ax3.set_ylabel('Pressure / mbar', fontsize=axis_fontsize, rotation=270, labelpad=14)
ax3.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax3.yaxis.set_major_locator(MultipleLocator(0.2e-9))
ax3.yaxis.set_minor_locator(MultipleLocator(0.1e-9))
ax3.set_ylim(0.9e-10, 1e-9)

ax2.text(0.25, 1.18, '(a)', ha='center', fontsize=8, font='Times New Roman')
ax2.text(0.25, -0.53, '(b)', ha='center', fontsize=8, font='Times New Roman')

fig1.savefig('Neutrals-Vs-Ions_Plot_v5.pdf', dpi=300, bbox_inches='tight', format='pdf')