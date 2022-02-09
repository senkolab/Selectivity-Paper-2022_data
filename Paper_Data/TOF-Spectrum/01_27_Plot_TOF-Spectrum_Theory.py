# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:13:25 2020

@author: brend
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.optimize import curve_fit
from TOFSpectrumFunctions import *
from FitFunctions import *
from plotFunctions import *
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

plt.style.use('seaborn-colorblind')


#Atom parameters
lifetime = 8.36e-9
natural_rate = 1/lifetime
natural_rate_rads = natural_rate/(2*math.pi)
natural_rate_rads_mhz = natural_rate_rads*1e-6
#553 nm laser parameters
ionization_power = 1e-6
pulse_energy = 60e-6
ablation_waist = 98e-6
fluence = pulse_energy / (math.pi*ablation_waist**2)
fluence_cm2 = fluence/100**2
ionization_waist = 32e-6
#Plotting parameters
color_map = plt.get_cmap('PuBu_r')

color_spectrum = "tab:purple"
color_TOF = plt.rcParams['axes.prop_cycle'].by_key()['color'][2]

slow_atom_cutoff = 160e-6
pulse_time = 142.5*1e-6
time_TOF_cutoff = 145e-6

plot_TOF_spectrum = True #TOF spectrum data graph, 2D fit graph, Doppler shifts on both

plot_spectrum = True #TOF spectrum data graph fast atoms cutoff, 553 nm spectrum graph with fit
plot_ba133 = False
spectrum_isotopes_label = True

plot_TOF = True #TOF spectrum data graph, TOF data graph with fit
plot_velocity = True #Velocity distribution graph
plot_all = False #TOF spectrum data graph with Doppler shifts, 553 nm spectrum graph with fit, TOF data graph, Velocity distribution graph



data_file_12_18_01 = r"12-2020to01-2021\TOF-Spectrum_8uW_60uJ_5HzRep_5Samp_*"
data_file = data_file_12_18_01

(freqs_data, times_data, counts_calib_ave, counts_calib_var, counts_calibs) = GetData(data_file)
num_times, num_freqs = len(np.unique(times_data)), len(np.unique(freqs_data))
max_times, min_times, max_freqs, min_freqs = max(times_data), min(times_data), max(freqs_data), min(freqs_data)
dx = np.round(max_times - min_times, 2)/(num_times-1)
dy = np.round(max_freqs - min_freqs, 6)/(num_freqs-1)



#F1: TOF spectrum data
y, x = np.arange(min_freqs, max_freqs + 1e-6, dy)*1e12, np.arange(min_times, max_times + 1e-6, dx)*1e-6
time_data, freq_data = np.meshgrid(x, y)
counts_calib_ave = np.reshape(counts_calib_ave, (num_freqs, num_times))
counts_calib_var = np.reshape(counts_calib_var, (num_freqs, num_times))
counts_data = counts_calib_ave.copy()
counts_dataerr = counts_calib_var.copy()
counts_data *= np.mean(counts_calibs)
lowest_nonzero_counts = np.min(counts_data[counts_data > 0])
counts_dataerr *= np.mean(counts_calibs)
counts_data += lowest_nonzero_counts

freq_plotting_offset = 541.43*1e12#offset for frequency in graphing
freq_dataplot = np.round((freq_data - freq_plotting_offset)*1e-6, 1)#Move freq_data based on graphing offset
# print(freq_data)
freq_calib_offset_datasets = 40#Amount of calibration change between data sets
freq_calib_start_datasets = 3100#Where in frequency the two data sets are separated
freq_dataplot[freq_dataplot>freq_calib_start_datasets] += freq_calib_offset_datasets
time_dataplot = time_data.copy()*1e6

time_doppler_range_plot = np.arange(143, 181, 0.1)
time_doppler_range = (time_doppler_range_plot - 142)*1e-6

#Fit to 2D TOF spectrum
time_range_fit_2D, freq_range_fit_2D = np.meshgrid(x, y)
time_range_fit_2D -= pulse_time
cols_delete = [0, 1, 2, 3]
time_range_fit_2D = np.delete(time_range_fit_2D, cols_delete, 1)
freq_range_fit_2D = np.delete(freq_range_fit_2D, cols_delete, 1)
counts_dataForFit = np.delete(counts_data, cols_delete, 1)
freq0_init, freq0_lower, freq0_upper = 541.43297*1e12, 541.43295*1e12, 541.43298*1e12
sat_init, sat_lower, sat_upper = 24, 10, 35
amp_init, amp_lower, amp_upper = 1e-3, 1e-5, 1
angle_init, angle_lower, angle_upper = (1*math.pi/180, 0.2*math.pi/180, 2*math.pi/180)
temp_init, temp_lower, temp_upper = 20000, 2000, 100000
bg_init, bg_lower, bg_upper = 0.1, 0.01, 1
params_init = np.array([freq0_init, sat_init, angle_init, amp_init, temp_init, bg_init])
lower_bound_array = np.array([freq0_lower, sat_lower, angle_lower, amp_lower, temp_lower, bg_lower])
upper_bound_array = np.array([freq0_upper, sat_upper, angle_upper, amp_upper, temp_upper, bg_upper])
pbounds = np.array([lower_bound_array, upper_bound_array])
fit_variables = np.vstack((freq_range_fit_2D.ravel(), time_range_fit_2D.ravel()))
params, pcov = curve_fit(TOF_Spectrum_v3, fit_variables, counts_dataForFit.ravel(), p0=params_init, bounds=pbounds)
[freq_2D_fit, sat_2D_fit, angle_2D_fit, amp_2D_fit, temp_2D_fit, bg_2D_fit] = params


fig_TOF_spect = plt.figure(figsize=(10, 10))
ax1_TOF_spect = fig_TOF_spect.add_subplot(111)
divider = make_axes_locatable(ax1_TOF_spect)
cax = divider.append_axes('right', size='5%', pad=0.05)
cf = ax1_TOF_spect.pcolormesh(time_dataplot-pulse_time*1e6, freq_dataplot, counts_data, shading='nearest', \
                              norm=colors.LogNorm(vmin=counts_data.min(), vmax=counts_data.max()), cmap=color_map)
fig_TOF_spect.colorbar(cf, cax=cax, orientation='vertical')
ax1_TOF_spect.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
ax1_TOF_spect.set_ylabel('Frequency (MHz + %3.2f THz)'%(freq_plotting_offset*1e-12),fontsize=30, labelpad=23)
ax1_TOF_spect.tick_params(axis='both', which='major', labelsize=20)
ax1_TOF_spect.set(xlim=(np.min(time_dataplot)-pulse_time*1e6, np.max(time_dataplot)-pulse_time*1e6), \
                  ylim=(np.min(freq_dataplot), np.max(freq_dataplot)))
dyd = np.round(dy*1e6, 1)#For rectangle patch on colorplot
rect = patches.Rectangle((np.min(time_dataplot)-pulse_time*1e6, freq_calib_start_datasets + dyd/2), \
                         np.max(time_dataplot)-pulse_time*1e6+3,freq_calib_offset_datasets,linewidth=0,facecolor='w', hatch='x')
ax1_TOF_spect.add_patch(rect)
ax1_TOF_spect = PlotDopplerShifts(time_doppler_range_plot-pulse_time*1e6, time_doppler_range, \
                                  freq_plotting_offset, ax1_TOF_spect, angle_2D_fit, True)
ax1_TOF_spect.text(175-pulse_time*1e6, (freq_2D_fit-freq_plotting_offset)*1e-6 - 50, \
                   f"$^{{138}}\mathrm{{Ba}}$", bbox=dict(facecolor='w', edgecolor='black'), \
                       color='black', fontsize=15)
ax1_TOF_spect.text(170-pulse_time*1e6, (freq_2D_fit-freq_plotting_offset)*1e-6 - 50 + Ba137_32.freq_Off*1e-6, \
                   f"$^{{137}}\mathrm{{Ba}} (F'=3/2)$", bbox=dict(facecolor='w', edgecolor='black'), \
                       color='black', fontsize=15)
# l6 = ax1_TOF_spect.vlines(0, np.min(freq_dataplot), np.max(freq_dataplot), \
    # color=color_TOF, label=r'Ablation pulse time', linestyle='dotted', linewidth=3)
# ax1_TOF_spect.legend(fontsize=15, loc=1, framealpha=1)
    

lowest_data_counts = np.min(counts_data)
xslice = slice(541.432880*1e12, 541.433800*1e12, 1e6)
yslice = slice(0.1e-6, 40e-6, 0.1e-6)
yy, xx = np.mgrid[xslice, yslice]
#zz = TOF_Spectrum((yy, xx), Ba138.TransitionFreq(), 20, 1, 1*math.pi/180)
zz = TOF_Spectrum_v3((yy, xx), freq_2D_fit, sat_2D_fit, angle_2D_fit, Amp=amp_2D_fit, Temp=temp_2D_fit)
amp_data_max = np.max(counts_data)-lowest_data_counts
amp_fit_max = np.max(zz)
scale_zz = amp_data_max/amp_fit_max
zz *= scale_zz
zz += lowest_data_counts
yy = np.round((yy-freq_plotting_offset)*1e-6, 1)
# xx = (xx + 142e-6)

fig_TOF_spect_theory = plt.figure(figsize=(10, 10))
ax2_TOF_spect_theory = fig_TOF_spect_theory.add_subplot(111)
divider = make_axes_locatable(ax2_TOF_spect_theory)
cax = divider.append_axes('right', size='5%', pad=0.05)
cf = ax2_TOF_spect_theory.pcolormesh(xx*1e6, yy, zz, norm=colors.LogNorm(vmin=zz.min(), vmax=zz.max()), shading='nearest', cmap=color_map)
fig_TOF_spect_theory.colorbar(cf, cax=cax, orientation='vertical')

ax2_TOF_spect_theory.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
ax2_TOF_spect_theory.tick_params(axis='both', which='major', labelsize=20)
ax2_TOF_spect_theory.set(xlim=(np.min(time_dataplot)-pulse_time*1e6, np.max(time_dataplot)-pulse_time*1e6), \
                  ylim=(np.min(freq_dataplot), np.max(freq_dataplot)))
ax2_TOF_spect_theory = PlotDopplerShifts(time_doppler_range_plot-pulse_time*1e6, time_doppler_range, \
                                  freq_plotting_offset, ax2_TOF_spect_theory, angle_2D_fit, True)
fit_text = (
    f"FreqFit: {freq_2D_fit*1e-12:.6f} THz, SFit: {sat_2D_fit:.0f}, "
    f"\nAngleFit: {angle_2D_fit*180/math.pi:.1f} Degrees, TempFit: {temp_2D_fit:.0f} K")
ax2_TOF_spect_theory.text(20, 3400, fit_text, fontsize=10, color='w')
# ax2_TOF_spect_theory.legend(handles = [l1, l2, l3, l4], labels=['Other Isotopes', '0 degrees', '1 degrees', '1.7 degrees'], fontsize=15, loc=1)
ax2_TOF_spect_theory.legend(fontsize=15, loc=1, framealpha=1)
plt.show()

# fig_TOF_spect.savefig('TOF-Spectrum_60uJ_supp_v2.pdf', dpi=100, bbox_inches='tight', format='pdf')
# fig_TOF_spect_theory.savefig('TOF-Spectrum_60uJ_theory_supp_v2.pdf', dpi=100, bbox_inches='tight', format='pdf')