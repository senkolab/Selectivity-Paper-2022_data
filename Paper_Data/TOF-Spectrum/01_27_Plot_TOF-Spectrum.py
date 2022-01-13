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
temp_init, temp_lower, temp_upper = 30000, 20000, 50000
bg_init, bg_lower, bg_upper = 0.1, 0.01, 1
params_init = np.array([freq0_init, sat_init, angle_init, amp_init, temp_init, bg_init])
lower_bound_array = np.array([freq0_lower, sat_lower, angle_lower, amp_lower, temp_lower, bg_lower])
upper_bound_array = np.array([freq0_upper, sat_upper, angle_upper, amp_upper, temp_upper, bg_upper])
pbounds = np.array([lower_bound_array, upper_bound_array])
fit_variables = np.vstack((freq_range_fit_2D.ravel(), time_range_fit_2D.ravel()))
params, pcov = curve_fit(TOF_Spectrum_v3, fit_variables, counts_dataForFit.ravel(), p0=params_init, bounds=pbounds)
[freq_2D_fit, sat_2D_fit, angle_2D_fit, amp_2D_fit, temp_2D_fit, bg_2D_fit] = params

for isotope in (Ba132, Ba133_g12, Ba133_g32, Ba134, Ba135_12, Ba135_32, \
                Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52, Ba138):
    if not plot_ba133 and isotope.m == 133:
        continue
    isotope.ChangeFreq(freq_2D_fit)

if plot_TOF_spectrum:
    fig_TOF_spect = plt.figure(figsize=(10, 10))
    ax1_TOF_spect = fig_TOF_spect.add_subplot(111)
    divider = make_axes_locatable(ax1_TOF_spect)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cf = ax1_TOF_spect.pcolormesh(time_dataplot-pulse_time*1e6, freq_dataplot, counts_data, shading='nearest', \
                                  norm=colors.LogNorm(vmin=counts_data.min(), vmax=counts_data.max()), cmap=color_map)
    fig_TOF_spect.colorbar(cf, cax=cax, orientation='vertical')
    ax1_TOF_spect.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
    ax1_TOF_spect.set_ylabel('Frequency (MHz) + %3.2f THz'%(freq_plotting_offset*1e-12),fontsize=30, labelpad=23)
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
    l5 = ax1_TOF_spect.vlines([slow_atom_cutoff*1e6-pulse_time*1e6], np.min(freq_dataplot), np.max(freq_dataplot), \
                              color=color_spectrum, label=r'$160\:\mu s$ slow cutoff', linestyle='dotted', linewidth=3)
    # l6 = ax1_TOF_spect.vlines(0, np.min(freq_dataplot), np.max(freq_dataplot), \
        # color=color_TOF, label=r'Ablation pulse time', linestyle='dotted', linewidth=3)
    ax1_TOF_spect.legend(fontsize=15, loc=1, framealpha=1)
    # ax1_TOF_spect.legend(handles = [l1, l2, l3, l4, l5, l6], \
    #            labels=[r'$1^\circ$ Doppler Shift', r'$0^\circ$ Doppler Shift', r'$1^\circ$ Doppler Shift', r'$1.7^\circ$ Doppler Shift'], \
    #                fontsize=12, loc=1, framealpha=1)










#F2: 553 nm spectrum by integrating over times 
time_data_slow = time_data[:, x > slow_atom_cutoff]
# time_data_slow = time_data_slow*1e6
freq_data_slow = freq_dataplot[:, x > slow_atom_cutoff]
freq_data_slowplot = freq_data_slow
counts_data_slow = counts_data[:, x > slow_atom_cutoff]
counts_error_data_slow = counts_dataerr[:, x > slow_atom_cutoff]

counts_spectrum_slow = np.sum(counts_data_slow, 1)
counts_error_data_sum = np.sum(counts_error_data_slow, 1)
counts_error_data_slow = np.mean(counts_error_data_slow, 1)
counts_error_data_slow *= counts_error_data_sum
freq_data_slow_plot = np.round((y - freq_plotting_offset)*1e-6, 1)
freq_data_slow_plot[freq_data_slow_plot>freq_calib_start_datasets] += freq_calib_offset_datasets

freq0_init, freq0_lower, freq0_upper = freq_2D_fit, 541.43295*1e12, 541.43298*1e12
sat_init, sat_lower, sat_upper = 10, 5, 20
amp_init, amp_lower, amp_upper = 100, 50, 250
bg_init, bg_lower, bg_upper = 3, 1, 4
params_init = np.array([freq0_init, sat_init, amp_init, bg_init])
lower_bound_array = np.array([freq0_lower, sat_lower, amp_lower, bg_lower])
upper_bound_array = np.array([freq0_upper, sat_upper, amp_upper, bg_upper])
pbounds = np.array([lower_bound_array, upper_bound_array])
params, pcov = curve_fit(SpectrumIonization, y, counts_spectrum_slow, p0=params_init, bounds=pbounds)
[freq_slow_fit, sat_slow_fit, amp_slow_fit, bg_slow_fit] = params
perr = np.sqrt(np.diag(pcov))

freq_range_fit = np.arange(min(y), max(y) + 100e6, 1e6)
counts_spectrum_fit = SpectrumIonization(freq_range_fit, freq_slow_fit, sat_slow_fit, amp_slow_fit, bg_slow_fit)
freq_range_fit_plot = np.round((freq_range_fit - freq_plotting_offset)*1e-6, 1)

counts_spectrum_fit_summed = np.zeros(freq_range_fit.size)
transition_strengths = np.array([])
isotope_freqs_indices = np.array([])
plotcolor = '0'
spectrum_scaling = 150

if plot_spectrum:
    fig_spectrum = plt.figure(figsize=(10, 10))
    ax2_spectrum = fig_spectrum.add_subplot(111)
    ax2_spectrum.errorbar(counts_spectrum_slow, freq_data_slow_plot, color=color_spectrum, \
                          xerr=counts_error_data_slow, fmt='o', alpha=0.5, label='Slow-atom spectrum')
    ax2_spectrum.plot(counts_spectrum_fit, freq_range_fit_plot, color=color_spectrum, label=f'Fit: s={sat_slow_fit:0.0f}')
    # textt = f"""
    # 553 nm Power: $8\: \mu W$
    # $^{{138}}\mathrm{{Ba}}$ Peak: {freq_slow_fit*1e-12:3.6f} THz
    # Saturation: x{sat_slow_fit:.0f}
    # """
    # ax2.text(30, 3100, textt, fontsize=12)
    # ax2.set_ylabel('Frequency (MHz) + %3.2f THz'%(freq_plotting_offset*1e-12), fontsize=30, labelpad=20)
    ax2_spectrum.yaxis.set_label_position("right")
    ax2_spectrum.yaxis.tick_right()
    ax2_spectrum.set_xlabel('Counts (photons)',fontsize=30, labelpad=23)
    ax2_spectrum.set_ylabel('Frequency (MHz) + %3.2f THz'%(freq_plotting_offset*1e-12),fontsize=30, labelpad=23)
    ax2_spectrum.tick_params(axis='both', which='major', labelsize=20)
    ax2_spectrum.set(xlim=(np.min(counts_spectrum_slow) - 2, np.max(counts_spectrum_slow) + 5), \
                     ylim=(np.min(freq_dataplot), np.max(freq_dataplot)))
    # ax2_spectrum.yaxis.set_ticklabels([])
    
    for isotope in (Ba132, Ba133_g12, Ba133_g32, Ba134, Ba135_12, Ba135_32, \
                    Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52, Ba138):
        if not plot_ba133 and isotope.m == 133:
            continue
        sat = sat_slow_fit
        fwhm_nat = natural_rate_rads_mhz*math.sqrt(1+sat)
        fwhm_doppler = CalculateDopplerBroadening(2000, isotope, 90-angle_2D_fit*180/math.pi)*1e-6
        fwhm = max(fwhm_nat, fwhm_doppler)
        if isotope.m in (133, 135, 137):
            transition_strength = float(isotope.AverageFTransitionStrength())
        else:
            transition_strength = float(isotope.AverageTransitionStrength())
        transition_strengths = np.append(transition_strengths, [transition_strength])    
        # print('Transition Strength of isotope %i: %f'%(isotope.m, transition_strength))
        xx = 2*(freq_range_fit - isotope.freq_Off - freq_slow_fit)*1e-6/fwhm
        counts_spectrum_theory = isotope.abund*transition_strength/(1 + np.square(xx))
        counts_spectrum_theory *= spectrum_scaling
        Index = np.argmax(counts_spectrum_theory)
        isotope_freqs_indices = np.append(isotope_freqs_indices, [Index])
        
        
        if spectrum_isotopes_label:
            if plotcolor == '0':
                plott = ax2_spectrum.plot(counts_spectrum_theory, (freq_range_fit-freq_plotting_offset)*1e-6, \
                                          linestyle='--', alpha=0.3, label='Isotope curves')
                spectrum_isotopes_label = False
            else:
                plott = ax2_spectrum.plot(counts_spectrum_theory, (freq_range_fit-freq_plotting_offset)*1e-6, \
                                          linestyle='--', color=plotcolor, alpha=0.3, label='Each isotope\'s peak')
                spectrum_isotopes_label = False
        else:
            if plotcolor == '0':
                plott = ax2_spectrum.plot(counts_spectrum_theory, (freq_range_fit-freq_plotting_offset)*1e-6, \
                                          linestyle='--', alpha=0.3)
            else:
                plott = ax2_spectrum.plot(counts_spectrum_theory, (freq_range_fit-freq_plotting_offset)*1e-6, \
                                          linestyle='--', color=plotcolor, alpha=0.3)
        plotcolor = plott[0].get_color()
        counts_spectrum_fit_summed += counts_spectrum_theory
        x_max = ax2_spectrum.dataLim.bounds[2] - ax2_spectrum.dataLim.bounds[0]
        ax2_spectrum.axhline(y=(freq_slow_fit - freq_plotting_offset)*1e-6 + isotope.freq_Off, \
                             xmin=0, xmax=(max(counts_spectrum_theory))/(x_max), color=plotcolor, alpha=0.3)
        
    i = 0
    for isotope in (Ba132, Ba133_g12, Ba133_g32, Ba134, Ba135_12, Ba135_32, \
                    Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52, Ba138):
        if not plot_ba133 and isotope.m == 133:
            continue
        isotope_text_offset_add = counts_spectrum_fit[int(isotope_freqs_indices[i])]
        isotope.Add_gvert(isotope_text_offset_add)
        ax2_spectrum.text(isotope.gvert+0.05, (freq_slow_fit - freq_plotting_offset)*1e-6 + isotope.freq_Off*1e-6, \
                          isotope.label(), rotation=25, fontsize=20)
        i += 1
    ax2_spectrum.legend(fontsize=15, loc=1, framealpha=1)
    
    
    
    
    






#F3: TOF data - integrate over frequencies
counts_TOF = counts_data[:, x > time_TOF_cutoff]
time_data_TOF = x[x>time_TOF_cutoff]
time_data_TOF_plot = time_data_TOF*1e6
time_data_TOF = time_data_TOF - pulse_time
counts_dataerr_TOF = counts_dataerr[:, x > time_TOF_cutoff]
counts_TOF = np.sum(counts_TOF, 0)
counts_dataerr_TOF_sum, counts_dataerr_TOF = np.sum(counts_dataerr_TOF, 0), np.mean(counts_dataerr_TOF, 0)
counts_dataerr_TOF *= counts_dataerr_TOF_sum
if plot_TOF:
    fig = plt.figure(figsize=(10, 10))
    ax2 = fig.add_subplot(111)
    # ax2.vlines(0, 0, 140, color=color_TOF, label=r'Ablation pulse time', linestyle='dotted', linewidth=2)
    ax2.errorbar(time_data_TOF_plot-pulse_time*1e6, counts_TOF, yerr=counts_dataerr_TOF, \
                 fmt='o', alpha=0.5, label=f'${fluence_cm2:0.2f}\:J/cm^2$')
    

temp_init, temp_lower, temp_upper = 35000, 25000, 45000
amp_init, amp_lower, amp_upper = 0.01, 0.0001, 1
bg_init, bg_lower, bg_upper = 0.1, 0.0001, 1
params_init = np.array([temp_init, amp_init, bg_init])
lower_bound_array = np.array([temp_lower, amp_lower, bg_lower])
upper_bound_array = np.array([temp_upper, amp_upper, bg_upper])
pbounds = np.array([lower_bound_array, upper_bound_array])
params, pcov = curve_fit(TOFDist, time_data_TOF, counts_TOF, p0=params_init, bounds=pbounds)
[temp_TOF_fit, amp_TOF_fit, bg_TOF_fit] = params
perr = np.sqrt(np.diag(pcov))
timerange_TOF_fit = np.arange(min(time_data_TOF), max(time_data_TOF) + 1e-6, 1e-8)
timerange_TOF_fit_plot = (timerange_TOF_fit + pulse_time)*1e6
counts_TOF_fit = TOFDist(timerange_TOF_fit, temp_TOF_fit, amp_TOF_fit, bg_TOF_fit)

if plot_TOF:
    # ax2.plot(timerange_TOF_fit_plot, counts_TOF_fit, color=color_TOF, label='Fit')
    ax2.legend(fontsize=15, loc=1, framealpha=1)
    ax2.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
    ax2.set_ylabel('Counts (photons)',fontsize=30, labelpad=23)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.set(xlim=(x[0]*1e6-pulse_time*1e6, np.max(time_data_TOF_plot)-pulse_time*1e6), \
            ylim=(np.min(counts_TOF)-5, np.max(counts_TOF)+10))
    plt.subplots_adjust(right=0.855)














#F4: fit scaled data to Maxwell-Boltzmann velocity distribution
if plot_velocity:
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(121)
counts_TOF_scaled = counts_TOF/time_data_TOF
countserr_TOF_scaled = counts_dataerr_TOF/time_data_TOF
counts_TOF_scaled_scaling = counts_TOF[0]/counts_TOF_scaled[0]
counts_TOF_scaled *= counts_TOF_scaled_scaling
countserr_TOF_scaled *= counts_TOF_scaled_scaling
if plot_velocity:
    plott = ax1.errorbar(14.6e-3/time_data_TOF, counts_TOF_scaled, yerr=countserr_TOF_scaled, \
                         fmt='o', alpha=0.5, label=f'${fluence_cm2:0.2f}\:J/cm^2$')

temp_init, temp_lower, temp_upper = 35000, 25000, 45000
amp_init, amp_lower, amp_upper = 1e-9, 1e-11, 1e-7
bg_init, bg_lower, bg_upper = 0.0001, 0.00000001, 0.01
params_init = np.array([temp_init, amp_init, bg_init])
lower_bound_array = np.array([temp_lower, amp_lower, bg_lower])
upper_bound_array = np.array([temp_upper, amp_upper, bg_upper])
pbounds = np.array([lower_bound_array, upper_bound_array])
params, pcov = curve_fit(VelocityDist, time_data_TOF, counts_TOF_scaled, p0=params_init, bounds=pbounds)
[temp_vel_fit, amp_vel_fit, bg_vel_fit] = params
perr = np.sqrt(np.diag(pcov))

time_vel_fit = np.arange(min(time_data_TOF), max(time_data_TOF) + 1e-6, 0.1e-6)
vel_fit_plot = 14.6e-3/time_vel_fit
time_vel_fit_plot = time_vel_fit*1e6
counts_vel_fit = VelocityDist(time_vel_fit, temp_vel_fit, amp_vel_fit, bg_vel_fit)

if plot_velocity:
    plotcolor = plott[0].get_color()
    ax1.plot(vel_fit_plot, counts_vel_fit, color=plotcolor, label=f'Fit: T={temp_vel_fit:0.0f} K')
    yboundd = ax1.get_ybound()
    trappable_cutoff = ax1.vlines([330], yboundd[0], yboundd[1], \
                              color='black', label=r'$330\:m/s$ trap depth', linestyle='dotted', linewidth=3)
    ax1.legend(fontsize=15, loc=4, framealpha=1)
    ax1.set_xlabel('Velocity (m/s)',fontsize=30, labelpad=20)
    ax1.set_ylabel('Intensity (arb. units)',fontsize=30, labelpad=23)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    DiffXRange = np.ptp(time_data_TOF)
    DiffYRange = np.ptp(counts_TOF_scaled)
    ax1.yaxis.set_ticklabels([])
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.set(xlim=(0, 4500), ylim=yboundd)
    # ax1.text(15, 30, f"""
    #          Temperature: {temp_vel_fit:.0f} K
    #          Amplitude: {amp_vel_fit:.2e}
    #          """, color='k', fontsize=12)










#F5: Final figure for publication
#F5G1: TOF spectrum data
color_spectrum = "tab:purple"
color_TOF = "fuchsia"
if plot_all:
    fig = plt.figure(figsize=(20, 20))
    ax1 = fig.add_subplot(221)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cf = ax1.pcolormesh(time_data*1e6, freq_dataplot, counts_data, shading='nearest', \
                        norm=colors.LogNorm(vmin=counts_data.min(), vmax=counts_data.max()), color_map=color_map)
    fig.colorbar(cf, cax=cax, orientation='vertical')
    ax1.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
    ax1.set_ylabel('Frequency (MHz) + %3.2f THz'%(freq_plotting_offset*1e-12),fontsize=30, labelpad=23)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set(xlim=(min(x)*1e6, max(x)*1e6), ylim=(np.min(freq_dataplot), np.max(freq_dataplot)))
    rect = patches.Rectangle((min(x)*1e6, freq_calib_start_datasets + dyd/2), \
                             max(x)*1e6,freq_calib_offset_datasets,linewidth=0,edgecolor='w',facecolor='w')
    #rect = patches.Rectangle((140,3100),50,500,linewidth=1,edgecolor='w',facecolor='w')
    ax1.add_patch(rect)
    l1, l2, l3, l4 = PlotDopplerShifts(time_doppler_range_plot, time_doppler_range, freq_plotting_offset, ax1, True)
    ax1.text(160, 3350, f"""
             Fit parameters: 
             Frequency $^{{138}}\mathrm{{Ba}}$: {freq_2D_fit*1e-12:.6f} THz
             Saturation: {sat_2D_fit:.0f}
             Doppler-shift angle: {angle_2D_fit*180/math.pi:.1f}$^\circ$
             Temperature: {temp_2D_fit:.0f} K
             """, color='w', fontsize=12)
    ax1.text(175, (freq_2D_fit-freq_plotting_offset)*1e-6 - 50, f"$^{{138}}\mathrm{{Ba}}$", color='r', fontsize=15)
    ax1.text(170, (freq_2D_fit-freq_plotting_offset)*1e-6 - 50 + Ba137_32.freq_Off*1e-6, \
             f"$^{{137}}\mathrm{{Ba}} (F'=3/2)$", color='r', fontsize=15)
    ax1.legend(handles = [l1, l2, l3, l4], \
               labels=[r'$1^\circ$ Doppler Shift', r'$0^\circ$ - $1.7^\circ$ Doppler Shift', r'$1^\circ$ Doppler Shift'], \
                   fontsize=15, loc=1, framealpha=1)
    ax1.vlines([slow_atom_cutoff*1e6], np.min(freq_dataplot), np.max(freq_dataplot), color=color_spectrum)
    ax1.vlines([pulse_time*1e6], np.min(freq_dataplot), np.max(freq_dataplot), color=color_TOF)
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    
    #F5G2: 553 nm spectrum - sum over slow atoms
    ax2 = fig.add_subplot(222)
    ax2.errorbar(counts_spectrum_slow, freq_data_slow_plot, xerr=counts_error_data_slow, fmt='o', alpha=0.5, color=color_spectrum)
    ax2.plot(counts_spectrum_fit, ySpectrumfitplot, color='k')
    textt = f"""
    Fit parameters:
    553 nm Power: $8\: \mu W$
    $^{{138}}\mathrm{{Ba}}$ Peak: {freq_slow_fit*1e-12:3.6f} THz
    Saturation: x{sat_slow_fit:.0f}
    """
    ax2.text(30, 3100, textt, fontsize=12)
    # ax2.set_ylabel('Frequency (MHz) + %3.2f THz'%(freq_plotting_offset*1e-12), fontsize=30, labelpad=20)
    ax2.yaxis.set_label_position("right")
    ax2.set_xlabel('Counts (photons)',fontsize=30, labelpad=23)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.set(xlim=(0, 55), ylim=(np.min(freq_dataplot), np.max(freq_dataplot)))
    ax2.yaxis.set_ticklabels([])
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')

    #F5G3: TOF data - summed over frequencies
    ax3 = fig.add_subplot(223)
    ax3.errorbar(time_data_TOF*1e6, counts_TOF, yerr=counts_dataerr_TOF, fmt='o', alpha=0.5, color=color_TOF)
    ax3.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
    ax3.set_ylabel('Counts (photons)',fontsize=30, labelpad=23)
    ax3.tick_params(axis='both', which='major', labelsize=20)
    ax3.set(xlim=(min(x)*1e6, max(x)*1e6), ylim=(5, 120))
    
    #F5G4: TOF data scaled by 1/t, to get velocity distribution
    ax4 = fig.add_subplot(224)
    ax4.errorbar(time_data_TOF*1e6, counts_TOF_scaled, yerr=countserr_TOF_scaled, fmt='o', alpha=0.5)
    ax4.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
    # ax4.set_ylabel('Counts (photons)',fontsize=30, labelpad=23)
    ax4.tick_params(axis='both', which='major', labelsize=20)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')