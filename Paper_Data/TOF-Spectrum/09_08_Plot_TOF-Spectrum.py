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
from scipy.integrate import quad
from TOFSpectrumFunctions import *
from FitFunctions import *
from plotFunctions import *
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib import rc
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
graph_edge_width = 0.5
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['axes.linewidth'] = graph_edge_width
# mpl.rcParams["font.family"] = "Times New Roman"
# rc('text', usetex=True)


plt.style.use('seaborn-colorblind')

#Plotting parameters
color_map = plt.get_cmap('PuBu_r')
color_spectrum = "tab:purple"
color_TOF = plt.rcParams['axes.prop_cycle'].by_key()['color'][2]

plot_alpha = 0.7
alpha_theory_spectrum = 0.3
markersize = 3
markeredgewidth = 0.5
linewidth = 1
elinewidth = 1
legend_fontsize = 7
axis_fontsize = 10
ticks_fontsize = 6
text_fontsize = 6
text_isotopes_fontsize = 6
graph_edge_width = 0.5
graph_tick_width = graph_edge_width

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
fig_width_in = fullwidth * inches_per_pt

#Atom parameters
lifetime = 8.36e-9
natural_rate = 1/lifetime
natural_rate_rads = natural_rate/(2*math.pi)
natural_rate_rads_mhz = natural_rate_rads*1e-6
#553 nm laser parameters
ionization_power = 1e-6
pulse_energy = 75e-6
ablation_waist = 98e-6
fluence = pulse_energy / (math.pi*ablation_waist**2)
fluence_cm2 = fluence/100**2
ionization_waist = 32e-6

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



data_file_09_08 = r"2021_09_08_TOF-Spectrum_75uJ\TOF-Spectrum_8uW_75uJ_5HzRep_5Samp_*"
data_file = data_file_09_08

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
freq0_init, freq0_lower, freq0_upper = 541.433045*1e12, 541.433*1e12, 541.43305*1e12
sat_init, sat_lower, sat_upper = 24, 10, 35
amp_init, amp_lower, amp_upper = 1e-3, 1e-5, 1
angle_init, angle_lower, angle_upper = (1*math.pi/180, 0.2*math.pi/180, 2*math.pi/180)
temp_init, temp_lower, temp_upper = 30000, 10000, 50000
bg_init, bg_lower, bg_upper = 0.05, 0.01, 0.1
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
    height_ratio1 = 0.95
    height_ratio2 = 0.6
    height_ratio = height_ratio1 + height_ratio2
    fig_height_in = (fig_width_in * height_ratio)/2
    fig_TOF_spect, ((ax1_TOF_spect, ax2_spectrum), (ax3_tof, ax4_vel)) = plt.subplots(2, 2, figsize=(fig_width_in, fig_height_in), \
                                                    dpi=300, gridspec_kw={'height_ratios': [0.95, 0.6]})
    divider = make_axes_locatable(ax1_TOF_spect)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cf = ax1_TOF_spect.pcolormesh(time_dataplot-pulse_time*1e6, freq_dataplot, counts_data, shading='nearest', \
                                  norm=colors.LogNorm(vmin=counts_data.min(), vmax=counts_data.max()), cmap=color_map)
    cb = fig_TOF_spect.colorbar(cf, cax=cax, orientation='vertical')
    ax1_TOF_spect = PlotDopplerShifts(time_doppler_range_plot-pulse_time*1e6, time_doppler_range, \
                                      freq_plotting_offset, ax1_TOF_spect, angle_2D_fit, True)
    ax1_TOF_spect.text(175-pulse_time*1e6, (freq_2D_fit-freq_plotting_offset)*1e-6 - 50, \
                       f"$^{{138}}\mathrm{{Ba}}$", bbox=dict(facecolor='w', edgecolor='black'), \
                           color='black', fontsize=text_fontsize, ha='center')
    ax1_TOF_spect.text(175-pulse_time*1e6, (freq_2D_fit-freq_plotting_offset)*1e-6 - 50 + Ba137_32.freq_Off*1e-6, \
                       f"$^{{137}}\mathrm{{Ba}}_b$", bbox=dict(facecolor='w', edgecolor='black'), \
                           color='black', fontsize=text_fontsize, ha='center')
    l5 = ax1_TOF_spect.vlines([slow_atom_cutoff*1e6-pulse_time*1e6], np.min(freq_dataplot), np.max(freq_dataplot), \
                              color=color_spectrum, label=r'Slow-atom cutoff', linestyle='dotted', lw=linewidth)
        
    cb.ax.tick_params(labelsize=ticks_fontsize, axis='both', which='major')
    cb.set_label('Counts / photon', rotation=270, labelpad=14)
    ax1_TOF_spect.set_xlabel('Time / $\mu s$', fontsize=axis_fontsize)
    ax1_TOF_spect.xaxis.set_label_position('top')
    ax1_TOF_spect.xaxis.set_ticks_position('top')
    ax1_TOF_spect.set_ylabel('Frequency / MHz + %3.2f THz'%(freq_plotting_offset*1e-12),fontsize=axis_fontsize)
    ax1_TOF_spect.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    ax1_TOF_spect.xaxis.set_major_locator(MultipleLocator(10))
    ax1_TOF_spect.xaxis.set_minor_locator(MultipleLocator(2))
    ax1_TOF_spect.yaxis.set_major_locator(MultipleLocator(200))
    ax1_TOF_spect.yaxis.set_minor_locator(MultipleLocator(50))
    ax1_TOF_spect.set(xlim=(np.min(time_dataplot)-pulse_time*1e6, np.max(time_dataplot)-pulse_time*1e6), \
                      ylim=(np.min(freq_dataplot), np.max(freq_dataplot)))
        
    ax1_TOF_spect.legend(fontsize=legend_fontsize, loc=1, framealpha=1)
    
    




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
# freq_data_slow_plot[freq_data_slow_plot>freq_calib_start_datasets] += freq_calib_offset_datasets

freq0_init, freq0_lower, freq0_upper = freq_2D_fit, 541.433*1e12, 541.4331*1e12
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
    # height_ratio = 0.95
    # fig_height_in = fig_width_in * height_ratio
    # fig_spectrum = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300)
    # ax2_spectrum = fig_spectrum.add_subplot(111)
    # height_ratio1 = 0.95
    # height_ratio2 = 0.6
    # height_ratio = height_ratio1 + height_ratio2
    # fig_height_in = fig_width_in * height_ratio
    # fig_Spect, (ax2_spectrum, ax4_vel) = plt.subplots(2, 1, figsize=(fig_width_in, fig_height_in), \
    #                                                 dpi=300, gridspec_kw={'height_ratios': [0.95, 0.6]})
    ax2_spectrum.errorbar(counts_spectrum_slow, freq_data_slow_plot, color=color_spectrum, \
                          xerr=counts_error_data_slow, fmt='o', alpha=plot_alpha, label='Slow-atom spectrum', \
                              ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
    ax2_spectrum.plot(counts_spectrum_fit, freq_range_fit_plot, color=color_spectrum, \
                      ms=markersize, mew=markeredgewidth, label=f'Fit: s={sat_slow_fit:0.0f}', lw=linewidth)
    # textt = f"""
    # 553 nm Power: $8\: \mu W$
    # $^{{138}}\mathrm{{Ba}}$ Peak: {freq_slow_fit*1e-12:3.6f} THz
    # Saturation: x{sat_slow_fit:.0f}
    # """
    
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
                                          linestyle='--', alpha=alpha_theory_spectrum, label='Isotope curves', \
                                              lw=linewidth)
                spectrum_isotopes_label = False
            else:
                plott = ax2_spectrum.plot(counts_spectrum_theory, (freq_range_fit-freq_plotting_offset)*1e-6, \
                                          linestyle='--', color=plotcolor, alpha=alpha_theory_spectrum, label='Each isotope\'s peak', \
                                              lw=linewidth)
                spectrum_isotopes_label = False
        else:
            if plotcolor == '0':
                plott = ax2_spectrum.plot(counts_spectrum_theory, (freq_range_fit-freq_plotting_offset)*1e-6, \
                                          linestyle='--', alpha=alpha_theory_spectrum, lw=linewidth)
            else:
                plott = ax2_spectrum.plot(counts_spectrum_theory, (freq_range_fit-freq_plotting_offset)*1e-6, \
                                          linestyle='--', color=plotcolor, alpha=alpha_theory_spectrum, \
                                              lw=linewidth)
        plotcolor = plott[0].get_color()
        counts_spectrum_fit_summed += counts_spectrum_theory
        x_max = ax2_spectrum.dataLim.bounds[2] - ax2_spectrum.dataLim.bounds[0]
        ax2_spectrum.axhline(y=(freq_slow_fit - freq_plotting_offset)*1e-6 + isotope.freq_Off, \
                             xmin=0, xmax=(max(counts_spectrum_theory))/(x_max), color=plotcolor, alpha=alpha_theory_spectrum)
        
    i = 0
    for isotope in (Ba132, Ba133_g12, Ba133_g32, Ba134, Ba135_12, Ba135_32, \
                    Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52, Ba138):
        if not plot_ba133 and isotope.m == 133:
            continue
        isotope_text_offset_add = counts_spectrum_fit[int(isotope_freqs_indices[i])]
        isotope.Add_gvert(isotope_text_offset_add)
        ax2_spectrum.text(isotope.gvert+0.05, (freq_slow_fit - freq_plotting_offset)*1e-6 + isotope.freq_Off*1e-6-20, \
                          isotope.label(), rotation=25, fontsize=text_isotopes_fontsize)
        i += 1
        
    ax2_spectrum.xaxis.set_label_position('top')
    ax2_spectrum.xaxis.set_ticks_position('top')
    ax2_spectrum.yaxis.set_label_position("right")
    ax2_spectrum.yaxis.set_ticks_position("right")
    ax2_spectrum.set_xlabel('Counts / photon',fontsize=axis_fontsize)
    ax2_spectrum.set_ylabel('Frequency / MHz + %3.2f THz'%(freq_plotting_offset*1e-12), \
                            fontsize=axis_fontsize, rotation=270, labelpad=14)
    ax2_spectrum.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    ax2_spectrum.xaxis.set_major_locator(MultipleLocator(20))
    ax2_spectrum.xaxis.set_minor_locator(MultipleLocator(5))
    ax2_spectrum.yaxis.set_major_locator(MultipleLocator(200))
    ax2_spectrum.yaxis.set_minor_locator(MultipleLocator(50))
    ax2_spectrum.set(xlim=(np.min(counts_spectrum_slow) - 2, np.max(counts_spectrum_slow) + 20), \
                     ylim=(np.min(freq_dataplot), np.max(freq_dataplot)+10))
        
    ax2_spectrum.legend(fontsize=legend_fontsize, loc=1, framealpha=1)





#F3: TOF data - integrate over frequencies
counts_TOF = counts_data[:, x > time_TOF_cutoff]
time_data_TOF = x[x>time_TOF_cutoff]
time_data_TOF_plot = time_data_TOF*1e6
time_data_TOF = time_data_TOF - pulse_time
counts_dataerr_TOF = counts_dataerr[:, x > time_TOF_cutoff]
counts_TOF = np.sum(counts_TOF, 0)
counts_dataerr_TOF_sum, counts_dataerr_TOF = np.sum(counts_dataerr_TOF, 0), np.mean(counts_dataerr_TOF, 0)
counts_dataerr_TOF *= counts_dataerr_TOF_sum
    
temp_init, temp_lower, temp_upper = 35000, 25000, 45000
amp_init, amp_lower, amp_upper = 0.01, 0.0001, 1
bg_init, bg_lower, bg_upper = 0.1, 0.001, 1
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
    # ax2.legend(fontsize=legend_fontsize, loc=1, framealpha=1)
    ax3_tof.errorbar(time_data_TOF_plot-pulse_time*1e6, counts_TOF, yerr=counts_dataerr_TOF, \
                 fmt='o', alpha=plot_alpha, label=f'${fluence_cm2:0.2f}\:J/cm^2$', \
                     ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)
        
    ax3_tof.set_xlabel('Time / $\mu s$',fontsize=axis_fontsize)
    ax3_tof.set_ylabel('Counts / photon',fontsize=axis_fontsize)
    ax3_tof.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    ax3_tof.xaxis.set_major_locator(MultipleLocator(10))
    ax3_tof.xaxis.set_minor_locator(MultipleLocator(2))
    ax3_tof.yaxis.set_major_locator(MultipleLocator(20))
    ax3_tof.yaxis.set_minor_locator(MultipleLocator(10))
    ax3_tof.set(xlim=(x[0]*1e6-pulse_time*1e6, np.max(time_data_TOF_plot)-pulse_time*1e6), \
            ylim=(np.min(counts_TOF)-5, np.max(counts_TOF)+13))
        
    pos = ax1_TOF_spect.get_position()
    pos2 = ax3_tof.get_position()
    ax3_tof.set_position([pos.x0, pos2.y0, 0.33, pos2.height])
    ax3_tof.text(17, 137, '(a)', ha='center', fontsize=10, font='Times New Roman')
    ax3_tof.text(17, -40, '(c)', ha='center', fontsize=10, font='Times New Roman')
    
    
    
    
    

#F4: fit scaled data to Maxwell-Boltzmann velocity distribution
# if plot_velocity:
    # height_ratio = 0.6
    # fig_height_in = fig_width_in * height_ratio
    # fig = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300)
    # ax4_vel = fig.add_subplot(111)
counts_TOF_scaled = counts_TOF/time_data_TOF
countserr_TOF_scaled = counts_dataerr_TOF/time_data_TOF
counts_TOF_scaled_scaling = counts_TOF[0]/counts_TOF_scaled[0]
counts_TOF_scaled *= counts_TOF_scaled_scaling
countserr_TOF_scaled *= counts_TOF_scaled_scaling
if plot_velocity:
    plott = ax4_vel.errorbar((14.6e-3/time_data_TOF)*1e-3, counts_TOF_scaled, yerr=countserr_TOF_scaled, \
                         fmt='o', alpha=plot_alpha, label=f'${fluence_cm2:0.2f}\:J/cm^2$', \
                             ms=markersize, mew=markeredgewidth, elinewidth=elinewidth)

temp_init, temp_lower, temp_upper = 35000, 25000, 45000
amp_init, amp_lower, amp_upper = 1e-9, 1e-11, 1e-7
bg_init, bg_lower, bg_upper = 0.01, 0.00001, 0.2
params_init = np.array([temp_init, amp_init, bg_init])
lower_bound_array = np.array([temp_lower, amp_lower, bg_lower])
upper_bound_array = np.array([temp_upper, amp_upper, bg_upper])
pbounds = np.array([lower_bound_array, upper_bound_array])
params, pcov = curve_fit(VelocityDist, time_data_TOF, counts_TOF_scaled, p0=params_init, bounds=pbounds)
[temp_vel_fit, amp_vel_fit, bg_vel_fit] = params
perr = np.sqrt(np.diag(pcov))
temp_vel_fit_rounded = np.round(temp_vel_fit, -3)
errtemp_vel_fit_rounded = np.round(perr[0], -2)

time_vel_fit = np.arange(min(time_data_TOF), max(time_data_TOF) + 1e-6, 1e-9)
vel_fit_plot = 14.6e-3/time_vel_fit
time_vel_fit_plot = time_vel_fit*1e6
counts_vel_fit = VelocityDist(time_vel_fit, temp_vel_fit, amp_vel_fit, bg_vel_fit)

if plot_velocity:
    plotcolor = plott[0].get_color()
    ax4_vel.plot(vel_fit_plot*1e-3, counts_vel_fit, color=plotcolor, 
                 label=f'Fit: T={temp_vel_fit_rounded:0.0f}({errtemp_vel_fit_rounded:0.0f}) K', \
             ms=markersize, mew=markeredgewidth, alpha=plot_alpha, linestyle='--')
    yboundd = ax4_vel.get_ybound()
    trappable_cutoff = ax4_vel.vlines([.330], yboundd[0], yboundd[1], \
                              color='black', label=r'$330\:m/s$ trap depth', linestyle='dotted', lw=linewidth)
    ax4_vel.legend(fontsize=legend_fontsize, loc=4, framealpha=1)
    ax4_vel.set_xlabel('Velocity / $km\cdot s^{-1}$',fontsize=axis_fontsize)
    ax4_vel.set_ylabel('Intensity / arb. units',fontsize=axis_fontsize, rotation=270, labelpad=14)
    ax4_vel.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    ax4_vel.xaxis.set_major_locator(MultipleLocator(1))
    ax4_vel.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax4_vel.yaxis.set_major_locator(MultipleLocator(20))
    ax4_vel.yaxis.set_minor_locator(MultipleLocator(5))
    DiffXRange = np.ptp(time_data_TOF)
    DiffYRange = np.ptp(counts_TOF_scaled)
    # ax4_vel.yaxis.set_ticklabels([])
    ax4_vel.yaxis.tick_right()
    ax4_vel.yaxis.set_label_position("right")
    ax4_vel.set(xlim=(0, 4.5), ylim=yboundd)
    # ax4_vel.text(15, 30, f"""
    #          Temperature: {temp_vel_fit:.0f} K
    #          Amplitude: {amp_vel_fit:.2e}
    #          """, color='k', fontsize=12)
    ax4_vel.text(2.5, 74, '(b)', ha='center', fontsize=10, font='Times New Roman')
    ax4_vel.text(2.5, -27, '(d)', ha='center', fontsize=10, font='Times New Roman')
    
peak_vel = vel_fit_plot[np.argmax(counts_vel_fit)]
counts_vel_fit_upper = VelocityDist(time_vel_fit, temp_vel_fit_rounded+errtemp_vel_fit_rounded, amp_vel_fit, bg_vel_fit)
peak_vel_upper = vel_fit_plot[np.argmax(counts_vel_fit_upper)]
counts_vel_fit_lower = VelocityDist(time_vel_fit, temp_vel_fit_rounded-errtemp_vel_fit_rounded, amp_vel_fit, bg_vel_fit)
peak_vel_lower = vel_fit_plot[np.argmax(counts_vel_fit_lower)]
    
# plt.show()
#fig_TOF_spect.savefig('TOF-Spectrum_Final_75uJ_v2.pdf', dpi=300, bbox_inches='tight', format='pdf')