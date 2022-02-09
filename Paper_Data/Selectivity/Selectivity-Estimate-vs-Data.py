# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 16:59:26 2021

@author: brend
"""
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Barium_Isotopes_data import *
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
# from matplotlib import rc
# rc('font', **{'family': 'serif', 'serif': ['Comic Sans']})
# rc('text', usetex=True)
# mpl.rcParams["font.family"] = "Times New Roman"
# mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams.update({'errorbar.capsize': 2})

plt.style.use('seaborn-colorblind')
# sns.set_palette('colorblind')
# sns.set()

lifetime_intermediate = 8.36e-9
freq_natural = 1/lifetime_intermediate
freq_natural_rad = freq_natural/(2*math.pi)

#Physical constants
amu_conv = 1.660539040e-27
h_bar = 1.0545718e-34
kelvin_conv = 273.15
c = 3e8
kB = 1.38064852e-23

plot_ba133 = False
isotopes_all = (Ba138, Ba132, Ba133_g12, Ba133_g32, Ba134, Ba135_12, Ba135_32, \
                    Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52)
isotopes_care = (Ba138, Ba137_32)
isotopes_other = (Ba132, Ba134, Ba135_12, Ba135_32, Ba135_52, Ba136, Ba137_12, Ba137_52)
# line_styles = ('-', '--', '-.', ':')

ba138_label = r'$^{138}\mathrm{Ba}^+$'
ba137_label = r'$^{137}\mathrm{Ba}^+$'
plotcolor_ba137 = 'tab:purple'
plotcolor_ba138 = 'tab:blue'
alpha_ba138 = 0.9
alpha_ba137 = 0.7
markersize = 3
markeredgewidth = 0.5
linewidth = 1
elinewidth = 1
legend_fontsize = 6.5
axis_fontsize = 10
ticks_fontsize = 6
graph_edge_width = 0.5
graph_tick_width = graph_edge_width

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
fig_width_in = columnwidth * inches_per_pt
height_ratio = 0.95
fig_height_in = fig_width_in * height_ratio

def spectrum_Gaussian(Isotope, freq, f0, sat, angle=1.0, a=1, bg=1):
    fwhm_nat = freq_natural_rad*math.sqrt(1+sat)
    fwhm_doppler = calc_doppler_broadening(Isotope, 2000)*1e-6
    fwhm = max(fwhm_nat, fwhm_doppler)
    transition_strength = float(Isotope.get_TransStrength())
        
    xx = 2*(freq - Isotope.freq_Off - f0)/fwhm
    # xx *= 1e-6 #Convert to MHz for display
    lineshape = a*Isotope.abund*transition_strength/(1 + np.square(xx))
    return lineshape

#Calculate the doppler broadening
def calc_doppler_broadening(Isotope, temp, angle=1):
    temp_abs = temp - kelvin_conv
    m = Isotope.m*amu_conv
    theta = angle*math.pi/180
    fwhm = math.sqrt(8*k_B*temp_abs*math.log(2)/(m*c**2))*Isotope.TransitionFreq()*math.cos(theta)
    return fwhm





fig_spectrum = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300)
ax1 = fig_spectrum.add_subplot(111)

freq_range = np.arange(-50e6, 520e6, 100e3)
lineshape_total = np.zeros(freq_range.size)

sat = 5
legend_entry = True
for Isotope in isotopes_other:
    if not plot_ba133 and Isotope.m == 133:
        continue
    lineshape_isotope = spectrum_Gaussian(Isotope, freq_range, 0, sat)
    # a=1/0.23899701287387767
    if legend_entry:
        plott = ax1.plot(freq_range*1e-6, lineshape_isotope, linestyle='--',\
                         alpha=0.3, color='k', label=f'Other isotopes', lw=linewidth)
        legend_entry = False
    else:
        plott = ax1.plot(freq_range*1e-6, lineshape_isotope, linestyle='--', alpha=0.3, color='k', lw=linewidth)
    lineshape_total += lineshape_isotope
    print(f"y max: {ax1.dataLim.bounds[3]}, isotope max: {max(lineshape_isotope)}")
    
for Isotope in isotopes_care:
    if not plot_ba133 and Isotope.m == 133:
        continue
    lineshape_isotope = spectrum_Gaussian(Isotope, freq_range, 0, sat)
    if Isotope.m == 138:
        lineshape_ba138 = lineshape_isotope
        plott = ax1.plot(freq_range*1e-6, lineshape_isotope, alpha=alpha_ba138, lw=linewidth, \
                         color=plotcolor_ba138, label=r'Isotope $^{%i}\mathrm{Ba}$'%(Isotope.m))
    elif Isotope.m == 137:
        lineshape_ba137 = lineshape_isotope
        plott = ax1.plot(freq_range*1e-6, lineshape_isotope, alpha=alpha_ba137, lw=linewidth, \
                         color=plotcolor_ba137, label=r'Isotope $^{%i}\mathrm{Ba}$'%(Isotope.m))
    lineshape_total += lineshape_isotope
    print(f"y max: {ax1.dataLim.bounds[3]}, isotope max: {max(lineshape_isotope)}")
    
ba137_range_min, ba137_range_max = Ba137_32.freq_Off*1e-6 - 50, Ba137_32.freq_Off*1e-6 + 50

ax1.plot(freq_range*1e-6, lineshape_total, color='k', label='Total lineshape', lw=linewidth, alpha=0.7)
# ax1.axhline(y=np.max(ba137_ratio), xmin=0, xmax=np.max(freq_range)*1e-6, color='k')
    
ax1.set(xlim=(-50, 510), ylim=(0, np.max(lineshape_total)))
# ax1.set_yscale('log')
ax1.set_xlabel('Frequency (MHz)',fontsize=axis_fontsize)
ax1.set_ylabel('Neut. fluor. (arb. units)',fontsize=axis_fontsize)
ax1.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax1.xaxis.set_major_locator(MultipleLocator(100))
ax1.xaxis.set_minor_locator(MultipleLocator(50))
ax1.yaxis.set_major_locator(MultipleLocator(0.1))
ax1.yaxis.set_minor_locator(MultipleLocator(0.05))
# ax1.yaxis.set_ticklabels([])

handles, labels = ax1.get_legend_handles_labels()
handles2 = [handles[1], handles[2], handles[0], handles[3]]
labels2 = [labels[1], labels[2], labels[0], labels[3]]

ax1.legend(fontsize=legend_fontsize, loc=3, framealpha=0.7, handles=handles2, labels=labels2)


sat2 = 10
lineshape_total_2 = np.zeros(freq_range.size)
for Isotope in isotopes_other:
    if not plot_ba133 and Isotope.m == 133:
        continue
    lineshape_isotope_2 = spectrum_Gaussian(Isotope, freq_range, 0, sat2)
    lineshape_total_2 += lineshape_isotope_2
    print(f"y max: {ax1.dataLim.bounds[3]}, isotope max: {max(lineshape_isotope)}")
    
for Isotope in isotopes_care:
    if not plot_ba133 and Isotope.m == 133:
        continue
    lineshape_isotope_2 = spectrum_Gaussian(Isotope, freq_range, 0, sat2)
    if Isotope.m == 137:
        lineshape_ba137_2 = lineshape_isotope_2
    lineshape_total_2 += lineshape_isotope_2
    print(f"y max: {ax1.dataLim.bounds[3]}, isotope max: {max(lineshape_isotope)}")
    
ba137_ratio = lineshape_ba137/lineshape_total
ba138_ratio = lineshape_ba138/lineshape_total
ba137_ratio_2 = lineshape_ba137_2/lineshape_total_2

ba137_exp_sat = 0.37
ba137_exp_sat_err = 0.05
ba137_exp_sat2 = 0.32
ba137_exp_sat2_err = 0.05
ba138_exp_sat = 0.85
ba138_exp_sat_err = 0.03
natural_ln_137 = math.log(Ba137_32.abund/(1-Ba137_32.abund))
natural_ln_138 = math.log(Ba138.abund/(1-Ba138.abund))
exp_ln_137_sat = math.log(ba137_exp_sat/(1-ba137_exp_sat))
exp_ln_137_sat2 = math.log(ba137_exp_sat2/(1-ba137_exp_sat2))
exp_ln_138_sat = math.log(ba138_exp_sat/(1-ba138_exp_sat))
ba137_enhancement_sat = exp_ln_137_sat - natural_ln_137
ba137_enhancement_sat2 = exp_ln_137_sat2 - natural_ln_137
ba138_enhancement_sat = exp_ln_138_sat - natural_ln_138

th_sel_label = r'$p_{sel}^{th}$'
exp_sel_label = r'$p_{sel}^{exp}$'

ax2 = ax1.twinx()
ax2.plot(freq_range*1e-6, ba137_ratio, linestyle='--', color=plotcolor_ba137, \
         label=r'%s %s, $s_1$'%(ba137_label, th_sel_label), lw=linewidth, alpha=alpha_ba137)
ax2.errorbar(Ba137_32.freq_Off*1e-6, ba137_exp_sat, yerr=ba137_exp_sat_err, alpha=alpha_ba137,\
             fmt='o', color=plotcolor_ba137, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth, \
                 label=r'%s %s, $s_1$'%(ba137_label, exp_sel_label))
ax2.plot(freq_range*1e-6, ba138_ratio, linestyle='--', color=plotcolor_ba138, \
         label=r'%s %s, $s_1$'%(ba138_label, th_sel_label), lw=linewidth, alpha=alpha_ba138)
ax2.errorbar(Ba138.freq_Off*1e-6, ba138_exp_sat, yerr=ba138_exp_sat_err, alpha=alpha_ba138,\
             fmt='d', color=plotcolor_ba138, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth, \
                 label=r'%s %s, $s_1$'%(ba138_label, exp_sel_label))
ax2.plot(freq_range*1e-6, ba137_ratio_2, linestyle='-.', color=plotcolor_ba137, \
         label=r'%s %s, $s_2$'%(ba137_label, th_sel_label), lw=linewidth, alpha=alpha_ba137)
ax2.errorbar(Ba137_32.freq_Off*1e-6, ba137_exp_sat2, yerr=ba137_exp_sat2_err, alpha=alpha_ba137,\
             fmt='s', color=plotcolor_ba137, ms=markersize, mew=markeredgewidth, elinewidth=elinewidth, \
                 label=r'%s %s, $s_2$'%(ba137_label, exp_sel_label))
ax2.axhline(y=Ba138.abund, xmin=0, xmax=np.max(freq_range)*1e-6, color=plotcolor_ba138, \
            linestyle=':', label=r'$^{138}\mathrm{Ba}^+$ $p_{nat}$', alpha=alpha_ba138)
ax2.axhline(y=Ba137_32.abund, xmin=0, xmax=np.max(freq_range)*1e-6, color=plotcolor_ba137, \
            linestyle=':', label=r'$^{137}\mathrm{Ba}^+$ $p_{nat}$', alpha=alpha_ba137)
ax2.arrow(Ba137_32.freq_Off*1e-6+20, max(ba137_ratio_2), 0, max(ba137_ratio)-max(ba137_ratio_2), \
          length_includes_head=True, head_width=5, head_length=0.01, fc='k', linewidth=linewidth/2)
ax2.text(Ba137_32.freq_Off*1e-6+78, max(ba137_ratio_2)+0.03, 'Lower\nlinewidth', ha='center', fontsize=6)
ax2.arrow(Ba137_32.freq_Off*1e-6+20, ba137_exp_sat2, 0, ba137_exp_sat-ba137_exp_sat2, \
          length_includes_head=True, head_width=5, head_length=0.01, fc='k', linewidth=linewidth/2)

handles, labels = ax2.get_legend_handles_labels()
handles2 = [handles[3], handles[0], handles[2], handles[1], handles[4], handles[5], handles[7], handles[6]]
labels2 = [labels[3], labels[0], labels[2], labels[1], labels[4], labels[5], labels[7], labels[6]]
# handles2 = handles
# labels2 = labels
    
ax2.set_ylabel('Loading Probability',fontsize=axis_fontsize, rotation=270, labelpad=14)
ax2.tick_params(axis='both', which='major', labelsize=ticks_fontsize, width=graph_tick_width)
ax2.yaxis.set_major_locator(MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
ax2.set(ylim=(0, 1.09))
ax2.legend(fontsize=legend_fontsize, loc=1, framealpha=0.7, ncol=2, handles=handles2, labels=labels2)

plt.show()
fig_spectrum.savefig('Selectivity_Ba137_Ba138_v5.pdf', dpi=300, bbox_inches='tight', format='pdf')