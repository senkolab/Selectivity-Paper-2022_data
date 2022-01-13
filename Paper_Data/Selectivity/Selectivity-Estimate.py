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
    

fig_spectrum = plt.figure(figsize=(10, 10))
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
                         alpha=0.1, color='k', label=f'Other isotopes')
        legend_entry = False
    else:
        plott = ax1.plot(freq_range*1e-6, lineshape_isotope, linestyle='--', alpha=0.1, color='k')
    lineshape_total += lineshape_isotope
    print(f"y max: {ax1.dataLim.bounds[3]}, isotope max: {max(lineshape_isotope)}")
    
for Isotope in isotopes_care:
    if not plot_ba133 and Isotope.m == 133:
        continue
    lineshape_isotope = spectrum_Gaussian(Isotope, freq_range, 0, sat)
    plott = ax1.plot(freq_range*1e-6, lineshape_isotope, alpha=1, label=f'Isotope Ba{Isotope.m}')
    lineshape_total += lineshape_isotope
    print(f"y max: {ax1.dataLim.bounds[3]}, isotope max: {max(lineshape_isotope)}")
    
ba137_range_min, ba137_range_max = Ba137_32.freq_Off*1e-6 - 50, Ba137_32.freq_Off*1e-6 + 50

ax1.plot(freq_range*1e-6, lineshape_total, color='k', label='Total lineshape')
# ax1.axhline(y=np.max(ba137_ratio), xmin=0, xmax=np.max(freq_range)*1e-6, color='k')
    
ax1.set(xlim=(-50, 510), ylim=(0, np.max(lineshape_total)))
# ax1.set_yscale('log')
ax1.set_xlabel('Frequency / MHz',fontsize=23, labelpad=23)
ax1.set_ylabel('Neutral fluorescence / arb. units',fontsize=23, labelpad=23)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.legend(fontsize=15, loc=3, framealpha=0.5)


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
    lineshape_total_2 += lineshape_isotope_2
    print(f"y max: {ax1.dataLim.bounds[3]}, isotope max: {max(lineshape_isotope)}")
    
ba137_ratio = lineshape_isotope/lineshape_total
ba137_ratio_2 = lineshape_isotope_2/lineshape_total_2
    
ax2 = ax1.twinx()
ax2.plot(freq_range*1e-6, ba137_ratio, linestyle='-.', label=f'Ba137 selectivity s={sat}')
ax2.plot(freq_range*1e-6, ba137_ratio_2, linestyle='-.', label=f'Ba137 selectivity s={sat2}')
ax2.set_ylabel('Loading Probability',fontsize=23, labelpad=23)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set(ylim=(0, 1))
ax2.legend(fontsize=15, loc=1, framealpha=1)