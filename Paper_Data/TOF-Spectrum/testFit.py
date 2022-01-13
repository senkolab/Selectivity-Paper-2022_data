# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 13:13:17 2021

@author: brend
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.optimize import curve_fit
from TOFSpectrumFunctions import *
from BariumIsotopesData import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

Lifetime = 8.36e-9
FreqNat = 1/Lifetime
FreqNatRad = FreqNat/(2*math.pi)
#FreqNatRad *= 1e-6
c = 3e8
kB = 1.38064852e-23
cmap = plt.get_cmap('PuBu_r')

def DopplerShift(Angle, Time, d = 17e-3):
    return 1 + d*math.sin(Angle)/(c*Time)

def VelocityDist(Time, Temp, Amp, BG, d=17e-3, mass=138*1.66054e-27):
    x = d**3*math.e**(-mass*d**2/(2*Time**2*kB*Temp))/(Time**3)
    x *= Amp
    x += BG
    return x

def TOF_Spectrum(X, F0, Saturation, Angle, AmpMax=1, Temp=100000, d=17e-3):
    Freq, Time = X
    #TransFreq = 541.43300e12
    #s = SaturationLevel(Power, Lifetime, TransFreq, BeamWaist)
    FWHM = FreqNatRad*math.sqrt(1+Saturation)
    LineshapeTotal = 0

    i = 0
    for Isotope in (Ba132, Ba134, Ba135_12, Ba135_32, \
                Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52, Ba138):
        TransStrength = get_TransStrength(Isotope)
        DopplerS = DopplerShift(Angle, Time)
        FreqPeak = (Isotope.freq_Off+ F0)*DopplerS
        x = 2*(Freq - FreqPeak)/FWHM
        Lineshape = Isotope.abund*TransStrength/(1 + np.square(x))
        SpeedDistr = d**3*math.e**(-Isotope.mKg*d**2/(2*Time**2*kB*Temp))/(Time**3)
        Lineshape *= SpeedDistr
        LineshapeTotal += Lineshape
        i += 1
    
    LineshapeTotal *= AmpMax
    return LineshapeTotal

offset = 541.43*1e12
y, x = np.arange(541.432880*1e12, 541.433800*1e12 + 1, 1e6), np.arange(0.1e-6, 38e-6, 0.1e-6)
xslice = slice(541.432880*1e12, 541.433800*1e12, 1e6)
yslice = slice(0.1e-6, 38e-6, 0.1e-6)
yy, xx = np.mgrid[xslice, yslice]
Freq0 = 541.43297*1e12
Saturation0 = 20
Angle0 = 1/180*math.pi
Temperature0 = 2000




zz = TOF_Spectrum((yy, xx), Freq0, Saturation0, Angle0, Temp=Temperature0)
zz *= 1e-4
zz += 0.2
yy = np.round((yy-offset)*1e-6, 1)
xx = (xx + 142e-6)


#F1G1: theoretical TOF spectrum 2D
fig = plt.figure(figsize=(30, 10))
ax1 = fig.add_subplot(131)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cf = ax1.pcolormesh(xx*1e6, yy, zz, norm=colors.LogNorm(vmin=zz.min(), vmax=zz.max()), shading='nearest', cmap=cmap)
fig.colorbar(cf, cax=cax, orientation='vertical')



#F1G2: Integrated over frequencies
AblationSignalTime = 142*1e-6
AblationSignalTimeCut = AblationSignalTime + 3e-6
x += AblationSignalTime
ZTOFsummed = zz[:, x > AblationSignalTimeCut]
xTOFsummed = x[x > AblationSignalTimeCut]
ZTOFsummed = np.sum(ZTOFsummed, 0)
AmpSummed = 1e-6
ZTOFsummed *= AmpSummed
ax2 = fig.add_subplot(132)
ax2.plot(xTOFsummed*1e6, ZTOFsummed, alpha=0.5)
ax2.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
ax2.set_ylabel('Counts (photons)',fontsize=30, labelpad=23)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set(xlim=(min(x)*1e6, max(x)*1e6), ylim=(min(ZTOFsummed), max(ZTOFsummed)))



#F1G3: Scaled by 1/time
xActualTimes = xTOFsummed - AblationSignalTime
ZVelocity = ZTOFsummed/xActualTimes
VelocityDataScale = ZTOFsummed[0]/ZVelocity[0]
ZVelocity *= VelocityDataScale
ax3 = fig.add_subplot(133)
ax3.plot(xActualTimes*1e6, ZVelocity, alpha=0.5)
ax3.set_xlabel('Time ($\mu s$)',fontsize=30, labelpad=20)
ax3.set_ylabel('Counts (photons)',fontsize=30, labelpad=23)
ax3.tick_params(axis='both', which='major', labelsize=20)


ZVelocityMB = VelocityDist(xActualTimes, Temperature0, 10, 0)
Amp = 1e-10
ZVelocityMB *= Amp
ax3.plot(xActualTimes*1e6, ZVelocityMB, alpha=0.5)