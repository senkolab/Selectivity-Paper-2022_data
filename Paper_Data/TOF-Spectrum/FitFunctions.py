# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:28:38 2020

@author: brend
"""
import math
import numpy as np

Lifetime = 8.36e-9
FreqNat = 1/Lifetime
FreqNatRad = FreqNat/(2*math.pi)
#FreqNatRad *= 1e-6
c = 3e8
kB = 1.38064852e-23

#import all isotope classes
from BariumIsotopesData import *
from TOFSpectrumFunctions import *


def SpectrumIonization(Freq, F0, Saturation, Amp, BG):
    FWHM = FreqNatRad*math.sqrt(1+Saturation)
    Lineshape = 0
    i = 0
    for Isotope in (Ba132, Ba134, Ba135_12, Ba135_32, \
                Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52, Ba138):
        TransStrength = get_TransStrength(Isotope)
        x = 2*(Freq - Isotope.freq_Off - F0)/FWHM
        Lineshape += Isotope.abund*TransStrength/(1 + np.square(x))
        i += 1
    Lineshape *= Amp
    Lineshape += BG
    return Lineshape

def DopplerShift(Angle, Time, d = 14.6e-3):
    return 1 + d*math.sin(Angle)/(c*Time)

#Calculate peak freq. for an isotope based on doppler shift
def DopplerShiftFrequency(Time, Angle, Isotope, d=14.6e-3):
    Doppler = 1 + d*math.sin(Angle)/(c*Time)
    return (Isotope.TransitionFreq())*Doppler

def VelocityDist(Time, Temp, Amp, BG, d=14.6e-3, mass=137.327*1.66054e-27):
    x = d**3*math.e**(-mass*d**2/(2*Time**2*kB*Temp))/Time**3
    x *= Amp
    x += BG
    return x

def TOFDist(Time, Temp, Amp, BG, d=14.6e-3, mass=137.327*1.66054e-27):
    x = d**3*math.e**(-mass*d**2/(2*Time**2*kB*Temp))/Time**2
    x *= Amp
    x += BG
    return x

def TOF_Spectrum_v3(X, F0, Saturation, Angle, Amp=1, Temp=100000, Background=0, d=14.6e-3):
    Freq, Time = X
    #TransFreq = 541.43300e12
    #s = SaturationLevel(Power, Lifetime, TransFreq, BeamWaist)
    FWHM = FreqNatRad*math.sqrt(1+Saturation)
    LineshapeTotal = 0

    for i, Isotope in enumerate([Ba132, Ba134, Ba135_12, Ba135_32, \
                Ba135_52, Ba136, Ba137_12, Ba137_32, Ba137_52, Ba138]):
        TransStrength = get_TransStrength(Isotope)
        DopplerS = DopplerShift(Angle, Time)
        FreqPeak = (Isotope.freq_Off+ F0)*DopplerS
        x = 2*(Freq - FreqPeak)/FWHM
        Lineshape = Isotope.abund*TransStrength/(1 + np.square(x))
        SpeedDistr = d**3*math.e**(-Isotope.mKg*d**2/(2*Time**2*kB*Temp))/Time**3
        SpeedDistr *= Time#Scale by time, since faster atoms fluoresce more
        Lineshape *= SpeedDistr
        LineshapeTotal += Lineshape
    
    LineshapeTotal *= Amp
    LineshapeTotal += Background
    return LineshapeTotal

def TOF_Spectrum_v2(X, F0, Saturation, Angle, AmpMax=1, Temp=100000, d=14.6e-3):
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
        SpeedDistr = d**2*math.e**(-Isotope.mKg*d**2/(2*Time**2*kB*Temp))/Time
        Lineshape *= SpeedDistr
        LineshapeTotal += Lineshape
        i += 1
    
    LineshapeTotal *= AmpMax
    return LineshapeTotal