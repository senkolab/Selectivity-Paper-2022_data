# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 14:55:40 2020

@author: brend
"""

import math
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j

#Physical constants
amu_Conv = 1.660539040e-27
c = 2.99792458e8
k_B = 1.38064852e-23
h_bar = 1.0545718e-34
kelvin_Conv = 273.15

OffsetAll = 3

#Class for different isotopess 553 nm transitions
class Transition:
    def __init__(self, Isotope, AtomicNumbersFine, Freq_Offset=0, Freq=541.4329619e12, HyperFine = [0, 0], ZeemanNumbers = [0, 0], Abundance = 1, GraphVertOff = 0, alt_label=''):
        self.m = Isotope
        self.mKg = Isotope*1.66054e-27
        self.I = AtomicNumbersFine[0]
        self.J = AtomicNumbersFine[1]
        self.Jp = AtomicNumbersFine[2]
        self.F = HyperFine[0]
        self.Fp = HyperFine[1]
        self.m_F = ZeemanNumbers[0]
        self.m_Fp = ZeemanNumbers[1]
        self.freq = Freq
        self.freqTHz = Freq
        self.freq_Off = Freq_Offset
        self.abund = Abundance
        self.gvert = GraphVertOff + OffsetAll
        self.alt_label = alt_label
    def S_FF(self, m_F):
        Coeff = math.sqrt((2*self.Fp + 1)*(2*self.F + 1)*(2*self.J+1))
        SixJ = wigner_6j(self.J, self.Jp, 1, self.Fp, self.F, self.I)
        return Coeff*SixJ
    def TransitionStrength(self, m_F, m_Fp, F=False, Fp=False):
        if not F:
            F = self.F
        if not Fp:
            Fp = self.Fp
        Coeff = self.S_FF(m_F)*(-1)**(self.J+self.I-m_F)
        ThreeJ = wigner_3j(self.Fp, 1, self.F, m_Fp, (m_F - m_Fp), -m_F)
        return abs(Coeff*ThreeJ)    
    def AverageFTransitionStrength(self):
        Strengths = np.array([])
        for m_F in np.arange(-self.F, self.F+1, 1):
            for m_Fp in np.arange(-self.Fp, self.Fp, 1):
                Strength = self.TransitionStrength(m_F, m_Fp).evalf()
                Strengths = np.append(Strengths, [Strength])
        #Strengths = np.ma.masked_equal(Strengths,0)
        Strengths = np.square(Strengths)
        # if self.m == 137 and self.alt_label == 'b':
        #     print(Strengths)
        Average_Strength = np.mean(Strengths)
        return Average_Strength    
    def AverageTransitionStrength(self):
        Fs = np.arange(abs(min(self.I - self.J, self.J - self.I)), self.I + self.J + 0.1, 1)
        Fps = np.arange(abs(min(self.I - self.Jp, self.Jp - self.I)), self.I + self.Jp + 0.1, 1)
        Strengths = np.array([])
        for F in Fs:
            self.F = F
            for Fp in Fps:
                self.Fp = Fp
                Strength = self.AverageFTransitionStrength()
                Strengths = np.append(Strengths, [Strength])
        Average_Strength = np.mean(Strengths)
        return Average_Strength    
    def TransitionFreq(self):
        return self.freq + self.freq_Off    
    def ChangeHyperFine(self, F, Fp):
        self.F = F
        self.Fp = Fp        
    def Add_gvert(self, Add):
        self.gvert += Add
    def label(self):
        return f'.$^{{{self.m}}}\mathrm{{Ba}}_{{{self.alt_label}}}$'
    def ChangeFreq(self, peak_freq):
        self.freq = peak_freq
        self.freqTHz = self.freq
        
#Get the average transition strength of an isotope
def get_TransStrength(Isotope):
    if Isotope.m in (133, 135, 137):
        TransStrength = float(Isotope.AverageFTransitionStrength())
    else:
        TransStrength = float(Isotope.AverageTransitionStrength())
    return TransStrength
    
#Calculate the doppler broadening
def CalculateDopplerBroadening(Temp, Isotope, Angle):
    T = Temp - kelvin_Conv
    m = Isotope.m*amu_Conv
    theta = Angle*math.pi/180
    FWHM = math.sqrt(2*k_B*T/m)*(Isotope.TransitionFreq()*math.cos(theta)/c)
    return FWHM

#Function to get all data from data files
def GetData(Path):
    data_files = glob.glob(Path)#Grab the data path
    data_files = np.asarray(data_files)#Turn the files into an array of files
    data_Avgs = np.array([])#Setup all of the arrays to save data in
    data_StDevs = np.array([])
    dataFreqs = np.array([])
    dataTimes = np.array([])
    dataAve_Cal = np.array([])
    dataVar_Cal = np.array([])
    data_Calibs = np.array([])
    for i in range(0, data_files.size):#Go through the data files
        data = pd.read_csv(data_files[i], delimiter = "[\],\[\t]", engine='python',header=None)
        data = np.asarray(data)
        indicesBad = np.where(np.logical_not(np.isnan(data))[0] == False)#Delete NaN data
        data = np.delete(data, indicesBad, 1)
        dataFreq = data[:, 1]
        dataFreq = np.round(dataFreq, 6)
        dataTime = data[:, 2]
        dataAve = data[:, 3]
        dataVar = data[:, 5]
        dataRaw = data[:, 6:]
        calibFreq = dataFreq[0]
        for j in range(1, len(dataFreq), 2):#Go through and do calibration scaling
            CalibrationScale1 = dataAve[j-1]
            VarCalib1 = dataVar[j-1]
            CalibrationScale2 = dataAve[j+1]
            VarCalib2 = dataVar[j+1]
            Calibration = (CalibrationScale1 + CalibrationScale2)/2#Calibrate scale is average of surrounding calibrations
            # VarCalib = 0.25*VarCalib1 + 0.25*VarCalib2 + 
            #print("Frequency: %f THz, Time: %i us, Counts: %f, Calibration Scaling: %f"%(dataFreq[j], dataTime[j], dataAve[j], Calibration))
            dataFreqs = np.append(dataFreqs, dataFreq[j])
            dataTimes = np.append(dataTimes, dataTime[j])
            dataAve = np.mean(dataRaw, axis=1)
            dataStd = np.std(dataRaw, axis=1)
            dataAve_Cal = np.append(dataAve_Cal, dataAve[j]/Calibration)#Calibrate the data
            #dataAve_UnCal = np.append(dataAve_UnCal, dataAve[j])
            if j == 1:#For the first data, add first 2 calibrations
                data_Calibs = np.append(data_Calibs, [CalibrationScale1, CalibrationScale2])
            else:
                data_Calibs = np.append(data_Calibs, CalibrationScale2)
            dataVar_Cal = np.append(dataVar_Cal, dataVar[j]/Calibration)
    return (dataFreqs, dataTimes, dataAve_Cal, dataVar_Cal, data_Calibs)