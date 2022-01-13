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
    #Get the average transition strength of an isotope
    def get_TransStrength(self):
        if self.m in (133, 135, 137):
            TransStrength = float(self.AverageFTransitionStrength())
        else:
            TransStrength = float(self.AverageTransitionStrength())
        return TransStrength
    
