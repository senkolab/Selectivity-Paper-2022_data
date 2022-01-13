# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 16:49:42 2020

@author: Ba138_Ion
"""
from Spectrum_functions import *


AtomNums = [0, 0, 1]
Ba132 = Transition(132, AtomNums, Freq_Offset=167.9e6, Abundance=0.001, GraphVertOff=1)
AtomNums = [1/2, 0, 1]
Ba133_g12 = Transition(133, AtomNums, Freq_Offset=-23.3e6, HyperFine=[1/2, 1/2], Abundance=0.01, alt_label='a')
AtomNums = [1/2, 0, 1]
Ba133_g32 = Transition(133, AtomNums, Freq_Offset=386.65e6, HyperFine=[1/2, 3/2], Abundance=0.01, alt_label='b')
AtomNums = [0, 0, 1]
Ba134 = Transition(134, AtomNums, Freq_Offset=142.8e6, Abundance=0.0242, GraphVertOff=1)
AtomNums = [3/2, 0, 1]
Ba135_12 = Transition(135, AtomNums, Freq_Offset=547.3e6, HyperFine=[3/2, 1/2], Abundance=0.0659, GraphVertOff=6, alt_label='a')
AtomNums = [3/2, 0, 1]
Ba135_32 = Transition(135, AtomNums, Freq_Offset=326.7e6, HyperFine=[3/2, 3/2], Abundance=0.0659, GraphVertOff=-0.03, alt_label='b')
AtomNums = [3/2, 0, 1]
Ba135_52 = Transition(135, AtomNums, Freq_Offset=121.6e6, HyperFine=[3/2, 5/2], Abundance=0.0659, GraphVertOff=9, alt_label='c')
AtomNums = [0, 0, 1]
Ba136 = Transition(136, AtomNums, Freq_Offset=128.02e6, Abundance=0.0785, GraphVertOff=3)
AtomNums = [3/2, 0, 1]
Ba137_12 = Transition(137, AtomNums, Freq_Offset=549.47e6, HyperFine=[3/2, 1/2], Abundance=0.1123, GraphVertOff=0, alt_label='a')
AtomNums = [3/2, 0, 1]
Ba137_32 = Transition(137, AtomNums, Freq_Offset=274.56e6, HyperFine=[3/2, 3/2], Abundance=0.1123, GraphVertOff=-0.018, alt_label='b')
AtomNums = [3/2, 0, 1]
Ba137_52 = Transition(137, AtomNums, Freq_Offset=63.43e6, HyperFine=[3/2, 5/2], Abundance=0.1123, GraphVertOff=10, alt_label='c')
AtomNums = [0, 0, 1]
Ba138 = Transition(138, AtomNums, Abundance=0.717, GraphVertOff=2)