# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 17:19:06 2021

@author: brend
"""

from BariumIsotopesData import *
from FitFunctions import *


def PlotDopplerShifts(xDopplerplot, xDoppler, offset, ax1, angle, LegendEntry, plot_linewidth=1):
    for i, Isotope in enumerate([Ba134, Ba136, Ba135_52, Ba137_52, Ba135_32, Ba135_12, Ba137_12]):
        # yDoppler = DopplerShiftFrequency(xDoppler, 1*math.pi/180, Isotope)
        yDoppler = DopplerShiftFrequency(xDoppler, angle, Isotope)
        yDoppler -= offset
        yDoppler *= 1e-6
        colorr = "plum"
        if LegendEntry:
            if i == 0:
                l1, = ax1.plot(xDopplerplot, yDoppler, color=colorr, label=f'Doppler shift ${angle*180/math.pi:0.1f}^\circ$', \
                               lw=plot_linewidth)
            else:
                ax1.plot(xDopplerplot, yDoppler, color=colorr, lw=plot_linewidth)
        else:
            l1 = 0
            ax1.plot(xDopplerplot, yDoppler, color='k', lw=plot_linewidth)
    for i, Isotope in enumerate([Ba137_32, Ba138]):
        colorr = 'tab:purple'
        yDoppler = DopplerShiftFrequency(xDoppler, 0*math.pi/180, Isotope)
        yDoppler -= offset
        yDoppler *= 1e-6
        if LegendEntry:
            if i == 0:
                l2, = ax1.plot(xDopplerplot, yDoppler, color=colorr, linestyle='--', label=r'Doppler shift $0^\circ$', \
                               lw=plot_linewidth)
            else:
                ax1.plot(xDopplerplot, yDoppler, color=colorr, linestyle='--', lw=plot_linewidth)
        else:
            l2 = 0
            ax1.plot(xDopplerplot, yDoppler, color=colorr, linestyle='--', lw=plot_linewidth)
        yDoppler = DopplerShiftFrequency(xDoppler, angle, Isotope)
        yDoppler -= offset
        yDoppler *= 1e-6
        if LegendEntry:
            if i == 0:
                l3, = ax1.plot(xDopplerplot, yDoppler, color=colorr, label=f'Doppler shift ${angle*180/math.pi:0.1f}^\circ$', \
                               lw=plot_linewidth)
            else:
                ax1.plot(xDopplerplot, yDoppler, color=colorr, lw=plot_linewidth)
        else:
            l3 = 0
            ax1.plot(xDopplerplot, yDoppler, color=colorr, lw=plot_linewidth)
        yDoppler = DopplerShiftFrequency(xDoppler, 1.4*math.pi/180, Isotope)
        yDoppler -= offset
        yDoppler *= 1e-6
        if LegendEntry:
            if i == 0:
                l4, = ax1.plot(xDopplerplot, yDoppler, color=colorr, linestyle='-.', label=r'Doppler shift $1.4^\circ$', \
                               lw=plot_linewidth)
            else:
                ax1.plot(xDopplerplot, yDoppler, color=colorr, linestyle='-.', lw=plot_linewidth)
        else:
            l4 = 0
            ax1.plot(xDopplerplot, yDoppler, color=colorr, linestyle='-.', lw=plot_linewidth)
    return ax1