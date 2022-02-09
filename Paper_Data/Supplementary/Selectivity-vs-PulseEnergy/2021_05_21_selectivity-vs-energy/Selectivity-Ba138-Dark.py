# -*- coding: utf-8 -*-
"""
Created on Sat May 22 16:05:18 2021

@author: brend
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import FormatStrFormatter

columnwidth = 225.8775
fullwidth = 469.75502
inches_per_pt = 1/72.27
golden_ratio = (5**.5 - 1) / 2
fig_width_in = fullwidth * inches_per_pt
height_ratio = golden_ratio
fig_height_in = fig_width_in * height_ratio

energies = np.array([82, 95, 115, 140])
fluences = energies*1e-6/(math.pi*(98*1e-4)**2)

ba138_selectivity_82uJ = 27/27
known_ba_selectivity_82uJ = 27/27

ba138_selectivity_95uJ = 32/33
known_ba_selectivity_95uJ = 33/33

ba138_selectivity_115uJ = 32/42
known_ba_selectivity_115uJ = 33/42

ba138_selectivity_140uJ = 20/32
known_ba_selectivity_140uJ = 21/32

selectivity_ba138 = [ba138_selectivity_82uJ, ba138_selectivity_95uJ, ba138_selectivity_115uJ, ba138_selectivity_140uJ]
selectivity_barium_known = [known_ba_selectivity_82uJ, known_ba_selectivity_95uJ, known_ba_selectivity_115uJ, known_ba_selectivity_140uJ]

fig = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300)
ax1 = fig.add_subplot(111)

ax1.plot(fluences, selectivity_ba138, label='Ba138 selectivity', linestyle='', marker='o', markersize=3, alpha=0.6, mew=0.5)
ax1.plot(fluences, selectivity_barium_known, label='Overall barium selectivity less Ba135', linestyle='', marker='o', markersize=3, alpha=0.6, mew=0.5)

ax1.set_xlabel(r"Fluence ($J\: cm^{-2}$)",fontsize=10)
ax1.set_ylabel(f'Selectivity', fontsize=10)
ax1.tick_params(axis='both', which='major', labelsize=6, width=0.5)
# ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.legend(loc='lower left', fontsize=8)

fig.savefig('Selectivity-vs-Fluence_v2.pdf', dpi=300, bbox_inches='tight', format='pdf')