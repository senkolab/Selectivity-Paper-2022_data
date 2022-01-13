# -*- coding: utf-8 -*-
"""
Created on Sat May 22 16:05:18 2021

@author: brend
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import FormatStrFormatter

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

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(111)

ax1.plot(fluences, selectivity_ba138, label='Ba138 selectivity', linestyle='', marker='o', markersize=13, alpha=0.6)
ax1.plot(fluences, selectivity_barium_known, label='Overall barium selectivity less Ba135', linestyle='', marker='o', markersize=13, alpha=0.6)

ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel(r"Fluence / $J \cdot cm^{-2}$",fontsize=30)
ax1.set_ylabel(f'Selectivity', fontsize=30)
ax1.tick_params(axis='both', which='major', labelsize=20)
# ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.legend(loc='lower left', fontsize=20)