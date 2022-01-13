# -*- coding: utf-8 -*-
"""
Created on Fri May 22 16:41:22 2020

@author: brend
"""

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import os
import math

with open('Use-1-spot.csv', 'r') as myFile:
    data = list(csv.reader(myFile, delimiter=','))
del data[0]
    
time = np.array([])
counts = np.array([])
i = 1;
for points in data:
    if i<1000 and i%5 == 0:
        time = np.append(time, float(points[0]))
        counts = np.append(counts, float(points[1]))
    else:
        pass
    i+=1
    
time = (time - time[0])*2
#plt.figure(figsize=(20,10))
plt.plot(time,counts, marker='o', linestyle='None', mfc='#1dca53', mec='#1dca53')
plt.text(11700, 5000, "2 Hz, 50 uJ Pulse Energy")
plt.xlabel("Pulse Count")
plt.ylabel("Counts (a.u.)")
plt.savefig('Pitting_02-12-2020_v2.png')
plt.show()