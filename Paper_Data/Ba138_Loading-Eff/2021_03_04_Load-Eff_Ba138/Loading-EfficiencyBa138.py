# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 15:50:57 2021

@author: brend
"""

import math
import numpy as np
import matplotlib.pyplot as plt
#print(plt.style.available)
plt.style.use('seaborn-colorblind')
# plt.rcParams['figure.facecolor'] = '1'
import glob
import pandas as pd
import re

def prepareData(data_files): #Inputs - data_files: string (ex:"Load-Effic_Ba138_140uJ_80uW.*")
    #Build sanitized data file and get metadata and header
    new_file = data_files[:-1] + '_sanitized'
    data_files = glob.glob(data_files)
    if not glob.glob(new_file + '.txt'): #Check if the data file has already been sanitized
        new_file = sanitizeData(data_files, new_file)
    else:
        new_file += '.txt'
    meta_data, Header = getMetaData(data_files[0])
    Header = [i.replace('\n', '') for i in Header]
    return new_file, meta_data, Header #Returns - new_file: string (full filename), meta_data: list, Header: list

def sanitizeData(data_files, new_file): #Inputs - data_files: list of strings, new_file: string
    #Sanitize the raw data file into a nice new file
    with open(new_file + '.txt', 'w+') as new: #Write new file
        for i, data_file in enumerate(data_files): #Go through all data files
            with open(data_file, 'r') as file:
                for line in file:
                    Words = line.split("\t")
                    if Words[0].isdigit(): #Only write lines with digits for the first character
                        new.write(line)
    new_file_full = glob.glob(new_file + '.txt')[0]
    return  new_file_full #Returns - new_file_full: string (full filename)

def getMetaData(data_file): #Inputs - data_file: string (full filename)
    #Pull the metadata from the main data file
    Metadata = ''
    Header = ''
    with open(data_file, 'r') as file:
        for line in file: #Check each line, looking for tag Metadata
            Words = line.split(",")
            # print(Words)
            if Words[0] == 'Metadata':
                Metadata = Words[1:] #Skip the Metadata tag part of the line
                break
        for line in file: #Check each line, looking for the tag h for Header
            Words = line.split("\t")
            if Words[0][0] == 'h':
                Header = Words[1:]
                Header.insert(0, Words[0][1:]) #Take out the h tag
                break
    return Metadata, Header #Returns - Metadata: list, Header: list

def getFileNameData(data_file): #Inputs - data_file: string (full filename)
    #Pull data from filename
    Isotope, PulseEnergy, LaserPower = re.split("_", re.sub("[^0-9|^_]", "", new_file))[1:4]
    Isotope, PulseEnergy, LaserPower = int(Isotope), int(PulseEnergy), int(LaserPower)
    return  Isotope, PulseEnergy, LaserPower #Returns - Isotope: int, PulseEnergy: int, LaserPower: int

def getSanitizedData(data_file): #Inputs - data_file: string (full filename)
    #Extract all data from the sanitized data file
    raw_data = np.asarray(pd.read_csv(data_file, delimiter = "[\],\[\t]", engine='python',header=None))
    return raw_data #Returns - raw_data: np 2D array

def efficiencyWriteString(PulseEnergy, LaserPower, CoolingFreq, RepumpFreq, \
                        WindowStart, WindowWidth, DataPoints, LoadingEff): 
    #Inputs - PulseEnergy: int, LaserPower: int, CoolingFreq: float, RepumpFreq: float, 
    #WindowStart: float, WindowWidth: float, DataPoints: int, LoadingEff: float
    #Create nicely formatted string for putting parameter data in plot
    WriteString = f"""
PulseEnergy: {PulseEnergy}
553nmPower: {LaserPower}
CoolingFreq: {CoolingFreq}
RepumpFreq: {RepumpFreq}
WindowStart:{WindowStart}
WindowWidth:{WindowWidth}
DataPoints:{raw_data.shape[0]}
Loading Efficiency:{LoadingEff*100:0.1f}%"""
    return WriteString #Returns - WriteString: f-string

def analyzeLoading(raw_data, Verbose=False): #Inputs - raw_data: np 2D array, Verbose: bool
    #Analyze the loading data to find the number of ions trapped each time, and the loading efficiency
    BG = raw_data[:,5]
    BGStd = raw_data[:,6]
    IonCountsRaw = raw_data[:,2]
    IonCounts = IonCountsRaw.copy() - BG #Total Brightness of the ions
    #Find the average brightness of a single ion
    SingleAve, SingleStd = findSingleIonAverage(IonCounts, np.mean(BGStd), Verbose=Verbose)
    ThresholdSingle = SingleAve/6 #Using SingleStd sometimes fails, because of drifting signal
    IonsTrapped = countIons(IonCounts, SingleAve, ThresholdSingle, BGStd, Verbose=Verbose) #Count the number of ions trapped
    Attempts = raw_data[:,1] 
    IonsTrappedExpanded = expandLoadingData(Attempts, IonsTrapped) #Add in data points when none trapped
    TempIonsTrappedExpanded = IonsTrappedExpanded.copy()
    TempIonsTrappedExpanded[TempIonsTrappedExpanded < 0] = 0 #Only Ba138, for loading efficiency
    LoadingEfficiency = np.mean(TempIonsTrappedExpanded)
    if Verbose:
        print(f'Inside analyzeLoading() function')
        print(f'BGave:{np.mean(BG)}, BGStdave:{np.mean(BGStd)}')
        print(f'IonCounts:{IonCounts}')
        print(f'IonsTrapped:{IonsTrapped.astype(np.int32)}')
        print(f'IonsTrappedExpanded:{IonsTrappedExpanded.astype(np.int32)}')
        print(f'LoadingEfficiency:{LoadingEfficiency}')
        print('\n')
    return IonsTrappedExpanded, LoadingEfficiency #Returns - IonsTrappedExpanded: np 1D array, LoadingEfficiency: float

def findSingleIonAverage(Counts, BGStd, ThresholdFactor=10, Verbose=False): 
    #Inputs - Counts: np 1D array, BGStd: float, Verbose: bool
    #Figure out how bright an average single ion is
    Threshold = ThresholdFactor*BGStd
    CountsNonBG = Counts[Counts > Threshold] #Get instances where counts higher than BG
    SmallestNonBG = CountsNonBG[np.argpartition(CountsNonBG, 3)[:3]]
    SmallestNonBGAve = np.mean(SmallestNonBG) #Average the smallest 3 non-BG instances
    CountsSingle = CountsNonBG[CountsNonBG < 1.5*SmallestNonBGAve] #Get only the single ion instances
    SingleAve = np.mean(CountsSingle)
    SingleStd = np.std(CountsSingle)
    if Verbose:
        print(f'Inside findSingleIonAverage() function')
        print(f'Threshold:{Threshold}')
        print(f'SmallestNonBG:{SmallestNonBG}')
        print(f'CountsSingle:{CountsSingle}')
        print(f'SingleAve:{SingleAve}')
        print(f'SingleStd:{SingleStd}')
        fig = plt.figure(figsize=(20,10)) #Plot the ion counts and brightness of numbers of ions
        ax1 = fig.add_subplot(111)
        ax1.plot(CountsNonBG)
        for i in range(1, 5): 
            IonsCountBase = i*SingleAve
            ax1.plot([0, CountsNonBG.size], [IonsCountBase, IonsCountBase])
        print('\n')
    return SingleAve, SingleStd #Returns - SingleAve: float, SingleStd: float

def countIons(Counts, SingleAve, ThresholdSingle, BGstd, Verbose=False):
    #Inputs - Counts: np 2D array, SingleAve: float, ThresholdSingle: float, BGstd: float, Verbose: bool
    #Count number of ions trapped for each instance
    IonsTrapped = np.zeros(Counts.shape)
    NumberNotFoundYet = np.ones(Counts.shape) #Keep track of whether an instance has been assigned a number of ions yet
    for i in range(1, 5):
        #Check if assigned number yet, and if counts are within the range for this number of ions
        Condition = np.array([(SingleAve - ThresholdSingle)*i < Counts, (SingleAve + ThresholdSingle)*i > Counts, NumberNotFoundYet])
        ConditionTrue = np.all(Condition, axis=0)
        IonsTrappedTemp = ConditionTrue*i
        NumberNotFoundYet -= ConditionTrue #If assigned, mark it assigned
        IonsTrapped += IonsTrappedTemp
        if Verbose:
            if i == 1:
                print(f'Inside countIons() function')
            print(f'i: {i}')
            print(f'IonsTrappedTemp:  {IonsTrappedTemp.astype(np.int32)}')
            print(f'NumberNotFoundYet:{NumberNotFoundYet.astype(np.int32)}')
    ConditionOtherIsotope = np.array([Counts > 3*BGstd, NumberNotFoundYet])#Check if assigned number yet, and if over background
    ConditionOtherTrue = np.all(ConditionOtherIsotope, axis=0)
    IonsTrappedTemp = ConditionOtherTrue*-1
    NumberNotFoundYet -= ConditionTrue #If assigned, mark it assigned
    IonsTrapped += IonsTrappedTemp
    if Verbose:
        print(f'IonsTrappedTemp:  {IonsTrappedTemp.astype(np.int32)}')
        print(f'NumberNotFoundYet:{NumberNotFoundYet.astype(np.int32)}')
        print(f'\n')
    return IonsTrapped #Returns - IonsTrapped: np 1D array

def expandLoadingData(Attempts, IonsTrapped): #Inputs - Attempts: np 1D array, IonsTrapped: np 1D array
    #Expand the loading data to include instances where no ions trapped
    IonsTrappedExpanded = np.array([])
    for i, Attemptss in enumerate(Attempts):
        IonsTrappedExpanded = np.concatenate((IonsTrappedExpanded, np.zeros(int(Attemptss-1)))) #Add zeros in trapped array
        IonsTrappedExpanded = np.append(IonsTrappedExpanded, IonsTrapped[i])
    return IonsTrappedExpanded #Returns - IonsTrappedExpanded: np 1D array

def plotEfficiencyVsPulse(raw_data, WriteString, PositionText): #Inputs - raw_data: np 2D array, WriteString: str, PositionText: list
    #Plot the loading efficiency and neutral fluorescence vs pulse number
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    
    ax1.plot(IonsTrapped, linestyle='', marker='o') #Plot number of ions trapped
    ax1.text(PositionText[0], PositionText[1], WriteString, fontsize=12) #Write parameter data in figure
    ax2.plot(np.cumsum(raw_data[:,1])-1, raw_data[:,4], color='r') #Plot neutral fluorescence
    
    ax2.set_ylabel(f'Counts (a.u.)', fontsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel(f'Pulse Number',fontsize=30)
    ax1.set_ylabel('Num. Trapped (-1 rep. other iso.)', fontsize=30)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    return fig, ax1, ax2




















data_files = "Load-Effic_Ba138_140uJ_80uW_No405*"
new_file, meta_data, Header = prepareData(data_files)

IonizationFreq, CoolingFreq, RepumpFreq, WindowStart, WindowWidth = [float(i.split(':')[1]) for i in meta_data]
Isotope, PulseEnergy, LaserPower = getFileNameData(new_file)

raw_data = getSanitizedData(new_file)
IonsTrapped, LoadingEff = analyzeLoading(raw_data, Verbose=True)
# print(IonsTrapped)
WriteString = efficiencyWriteString(PulseEnergy, LaserPower, CoolingFreq, RepumpFreq, WindowStart, WindowWidth, raw_data.shape[0], LoadingEff)

fig, ax1, ax2 = plotEfficiencyVsPulse(raw_data, WriteString, [0, -1])





data_files = "Load-Effic_Ba138_140uJ_80uW.*"
new_file, meta_data, Header = prepareData(data_files)

IonizationFreq, CoolingFreq, RepumpFreq, WindowStart, WindowWidth = [float(i.split(':')[1]) for i in meta_data]
Isotope, PulseEnergy, LaserPower = getFileNameData(new_file)

raw_data = getSanitizedData(new_file)
IonsTrapped, LoadingEff = analyzeLoading(raw_data, Verbose=True)
# print(IonsTrapped)
WriteString = efficiencyWriteString(PulseEnergy, LaserPower, CoolingFreq, RepumpFreq, WindowStart, WindowWidth, raw_data.shape[0], LoadingEff)

fig, ax1, ax2 = plotEfficiencyVsPulse(raw_data, WriteString, [20, 2.5])





data_files = "Load-Effic_Ba138_140uJ_8uW.*"
new_file, meta_data, Header = prepareData(data_files)

IonizationFreq, CoolingFreq, RepumpFreq, WindowStart, WindowWidth = [float(i.split(':')[1]) for i in meta_data]
Isotope, PulseEnergy, LaserPower = getFileNameData(new_file)

raw_data = getSanitizedData(new_file)
IonsTrapped, LoadingEff = analyzeLoading(raw_data, Verbose=True)
# print(IonsTrapped)
WriteString = efficiencyWriteString(PulseEnergy, LaserPower, CoolingFreq, RepumpFreq, WindowStart, WindowWidth, raw_data.shape[0], LoadingEff)

fig, ax1, ax2 = plotEfficiencyVsPulse(raw_data, WriteString, [0, 2.1])
