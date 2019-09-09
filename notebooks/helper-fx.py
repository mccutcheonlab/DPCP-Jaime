# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 12:19:22 2019

@author: jmc010
"""

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

# import JM_general_functions as jmf
# import JM_custom_figs as jmfig

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

plt.style.use('seaborn-muted')

import os
import timeit

tic = timeit.default_timer()

datafolder = '..\\data\\'

class Rat(object):
    
    nRats = 0
    nSessions = 0
    
    def __init__(self, data):      
        self.rat = data
        self.sessions = {}
        
        Rat.nRats += 1
                
    def loadsession(self, data, header):
        self.session = 's'+str(data[3]) #should reference column of data with session number
        self.sessions[self.session] = Session(data, header, self.rat, self.session)
       
        Rat.nSessions += 1
        
class Session(object):
    
    def __init__(self, data, header, rat, session):
        self.hrow = {}
        for idx, col in enumerate(header):
            self.hrow[col] = data[idx]
        self.medfile = datafolder + self.hrow['medfile']
        self.distractorpresent = self.hrow['distractorpresent']
        self.condition = self.hrow['condition']
        self.rat = str(rat)
        self.session = session
        
    def extractlicks(self):
        self.licks, self.offset = medfilereader(self.medfile,
                                                    varsToExtract = ['e', 'f'],
                                                    remove_var_header = True)
        if len(self.licks) == len(self.offset)+1:
            self.licks = self.licks[:-1]
            
    def extractdistractors(self):
        self.distractors_calc = calcDistractors(self.lickData['licks'])
        
        if self.distractorpresent == 1:
            self.distractors, self.distractortype, self.distractedOrNot = jmf.medfilereader(self.medfile,
                                                                                        varsToExtract = ['i', 'j', 'k'],
                                                                                        remove_var_header = True)
            if len(self.distractors_calc) != len(self.distractors):
                print('Calculated distractors (n={}) do not match TTLs from data file (n={})!'.format(len(self.distractors_calc), len(self.distractors)))
            self.distractorstatus = 'Distractor'
        else:
            self.distractors = self.distractors_calc
            self.distractorstatus = 'Simulated distractor'

def medfilereader(filename, varsToExtract = 'all',
                  sessionToExtract = 1,
                  verbose = False,
                  remove_var_header = False):
    if varsToExtract == 'all':
        numVarsToExtract = np.arange(0,26)
    else:
        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
    
    f = open(filename, 'r')
    f.seek(0)
    filerows = f.readlines()[8:]
    datarows = [isnumeric(x) for x in filerows]
    matches = [i for i,x in enumerate(datarows) if x == 0.3]
    if sessionToExtract > len(matches):
        print('Session ' + str(sessionToExtract) + ' does not exist.')
    if verbose == True:
        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
        print('Analyzing session ' + str(sessionToExtract))
    
    varstart = matches[sessionToExtract - 1]
    medvars = [[] for n in range(26)]
    
    k = int(varstart + 27)
    for i in range(26):
        medvarsN = int(datarows[varstart + i + 1])
        
        medvars[i] = datarows[k:k + int(medvarsN)]
        k = k + medvarsN
        
    if remove_var_header == True:
        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
    else:
        varsToReturn = [medvars[i] for i in numVarsToExtract]

    if np.shape(varsToReturn)[0] == 1:
        varsToReturn = varsToReturn[0]
    return varsToReturn

def metafilereader(filename):
    
    f = open(filename, 'r')
    f.seek(0)
    header = f.readlines()[0]
    f.seek(0)
    filerows = f.readlines()[1:]
    
    tablerows = []
    
    for i in filerows:
        tablerows.append(i.split('\t'))
        
    header = header.split('\t')
    # need to find a way to strip end of line \n from last column - work-around is to add extra dummy column at end of metafile
    return tablerows, header

def isnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')
def lickCalc(licks, offset = [], burstThreshold = 0.25, runThreshold = 10, 
             binsize=60, histDensity = False, adjustforlonglicks='none'):
    # makes dictionary of data relating to licks and bursts
    if type(licks) != np.ndarray or type(offset) != np.ndarray:
        try:
            licks = np.array(licks)
            offset = np.array(offset)
        except:
            print('Licks and offsets need to be arrays and unable to easily convert.')
            return

    lickData = {}

    if len(offset) > 0:        
        lickData['licklength'] = offset - licks[:len(offset)]
        lickData['longlicks'] = [x for x in lickData['licklength'] if x > 0.3]
    else:
        lickData['licklength'] = []
        lickData['longlicks'] = []
    
    if adjustforlonglicks != 'none':
        if len(lickData['longlicks']) == 0:
            print('No long licks to adjust for.')
        else:
            lickData['median_ll'] = np.median(lickData['licklength'])
            lickData['licks_adj'] = int(np.sum(lickData['licklength'])/lickData['median_ll'])
            if adjustforlonglicks == 'interpolate':
                licks_new = []
                for l, off in zip(licks, offset):
                    x = l
                    while x < off - lickData['median_ll']:
                        licks_new.append(x)
                        x = x + lickData['median_ll']
                licks = licks_new
        
    lickData['licks'] = licks
    lickData['ilis'] = np.diff(np.concatenate([[0], licks]))
    lickData['shilis'] = [x for x in lickData['ilis'] if x < burstThreshold]
    lickData['freq'] = 1/np.mean([x for x in lickData['ilis'] if x < burstThreshold])
    lickData['total'] = len(licks)
    
    # Calculates start, end, number of licks and time for each BURST 
    lickData['bStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]  
    lickData['bInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]
    lickData['bEnd'] = [lickData['licks'][i-1] for i in lickData['bInd'][1:]]
    lickData['bEnd'].append(lickData['licks'][-1])

    lickData['bLicks'] = np.diff(lickData['bInd'] + [len(lickData['licks'])])    
    lickData['bTime'] = np.subtract(lickData['bEnd'], lickData['bStart'])
    lickData['bNum'] = len(lickData['bStart'])
    if lickData['bNum'] > 0:
        lickData['bMean'] = np.nanmean(lickData['bLicks'])
        lickData['bMean-first3'] = np.nanmean(lickData['bLicks'][:3])
    else:
        lickData['bMean'] = 0
        lickData['bMean-first3'] = 0
    
    lickData['bILIs'] = [x for x in lickData['ilis'] if x > burstThreshold]

    # Calculates start, end, number of licks and time for each RUN
    lickData['rStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]  
    lickData['rInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]
    lickData['rEnd'] = [lickData['licks'][i-1] for i in lickData['rInd'][1:]]
    lickData['rEnd'].append(lickData['licks'][-1])

    lickData['rLicks'] = np.diff(lickData['rInd'] + [len(lickData['licks'])])    
    lickData['rTime'] = np.subtract(lickData['rEnd'], lickData['rStart'])
    lickData['rNum'] = len(lickData['rStart'])

    lickData['rILIs'] = [x for x in lickData['ilis'] if x > runThreshold]
    try:
        lickData['hist'] = np.histogram(lickData['licks'][1:], 
                                    range=(0, 3600), bins=int((3600/binsize)),
                                          density=histDensity)[0]
    except TypeError:
        print('Problem making histograms of lick data')
        
    return lickData    

def distractedOrNot(distractors, licks, delay=1):   
    firstlick = []
    distractedArray = []

    for d in distractors:               
        try:
            firstlick.append([i-d for i in licks if (i > d)][0])
        except IndexError:
            firstlick.append(np.NaN)

    distractedArray = np.array([i>delay for i in firstlick], dtype=bool)
    
    if np.isnan(firstlick)[-1] == 1: 
        distractedArray[-1] = True
    
    return firstlick, distractedArray

def calcDistractors(licks):
#    disStart = [licks[idx-2] for idx, val in enumerate(licks) if (val - licks[idx-2]) < 1]
#    disEnd = [val for idx, val in enumerate(licks) if (val - licks[idx-2]) < 1]
#    
#    distractors = [val for idx, val in enumerate(disEnd) if disStart[idx] - disEnd[idx-1] > 1]
    d = []

    i=2
    while i < len(licks):
        if licks[i] - licks[i-2] < 1:
            d.append(licks[i])
            i += 1
            try:
                while licks[i] - licks[i-1] < 1:
                    i += 1
            except IndexError:
                pass
        else:
            i += 1
    
    distractors = d
    
    return distractors
        
metafile = '..\data\DPCP1Masterfile.txt'
#metafile = 'R:\\DA_and_Reward\\kp259\THPH1\\THPH1 Scripts_170616\\thph1-forMatPy.txt'
metafileData, metafileHeader = metafilereader(metafile)

exptsuffix = ''
includecol = 11

rats = {}

for i in metafileData:
    if int(i[includecol]) == 1:
        rowrat = str(i[1])
        if rowrat not in rats:
            rats[rowrat] = Rat(rowrat)
        rats[rowrat].loadsession(i, metafileHeader)
        
for i in rats:
    for j in rats[i].sessions:
#        print('Analysing rat ' + i + ' in session ' + j)
        x = rats[i].sessions[j]
        
        x.extractlicks()


        x.lickData = lickCalc(x.licks,
                          offset = x.offset,
                          burstThreshold = 0.50)
        
        x.extractdistractors()
        x.firstlick, x.distractedArray = distractedOrNot(x.distractors, x.lickData['licks'])
        
