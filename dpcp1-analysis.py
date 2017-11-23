# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:39:30 2017
Analysis for dpcp1 experiment

@author: James Rig
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

import JM_general_functions as jmf
import JM_custom_figs as jmfig

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

plt.style.use('seaborn-muted')

import os
import timeit

tic = timeit.default_timer()

datafolder = 'R:\\DA_and_Reward\\jem64\\1704_DPCP\\DPCP_Data\\'

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
        self.licks, self.offset = jmf.medfilereader(self.medfile,
                                                    varsToExtract = ['e', 'f'],
                                                    remove_var_header = True)
        if len(self.licks) == len(self.offset)+1:
            self.licks = self.licks[:-1]
            
    def extractdistractors(self):
        self.distractors_calc = jmf.calcDistractors(self.lickData['licks'])
        
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
        
metafile = 'R:\\DA_and_Reward\\jem64\\1704_DPCP\\DPCP1Masterfile.txt'
#metafile = 'R:\\DA_and_Reward\\kp259\THPH1\\THPH1 Scripts_170616\\thph1-forMatPy.txt'
metafileData, metafileHeader = jmf.metafilereader(metafile)

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


        x.lickData = jmf.lickCalc(x.licks,
                          offset = x.offset,
                          burstThreshold = 0.50)
        
        x.extractdistractors()
        x.firstlick, x.distractedArray = jmf.distractedOrNot(x.distractors, x.lickData['licks'])
        
i = 'dpcp1.1'
j = 's2'

x = rats[i].sessions[j]

cumsumBL = []
cumsumBLsal = []
cumsumBLpcp = []

cumsumDIS= []
cumsumDISsal = []
cumsumDISpcp = []

for i in rats:
    j = 's3'
    x = rats[i].sessions[j]
    cumsumBL.append(x.firstlick)
    if x.condition == 'SAL':
        cumsumBLsal.append(x.firstlick)
    elif x.condition == 'PCP':
        cumsumBLpcp.append(x.firstlick)
        
    j = 's4'
    x = rats[i].sessions[j]
    cumsumDIS.append(x.firstlick)
    if x.condition == 'SAL':
        cumsumDISsal.append(x.firstlick)
    elif x.condition == 'PCP':
        cumsumDISpcp.append(x.firstlick)
    
groupFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
gs1 = gridspec.GridSpec(5, 2)
gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
plt.suptitle('Group data')
    
ax = plt.subplot(gs1[0, 0])
for x in cumsumBL:
    jmfig.cumulativelickFig(ax, x, color='grey')
avg = [item for rat in cumsumBL for item in rat]
jmfig.cumulativelickFig(ax, avg, color='k')
ax.set_title('Baseline')
    
    
ax = plt.subplot(gs1[0, 1])
for x in cumsumDIS:
    jmfig.cumulativelickFig(ax, x, color='grey')
avg = [item for rat in cumsumDIS for item in rat]
jmfig.cumulativelickFig(ax, avg, color='k')
ax.set_title('Distraction Day')

ax = plt.subplot(gs1[1, 0])
for x in cumsumDISsal:
    jmfig.cumulativelickFig(ax, x, color='grey')
avgDISsal = [item for rat in cumsumDISsal for item in rat]
jmfig.cumulativelickFig(ax, avgDISsal, color='k')
ax.set_title('Distraction Day - Saline')
    
    
ax = plt.subplot(gs1[1, 1])
for x in cumsumDISpcp:
    jmfig.cumulativelickFig(ax, x, color='grey')
avgDISpcp = [item for rat in cumsumDISpcp for item in rat]
jmfig.cumulativelickFig(ax, avgDISpcp, color='k')
ax.set_title('Distraction Day - PCP')
    
ks_2samp(avgDISsal, avgDISpcp)
        

