# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:39:30 2017
Analysis for dpcp1 experiment

@author: James Rig
"""
import sys
sys.path.insert(0,'C:\\Github\\functions-and-figures\\')

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import scipy.optimize as opt
import scipy.stats as stats

import JM_general_functions as jmf
import JM_custom_figs as jmfig

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

import pandas as pd

plt.style.use('seaborn-muted')

import os
import timeit



tic = timeit.default_timer()

#userhome = os.path.expanduser('~')
datafolder = 'data\\'

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
    
    def __init__(self, sessionID, metafiledata, hrows, datafolder):
        self.sessionID = sessionID
        self.rat = metafiledata[hrows['rat']].replace('.', '-')
        self.session = metafiledata[hrows['session']]
        self.medfile = datafolder + metafiledata[hrows['medfile']]
        self.distractorpresent = metafiledata[hrows['distractorpresent']]
        self.condition = metafiledata[hrows['condition']]
 
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
            
    def calculate_pdp_prob(self):
        self.pdp_xdata = np.sort(self.firstlick)
        #self.pdp_xdata = np.log(np.sort(self.firstlick))
        self.pdp_ydata = [1-i/len(self.pdp_xdata) for i,val in enumerate(self.pdp_xdata)]

    def fit_singleexp(self):
        x0=[0.5]
        try:
            self.fit1=opt.curve_fit(singleexp, self.pdp_xdata, self.pdp_ydata, x0)
            self.sexp_alpha=self.fit1[0][0]
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                    self.pdp_ydata, singleexp(self.pdp_xdata,
                                     self.sexp_alpha))
            self.sexp_rsq=r_value**2
        except:
            print('could not fit single exp')
            self.sexp_alpha=0
            self.sexp_rsq=0
#  

    def fit_doubleexp(self):
        x0=np.array([0.5, 1, 1])
        try:
            self.fit=opt.curve_fit(doubleexp, self.pdp_xdata, self.pdp_ydata, x0)
            self.dexp_alpha=self.fit[0][0]
            self.dexp_beta=self.fit[0][1]
            self.dexp_tau=self.fit[0][2]
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                    self.pdp_ydata, doubleexp(self.pdp_xdata,
                                     self.dexp_alpha,
                                     self.dexp_beta,
                                     self.dexp_tau))
            self.dexp_rsq=r_value**2
            print('Parameters found')
        except:
            print('Optimal parameters not found for', self.rat, 'session', self.session)
            self.dexp_alpha=0
            self.dexp_beta=0
            self.dexp_tau=0
            self.dexp_rsq=0


def singleexp(t, alpha):
    return np.exp(-alpha*t)

def doubleexp(t, alpha, beta, tau):
    return alpha*(np.exp(-beta*t)) + (1-alpha)*np.exp(-tau*t)

def metafile2sessions(metafile, datafolder):
#    jmf.metafilemaker(xlfile, metafile, sheetname=sheetname, fileformat='txt')
    rows, header = jmf.metafilereader(metafile)
    
    hrows = {}
    for idx, field in enumerate(header):
        hrows[field] = idx
    
    sessions = {}
    
    for row in rows:
        sessionID = row[hrows['rat']].replace('.','-') + '_' + row[hrows['session']]
        sessions[sessionID] = Session(sessionID, row, hrows, datafolder)
    
    return sessions   


     
metafile = 'data\\DPCP1Masterfile.txt'
datafolder='data\\'

sessions = metafile2sessions(metafile, datafolder)

for session in sessions:
      
    x = sessions[session]
    x.extractlicks()
    x.lickData = jmf.lickCalc(x.licks,
                          offset = x.offset,
                          burstThreshold = 0.50)
    x.extractdistractors()
    x.firstlick, x.distractedArray = jmf.distractedOrNot(x.distractors, x.lickData['licks'])
    x.calculate_pdp_prob()
    x.fit_singleexp()
    x.fit_doubleexp()
    



df = pd.DataFrame([sessions[s].rat for s in sessions], columns=['rat'])
    
#exptsuffix = ''
#includecol = 11

#rats = {}
#
#for i in metafileData:
#    if int(i[includecol]) == 1:
#        rowrat = str(i[1])
#        if rowrat not in rats:
#            rats[rowrat] = Rat(rowrat)
#        rats[rowrat].loadsession(i, metafileHeader)
        
#for i in rats:
#    for j in rats[i].sessions:
##        print('Analysing rat ' + i + ' in session ' + j)
#        x = rats[i].sessions[j]
#        
#        x.extractlicks()
#
#
#        x.lickData = jmf.lickCalc(x.licks,
#                          offset = x.offset,
#                          burstThreshold = 0.50)
#        
#        x.extractdistractors()
#        x.firstlick, x.distractedArray = jmf.distractedOrNot(x.distractors, x.lickData['licks'])
#        
#        x.calculate_pdp_prob()
#        x.fit_doubleexp()
#        
#        
#i = 'dpcp1.16'
#j = 's4'
#
#session = rats[i].sessions[j]
#
#j='s3'
#pd.DataFrame([s.rat for s in rats])

#for i in rats:
#    j = 's3'
#    x = rats[i].sessions[j]
#    cumsumBL.append(x.firstlick)
#    if x.condition == 'SAL':
#        cumsumBLsal.append(x.firstlick)
#    elif x.condition == 'PCP':
#        cumsumBLpcp.append(x.firstlick)
#        
#    j = 's4'
#    x = rats[i].sessions[j]
#    cumsumDIS.append(x.firstlick)
#    if x.condition == 'SAL':
#        cumsumDISsal.append(x.firstlick)
#    elif x.condition == 'PCP':
#        cumsumDISpcp.append(x.firstlick)
#    
#
#cumsumBL = []
#cumsumBLsal = []
#cumsumBLpcp = []
#
#cumsumDIS= []
#cumsumDISsal = []
#cumsumDISpcp = []

#for i in rats:
#    j = 's3'
#    x = rats[i].sessions[j]
#    cumsumBL.append(x.firstlick)
#    if x.condition == 'SAL':
#        cumsumBLsal.append(x.firstlick)
#    elif x.condition == 'PCP':
#        cumsumBLpcp.append(x.firstlick)
#        
#    j = 's4'
#    x = rats[i].sessions[j]
#    cumsumDIS.append(x.firstlick)
#    if x.condition == 'SAL':
#        cumsumDISsal.append(x.firstlick)
#    elif x.condition == 'PCP':
#        cumsumDISpcp.append(x.firstlick)
#    
#groupFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
#gs1 = gridspec.GridSpec(5, 2)
#gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
#plt.suptitle('Group data')
#    
#ax = plt.subplot(gs1[0, 0])
#for x in cumsumBL:
#    jmfig.cumulativelickFig(ax, x, color='grey')
#avg = [item for rat in cumsumBL for item in rat]
#jmfig.cumulativelickFig(ax, avg, color='k')
#ax.set_title('Baseline')
#    
#    
#ax = plt.subplot(gs1[0, 1])
#for x in cumsumDIS:
#    jmfig.cumulativelickFig(ax, x, color='grey')
#avg = [item for rat in cumsumDIS for item in rat]
#jmfig.cumulativelickFig(ax, avg, color='k')
#ax.set_title('Distraction Day')
#
#ax = plt.subplot(gs1[1, 0])
#for x in cumsumDISsal:
#    jmfig.cumulativelickFig(ax, x, color='grey')
#avgDISsal = [item for rat in cumsumDISsal for item in rat]
#jmfig.cumulativelickFig(ax, avgDISsal, color='k')
#ax.set_title('Distraction Day - Saline')
#    
#    
#ax = plt.subplot(gs1[1, 1])
#for x in cumsumDISpcp:
#    jmfig.cumulativelickFig(ax, x, color='grey')
#avgDISpcp = [item for rat in cumsumDISpcp for item in rat]
#jmfig.cumulativelickFig(ax, avgDISpcp, color='k')
#ax.set_title('Distraction Day - PCP')
#    
##ks_2samp(avgDISsal, avgDISpcp)
#        
#
