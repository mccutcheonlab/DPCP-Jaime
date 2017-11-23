# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 08:20:30 2017

@author: James Rig
"""

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(6, 2))

#fig, ax = plt.figure(figsize=(6, 2.5), dpi=600)
#gs1 = gridspec.GridSpec(1, 2)
#gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
    
for x in cumsumBL:
    jmfig.cumulativelickFig(ax[0], x, color='grey')
avg = [item for rat in cumsumBL for item in rat]
jmfig.cumulativelickFig(ax[0], avg, color='k')
ax[0].set_title('Baseline')
    

for x in cumsumDIS:
    jmfig.cumulativelickFig(ax[1], x, color='grey')
avg = [item for rat in cumsumDIS for item in rat]
jmfig.cumulativelickFig(ax[1], avg, color='blue')
ax[1].set_title('Distraction Day')

ax[0].set_ylabel('Cumulative probability')

fig.text(0.50, -0.15, 'Lick pauses (s)', ha='center')
plt.savefig('R:\\DA_and_Reward\\jem64\\1704_DPCP\\figures\\lick-pauses.eps')

