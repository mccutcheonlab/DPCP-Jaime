# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:20:26 2019

@author: admin
"""

import scipy.optimize as opt
import scipy.stats as stats

session = rats[i].sessions[j]



def calculate_pdp_prob(pdps):
    
    
    xdata = np.sort(pdps)
    ydata = [1-i/len(xdata) for i,val in enumerate(xdata)]
    
    return xdata, ydata

def doubleexp(t, alpha, beta, tau):
    return alpha*(np.exp(-beta*t)) + (1-alpha)*np.exp(-tau*t)

pdps = session.firstlick

xdata, ydata = calculate_pdp_prob(session.firstlick)

f, ax = plt.subplots()
ax.scatter(xdata, ydata, color='none', edgecolors='grey')


#alpha=0.1
#beta=10
#tau=5
#
#
#ax.plot(xdata, doubleexp(xdata, alpha, beta, tau), color='g')

x0=np.array([0.5, 1, 1])
fit=opt.curve_fit(doubleexp, xdata, ydata, x0)

alpha=fit[0][0]
beta=fit[0][1]
tau=fit[0][2]

print(alpha, beta, tau)

ax.plot(xdata, doubleexp(xdata, alpha, beta, tau), color='orange')
ax.set_xscale('log')


slope, intercept, r_value, p_value, std_err = stats.linregress(ydata, doubleexp(xdata, alpha, beta, tau))

print(r_value**2)



def fit_weibull(xdata, ydata):
    x0=np.array([0.1, 1])
    fit=opt.curve_fit(weib_davis, xdata, ydata, x0)
    alpha=fit[0][0]
    beta=fit[0][1]
    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata, weib_davis(xdata, alpha, beta))
    r_squared=r_value**2
    
    return alpha, beta, r_squared

