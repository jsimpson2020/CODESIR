#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 15:06:55 2020

@author: jsimpson
"""


# uses the euler method to numerically solve coupled first order ODE's
# for the SIR ODE model using the euler method
# S = Susceptible, I = Infected, R = Removed (recovered/deceased/immune)
# dS/dt = -beta*S*I/N, dI/dt = beta*S*I-gamma*I, dR/dt=gamma*I
# where S+I+R=constant=total population
# https://www.davidketcheson.info/2020/03/19/SIR_Estimating_parameters.html
# https://www.davidketcheson.info/2020/03/19/SIR_predictions.html
# basic assumption beta = 0.25 and gamma = 0.05

# importing packages
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
import sys

def mixed_quarantine_comparison(spacing, psd, q1, r):
    #mitigation factor is q
    #spacing is the spacing between trial runs in days
    #q1 and r are varying mitigation factors
    #q1 is the initial quarantine level (between 0 and 1)
    #r is the q value when quarantine is lifted (between 0 and 1)
    #psd is the ratio of people participating in quarantine (between 0 and 1)
    #this version turns has both groups quarantining, but one group stops
    #after j > some number of days determined by the spacing factor
    
    #initializing parameters
    # time step (day)
    dt = 0.01
    # average number of people that come within infection range of an infected individual
    # per day
    # ASSUMES NO SOCIAL DISTANCING OR MITIGATION
    beta = 0.25
    # remove probability per day = 1/(recovery time)
    # approximately 16 days (varies from 14 days to a month)
    gamma = 0.05
    
    # total population
    N = 1E6
    # total time (days)
    tottime = 365
    # initial percent removed (immune)
    pr = 0.0
    # initial percent infected
    pi = 0.01
    # initial percent susceptible
    ps = 1-pr-pi
    # percent isolators (social distance)
    #psd = 0.5
    
    niter = int(math.ceil(tottime/dt))
    t = np.arange(0, tottime, dt)   
    S = np.zeros(niter)
    S1 = np.zeros(niter)
    S2 = np.zeros(niter)
    I = np.zeros(niter)
    
     
    S[0] = ps*N
    S1[0]= psd*S[0]
    S2[0]= S[0]-S1[0]
    I[0] = pi*N
    sum = 0
    
#q is quarantine level for group that remains in quarantine after restrictions lifted
#k is quarantine level for group that leaves quarantine after restrictions lifted
    for n in range(8):
        if n <= 5:
            for j in range(niter-1):
                if j <= (spacing * (n+1) / dt):
                    q = q1
                    k = q1
                else:
                    q = q1
                    k = r
                dS1dt=-q*beta/N*S1[j]*I[j]
                dS2dt=-k*beta/N*S2[j]*I[j]
                dIdt=-dS1dt-dS2dt-gamma*I[j]
                S1[j+1] = S1[j] + dt*dS1dt 
                S2[j+1] = S2[j] + dt*dS2dt
                I[j+1] = I[j] + dt*dIdt  
        
            S=S1+S2   
            plt.plot(t, I)
        elif n == 6:
            for k in range(niter-1):
                dSdt=-(beta/N)*S[k]*I[k]
                dIdt=(beta/N)*S[k]*I[k]-gamma*I[k]
                sum = sum + I[k]*dt
                S[k+1] = S[k] + dt*dSdt 
                I[k+1] = I[k] + dt*dIdt
            plt.plot(t,I)
        else:
            for l in range(niter-1):
                dS1dt=-q1*beta/N*S1[l]*I[l]
                dS2dt=-beta/N*S2[l]*I[l]
                dIdt=-dS1dt-dS2dt-gamma*I[l]
                S1[l+1] = S1[l] + dt*dS1dt 
                S2[l+1] = S2[l] + dt*dS2dt
                I[l+1] = I[l] + dt*dIdt
            plt.plot(t,I)
    
    plt.legend([str(int(spacing)) + ' days of q = ' + str(float(q1)), 
                str(int(2 * spacing)) + ' days of q = ' + str(float(q1)), 
                str(int(3 * spacing)) + ' days of q = ' + str(float(q1)), 
                str(int(4 * spacing)) + ' days of q = ' + str(float(q1)),
                str(int(5 * spacing)) + ' days of q = ' + str(float(q1)),
                str(int(6 * spacing)) + ' days of q = ' + str(float(q1)),
                'no change', 'permanent q = ' + str(float(q1))],
                loc='upper right')
    plt.title('Mixed SIR model for removing restrictions with \n q = ' + str(float(q1))
              + ', \u03B2 = ' + str(beta) + ', \u03B3 = ' + str(gamma)
              + ', % isolators = '+ str(round(float(psd*100),3)))
    plt.xlabel('days elapsed since 1 percent of the population became infected')
    plt.ylabel('population')
   
    return plt.show()

#for m in range(10):
mixed_quarantine_comparison(21, 0.9, 0.5, 1)
