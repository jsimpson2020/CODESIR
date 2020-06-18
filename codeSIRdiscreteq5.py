#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:55:50 2020

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

def quarantine_comparison(spacing, q1, r):
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
    N = 328.2E6
    # total time (days)
    tottime = 200
    # initial percent removed (immune)
    pr = 0.0
    # initial percent infected
    pi = 0.01
    # initial percent susceptible
    ps = 1-pr-pi
    
    niter = int(math.ceil(tottime/dt))
    t = np.arange(0, tottime, dt)   
    S = np.zeros(niter)
    I = np.zeros(niter)
    
     
    S[0] = ps*N
    I[0] = pi*N
    sum = 0
    
    #mitigation factor is q
    #First_change is the day that restrictions are lifted
    #spacing is the spacing between trial runs in days
    #q1, q2, and r are varying mitigation factors
    #q1 is the initial quarantine level (between 0 and 1)
    #q2 is the return to quarantine levle (between 0 and 1)
    #r is the q value when quarantine is lifted (between 0 and 1)

    for n in range(8):
        if n <= 5:
            for j in range(niter-1):
                if j <= (spacing * (n+1) / dt):
                    q = q1
                else:
                    q = r
                dSdt=-(q*beta/N)*S[j]*I[j]
                dIdt=(q*beta/N)*S[j]*I[j]-gamma*I[j]
                sum = sum + I[j]*dt
                S[j+1] = S[j] + dt*dSdt 
                I[j+1] = I[j] + dt*dIdt
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
                dSdt=-(q1*beta/N)*S[l]*I[l]
                dIdt=(q1*beta/N)*S[l]*I[l]-gamma*I[l]
                sum = sum + I[l]*dt
                S[l+1] = S[l] + dt*dSdt 
                I[l+1] = I[l] + dt*dIdt
            plt.plot(t,I)
    
    plt.legend([str(int(spacing)) + ' days of q = ' + str(float(q1)), 
                str(int(2 * spacing)) + ' days of q = ' + str(float(q1)), 
                str(int(3 * spacing)) + ' days of q = ' + str(float(q1)), 
                str(int(4 * spacing)) + ' days of q = ' + str(float(q1)),
                str(int(5 * spacing)) + ' days of q = ' + str(float(q1)),
                str(int(6 * spacing)) + ' days of q = ' + str(float(q1)),
                'no change', 'permanent q = ' + str(float(q1))])
    plt.title(' SIR model for removing restrictions with \n q = ' + str(float(q1))
              + ', \u03B2 = ' + str(beta) + ', \u03B3 = ' + str(gamma))
    plt.xlabel('days elapsed since 1 percent of the population became infected')
    plt.ylabel('population')
   
    return plt.show()

for m in range(10):
    quarantine_comparison(14, (m+6)/20, 1)
