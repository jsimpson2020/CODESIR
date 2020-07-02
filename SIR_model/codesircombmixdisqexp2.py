#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 15:20:59 2020

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

def mixed_quarantine_comparison(day_lifted, spacing, psd, q1, r, q2):
    #mitigation factor is q
    #spacing is the spacing between trial runs in days
    #q1 and r are varying mitigation factors
    #q1 is the initial quarantine level for isolators (between 0 and 1)
    #r is the q value when quarantine is lifted (between 0 and 1)
    #q2 is the initial quarantine level for nonisolators (between 0 and 1)
    #psd is the ratio of people participating in quarantine (between 0 and 1)
    
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
    tottime = 250
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
    
    #First_change is the day that restrictions are lifted
    first_change = day_lifted / dt
    
#q is quarantine factor for coorperative group
#k is quarantine factor for group that refuses to go back to quarantine
    for n in range(6):
        if n <= 4:
            second_change = first_change + (spacing * n / dt)
            for j in range(niter-1):
                if j <= first_change:
                    q = q1
                    k = q2
                elif j <= second_change and j > first_change:
                    q = r
                    k = r
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
        else:
            for k in range(niter-1):
                dSdt=-(beta/N)*S[k]*I[k]
                dIdt=(beta/N)*S[k]*I[k]-gamma*I[k]
                sum = sum + I[k]*dt
                S[k+1] = S[k] + dt*dSdt 
                I[k+1] = I[k] + dt*dIdt
            plt.plot(t,I)
    plt.legend(['no quarantine lift',
                str(int(spacing)) + ' days relaxed', 
                str(int(2 * spacing)) + ' days relaxed', 
                str(int(3 * spacing)) + ' days relaxed', 
                str(int(4 * spacing)) + ' days relaxed',
                'no quarantine'], loc='upper right')
    plt.title('Mixed SIR model lifting restrictions after ' + str(int(day_lifted))
              + ' days of quarantine \n and replacing restrictions on ' + str(int(spacing))
              + ' day intervals (I only) \n q1 = ' + str(float(q1)) + ', r = ' + str(float(r))
              + ', q2 = ' + str(float(q2)) + ', \u03B2 = ' + str(beta) + ', \u03B3 = ' +
              str(gamma) + ', psd = ' + str(float(psd)*100) + '%')
    plt.xlabel('days elapsed since 1 percent of the population became infected')
    plt.ylabel('population')
   
    return plt.show()

#for m in range(10):
mixed_quarantine_comparison(77, 5, 0.9, 0.3, 1, 0.7)