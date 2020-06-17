#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Jun 16 13:55:06 2020

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
import random

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
first_change = float(input("Enter day quarantine lifted (integer between 1 and 100): ")) / dt
#spacing is the spacing between trial runs in days
spacing = float(input("Enter comparison spacing (integer between 1 and 10): "))
#q1, q2, and r are varying mitigation factors
#q1 is the initial quarantine level (between 0 and 1)
#q2 is the return to quarantine levle (between 0 and 1)
#r is the q value when quarantine is lifted (between 0 and 1)
q1 = float(input("Enter 1st quarantine value (between 0 and 1): "))
q2 = float(input("Enter 2nd quarantine value (between 0 and 1): "))
r = float(input("Enter relaxation q value (bigger than quarantine values up to 1): "))

for n in range(5):
    second_change = first_change + (spacing * n / dt)
    qlist = []
    
    for i in range(niter-1):
        if i <= first_change: 
            q = q1
        elif i <= second_change and i > first_change:
            q = 1
        else:
            q = q2
        qlist.append(q)

    for j in range(niter-1):
        dSdt=-(qlist[j]*beta/N)*S[j]*I[j]
        dIdt=(qlist[j]*beta/N)*S[j]*I[j]-gamma*I[j]
        sum = sum + I[j]*dt
        S[j+1] = S[j] + dt*dSdt 
        I[j+1] = I[j] + dt*dIdt

    # plotting
    plt.plot(t, I)
    plt.legend(['no explicit quarantine lift', 
                str(int(spacing)) + ' days relaxed', 
                str(int(2 * spacing)) + ' days relaxed', 
                str(int(3 * spacing)) + ' days relaxed', 
                str(int(4 * spacing)) + ' days relaxed'])
    plt.title('SIR ode model for lifting and replacing quarantines')
    plt.xlabel('days elapsed since 1 percent of the population became infected')
    plt.ylabel('population')

plt.show()
