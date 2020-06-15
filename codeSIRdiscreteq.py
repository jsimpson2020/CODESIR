#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:27:56 2020

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
N = 7.7E9
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


#mitigation factor
qlist = []

for i in range(niter-1):
    #multiply i value by 100 to get the day value
    #Add conditionals if necessary
    if i <=2500: 
        q=0.6
    elif i<=7500 and i>2500:
        q=0.8
    else:
        q=1
    qlist.append(q)

for j in range(niter-1):
    dSdt=-(qlist[j]*beta/N)*S[j]*I[j]
    dIdt=(qlist[j]*beta/N)*S[j]*I[j]-gamma*I[j]
    sum = sum + I[j]*dt
    S[j+1] = S[j] + dt*dSdt 
    I[j+1] = I[j] + dt*dIdt  

   
       
R = N-S-I

# plotting
fig = plt.figure()
plt.plot(t, S, 'k', label = 'susceptible')
plt.plot(t, I, 'm', label = 'infected')
plt.plot(t, R, 'b', label = 'recovered')
plt.plot(t, R+S+I, 'y', label = 'total')
plt.gca().legend(('susceptible','infected','recovered','total'))
plt.title('simple SIR ode model with beta = 0.25, gamma = 0.05')
plt.xlabel('days elapsed since 1 percent of the population became infected')
plt.ylabel('population')
plt.show()
