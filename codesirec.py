#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:02:00 2020

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
# https://towardsdatascience.com/infectious-disease-modelling-beyond-the-basic-sir-model-216369c584c4
# https://idmod.org/docs/malaria/model-sir.html

# importing packages
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
import sys


#mitigation factor is q
#psd is the ratio of extreme social distancers (between 0 and 1)
#pss is the ratio of super spreaders (between 0 and 1)
#avoid pss+psd > 1

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
N = 40E6
# total time (days)
tottime = 2000
# initial percent removed (immune)
pr = 0.0
# initial percent infected
pi = 0.0001
# initial percent susceptible
ps = 1-pr-pi
# percent extreme isolators (social distance)
psd = 0.2
# percent super spreaders
pss = 0.2
# day restrictions lifted
day_lifted = 100
# rate of immunity loss (ratio of recovered people who lose immunity each day)
r = 0.001
# death rate
d = 0.06
# rate at which people die (1/ time to die after infection)
rho = 1/20 #20 days to die?

niter = int(math.ceil(tottime/dt))
t = np.arange(0, tottime, dt)   
S = np.zeros(niter)
R = np.zeros(niter)
I = np.zeros(niter)
D = np.zeros(niter)
    
     
S[0] = ps*N #total population
I[0] = pi*N
R[0] = pr*N
D[0] = d*R[0]

    

for j in range(niter-1):
    dSdt=-beta/N*S[j]*I[j] + r*R[j]
    S[j+1] = S[j] + dt*dSdt
    dIdt=beta/N*S[j]*I[j]-(1-d)*gamma*I[j]-d*rho*I[j]
    dRdt=(1-d)*gamma*I[j] - r*R[j]
    dDdt=d*rho*I[j]
    S[j+1] = S[j] + dt*dSdt
    I[j+1] = I[j] + dt*dIdt
    R[j+1] = R[j] + dt*dRdt
    D[j+1] = D[j] + dt*dDdt
plt.plot(t, I)  
