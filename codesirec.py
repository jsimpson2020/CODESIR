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
tottime = 730
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
S1 = np.zeros(niter)
R1 = np.zeros(niter)
I1 = np.zeros(niter)
D1 = np.zeros(niter)
  
S2 = np.zeros(niter)
R2 = np.zeros(niter)
I2 = np.zeros(niter)
D2 = np.zeros(niter)
    
     
S1[0] = ps*N #total population
I1[0] = pi*N
R1[0] = pr*N
D1[0] = d*R1[0]

S2[0] = ps*N #total population
I2[0] = pi*N
R2[0] = pr*N
D2[0] = d*R2[0]

    

for j in range(niter-1):
    dI1dt=beta/N*S1[j]*I1[j]-(1-d)*gamma*I1[j]-d*rho*I1[j]
    I1[j+1] = I1[j] + dt*dI1dt
    dS1dt=-beta/N*S1[j]*I1[j+1] + r*R1[j]
    dR1dt=(1-d)*gamma*I1[j+1] - r*R1[j]
    dD1dt=d*rho*I1[j+1]
    S1[j+1] = S1[j] + dt*dS1dt
    R1[j+1] = R1[j] + dt*dR1dt
    D1[j+1] = D1[j] + dt*dD1dt  

    dS2dt=-beta/N*S2[j]*I2[j] + r*R2[j]
    dI2dt=beta/N*S2[j]*I2[j]-(1-d)*gamma*I2[j]-d*rho*I2[j]
    dR2dt=(1-d)*gamma*I2[j] - r*R2[j]
    dD2dt=d*rho*I2[j]
    S2[j+1] = S2[j] + dt*dS2dt
    I2[j+1] = I2[j] + dt*dI2dt
    R2[j+1] = R2[j] + dt*dR2dt
    D2[j+1] = D2[j] + dt*dD2dt
    
I=(I1-I2)/I1*100    
plt.plot(t, I)
plt.title('SIRS model test difference between Euler-Cromer and Euler methods \n (Infected Population)')
plt.xlabel('days elapsed since 0.01 percent of the population became infected')
plt.ylabel('% difference in population output')
