#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:58:43 2020

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
tottime = 500
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
day_lifted = 250
# death rate
d = 0.06
# rate at which people die (1/ time to die after infection)
rho = 1/20 #20 days to die?

niter = int(math.ceil(tottime/dt))
t = np.arange(0, tottime, dt)   
S = np.zeros(niter)
S1 = np.zeros(niter)
S2 = np.zeros(niter)
S3 = np.zeros(niter)
R = np.zeros(niter)
I = np.zeros(niter)
D = np.zeros(niter)
    
     
S[0] = ps*N #total population
S1[0]= psd*S[0] #extreme social distancers
S2[0]= pss*S[0] #super spreaders
S3[0]= S[0]-S1[0]-S2[0] #average people
I[0] = pi*N
R[0] = pr*N
D[0] = d*R[0]
    
for j in range(niter-1):
    if j <= day_lifted / dt:
        q = 0 #extreme isolation
        k = 1.5 #extreme spreading
        p = 0.5 #normal restricted spreading
    else:
        q = 0.5 #still some reduction
        k = 1.5 #still extreme spreading
        p = 1 #resumed normal behavior
    dS1dt=-q*beta/N*S1[j]*I[j]
    dS2dt=-k*beta/N*S2[j]*I[j]
    dS3dt=-p*beta/N*S3[j]*I[j]
    dIdt=q*beta/N*S1[j]*I[j]+k*beta/N*S2[j]*I[j]+p*beta/N*S3[j]*I[j]-(1-d)*gamma*I[j]-d*rho*I[j]
    dRdt=(1-d)*gamma*I[j]
    dDdt=d*rho*I[j]
    S1[j+1] = S1[j] + dt*dS1dt
    S2[j+1] = S2[j] + dt*dS2dt
    S3[j+1] = S3[j] + dt*dS3dt
    I[j+1] = I[j] + dt*dIdt
    R[j+1] = R[j] + dt*dRdt
    D[j+1] = D[j] + dt*dDdt
        
S=S1+S2+S3
N=S+I+R+D
plt.plot(t, S)   
plt.plot(t, S3)
plt.plot(t, S1)
plt.plot(t, S2)
plt.plot(t, I)
plt.plot(t, R)
plt.plot(t, D)
plt.plot(t, N) #check that the total population is constant
    
plt.legend(['Total susceptible', str(round(((1-pss-psd)*100),3)) + '% "average" people',
           str(round((psd*100),3)) + '% extreme isolators',
           str(round((pss*100),3)) + '% super spreaders', 'Infected', 'Recovered',
           'Dead', 'Total'], loc='upper right')
plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days')
plt.xlabel('days elapsed since 0.01 percent of the population became infected')
plt.ylabel('population')
#plt.ylim(0,2E7)
#plt.xlim(0,500)
plt.show()
