#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:36:32 2020

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
tottime = 16000
# initial percent removed (immune)
pr = 0.0
# initial percent infected
pi = 0.0001
# initial percent susceptible
ps = 1-pr-pi
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

rvar = 0.01 #varied rate which recovered individuals return to the susceptible
# population due to loss of immunity.
q = 1
L=[]

    
for n in range(5):
    for j in range(niter-1):
        dSdt=-q*beta/N*S[j]*I[j] + (rvar+0.01*(n+0))*R[j]
        dIdt=q*beta/N*S[j]*I[j]-(1-d)*gamma*I[j]-d*rho*I[j]
        dRdt=(1-d)*gamma*I[j] - (rvar+0.01*(n+0))*R[j]
        dDdt=d*rho*I[j]
        S[j+1] = S[j] + dt*dSdt
        I[j+1] = I[j] + dt*dIdt
        R[j+1] = R[j] + dt*dRdt
        D[j+1] = D[j] + dt*dDdt
        
    lnI=np.log(I)
    m=np.polyfit(t[800000:], lnI[800000:], 1)
    L.append(m[0])
    
    plt.plot(t, lnI)    
plt.legend(['r = 0.01, m = ' + str(L[0]), 'r = 0.02, m = ' + str(L[1]),
            'r = 0.03, m = ' + str(L[2]), 'r = 0.04, m = ' + str(L[3]),
            'r = 0.05, m = ' + str(L[4])], loc='lower left')
plt.title('SIRS model with varied resusceptibility rates \n (Infected Population semi-log plot)')
plt.xlabel('days elapsed since 0.01 percent of the population became infected')
plt.ylabel('natural log of population')
#plt.xlim(250,2000)
plt.ylim(-40, 20)  
plt.show()
        



