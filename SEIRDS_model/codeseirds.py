#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 15:22:19 2020

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

#initializing parameters
# time step (day)
dt = 0.01
# average number of people that come within infection range of an infected individual
# per day
# ASSUMES NO SOCIAL DISTANCING OR MITIGATION
beta = 0.25
# remove probability per day = 1/(recovery time)
# approximately 16 days (varies from 14 days to a month)
gamma = 0.07 # approximately 14 days infected after exposure

# total population
N = 5E4
# total time (days)
tottime = 365
# initial percent removed (immune)
pr = 0.0
# initial percent exposed
pe = 0.0
# initial percemt infected
pi = 0.0001
# initial percent susceptible
ps = 1-pr-pe-pi
# rate at which exposed become infected (delta)
delta = 0.20 # 5 days exposed
# death rate alpha
alpha = 0.05 #5%
# rate at which people die
rho = 0.05 #20 days infected before death

niter = int(math.ceil(tottime/dt))
t = np.arange(0, tottime, dt)   
S = np.zeros(niter)
R = np.zeros(niter)
I = np.zeros(niter)
E = np.zeros(niter)
D = np.zeros(niter)
    
     
S[0] = ps*N #total population
I[0] = pi*N
R[0] = pr*N
E[0] = pe*N
D[0] = 0
 

tau = 0.01 #varied rate which recovered individuals return to the susceptible
# population due to loss of immunity.
q = 1

for j in range(niter-1):
    dSdt = -q*beta/N*S[j]*I[j] + tau*R[j]
    dEdt = q*beta/N*S[j]*I[j] - delta*E[j]
    dIdt = delta*E[j] - (1-alpha)*gamma*I[j] - alpha*rho*I[j]
    dRdt = (1-alpha)*gamma*I[j] - tau*R[j]
    dDdt = alpha*rho*I[j]
    S[j+1] = S[j] + dt*dSdt
    E[j+1] = E[j] + dt*dEdt
    I[j+1] = I[j] + dt*dIdt
    R[j+1] = R[j] + dt*dRdt
    D[j+1] = D[j] + dt*dDdt
N=S+E+I+R+D
plt.plot(t, S)
plt.plot(t, E)
plt.plot(t, I)
plt.plot(t, R)
plt.plot(t, D)
plt.plot(t, N)
plt.legend(['susceptible', 'exposed', 'infected', 'recovered', 'dead', 'total'],
           loc='upper right')
plt.title('SEIRDS model')
plt.xlabel('days elapsed since 0.01 percent of the population became infected')
plt.ylabel('population')
#plt.xlim(250,2000)
#plt.ylim(0, 2E6)  
plt.show()
        