#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 12:40:35 2020

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
x = 10
dt1 = 0.0001
dt2 = dt1*x
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
# death rate
d = 0.05
# rate at which people die (1/ time to die after infection)
rho = 1/20 #20 days to die?

niter1 = int(math.ceil(tottime/dt1))
t1 = np.arange(0, tottime, dt1)   
S1 = np.zeros(niter1)
R1 = np.zeros(niter1)
I1 = np.zeros(niter1)
D1 = np.zeros(niter1)
    
     
S1[0] = ps*N #total population
I1[0] = pi*N
R1[0] = pr*N
D1[0] = d*R1[0]

rvar = 0.001 #varied rate which recovered individuals return to the susceptible
# population due to loss of immunity.
q = 1


    

for j in range(niter1-1):
    dS1dt=-q*beta/N*S1[j]*I1[j] + rvar*R1[j]
    dI1dt=q*beta/N*S1[j]*I1[j]-(1-d)*gamma*I1[j]-d*rho*I1[j]
    dR1dt=(1-d)*gamma*I1[j] - rvar*R1[j]
    dD1dt=d*rho*I1[j]
    S1[j+1] = S1[j] + dt1*dS1dt
    I1[j+1] = I1[j] + dt1*dI1dt
    R1[j+1] = R1[j] + dt1*dR1dt
    D1[j+1] = D1[j] + dt1*dD1dt



niter2 = int(math.ceil(tottime/dt2))
t2 = np.arange(0, tottime, dt2)   
S2 = np.zeros(niter2)
R2 = np.zeros(niter2)
I2 = np.zeros(niter2)
D2 = np.zeros(niter2)
    
     
S2[0] = ps*N #total population
I2[0] = pi*N
R2[0] = pr*N
D2[0] = d*R2[0]

for j in range(niter2-1):
    dS2dt=-q*beta/N*S2[j]*I2[j] + rvar*R2[j]
    dI2dt=q*beta/N*S2[j]*I2[j]-(1-d)*gamma*I2[j]-d*rho*I2[j]
    dR2dt=(1-d)*gamma*I2[j] - rvar*R2[j]
    dD2dt=d*rho*I2[j]
    S2[j+1] = S2[j] + dt2*dS2dt
    I2[j+1] = I2[j] + dt2*dI2dt
    R2[j+1] = R2[j] + dt2*dR2dt
    D2[j+1] = D2[j] + dt2*dD2dt

I=np.zeros(niter2)
I[0]=0
    
for i in range(niter2-1):
    I[i+1]=(I1[x*(i+1)]-I2[i+1])/I1[x*(i+1)]*100
    

    
plt.plot(t2, I)    
#plt.plot(t1, I1)
#plt.plot(t2, I2)  
#plt.legend(['r = 1','r = 0.1', 'r = 0.01', 'r = 0.001', 'r = 0.0001',
#            'r = 0.00001'], loc='upper left')
plt.title('SIRS model test difference in dt (dt1 = ' + str(dt1) +
          ', dt2 = ' + str(dt2) + ')  \n (I1-I2) vs t (Infected Population)')
plt.xlabel('days elapsed since 0.01 percent of the population became infected')
plt.ylabel('% difference in population output')
#plt.xlim(600,730)
#plt.ylim(0, 3E5)  
        