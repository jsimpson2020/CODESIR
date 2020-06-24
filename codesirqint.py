#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:05:37 2020

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
tottime = 365
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
day_lifted = 50
# rate of immunity loss (ratio of recovered people who lose immunity each day)
r = 0.001
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
R1 = np.zeros(niter)
R2 = np.zeros(niter)
R3 = np.zeros(niter)
I = np.zeros(niter)
I1 = np.zeros(niter)
I2 = np.zeros(niter)
I3 = np.zeros(niter)
D = np.zeros(niter)
D1 = np.zeros(niter)
D2 = np.zeros(niter)
D3 = np.zeros(niter)

I_T = np.zeros(niter)
I_T1 = np.zeros(niter)
I_T2 = np.zeros(niter)
I_T3 = np.zeros(niter)


    
     
S[0] = ps*N #total population
S1[0]= psd*S[0] #extreme social distancers
S2[0]= pss*S[0] #super spreaders
S3[0]= S[0]-S1[0]-S2[0] #average people
I[0] = pi*N
I1[0] = psd*I[0] #extreme social distancers
I2[0] = pss*I[0] #super spreaders
I3[0] = I[0]-I1[0]-I2[0] #average people
R[0] = pr*N
R1[0] = psd*R[0] #extreme social distancers
R2[0] = pss*R[0] #super spreaders
R3[0] = R[0]-R1[0]-R2[0] #average people
D[0] = d*R[0]
D1[0] = d*R1[0]
D2[0] = d*R2[0]
D3[0] = d*R3[0]

I_T[0] = I[0]
I_T1[0] = I1[0]
I_T2[0] = I2[0]
I_T3[0] = I3[0]


    
for n in range(20):
    for j in range(niter-1):
        if j <= (day_lifted + 10*n) / dt:
            q = 0 #extreme isolation
            k = 1 #extreme spreading
            p = 0.5 #normal restricted spreading
        else:
            q = 0.5 #still some reduction
            k = 1 #still extreme spreading
            p = 0.9 #resumed mostly normal behavior
        dS1dt=-q*beta/N*S1[j]*I[j] + r*R1[j]
        dS2dt=-k*beta/N*S2[j]*I[j] + r*R2[j]
        dS3dt=-p*beta/N*S3[j]*I[j] + r*R3[j]
        dIdt=q*beta/N*S1[j]*I[j]+k*beta/N*S2[j]*I[j]+p*beta/N*S3[j]*I[j]-(1-d)*gamma*I[j]-d*rho*I[j]
        dI1dt=q*beta/N*S1[j]*I[j]-(1-d)*gamma*I1[j]-d*rho*I1[j]
        dI2dt=k*beta/N*S2[j]*I[j]-(1-d)*gamma*I2[j]-d*rho*I2[j]
        dI3dt=p*beta/N*S3[j]*I[j]-(1-d)*gamma*I3[j]-d*rho*I3[j]
        dRdt=(1-d)*gamma*I[j] - r*R[j]
        dR1dt=(1-d)*gamma*I1[j] - r*R1[j]
        dR2dt=(1-d)*gamma*I2[j] - r*R2[j]
        dR3dt=(1-d)*gamma*I3[j] - r*R3[j]
        dDdt=d*rho*I[j]
        dD1dt=d*rho*I1[j]
        dD2dt=d*rho*I2[j]
        dD3dt=d*rho*I3[j]
        S1[j+1] = S1[j] + dt*dS1dt
        S2[j+1] = S2[j] + dt*dS2dt
        S3[j+1] = S3[j] + dt*dS3dt
        I[j+1] = I[j] + dt*dIdt
        I1[j+1] = I1[j] +dt*dI1dt
        I2[j+1] = I2[j] +dt*dI2dt
        I3[j+1] = I3[j] +dt*dI3dt
        R[j+1] = R[j] + dt*dRdt
        R1[j+1] = R1[j] + dt*dR1dt
        R2[j+1] = R2[j] + dt*dR2dt
        R3[j+1] = R3[j] + dt*dR3dt
        D[j+1] = D[j] + dt*dDdt
        D1[j+1] = D1[j] + dt*dD1dt
        D2[j+1] = D2[j] + dt*dD2dt
        D3[j+1] = D3[j] + dt*dD3dt
        I_T[j+1] = I_T[j] +I[j]*dt
        I_T1[j+1] = I_T1[j] +I1[j]*dt
        I_T2[j+1] = I_T2[j] +I2[j]*dt
        I_T3[j+1] = I_T3[j] +I3[j]*dt
    I_T = I_T*((1-d)*gamma+d*rho) #multiplied by gamma to adjust for overcountring
    I_T1 = I_T1*((1-d)*gamma+d*rho) #multiplied by gamma to adjust for overcountring
    I_T2 = I_T2*((1-d)*gamma+d*rho) #multiplied by gamma to adjust for overcountring
    I_T3 = I_T3*((1-d)*gamma+d*rho) #multiplied by gamma to adjust for overcountring
    L = [max(I_T).round(0), max(I_T1).round(0), max(I_T2).round(0), max(I_T3).round(0),
         max(D).round(0)]
    print(L)
    plt.plot(t, I_T, 'k')
    plt.plot(t, I_T1, 'b')
    plt.plot(t, I_T2, 'g')
    plt.plot(t, I_T3, 'r')
plt.show()