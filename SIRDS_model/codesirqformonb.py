#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 10:27:02 2020

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

    
for j in range(niter-1):
    if j <= day_lifted / dt:
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
        
S=S1+S2+S3
N=S+I+R+D
T1=S1+I1+R1+D1
T2=S2+I2+R2+D2
T3=S3+I3+R3+D3


for n in range(1,2):
    if n == 0:
        plt.plot(t, S)   
        plt.plot(t, S3)
        plt.plot(t, S1)
        plt.plot(t, S2)
        plt.legend(['Total susceptible', str(round(((1-pss-psd)*100),3)) + '% "average" people',
                    str(round((psd*100),3)) + '% extreme isolators',
                    str(round((pss*100),3)) + '% super spreaders'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n (Susceptible population)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()
    elif n == 1:
        plt.plot(t, I)
        plt.plot(t, I3)
        plt.plot(t, I1)
        plt.plot(t, I2)
        plt.legend(['Total infected', str(round(((1-pss-psd)*100),3)) + '% "average" people',
                    str(round((psd*100),3)) + '% extreme isolators',
                    str(round((pss*100),3)) + '% super spreaders'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n (Infected population)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()    
    elif n == 2:
        plt.plot(t, R)
        plt.plot(t, R3)
        plt.plot(t, R1)
        plt.plot(t, R2)
        plt.legend(['Total recovered', str(round(((1-pss-psd)*100),3)) + '% "average" people',
                    str(round((psd*100),3)) + '% extreme isolators',
                    str(round((pss*100),3)) + '% super spreaders'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n (Recovered population)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()
    elif n == 3:
        plt.plot(t, D)
        plt.plot(t, D3)
        plt.plot(t, D1)
        plt.plot(t, D2)
        plt.legend(['Total dead', str(round(((1-pss-psd)*100),3)) + '% "average" people',
                    str(round((psd*100),3)) + '% extreme isolators',
                    str(round((pss*100),3)) + '% super spreaders'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n (Dead population)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()
    elif n == 4:
        plt.plot(t, S1)
        plt.plot(t, I1)
        plt.plot(t, R1)
        plt.plot(t, D1)
        plt.plot(t, T1)
        plt.legend(['Susceptible', 'Infected', 'Recovered', 'Dead', 'Total'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n (Extreme isolators)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()    
    elif n == 5:
        plt.plot(t, S2)
        plt.plot(t, I2)
        plt.plot(t, R2)
        plt.plot(t, D2)
        plt.plot(t, T2)
        plt.legend(['Susceptible', 'Infected', 'Recovered', 'Dead', 'Total'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n (Super spreaders)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()
    elif n == 6:
        plt.plot(t, S3)
        plt.plot(t, I3)
        plt.plot(t, R3)
        plt.plot(t, D3)
        plt.plot(t, T3)
        plt.legend(['Susceptible', 'Infected', 'Recovered', 'Dead', 'Total'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n ("Average" people)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()    
    else:
        plt.plot(t, S)
        plt.plot(t, I)
        plt.plot(t, R)
        plt.plot(t, D)
        plt.plot(t, N)
        plt.legend(['Susceptible', 'Infected', 'Recovered', 'Dead', 'Total'], loc='upper right')
        plt.title('Mixed SIR model for lifting restrictions after ' + str(day_lifted)+ ' days'
                  + '\n (Total population)')
        plt.xlabel('days elapsed since 0.01 percent of the population became infected')
        plt.ylabel('population')
        plt.show()    

#plt.ylim(0,1E6)
#plt.xlim(0,250)'''
