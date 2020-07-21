#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 14:18:33 2020

@author: jsimpson
"""

''' This document contains functions for each model that we will be using'''

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

### Basic SIR model from Brauer 2008
def BrauerSIR(initial_population=100000, length_of_infection=20,
              avg_contacts_per_day=0.25, total_time=365,
              initial_percent_recovered=0.0, initial_percent_infected=1.0,
              plot=False, parameters=False):
    ''' Uses Euler's method to solve a simple SIR model returns plot '''

    dt = 0.01
    
    lmbda = avg_contacts_per_day
    
    alpha = 1/length_of_infection


    niter = int(math.ceil(total_time/dt))
    t = np.arange(0, total_time, dt)   
    S = np.zeros(niter)
    I = np.zeros(niter)
    R = np.zeros(niter)
    
    pi = initial_percent_infected/100
    pr = initial_percent_recovered/100
    
    ps = 1 - pr - pi
    
    S[0] = ps*initial_population
    I[0] = pi*initial_population
    R[0] = pr*initial_population
    
    beta = lmbda/initial_population
    
    R0 = beta*S[0]/alpha

    for j in range(niter-1):
        dSdt = -beta*S[j]*I[j]
        dIdt = beta*S[j]*I[j] - alpha*I[j]
        dRdt = alpha*I[j]
        S[j+1] = S[j] + dt*dSdt
        I[j+1] = I[j] + dt*dIdt
        R[j+1] = R[j] + dt*dRdt
        if I[j] < 1.0:
            I[j+1]=0
    N = S + I + R
    
    if plot == True:
        plt.plot(t, S, 'k', label = 'susceptible')
        plt.plot(t, I, 'm', label = 'infected')
        plt.plot(t, R, 'b', label = 'recovered')
        plt.plot(t, N, 'y', label = 'total')
        plt.gca().legend(('susceptible','infected','recovered','total'),loc='right')
        plt.title('Brauer 2008 SIR model R0 = ' + str(R0.round(4)))
        plt.xlabel('days since ' + str(initial_percent_infected) + 
                   ' percent of population infected')
        plt.ylabel('population')
        plt.show()
        
    if parameters == True:
        print('Total time = ' + str(total_time) + ' days',
              'Average contacts with an infective sufficient to transmit disease per day = '
              + str(lmbda) + ' people',
              'Length of infection = ' + str(length_of_infection) + ' days',
              'Initial population = ' + str(initial_population) + ' people',
              'Initial percent recovered = ' +str(initial_percent_recovered) + '%',
              'Initial percent infected = ' + str(initial_percent_infected) + '%',
              sep = "\n")
    
    df = pd.DataFrame({'days':t, 'susceptible':S, 'infected':I, 'recovered':R, 'total':N})

    return df

### Basic SEIR model from Brauer 2008
def BrauerSEIR(initial_population=100000, length_of_infection=20,
               length_of_exposure=5, avg_contacts_per_day=0.25, total_time=365,
              initial_percent_recovered=0.0, initial_percent_infected=1.0,
              initial_percent_exposed=0.0, plot=False, parameters=False):
    ''' Uses Euler's method to solve a simple SEIR model '''

    dt = 0.01
    
    lmbda = avg_contacts_per_day
    
    alpha = 1/length_of_infection
    
    kappa = 1/length_of_exposure

    niter = int(math.ceil(total_time/dt))
    t = np.arange(0, total_time, dt)   
    S = np.zeros(niter)
    E = np.zeros(niter)
    I = np.zeros(niter)
    R = np.zeros(niter)
    
    pi = initial_percent_infected/100
    pe = initial_percent_exposed/100
    pr = initial_percent_recovered/100
    
    ps = 1 - pr - pi - pe
    
    S[0] = ps*initial_population
    E[0] = pe*initial_population
    I[0] = pi*initial_population
    R[0] = pr*initial_population
    
    beta = lmbda/initial_population
    
    R0 = beta*S[0]/alpha

    for j in range(niter-1):
        dSdt=-(beta)*S[j]*I[j]
        dEdt=beta*S[j]*I[j]-kappa*E[j]
        dIdt=kappa*E[j]-alpha*I[j]
        dRdt=alpha*I[j]
        S[j+1] = S[j] + dt*dSdt
        E[j+1] = E[j] + dt*dEdt
        I[j+1] = I[j] + dt*dIdt
        R[j+1] = R[j] + dt*dRdt
        if I[j] < 1.0:
            I[j+1]=0
    N=S+E+I+R
    
    if plot == True:
        plt.plot(t, S, 'k', label = 'susceptible')
        plt.plot(t, E, 'g', label = 'exposed')
        plt.plot(t, I, 'm', label = 'infected')
        plt.plot(t, R, 'b', label = 'recovered')
        plt.plot(t, N, 'y', label = 'total')
        plt.gca().legend(('susceptible', 'exposed', 'infected','recovered','total'),loc='right')
        plt.title('Brauer 2008 SEIR model R0 = ' + str(R0.round(4)))
        plt.xlabel('days since ' + str(initial_percent_infected) + 
                   ' percent of population infected')
        plt.ylabel('population')
        plt.show()
        
    if parameters == True:
        print('Total time = ' + str(total_time) + ' days',
              'Average contacts with an infective sufficient to transmit disease per day = '
              + str(lmbda) + ' people',
              'Length of infection = ' + str(length_of_infection) + ' days',
              'Length of exposure = ' + str(length_of_exposure) + ' days',
              'Initial population = ' + str(initial_population) + ' people',
              'Initial percent recovered = ' +str(initial_percent_recovered) + '%',
              'Initial percent infected = ' + str(initial_percent_infected) + '%',
              'Initial percent exposed = ' + str(initial_percent_exposed) + '%',
              sep = "\n")
        
    df = pd.DataFrame({'days':t, 'susceptible':S, 'exposed':E, 'infected':I,
                       'recovered':R, 'total':N})
    
    return df

### Basic SIS model from Brauer 2008
def BrauerSIS(initial_population=100000, length_of_infection=20,
              avg_contacts_per_day=0.25, total_time=365,
              initial_percent_recovered=0.0, initial_percent_infected=1.0,
              fraction_recovered=1, plot=False, parameters=False):
    ''' Uses Euler's method to solve a simple SIS model returns plot '''

    dt = 0.01
    
    lmbda = avg_contacts_per_day
    
    alpha = 1/length_of_infection
    
    f = fraction_recovered

    niter = int(math.ceil(total_time/dt))
    t = np.arange(0, total_time, dt)   
    S = np.zeros(niter)
    I = np.zeros(niter)
    R = np.zeros(niter)
    
    pi = initial_percent_infected/100
    pr = initial_percent_recovered/100
    
    ps = 1 - pr - pi
    
    S[0] = ps*initial_population
    I[0] = pi*initial_population
    R[0] = pr*initial_population
    
    beta = lmbda/initial_population
    
    R0 = beta*S[0]/alpha

    for j in range(niter-1):
        dSdt = -beta*S[j]*I[j] + f*alpha*I[j]
        dIdt = beta*S[j]*I[j] - alpha*I[j]
        dRdt = (1-f)*alpha*I[j]
        S[j+1] = S[j] + dt*dSdt
        I[j+1] = I[j] + dt*dIdt
        R[j+1] = R[j] + dt*dRdt
        if I[j] < 1.0:
            I[j+1]=0
    N = S + I + R
    
    if plot == True:
        plt.plot(t, S, 'k', label = 'susceptible')
        plt.plot(t, I, 'm', label = 'infected')
        plt.plot(t, R, 'b', label = 'recovered')
        plt.plot(t, N, 'y', label = 'total')
        plt.gca().legend(('susceptible','infected','recovered','total'),loc='right')
        plt.title('Brauer 2008 SIS model R0 = ' + str(R0.round(4)))
        plt.xlabel('days since ' + str(initial_percent_infected) + 
                   ' percent of population infected')
        plt.ylabel('population')
        plt.show()
        
    if parameters == True:
        print('Total time = ' + str(total_time) + ' days',
              'Average contacts with an infective sufficient to transmit disease per day = '
              + str(lmbda) + ' people',
              'Length of infection = ' + str(length_of_infection) + ' days',
              'Fraction of population who recover  = ' + str(fraction_recovered),
              'Initial population = ' + str(initial_population) + ' people',
              'Initial percent recovered = ' +str(initial_percent_recovered) + '%',
              'Initial percent infected = ' + str(initial_percent_infected) + '%',
              sep = "\n")
    
    df = pd.DataFrame({'days':t, 'susceptible':S, 'infected':I, 'recovered':R, 'total':N})

    return df

### SIR model with birth and death from Brauer 2008
def BrauerSIRbirthdeath(initial_population=100000, length_of_infection=20,
                        avg_contacts_per_day=0.25, total_time=365,
                        initial_percent_recovered=0.0, initial_percent_infected=1.0,
                        fraction_recovered=1.0, deaths_per_100000_per_day=2.4, 
                        births_per_100000_per_day=3.3, xmin=0, xmax=365, ymin=0,
                        ymax=105000, plot=False, parameters=False):
        
    dt = 0.01
    
    lmbda = avg_contacts_per_day
    
    alpha = 1/length_of_infection
    
    mu = deaths_per_100000_per_day/100000
    
    nu = births_per_100000_per_day/100000
    
    f = fraction_recovered
    
    niter = int(math.ceil(total_time/dt))
    t = np.arange(0, total_time, dt)   
    S = np.zeros(niter)
    I = np.zeros(niter)
    R = np.zeros(niter)
    N = np.zeros(niter)
    
    pi = initial_percent_infected/100
    pr = initial_percent_recovered/100
    ps = 1-pi-pr

    S[0] = ps*initial_population
    I[0] = pi*initial_population
    R[0] = pr*initial_population
    N[0] = initial_population
    
    for j in range(niter-1):
        dSdt=nu*N[j]-lmbda/N[j]*S[j]*I[j]-mu*S[j]
        dIdt=lmbda/N[j]*S[j]*I[j]-mu*I[j]-alpha*I[j]
        dNdt=nu*N[j]-(1-f)*alpha*I[j]-mu*N[j]
        S[j+1] = S[j] + dt*dSdt 
        I[j+1] = I[j] + dt*dIdt
        N[j+1] = N[j] + dt*dNdt
        if I[j] < 1.0:
            I[j+1]=0
    R = N-S-I
    
    R0 = lmbda/(mu + alpha)
    
    R0approx = round(R0, 3)
    
    if plot == True:
        plt.plot(t, S, 'k')
        plt.plot(t, I, 'm')
        plt.plot(t, R, 'b')
        plt.plot(t, N, 'y')
        plt.gca().legend(('susceptible','infected','recovered','total'),loc='right')
        plt.title('Brauer 2008 SIR model with birth and death, R0 = ' + str(R0approx))
        plt.xlabel('days since ' + str(initial_percent_infected) +
                   ' percent of population infected')
        plt.ylabel('population')
        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)
        plt.show()
    
    if parameters == True:
        print('Total time = ' + str(total_time) + ' days',
              'Average contacts with an infective sufficient to transmit disease per day = '
              + str(lmbda) + ' people',
              'Length of infection = ' + str(length_of_infection) + ' days',
              'Deaths per 100000 per day = ' + str(deaths_per_100000_per_day) + ' deaths',
              'Births per 100000 per day = ' + str(births_per_100000_per_day) + ' births',
              'Fraction of population who recover = ' + str(fraction_recovered),
              'Initial population = ' + str(initial_population) + ' people',
              'Initial percent recovered = ' +str(initial_percent_recovered) + '%',
              'Initial percent infected = ' + str(initial_percent_infected) + '%',
              'Total population at the end of 1st year = ' + str(N[36600].round(0)),
              sep = "\n")
        
    df = pd.DataFrame({'days':t, 'susceptible':S, 'infected':I, 'recovered':R, 'total':N})
    
    return df

### SIR model with birth, death and temporary immunity from Brauer 2008
def BrauerSIRSbirthdeath(initial_population=100000, length_of_infection=20,
                        avg_contacts_per_day=0.25, total_time=365,
                        initial_percent_recovered=0.0, initial_percent_infected=1.0,
                        fraction_recovered=1.0, deaths_per_100000_per_day=2.4, 
                        births_per_100000_per_day=3.3, length_of_immunity=90,
                        xmin=0, xmax=365, ymin=0, ymax=105000, plot=False,
                        parameters=False):
    
    dt = 0.01
    
    lmbda = avg_contacts_per_day
    
    alpha = 1/length_of_infection
    
    mu = deaths_per_100000_per_day/100000
    
    nu = births_per_100000_per_day/100000
    
    f = fraction_recovered
    
    r = 1/length_of_immunity
    
    niter = int(math.ceil(total_time/dt))
    t = np.arange(0, total_time, dt)   
    S = np.zeros(niter)
    I = np.zeros(niter)
    R = np.zeros(niter)
    N = np.zeros(niter)
    
    pi = initial_percent_infected/100
    pr = initial_percent_recovered/100
    ps = 1-pi-pr

    S[0] = ps*initial_population
    I[0] = pi*initial_population
    R[0] = pr*initial_population
    N[0] = initial_population
    
    for j in range(niter-1):
        dSdt=nu*N[j]-lmbda/N[j]*S[j]*I[j]-mu*S[j]+r*R[j]
        dIdt=lmbda/N[j]*S[j]*I[j]-mu*I[j]-alpha*I[j]
        dRdt=f*alpha*I[j]-mu*R[j]-r*R[j]
        dNdt=nu*N[j]-(1-f)*alpha*I[j]-mu*N[j]
        S[j+1] = S[j] + dt*dSdt 
        I[j+1] = I[j] + dt*dIdt
        R[j+1] = R[j] + dt*dRdt
        N[j+1] = N[j] + dt*dNdt
        if I[j] < 1.0:
            I[j+1]=0
    
    R0 = lmbda/(mu + alpha)
    
    R0approx = round(R0, 3)
    
    if plot == True:
        plt.plot(t, S, 'k')
        plt.plot(t, I, 'm')
        plt.plot(t, R, 'b')
        plt.plot(t, N, 'y')
        plt.gca().legend(('susceptible','infected','recovered','total'),loc='right')
        plt.title('SIR model with birth, death and temporary immunity \n R0 = ' + str(R0approx))
        plt.xlabel('days since ' + str(initial_percent_infected) +
                   ' percent of population infected')
        plt.ylabel('population')
        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)
        plt.show()
    
    if parameters == True:
        print('Total time = ' + str(total_time) + ' days',
              'Average contacts with an infective sufficient to transmit disease per day = '
              + str(lmbda) + ' people',
              'Length of infection = ' + str(length_of_infection) + ' days',
              'Deaths per 100000 per day = ' + str(deaths_per_100000_per_day) + ' deaths',
              'Births per 100000 per day = ' + str(births_per_100000_per_day) + ' births',
              'Fraction of population who recover with temporary immunity = '
              + str(fraction_recovered),
              'Length of immunity = ' + str(length_of_immunity) + ' days',
              'Initial population = ' + str(initial_population) + ' people',
              'Initial percent recovered = ' +str(initial_percent_recovered) + '%',
              'Initial percent infected = ' + str(initial_percent_infected) + '%',
              sep = "\n")
    
    df = pd.DataFrame({'days':t, 'susceptible':S, 'infected':I, 'recovered':R, 'total':N})
    
    return df
