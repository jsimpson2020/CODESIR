#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:27:13 2020

@author: jsimpson
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import random

def modifiedBrauerSIRSbd(initial_population=100000, length_of_infection=20,
                        avg_contacts_per_day=0.25, total_time=365,
                        initial_percent_recovered=0.0, initial_percent_infected=1.0,
                        fraction_recovered=1.0, deaths_per_100000_per_day=2.4, 
                        births_per_100000_per_day=3.3, length_of_immunity=90,
                        dt = 0.01, xmin=0, xmax=365, ymin=0, ymax=105000, plot=False,
                        parameters=False):
    
    random_number = random.random()
    
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
        if I[j] > 0:
            dSdt=nu*N[j]-lmbda/N[j]*S[j]*I[j]-mu*S[j]+r*R[j]
            dIdt=lmbda/N[j]*S[j]*I[j]-mu*I[j]-alpha*I[j]
            dRdt=f*alpha*I[j]-mu*R[j]-r*R[j]
            dNdt=nu*N[j]-(1-f)*alpha*I[j]-mu*N[j]
            S[j+1] = S[j] + dt*dSdt 
            I[j+1] = I[j] + dt*dIdt
            R[j+1] = R[j] + dt*dRdt
            N[j+1] = N[j] + dt*dNdt
            if random_number < 1/I[j]:
                R[j+1] = R[j] + dt*dRdt + I[j]
                I[j+1] = 0
        else:
            dSdt=nu*N[j]-mu*S[j]+r*R[j]
            dRdt=-mu*R[j]-r*R[j]
            dNdt=nu*N[j]-mu*N[j]
            S[j+1] = S[j] + dt*dSdt 
            R[j+1] = R[j] + dt*dRdt
            N[j+1] = N[j] + dt*dNdt
            
    
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
              'R0 = ' + str(R0approx), 'Random number = ' + str(random_number),
              sep = "\n")
    
    df = pd.DataFrame({'days':t, 'susceptible':S, 'infected':I, 'recovered':R, 'total':N})
    
    return df

modifiedBrauerSIRSbd(plot=True, length_of_immunity=1000)
