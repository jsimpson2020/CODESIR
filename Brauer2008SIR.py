#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 12:42:06 2020

@author: jsimpson
"""
#Uses Brauer 2008 Kermack-McKendrick model:

#S(t) is the number of individuals who are susceptible to the disease, that is,
# who are not (yet) infected at time t.
#I(t) is the number of infected individuals, assumed infectious and able to spread
# the disease by contact with susceptibles.
#R(t) is the number of individuals who have been infected and then removed from the
#possibility of being infected again or of spreading infection.

#Diff eq's:
    #S'=-betaSI
    #I'=betaSI-alphaI
    #R'=alphaI
    
#Model assumptions:
    #(1) An average member of the population makes contact sufficient to transmit
    # infection with betaN others per unit time, where N represents total
    # population size (mass action incidence).
    #(2) Infectives leave the infective class at a rate alphaI per unit time.
    #(3) There is no entry into or departure from the population, except possibly
    # through death from the disease.
    
# R0 is the basic reproduction number, which is betaS(0)/alpha

# importing packages
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math

#initializing parameters
# time step (day)
dt = 0.01
# fraction of total population that an average member of the population makes contact
# with sufficient to transmit infection per unit time.
# units are [1/N][1/t] let's choose N=individual, t=days
beta = 0.000595
# length of infective period is 1/alpha -> alpha is 1/(length of infective period)
# units are [1/t] let's choose t=days
alpha = 0.091


# total time (days)
tottime = 135

#Initializing arrays
niter = int(math.ceil(tottime/dt))
t = np.arange(0, tottime, dt)   
S = np.zeros(niter)
I = np.zeros(niter)

#Initial population values 
S[0] = 254
I[0] = 7
N = S[0] + I[0]

R0 = beta*S[0]/alpha


for j in range(niter-1):
    dSdt=-(beta)*S[j]*I[j]
    dIdt=(beta)*S[j]*I[j]-alpha*I[j]
    S[j+1] = S[j] + dt*dSdt 
    I[j+1] = I[j] + dt*dIdt  
   
       
R = N-S-I



 #plotting
plt.plot(t, S, 'k', label = 'susceptible')
plt.plot(t, I, 'm', label = 'infected')
plt.plot(t, R, 'b', label = 'recovered')
plt.plot(t, R+S+I, 'y', label = 'total')
plt.gca().legend(('susceptible','infected','recovered','total'))
plt.title('Brauer 2008 SIR model for Eyam Plague')
plt.xlabel('days since beginning of outbreak')
plt.ylabel('population')
plt.show()


 #specific plotting for Brauer example
plt.plot(S, I, 'k')
plt.title('Brauer 2008 SIR model for Eyam Plague')
plt.xlabel('S(t)')
plt.ylabel('I(t)')
plt.xlim(0, 260)
plt.ylim(0, 32)
plt.xticks(np.arange(0, 270, step=10), 
           [' ',' ',' ',' ',' ', 50,' ',' ',' ',' ',100,' ',' ',' ',' ',150,
            ' ',' ',' ',' ', 200,' ',' ',' ',' ', 250, ' '])
plt.yticks(np.arange(0, 32, step=1),
           [0,' ',' ',' ',' ',5,' ',' ',' ',' ',10,' ',' ',' ',' ',15,
            ' ',' ',' ',' ',20,' ',' ',' ',' ',25,' ',' ',' ',' ',30,' ', ' '])
plt.show()

