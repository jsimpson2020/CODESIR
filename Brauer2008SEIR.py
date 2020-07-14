#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 10:41:03 2020

@author: jsimpson
"""


#Uses Brauer model:

#S(t) is the number of individuals who are susceptible to the disease, that is,
# who are not (yet) infected at time t.
#I(t) is the number of infected individuals, assumed infectious and able to spread
# the disease by contact with susceptibles.
#R(t) is the number of individuals who have been infected and then removed from the
#possibility of being infected again or of spreading infection.
#E(t) is the number of individuals who have been exposed, but are not infectious yet

#Diff eq's:
    #S'=-betaSI
    #E'=betaSI-kappaE
    #I'=kappaE-alphaI
    #R'=alphaI
    
#Model assumptions:
    #(1) An average member of the population makes contact sufficient to transmit
    # infection with betaN others per unit time, where N represents total
    # population size (mass action incidence).
    #(2) Infectives leave the infective class at a rate alphaI per unit time.
    #(3) There is no entry into or departure from the population, except possibly
    # through death from the disease.
    #(4) Incorporate an exponentially distributed exposed period with mean exposed
    # period 1/kappa
    
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
# length of exposed period is 1/kappa -> kappa is 1/(length of exposed period)
# units are [1/t] let's choose t=days
kappa = 0.2


# total time (days)
tottime = 135

#Initializing arrays
niter = int(math.ceil(tottime/dt))
t = np.arange(0, tottime, dt)   
S = np.zeros(niter)
E = np.zeros(niter)
I = np.zeros(niter)
R = np.zeros(niter)

#Initial population values 
S[0] = 247
E[0] = 7
I[0] = 7
R[0] = 0
N = S[0] + I[0] + E[0] + R[0]

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
   
       
N=S+E+I+R



 #plotting
plt.plot(t, S, 'k', label = 'susceptible')
plt.plot(t, E, 'g', label = 'exposed')
plt.plot(t, I, 'm', label = 'infected')
plt.plot(t, R, 'b', label = 'recovered')
plt.plot(t, R+S+I, 'y', label = 'total')
plt.gca().legend(('susceptible', 'exposed', 'infected','recovered','total'))
plt.title('Brauer 2008 SIR model for Eyam Plague with exposure')
plt.xlabel('days since beginning of outbreak')
plt.ylabel('population')
plt.show()