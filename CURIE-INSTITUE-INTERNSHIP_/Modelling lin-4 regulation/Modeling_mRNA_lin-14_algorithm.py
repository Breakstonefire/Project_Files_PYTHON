# # # # # # # # Importing libraries # # # # # # # #

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# # # # # # # # Definition of the main terms in the chemical equations # # # # # # # #

def constant(A,N):
    y = []
    for i in range(N):
        y.append(A)
    return y

def triangle(y0,y1,f,phi,t0,t1,N):
    # correction de la phase pour que le signal évolue comme un sinus lorsqu'on lui renseigne exactement les mêmes paramètres
    phi = phi*np.pi/180 # conversion préliminaire de la phase en radians
    phi = phi+np.pi/2

    A = (y1-y0)/2
    B = A+y0
    t = np.linspace(t0,t1,N,True)
    return A*signal.sawtooth(2*np.pi*f*t+phi,0.5)+B

def sinus(y0,y1,f,phi,t0,t1,N):
    phi = phi*np.pi/180 # conversion préliminaire de la phase en radians
    
    A = (y1-y0)/2
    B = B = A+y0
    t = np.linspace(t0,t1,N,True)
    return A*np.sin(2*np.pi*f*t+phi)+B

def C_Processing(alpha,beta,gamma,C_t,dt):
    n_lin14             = C_t[0]
    n_lin4_unActivated  = C_t[1]
    n_lin4_Activated    = C_t[2]
    n_Activator_free    = C_t[3]
    n_miR4_free         = C_t[4]
    n_mR14_free         = C_t[5]
    n_Pr14              = C_t[6]

    n_lin4_Activated_dt     = (1-gamma*dt)*n_lin4_Activated                     + gamma*n_lin4_unActivated*n_Activator_free*dt
    n_lin4_unActivated_dt   = (1-gamma*n_Activator_free*dt)*n_lin4_unActivated  + gamma*n_lin4_Activated*dt
    n_miR4_free_dt          = (1-beta*dt-gamma*n_mR14_free*dt)*n_miR4_free      + alpha*n_lin4_unActivated*dt - alpha*n_lin4_Activated*dt
    n_mR14_free_dt          = (1-beta*dt-gamma*n_miR4_free*dt)*n_mR14_free      + alpha*n_lin14*dt
    n_Pr14_dt               = (1-beta*dt)*n_Pr14                                + alpha*n_mR14_free*dt

    C_t_dt = [n_lin14,n_lin4_unActivated_dt,n_lin4_Activated_dt,n_Activator_free,n_miR4_free_dt,n_mR14_free_dt,n_Pr14_dt]
    return C_t_dt

# # # # # # # # Setting the global parameters for the Activator oscillations # # # # # # # #

y0  = 0
y1  = 2
f   = 1
phi = 0
t0  = 0
t1  = 4
N   = 1000

alpha = 1.1
beta = 1.2
gamma = 0.9

# # # # # # # # Computing the useful constants # # # # # # # #

t = t0
dt = (t1-t0)/(N-1)
L_t = np.linspace(t0,t1,N,True)
cpt = 0 # variable counting the number of time the while loop is visited (see below)

# # # # # # # # Generating all the activator concentration values in time # # # # # # # #

L_Activator_free = []
tri = triangle(y0,y1,f,phi,t0,t1,N)
sin = sinus(y0,y1,f,phi,t0,t1,N)
Magnitude = y1-y0
L_Activator_free = np.multiply(tri,sin/Magnitude)

# # # # # # # # Declaring the inital values of concentrations # # # # # # # #

n_lin14             = 7
n_lin4_unActivated  = 1
n_lin4_Activated    = 0
n_Activator_free    = L_Activator_free[0]
n_miR4_free         = 0
n_mR14_free         = 0
n_Pr14              = 0

C_t = [n_lin14,n_lin4_unActivated,n_lin4_Activated,n_Activator_free,n_miR4_free,n_mR14_free,n_Pr14]

# # # # # # # # Creating the result empty lists # # # # # # # #

L_lin14             = [n_lin14]
L_lin4_unActivated  = [n_lin4_unActivated]
L_lin4_Activated    = [n_lin4_Activated]
# L_Activator_free is already created above
L_miR4_free         = [n_miR4_free]
L_mR14_free         = [n_mR14_free]
L_Pr14              = [n_Pr14]

# # # # # # # # Launching the main computations # # # # # # # #

while (t<=t1-dt):
    cpt = cpt + 1
    C_t = C_Processing(alpha,beta,gamma,C_t,dt)
    # mise à jour des valeurs des concentrations
    n_lin14             = C_t[0]
    n_lin4_unActivated  = C_t[1]
    n_lin4_Activated    = C_t[2]
    n_Activator_free    = L_Activator_free[cpt]
    C_t[3]              = n_Activator_free
    n_miR4_free         = C_t[4]
    n_mR14_free         = C_t[5]
    n_Pr14              = C_t[6]

    L_lin14.append(n_lin14)
    L_lin4_unActivated.append(n_lin4_unActivated)
    L_lin4_Activated.append(n_lin4_Activated)
    # L_Activator_free is already filled in
    L_miR4_free.append(n_miR4_free)
    L_mR14_free.append(n_mR14_free)
    L_Pr14.append(n_Pr14)
    
    t = t + dt

plt.subplot(221)
plt.title("Activator: sine*triangle")
plt.plot (L_t,np.array(L_Activator_free)/2, c = 'r') #, label = 'L_Activator_free')
plt.xlabel('Scaled larval stages')
plt.legend()
plt.grid(True)

plt.subplot(222)
plt.title("C")
plt.plot (L_t,L_miR4_free, c = 'g', label = 'L_miR4_free')
plt.legend()
plt.grid(True)

plt.subplot(223)
plt.plot (L_t,L_mR14_free, c = 'b', label = 'L_mR14_free')
plt.legend()
plt.grid(True)

plt.subplot(224)
plt.title("C")
plt.plot (L_t,L_Pr14, c = 'y', label = 'L_Pr14')
plt.legend()
plt.grid(True)
plt.show()
