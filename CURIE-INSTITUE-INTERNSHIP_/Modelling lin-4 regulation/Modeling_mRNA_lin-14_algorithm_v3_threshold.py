# # # # # # # # Importing libraries # # # # # # # #

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# # # # # # # # Definition of the functions # # # # # # # #

def constant(y1,N):
    y = np.empty(0)
    for i in range(N):
        y = np.append(y,y1)
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

def juxtapose_two_signal_functions(f1,y0,y1,f,phi,t0,t1,n1,f2,cst,n2):
    # This function returns one list containing the values of the first function given in parameter and then the values of the second function given in parameter
    # f1 is typically an oscillating function as sinus or triangle
    # f2 is typically a constant or a slope
    # n1 and n2 are the number of points to compute for f1 and f2 respectively
    for i in range(n1):
        if(f1=="sinus"):
            sin = sinus(y0,y1,f,phi,t0,t1,n1)
    return

def verify_positive_value(x,epsilon):
    if x<epsilon:
        return 0
    else:
        return x

def C_Processing(L_alpha,L_beta,L_gamma,C_t,dt,nA_threshold):
    [alpha0 , alpha1] = L_alpha
    [beta0 , beta1 , beta2] = L_beta
    [gamma0 , gamma1] = L_gamma
    
    [nmiR4 , nmR14 , ncomplex , nP14 , nA] =  C_t
    real_nA = nA
    if (nA<nA_threshold):
        nA=0
    nP14_dt = nP14*(1-beta2*dt) + alpha1*nmR14*dt
    ncomplex_dt = ncomplex*(1-gamma1*dt) + gamma0*nmR14*nmiR4*dt
    nmR14_dt = nmR14*(1-beta1*dt-gamma0*nmiR4*dt) + alpha0*nA*dt
    nmiR4_dt = nmiR4*(1-beta0*dt-gamma0*nmR14*dt) + alpha0*nA*dt + gamma1*ncomplex*dt
    nA_dt = real_nA
    
    C_t_dt = [nmiR4_dt , nmR14_dt , ncomplex_dt , nP14_dt , nA_dt]
    
    return C_t_dt

# # # # # # # # Setting the global parameters for the Activator oscillations # # # # # # # #

y0  = 0     # minimum value of the activator signal
y1  = 1     # maximum value of the activator signal 
f   = 1     # frequency of the activator signal during each larval stage
phi = -90   # phase of the signal in DEGREES, the phase conversion from degrees to radians will be made in the processes
t0  = 0     # starting time for the modelling
t1  = 10    # ending time for the modelling
n_T = 4     # number of periods of the Activator signal to compute and show
n_osc = 100 # number of points in the plot of one oscillation
n_cst = 300 # number of points in the plot of the constant part of each period of the activator signal
n_post_oscillation = 0 # number of points to add at the end of the activator signal to plot the system evolution even after all the oscillations
N = n_T*(n_osc+n_cst) + n_post_oscillation

# # # # # # # # Declaring the chemical reaction constants # # # # # # # #

alpha0  = 20
beta0   = 0.2
gamma0  = 0.1
gamma1  = 8

alpha1  = 5
beta1   = 0.01
beta2   = 1

# # # # # # # # Generating all the activator concentration values in time # # # # # # # #

sin = sinus(y0,y1,f,phi,0,1/f,n_osc)
F_oscillation = sin
y_cst = F_oscillation[len(F_oscillation)-1]
F_constant = constant(y_cst,n_cst)
Pattern_function = np.append(F_oscillation,F_constant)
L_Activator = np.empty(0)


for i in range(n_T):
    # This for loop adds n_T time the pattern function to the Activator list
    L_Activator = np.append(L_Activator,Pattern_function)

# This line creates the constant value that will be taken by the activator after the oscillations
L_Activator_free_post_oscillation = constant(L_Activator[len(L_Activator)-1],n_post_oscillation)
# This line adds the constant function at the end of the Activator list
L_Activator = np.append(L_Activator,L_Activator_free_post_oscillation)

# # # # # # # # Declaring the inital values of concentrations # # # # # # # #

nA      = L_Activator[0]
nA_threshold = y1*30/100
nmiR4   = 0
nmR14   = 50
ncomplex= 0
nP14    = 250

# # # # # # # # Computing the useful constants and registering the first lists # # # # # # # #

t = t0
dt = (t1-t0)/(N-1)
if (dt>=1):
    print("Be careful to verify that the degradation rate * dt are lower than 1.")
T = 1/f
time_span = min(t1,n_T*T)
L_t = np.linspace(t0,t1,N,True)
cpt = 0 # variable counting the number of time the while loop is visited (see below)

L_alpha = [alpha0   ,   alpha1]
L_beta  = [beta0    ,   beta1   ,   beta2]
L_gamma = [gamma0   ,   gamma1]
C_t     = [nmiR4    ,   nmR14   ,   ncomplex , nP14 , nA]

# # # # # # # # Creating the result lists # # # # # # # #

# L_Activator_free is already created above
L_miR4  = [nmiR4]
L_mR14  = [nmR14]
L_complex = [ncomplex]
L_Pr14  = [nP14]

# # # # # # # # Launching the main computations # # # # # # # #

while (cpt<N-1 or t<=t1-dt):

    cpt = cpt + 1
    t = t + dt
    
    C_t = C_Processing(L_alpha,L_beta,L_gamma,C_t,dt,nA_threshold)

    # refreshing the current concentration values
    [nmiR4 , nmR14 , ncomplex , nP14 , nA] =  C_t
    nA      = L_Activator[cpt]
    C_t[4]  = nA
    nmiR4   = C_t[0]
    nmR14   = C_t[1]
    ncomplex= C_t[2]
    nP14   = C_t[3]
    
    # L_Activator_free is already filled in
    L_miR4.append(nmiR4)
    L_mR14.append(nmR14)
    L_complex.append(ncomplex)
    L_Pr14.append(nP14)

plt.subplot(231)
title = "Activator signal pattern"
plt.title(title)
m = min(len(L_t),len(Pattern_function))
Lt_norm = np.array(L_t[0:m])/max(L_t[0:m])
plt.plot(Lt_norm, Pattern_function[0:m] , c = 'r',label = 'Pattern function')
plt.plot(Lt_norm,constant(nA_threshold,m), label = 'Activator threshold')
plt.legend()
plt.xlabel('One larval stage')
plt.ylabel('Norm. conc.')
plt.ylim(-0.1,1.1)
plt.grid(True)

plt.subplot(232)
title = "Activator = sine with threshold"
plt.title(title)
m = min(len(L_t),len(L_Activator))
Lt_norm = 4*np.array(L_t[0:m])/max(L_t[0:m])
plt.plot (Lt_norm , L_Activator[0:m] , c = 'r')# , label = 'L_Activator')
plt.legend()
plt.xlabel('Scaled larval stages')
plt.ylabel('Norm. conc.')
plt.grid(True)

plt.subplot(233)
plt.title("miRNA lin-4 concentration in time")
plt.plot (L_t,L_miR4, c = 'g', label = 'L_miR4')
plt.legend()
plt.grid(True)

plt.subplot(234)
plt.title("mRNA lin-14 concentration in time")
plt.plot (L_t,L_mR14, c = 'b', label = 'L_mR14')
plt.legend()
plt.grid(True)

plt.subplot(235)
plt.title("LIN-14 concentration in time")
plt.plot (L_t,L_Pr14, c = 'y', label = 'L_Pr14')
plt.legend()
plt.grid(True)

plt.subplot(236)
plt.title("Constants - Parameters")
plt.xlim(0,10)
plt.ylim(-2,10)
plt.text(1,7,"Activator thresh.="+str(nA_threshold))
plt.text(1,5,"mRNA lin-14="+str(50))
plt.text(1,3,"alpha0  ="+str(alpha0))
plt.text(1,1,"alpha1  ="+str(alpha1))
plt.text(6,7,"beta0   ="+str(beta0))
plt.text(6,5,"beta1   ="+str(beta1))
plt.text(6,3,"beta2   ="+str(beta2))
plt.text(6,1,"gamma0  ="+str(gamma0))
plt.text(6,-1,"gamma1  ="+str(gamma1))
plt.plot()
plt.legend()
plt.grid(True)
plt.show()
