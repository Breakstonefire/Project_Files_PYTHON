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

def verify_positive_value(x,epsilon):
    if x<epsilon:
        return 0
    else:
        return x

def gaussian_part_of_solution(D,t0,t1,Nt,x0,x1,xmean,Nx):
    coef = 1/sqrt(4*pi*D)
    x = np.linspace(x0,x1,Nx,True)
    t = np.linspace(t0,t1,Nt,True)
    u = np.zeros(shape(Nt,Nx))
    for i in range(Nx):
        xi = x[i]
        ut = np.empty(0)
        for j in range(Nt):
            tj = t[j]
            uij = np.exp(-((xi-xmean)**2)/(4*D*tj))/sqrt(tj)
            ut = np.append(ut,uij)
        u[i] = ut
    u = u*coef
    return u

def Ch_t_dt(D,t0,t1,Nt,x0,x1,Nx,Ch_t,Sigma,i):
    dx = (x1-x0)/(Nx-1)
    dt = (t1-t0)/(Nt-1)
    coef = D*dt/(dx*dx)
    sigma_x = Sigma[i]
    Ch_x = Ch_t[i]

    if (i+1>=Nx):
        Ch_x_plus_dx = Ch_x
    else:
        Ch_x_plus_dx = Ch_t[i+1]
        
    if (i-1<0):
        Ch_x_minus_dx = Ch_x
    else:
        Ch_x_minus_dx = Ch_t[i-1]

    if Ch_x_plus_dx<0:
        Ch_x_plus_dx=0
    if Ch_x<0:
        Ch_x = 0
    if Ch_x_minus_dx<0:
        Ch_x_minus_dx = 0
    
    Ch_t_plus_dt = coef*(Ch_x_plus_dx -2*Ch_x + Ch_x_minus_dx)+ Ch_x +sigma_x*dt
    if Ch_t_plus_dt<0:
        return 0
    return Ch_t_plus_dt

def Cell_step(Ncell,x0,x1,Nx,N_free_left,N_free_right):
    #be careful to have at least Ncell<10*Nx to have a correct result to plot (in terms of visualiation).
    Delta = (Nx-N_free_left-N_free_right-1)*dx/(Ncell-1)
    return Delta

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
phi = -90   # phase of the signal in DEGREES, the phase conversion from degrees to radians will be made in the computations
t0  = 0     # starting time for the modelling
t1  = 10    # ending time for the modelling
n_T = 4     # number of periods of the Activator signal to compute and show
n_osc = 100 # number of points in the plot of one oscillation
n_cst = 300 # number of points in the plot of the constant part of each period of the activator signal
n_post_oscillation = 1000 # number of points to add at the end of the activator signal to plot the system evolution even after all the oscillations
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
y_cst = F_oscillation[len(F_oscillation)-1] #takes the last value of Foscillation
F_constant = constant(y_cst,n_cst)
Pattern_function_85 = np.append(F_oscillation,F_constant)
L_NHR85 = np.empty(0)
L_NHR23 = np.empty(0)
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
nA_threshold = y1*70/100
nmiR4   = 0
nmR14   = 50
ncomplex= 0
nP14    = 0

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





## beginning line 1 of plots

plt.subplot(341)
title = "Hormone Concentration pattern"
plt.title(title)
m = min(len(L_t),len(Pattern_function))
plt.plot(L_t[0:m], Pattern_function[0:m] , c = 'm',label = 'Pattern function')
plt.plot(L_t[0:m],constant(nA_threshold,m),label = 'Activator threshold')
plt.legend()
plt.ylim(-0.1,1.1)
plt.grid(True)

plt.subplot(342)
title = "Analytical solution without sink to the edges"
plt.title(title)
m = min(len(L_t),len(L_Activator))
plt.plot (L_t[0:m] , L_Activator[0:m] , c = 'r' , label = 'L_Activator')
plt.legend()
plt.grid(True)

plt.subplot(343)
plt.title("Hormone Creation-Degradation term")
plt.plot (L_t,L_miR4, c = 'g', label = 'L_miR4')
plt.legend()
plt.grid(True)

plt.subplot(344)
plt.title("Numerical solution with sink to the edges")
plt.plot (L_t,L_mR14, c = 'b', label = 'L_mR14')
plt.legend()
plt.grid(True)

## beginning line 2 of plots

plt.subplot(345)
plt.title("Nuclear Hormone Receptors concentrations")
plt.plot (L_t,L_Pr14, c = 'y', label = 'L_Pr14')
plt.legend()
plt.grid(True)

plt.subplot(346)
plt.title("Activator resulting concentration")
plt.xlim(0,10)
plt.ylim(-2,10)
plt.text(1,7,"Activator thresh.="+str(nA_threshold))
plt.text(1,5,"mRNA lin-14="+str(50))
plt.plot()
plt.legend()
plt.grid(True)

plt.subplot(347)
plt.title("Unactivated/Activated Activator concentrations - in time")
plt.xlim(0,10)
plt.ylim(-2,10)
plt.text(1,7,"Activator thresh.="+str(nA_threshold))
plt.text(1,5,"mRNA lin-14="+str(50))
plt.plot()
plt.legend()
plt.grid(True)

plt.subplot(348)
plt.title("Mean Unactivated/Activated Activator concentrations - in space")
plt.xlim(0,10)
plt.ylim(-2,10)
plt.text(1,7,"Activator thresh.="+str(nA_threshold))
plt.text(1,5,"mRNA lin-14="+str(50))
plt.plot()
plt.legend()
plt.grid(True)

## beginning line 3 of plots

plt.subplot(349)
plt.title("Unactivated/Activated Activator concentrations - in time")
plt.xlim(0,10)
plt.ylim(-2,10)
plt.text(1,7,"Activator thresh.="+str(nA_threshold))
plt.text(1,5,"mRNA lin-14="+str(50))
plt.plot()
plt.legend()
plt.grid(True)

plt.subplot(3410)
plt.title("Activated Activator concentrations - in time")
plt.xlim(0,10)
plt.ylim(-2,10)
plt.text(1,7,"Activator thresh.="+str(nA_threshold))
plt.text(1,5,"mRNA lin-14="+str(50))
plt.plot()
plt.legend()
plt.grid(True)

plt.subplot(3411)
plt.title("Characteristic timings from Activated Activator concentration - in space")
plt.xlim(0,10)
plt.ylim(-2,10)
plt.text(1,7,"Activator thresh.="+str(nA_threshold))
plt.text(1,5,"mRNA lin-14="+str(50))
plt.plot()
plt.legend()
plt.grid(True)



plt.show()
