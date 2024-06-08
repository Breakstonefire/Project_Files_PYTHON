# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing libraries # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from math import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Definition of the functions # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def constant(y1,N):
    y = np.empty(0)
    for i in range(N):
        y = np.append(y,y1)
    return y

def sinus(y0,y1,f,phi,t0,t1,dt):
    N = 1+int(round((t1-t0)/dt))
    phi = phi*np.pi/180 # conversion préliminaire de la phase en radians

    A = (y1-y0)/2
    B = A+y0
    t = np.linspace(t0,t1,N,True)
    return A*np.sin(2*np.pi*f*t+phi)+B

def compute_truncated_sin(s,f,phi,nT,dt):
    T = 1/f
    w = 2*np.pi*f
    t0 = (asin(s)-phi)/w
    if t0<0:
        while(t0<0):
            t0+=T
    t1 = t0+nT*T
    N = 1+int(round((t1-t0)/dt))
    L_t = np.linspace(t0,t1,N,True)
    sin = np.sin(w*L_t+phi)
    return [L_t,sin]

def Activator(A1,A2,f23,f85,delta_phi,t0,t1,dt):
    N = 1+int(round((t1-t0)/dt))
    t =  np.linspace(t0,t1,N,True)
    w23 = 2*np.pi*f23
    w85 = 2*np.pi*f85
    s = (A1*A2/4)*(1+np.sin(w23*t))*(1+np.sin(w85*t+delta_phi))
    return s

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Computing the analytical solution for the miRNA lin4 evolution during the non-null activator signal # # # # # # # #

## DURING THE ACTIVATION PEAK:
## CONDITIONS to use this solution:
# The lin4 anaytical solution is written for this chemical equation : dlin4/dt = alpha_lin4*Activator(t) - beta_lin4*lin4
# Where Activator(t) = 0.5*A1*(sin(2pi*f*t-pi/2)+1)*0.5*A2*(sin(2pi*f*t-pi/2+delta_phi)+1)
# Compute the result only during the moment NHR23 and NHR85 sin waves are superposing and only beetween twonull values of the resulting solution
# Between two null values of the analytical solution, it corresponds to the peak of {NHR23 ; NHR85} we are looking to implement.

def Sol_Homog_Particular(A1,A2,alpha_lin4,beta_lin4,lin4_t0,f,phi0,delta_phi,t0,t1,N):
    a = alpha_lin4
    b = beta_lin4
    b2 = pow(beta_lin4,2)
    w = 2*np.pi*f
    w2 = pow(w,2)
    cphi = np.cos(delta_phi)
    sphi = np.sin(delta_phi)

    t       = np.linspace(t0,t1,N,True)

    sft     = np.sin(w*t+phi0)
    s2ft    = np.sin(2*w*t+2*phi0)
    cft     = np.cos(w*t+phi0)
    c2ft    = np.cos(2*w*t+2*phi0)
    exp_plus_t  = np.exp(b*t)
    exp_minus_t = np.exp(-b*t)
    
    sft0     = np.sin(w*t0+phi0)
    s2ft0    = np.sin(2*w*t0+2*phi0)
    cft0     = np.cos(w*t0+phi0)
    c2ft0    = np.cos(2*w*t0+2*phi0)
    exp_plus_t0  = np.exp(b*t0)
    exp_minus_t0 = np.exp(-b*t0)

#   Sol_Homogeneous corresponds to the Homogeneous equation solution computed beetween t0 and t1.
    Sol_Homogeneous = lin4_t0*(exp_minus_t - exp_minus_t0)

#   Sol_0f corresponds to the constant part of the particular equation solution computed beetween t0 and t1.
    constant = (a*np.sqrt(A1*A2)/4)*(1+0.5*cphi)
    Sol_Particular_0f = constant*(exp_plus_t - exp_plus_t0)
#   Sol_1f corresponds to the sine and cosine terms depending on f in the particular equation solution computed beetween t0 and t1.
    constant = (a*np.sqrt(A1*A2)/4)/(b2+w2)  
    Sol_Particular_1f = constant*( exp_plus_t*((b-w)*sft-(b+w)*cft) - exp_plus_t0*((b-w)*sft0-(b+w)*cft0) )
#   Sol_2f corresponds to the sine and cosine terms depending on 2f in the particular equation solution computed beetween t0 and t1.
    constant = (-a*np.sqrt(A1*A2)/2)/(b2+4*w2)
    Sol_Particular_2f = constant*( exp_plus_t*(b*s2ft-w*c2ft) - exp_plus_t0*(b*s2ft0-w*c2ft0))
    Whole_solution = lin4_t0 + Sol_Homogeneous + Sol_Particular_0f + Sol_Particular_1f + Sol_Particular_2f

    for i in range(len(Whole_solution)):
        if Whole_solution[i]<0:
            Whole_solution[i]=0
    return Whole_solution

## AFTER THE ACTIVATION PEAK :
## CONDITIONS to use this solution:
# The lin4 anaytical solution is now written for this chemical equation : dlin4/dt = -beta_lin4*lin4
# The solution is a classic decreasing exponential

def Clin4_analytic_after_activator_peak(lin4_t0,beta_lin4,t0,t1,N):
    if N==0:
        return np.empty(0)
    t = np.linspace(t0,t1,N,True)
    sol = lin4_t0*np.exp(-beta_lin4*t)
    return sol

def C_Processing(L_alpha,L_beta,gamma_lin4_lin14,C_t,dt):
    [alpha_lin4 , alpha_lin14   , alpha_LIN14  ] = L_alpha
    [beta_lin4  , beta_lin14    , beta_LIN14   ] = L_beta
    [lin4       , lin14         , LIN14     , A] = C_t

    lin4_dt  = lin4  * (1-beta_lin4 *dt)+ alpha_lin4*A*dt
    lin14_dt = lin14 * (1-beta_lin14*dt - gamma_lin4_lin14*lin4*dt) + alpha_lin14*A*dt
    LIN14_dt = LIN14 * (1-beta_LIN14*dt)+ alpha_LIN14*lin14*dt
    
    if(lin4_dt<0):
        lin4_dt=0
    if(lin14_dt<0):
        lin14_dt=0
    if(LIN14_dt<0):
        LIN14_dt=0
    C_t_dt   = [lin4_dt , lin14_dt , LIN14_dt , A]

    return C_t_dt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Setting the global parameters for the Activator oscillations  # # # # # # # # # # # # #

NHR23_min = 0                   # minimum value of the NHR23 signal
NHR23_max = 10000                # maximum value of the NHR23 signal
A23 = NHR23_max-NHR23_min
f23 = 1                         # frequency of the NHR23 signal during each larval stage
T23 = 1/f23
t0_23 = 0
t1_23 = t0_23 + T23
phi23 = (-np.pi/2)              # phase of the NHR23 signal in RADIANS

delta_phi = (6.13-2.37)*180/np.pi                # enter the phase shift in DEGREES in this line (it will be converted in radians after)
delta_phi*=-np.pi/180

NHR85_min = 0                   # minimum value of the NHR85 signal
NHR85_max = 10000                # maximum value of the NHR85 signal
A85 = NHR85_max-NHR85_min
phi85 = (-np.pi/2)              # phase of the NHR85 signal in RADIANS
f85 = 1                         # frequency of the NHR85 signal during each larval stage
T85 = 1/f85
t0_85 = T85*((3/4)-phi85-delta_phi)
t1_85 = t0_85 + 85

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Declaring the chemical reaction constants # # # # # # # # # # # # # # # # # # # # # # #

# ASSUMPTIONS :
### (1)
# LIN14 protein is 539 units (=3 nucleotides) long so the mRNA lin-14 is approximately 3*539 = 1617 nucleotides long : https://www.uniprot.org/uniprot/Q21446
# lin-4 miRNA is only 22 nucleotides long
#   => one miRNA lin4 is read 1617/22= (~) 80 times faster than mRNA lin14
#   => the natural production rate of lin4 and lin14 are related the same way:
#   => alpha_lin4 = 1/tau_lin4 = 80/tau_lin14 = 80*alpha_lin14
#   => alpha_lin4 = 80*alpha_lin14

### (2)
# tau = half-life time of lin4 mRNA =~ 20.7h
#   => dlin14/dt =~ (after the activator has disappeared) -beta_lin14*lin14
#   => beta_lin14 = ln(2)/tau
#   => beta_lin14 =~ 10^-5 sec^-1

### (3)
# beta_lin4 << beta_lin14

# Known parameters due to measurements - reasonable assumptions
tau_lin4 = 21*3600                      # corresponds to 21 hours in sec
tau_lin14 = 5*3600                      # corresponds to 2 hours in sec

tau_lin4/=10*3600                       # here we divide by the duration in seconds of one larval stage (We do the assumption that a larval stage at 25°C has always the same duration which is 10*3600 sec = 10 hours)
tau_lin14/=10*3600                      # here we divide by the duration in seconds of one larval stage (We do the assumption that a larval stage at 25°C has always the same duration which is 10*3600 sec = 10 hours)

beta_lin4           = log(2)/tau_lin4   # (2) mRNAlin4 has a half-time period equals to ~21h => beta_lin14 =~ 10^-5 sec^-1 (at 25°C)
beta_lin14          = log(2)/tau_lin14  # (3) the lin4 miRNA is relatively less degraded compared to lin14 mRNA

# Deduced parameters due to relations with other parameters
alpha_lin14         = 0.2*beta_lin14      # (1) this assumption comes from the ratio of gene's length of lin4 and lin14
gamma_lin4_lin14    = 0#beta_lin4*beta_lin14/700

# Unknown parameters
alpha_lin4          = alpha_lin14
alpha_LIN14         = alpha_lin14/1000 #at beta_LIN14 constant, diminishing alpha_LIN14 causes a fall of the amplitude bumps, it dampens the oscillations of LIN-14
beta_LIN14          = beta_lin14


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Declaring the inital values of concentrations # # # # # # # #

# ASSUMPTIONS :
### (1)
# N(lin14 ; t=0) =~0
# N(lin4  ; t=0) =~0

### (2)
# LIN14(t=0)>>1 , very bright when the embryo hatches and then decreasing very fast at 9h in L1
# Number of cells in the C. elegans just after hatching is 558
# N(LIN-14 ; t=0) = ?
# (see this article) -> Caenorhabditis Elegans: Development from the Perspective of the Individual Cell - https://www.ncbi.nlm.nih.gov/books/NBK26861/#:~:text=Almost%20Perfectly%20Predictable-,C.,worm%20inside%20the%20egg%20shell.

n_lin4  = 0         # (1)
n_lin14 = 3000       # (1)
n_LIN14 = 1000       # (2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Generating all the activator concentration values in time # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # #  Computing the Activator signal # # # # # # # # # # # # # # #

n_T = 4                     # number of periods of the whole Activator signal to compute and show (oscillating part + constant part = 1 period)
n_post_oscillation = 0

dt = 1/(10000*min(f23,f85))

[t23,sin23] = compute_truncated_sin(-1,f23,phi23,1,dt)
[t85,sin85] = compute_truncated_sin(-1,f85,phi85+delta_phi,1,dt)

sin23 = A23*(1+sin23)/2
t0_23 = t23[0]
t1_23 = t23[-1]

sin85 = A85*(1+sin85)/2
t0_85 = t85[0]
t1_85 = t85[-1]

t85_i = t85[-1]
t23_i = t23[0]

N_flat23 = 0
while(t85_i>t1_23 and 0<=N_flat23<=len(t85)):
    N_flat23+=1
    t85_i = t85[len(t85)-N_flat23]
N_flat85 = 0
while(t23_i<t0_85 and 0<=N_flat85<=len(t23)):
    N_flat85+=1
    t23_i = t23[N_flat85]
if(N_flat23!=N_flat85):
    diff = abs(N_flat23-N_flat85)
    print("WARNING in 'Computing the Activator signal' part")
    print("Be careful the flat parts of the two NHRs haven't the same length, {",N_flat23,";",N_flat85,"}")
    if(N_flat23<N_flat85):
        message = "NHR-23"
        N_flat23 = N_flat85
    else:
        message = "NHR-85"
        N_flat85=N_flat23
    print("WARNING in 'Computing the Activator signal' part")
    print("The flat part of "+message+" has been completed with",diff,"point(s).")
Flat_signal_23 = constant(0,N_flat23)
Pattern_function_23 = np.append(sin23,Flat_signal_23)
Flat_signal_85 = constant(0,N_flat85)
Pattern_function_85 = np.append(Flat_signal_85,sin85)

N_pattern = len(Pattern_function_23)
t_pattern = np.linspace(t0_23,t1_85,N_pattern,True)
Pattern_function = Pattern_function_23 * Pattern_function_85

N_whole_signal = n_T*N_pattern + n_post_oscillation
Lt_whole_signal = np.linspace(t0_23,n_T*t1_85,N_whole_signal,True)

L_Activator = np.empty(0)
Amplitude_correction_term = 1/np.sqrt(A23*A85)
Pattern_function*=Amplitude_correction_term

for i in range(n_T):
    # This 'for' loop adds n_T time the pattern function to the whole Activator list
    L_Activator = np.append(L_Activator,Pattern_function)

# This line creates the constant value that will be taken by the activator after the oscillations
L_Activator_post_oscillation = constant(L_Activator[len(L_Activator)-1],n_post_oscillation)

# This line adds the constant function at the end of the Activator list
L_Activator_with_post_oscillation = np.append(L_Activator,L_Activator_post_oscillation)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Registering the first lists # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

L_alpha = [alpha_lin4   , alpha_lin14   , alpha_LIN14   ]
L_beta  = [beta_lin4    , beta_lin14    , beta_LIN14    ]
C_t     = [n_lin4       , n_lin14       , n_LIN14       , L_Activator[0]]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Creating the result lists # # # # # # # # # # # # # # # # # # #
nA              = L_Activator[0]
lin4_t0         = n_lin4
lin4_numeric    = [lin4_t0]
lin14_numeric   = [n_lin14]
LIN14_numeric   = [n_LIN14]
lin4_analytic   = np.empty(0)

t0_analytic = t0_23             # temps de début d'un pattern
t1_analytic = t0_85 - dt        # temps de fin de la partie où le signal est plat
t2_analytic = t1_analytic + dt  # temps de début de la partie où le signal est oscillant
t3_analytic = t1_23             # temps de fin de la partie où le signal est oscillant
t4_analytic = t3_analytic + dt  # temps de début de la partie où le signal est plat après la première oscillation
t5_analytic = t1_85             # temps de fin d'un pattern

N_analytic_before = N_flat85
N_analytic_during = N_pattern - N_flat23 - N_flat85
N_analytic_after  = N_flat23

deltaT1 = t1_analytic - t0_analytic
deltaT2 = t3_analytic - t2_analytic
deltaT3 = t5_analytic - t4_analytic

lin4_analytic = Clin4_analytic_after_activator_peak(lin4_t0,beta_lin4,t0_analytic,t1_analytic,N_analytic_before)
lin4_t2analytic = lin4_analytic[-1]
lin4_analytic = np.append(lin4_analytic,Sol_Homog_Particular(NHR23_max-NHR23_min,NHR85_max-NHR85_min,alpha_lin4,beta_lin4,lin4_t2analytic,f23,phi23,delta_phi,t2_analytic,t3_analytic,N_analytic_during))
lin4_t4analytic = lin4_analytic[-1]
lin4_analytic = np.append(lin4_analytic,Clin4_analytic_after_activator_peak(lin4_t4analytic,beta_lin4,dt,dt+deltaT3,N_analytic_after))

for i in range(n_T-1):
    t0_analytic = t5_analytic + dt          # temps de début d'un pattern
    t1_analytic = t0_analytic + deltaT1     # temps de fin de la partie où le signal est plat
    t2_analytic = t1_analytic + dt          # temps de début de la partie où le signal est oscillant
    t3_analytic = t2_analytic + deltaT2     # temps de fin de la partie où le signal est oscillant
    t4_analytic = t3_analytic + dt          # temps de début de la partie où le signal est plat après la première oscillation
    t5_analytic = t4_analytic + deltaT3     # temps de fin d'un pattern

    lin4_t0         = lin4_analytic[-1]
    lin4_analytic   = np.append(lin4_analytic,Clin4_analytic_after_activator_peak(lin4_t0,beta_lin4,dt,dt+deltaT1,N_analytic_before))
    lin4_t0         = lin4_analytic[-1]
    lin4_analytic   = np.append(lin4_analytic,Sol_Homog_Particular(NHR23_max-NHR23_min,NHR85_max-NHR85_min,alpha_lin4,beta_lin4,lin4_t0,f23,phi23,delta_phi,t0_85,t1_23,N_analytic_during))
    lin4_t0         = lin4_analytic[-1]
    lin4_analytic   = np.append(lin4_analytic,Clin4_analytic_after_activator_peak(lin4_t0,beta_lin4,dt,dt+deltaT3,N_analytic_after))

lin4_t0         = lin4_analytic[-1]
lin4_analytic   = np.append(lin4_analytic,Clin4_analytic_after_activator_peak(lin4_t0,beta_lin4,0,(n_post_oscillation-1)*dt,n_post_oscillation))

N_lin4_analytic_pattern = n_T*(N_analytic_before+N_analytic_during+N_analytic_after)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Launching the main computations # # # # # # # # # # # # # # # #

t = t0_23
cpt = 0 # variable counting the number of time the while loop is visited (see below)
cpt_max = len(L_Activator_with_post_oscillation)
integrity_boolean = 0
while (cpt<cpt_max-1):
    C_t = C_Processing(L_alpha,L_beta,gamma_lin4_lin14,C_t,dt)
    [lin4_dt , lin14_dt , LIN14_dt , nA] = C_t # refreshing the current concentration values
    nA     = L_Activator_with_post_oscillation[cpt] #refreshing the activator concentration value
    C_t[3] = nA

    lin4_numeric.append(lin4_dt)
    lin14_numeric.append(lin14_dt)
    LIN14_numeric.append(LIN14_dt)

    if(1-beta_lin14*dt - gamma_lin4_lin14*lin4_dt*dt<0 and integrity_boolean == 0):
        print(" ")
        print("WARNING")
        print("Be careful, beta_lin14 and gamma_lin4_lin14 are such that lin14 can be totally consumed.")
        print(" ")
        integrity_boolean = 1
    
    cpt+=1
    t+=dt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Color code for the plotting part  # # # # # # # # # # # # # # #

# NHR23         - cyan - c
# NHR85         - magenta - m
# Activator     - green - g
# lin4 numeric  - red - r
# lin4 analytic - blue - b
# lin14 numeric - orange = [255,165,0] (HAS TO BE SCALED BETWEEN 0 AND 1 => DIVIDE BY 255: [1,165/255,0])
# LIN14 numeric - yellow = y

# b: blue.
# g: green.
# r: red.
# c: cyan.
# m: magenta.
# y: yellow.
# k: black.
# w: white.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Plotting the results  # # # # # # # # # # # # # # # # # # # # #

## Rescaling the time axis :
rescaling_factor= (n_T*N_pattern+n_post_oscillation)/((N_pattern-1)*(Lt_whole_signal[-1]-Lt_whole_signal[0]))
Lt_whole_signal = (rescaling_factor+np.zeros(len(Lt_whole_signal)))*Lt_whole_signal

#graphe 1
#NHR23 et NHR85 et activateur sur une période

plt.subplot(231)
plt.title("Activator signal pattern")
plt.plot(t_pattern/max(t_pattern),Pattern_function_23/max(Pattern_function_23)  , c = 'k', label = 'NHR85 Pattern')
plt.plot(t_pattern/max(t_pattern),Pattern_function_85/(2.5*max(Pattern_function_85))  , c = 'r', label = 'NHR23 Pattern')
plt.plot(t_pattern/max(t_pattern),Pattern_function/(2.5*sqrt(max(Pattern_function_23)*max(Pattern_function_85))), c = 'm', label = 'Activator Pattern')
plt.xlabel('Scaled larval stage')
plt.ylabel("Norm. conc.")
#plt.ylabel('A(t)')
plt.legend()
plt.grid(True)

#graphe 2
#tout le signal sur toutes les périodes rescale (une unité de temps = un larval stage)

plt.subplot(232)
plt.title("Activator = product of concentrations")
plt.plot(Lt_whole_signal,L_Activator_with_post_oscillation/max(L_Activator_with_post_oscillation) , color = 'm')#[0,1,0], label = 'Whole activator signal in time')
plt.xlabel('Scaled larval stages')
plt.ylabel("Norm. conc.")
plt.legend()
plt.grid(True)

#graphe 3
#lin4 numeric et analytic sur toutes les périodes

plt.subplot(233)
#plt.title("Normalized lin4 miRNA concentration")
plt.plot(Lt_whole_signal,lin4_analytic/max(lin4_analytic) , c = 'b', label = 'analytical solution')
plt.plot(Lt_whole_signal,lin4_numeric/max(lin4_numeric)  , c = 'r', label = 'numerical solution')
plt.xlabel('Scaled larval stages')
plt.ylabel('lin4(t)')
plt.legend()
plt.grid(True)

#graphe 4
#lin 14 numeric

plt.subplot(234)
#plt.title("Normalized lin14 concentration in time")
plt.plot(Lt_whole_signal,np.array(lin14_numeric)/max(lin14_numeric), color = [1,165/255,0], label = 'numerical solution')
plt.xlabel('Scaled larval stages')
plt.ylabel('lin14(t)')
plt.legend()
plt.grid(True)
#graphe 5
#LIN14 numeric

plt.subplot(235)
#plt.title("Normalized LIN14 concentration in time")
plt.plot (Lt_whole_signal,np.array(LIN14_numeric)/max(LIN14_numeric), c = 'y', label = 'numerical solution')
plt.xlabel('Scaled larval stages')
plt.ylabel('LIN14(t)')
plt.legend()
plt.grid(True)

plt.show()
