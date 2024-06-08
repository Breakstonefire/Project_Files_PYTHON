# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing libraries # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import cm
from numpy import *
from scipy import signal
from math import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Definition of the functions # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def constant(y1,N):
    y = np.empty(0)
    for i in range(N):
        y = np.append(y,y1)
    return y

def Sigma_Computation(X,Y,boolean_interpolation):
    ymax = max(Y)
    ym2 = ymax/2
    N = len(Y)
    ind_max = 0
    test_left = 0
    test_right = 0
    for i in range(N-1):
        yi = Y[i]
        yii = Y[i+1]
        if (yi<=ym2 and yii>ym2):
            ileft = i
            test_left = 1
        if (yi==ymax):
            ind_max = i
        if (yi>ym2 and yii<=ym2):
            iright = i+1
            test_right = 1
    if(test_left==0 and test_right==0):
        return 0
    if(test_left==0):
        xil = X[ind_max]
        print(" ")
        print("WARNING")
        print("Be careful, the standard deviation can't be computed by the algorithm.")
        print(" ")
    elif(test_left==0):
        xiir = X[ind_max]
        print(" ")
        print("WARNING")
        print("Be careful, the standard deviation can't be computed by the algorithm.")
        print(" ")
    else:
        xil = X[ileft]
        xiir = X[iright+1]
    if(boolean_interpolation==0):
        dx = xiir-xil
    else:
        yil = Y[ileft]
        xiil = X[ileft+1]
        yiil = Y[ileft+1]
        xir = X[iright]
        yir = Y[iright]
        yiir = Y[iright+1]
        if(yi==yii):
            return 0
        xleft_interpolated = linear_interpolation(xil,yil,xiil,yiil,ymax)
        xright_interpolated = linear_interpolation(xir,yir,xiir,yiir,ymax)
        dx = xright_interpolated-xleft_interpolated    
    coeff = 1/(sqrt(2*np.log(2)))
    sigma = coeff*dx
    return sigma

def linear_interpolation(xi,yi,xii,yii,ymax):
    dx = xii-xi
    dy = yii-yi
    a = dy/dx
    b = yi-a*xi
    x_interpolated = (0.5*ymax-b)/a
    if a==0:
        x_interpolated = (0.5*ymax-b)
        x_interpolated*= dx
        x_interpolated/=dy
    return x_interpolated

def gaussian_part_of_solution(D,t0,t1,Nt,x0,x1,xmean,Nx):
    coef = 1/sqrt(4*pi*D)
    x = np.linspace(x0,x1,Nx,True)
    t = np.linspace(t0,t1,Nt,True)
    u = [[]]
    Sigma = [0 for i in range(Nt)]
    for k in range(Nt):
        u.append([0 for u in range(Nx)])
    for j in range(Nt):
        tj = t[j]
        ux = np.empty(0) #creating the list containing the values of concentration for one given t at every 'x'
        for i in range(Nx):
            xi = x[i]
            if(tj==0):
                uji = 0
            else:
                a = (xi-xmean)
                a/= 2*sqrt(D*tj)
                a*=a
                uji = np.exp(-a)    #computing the value of the 'u(xi,tj)' concentration
                uji/= sqrt(tj)      #computing the value of the 'u(xi,tj)' concentration
            ux = np.append(ux,uji)  # adding the value to the ut list
        ux*=coef
        u[j] = ux # adding the list of all concentration at one given time for every 'x' to the 2D list named 'u' # multiplying the 'u' list by the analytical coefficient in front of the varibale part of the solution
    Sigma = np.sqrt(2*D*t)
    return u,Sigma

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Setting the global parameters for the time and space coordinates  # # # # # # # # # # #

# ASSUMPTIONS :
### (1)
# D = kB*T/(6*pi*eta*R)
# kB = 1.380649 × 10-23 m2 kg s-2 K-1
# T belongs to this range of temperature: [15:25]°C => [288.15 : 303.15]°K
# eta is approximated by the viscosity of the cytoplasm, since cytoplasm is composed mainly of water, the dynamic viscosity of water can be taken for consideration which is, 8.90 × 10−4 Pa*s at about 25°C.
# R is the radius of a hormone wich can be approximated by a value belonging to [1.8 : 3.5]nm of diameter so ~1nm radius.

# Dtheoretical = 1.9 * 10^-10 ~=~ 2e-10 m²/s = 2e8 nm²/s

D = 2e-10 #ENTER THE VALUE IN m²/s

# Known parameters due to measurements - reasonable assumptions
# Unknown parameters
# Deduced parameters due to relations with other parameters

#### Space parameters - ENTER THE VALUES IN nm for computations safety and non-divide-by-zero operations ####
Nx      = 41                    # number of space points between x0 and x1, the x boundaries on the x axis
x0      = 0                     # The C.elegans is ~200 µm long in L1 growing to 1mm in adulthood, so here for the modelling we choose a L1 worm.
x1      = 2e-4                  # The C.elegans is ~200 µm long in L1 growing to 1mm in adulthood, so here for the modelling we choose a L1 worm.
x = np.linspace(x0,x1,Nx,True)
middle  = int((Nx-1)/2)
dx = (x1-x0)/(Nx-1)

#### Time parameters ####

t0      = 0                             # The initial time value corresponds here tp the begining of the L1 and will be set to 0
t1      = t0+60                         # The last time value corresponds to an "infinite" time value or at least to the duration of the egg to hatch, so approximately 9H.
                                        # The hatching happends ~9H after the egg formation => t1 ~=~ 9*3600 = 32400
Nt      = 1+int((Nx-1)*(Nx-1)*(t1-t0)/(0.92*90))             # number of time points between t0 and t1, the time boundaries on the t axis
dt      = (t1-t0)/(Nt-1)
t       = np.linspace(t0,t1,Nt,True)

#### Source and sink term profile ####
s0 = 0
s1 = 0             # The source term has to have the same time unit as D (the diffusion constant)
source_term = constant(s0,Nx)
source_term[middle] += s1

#### Initial concentration profile ####
U_middle = 100000
U_edges = 0

#### Computing useful constants for the next computations ####

x_middle = x[middle]
dx2 = dx**2
alpha = D*dt*(Nx-1)*(Nx-1)*pow(x1-x0,-2) #D*dt/dx2

boolean_interpolation = 1

#### Warning messages ####

if(Nx>10*sqrt(Nt/(t1-t0))):
    print(" ")
    print("WARNING")
    print("Be careful, Nx<=10*sqrt(Nt/(t1-t0)) isn't respected.")
    print(" ")

if alpha>0.5:
    print(" ")
    print("WARNING")
    print("Be careful, alpha = D*dt/dxdx",alpha,"is higher than 1/2 - Stability condition.")
    alpha = 0.4
    dt = alpha*dx*dx/D
    print("Now dt =",dt,",so alpha = 0.4.")
else:
    print(" ")
    print("Stability condition validated : alpha = D*dt/dxdx =",alpha,"is lower than 1/2.")
    print(" ")
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Creating the result lists # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

U               = np.zeros(Nx)
U_middle_t          = np.zeros(Nt)
#U_middle_t_analytic = U_middle/np.sqrt(np.pi*D*t)
laplacien_dx2   = np.zeros(Nx)
Flux            = np.zeros(Nx)
Flux_middle_t       = np.zeros(Nt)
Flux_edges_t        = np.zeros(Nt)
Flux_edges_left_t   = np.zeros(Nt)
Flux_edges_right_t  = np.zeros(Nt)

U[0]  = U_edges
U[middle] = U_middle
U[Nx-1]= U_edges

U_initial = list(U)

Sigma_numeric = []

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Launching the main computations # # # # # # # # # # # # # # # #
xmean = (x0+x1)/2
Gauss_solution_X_T , Sigma_analytic = gaussian_part_of_solution(D,t0,t1,Nt,x0,x1,xmean,Nx)

U_all_x = [[0 for u in range(Nx)] for i in range(Nt)]
U_all_x[0] = U_initial
Sigma_numeric.append(0)
for i in range(1,Nt):
    for k in range(1,Nx-1):
        laplacien_dx2[k] = U[k+1]-2*U[k]+U[k-1]
        Flux[k] = -D*(U[k+1]-U[k])/dx
    for k in range(1,Nx-1):
        U[k] = U[k] + alpha*laplacien_dx2[k] + source_term[k]*dt
    for k in range(0,len(U),1):
        uk = U[k]
        U_all_x[i][k] = uk
    sigma = Sigma_Computation(x,U,boolean_interpolation)
    if(sigma<0):
        sigma = Sigma_numeric[len(Sigma_numeric)-1]
        print(" ")
        print("WARNING")
        print("Be careful, the standard deviation is negative so it has been corrected to the previous positive value.")
        print(" ")
    Sigma_numeric.append(sigma)
    U_middle_t[i]           = U[middle]
    Flux_middle_t[i]        = Flux[middle]
    Flux[0]                 = -D*(U[1]-U[0])/dx
    Flux[Nx-1]              = -D*(U[Nx-1]-U[Nx-2])/dx
    Flux_edges_left_t[i]    = Flux[0]
    Flux_edges_right_t[i]   = Flux[Nx-1]

###################################################################################################
#### THE FOLLOWING BLOCK CAN BE COPY-PASTE IN EVERY CODE TO PLOT K curves composed of N VALUES ####
#RED TO BLACK

# K =     # value corresponding to the number of 'y(x)' curves to plot
# N =     # value corresponding to the number of points to plot for each curve (= len(x))

Lcolors = [[]]
for i in range(1276):
    Lcolors.append([0,0,0])

# Theoretically the number of triplet-color-values in Lcolors is equal to 5*(2^8-1) = 1275
# INCREASING THE VALUE
index = 0
for green in range(0,256,1):
    rgb = [255,green,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for red in range(254,-1,-1):
    rgb = [red,255,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# INCREASING THE VALUE
for blue in range(1,256,1):
    rgb = [0,255,blue]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for green in range(254,-1,-1):
    rgb = [0,green,255]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for blue in range(254,0,-1):
    rgb = [0,0,blue]
    Lcolors[index] = rgb
    index+=1
# DELETING THE DOUBLE VALUES IN THE END
for i in range(len(Lcolors)-1,0,-1):
    l1 = Lcolors[i]
    l2 = Lcolors[i-1]
    if l1==l2:
        print(" ")
        print("WARNING when computing the colors for plots: doublon at index",i)
        print(" ")
        del Lcolors[i]
del Lcolors[-1]
len_col = len(Lcolors)
# Theoretically the number of triplet-color-values in Lcolors is equal to 5*(2^8-1) = 1275

# Then we normalize the values of Lcolors from [0:255] to [0:1]
for i in range(len(Lcolors)):
    for j in range(len(Lcolors[0])):
        Lcolors[i][j]/=255

Lcolors.reverse()

#### END OF THE GENERIC BLOCK ####
##################################
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Plotting the results  # # # # # # # # # # # # # # # # # # # # #

# Create a figure
fig = plt.figure()

##############################################################################

ax = fig.add_subplot(331) # Add a subplot
plt.subplot(331)
title = "Source term shape in space"
plt.title(title)
plt.plot (x,source_term, c = 'g', label = 'Source term')
plt.xlabel('Worm axis (m)')
plt.ylabel('y')
plt.legend()
plt.grid(True)

##############################################################################

ax = fig.add_subplot(332)
plt.subplot(332)
title = "Initial concentration in space"
plt.title(title)
#plt.plot (x,U_initial, c = 'k', label = 'Initial concentration')
plt.vlines(x=x_middle, ymin=0, ymax=U_middle, color = 'k', label = 'Initial concentration', linestyle='-')
plt.xlim([x0,x1])
plt.xlabel('Worm axis (m)')
plt.ylabel('C(x;t=0)')
plt.legend()
plt.grid(True)

##############################################################################

ax = fig.add_subplot(333)
plt.subplot(333)
title = "Concentration in the center in time"
plt.title(title)
plt.plot (t,U_middle_t,             color = 'r', label = 'C_numeric')
#plt.plot (t,U_middle_t_analytic,    color = 'b', label = 'C_analytic')
plt.xlabel('Time axis (s)')
plt.ylabel('C(center;t)')
plt.legend()
plt.grid(True)

##############################################################################

ax = fig.add_subplot(334)
plt.subplot(334)
title = "Concentrations at stationary state - Numeric"
colorbar_legend = "Rising time: BLACK -> RED"

plt.title(title)
K = len(U_all_x)
for k in range(K):
    yk = U_all_x[k]
    if(K<len_col):
        ktild = int(k*len_col/K)
    elif(K==len_col):
        ktild = k
    elif(K>len_col):
        ktild = int(k*len_col/K)
    ck = Lcolors[ktild]
    ymax = max(yk)
    if ymax==0:
        if k==0:
            plt.vlines(x=x_middle, ymin=0, ymax=1, color=[0,0,0], linestyle='-')
        else:
            plt.vlines(x=x_middle, ymin=0, ymax=1, color=ck, linestyle='-')
    else:
        for i in range(len(yk)):
            yk[i]/=ymax
        plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
plt.vlines(x=x_middle, ymin=0, ymax=1, color=Lcolors[0], linestyle='-')
plt.xlabel('Worm axis (m)')
plt.ylabel('Norm_C(x,t=∞)')
plt.legend()
plt.grid(True)
c = np.array([t0,t1])
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
cmap.set_array([])
cbar = plt.colorbar(cmap)
cbar.set_label(colorbar_legend, rotation=90)

##############################################################################

ax = fig.add_subplot(335)
title           = "Concentrations at stationary state - Analytic"
subtitle        = "Gaussian shape"
colorbar_legend = "Rising time: BLACK -> RED"

plt.title(title)
plt.suptitle(subtitle)
K = len(Gauss_solution_X_T)
ylim_min = 0
for k in range(K):
    yk = Gauss_solution_X_T[k]
    if(K<len_col):
        ktild = int(k*len_col/K)
    elif(K==len_col):
        ktild = k
    elif(K>len_col):
        ktild = int(k*len_col/K)
    ck = Lcolors[ktild]
    ymax = max(yk)
    if ymax==0:
        if k==0:
            plt.vlines(x=x_middle, ymin=0, ymax=1, color=[0,0,0], linestyle='-')
        else:
            #plt.vlines(x=x_middle, ymin=0, ymax=1, color=ck, linestyle='-')
            print("WARNING")
            print("Be careful, the analytical solution at time",t[min(k,K-2)],"cannot be correctly computed.")
    else:
        for i in range(len(yk)):
            yk[i]/=ymax
        #ylim_min = min(yk)
        plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
plt.xlabel('Worm axis (m)')
plt.ylabel('Norm_C(x,t=∞)')
plt.ylim([ylim_min,1])
plt.legend()
plt.grid(True)
c = np.array([t0,t1])
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
cmap.set_array([])
cbar = plt.colorbar(cmap)
cbar.set_label(colorbar_legend, rotation=90)

##############################################################################

ax = fig.add_subplot(336)
plt.subplot(336)
title = "Standard deviations in time"
subtitle = ""
plt.title(title)
plt.suptitle(subtitle)
plt.plot(t,Sigma_analytic, c = 'b', marker = '', label = 'analytical solution')
plt.plot(t,Sigma_numeric , c = 'r', marker = '', label = 'numerical solution')
plt.xlabel('t')
plt.ylabel('Sigma(m)')
plt.legend()
plt.grid(True)

##############################################################################

ax = fig.add_subplot(337)
plt.subplot(337)
title = "Flux in space at stationary state"
plt.title(title)
plt.plot (x[1:len(x)-2],Flux[1:len(x)-2], c = 'b', label = 'Flux')
plt.xlabel('Worm axis (m)')
plt.ylabel('Flux(m/s)')
plt.legend()
plt.grid(True)


##############################################################################

plt.show()
