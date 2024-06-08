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

def gaussian_part_of_solution(D,t0,t1,Nt,x0,x1,xmean,Nx):
    coef = 1/sqrt(4*pi*D)
    x = np.linspace(x0,x1,Nx,True)
    t = np.linspace(t0,t1,Nt,True)
    u = [[]]
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
                uji = np.exp(-((xi-xmean)**2)/(4*D*tj))/sqrt(tj) #computaing the value of the 'u(xi,tj)' concentration
            ux = np.append(ux,uji) # adding the value to the ut list
        ux*=coef
        u[j] = ux # adding the list of all concentration at one given time for every 'x' to the 2D list named 'u' # multiplying the 'u' list by the analytical coefficient in front of the varibale part of the solution
    return u

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Setting the global parameters for the time and space coordinates  # # # # # # # # # # #

D = 500

Nx       = 101          # number of space points between x0 and x1, the x boundaries on the x axis
x0      = -10           #
x1      = 10            #
x = np.linspace(x0,x1,Nx,True)
dx = x[1]-x[0]

Nt      = 501           # number of time points between t0 and t1, the time boundaries on the t axis
t0      = 0
dt      = 1e-5
t1      = (Nt-1)*dt
t       = np.linspace(t0,t1,Nt,True)

#### Computing useful constants for hte next computations ####

middle  = int((Nx-1)/2) # 
dx2 = dx**2
alpha = D*dt/dx2

if alpha>0.5:
    print("Be careful, alpha = D*dt/dxdx is higher than 1/2 - Stability condition.")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Creating the result lists # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

U               = np.zeros(Nx)
U_middle_t          = np.zeros(Nt)
laplacien_dx2   = np.zeros(Nx)
Flux            = np.zeros(Nx)
Flux_middle_t       = np.zeros(Nt)
Flux_edges_t        = np.zeros(Nt)
Flux_edges_left_t   = np.zeros(Nt)
Flux_edges_right_t  = np.zeros(Nt)

s0 = 0
s1 = 0
source_term = constant(s0,Nx)
source_term[middle] += s1

U_middle = 10
U_edges = 0

U[0]  = U_edges
U[middle] = U_middle
U[Nx-1]= U_edges

U_initial = list(U)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Launching the main computations # # # # # # # # # # # # # # # #
xmean = (x0+x1)/2
Gauss_solution_X_T = gaussian_part_of_solution(D,t0,t1,Nt,x0,x1,xmean,Nx)

U_all_x = [[]]
for k in range(Nt):
    U_all_x.append([0 for u in range(Nx)])
for i in range(Nt):
    for k in range(1,Nx-1):
        u = U[k+1]-2*U[k]+U[k-1]
        laplacien_dx2[k] = u
        Flux[k]     = -D*(U[k+1]-U[k])/dx
    for k in range(1,Nx-1):
        U[k] += alpha*laplacien_dx2[k]+dt*source_term[k]
    U_all_x[i] = U
    U_middle_t[i]           = U[middle]
    Flux_middle_t[i]        = Flux[middle]
    Flux[0]                 = -D*(U[1]-U[0])/dx
    Flux[Nx-1]              = -D*(U[Nx-1]-U[Nx-2])/dx
    Flux_edges_left_t[i]    = Flux[0]
    Flux_edges_right_t[i]   = Flux[Nx-1]

###################################################################################################
#### THE FOLLOWING BLOCK CAN BE COPY-PASTE IN EVERY CODE TO PLOT K curves composed of N VALUES ####

Lcolors = [[]]
for i in range(1531):
    Lcolors.append([0,0,0])

# Theoretically the number of triplet-color-values in Lcolors is equal to 6*(2^8-1) = 1530
# DECREASING THE VALUE
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
# DECREASING THE VALUE
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
for red in range(1,256,1):
    rgb = [red,0,255]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for blue in range(254,0,-1):
    rgb = [255,0,blue]
    Lcolors[index] = rgb
    index+=1
# DELETING THE DOUBLE VALUES IN THE END
for i in range(len(Lcolors)-1,0,-1):
    l1 = Lcolors[i]
    l2 = Lcolors[i-1]
    if l1==l2:
        print("doublon at index",i)
        del Lcolors[i]
del Lcolors[-1]
len_col = len(Lcolors)
# Theoretically the number of trilpet-color-values in Lcolors is equal to 6*(2^8-1) = 1530

# Then we normalize the values of Lcolors from [0:255] to [0:1]
for i in range(len(Lcolors)):
    for j in range(len(Lcolors[0])):
        Lcolors[i][j]/=255

#### END OF THE GENERIC BLOCK ####
##################################
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Plotting the results  # # # # # # # # # # # # # # # # # # # # #

# Create a figure
fig = plt.figure()

##############################################################################
# Add a subplot
ax = fig.add_subplot(331)
plt.subplot(331)
title = "Source term shape in space"
plt.title(title)
plt.plot (x,source_term, c = 'g', label = 'Source term')
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.grid(True)

##############################################################################
ax = fig.add_subplot(332)
plt.subplot(332)
title = "Initial conditions at boundaries"
plt.title(title)
plt.plot (x,U_initial, c = 'k', label = 'Initial concentration')
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.grid(True)

ax = fig.add_subplot(333)

plt.subplot(333)
title = "Concentration in middle-space in time"
plt.title(title)
plt.plot (t,U_middle_t, color = 'b', label = 'Concentration middle-space')
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.grid(True)

##############################################################################
ax = fig.add_subplot(334)
plt.subplot(334)
title = "Concentrations at stationary state - Numeric"
plt.title(title)
#plt.plot (x,U, color = (0.6,0.2,0.9) , marker = '', label = 'numerical solution')
K = len(U_all_x)
for k in range(K):
    yk = U_all_x[k]
    if(K<len_col):
        ktild = int(k*len_col/K)
    elif(K==len_col):
        ktild = k
    elif(K>len_col):
        if(k<len_col):
            ktild = k
        if(k>=len_col):
            ktild = int(10*((k/len_col)-int(k/len_col)))
    ck = Lcolors[ktild]
    ymax = max(yk)
    if ymax==0:
        ymax=1
    for i in range(len(yk)):
        yk[i]/=ymax
    plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.grid(True)

##############################################################################
ax = fig.add_subplot(335)
c = np.arange(0, 2)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
cmap.set_array([])
fig.colorbar(cmap, ticks=c)

#fig = plt.subplot(335)
title = "Concentrations at stationary state - Analytic"
subtitle = "Gaussian shape"
plt.title(title)
plt.suptitle(subtitle)
K = len(Gauss_solution_X_T)
for k in range(K):
    yk = Gauss_solution_X_T[k]
    if(K<len_col):
        ktild = int(k*len_col/K)
    elif(K==len_col):
        ktild = k
    elif(K>len_col):
        if(k<len_col):
            ktild = k
        if(k>=len_col):
            ktild = int(10*((k/len_col)-int(k/len_col)))
    ck = Lcolors[ktild]
    ymax = max(yk)
    if ymax==0:
        ymax=1
    for i in range(len(yk)):
        yk[i]/=ymax
    plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
plt.xlabel('x')
plt.ylabel('Normalized concentration')
plt.legend()
plt.grid(True)
##############################################################################
##
##
##plt.subplot(335)
##title = "Concentrations at stationary state - Analytic"
##subtitle = "Gaussian shape"
##plt.title(title)
##plt.suptitle(subtitle)
##K = len(Gauss_solution_X_T)
##for k in range(K):
##    yk = Gauss_solution_X_T[k]
##    if(K<len_col):
##        ktild = int(k*len_col/K)
##    elif(K==len_col):
##        ktild = k
##    elif(K>len_col):
##        if(k<len_col):
##            ktild = k
##        if(k>=len_col):
##            ktild = int(10*((k/len_col)-int(k/len_col)))
##    ck = Lcolors[ktild]
##    ymax = max(yk)
##    if ymax==0:
##        ymax=1
##    for i in range(len(yk)):
##        yk[i]/=ymax
##    plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
##yk = Gauss_solution_X_T[int(len(Gauss_solution_X_T)/2-1)]
##yk = np.array(yk)
##print(yk)
##colors = plt.cm.hsv(yk / float(max(yk)))
##sm = plt.cm.ScalarMappable(cmap=plt.cm.hsv, norm=plt.Normalize(vmin=0, vmax=1))
##sm._A = []
##plt.colorbar(sm)
##plt.bar(range(len(yk)), yk, color = colors)
##
##plt.xlabel('x')
##plt.ylabel('Normalized concentration')
##plt.legend()
##plt.grid(True)

##############################################################################
##############################################################################
##cmap = mpl.cm.get_cmap('jet', K)
##plt.subplot(336)
##title = "Concentrations at stationary state - Analytic"
##subtitle = "Gaussian shape"
##plt.title(title)
##plt.suptitle(subtitle)
##K = len(Gauss_solution_X_T)
##for k in range(K):
##    yk = Gauss_solution_X_T[k]
##    if(K<len_col):
##        ktild = int(k*len_col/K)
##    elif(K==len_col):
##        ktild = k
##    elif(K>len_col):
##        if(k<len_col):
##            ktild = k
##        if(k>=len_col):
##            ktild = int(10*((k/len_col)-int(k/len_col)))
##    ck = Lcolors[ktild]
##    ymax = max(yk)
##    if ymax==0:
##        ymax=1
##    for i in range(len(yk)):
##        yk[i]/=ymax
##    plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
##plt.xlabel('x')
##plt.ylabel('Normalized concentration')
##plt.legend()
##plt.grid(True)

##c = np.arange(1, 100)
##norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
##cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
##cmap.set_array([])
##fig, ax = plt.subplots(dpi=100)
##fig.colorbar(cmap, ticks=c)
##plt.show()

##min, max = (-40, 30)
##step = 10
### Setting up a colormap that's a simple transtion
##mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
### Using contourf to provide my colorbar info, then clearing the figure
##Z = [[0,0],[0,0]]
##levels = range(min,max+step,step)
##CS3 = plt.contourf(Z, levels, cmap=mymap)
##plt.clf()
###plt.colorbar(ax,CS3) # using the colorbar info I got from contourf

##############################################################################
##############################################################################
##############################################################################

##plt.subplot(337)
##title = "Flux in space at stationary state"
##plt.title(title)
##plt.plot (x[1:len(x)-2],Flux[1:len(x)-2], c = 'b', label = 'Flux')
##plt.xlabel('')
##plt.ylabel('')
##plt.legend()
##plt.grid(True)
##
##plt.subplot(338)
##title = "Flux in time"
##plt.title(title)
##plt.plot (t,Flux_middle_t, c = 'b', label = 'Flux middle-space')
##plt.xlabel('')
##plt.ylabel('')
##plt.legend()
##plt.grid(True)
##
##plt.subplot(339)
##title = "Flux in time"
##plt.title(title)
##plt.plot (t,Flux_edges_left_t, c = 'c', label = 'Flux left edge')
##plt.plot (t,Flux_edges_right_t, c = 'm', label = 'Flux right edge')
##plt.xlabel('')
##plt.ylabel('')
##plt.legend()
##plt.grid(True)
##

plt.show()
