# This code is describing how a hormone diffuses in a worm.
# We compute with the 1D diffusion equation the values of the concentration in time and spac.
# We plot the results in a layout where different plots are showed (flux in time, flux in space, 
# Here the viscosity of the nucleus which is different from the one of the cytoplasm where hromones diffuse too is taken into account.
# At the end, it seems that the hormones are having an x² shape in terms of concentration at stationnary state (at the end of the model)
# This shape seems to be reached at only one minute after diffusion started.
# We also compute the error compared to an "x²" curve (with the good coefficients).

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing libraries # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
import numpy as np
# import math
# import subprocess
# from matplotlib import cm
# from numpy import *
# from scipy import signal
from math import sqrt


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Definition of the functions # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def constant(y1,N): # This function is useful for creating an array with constant values
    y = np.empty(0)
    for i in range(N):
        y = np.append(y,y1)
    return y

def derivative(X,Y): # This function is useful to have the derivative values of a function in an array
    N = len(Y)
    dY = []
    for i in range(N-1):
        xi = X[i]
        xii = X[i+1]
        yi = Y[i]
        yii = Y[i+1]
        dx = xii-xi
        dy = yii-yi
        dY.append(dy/dx)
    return dY

def LeastSquaresMethod(Ytheor,Yexp): # This function is useful to compute the error of a curve compared to another theoretical curve with the least square method
    N = len(Ytheor)
    eps = 0
    for i in range(N):
        dy = Ytheor[i]-Yexp[i]
        eps+=dy*dy
    return eps

def Sigma_Computation(X,Y,boolean_interpolation):
    ymax = max(Y)           # we take back the value of the highest peak of the curve
    ym2 = ymax/2            # we are looking for sigma (which is the half-height gap or the "the standard deviation" for gaussians for example) so we need to know where the curve is equal to the half value of the peak
    N = len(Y)
    ind_max = 0             # this variable corresponds to the index of the maximum value of the curve in the Y array
    test_left = 0           # this boolean will say whether we found a value for the left x that is needed in the sigma computation 
    test_right = 0          # this boolean will say whether we found a value for the right x that is needed in the sigma computation
    for i in range(N-1):    # we look at every Y values to find the indexes that are needed to identify the x-associated values for the sigma computation
        yi = Y[i]
        yii = Y[i+1]
        if (yi<=ym2 and yii>ym2):   # this condition ensures that we are currently observing the y value that is at the left of the peak and that is the most close to the half-peak value
            ileft = i               # we then save the value of the index to know after the x value associated
            test_left = 1
        if (yi==ymax):              # this condition ensures that we are currently observing the y value that is the exact half-peak value
            ind_max = i             # we then save the value of the index to know after the x value associated
        if (yi>ym2 and yii<=ym2):   # this condition ensures that we are currently observing the y value that is at the right of the peak and that is the most close to the half-peak value
            iright = i+1            # we then save the value of the index to know after the x value associated
            test_right = 1
    if(test_left==0 and test_right==0):
        return 0
    if(test_left==0):
        xil = X[ind_max]
##        print(" ")
##        print("WARNING")
##        print("Be careful, the standard deviation can't be computed by the algorithm.")
##        print(" ")
    elif(test_right==0):
        xiir = X[ind_max]
##        print(" ")
##        print("WARNING")
##        print("Be careful, the standard deviation can't be computed by the algorithm.")
##        print(" ")
    else:
        xil = X[ileft]
        xiir = X[iright]
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

def linear_interpolation(xi,yi,xii,yii,ymax): # This function computes standard deviation of a quasi-gaussian (or -x² shape) curve using a linear interpolation of the curve at the given points given in parameters
    dx = xii-xi
    dy = yii-yi
    a = dy/dx
    b = yi-a*xi
    if a==0:
##        print(" ")
##        print("WARNING")
##        print("Be careful, the standard deviation can't be computed by the algorithm.")
##        print(" ")
        x_interpolated = 0
    else:
        x_interpolated = (0.5*ymax-b)/a
    return x_interpolated

def invert_color(Color): # This function inverts the color of a chosen curve whose color is already set.
                            # It is useful when you want to put a black graph on a black font instead of a white graph on a black font and the inverse is true too.
    Color_inverted = np.array(list(Color))
    if(max(Color_inverted)>1):
        for i in range(len(Color_inverted)):
            Color_inverted/=255
    for i in range(len(Color_inverted)):
        x = Color_inverted[i]
        if(x<0.5):
            delta = 0.5-x
            x = 0.5+delta
        elif(x>0.5):
            delta = x-0.5
            x = 0.5-delta
        Color_inverted[i] = x
    return Color_inverted

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Setting the global parameters for the time and space coordinates  # # # # # # # # # # #

# ASSUMPTIONS :
### (1)
# D = kB*T/(6*pi*eta*R)
# kB = 1.380649 × 10-23 m2 kg s-2 K-1
# T belongs to this range of temperature: [15:25]°C => [288.15 : 303.15]°K
# eta is approximated by the viscosity of the cytoplasm, since cytoplasm is composed mainly of water, the dynamic viscosity of water can be taken for consideration which is, 8.90 × 10−4 Pa*s at about 25°C.
# R is the radius of a hormone wich can be approximated by a value belonging to [1.8 : 3.5]nm of diameter so ~1nm radius.

Nucleus_size = 6*1e-6 #the nucleus diameter in mammalians is approx. 6µm

# Dtheoretical = 1.9 * 10^-10 ~= 2e-10 m²/s = 2e8 nm²/s
D_cytoplasm = 1.9e-10 #ENTER THE VALUE IN m²/s
D_nucleus = 2*D_cytoplasm
D = max(D_cytoplasm,D_nucleus)

# Known parameters due to measurements - reasonable assumptions
# Unknown parameters
# Deduced parameters due to relations with other parameters

#### Space parameters - ENTER THE VALUES IN nm for computations safety and non-divide-by-zero operations ####
Nx      = 300                       # number of space points between x0 and x1, the x boundaries on the x axis
x0      = 0                         # The C.elegans is ~200 µm long in L1 growing to 1mm in adulthood, so here for the modelling we choose a L1 worm.
x1      = 2e-4                      # The C.elegans is ~200 µm long in L1 growing to 1mm in adulthood, so here for the modelling we choose a L1 worm.
dx = (x1-x0)/(Nx-1)                 # this is the step parameter in space
Nx_in_the_nucleus = int(Nx*Nucleus_size/(x1-x0))
x = np.linspace(x0,x1,Nx,True)

#### Time parameters ####
t0 = 0
t1 = t0 + 20

Nt = 1 + int(2*D*(t1-t0)/(dx*dx))
dt      = (t1-t0)/(Nt-1)
t       = np.linspace(t0,t1,Nt,True)

#### Source and sink term profile ####
s0 = 0                              # COMMENT
s1 = 10000                          # The source term has to have the same time unit as D (the diffusion constant) times the unit of the "amount" or "concentration" of hormones
source_term = constant(s0,Nx)

#### Worm - cell - nucleus parameters
Ncell = 20
i_step = int((Nx-1)/(Ncell+1))      # this is the value that correspons to the number of index to jump to go from a cell to the the next one (the distance between two cells is obviously i_step*dx)
#print(i_step)

#### Adding the source terms of every cell to the source_term array and computing the array of Diffusion constant in space
Diffusion_array = np.ones(Nx)
Diffusion_array*= D_cytoplasm

j_step = 1 + int(Nucleus_size/dx)
#print(j_step)
j_step_corr = -int((1+j_step)/2)
for i in range(Ncell):
    source_term[i_step+i*i_step]+= s1
    for j in range(1,j_step+1):
        Diffusion_array[i_step+i*i_step+j_step_corr+j] = D_nucleus

#### Choosing the rate of curves to plot for the color plot (less curves to plot => lower time to charge and show the result)
### You can either choose directly the number of curves to show in general [0 ; Nt] or choose the number of curves that you don't show for one curve you show (this is cpt_max)
Ncurves_to_show = 250   # If you don't want the plotting part to take too much time after all the results are computed, you can set a reasonable value of 100 to 300 curves to show in total.
                        # This still allows you to see the evolution of the concentration in space and time during all he simulation.
cpt_max = int(Nt/Ncurves_to_show) # This corresponds to the number of curve computed to show only one curve among them (the last one from the 'cpt_max' curves will be shown)

#### Initial concentration profile ####
U_middle = 0
U_edges  = 0

#### Computing useful constants for the next computations ####
middle  = int((Nx-1)/2)             # computing the index corresponding to the middle of the x axis
x_middle = x[middle]                # computing the x value corresponding to the middle of the x axis
alpha = D*dt*pow((Nx-1)/(x1-x0),2)  # alpha = D*dt/dx² = D*dt*(Nx-1)²/(x1-x0)²

boolean_interpolation = 0           # this boolean tells whether there will be a linear interpolation in the computation of every standard deviation values in the function "Sigma_Computation"

## Plots parameters
# Setting the different sizes for all plots
DefaultFigureWidth    = 5
DefaultFigureHeight   = 3
DefaultDPIResolution  = 600 # 800 is enough for detailed plots on a PC screen
DefaultFaceColor      = [1 , 1 , 1] # Spyder Black Background color is [25/100 , 35/100 , 45/100]
DefaultEdgeColor      = [0.1 , 0.1 , 0.1] # Spyder Black Background color is [25/100 , 35/100 , 45/100]
DefaultLayout         = 'constrained' # 'constrained', 'compressed', 'tight'

#### Warning messages ####

if(Nx_in_the_nucleus<10):
    print(" ")
    print("WARNING - 1")
    print("Be careful, you discretized the nuclei with less than 10 points (only ",Nx_in_the_nucleus," points).")
    print("Increase the number of points in x-axis to ",int(10*(x1-x0)/Nucleus_size)," points.")
    print(" ")
else:
    print(" ")
    print("The nuclei is discretized with ",Nx_in_the_nucleus," points.")
    print(" ")

if alpha>0.5:
    print(" ")
    print("WARNING - 2")
    print("Be careful, alpha = D*dt/dxdx",alpha,"is higher than 1/2 in the nucleus - Stability condition.")
    alpha = 0.49
    dt = alpha*dx*dx/D
    print("Now dt =",dt,",so alpha = 0.49")
else:
    print(" ")
    print("Stability condition validated : alpha = D*dt/dxdx =",alpha,"is lower than 1/2.")
    print(" ")

tcomputing = 0.007*(t1-t0)*Nx
tc_h = int(tcomputing/3600)
tc_mn = int((tcomputing-tc_h*3600)/60)
tc_sec = int(tcomputing - 60*tc_mn - 3600*tc_h)

tshowing = 150*1e-7*tcomputing*Nt/cpt_max
ts_h = int(tshowing/3600)
ts_mn = int((tshowing-ts_h*3600)/60)
ts_sec = int(tshowing - 60*ts_mn - 3600*ts_h)
print(" ")
print("The process will take approximately:")
print(" - ",tc_h,"h ",tc_mn,"mn ",tc_sec," sec to compute the results.")
print(" - ",ts_h,"h ",ts_mn,"mn ",ts_sec," sec to show the results.")
print(" ")

str_continue = "yes"
# str_continue = str(input("Do you still want to continue? (y/n)"))

if(str_continue=="n"):
    exit()
else:
    print(" ")
    print("Starting the process.")
    print(" ")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
LSM                 = []
LSM_time            = []

U[0]  = U_edges
U[middle] = U_middle
U[Nx-1]= U_edges

U_initial = list(U)

Sigma_numeric = []

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Launching the main computations # # # # # # # # # # # # # # # #

xmean = (x0+x1)/2

print(" ")
print("1 - Computing the numerical solution.")
print(" ")
U_all_x = [[0 for u in range(Nx)] for i in range(Nt)]
U_all_x[0] = U_initial
Sigma_numeric.append(0)
sigma = 0.0001
for i in range(1,Nt):
    for k in range(1,Nx-1):
        laplacien_dx2[k] = U[k+1]-2*U[k]+U[k-1]
        D = Diffusion_array[k-1]
        Flux[k] = -D*(U[k+1]-U[k])/dx
    for k in range(1,Nx-1):
        D = Diffusion_array[k-1]
        alpha = D*dt*pow((Nx-1)/(x1-x0),2)
        U[k]+=alpha*laplacien_dx2[k] + source_term[k]*dt
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

###########################################################################################################
#### THE FOLLOWING BLOCK CAN BE COPY-PASTE IN EVERY CODE TO PLOT K colored curves composed of N VALUES ####
#RED TO BLACK

# K =     # value corresponding to the number of 'y(x)' curves to plot
# N =     # value corresponding to the number of points to plot for each curve (= len(x))

Lcolors = [[]]
for i in range(1276):
    Lcolors.append([0,0,0])

# Theoretically the number of triplet-color-values in Lcolors is equal to 5*(2^8-1) = 1275

# Color is RED
# INCREASING THE VALUE of GREEN
index = 0
for green in range(0,256,1):
    rgb = [255,green,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)

# Color is YELLOW
# DECREASING THE VALUE of RED
for red in range(254,-1,-1):
    rgb = [red,255,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)

# Color is GREEN
# INCREASING THE VALUE of BLUE
for blue in range(1,256,1):
    rgb = [0,255,blue]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)

# Color is CYAN
# DECREASING THE VALUE of GREEN
for green in range(254,-1,-1):
    rgb = [0,green,255]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)

# Color is BLUE
# DECREASING THE VALUE of BLUE
for blue in range(254,0,-1):
    rgb = [0,0,blue]
    Lcolors[index] = rgb
    index+=1

# Color is BLACK
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

#################################################################################
##Parameters for the plots

µm_scaling = 1
normalize = 1

#Enter here the colors you want to have with or without any 'color invert' tool after to put the figures in some slides and/or presentations
color_fig1 = [0.5,0.5,0.5] # default color - grey
color_fig2 = [1,1,0] #corresponds to yellow
color_fig3 = [1,1,0] #corresponds to yellow
color_fig4 = [0.5,0.5,0.5] # default color - grey
color_fig5 = [0.5,0.5,0.5] # default color - grey
color_fig6 = [0.5,0.5,0.5] # default color - grey
color_fig7 = [0.5,0.5,0.5] # default color - grey
color_fig8 = [0.5,0.5,0.5] # default color - grey
color_fig9 = [0.5,0.5,0.5] # default color - grey

will_the_colors_be_inverted_after = 0 #enter 0 for no or 1 for yes if you plan to invert the colors of the figures after

###############################################################################
##Eventually Normalizing / Scaling / Changing colors / 

space_unit = "(m)"
if (µm_scaling==1):
    print(" ")
    print("2 - Scaling the results.")
    print(" ")
    space_unit = "(μm)"
    x = x*1e6
    x0 = x0*1e6
    x1 = x1*1e6
    x_middle = x_middle*1e6
    Sigma_numeric = np.array(Sigma_numeric)*1e6

if (normalize==1):
    print(" ")
    print("3 - Normalizing the results.")
    print(" ")
    source_term = source_term/max(source_term)

if (will_the_colors_be_inverted_after==1):
    print(" ")
    print("4 - Inverting the colors of the plots.")
    print(" ")
    color_fig1 = invert_color(color_fig1)
    color_fig2 = invert_color(color_fig2)
    color_fig3 = invert_color(color_fig3)
    color_fig4 = invert_color(color_fig4)
    color_fig5 = invert_color(color_fig5)
    color_fig6 = invert_color(color_fig6)
    color_fig7 = invert_color(color_fig7)
    color_fig8 = invert_color(color_fig8)
    color_fig9 = invert_color(color_fig9)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Plotting the results  # # # # # # # # # # # # # # # # # # # # #

print(" ")
print("5 - Plotting the results.")
print(" ")

# Create a figure
fig = plt.figure(figsize = (DefaultFigureWidth , DefaultFigureHeight) , dpi = DefaultDPIResolution , facecolor = DefaultFaceColor , edgecolor = DefaultEdgeColor , layout = DefaultLayout)

##############################################################################

ax = fig.add_subplot(331) # Add a subplot
plt.subplot(331)
title = "Source distribution in space"
#plt.title(title)
plt.plot(x,source_term, color = color_fig1, label = 'Source term')
plt.xlabel('Worm axis'+space_unit)
plt.ylabel('Source term')
plt.legend(loc='upper right')
plt.grid(True)

################################################################################

ax = fig.add_subplot(332)
plt.subplot(332)
title = "Initial concentration in space"
#plt.title(title)
plt.plot (x,U_initial, color = color_fig2, label = 'Initial concentration')
#plt.vlines(x=x_middle, ymin=0, ymax=U_middle, color = 'k', label = 'C(t=0)', linestyle='-')
plt.xlim([x0,x1])
plt.xlabel('Worm axis'+space_unit)
plt.ylabel('Hormone concentration')
plt.legend(loc='upper right')
plt.grid(True)

##############################################################################

ax = fig.add_subplot(333)
plt.subplot(333)
title = "Norm. hormone conc. in the center in time"
#plt.title(title)
plt.plot (t, np.array(U_middle_t)/max(U_middle_t), color = color_fig3, label = '[H](x=middle;t)')
plt.xlabel('t(s)')
plt.ylabel('Norm. hormone conc.')
plt.legend(loc='lower right')
plt.grid(True)

##############################################################################

xsquare = x
x1_new = x1
ysquare = [-4*u*u*pow(x1_new,-2)+(4*u/x1_new) for u in xsquare]

ax = fig.add_subplot(334)
plt.subplot(334)
title = "Numerical solution"
colorbar_legend = "Rising time: BLACK -> RED"

#plt.title(title)
K = len(U_all_x)
cpt = 0
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
    if cpt==cpt_max:
        if ymax==0:
            if k==0:
                print("ymax = 0 and k = 0")
                #plt.vlines(x=x_middle, ymin=0, ymax=1, color=[0,0,0], linestyle='-')
            else:
                plt.vlines(x=x_middle, ymin=0, ymax=1, color=ck, linestyle='-')
        else:
            for i in range(len(yk)):
                yk[i]/=ymax
            plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
            LSM.append(LeastSquaresMethod(ysquare,yk))
            LSM_time.append(t[k])
    if cpt==cpt_max:
        cpt=-1
    cpt+=1
#plt.vlines(x=x_middle, ymin=0, ymax=1, color=Lcolors[0], linestyle='-')

#plt.plot(xsquare,np.array(ysquare)/max(ysquare),color = 'r', marker = '.')

plt.xlabel('Worm axis'+space_unit)
plt.ylabel('Norm. hormone conc.')
plt.legend(loc='lower right')
plt.grid(True)
c = np.array([t0,t1])
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
color_map_to_use = (mpl.colors.ListedColormap(Lcolors).with_extremes(over=Lcolors[0], under=Lcolors[len(Lcolors)-1]))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=color_map_to_use)
cmap.set_array([])
cbar = plt.colorbar(cmap)
cbar.set_label(colorbar_legend, rotation=90)

################################################################################

ax = fig.add_subplot(335)
plt.subplot(335)
title = "Standard deviations in time"
#subtitle = 'interpolation='+str(boolean_interpolation)
if(boolean_interpolation==0):
    label_str = 'numerical solution - no interpolation'
else:
    label_str = 'interpolated numerical solution'
#plt.title(title)
#plt.suptitle(subtitle)
plt.plot(t,np.array(Sigma_numeric), color = color_fig5, marker = '', label = label_str)
plt.xlabel('t(s)')
plt.ylabel('Sigma'+space_unit)
plt.legend(loc='lower right')
plt.grid(True)

################################################################################

ax = fig.add_subplot(336)
plt.subplot(336)
title = "Flux in space at stationary state"
#plt.title(title)
plt.plot (x[1:min(len(x),len(Flux))-1],Flux[1:min(len(x),len(Flux))-1], color = color_fig6, label = 'Flux')
plt.xlabel('Worm axis'+space_unit)
plt.ylabel('Flux(m/s)')
plt.legend(loc='lower right')
plt.grid(True)

################################################################################

ax = fig.add_subplot(337)
plt.subplot(337)

title = "Least Squares Method"
#plt.title(title)
plt.plot (LSM_time,LSM, color = color_fig7, label = 'LSM')
plt.xlabel('t(s)')
plt.ylabel('Error[numeric,x²](t)')
plt.legend(loc='upper right')
plt.grid(True)

################################################################################

plt.show()
#plt.savefig('fig_to_save.png')
