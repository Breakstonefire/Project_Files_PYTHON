from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from math import *
import numpy as np
import matplotlib.pyplot as plt

def Rcall(r):
    s = r**2
    return s

def Tcall(Rcall):
    s = 1-Rcall
    return s

def coeff(Rcall):
    s = 4*Rcall/(Tcall(Rcall)**2)
    return s

def sin2(x):
    s = sin(x)*sin(x)
    return s

#def thetaprime(theta,n):
#    s = arcsin(sin(theta)/n)
#    return s

def DeltaPhi(lbd1,lbd):
    s = 2*pi*lbd1/lbd
    return s

def T3d(lbdmin,lbdmax,lbdstep,lbd1,Rmin,Rmax,Rstep):
    fig = figure(1)
    ax = Axes3D(fig)

    lbdmin/=10**9
    lbdmax/=10**9
    lbdstep/=10**9
    lbd1/=10**9

    Lambda = np.arange(lbdmin , lbdmax , lbdstep)
    Rcalls = np.arange(Rmin , Rmax , Rstep)
    Lambda, Rcalls = np.meshgrid(Lambda, Rcalls)
    T = 1/(1+coeff(Rcalls)*np.sin(DeltaPhi(lbd1,Lambda)/2)*np.sin(DeltaPhi(lbd1,Lambda)/2))
    
    ax.plot_surface(Lambda, Rcalls, T, rstride=1, cstride=1, cmap='hot')
    show()
    return 0

def T2d(lbdmin,lbdmax,lbdstep,lbd1,ListeRcalls):
    fig1 = figure(1)
    lbdmin/=10**9
    lbdmax/=10**9
    lbdstep/=10**9
    lbd1/=10**9

    Lambda = []
    ListeDesListesDeT = []

    for j in range(int((lbdmax-lbdmin)/lbdstep)):
        Lambda.append(lbdmin+j*lbdstep)

    for n in range(len(ListeRcalls)):
        ListeT1 = []
        for k in range(len(Lambda)):
            valT = 1/(1+coeff(ListeRcalls[n])*float(sin2(DeltaPhi(lbd1,Lambda[k])/2)))
            ListeT1.append(valT)
        ListeDesListesDeT.append(ListeT1)

    for n in range(len(ListeRcalls)):
        plt.plot(Lambda,ListeDesListesDeT[n],".--",linewidth = 1,label=ListeRcalls[n]) #scatter affiche les points sans les relier, 'plot' les relie
    plt.grid()
    title("T(R,Lambda)")
    xlabel("Lambda(nm)")
    ylabel("R")
    plt.legend()
    plt.show()
    return 0

T2d(100,800,1,450,[0,0.1,0.5,0.8,0.95,0.99])
T3d(100,800,2,450,0,0.99,0.02)
