#graphe 3D de la probabilité de passage de l'état de basse énergie (ms=-1/2) à l'état de haute énergie(ms=+1/2) en fonction
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from math import *
import numpy as np
import matplotlib.pyplot as plt

def coeffs(w0,w1,w):
    wL = 2*w0
    dw = wL - w
    
    r = dw/(2*w1)
    s = 1/r
    r2 = r**2
    s2 = s**2
    term = np.sqrt(1+(s**2))
    c0 = (-w1 + 1/np.sqrt((1/(w1**2))+(2/dw)**2))**2
    c1 = 1/(-1+(1/term))
    c2 = w1*np.sqrt(1+(dw/2)**2)
    c3 = -c1
    c4 = c2
    #les ci apparaissent de la manière suivante dans la formule
    #de la probabilité: Proba = c0*(c1*cos(c2*t)+c3*sin(c4*t))**2
    Lc=[]
    Lc.append(c0)
    Lc.append(c1)
    Lc.append(c2)
    Lc.append(c3)
    Lc.append(c4)
    return Lc

def cosinus(X):
    S = np.cos(X)
    return S

def sinus(X):
    S = np.sin(X)
    return S

def Probapvm(w0,w1,wmin,wmax,wstep,tmin,tmax,tstep):
    fig1 = figure(1)
    ax = Axes3D(fig1)
    Lw = np.arange(wmin,wmax,wstep)
    Lt = np.arange(tmin,tmax,tstep)
    Lw,Lt = np.meshgrid(Lw,Lt)
    C = coeffs(w0,w1,Lw)
    C0,C1,C2,C3,C4 = C[0],C[1],C[2],C[3],C[4]
    Ppvm = C0*(C1*cosinus(C2*Lt)+C3*sinus(C4*Lt))**2
    print(type(Ppvm))
    print(len(Ppvm))
    print(len(Ppvm[0]))
    print(len(Ppvm[1]))
    ax.set_title("Porbabilité état haut vers état bas")
    ax.set_xlabel("Axe des w (GHz*precision choisie en entrée)")
    ax.set_ylabel("Axe des temps (s)")
    ax.set_zlabel("Probabilité de changement d'état")
    ax.plot_surface(Lw, Lt, Ppvm, rstride=1, cstride=1, cmap='Blues')
    #show()
    return None

##constantes
g = 2.002319314
esurm = (1.0602/9.1094)*10**(12)
B0 = 1


##paramètres d'entrée
wmin = 100
wmax = 140
wstep = 0.2
tmin = 0
#nT est le nombre de périodes à afficher (à remplir)
nT = 1
puissance_precision = 0
precision = 10**puissance_precision

    
##grandeurs calculées
w0 = g*esurm*B0*precision/(4*10**9) #rapporté en GHz
w1 = w0/1000
wmin*= precision
#calcul de la période des oscillations à partir de la pulsation du sinus à w = wL 
T = pi/w1
wmax*= precision
tmax = tmin + nT*T
tstep = T/500

Probapvm(w0,w1,wmin,wmax,wstep,tmin,tmax,tstep)
