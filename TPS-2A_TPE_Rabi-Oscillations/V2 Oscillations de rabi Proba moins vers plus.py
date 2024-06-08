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
    r2 = r**2
    denom = 1+r2
    c0 = 1/denom
    c1 = w1*np.sqrt(denom)
    Lc=[]
    Lc.append(c0)
    Lc.append(c1)
    return Lc

def sin2(X):
    S = np.sin(X)
    S*=S
    return S

def control():
    a=0
    a = int(input("Passer à la suite? (0 ou 1)"))
    if (a==0):
        return control()
    else:
        return None

def Probamvp(w0,w1,wmin,wmax,wstep,tmin,tmax,tstep,precision):
    fig1 = figure(1)
    ax = Axes3D(fig1)
    Lw = np.arange(wmin,wmax,wstep)
    Lt = np.arange(tmin,tmax,tstep)
    Lwbis,Ltbis = np.meshgrid(Lw,Lt)
    
    C = coeffs(w0,w1,Lwbis)
    C0 = C[0]
    C1 = C[1]
    Pmvp3D = C0 * sin2(C1*Ltbis)
    
    cwL = coeffs(w0,w1,2*w0)
    c0 = cwL[0]
    c1 = cwL[1]
    Pmvp2D = c0*sin2(c1*Lt)
    Pmoy2D = []
    Pmoy2D.append(Pmvp2D[0])
    for i in range(1,len(Lt)):
        s=0
        s+=Pmoy2D[i-1]
        s+=Pmvp2D[i]
        Pmoy2D.append(s)
    for i in range(0,len(Lt)):
        Pmoy2D[i]/=(i+1)

    print(type(Lwbis),type(Ltbis),type(Pmvp3D))



    ax.set_title("Porbabilité état bas vers état haut")
    precision = str(precision)
    wlabel = "Axe des w ("+precision+" GHz)"
    ax.set_xlabel(wlabel)
    ax.set_ylabel("Axe des temps (s)")
    ax.set_zlabel("Probabilité de changement d'état")
    ax.plot_surface(Lwbis, Ltbis, Pmvp3D, rstride=1, cstride=1, cmap='Blues', label="P(w;t)")
    show()

    plot(Lt,Pmvp2D,color="blue",label="P(w=wL;t)")
    plot(Lt,Pmoy2D,color="red",label="Pmoyen(w=wL;t)")
    plt.title("Probabilité état bas vers état haut (en bleu)")
    legend(loc='upper left')
    show()
    return 


##constantes
g = 2.002319314
esurm = (1.0602/9.1094)*10**(12)
B0 = 1


##paramètres d'entrée
wmin = 110
wmax = 120
wstep = 2
tmin = 0
#nT est le nombre de périodes à afficher (à remplir)
nT = 20
#plus puissance_precision augmente, plus on va pouvoir calculer à des intervalles faibles sur l'axe des w
## En augmentant de un 'puissance_precision' et en multipliant par 10 wstep, on se retrouve avec le même graphe qu'avant mais
## les w seront en 'precision*GHz' au lieu d'être simplement en GHz si on avait choisi precision = à 0
puissance_precision = 0
precision = 10**puissance_precision


##grandeurs calculées
w0 = g*esurm*B0*precision/(4*10**9) #rapporté en GHz
w1 = w0/1000
wmin*= precision
#calcul de la période des oscillations à partir de la pulsation du sinus carré à w = wL
T = pi/w1
wmax*= precision
tmax = tmin + nT*T
tstep = T/1000
Ncalculs = (wmax-wmin)*(tmax-tmin)/(wstep*tstep)
Ncalculs = int(Ncalculs)

print("Nombre de points calculés = ",Ncalculs)
if (Ncalculs >= 500000):
    print("Le graphe n'est plus fluide pour un nombre de points à afficher dépassant 500 000.")
    go_on = str(input("Voulez vous continuer: (oui ou non)"))
    if (go_on == "oui"):
        Probamvp(w0,w1,wmin,wmax,wstep,tmin,tmax,tstep,precision)
else:
    Probamvp(w0,w1,wmin,wmax,wstep,tmin,tmax,tstep,precision)

