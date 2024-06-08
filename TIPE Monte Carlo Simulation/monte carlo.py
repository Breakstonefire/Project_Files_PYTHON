from math import *
import numpy as np
import random
import statistics
import matplotlib.pyplot as plt
from random import gauss


##    x est un double
##    n est l'entier indiquant le nombre de chiffres post-virgule à retenir pour x écrit en écriture scientifique
def arrondir(x,n):
    if x==0:
        return 0
    i=0
    a=x
    while a>=10:
        i+=1
        a=x/(10**i)
    x/=10**i #abaisse y jusqu'à un float dans [1,10[
    j=0
    a=x
    while abs(a)<1:
        j+=1
        a=x*(10**j)
    x*=10**j #augmente y jusqu'à un float dans [1,10[
    for k in range(0,n,1):
        x*=10 #augmente y jusqu'à un float>=1 dont la partie entière et le 1er chiffre post-virgule décrivent l'arrondi finale
    x=(int)(x+0.5) #arrondi de y à l'entier le plus proche
    x*=10**(i-j-n)
    return x


##    tau20= temps de vol aller-retour entre voiture et obstacle
##    taubar0= 2*v*T/c où v=vitesse voiture;T=période de mesure;c=célérité du son dans l'air (c'est le "tau" dans la théorie)
##    utau1= incertitude sur tau1 qui suit environ la loi affine: utau1=0.003*tau1(microsec)+20(microsec) et environ égale à 60microsec pour un obstacle à 3m 
##    utau2= incertitude sur tau2 (la même que celle de tau1)
##    utaubar = incertitude sur taubar
##    teta2 = angle formé par l'axe reliant la voiture à l'obstacle et l'axe de direction de la voiture
##    N = nombre d'appels de la fonction montecarlo à teta2 fixé
##    n = nombre de chiffres post-virgule à afficher après le calcul des différentes variables
def montecarlo(tau20,utau2,utau1,taubar0,utaubar,teta2,N,n):
    teta2rad=teta2*pi/180
    tau10=(tau20**2+taubar0**2+2*tau20*taubar0*cos(teta2rad))**0.5
    T1,T2,T,Lteta,Ltetahist=([],[],[],[],[])
    random.seed(None,2) #réinitialisation de la fonction gaussienne de python au cas où montecarlo est appelé successivement
    for i in range (N):
        tau2=random.gauss(tau20,utau2)
        T2.append(tau2) #ajout d'une valeur proche de tau2 dans T2
    random.seed(None,2)
    for i in range (N):
        tau1=random.gauss(tau10,utau1)
        T1.append(tau1) #même ajout pour T1
    random.seed(None,2)
    for i in range (N):
        taubar=random.gauss(taubar0,utaubar)
        T.append(taubar) #même ajout pour T
        numerateur=T1[i]**2-T2[i]**2-T[i]**2
        denominateur=2*T2[i]*T[i]
        costeta=numerateur/denominateur
        Lteta.append(costeta)
        Ltetahist.append(int(1000*costeta))
    moy=arrondir(statistics.mean(Lteta),n) #moyenne arrondie à n chiffres post-virgule des teta obtenus
    sigma=arrondir(statistics.stdev(Lteta,moy),n) #écart type sur teta arrondi à n chiffres post-virgule
    if moy>1:
        reporteràlamoy=0
    if moy-sigma>1:
        Ecartteta2=0
    if moy<-1:
        reporteràlamoy=180
    if moy-sigma<-1:
        Ecartteta2=180
    if  -1<=moy<=1:
        reporteràlamoy=acos(moy)*180/pi
    if -1<=moy-sigma<=1:
        Ecartteta2=acos(moy-sigma)*180/pi
    Ecartteta2-=reporteràlamoy
    Ecartteta2=arrondir(Ecartteta2,n) #calcul de l'écart à la moyenne de teta2 arrondi à n chiffres post-virgule
##    titre="Histogramme des mesures, E= ",moy,"μs, U= ",sigma," Estimation de l'incertitude sur teta= ",Ecartteta2
##    plt.hist(Ltetahist,bins = [min(Ltetahist)+pas for pas in range(0,max(Ltetahist)-min(Ltetahist)+1,1)]) #tracé de l'histogramme
##    plt.xlabel('cos(teta) obtenus')
##    plt.ylabel("Nombre d'apparition des cos(teta)")
##    plt.title(titre)
##    plt.show()
    return moy,sigma,Ecartteta2

def UtetaFctionTeta(tau20,utau2,utau1,taubar0,utaubar,teta2min,teta2max,pas,N,n):
    Luteta=[]
    Lteta=[]
    teta=teta2min
    while teta2min<=teta<=teta2max:
        uteta=montecarlo(tau20,utau2,utau1,taubar0,utaubar,teta,N,n)[2]
        Luteta.append(uteta)
        Lteta.append(teta)
        teta+=pas
    plt.plot(Lteta,Luteta) #tracé du graphe
    plt.xlabel('teta')
    plt.ylabel('Uteta')
    titre="Graphe de U(teta)=f(teta)"
    plt.title(titre)
    plt.show()
    return Luteta
