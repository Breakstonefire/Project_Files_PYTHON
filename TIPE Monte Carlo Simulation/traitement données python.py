from math import *
#import matplotlib.pyplot as plt
#import numpy as np
import statistics

def répartition(c):
    c.replace(',0',',')
    virg=','
    lc=len(c)
    L=[]
    a=""
    for i in range(lc):
        if c[i]!=virg:
            a+=c[i]
        else:
            b=int(a)
            if b<23500:
                L.append(b)
            a=""
    l1=len(L)
    moy=(statistics.mean(L))
    sigma=(statistics.pstdev(L))
    var=sigma**2
    maxi=max(L)
    b=(2*pi*var)**-0.5
    xgauss=[]
    ygauss=[]
    for i in range(min(L)-100,maxi+100,1):
        a=((i-moy)**2)/(2*var)
        y=b*exp(-a)*3*l1 #4*l car l'aire occupée par l'histogramme vaut l fois le pas ici choisi égal à 4
        xgauss.append(i)
        ygauss.append(y)
    plt.plot(xgauss,ygauss)
    plt.hist(L,bins = [min(L)+pas for pas in range(0,max(L)-min(L)+1,4)])
    plt.xlabel('Temps de vol')
    plt.ylabel("Nombre d'apparition des temps de vol")
    titre="Histogramme des mesures et modèle gaussien, E= ",int(moy),"μs, U= ",int(sigma)
    plt.title(titre)
    plt.show()
    return moy,sigma,l1

