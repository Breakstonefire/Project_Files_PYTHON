#tracé d'histogramme
from math import *
#import hist
import statistics
import numpy as np
import matplotlib.pyplot as plt
import random
from random import gauss

def histog(c):
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
    moy=int(statistics.mean(L))
    sigma=int(statistics.pstdev(L))
    var=sigma**2
    xgauss=[]
    ygauss=[]
    maxi=0
    l=len(L)
    for i in range(l):
        if L[i]>maxi and L[i]!=23500:
            maxi=L[i]
    b=(var*2*pi)**-0.5
    for i in range(min(L)-1000,maxi+1000,1):
        a=((i-moy)**2)/(2*var)
        y=b*exp(-a)*4*l
        xgauss.append(i)
        ygauss.append(y)
    plt.plot(xgauss,ygauss)
    plt.hist(L,bins = [min(L)+pas-1 for pas in range(0,max(L)-min(L)+1,4)])#, range = (min(L)-1, max(L)+1), bins = 1000000, color = 'red',edgecolor = 'blue')
    plt.xlabel('Temps de vol')
    plt.ylabel("Nombre d'apparition des temps de vol")
    titre="Histogramme des mesures et modèle gaussien, E= ",moy,"μs, U= ",sigma
    plt.title(titre)
    plt.show()
    return moy,sigma


