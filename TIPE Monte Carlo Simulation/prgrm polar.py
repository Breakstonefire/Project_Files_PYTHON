#programme de traçage de courbes en polaire
from math import *
import numpy as np
import matplotlib.pyplot as plt

def polar(Lg,Ld):
    l=len(Lg)
    b="oui"
    a=b
    Lrg=[]
    Lrd=[]
    Lmug=[]
    Lmud=[]
    h=19.55
    LAd=[]
    LAg=[]
    Ltetag=[]
    Ltetad=[]
    
    
    for i in range(len(Lg)):
        if Lg[i][0]<=400:
            Lrg.append(Lg[i][0])
            LAg.append(180*asin(h/Lg[i][0])/pi)
            Lmug.append(Lg[i][1])
            Ltetag.append(Lg[i][1]-180*asin(h/Lg[i][0])/pi)
    for i in range(len(Ld)):
        if Ld[i][0]<=400:
            Lrd.append(Ld[i][0])
            LAd.append(180*asin(h/Ld[i][0])/pi)
            Lmud.append(Ld[i][1])
            Ltetad.append(Ld[i][1]-180*asin(h/Ld[i][0])/pi)
                
    if a==b:
        for i in range (len(Ltetag)):
            Ltetag[i]=+Ltetag[i]*pi/180
        for i in range (len(Ltetad)):
            Ltetad[i]=-Ltetad[i]*pi/180
        ax=plt.subplot(projection='polar')
        ax.scatter(Ltetag,Lrg)
        ax.scatter(Ltetad,Lrd)
        #plt.figure(1)
        #plt.polar(Ltetag,Lrg, "b" ) # Courbe polaire en bleu
        #plt.polar(Ltetad,Lrd, "b" ) # Courbe polaire en bleu
        plt.show()
    return 
