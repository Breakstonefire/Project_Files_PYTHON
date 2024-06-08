from math import *
import matplotlib.pyplot as plt
import numpy as np

def répartangle(L,d,ecart):
    l=len(L)
    S=0
    n=0
    for i in range(l):
        term=L[i][0]
        if d-10<=term<=d+10:
            S+=term
            n+=1
    Moy=S/n
    Lang=[]
    for i in range(l-1):
        if L[i][0]>Moy+ecart and Moy+ecart>=L[i+1][0]>=Moy-ecart:
            if L[i][1]<90:
                Lang.append(180-L[i][1])
            else:
                Lang.append(L[i+1][1])
        elif L[i+1][0]>Moy+ecart and Moy+ecart>=L[i][0]>=Moy-ecart:
            if L[i+1][1]<90:
                Lang.append(180-L[i+1][1])
            else:
                Lang.append(L[i+1][1])
    l1=len(Lang)
    L2=[]
    X=[]
    E1=0 #déclare l'espérance à élever au carré
    E2=0 #déclare l'espérance du carré
    Lproba=[] #déclare la liste des probabilités
    for k in range (l2):
        P=L2[k][1]/l1
        E1+=L2[k][0]*P
        Lproba.append(P) #ordonnée du graphe
        E2+=(L2[k][0]**2)*P
    V=E2-(E1**2)
    Sigma=(V**0.5)
    #commentaire='(Moyenne;écart type)= '+str(E1)+';'+str(Sigma)
    #titre=commentaire
    
    return l1,E1,V,Sigma
