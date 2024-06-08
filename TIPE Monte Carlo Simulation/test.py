from math import *
import matplotlib.pyplot as plt
import numpy as np
def affichage():
    L=[[100,5862.074529667149, 31.85693330043641],[120,7060.191111111111, 51.70881215345202],[140,8206.19776536313, 44.115009613720915],[160,9418.848314606741, 52.094129610312],[180,10648.675830469645, 51.60415528904738],[200,11844.612826603325, 50.45658256203339],[220,13743.405479452054, 70.6951741231273],[240,15229.6327014218, 67.966528589117],[260,16116.20216606498, 69.6306388516925],[280,17633.418367346938, 75.8205372495736],[300,17546.32754342432, 62.2912921682643],[320,18087.01765447667, 78.15116067378205],[340,20113.734666666667, 93.8409273908859],[360,21054.43918053777, 89.6243448676343],[380,21999.72785622593, 72.46237840675973]]
    titre="Ecart type en fonction du temps de vol"  
    #plt.scatter(Ldist,Lectypang)
    D=[]
    E=[]
    U=[]
    V=[]
    l=len(L)
    D2=[]
    for i in range(l):
        D.append(L[i][0])
        D2.append(L[i][0]**2)
        E.append(L[i][1])
        U.append(L[i][2])
        V.append(L[i][2]**2)
    #plt.scatter(X,Y1,'b')
    plt.scatter(E,U)
    plt.legend()
    plt.title(titre)
    plt.xlabel('Moyenne des tof')
    plt.ylabel('Ecart type')
    plt.grid()
    plt.show()
