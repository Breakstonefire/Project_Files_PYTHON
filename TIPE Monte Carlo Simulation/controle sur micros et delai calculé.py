from math import *
import matplotlib.pyplot as plt
import numpy as np

def affichage(L):
    X=[]
    Y=[]
    for i in range(len(L)-1):
        X.append(i)
        Y.append((L[i+1][0]-L[i][0])/1000) #calcule le délai en ms
    plt.plot(X,Y)
    plt.legend()
    #plt.title(titre)
    plt.xlabel('')
    plt.ylabel('')
    plt.grid()
    plt.show()
