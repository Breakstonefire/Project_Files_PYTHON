from math import *
import numpy as np
import matplotlib.pyplot as plt

def incert(lim,pas,coef):
    X=[]
    Y=[]
    i=35
    while(35>=i>=lim):
        X.append(i)
        Y.append(coef/sin(i*pi/180))
        i+=pas
    titre="u(teta)=f(teta)"
    plt.plot(X,Y)
    plt.legend()
    plt.title(titre)
    plt.xlabel('Teta(°)')
    plt.ylabel('u(teta) (°)')
    plt.grid()
    plt.show()
