from math import *
import matplotlib.pyplot as plt
import numpy as np

def affichage():
    Lectdelaims=[8.087327499506298,6.082501802,2.2116965140926643,7.946495459172216,8.03691153251735,7.955032260481797,7.9224411813227125,7.990362790303153,7.980945723658466]
    Ldelai=[200,160,140,120,100,80,70,60,50]
    for i in range(len(Ldelai)):
        Lectdelaims[i]**=2
        Ldelai[i]**=2
    plt.plot(Ldelai,Lectdelaims)
    plt.legend()
    #plt.title(titre)
    plt.xlabel('')
    plt.ylabel('')
    plt.grid()
    plt.show()
