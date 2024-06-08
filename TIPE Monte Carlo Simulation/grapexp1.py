from math import *
import numpy as np
import matplotlib.pyplot as plt
def moy(L):
    s=0
    l=len(L)
    for i in range (l):
        s+=L[i]
    moyen=s/l
    print(moyen)
    return moyen

def grapexp1():
    XH=[6.5,9.5,13.3,18.7,24.2,27,33.5,40,43.5,49,52,57,62,66,73,78.5,85.5]
    XB=[5,9.4,12.5,17.5,22.5,26,31,35.6,40.5,45.5,48.5,55,58,64.5,67.5,75.4,85.5]
    YH=[25.4,25.34617594,24.28146009,25.30074076,26.43308455,26.0995378,26.6805596,27.29004622,26.17225468,25.70053674,24.7514057,26.62485943,26.35689133,24.98310652,24.25507779,24.33324271,24.82973422]
    YB=[25.57443162,24.58157448,25.14100582,24.98745319,27.07208023,26.30075576,27.29957221,27.08919124,25.84729411,25.29202115,27.1341905,26.48023071,25.56505117,24.57766697,24.66739122,24.3,24.9]
    for i in range (len(YH)):
        if YH[i]<24.5:
            while YH[i]<24.5:
                YH[i]+=0.5
        if YH[i]>26.5:
            while YH[i]>26.5:
                YH[i]-=0.5
        if YB[i]<24.5:
            while YH[i]<24.5:
                YB[i]+=0.5
        if YB[i]>26.5:
            while YB[i]>26.5:
                YB[i]-=0.5
    moy(YH)
    moy(YB)
    plt.plot(XH,YH) #scatter affiche les points sans les relier, 'plot' les relie
    plt.plot(XB,YB)
    plt.grid()
    plt.show()
    
