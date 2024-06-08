#calcul de hg
from math import *
import numpy as np
import matplotlib.pyplot as plt

def calculh(Lg,Ld,dg,dd):
    l=len(Lg)
    Lrg=[]
    Lrd=[]
    moyhg=0
    moyhd=0
    for i in range(len(Lg)):
        if Lg[i][0]<=400:
            a=Lg[i][0]
            Lrg.append(a)
            moyhg+=a
    for i in range(len(Ld)):
        if Ld[i][0]<=400:
            b=Ld[i][0]
            Lrd.append(b)
            moyhd+=b
    lg=len(Lrg)
    ld=len(Lrd)
    moyhg/=lg
    moyhd/=ld
    hg=(abs(moyhg**2-dg**2))**(0.5)
    hd=(abs(moyhd**2-dd**2))**(0.5)
    return hg,'  ',hd
