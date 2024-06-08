#tracé d'histogramme
from math import *
import hist
import statistics
import numpy as np
import matplotlib.pyplot as plt

def Em(Lind,Et):
    em = 0
    for i in range(len(Lind)):
        ind = Lind[i]
        em+=Et[ind]
    em/=1+len(L)
    return em

def F_Elg(Titles,i):
    name_fic = Titles[i]
    f = np. genfromtxt(name_fic,delimiter=';',usecols=np.arange(0,2))
    print(f)
    F = []
    Elg = []
    for i in range(1,len(f)):
        force = f[i][1]
        elong = f[i][0]
        F.append(force)
        Elg.append(elong)
    return F,Elg


Et = [220 , 2370 , 2360 , 2380 , 2200 , 2360 , 1850 , 2110 , 2120]
eps = [0.9 , 3.7 , 3.6 , 3.7 , 3.7 , 3.8 , 3.9 , 3.9 , 3.9 ]
Titles = ["PLA_10_rect1.TRA","PLA_10_rect2_pas_casse.TRA","PLA_10_rect3.TRA","PLA_50_grid_1_bizzar.TRA","PLA_50_grid_2.TRA","PLA_50_grid_3.TRA","PLA_100_rect_1.TRA","PLA_100_rect2.TRA","PLA_100_rect3.TRA"]
Lcolors = ["b","y","r"]
#f = np. genfromtxt('PLA_10_rect1.TRA',delimiter=';',usecols=np.arange(0,2))
F = []
Elg = []


L = [6,7,8]


for k in range(len(L)):
    col = Lcolors[k]
    F,Elg = F_Elg(Titles,L[k])
    plt.plot(Elg,F,col)
plt.xlabel("Relative Elongation(%)")
plt.ylabel("Force(MPa)")
Emoyen = round(Em(L,Et))
eps_moyen = round(Em(L,eps),1)
titre="Mechanical Uniaxial Tensile Test on PLA - 100%"
subtitle = 'Rectilinear Pattern ; E='+str(Emoyen)+' ; eps='+str(eps_moyen)+'%'
plt.suptitle(titre)
plt.title(subtitle)
plt.show()


