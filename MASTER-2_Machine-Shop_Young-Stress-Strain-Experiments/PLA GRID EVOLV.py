#tracé d'histogramme
from math import *
import hist
import statistics
import numpy as np
import matplotlib.pyplot as plt
import os

#Astuces pour obtenir rapidement les titres de tous les fichiers d'un dossier:
#entrer import os (dans le shell)
#entrer os.listdir()
#copier le résultat qui contient tous les noms des fichiers (et dossiers)

Titles = os.listdir()
Lcolors = ["b-","y-","r-"]
#f = np. genfromtxt('PLA_10_rect1.TRA',delimiter=';',usecols=np.arange(0,2))
#usecols permet de sélectionner les colonnes de données, en mettant (0,2) on prend la colonne d'indice 0 à 1

pattern = "Star"
Type_of_Polymer = "HIPS"
filling_rate = 50

def f(L):
    Li = list(L)
    maxi = max(L)
    for i in range(len(L)):
        Li[i]/=maxi
    return Li


Rate_PLA_rect = [10,50,100]
E_PLA_rect = [1577/1577,1520/1577,1519/1577] 
eps_PLA_rect = [2.4,2.5,2.9]
eps_PLA_rect = f(eps_PLA_rect)

Rate_PLA_square = [25,50]
E_PLA_square = [602.1/604,604/604]
eps_PLA_square = [2,1.9]
eps_PLA_square = f(eps_PLA_square)

Rate_PLA_star = [25,50]
E_PLA_star = [623.2/623.2,622.1/623.2]
eps_PLA_star = [2.2,2.6]
eps_PLA_star = f(eps_PLA_star)

Rate_HIPS_cube = [25,50]
E_HIPS_cube = [323/323.4,323.4/323.4]
eps_HIPS_cube = [10.3/10.3,10.1/10.3]
eps_HIPS_cube = f(eps_HIPS_cube)

Rate_HIPS_star = [25,50]
E_HIPS_star = [305.3/322.8,322.8/322.8]
eps_HIPS_star = [9.7,9.2]
eps_HIPS_star = f(eps_HIPS_star)

##plt.plot(Rate_PLA_rect,E_PLA_rect,'b',      label = "PLA - Rectilinear")
##plt.plot(Rate_PLA_square,E_PLA_square,'y',  label = "PLA - Square")
##plt.plot(Rate_PLA_star,E_PLA_star,'r',      label = "PLA - Star")
##plt.plot(Rate_HIPS_cube,E_HIPS_cube,'g',    label = "HIPS - Cubic")
##plt.plot(Rate_HIPS_star,E_HIPS_star,'k',    label = "HIPS - Star")
##plt.xlabel("Filling rate(%)")
##plt.ylabel("Young Modulus Ratio(-)")

plt.plot(Rate_PLA_rect,eps_PLA_rect,'b',      label = "PLA - Rectilinear")
plt.plot(Rate_PLA_square,eps_PLA_square,'y',  label = "PLA - Square")
plt.plot(Rate_PLA_star,eps_PLA_star,'r',      label = "PLA - Star")
plt.plot(Rate_HIPS_cube,eps_HIPS_cube,'g',    label = "HIPS - Cubic")
plt.plot(Rate_HIPS_star,eps_HIPS_star,'k',    label = "HIPS - Star")
plt.xlabel("Filling rate(%)")
plt.ylabel("Relative elongation ratio(-)")

titre = 'Mechanical Uniaxial Tensile Test - varying filling rates'
subtitle = 'Study on relative elongation'
plt.suptitle(titre)
plt.title(subtitle)
plt.legend()
plt.grid(True)
plt.show()
