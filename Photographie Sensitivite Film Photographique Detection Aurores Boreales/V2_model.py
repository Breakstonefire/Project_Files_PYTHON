#### ATTENTION INCOMPREHENSION DANS LE ROLE DE L'ISO : EN FAIT
# L'ISO EST PLUTOT LA SENSIBILITE DU FILM PHOTOGRAPHIQUE DONC SELON WIKIPEDIA
# C'EST LE RAPPORT H0/H IL FAUT FAIRE COMME SI ISO = H0*cste_d_offset DANS
# LE PROGRAMME PRECEDENT ET DANS CELUI CI

# Ce programme permet (à partir du précédent, donc de la v1) de tracer la courbe du rapport signal sur bruit en fonction des trois paramètres de l'expérience:
# 1. Le temps d'exposition totale allant de 1/3200 seconde à 15 secondes
# 2. Le diamètre d'ouverture du diaphragme (àdistance focale fixée et constante)
# 3. La sensibilité ISO
# Le graphe final peut en réalité être représenté par trois graphes séparés donnant l'allure moyenne de la dépendance de chacun des paramètres en fixant les deux autres
# Ou bien être un ensemble de graphes 3D qui montrent la dép. du RSB aux deux autres paramètres (pour toutes leurs valeurs) en fixant le dernier paramètre (de préférence celui qui a le moins de valeurs possibles)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing libraries # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from math import *
import random
from pylab import *
import random
from mpl_toolkits.mplot3d import Axes3D

#### EXCELLENT COURS SUR LA SENSITOMETRIE ####:
# http://www.photographeaparis.fr/tech-photo/index.htm

#### Caractéristiques de l'appareil photo ###
# https://www.canon.fr/for_home/product_finder/cameras/digital_camera/powershot/powershot_sx170_is/specification.html

#### Lien vers la règle NPF pour un calcul du temps d'exposition maximal théorique permettant d'éviter un filé de source lumineuse sur l'image lors de l'acquisition :
# http://web.archive.org/web/20200220123345/https://www.sahavre.fr/tutoriels/astrophoto/34-regle-npf-temps-de-pose-pour-eviter-le-file-d-etoiles

#### Fonctions du calcul principal ####

##I_Signal_Pur_1(I_Signal_Pur_0,alpha,L)
##I_lumineuse_1(I_Signal_Pur_1,I_Bruit_Photon_Signal,I_Bruit_Pollution_Lumineuse)
##I_lumineuse_2(I_lumineuse_1,Dmin,Dmax,Hc,Techantillonnage,Dod,f,lambda_do)
##I_electrique_1(I_Saturation,Temperature,Techantillonnage,I_lumineuse_2,Surface_Capteur,U_seuil_photodiode_passante,I_Bruit_Courant_Obscurite)
##I_electrique_3(Techantillonnage,constante_conversion,ISO,I_electrique_1)

def I_Signal_Pur_1(I_Signal_Pur_0,alpha,L):
    return I_Signal_Pur_0*np.exp(-alpha*L)

def I_lumineuse_1(I_Signal_Pur_1,I_Bruit_Photon_Signal,I_Bruit_Pollution_Lumineuse):
    return I_Signal_Pur_1 + I_Bruit_Photon_Signal + I_Bruit_Pollution_Lumineuse

def I_lumineuse_2(I_lumineuse_1,Dmin,Dmax,Hc,Techantillonnage,Dod,f,lambda_do):
    # Cette fonction est une fonction logistique utilisée pour la sensitivité d'un film photographique lorsqu'on trace sa caractéristique (absorbance (ou densité optique) en fonction du log en base 10 de l'exposition (ou la lumination)). La page wikipedia "Fonction logistique (Verhulst)" l'explique parfaitement.
    Lumination_ratio = Hc*f/(I_lumineuse_1*Techantillonnage*Dod)
    exposant = Dmin+(Dmax-Dmin)/(1+pow(Lumination_ratio,lambda_do/log(10)))
    return I_lumineuse_1/pow(10,exposant)

def I_electrique_1(I_Saturation,Temperature,Techantillonnage,I_lumineuse_2,Surface_Capteur,U_seuil_photodiode_passante,I_Bruit_Courant_Obscurite):
    kB = 1.380649/pow(10,23)
    e = 1.602/pow(10,19)
    return I_Saturation*((1/(kB*Temperature))*(Techantillonnage*I_lumineuse_2*Surface_Capteur - e*U_seuil_photodiode_passante)-1) + I_Bruit_Courant_Obscurite

def I_electrique_3(Techantillonnage,constante_conversion,ISO,I_electrique_1):
    return Techantillonnage*constante_conversion*ISO*I_electrique_1

###################################
#### Fonctions complémentaires ####

def focal(L,grandissement):
    return L/(1-(1/grandissement))

def factoriel(k):
    if (k<=1):
        return 1
    else:
        return factoriel(k-1)

def fonction_masse_poisson(Lambda,N_points,Nacquisition):
    prod_k = 1
    Lx_fct_masse = []
    L_fct_masse = []
    for k in range(N_points):
        s = Nacquisition*np.exp(-Lambda)*pow(Lambda,k)/prod_k
        if (k>1):
            prod_k*=k
        L_fct_masse.append(s)
        Lx_fct_masse.append(Lambda*k)
    return Lx_fct_masse,L_fct_masse

# fonctions de modélisation d'une loi de poisson récupérées en ligne sur ce site : https://python.jpvweb.com/python/mesrecettespython/doku.php?id=loi_poisson
def poisson(k,m):
    """poisson(k,m): donne la probabilité d'avoir k évènements distribués selon une loi de Poisson de paramètre m"""
    p=e**(-m)
    for i in range(0,k):
        p*=m/k
        k-=1
    return p

def hpoisson(m,nb=0):
    """hpoisson(m,nb=0): Génération de valeurs tirées au hasard selon une distribution de Poisson"""
    def _hpoisson(m):
        ph = random.random()
        k = 0
        pc = poisson(k,m)
        while ph>pc:
            k+=1
            pc+=poisson(k,m)
        return k
 
    if nb==0:
        return _hpoisson(m)
    else:
        R=[]
        for j in range(0,nb):
            R.append(_hpoisson(m))
        return R

#print(hpoisson(100,13))

def Ecriture_scientifique(x):
    if (x<0):
        x=abs(x)
    if (x==0):
        return 0,0
    if (x<1):
        cpt = -1
        x*=10
        while(x<1):
            x*=10
            cpt-=1
        return x,cpt
    if (x==1):
        return x,0    
    if (x>1):
        cpt = 1
        x/=10
        while(x>1):
            x/=10
            cpt+=1
        return x,cpt

def Bins_and_range(L,Y_moyen):
    range_ = min(L),max(L)
    bins = int((max(L)-min(L))/Y_moyen)
    return bins,range_

#### Listes et constantes ####

focal_min = 5/1000
focal_max = 80/1000
L_ISO = [100]#,200,400,800,1600]
L_Texposition = [1/3200,1/2500,1/2000,1/1600,1/1250,1/1000,1/800,1/640,1/500,1/400,1/320,1/250,1/200,1/160,1/125,1/100,1/80,1/60,1/50,1/40,1/30,1/25,1/20,1/15,1/13,1/10,1/8,1/6,1/5,1/4,0.3,0.4,0.5,0.6,0.8,1,1.3,1.6,2,2.5,3.2,4,5,6,8,10,13,15]
L_Diametre_ouverture_diaphragme = [1/3.5,1/4,1/4.5,1/5,1/5.6,1/6.3,1/7.1,1/8]

##L_ISO = [100,200,400,800,1600,3200]
##L_Texposition = [1/3200,1/2500,1/2000,1/1600,1/1500]
##L_Diametre_ouverture_diaphragme = [1/4.5,1/4,1/3.5,]

####################
#### PARAMETRES ####
# LES PARAMETRES SONT DONNES DANS L'ORDRE CHRONOLOGIQUE LORSQU'ON POSE LE SYSTEME GLOBAL SUR PAPIER ET QU'ON L'ETUDIE DEPUIS LE DEBUT.
#grandissement = 1 # le grandissement est la valeur du zoom donc >=1
I_Signal_Pur_0 = 0.20                                   # ESTIMABLE C'est l'intensité du signal pur sous une aurore boréale environ aussi lumineuse que la moitié de la lune.
lambda_lumiere = 550/pow(10,9)                          # ESTIMABLE On prend la longueur d'onde du vert (~550nm) car ce sont les aurores courantes entre 100 et 300km d'altitude. La longueur d'onde varie car les couleurs des aurores dépendent de l'altitude et de la longueur d'onde.
h_cst_Planck = 6.62607015/pow(10,34)                    # CSTE PHY
c_celerite_lumiere = 299792458                          # CSTE PHY vitesse de la lumiere dans le vide (m/s)
n_indice_refraction_air = 1.000293                      # ESTIMABLE adimensionné
v_lumiere = c_celerite_lumiere/n_indice_refraction_air  # ESTIMABLE
alpha = 0.00001                                         # ESTIMABLE coefficient d'absorption (m^-1) du milieu dans lequel se propage l'onde électromagnétique depuis l'objet jusqu'à l'objectif. Lien vers une page web utile pour les valeurs de l'air dans toute l'atmosphère : http://ressources.univ-lemans.fr/AccesLibre/UM/Pedago/physique/02/optigeo/loidebeer.html
L = 100*pow(10,3)                                       # ESTIMABLE Distance entre l'aurore boréale et le capteur photographique. Les aurores sont à une altitude variable allant de 100 à 600km d'altitude (se référer à la couleur pour connaitre l'altitude moyenne)
f = focal_min                                           # PARAMETRE focal(L,grandissement) # la valeur de la focal dépend du zoom
Texposition = 1                                         # PARAMETRE temps (sec) pendant lequel le(/les) capteurs restent ouverts et aptent de la lumière et sont donc excités et envoient un signal électrique au reste du circuit électronique
Techantillonnage = min(L_Texposition)                   # ESTIMABLE c'est le temps le plus court durant lequel l'appareil peut prendre une image (a priori il est <= au temps d'obturation minimal = 1/3200 sec, on prendra donc cette valeur en référence)
Nacquisition = int(Texposition/Techantillonnage)        # ESTIMABLE c'est le nombre de fois où l'appareil va acquérir une image (peu lumineuse ou non) pendant tout le temps d'exposition "Texposition" en prenant une image tous les "Techantillonnage"

print(Nacquisition) #ligne de contrôle de la valeur du nombre de points qui seront calculés lors d'une acquisition

#### bruits et photocourants parasites :
# Quelques liens vers les bruits présents dans une photodiode : https://fr.lambdageeks.com/photo-detector-photo-diode/
# https://media4.obspm.fr/public/ressources_lu/pages_analyser/bruit-photons_impression.html

#### Point sur les intensités lumineuses de référence de nuit :
#   Éclairement lumineux    Exemple
#   0,25 lux	            Pleine lune par une nuit claire
#   0,01 lux	            Quartier de lune
#   0,002 lux	            Ciel étoilé sans lune

Dmin = 0                                                # ESTIMABLE On suppose que lorsque l'intensité lumineuse incidente tend vers 0 (donc la lumination de même), on a une absorbance (ou donc une densité optique) du capteur qui tend vers 0 puisqu'aucun photon n'est absorbé (hormis ceux du bruit mais le Dmin est exprimé dans un cas idéal d'exposition lumineuse).
Dmax = 1                                                # INCONNUE C'est la valeur maximale d'absorbance que le capteur adopte lorsqu'il est soumis à des intensités lumineus très élevées.
Hc = 100000                                             # INCONNUE
Dod = 0.01                                              # INCONNUE 
lambda_do = 1                                           # INCONNUE
I_Saturation = 10/pow(10,9)                             # Pour la valeur d'application numérique on se base sur les transistors  à jonctions polarisées (PN et/ou NP) où on trouve des Isat de l'ordre de 10 nA), on peut se référer à ce lien : http://res-nlp.univ-lemans.fr/NLP_C_M15_G01/co/Contenu_c4.html
Temperature = 273.15 + 2                                # Température en degrés Kelvin ici on prend 2°C pour une nuit d'hiver classique en Novembre en Norvège
Surface_Capteur = 4.363117258846936e-12                 # Surface (en m²) élémentaire d'un capteur unitaire du capteur CCD de l'appareil photo.
U_seuil_photodiode_passante = 1
constante_conversion = 1                                # Cette constante est une constante qui permet de conserver l'unité de l'intensité électrique lorsqu'on multiplie l'intensité électrique en sortie de capteur par le gain après être passé dans l'amplificateur analogique. étant une simple constante permettant ce respect des unités, on la prend égale à 1.
ISO = L_ISO[0]                                          # Les valeurs de l'ISO sont dans la liste "L_ISO" et caractérisent la sensibilité du film photographique composé par tout les capteurs du CCD.

I_BPS_moyen = I_Signal_Pur_0                            # Intensité moyenne (W/m²) du bruit de photon du signal. On prend ici le signal pur lui même (ou la racine carrée du signal pur) pour avoir une approximation correcte du bruit
I_BPL_moyen = 0.01                                      # Intensité moyenne (W/m²) du bruit de pollution lumineuse. On prend l'équivalent d'un éclairement de nuit sous un quartier de lune.
I_BCO_moyen = I_Saturation                              # Intensité électrqiue moyenne du bruit dû au courant d'obscurité, qui est un courant résiduel et parasite de notre signal électrique 

n_sigma = 0 #coefficient mulitplicateur devant les différents écarts types des bruits lorsque l'on affiche les histogrammes, ce nombre indique jusqu'à combien d'écarts types on affiche l'histogramme.
n_sigma+=1

#### FIN DES PARAMETRES ####
############################

# Calcul des Intensités des différents bruits selon leurs processus de poisson respectifs
I_BPS_moyen_ecriture_scientifique,power_I_BPS_moyen = Ecriture_scientifique(I_BPS_moyen)
I_BPL_moyen_ecriture_scientifique,power_I_BPL_moyen = Ecriture_scientifique(I_BPL_moyen)
I_BCO_moyen_ecriture_scientifique,power_I_BCO_moyen = Ecriture_scientifique(I_BCO_moyen)

L_I_Bruit_Photon_Signal = hpoisson(I_BPS_moyen_ecriture_scientifique,Nacquisition)
for i in range(len(L_I_Bruit_Photon_Signal)):
    L_I_Bruit_Photon_Signal[i]*=pow(10,power_I_BPS_moyen)

L_I_Bruit_Pollution_Lumineuse = hpoisson(I_BPL_moyen_ecriture_scientifique,Nacquisition)
for i in range(len(L_I_Bruit_Pollution_Lumineuse)):
    L_I_Bruit_Pollution_Lumineuse[i]*=pow(10,power_I_BPL_moyen)

L_I_Bruit_Courant_Obscurite = hpoisson(I_BCO_moyen_ecriture_scientifique,Nacquisition)
for i in range(len(L_I_Bruit_Courant_Obscurite)):
    L_I_Bruit_Courant_Obscurite[i]*=pow(10,power_I_BCO_moyen)

L_IBPS_cumul = []
L_IBPL_cumul = []
L_IBCO_cumul = []
L_I_elec_3 = []
L_I_elec_3_cumul = []
L_I_elec_3_non_bruitee = []
L_I_elec_3_non_bruitee_cumul = []
Lt = []

Ielec_bruitee_cumul = 0
Ielec_non_bruitee_cumul = 0
IBPS_cumul = 0
IBPL_cumul = 0
IBCO_cumul = 0

L_RSB = [[[] for u in L_Diametre_ouverture_diaphragme] for v in L_ISO]

for i_ISO in range(0,len(L_ISO)):
    ISO = L_ISO[i_ISO]
    print(" ")
    print("ISO = ", ISO)
    for i_Diametre_ouverture_diaphragme in range(0,len(L_Diametre_ouverture_diaphragme)):
        Dod = L_Diametre_ouverture_diaphragme[i_Diametre_ouverture_diaphragme]
        print(" ")
        print("Diamètre d'Ouverture du Diapragme = ", Dod)
        for i_Texposition in range(0,len(L_Texposition)):
            Texposition = L_Texposition[i_Texposition]

            # Début de la section permettant de calculer le RSB selon les différentes valeurs des paramètres, en fin de boucle. 
            RSB = 0
            Nacquisition = int(Texposition/Techantillonnage)

            # Section de génération des bruits
            L_I_Bruit_Photon_Signal = hpoisson(I_BPS_moyen_ecriture_scientifique,Nacquisition)
            for i in range(len(L_I_Bruit_Photon_Signal)):
                L_I_Bruit_Photon_Signal[i]*=pow(10,power_I_BPS_moyen)
            L_I_Bruit_Pollution_Lumineuse = hpoisson(I_BPL_moyen_ecriture_scientifique,Nacquisition)
            for i in range(len(L_I_Bruit_Pollution_Lumineuse)):
                L_I_Bruit_Pollution_Lumineuse[i]*=pow(10,power_I_BPL_moyen)
            L_I_Bruit_Courant_Obscurite = hpoisson(I_BCO_moyen_ecriture_scientifique,Nacquisition)
            for i in range(len(L_I_Bruit_Courant_Obscurite)):
                L_I_Bruit_Courant_Obscurite[i]*=pow(10,power_I_BCO_moyen)

            # Section de calcul du signal bruité et non-bruité, enfin du RSB.
            for i in range(Nacquisition):
                I_Bruit_Photon_Signal = L_I_Bruit_Photon_Signal[i]
                I_Bruit_Pollution_Lumineuse = L_I_Bruit_Pollution_Lumineuse[i]
                I_Bruit_Courant_Obscurite = L_I_Bruit_Courant_Obscurite[i]

                IBPS_cumul+= I_Bruit_Photon_Signal
                IBPL_cumul+= I_Bruit_Pollution_Lumineuse
                IBCO_cumul+= I_Bruit_Courant_Obscurite

                # Calcul de l'intensité électrique résultante à un instant donné en prenant en compte les bruits
                I_SPUR_1 = I_Signal_Pur_1(I_Signal_Pur_0,alpha,L) #Intensitée lumineuse du singal pur après avoir traversé l'atmosphère
                I_LUM_1 = I_lumineuse_1(I_SPUR_1,I_Bruit_Photon_Signal,I_Bruit_Pollution_Lumineuse)
                I_LUM_2 = I_lumineuse_2(I_LUM_1,Dmin,Dmax,Hc,Techantillonnage,Dod,f,lambda_do)
                I_ELEC_1 = I_electrique_1(I_Saturation,Temperature,Techantillonnage,I_LUM_2,Surface_Capteur,U_seuil_photodiode_passante,I_Bruit_Courant_Obscurite)
                I_ELEC_3 = I_electrique_3(Techantillonnage,constante_conversion,ISO,I_ELEC_1)
                Ielec_bruitee_cumul+=I_ELEC_3

                # Calcul de l'intensité électrique résultante à un instant donné non-bruitée
                I_LUM_2 = I_lumineuse_2(I_SPUR_1,Dmin,Dmax,Hc,Techantillonnage,Dod,f,lambda_do)
                I_ELEC_1 = I_electrique_1(I_Saturation,Temperature,Techantillonnage,I_LUM_2,Surface_Capteur,U_seuil_photodiode_passante,0)
                I_ELEC_3 = I_electrique_3(Techantillonnage,constante_conversion,ISO,I_ELEC_1)
                Ielec_non_bruitee_cumul+= I_ELEC_3
            
            RSB = Ielec_non_bruitee_cumul/(Ielec_bruitee_cumul-Ielec_non_bruitee_cumul)
            L_RSB[i_ISO][i_Diametre_ouverture_diaphragme].append(RSB)

##########################################################################################################################################

fig1 = figure(1)
ax = Axes3D(fig1)
Lx = np.array(L_Texposition)
Ly = np.array(L_Diametre_ouverture_diaphragme)
Lx,Ly = np.meshgrid(Lx,Ly)
Matrix_2D = np.array(L_RSB[0])
ax.set_title("Carte du RSB selon {Texpo,Diam.ouv.diaph.}")
ax.set_xlabel("Texpo(s)")
ax.set_ylabel("Diamètre d'ouverture du diaphragme(m)")
ax.set_zlabel("RSB(ISO=L_ISO[0])")
ax.scatter(Lx, Ly, Matrix_2D, c = 'r', marker = 'o')
show()

### VERSION V0 QUI MARCHE:
##
##fig1 = figure(1)
##ax = Axes3D(fig1)
##Lx = np.array(L_Texposition)
##Ly = np.array(L_Diametre_ouverture_diaphragme)
##Lx,Ly = np.meshgrid(Lx,Ly)
##Matrix_2D = np.array(L_RSB[0])
##ax.set_title("Carte du RSB selon {Texpo,Diam.ouv.diaph.}")
##ax.set_xlabel("Texpo(s)")
##ax.set_ylabel("Diamètre d'ouverture du diaphragme(m)")
##ax.set_zlabel("RSB(ISO=L_ISO[0])")
##ax.plot_surface(Lx, Ly, Matrix_2D, rstride=1, cstride=1, cmap='Blues')
##show()

