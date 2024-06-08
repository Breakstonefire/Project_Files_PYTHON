# DANS CE PROGRAMME L'ISO EST BIEN DEFINIT COMME SUIT:
# L'ISO EST LA SENSIBILITE DU FILM PHOTOGRAPHIQUE DONC SELON WIKIPEDIA
# C'EST LE RAPPORT H0/H IL FAUT FAIRE COMME SI ISO = H0*cste_d_offset DANS CE PROGRAMME

# Ce programme permet de tracer la courbe du rapport signal sur bruit en fonction des trois paramètres de l'expérience:
# 1. Le temps d'exposition totale allant de 1/3200 seconde à 15 secondes
# 2. Le diamètre d'ouverture du diaphragme (à distance focale fixée et constante ; pour des observations d'objets lointains à quasi-infini, la focale est maximale)
# 3. La sensibilité ISO (variant de 100 à 1600 pour la CANON SX 170)
# Le graphe final peut en réalité être représenté par trois graphes séparés donnant l'allure moyenne de la dépendance de chacun des paramètres en fixant les deux autres
# Ou bien être un ensemble de graphes 3D qui montrent la dép. du RSB aux deux autres paramètres (pour toutes leurs valeurs) en fixant le dernier paramètre (de préférence celui qui a le moins de valeurs possibles)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing libraries # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import division
from random import random
from numpy import *
from math import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from mpl_toolkits import axes_grid1
from matplotlib import cm

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
import pygame

#### EXCELLENT COURS SUR LA SENSITOMETRIE ####:
# http://www.photographeaparis.fr/tech-photo/index.htm
#### Caractéristiques de l'appareil photo ###
# https://www.canon.fr/for_home/product_finder/cameras/digital_camera/powershot/powershot_sx170_is/specification.html
#### Lien vers la règle NPF pour un calcul du temps d'exposition maximal théorique permettant d'éviter un filé de source lumineuse sur l'image lors de l'acquisition :
# http://web.archive.org/web/20200220123345/https://www.sahavre.fr/tutoriels/astrophoto/34-regle-npf-temps-de-pose-pour-eviter-le-file-d-etoiles

#### Fonctions du calcul principal ####
def I_Signal_Pur_1(I_Signal_Pur_0,alpha,L):
    #print("I_Signal_Pur_1 OKAY")
    return I_Signal_Pur_0*np.exp(-alpha*L)
def I_lumineuse_1(I_Signal_Pur_1,I_Bruit_Photon_Signal,I_Bruit_Pollution_Lumineuse):
    #print("I_lumineuse_1 OKAY")
    return I_Signal_Pur_1 + I_Bruit_Photon_Signal + I_Bruit_Pollution_Lumineuse
def I_lumineuse_2(I_lumineuse_1,Dod,maxi_Dod):
    return I_lumineuse_1*Dod/maxi_Dod
def I_lumineuse_3(I_lumineuse_2,Dmin,Dmax,Hc,Techantillonnage,Dod,f,lambda_do): # Hc EST EN FAIT L'ISO DANS LA SUITE DES CALCULS DU PROGRAMME
    # Cette fonction est une fonction logistique utilisée pour la sensitivité d'un film photographique lorsqu'on trace sa caractéristique (absorbance (ou densité optique) en fonction du log en base 10 de l'exposition (ou la lumination)). La page wikipedia "Fonction logistique (Verhulst)" l'explique parfaitement.
    H = f/(I_lumineuse_2*Techantillonnage*Dod) #calcul de la lumination ou aussi appelée exposition
    Lumination_ratio = Hc/H
    exposant = Dmin+(Dmax-Dmin)/(1+pow(Lumination_ratio,lambda_do/log(10)))
    #print("I_lumineuse_2 OKAY")
    return I_lumineuse_2/pow(10,exposant)
def I_electrique_1(I_Saturation,Temperature,Techantillonnage,I_lumineuse_3,Surface_Capteur,U_seuil_photodiode_passante,I_Bruit_Courant_Obscurite):
    kB = 1.380649/pow(10,23)
    e = 1.602/pow(10,19)
    I_elec_1_non_bruite = I_Saturation*((1/(kB*Temperature))*(Techantillonnage*I_lumineuse_3*Surface_Capteur - e*U_seuil_photodiode_passante)-1)
    #print("I_electrique_1 OKAY")
    if (I_elec_1_non_bruite<=0):
        return I_Bruit_Courant_Obscurite
    else:
        return I_elec_1_non_bruite + I_Bruit_Courant_Obscurite
def I_electrique_2(Techantillonnage,constante_conversion,I_electrique_1):
    #print("I_electrique_3 OKAY")
    return Techantillonnage*constante_conversion*I_electrique_1

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
def _hpoisson(m):
    ph = random.random()
    k = 0
    pc = poisson(k,m)
    while ph>pc:
        k+=1
        pc+=poisson(k,m)
    return k
def hpoisson(m,nb=0):
    """hpoisson(m,nb=0): Génération de valeurs tirées au hasard selon une distribution de Poisson"""
    if nb==0:
        return _hpoisson(m)
    else:
        R=[]
        for j in range(0,nb):
            R.append(_hpoisson(m))
        return R
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
    if (x>1 and x<10):
        return x,0
    if (x>=10):
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
def reduction_nombre_elems_liste(L,facteur_reduction):
    length = len(L)
    L_reduite = []
    if (length/facteur_reduction<1):
        L_reduite.append(L[0])
    else :
        cpt = 0
        delta_increment = 1+int(length/int(length/facteur_reduction))
        for i in range(length):
            if i==0:
                L_reduite.append(L[0])
            if (cpt==delta_increment):
                L_reduite.append(L[i])
                cpt = 0
            cpt+=1
    return L_reduite
def mean_z_Matrix_2D(Matrix_2D):
    sum_z = 0
    for iy in range(len(Matrix_2D)):
        for ix in range(len(Matrix_2D[0])):
            sum_z+=Matrix_2D[iy][ix]
    sum_z/=len(Matrix_2D)
    sum_z/=len(Matrix_2D[0])
    return sum_z
def find_zmin_zmax_in_3D_matrix(L_RSB,parametre_de_la_methode):
    length = len(L_RSB)
    zmin,zmax = 0,0
    if(parametre_de_la_methode=="z_extremes"):
        L_zmin = []
        L_zmax = []
        for i in range(length):
            Matrix_2D = L_RSB[i]
            RSB_max_ISO_fixe = max([max(Matrix_2D[u]) for u in range(len(Matrix_2D))])
            RSB_min_ISO_fixe = min([min(Matrix_2D[u]) for u in range(len(Matrix_2D))])
            L_zmax.append(RSB_max_ISO_fixe)
            L_zmin.append(RSB_min_ISO_fixe)
        zmin = min(L_zmin)
        zmax = max(L_zmax)
    else:
        coefficient_de_marge_supplementaire_au_dessus_en_dessous = parametre_de_la_methode
        L_zmean = []
        for i in range(length):
            Matrix_2D = L_RSB[i]
            zm_i = mean_z_Matrix_2D(Matrix_2D)
            L_zmean.append(zm_i)
        zmin = min(L_zmean)/coefficient_de_marge_supplementaire_au_dessus_en_dessous
        zmax = max(L_zmean)*coefficient_de_marge_supplementaire_au_dessus_en_dessous
    return zmin,zmax
def ajouter_2_listes(L1,L2):
    L = [0 for u in range(len(L2))]
    for i in range(len(L2)):
        L[i] = L1[i]+L2[i]
        if len(L)!=len(L1):
            print("WARNING - SHAPE MISMATCH in ajouter_2_listes")
    return L
def ajouter_2_matrices_dim_2(M1_2,M2_2):
    M_2 = [[0 for u in range(len(M2_2[0]))] for v in range(len(M2_2))]
    for i in range(len(M2_2)):
        M_2[i] = ajouter_2_listes(M1_2[i],M2_2[i])
        if len(M_2)!=len(M1_2):
            print("WARNING 2.0 - SHAPE MISMATCH in ajouter_2_matrices_dim_2")
        if len(M_2[0])!=len(M1_2[0]):
            print("WARNING 2.1 - SHAPE MISMATCH in ajouter_2_matrices_dim_2")
    return M_2
def ajouter_2_matrices_dim_3(M1_3,M2_3):
    M_3 = [[[0 for u in range(len(M2_3[0][0]))] for v in range(len(M2_3[0]))] for w in range(len(M2_3))]
    for i in range(len(M2_3)):
        M_3[i] = ajouter_2_matrices_dim_2(M1_3[i],M2_3[i])
        if len(M_3)!=len(M1_3):
            print("WARNING 3.0 - SHAPE MISMATCH in ajouter_2_matrices_dim_3")
        if len(M_3[0])!=len(M1_3[0]):
            print("WARNING 3.1 - SHAPE MISMATCH in ajouter_2_matrices_dim_3")
        if len(M_3[0][0])!=len(M1_3[0][0]):
            print("WARNING 3.2 - SHAPE MISMATCH in ajouter_2_matrices_dim_3")
    return M_3
def multiplier_par_cste_matrice_dim_3(M_3,cste):
    for i in range(len(M_3)):
        for j in range(len(M_3[0])):
            for k in range(len(M_3[0][0])):
                M_3[i][j][k]*=cste
    return M_3
def Moyennage_Maps_3D(List_of_RSB_maps_for_every_ISO):
    N = len(List_of_RSB_maps_for_every_ISO)
    L_RSB_moyennes = List_of_RSB_maps_for_every_ISO[0] #[[[0 for u in range(len(L_RSB[0][0]))] for v in range(len(L_RSB[0]))] for w in range(len(L_RSB))]
    for i in range(1,N):
        M2_3 = List_of_RSB_maps_for_every_ISO[i]
        L_RSB_moyennes = ajouter_2_matrices_dim_3(L_RSB_moyennes,M2_3)
    L_RSB_moyennes = multiplier_par_cste_matrice_dim_3(L_RSB_moyennes,1/N)
    return L_RSB_moyennes
def add_colorbar(im, aspect_parameter=30, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect_parameter)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)
def sum_xy(Lx,Ly):
    somme = 0
    for i in range(len(Lx)):
        xi = Lx[i]
        yi = Ly[i]
        somme+=xi*yi
    return somme
def fit_affine(Lx,Ly):
    over_n = 1/len(Lx)
    sx = sum(Lx)
    sy = sum(Ly)
    sx2 = sum_xy(Lx,Lx)
    #sy2 = sum_xy(Ly,Ly)
    sxy = sum_xy(Lx,Ly)
    xmean = sx*over_n
    ymean = sy*over_n
    covariance_xy = over_n*sxy-xmean*ymean
    variance_x = over_n*sx2-pow(xmean,2)
    a = covariance_xy/variance_x
    b = ymean-a*xmean
    return a,b
def fit_RSB_map(amin_x,amax_x,nax,amin_y,amax_y,nay,bmin,bmax,nb,Lx,Ly,L_RSB):
    len_x = len(Lx)
    len_y = len(Ly)
    Matrix_2D_fit = [[0 for u in range(len_x)] for v in range(len_y)]
    L_error = [[[0 for ib in range(nb)] for i_ay in range(nay)] for i_ax in range(nax)]
    for i_ax in range(nax):
        ax = amin_x+(amax_x-amin_x)*i_ax/(nax-1)
        for i_ay in range(nay):
            ay = amin_y+(amax_y-amin_y)*i_ay/(nay-1)
            for i_b in range(nb):
                b = bmin+(bmax-bmin)*i_b/(nb-1)
                error = 0
                for i in range(len_y):
                    y = Ly[i]
                    for j in range(len_x):
                        x = Lx[j]
                        RSB = L_RSB[i][j]
                        RSB_fit = ax*x+ay*y+b
                        e = pow(abs(RSB_fit-RSB),2)
                        error+= e
                L_error[i_ax][i_ay][i_b] = error
    error_min = L_error[0][0][0]
    iax_fit = 0
    iay_fit = 0
    ib_fit = 0
    for i_ax in range(nax):
        for i_ay in range(nay):
            for i_b in range(nb):
                error = L_error[i_ax][i_ay][i_b]
                if (error<error_min):
                    error_min = error
                    iax_fit = i_ax
                    iay_fit = i_ay
                    ib_fit = i_b
    ax_fit = amin_x+(amax_x-amin_x)*iax_fit/(nax-1)
    ay_fit = amin_y+(amax_y-amin_y)*iay_fit/(nay-1)
    b_fit = bmin+(bmax-bmin)*ib_fit/(nb-1)
    for i in range(len_y):
        y = Ly[i]
        for j in range(len_x):
            x = Lx[j]
            Matrix_2D_fit[i][j] = ax_fit*x+ay_fit*y+b_fit
    print("ax=",ax_fit)
    print("ay=",ay_fit)
    print("b=",b_fit)
    return Matrix_2D_fit,[ax_fit,ay_fit,b_fit]
def moyenne_liste(L):
    return sum(L)/len(L)
########################################################################################################
#### DEBUT DU PROGRAMME - PARAMETRES - CALCUL DES INTENSITES - CALCUL DES FITS DE PLANS - AFFICHAGE ####
#### Listes et constantes ####
MODE = "all"
#MODE = "fast"
# MODE = "very_fast" #"all" le mode fast permet de calculer la pile de maps de rsb pour un ensemble réduit de paramètre contrairement au mode "all" où toutes les combinaisons de paramètrs sont prises en compte

#parametre_de_la_methode = 1+0.0001 # Ce coefficient (>=1) permet d'abaisser le zmin autant que l'on fait augmenter le zmax, on divise le zmin précédent par ce coeff et on multiplie le zmax précédent par ce coeff.
parametre_de_la_methode = "z_extremes" # Si ce parametre vaut "z_extremes" alors on cherchera le zmin;zmax (ou le RSBminRSBmax) dans tous les plans pour tout ISO lors de l'affichage final des résultats

Nb_realisations_simulation = 100 # nombre de fois où l'on simule les maps du RSB en fonction de l'ISO, du temps d'exposition et du diamètre d'ouverture du diaphragme, le RSB de chaque map 3D est alors moyenné en tout point sur ce nombre de réalisation
nax = 100 #nombre de fois où on calcule l'erreur entre le fit et le plan expérimental à b et ay fixés dans l'équation de plan du fit : P_fit(x,y) = ax*x+ay*y+b
nay = 100 #nombre de fois où on calcule l'erreur entre le fit et le plan expérimental à ax et b fixés dans l'équation de plan du fit : P_fit(x,y) = ax*x+ay*y+b
nb = 100 #nombre de fois où on calcule l'erreur entre le fit et le plan expérimental à ax et ay fixés dans l'équation de plan du fit : P_fit(x,y) = ax*x+ay*y+b
# => Le nombre d'erreurs calculées est donc nax*nay*nb

affichage_sur_fenetres_multiples = "oui"
affichage_sur_fenetres_multiples = "non"

##############################
L_ISO = [100,200,400,800,1600]
L_Texposition = [1/3200,1/2500,1/2000,1/1600,1/1250,1/1000,1/800,1/640,1/500,1/400,1/320,1/250,1/200,1/160,1/125,1/100,1/80,1/60,1/50,1/40,1/30,1/25,1/20,1/15,1/13,1/10,1/8,1/6,1/5,1/4,0.3,0.4,0.5,0.6,0.8,1,1.3,1.6,2,2.5,3.2,4,5,6,8,10,13,15]
L_Diametre_ouverture_diaphragme = [1/3.5,1/4,1/4.5,1/5,1/5.6,1/6.3,1/7.1,1/8]

facteur_reduction = 1
L_Texposition = reduction_nombre_elems_liste(L_Texposition,facteur_reduction)
L_Texposition = [1,1.3,1.6,2,2.5,3.2,4,5,6,8,10,13,15]

if MODE=="all":
    print("Go take a coffee you have time...")
if MODE=="fast":
    L_ISO = [100,1600]
    facteur_reduction = 10
    #L_Texposition = [1/6,1/5,1/4,0.3,0.4,0.5,0.6,0.8,1,1.3,1.6,2,2.5,3.2,4,5,6,8,10,13,15]
    L_Texposition = reduction_nombre_elems_liste(L_Texposition,facteur_reduction)
    facteur_reduction = 4
    L_Diametre_ouverture_diaphragme = reduction_nombre_elems_liste(L_Diametre_ouverture_diaphragme,facteur_reduction)
if MODE=="very_fast":
    L_ISO = [100]#,1600]
    facteur_reduction = 16
    #L_Texposition = [1/6,1/5,1/4,0.3,0.4,0.5,0.6,0.8,1,1.3,1.6,2,2.5,3.2,4,5,6,8,10,13,15]
    L_Texposition = reduction_nombre_elems_liste(L_Texposition,facteur_reduction)
    facteur_reduction = 3
    L_Diametre_ouverture_diaphragme = reduction_nombre_elems_liste(L_Diametre_ouverture_diaphragme,facteur_reduction)

maxi_Dod = max(L_Diametre_ouverture_diaphragme)

##############################
#### DEBUT DES PARAMETRES ####
# LES PARAMETRES SONT DONNES DANS L'ORDRE CHRONOLOGIQUE LORSQU'ON POSE LE SYSTEME GLOBAL SUR PAPIER ET QU'ON L'ETUDIE DEPUIS LE DEBUT.

#grandissement = 1 # le grandissement est la valeur du zoom donc >=1
I_Signal_Pur_0 = 0.10                                   # ESTIMABLE C'est l'intensité du signal pur sous une aurore boréale environ aussi lumineuse que la moitié de la lune.
lambda_lumiere = 550/pow(10,9)                          # ESTIMABLE On prend la longueur d'onde du vert (~550nm) car ce sont les aurores courantes entre 100 et 300km d'altitude. La longueur d'onde varie car les couleurs des aurores dépendent de l'altitude et de la longueur d'onde.
h_cst_Planck = 6.62607015/pow(10,34)                    # CSTE PHY
c_celerite_lumiere = 299792458                          # CSTE PHY vitesse de la lumiere dans le vide (m/s)
n_indice_refraction_air = 1.000293                      # ESTIMABLE adimensionné
v_lumiere = c_celerite_lumiere/n_indice_refraction_air  # ESTIMABLE
alpha = 0.00001                                         # ESTIMABLE coefficient d'absorption (m^-1) du milieu dans lequel se propage l'onde électromagnétique depuis l'objet jusqu'à l'objectif. Lien vers une page web utile pour les valeurs de l'air dans toute l'atmosphère : http://ressources.univ-lemans.fr/AccesLibre/UM/Pedago/physique/02/optigeo/loidebeer.html
L = 100*pow(10,3)                                       # ESTIMABLE Distance entre l'aurore boréale et le capteur photographique. Les aurores sont à une altitude variable allant de 100 à 600km d'altitude (se référer à la couleur pour connaitre l'altitude moyenne)
focal_min = 5/1000
focal_max = 80/1000
f = focal_max                                           # PARAMETRE focal(L,grandissement) # la valeur de la focal dépend du zoom
Texposition = 1                                         # PARAMETRE temps (sec) pendant lequel le(/les) capteurs restent ouverts et aptent de la lumière et sont donc excités et envoient un signal électrique au reste du circuit électronique
Techantillonnage = min(L_Texposition)                   # ESTIMABLE c'est le temps le plus court durant lequel l'appareil peut prendre une image (a priori il est <= au temps d'obturation minimal = 1/3200 sec, on prendra donc cette valeur en référence)
Nacquisition = int(Texposition/Techantillonnage)        # ESTIMABLE c'est le nombre de fois où l'appareil va acquérir une image (peu lumineuse ou non) pendant tout le temps d'exposition "Texposition" en prenant une image tous les "Techantillonnage"
Dmin = 0                                                # ESTIMABLE On suppose que lorsque l'intensité lumineuse incidente tend vers 0 (donc la lumination de même), on a une absorbance (ou donc une densité optique) du capteur qui tend vers 0 puisqu'aucun photon n'est absorbé (hormis ceux du bruit mais le Dmin est exprimé dans un cas idéal d'exposition lumineuse).
Dmax = 1                                                # INCONNUE C'est la valeur maximale d'absorbance que le capteur adopte lorsqu'il est soumis à des intensités lumineus très élevées.
# Hc est calculé à partir des paramètres précédents et de ceux qui suivent (cf ci-dessous pour les détails du calcul
Dod = 0.01                                              # INCONNUE 
lambda_do = 1                                           # INCONNUE
I_Saturation = 10/pow(10,9)                             # Pour la valeur d'application numérique on se base sur les transistors  à jonctions polarisées (PN et/ou NP) où on trouve des Isat de l'ordre de 10 nA), on peut se référer à ce lien : http://res-nlp.univ-lemans.fr/NLP_C_M15_G01/co/Contenu_c4.html
Temperature = 273.15 + 2                                # Température en degrés Kelvin ici on prend 2°C pour une nuit d'hiver classique en Novembre en Norvège
Surface_Capteur = 4.363117258846936e-12                 # Surface (en m²) élémentaire d'un capteur unitaire du capteur CCD de l'appareil photo.
U_seuil_photodiode_passante = 1
constante_conversion = 1                                # Cette constante est une constante qui permet de conserver l'unité de l'intensité électrique lorsqu'on multiplie l'intensité électrique en sortie de capteur par le gain après être passé dans l'amplificateur analogique. étant une simple constante permettant ce respect des unités, on la prend égale à 1.
ISO = L_ISO[0]                                          # Les valeurs de l'ISO sont dans la liste "L_ISO" et caractérisent la sensibilité du film photographique composé par tout les capteurs du CCD.
# Section d'explication de la valeur (ordre de grandeur de Hc) : d'après les formules qui suivent, on peut exprimer Hc (=H0) en fonction des autres paramètres.
#1)a) ISO = H0 / H
#1)b) H0  = ISO * H =~ ISO * Hthéorique_idéal
#2)a) Hthéorique_idéal = ISO * Texpo * Dod / f (ISO mal dimensionné pour le problème si on garde cette formule telle quelle)
#2)b) Hthéorique_idéal = Imoyenne * ISO * Texpo * Dod / f (avec Imoyenne en W/m² une intenisté moyenne sur le temps d'exposition Texpo)
#3)   H0 =~ ISO * Hthéorique_idéal = ISO * Imoyenne * ISO * Texpo * Dod /f
#=> H0 =~ Imoyenne * ISO² * Texpo * Dod / f
#	 (W/m²)   * (-)  * (sec) * (m) / (m) = W/m² * sec = lux*sec (c'est bien l'intensité de la lumination)
Hc = I_Signal_Pur_1(I_Signal_Pur_0,alpha,L)*pow(moyenne_liste(L_ISO),2)*moyenne_liste(L_Texposition)*moyenne_liste(L_Diametre_ouverture_diaphragme)/f  # ESTIMABLE C'est la valeur de référence de la lumination telle que l'image est sous éclairement "optimale" à un ISO donné
Hc_sur_ISO = I_Signal_Pur_1(I_Signal_Pur_0,alpha,L)*moyenne_liste(L_Texposition)*moyenne_liste(L_Diametre_ouverture_diaphragme)/f
L_Hc_fction_de_ISO = [Hc_sur_ISO*u for u in(L_ISO)]
# Explication supplémentaire : On considère le signal moyen après traversée de l'atmosphère (donc absorption partielle)
# on considère la valeur moyenne des ISO comme valeur de référence pour le caclul de H0 (ou Hc) de même pour Dod
print("Hc(ISO=",L_ISO[0],")=",L_Hc_fction_de_ISO[0])

#### bruits et photocourants parasites :
# Quelques liens vers les bruits présents dans une photodiode : https://fr.lambdageeks.com/photo-detector-photo-diode/
# https://media4.obspm.fr/public/ressources_lu/pages_analyser/bruit-photons_impression.html
#### Point sur les intensités lumineuses de référence de nuit :
#   Éclairement lumineux    Exemple
#   0,25 lux	            Pleine lune par une nuit claire
#   0,01 lux	            Quartier de lune
#   0,002 lux	            Ciel étoilé sans lune
I_BPS_moyen = I_Signal_Pur_0                            # Intensité moyenne (W/m²) du bruit de photon du signal. On prend ici le signal pur lui même (ou la racine carrée du signal pur) pour avoir une approximation correcte du bruit
I_BPL_moyen = 0.05                                      # Intensité moyenne (W/m²) du bruit de pollution lumineuse. On prend l'équivalent d'un éclairement de nuit sous un quartier de lune.
I_BCO_moyen = I_Saturation                              # Intensité électrqiue moyenne du bruit dû au courant d'obscurité, qui est un courant résiduel et parasite de notre signal électrique 

#L_Parametres_constantes = [focal_min,focal_max,L_ISO,L_Texposition,facteur_reduction,L_Diametre_ouverture_diaphragme,I_Signal_Pur_0,lambda_lumiere,h_cst_Planck,c_celerite_lumiere,n_indice_refraction_air,v_lumiere,alpha,L,f,Texposition,Techantillonnage,Nacquisition,Dmin,Dmax,Hc,Dod,lambda_do,I_Saturation,Temperature,Surface_Capteur,U_seuil_photodiode_passante,constante_conversion,ISO,I_BPS_moyen,I_BPL_moyen,I_BCO_moyen,n_sigma]
#### FIN DES PARAMETRES ####
############################

######################################
#### DEBUT DES CALCULS DES BRUITS ####
# Calcul des Intensités des différents bruits selon leurs processus de poisson respectifs
I_BPS_moyen_ecriture_scientifique,power_I_BPS_moyen = Ecriture_scientifique(I_BPS_moyen)
I_BPL_moyen_ecriture_scientifique,power_I_BPL_moyen = Ecriture_scientifique(I_BPL_moyen)
I_BCO_moyen_ecriture_scientifique,power_I_BCO_moyen = Ecriture_scientifique(I_BCO_moyen)
Ielec_bruitee_cumul = 0
Ielec_non_bruitee_cumul = 0
IBPS_cumul = 0
IBPL_cumul = 0
IBCO_cumul = 0

L_IBPS_cumul = []
L_IBPL_cumul = []
L_IBCO_cumul = []
L_I_elec_3 = []
L_I_elec_3_cumul = []
L_I_elec_3_non_bruitee = []
L_I_elec_3_non_bruitee_cumul = []
Lt = []
L_RSB = [[[0 for u in L_Texposition] for v in L_Diametre_ouverture_diaphragme] for w in L_ISO] #Cette liste contiendra tous les RSB pour tous les bruits simulés pour tous les paramètres de Texpo et Dod à ISO fixé
List_of_RSB_maps_for_every_ISO = [] #cette liste contiendra les différentes map L_RSB pour chacun des ISO
L_RSB_non_bruite = [[[0 for u in L_Texposition] for v in L_Diametre_ouverture_diaphragme] for w in L_ISO] #Cette liste est l'analogue de L_RSB mais pour un signal non-bruite de bout en bout (le RSB est calculé comme Signal/Bruit (où Bruit=constante=1)
List_of_RSB_maps_for_every_ISO_non_bruite = [] #cette liste contiendra toutes les maps de RSB pour un bruit unitaire constant (ie. un signal non bruite de bout en bout)

#List_of_RSB_maps_for_every_ISO = [[[[] for u in L_Diametre_ouverture_diaphragme] for v in L_ISO] for w in range(Nb_realisations_simulation)]
#L_Grandeurs_Calculees = [I_BPS_moyen_ecriture_scientifique,power_I_BPS_moyen,I_BPL_moyen_ecriture_scientifique,power_I_BPL_moyen,I_BCO_moyen_ecriture_scientifique,power_I_BCO_moyen,L_IBPS_cumul,L_IBPL_cumul,L_IBCO_cumul,L_I_elec_3,L_I_elec_3_cumul,L_I_elec_3_non_bruitee,L_I_elec_3_non_bruitee_cumul,Lt,Ielec_bruitee_cumul,Ielec_non_bruitee_cumul,IBPS_cumul,IBPL_cumul,IBCO_cumul,L_RSB]

#### FIN DES CALCULS DES BRUITS ####
####################################

###########################
#### DEBUT DES CALCULS ####
print("Lancement des calculs...")
pygame.init() #cette ligne permet de démarrer le chronomètre pour mesurer le temps que la section algorithmique ci-dessous met pour s'effectuer, le temps s'affiche avec l'écriture de la ligne qui suit : "print(pygame.time.get_ticks())"

#### SECTION DE CALCUL DE L'INTENSITE NON BRUITEE ####
len_ISO = len(L_ISO)
for i_ISO in range(len_ISO):
    ISO = L_ISO[i_ISO]
    Hc = L_Hc_fction_de_ISO[i_ISO]
    print("Hc=",Hc)
    I_Bruit_Photon_Signal = 0
    I_Bruit_Pollution_Lumineuse = 0
    I_Bruit_Courant_Obscurite = 0
    I_SPUR_1 = I_Signal_Pur_1(I_Signal_Pur_0,alpha,L) #Intensitée lumineuse du singal pur après avoir traversé l'atmosphère
    I_LUM_1 = I_lumineuse_1(I_SPUR_1,I_Bruit_Photon_Signal,I_Bruit_Pollution_Lumineuse)
    I_LUM_2 = I_lumineuse_2(I_LUM_1,Dod,maxi_Dod)
    for i_Diametre_ouverture_diaphragme in range(0,len(L_Diametre_ouverture_diaphragme)):
        Dod = L_Diametre_ouverture_diaphragme[i_Diametre_ouverture_diaphragme]
        I_LUM_3 = I_lumineuse_3(I_LUM_2,Dmin,Dmax,Hc,Techantillonnage,Dod,f,lambda_do)
        I_ELEC_1 = I_electrique_1(I_Saturation,Temperature,Techantillonnage,I_LUM_3,Surface_Capteur,U_seuil_photodiode_passante,I_Bruit_Courant_Obscurite)
        I_ELEC_2 = I_electrique_2(Techantillonnage,constante_conversion,I_ELEC_1)
        for i_Texposition in range(0,len(L_Texposition)):
            Texposition = L_Texposition[i_Texposition]
            Nacquisition = int(Texposition/Techantillonnage)
            Ielec_non_bruitee_cumul = I_ELEC_2*Nacquisition
            L_RSB_non_bruite[i_ISO][i_Diametre_ouverture_diaphragme][i_Texposition] = Ielec_non_bruitee_cumul

print("Fin des calculs.")
print("Calcul terminé en ",int(pygame.time.get_ticks()/1000)," secondes")
#### FIN DES CALCULS ####
#########################

###########################################################################################################
#### THE FOLLOWING BLOCK CAN BE COPY-PASTE IN EVERY CODE TO PLOT K colored curves composed of N VALUES ####
#RED TO BLACK
# K =     # value corresponding to the number of 'y(x)' curves to plot
# N =     # value corresponding to the number of points to plot for each curve (= len(x))
Lcolors = [[]]
for i in range(1276):
    Lcolors.append([0,0,0])
# Theoretically the number of triplet-color-values in Lcolors is equal to 5*(2^8-1) = 1275
# Color is RED
# INCREASING THE VALUE of GREEN
index = 0
for green in range(0,256,1):
    rgb = [255,green,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# Color is YELLOW
# DECREASING THE VALUE of RED
for red in range(254,-1,-1):
    rgb = [red,255,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# Color is GREEN
# INCREASING THE VALUE of BLUE
for blue in range(1,256,1):
    rgb = [0,255,blue]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# Color is CYAN
# DECREASING THE VALUE of GREEN
for green in range(254,-1,-1):
    rgb = [0,green,255]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# Color is BLUE
# DECREASING THE VALUE of BLUE
for blue in range(254,0,-1):
    rgb = [0,0,blue]
    Lcolors[index] = rgb
    index+=1
# Color is BLACK
# DELETING THE DOUBLE VALUES IN THE END
for i in range(len(Lcolors)-1,0,-1):
    l1 = Lcolors[i]
    l2 = Lcolors[i-1]
    if l1==l2:
        print(" ")
        print("WARNING when computing the colors for plots: doublon at index",i)
        print(" ")
        del Lcolors[i]
del Lcolors[-1]
len_col = len(Lcolors)
# Theoretically the number of triplet-color-values in Lcolors is equal to 5*(2^8-1) = 1275
# Then we normalize the values of Lcolors from [0:255] to [0:1]
for i in range(len(Lcolors)):
    for j in range(len(Lcolors[0])):
        Lcolors[i][j]/=255
Lcolors.reverse()
#### END OF THE GENERIC BLOCK ####
##################################

##############################
#### DEBUT DE L'AFFICHAGE ####
pygame.init() #cette ligne permet de démarrer le chronomètre pour mesurer le temps que la section algorithmique ci-dessous met pour s'effectuer, le temps s'affiche avec l'écriture de la ligne qui suit : "print(pygame.time.get_ticks())"
# Début de la section de remplissage des élémentes de la figure affichée
print("Démarrage de l'affichage des résultats")

len_x = len(L_RSB_non_bruite[0])
len_y = len(L_RSB_non_bruite[0][0])
Lx = np.array(log(L_Texposition))
dx = (Lx[len(Lx)-1]-Lx[0])
Ly = np.array(L_Diametre_ouverture_diaphragme)
dy = (Ly[len(Ly)-1]-Ly[0])
Lx,Ly = np.meshgrid(Lx,Ly)

if(affichage_sur_fenetres_multiples=="oui"):
    for i in range(len(L_RSB_non_bruite)):
        if(len(L_RSB_non_bruite)==1):
            i_color = 0
        else:
            i_color = int((i-len_RSB)*(len(Lcolors)-1)/(len(L_RSB_non_bruite)-1))
        color_for_plot = Lcolors[i_color]
        fig1 = figure(1)
        ax = Axes3D(figure(i+1))
        # A StrMethodFormatter is used automatically
        ax.zaxis.set_major_formatter('{x:.02f}')
        # Set the titles and legends
        ax.set_title("Map RSB_Non_Bruite=f(ISO,Texpo.,D.o.d.)")
        ax.set_xlabel("ln(Texpo(s))")
        ax.set_ylabel("Diamètre d'ouverture du diaphragme(m)")
        ax.set_zlabel("RSB(ISO)")
        # Customize the z axis.
        zmin,zmax = find_zmin_zmax_in_3D_matrix(L_RSB_non_bruite,"z_extremes")
        print("Signal PUR -> [zmin;zmax]=[",zmin,";",zmax,"]")
        ax.set_zlim(zmin,zmax)
        ax.zaxis.set_major_locator(LinearLocator(10))
        Matrix_2D = np.array(L_RSB_non_bruite[i-len_RSB]) # Plot the surface. Pour modifier le type de gradient de couleurs qu'on souhaite avoir, il faut remplacer le cm.Greys ou cm.coolwarm (par exemple) par cm."***" le truc qui vous intéresse sur le lien suivant: https://matplotlib.org/stable/tutorials/colors/colormaps.html
        ax.plot_surface(Lx, Ly, Matrix_2D, cmap=cm.coolwarm, linewidth=0, antialiased=False) #ax.plot_surface(Lx, Ly, Matrix_2D, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.scatter(Lx, Ly, Matrix_2D, color=color_for_plot, marker = 'o')
        ax.plot_surface(Lx,Ly,Matrix_2D, cmap=cm.Greys, linewidth=0, antialiased=False) #ax.plot_surface(Lx, Ly, Matrix_2D, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.scatter(Lx, Ly, Matrix_2D, color=color_for_plot, marker = 'o')
        # Add a color bar which maps the different values to colors.
        # fig1.colorbar(Matrix_2D, shrink=0.5, aspect=5)
        c = np.array([0,1])
        norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
        color_map_to_use = (mpl.colors.ListedColormap(Lcolors).with_extremes(over=Lcolors[0], under=Lcolors[len(Lcolors)-1]))
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=color_map_to_use)
        cmap.set_array([])
        # Set the colorbar
        colorbar_legend = "RSB LOW to HIGH ->->->"
        fig1.colorbar(cmap, location='right',orientation='vertical', label=colorbar_legend, aspect = 20 , pad = 0.05, shrink = 0.9, extendfrac='auto', spacing='uniform')
else:
#### FIGURE 2 ####
    fig2 = figure(2)
    ax = Axes3D(fig2)
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')
    # Set the titles and legends
    ax.set_title("Map RSB_Non_Bruite=f(ISO,Texpo,Dod)")
    ax.set_xlabel("ln(Texpo(s))")
    ax.set_ylabel("Diamètre d'ouverture du diaphragme(m)")
    ax.set_zlabel("RSB(ISO)")
    # Customize the z axis.
    zmin,zmax = find_zmin_zmax_in_3D_matrix(L_RSB_non_bruite,"z_extremes")
    print("Signal PUR -> [zmin;zmax]=[",zmin,";",zmax,"]")
    ax.set_zlim(zmin,zmax)
    ax.zaxis.set_major_locator(LinearLocator(10))
    for i in range(len(L_RSB_non_bruite)):
        if(len(L_RSB_non_bruite)==1):
            i_color = 0
        else:
            i_color = int(i*(len(Lcolors)-1)/(len(L_RSB_non_bruite)-1))
        color_for_plot = Lcolors[i_color]
        Matrix_2D = np.array(L_RSB_non_bruite[i]) # Plot the surface. Pour modifier le type de gradient de couleurs qu'on souhaite avoir, il faut remplacer le cm.Greys ou cm.coolwarm (par exemple) par cm."***" le truc qui vous intéresse sur le lien suivant: https://matplotlib.org/stable/tutorials/colors/colormaps.html
        ax.plot_surface(Lx, Ly, Matrix_2D, cmap=cm.coolwarm, linewidth=0, antialiased=False) #ax.plot_surface(Lx, Ly, Matrix_2D, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.scatter(Lx, Ly, Matrix_2D, color = color_for_plot, marker = 'o')
        ax.plot_surface(Lx, Ly, Matrix_2D, cmap=cm.Greys, linewidth=0, antialiased=False) #ax.plot_surface(Lx, Ly, Matrix_2D, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.scatter(Lx, Ly, Matrix_2D, color=color_for_plot, marker = 'o')
    # Add a color bar which maps the different values to colors.
    # fig2.colorbar(Matrix_2D, shrink=0.5, aspect=5)
    c = np.array([zmin,zmax])
    norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
    color_map_to_use = (mpl.colors.ListedColormap(Lcolors).with_extremes(over=Lcolors[0], under=Lcolors[len(Lcolors)-1]))
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=color_map_to_use)
    cmap.set_array([])
    # Set the colorbar
    colorbar_legend = "RSB LOW to HIGH ->->->"
    fig2.colorbar(cmap) #, location='right',orientation='vertical', label=colorbar_legend, aspect = 20 , pad = 0.05, shrink = 0.9, extendfrac='auto', spacing='uniform')

print("Affichage réalisé en ",int(pygame.time.get_ticks()/1000)," secondes")
show()
#### FIN DE L'AFFICHAGE ####
############################
