from math import *
print(' ')

#### HYPOTHESES de calcul ####
# 1. Le plan de focalisation (zone de l'espace où un objet paraît net mais =/= du plan focal objet ou image) de l'appareil est littéralement plan et pas courbé.
# 2. L'objet se déplace en ligne droite dans le plan de focalisation.
# 3. L'objet se déplace à vitesse constante.
# 4. L'objet a une taille (/largeur en 1D) correspondant exactement à la taille (/largeur) de la zone de l'espace décrite par le pixel d'intérêt, cette zone étant donc évidemment comprise dans le plan de focalisation.

# Nombre de pixels en format LARGE de photo (16 millions de pixels) : 4608 de large par 3456 de haut
NombrePixels_X    = 4608 # Nombre de pixels selon l'axe X de l'appareil photo lorsqu'il est posé à plat et regarde selon l'axe Y
NombrePixels_Z    = 3456 # Nombre de pixels selon l'axe Z de l'appareil photo lorsqu'il est posé à plat et regarde selon l'axe Y
delta_carre_pixel   = 1/100 / NombrePixels_X # approximation le temps d'avoir la valeur exacte du capteur de l'appareil
delta_theta_deg     = 100 # approximation avant la valeur exacte de l'angle d'ouverture de l'appareil photo <=> vision de l'appareil photo // vision d'un être humain ~130° en degrés
DistanceObjet        = 50 # Distance en mètres à l'objet qui est variable selon le cas de figure
VitesseObjet           = 10 # Vitesse de l'objet en mouvement en m/s

RecouvrementGlobal_pourcent = 50 # Pourcentage de recouvrement global du pixel d'intérêt 

# Calculs automatiques supplémentaires (conversions / variables intermédiaires / etc ...)
delta_theta_rad = pi * delta_theta_deg / 180 # conversion de l'angle d'ouverture de l'appareil photo en radians
RecouvrementGlobal_taux = 1 - (RecouvrementGlobal_pourcent / 100) # calcul du taux de recouvrement global de la zone de l'espace décrite par le pixel d'intérêt à partir du pourcentage de recrouvrement global.

# Calcul du temps d'exposition maximal / limite permettant toujours de s'affranchir du flou de mouvement sous les hyôthèses décrites au dessus.
TempsExpositionMaximal = 4 * tan(delta_theta_rad / 2) * DistanceObjet * RecouvrementGlobal_taux / (NombrePixels_X * VitesseObjet)

# Affichage des caractéristiques propres de l'appareil
print(' - - - - CARACTERISTIQUES DE L''APPAREIL - - - -')
message = "Nombre de pixels selon largeur  de l'image = {}".format(NombrePixels_X)
print(message)
message = "Nombre de pixels selon hauteur de l'image = {}".format(NombrePixels_Z)
print(message)
message = "Angle d'ouverture de l'appareil photo          = {}°".format(delta_theta_deg)
print(message)
print(' ')

# Affichage des paramètres de l'expérience
print(' - - - - PARAMETRES DE L''EXPERIENCE - - - -')
message = "Distance à l'objet = {} [m]".format(DistanceObjet)
print(message)
message = "Vitesse de l'objet = {} [m/s]".format(VitesseObjet)
print(message)
print(' ')

# Boucle for d'identification du dénominateur entier donnant la valeur la plus proche du temps d'exposition maximal à ne pas dépasser
#   L'identification se fait par le bas donc le temps d'exposition sous forme de fraction d'entiers donné est le plus haut temps d'exposition inférieur ou égal au temps d'exposition maximal
if TempsExpositionMaximal >= 1 :
    message = ''
else :
    for i in range(1, 20001) :
        if (1/(i + 1) <= TempsExpositionMaximal) and (TempsExpositionMaximal <= 1/i) :
            vitesse_obturation_reglage_appareil = i+1
            break

if TempsExpositionMaximal >= 0.3 :
    TempsExpositionMaximal_message = round(TempsExpositionMaximal , 2)
    message = "Le temps d'exposition maximal pour éviter tout flou de mouvement vaut {} [secondes]"\
              .format(TempsExpositionMaximal_message)
elif TempsExpositionMaximal >= 0.001 :
    TempsExpositionMaximal_message = round(TempsExpositionMaximal * 1000 , 3)
    message = "Le temps d'exposition maximal pour éviter tout flou de mouvement vaut {} [millisecondes] soit 1/{} [secondes]"\
              .format(TempsExpositionMaximal_message , vitesse_obturation_reglage_appareil)
else :
    TempsExpositionMaximal_message = round(TempsExpositionMaximal * 1000000 , 0)
    message = "Le temps d'exposition maximal pour éviter tout flou de mouvement vaut {} [microsecondes] soit 1/{} [secondes]"\
              .format(TempsExpositionMaximal_message , vitesse_obturation_reglage_appareil)

print(message)
## ## ##
