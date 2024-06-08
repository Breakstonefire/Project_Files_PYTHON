# librairies
import matplotlib.pyplot as plt
import numpy
import cv2

# Lien vers un site web contenant beaucoup d'opérations de base détaillées : https://python.plainenglish.io/image-processing-using-opencv-in-python-857c8cb21767

############################################################
#### COMMENTAIRE GENERAL ET DESCRIPTION DE CE QUE FAIT LE SCRIPT ####
##Dans ce programme, il est question de travailler sur une vidéo format mp4 prise par une caméra sur un drone.
##On peut sélectionner un intervalle continu d'images sur lesquelles on va calculer la valeur moyenne d'une des composantes RGB et ce pour chaque image.
##Les intensités moyennées de chaque image sont enfin affichées dans un plot.
##
##Coordonnées d'une image :
##
##Ligne numéro 0 en haut à gauche.
##Colonne numéro 0 tout à gauche.
##         __________________________
##        |                                 |     Colonne numéro *dernier* toute à droite.
##        |                                 |
##        |                                 |
##        |                                 |
##        |__________________________|
##
##Ligne numéro *dernier* tout en bas à gauche.
##
##########################################################################
#charger la vidéo dans la variable cap
cap = cv2.VideoCapture('video_fpv_crazyflie.mp4')
ret, frame = cap.read()

height, width = frame.shape[0:2]
N_pixels = height*width
print("Height=",height," and Width=",width)

N_frames_total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
print("N_frames = ",N_frames_total)

L_blue = []

# boucle infinie (ou finie réglée sur n itérations)
video_time_length = 104 # length of the video in seconds

t_start = 0 # ce paramètre permet de dire à quelle temps de la vidéo on souhaite lancer les calculs de moyennage d'intensité d'une image.
t_end = 103.9 # ATTENTION, un problème de "Nonetype Object" apparaît lorsqu'on inclu la denrière image n°2622 dans le calcul, ne comprenant pas comment résoudre cela, il est préférable de mettre un temps de fin de vidéo à traiter légèrement inférieur au réel...

frame_start = int(t_start*N_frames_total/video_time_length) # A NOTER que l'on peut également renseigner directement le numéro de frame voulue si on le souhaite plutôt que de donner un temps comem ci-dessus
frame_end = int(t_end*N_frames_total/video_time_length) # A NOTER que l'on peut également renseigner directement le numéro de frame voulue si on le souhaite plutôt que de donner un temps comem ci-dessus
print("Frames to compute from ",frame_start," to ",frame_end,".") # affichage du nombre de frames considérées pour les calculs
frame_current = frame_start # on initialise le numéro de la frame de départ (frame_start) dans la variable frame_current qui dictera dans la boucle while ci-dessosu à quelle frame on est.

N_frames_to_compute = frame_end - frame_start # affichage du nombre de frames considérées pour les calculs
print("N_frames_to_compute = ",N_frames_to_compute) # affichage du nombre de frames considérées pour les calculs

diminishing_explored_aera_Ratio = 13 # ce coefficient (compris entre 1 et minimum(height ; width) ) va venir diviser à la fois la valeur height et width des images de la vidéo, de sorte que l'on ne traite qu'un rectangle réduit (donc un nombre plus petit de pixels) sur chaque image, ce qui est moins couteux en calcul. Attention à ne pas trop élever ce nombre car on ne prendrait alors qu'un tout petit rectangle de l'image qui n'est plus caractéristique de l'intensité moyenne de chaque couleur.
new_height = int(height/diminishing_explored_aera_Ratio) # calcul de la hauteur du rectangle inscrit dans les images qui représente les pixels vraiment pris en compte pour la calcul
new_width = int(width/diminishing_explored_aera_Ratio)
N_pixels_new = new_height*new_width

while (True and frame_current<=frame_end):
        # stoquer l'image issue de la vidéo à l'instant t dans la variable "frame"
        ret, frame = cap.read()

        sum_blue = 0
        for i in range(new_height): #on parcourt toutes les lignes de l'image
                for j in range(new_width): #on parcourt toutes les colonnes de l'image
                        pixel = frame[i , j] # récupération des valeurs de rvb du pixel
                        blue = pixel[2]
                        sum_blue+=blue
                        
        sum_blue/=N_pixels_new # calcul de la valeur moyenne de bleue sur l'image
        L_blue.append(sum_blue)

        # print("N_Current_Frame =",frame_current) # L'affichage de cette phrase indiquant l'image en cours de traitement fait perdre beaucoup de temps à l'algorithme.
        frame_current+=1 # incrémentation

Lx = [u for u in range(frame_start , frame_end+1 )] # création de la liste contenant les numéros de frames traitées dans le programme

# Bloc d'affichage des résultats
plt.plot( Lx , L_blue , color='b' ,  label='I_blue mean value') # affichage de l'intensité moyenne de toutes les images traitées en fonction des images successives ( /en fonction du temps)
plt.xlabel('N_frame(t)')
plt.ylabel('Intensity(frame)')
plt.legend()
plt.grid(True)
plt.show()
 
#quiter le programme et fermer toutes les fenêtres ouvertes
cap.release()
cv2.destroyAllWindows()
