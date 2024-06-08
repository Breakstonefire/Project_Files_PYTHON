# librairies
import matplotlib.pyplot as plt
import numpy
import cv2

# Lien vers un site web contenant beaucoup d'opérations de base détaillées : https://python.plainenglish.io/image-processing-using-opencv-in-python-857c8cb21767

##########################################################################
#charger la vidéo dans la variable cap
cap = cv2.VideoCapture('video_fpv_crazyflie.mp4')
ret, frame = cap.read()

height, width = frame.shape[0:2]
N_pixels = height*width
print("Height=",height," and Width=",width)

N_frames_total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
print("N_frames = ",N_frames_total)

L_green = []

# boucle infinie (ou finie réglée sur n itérations)
video_time_length = 104 # length of the video in seconds

t_start = 0
t_end = 103.9 # ATTENTION, un problème de "Nonetype Object" apparaît lorsqu'on inclu la denrière image n°2622 dans le calcul, ne comprenant pas comment résoudre cela, il est préférable de mettre un temps de fin de vidéo à traiter légèrement inférieur au réel...

frame_start = int(t_start*N_frames_total/video_time_length)
frame_end = int(t_end*N_frames_total/video_time_length)
print("Frames to compute from ",frame_start," to ",frame_end,".")
frame_current = frame_start

N_frames_to_compute = 1 + frame_end - frame_start
print("N_frames_to_compute = ",N_frames_to_compute)

diminishing_explored_aera_Ratio = 13
new_height = int(height/diminishing_explored_aera_Ratio)
new_width = int(width/diminishing_explored_aera_Ratio)
N_pixels_new = new_height*new_width

while (True and frame_current<=frame_end):
        # stoquer l'image issue de la vidéo à l'instant t dans la variable "frame"
        ret, frame = cap.read()

        sum_green = 0
        for i in range(new_height): #on parcourt toutes les lignes de l'image
                for j in range(new_width): #on parcourt toutes les colonnes de l'image
                        pixel = frame[i , j] # récupération des valeurs de rvb du pixel
                        green = pixel[1]
                        sum_green+=green

        sum_green/=N_pixels_new # calcul de la valeur moyenne de vert sur l'image
        L_green.append(sum_green)

        # print("N_Current_Frame =",frame_current) # L'affichage de cette phrase indiquant l'image en cours de traitement fait perdre beaucoup de temps à l'algorithme.
        frame_current+=1

Lx = [u for u in range(frame_start , frame_end+1 )]

plt.plot( Lx , L_green , color='g' , label='I_green mean value')
plt.xlabel('N_frame(t)')
plt.ylabel('Intensity(frame)')
plt.legend()
plt.grid(True)
plt.show()
 
#quiter le programme et fermer toutes les fenêtres ouvertes
cap.release()
cv2.destroyAllWindows()
