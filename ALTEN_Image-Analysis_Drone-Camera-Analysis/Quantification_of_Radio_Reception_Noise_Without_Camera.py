# librairies
import matplotlib.pyplot as plt
import numpy
import cv2

# Lien vers un site web contenant beaucoup d'opérations de base détaillées : https://python.plainenglish.io/image-processing-using-opencv-in-python-857c8cb21767

def Identify_Dark_Current_Noise(frame):
        height, width = frame.shape[0:2]
        for i in range(height): #on parcourt toutes les lignes de l'image
                for j in range(width): #on parcourt toutes les colonnes de l'image
                        pixel = frame[i , j] # récupération des valeurs de rvb du pixel
                        red = pixel[0] # valeur de la composane bleue du pixel
                        green = pixel[1] # valeur de la composane verte du pixel
                        blue = pixel[2] # valeur de la composane rouge du pixel

                        sum_red+=red
                        sum_green+=green
                        sum_blue+=blue
        
##########################################################################
#charger la vidéo dans la variable cap
cap = cv2.VideoCapture('video_fpv_crazyflie.mp4')

ret, frame = cap.read()
#cv2.imshow('output', frame)
#cv2.imshow('Original Image', frame) 
#cv2.waitKey(0)

height, width = frame.shape[0:2]
N_pixels = height*width
print("Height=",height," and Width=",width)

L_blue = []
L_green = []
L_red = []


# boucle infinie (ou finie réglée sur n itérations)
N_frames_total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
print("N_frames = ",N_frames_total)

video_time_length = 104 # length of the video in seconds

t_start = 0
t_end = 104

frame_start = int(t_start*N_frames_total/video_time_length)
frame_end = int(t_end*N_frames_total/video_time_length)

Lx = [u for u in range(frame_start , frame_end+1 )]

N_frames_to_compute = frame_end - frame_start
print("N_frames_to_compute = ",N_frames_to_compute)

frame_current = frame_start
while (True and frame_current<=frame_end):
        # stoquer l'image issue de la vidéo à l'instant t dans la variable "frame"
        ret, frame = cap.read()
        # afficher le type des images contenues dans la vidéo (i.e. le type de l'objet "frame"), on attendra comme réponse : <class 'numpy.ndarray'>
        #print(type(frame))
        # afficher l'image contenue dans "frame"
        #cv2.imshow('output', frame)
        #cv2.waitKey(0)
        # quiter la boucle infinie lorqu'on appuie sur la touche 'q'

        sum_blue = 0
        sum_green = 0
        sum_red = 0
        for i in range(height): #on parcourt toutes les lignes de l'image
                for j in range(width): #on parcourt toutes les colonnes de l'image
                        [red,green,blue] = frame[i , j] # récupération des valeurs de rvb du pixel

                        sum_red+=red
                        sum_green+=green
                        sum_blue+=blue
                        
        sum_blue/=N_pixels # calcul de la valeur moyenne de bleue sur l'image
        sum_green/=N_pixels # calcul de la valeur moyenne de vert sur l'image
        sum_red/=N_pixels # calcul de la valeur moyenne de rouge sur l'image

        L_blue.append(sum_blue)
        L_green.append(sum_green)
        L_red.append(sum_red)

        print("N_Current_Frame =",frame_current)
        frame_current+=1
        if(cv2.waitKey(1) & 0xFF == ord('q')):
                break

plt.plot( Lx , L_red , color='r' , label='I_red mean value')
plt.plot( Lx , L_green , color='g' , label='I_green mean value')
plt.plot( Lx , L_blue , color='b' ,  label='I_blue mean value')
plt.xlabel('N_frame(t)')
plt.ylabel('Intensity(frame)')
plt.legend()
plt.grid(True)
plt.show()
 
#quiter le programme et fermer toutes les fenêtres ouvertes
cap.release()
cv2.destroyAllWindows()
