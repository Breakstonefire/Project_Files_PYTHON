# librairies
import matplotlib.pyplot as plt
import numpy
import cv2
import matplotlib.animation as animation

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

L_red = []
L_green = []
L_blue = []

# boucle infinie (ou finie réglée sur n itérations)
video_time_length = 104 # length of the video in seconds

t_start = 7
t_end = 10 # ATTENTION, un problème de "Nonetype Object" apparaît lorsqu'on inclu la denrière image n°2622 dans le calcul, ne comprenant pas comment résoudre cela, il est préférable de mettre un temps de fin de vidéo à traiter légèrement inférieur au réel...

frame_start = int(t_start*N_frames_total/video_time_length)
frame_end = int(t_end*N_frames_total/video_time_length)

frame_start = 0
frame_end = 2620

print("Frames to compute from ",frame_start," to ",frame_end,".")
frame_current = frame_start

N_frames_to_compute = frame_end - frame_start
print("N_frames_to_compute = ",N_frames_to_compute)

diminishing_explored_aera_Ratio = 1
new_height = int(height/diminishing_explored_aera_Ratio)
new_width = int(width/diminishing_explored_aera_Ratio)
N_pixels_new = new_height*new_width

while (True and frame_current<=frame_end):

        if(frame_current == 100):
                print("Frame 100")
        if(frame_current == 200):
                print("Frame 200")
        if(frame_current == 300):
                print("Frame 300")
        if(frame_current == 400):
                print("Frame 400")
        if(frame_current == 500):
                print("Frame 500")
        
        if(frame_current == 1000):
                print("Frame 1000")
        
        if(frame_current == 1500):
                print("Frame 1500")
        
        if(frame_current == 2000):
                print("Frame 2000")
        
        if(frame_current == 2500):
                print("Frame 2500")

        # stoquer l'image issue de la vidéo à l'instant t dans la variable "frame"
        ret, frame = cap.read()

        sum_red = 0
        sum_green = 0
        sum_blue = 0
        for i in range(new_height): #on parcourt toutes les lignes de l'image
                for j in range(new_width): #on parcourt toutes les colonnes de l'image
                        pixel = frame[i , j] # récupération des valeurs de rvb du pixel
                        red = pixel[0]
                        green = pixel[1]
                        blue = pixel[2]
                        sum_red+=red
                        sum_green+=green
                        sum_blue+=blue

        sum_red/=N_pixels_new # calcul de la valeur moyenne de rouge sur l'image
        sum_green/=N_pixels_new # calcul de la valeur moyenne de vert sur l'image
        sum_blue/=N_pixels_new # calcul de la valeur moyenne de bleue sur l'image
        L_red.append(sum_red)
        L_green.append(sum_green)
        L_blue.append(sum_blue)

        #print("N_Current_Frame =",frame_current) # L'affichage de cette phrase indiquant l'image en cours de traitement fait perdre beaucoup de temps à l'algorithme.
        frame_current+=1

Lx = [u for u in range(frame_start , frame_end+1 )]

with open("R_list.txt", 'w') as output:
    output.write("[")
    for row in L_red:
        output.write(str(row) + ',')
    output.write("]")

with open("G_list.txt", 'w') as output:
    output.write("[")
    for row in L_green:
        output.write(str(row) + ',')
    output.write("]")

with open("B_list.txt", 'w') as output:
    output.write("[")
    for row in L_blue:
        output.write(str(row) + ',')
    output.write("]")

with open("Lframes_list.txt", 'w') as output:
    output.write("[")
    for row in Lx:
        output.write(str(row) + ',')
    output.write("]")

# setting values to rows and column variables
rows = 1
columns = 1
# create figure
fig , ax1 = plt.subplots(nrows = rows, ncols = columns , layout = "constrained")
# showing plots
ax1.plot( Lx , L_red , color='r' , label='I_red mean value')
ax1.plot( Lx , L_green , color='g' , label='I_green mean value')
ax1.plot( Lx , L_blue , color='b' ,  label='I_blue mean value')
ax1.set_title("RGB mean intensity in time")
ax1.set_xlabel("Frame(t)")
ax1.set_ylabel("Intensity")
ax1.grid(True)
plt.show()

#quiter le programme et fermer toutes les fenêtres ouvertes
cap.release()
cv2.destroyAllWindows()
