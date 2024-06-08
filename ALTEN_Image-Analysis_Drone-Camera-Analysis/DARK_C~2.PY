# librairies
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import cv2
import os

import tkinter as tk
from tkinter import *
from tkinter import messagebox as mb

# Lien vers un site web contenant beaucoup d'opérations de base détaillées : https://python.plainenglish.io/image-processing-using-opencv-in-python-857c8cb21767

def video_to_frames(input_loc, output_loc):
    # function to exctrat all the frames from a video and to save them into a given location on the PC
    try:
        os.mkdir(output_loc)
    except OSError:
        pass
    # Log the time
    time_start = time.time()
    # Start capturing the feed
    cap = cv2.VideoCapture(input_loc)
    # Find the number of frames
    video_length = int(cap.get(cv2.CAP_PROP_FRAME_COUNT)) - 1
    print ("Number of frames: ", video_length)
    count = 0
    print ("Converting video..\n")
    # Start converting the video
    while cap.isOpened():
        # Extract the frame
        ret, frame = cap.read()
        if not ret:
            continue
        # Write the results back to output location.
        cv2.imwrite(output_loc + "/%#05d.jpg" % (count+1), frame)
        count = count + 1
        # If there are no more frames left
        if (count > (video_length-1)):
            # Log the time again
            time_end = time.time()
            # Release the feed
            cap.release()
            # Print stats
            print ("Done extracting frames.\n%d frames extracted" % count)
            print ("It took %d seconds forconversion." % (time_end-time_start))
            break

def save_one_frame(frame_number,name_to_save_string):
    try:
        os.mkdir(output_loc)
    except OSError:
        pass
    #charger la vidéo dans la variable cap
    cap = cv2.VideoCapture('video_fpv_crazyflie.mp4')
    ret, frame = cap.read()
    # Write the results back to output location.
    cv2.imwrite(name_to_save_string, frame)
    # Release the feed
    cap.release()

def No_Saturation_Verification(frame,Intensity_Saturation_Threshold,ratio_of_saturated_pixels_authorized):
    height, width = frame.shape[0:2]
    N_pixels = height*width
    N_pixels_sat_in_red = 0
    N_pixels_sat_in_green = 0
    N_pixels_sat_in_blue = 0
    for i in range(height): #on parcourt toutes les lignes de l'image
        for j in range(width): #on parcourt toutes les colonnes de l'image
            pixel = frame[i , j] # récupération des valeurs de rvb du pixel
            if(pixel[0]>=Intensity_Saturation_Threshold):
                N_pixels_sat_in_red+=1
            if(pixel[1]>=Intensity_Saturation_Threshold):
                N_pixels_sat_in_green+=1
            if(pixel[2]>=Intensity_Saturation_Threshold):
                N_pixels_sat_in_blue+=1
    if(N_pixels_sat_in_red > ratio_of_saturated_pixels_authorized*N_pixels):
        return False
    if(N_pixels_sat_in_green > ratio_of_saturated_pixels_authorized*N_pixels):
        return False
    if(N_pixels_sat_in_blue > ratio_of_saturated_pixels_authorized*N_pixels):
        return False
    return True

def find_mins_values(L,List_of_data_intervals_to_process_in_L,input_location,output_location,Intensity_Threshold):
    # Supprime les fichiers 'txt' et 'png' éventuellement déjà présents dûs aux précédentes exécutions de cet algorithme :
    for files in output_location: #os.listdir(output_location):
        if (files.endswith(".png") or files.endswith("txt")):
            os.remove(os.path.join(output_location, files))
    
    #charger la vidéo dans la variable cap et la première image dans la variable frame
    cap = cv2.VideoCapture('video_fpv_crazyflie.mp4')
    ret, frame = cap.read()
    # récupération des données caractériqtiques des images
    height, width = frame.shape[0:2]

    len_L = len(L) # stockage de la longueur de la liste des intensités données en paramètre 
    len_Lintervals = len(List_of_data_intervals_to_process_in_L) #stockage de la longueur de la liste des différents intervalles contenant les séquences d'images à analyser ci-après
    List_of_index = [] # création d'une liste permettant de garder en mémoire les numéros des images que l'on va enregistrer pour les traiter dans un autre programme éventuellement
    for i in range(len_Lintervals): # on parcourt tous les intervalles d'images à étudier
        List_of_current_frames_numbers_to_process = List_of_data_intervals_to_process_in_L[i] # on parcourt l'intervalle d'images renseigné par 'List_of_data_intervals_to_process_in_L'
        first_frame_number_to_process = List_of_current_frames_numbers_to_process[0] # on note le numéro de la première image de la séquence à traiter
        last_frame_number_to_process = List_of_current_frames_numbers_to_process[1] # pareil pour la dernière
        Intensity_min = min(L[first_frame_number_to_process : last_frame_number_to_process + 1]) # on ajoute +1 pour prendre en compte le dernier indice jusqu'auquel on cherche l'intensité minimale
        if(Intensity_min<=Intensity_Threshold): # on applique ici un critère de seuil maximal d'intensité moyenne de l'image à ne pas dépasser pour être considére comme une image contenant potentiellement uniquement du bruit de courant d'obscurité
            frame_number_of_min_Intensity = L.index(Intensity_min, first_frame_number_to_process, last_frame_number_to_process + 1) # on ajoute +1 pour prendre en compte le dernier indice jusqu'auquel on cherche l'intensité minimale
            cap.set(1,frame_number_of_min_Intensity) # on lit l'image en cours si le critère est validé
            ret, frame = cap.read() # stoquer l'image issue de la vidéo à l'instant t dans la variable "frame"
            if(No_Saturation_Verification(frame,Intensity_Saturation_Threshold,ratio_of_saturated_pixels_authorized)==True):
                cv2.imshow('frame_'+str(frame_number_of_min_Intensity)+'_to_show', frame) # on affiche l'image contenue dans "frame"
                res=mb.askquestion('Image selection.', 'Do you want to keep the image ?') # une fenêtre s'ouvre pour demander si l'image doit être gardée après inspection de l'image par l'utilisateur
                if (res == 'yes') :
                    cv2.imwrite(os.path.join(output_location , 'frame_'+str(frame_number_of_min_Intensity)+'.png'), frame) # on enregistre l'image dans le dossier renseigné avec le nom spécifié
                    List_of_index.append(frame_number_of_min_Intensity) # on ajoute le numéro de l'image fraichement traitée à la liste des numéros d'images gardées
                    print("Imean___saved =",Intensity_min) # on affiche à titre indicatif l'intensité moyenne de l'image en indiquant sur le script python en cours que c'est une image gardée
                else :
                    print("Imean_unsaved =",Intensity_min) # on affiche à titre indicatif l'intensité moyenne de l'image en indiquant sur le script python en cours que c'est une image rejetée
                cv2.destroyWindow('frame_'+str(frame_number_of_min_Intensity)+'_to_show') # on ferme enfin la fenêtre qui s'était ouverte avec l'image en cours de traitement qui était dedans
                cv2.waitKey(0) # cette commande doit venir après la précédente pour faire en sorte que la fenetre reste ouverte tant qu'on n'intéragit pas avec al boite de dialogue qu is'est ouverte précédemment
    with open(output_location+'\\'+'List_of_index.txt', 'w') as output: # on ouvre un fichier dont le nom est renseigné directement et dont la localisation est 'output_location' pour garder en mémoire la liste des images à traiter par la suite
        output.write(str(List_of_index)) # on écrit dans un fichier dont le nom est renseigné directement et dont la localisation est 'output_location' pour garder en mémoire la liste des images à traiter par la suite
    return List_of_index # la liste des numéros des images s'affiche également sur le prompteur python

def call():
    res=mb.askquestion('Image selection.', 'Do you want to keep the image ?')
    if res == 'yes' :
        root.destroy()
    else :
        mb.showinfo('Return', 'Following computations incoming...')
    return None

###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
# Paramètre permettant de pré-sélectionner les images candidates qui contiennent uniquement des couleurs dûes au bruti de courant d'obscurité de la caméra
Intensity_Threshold = 70
Intensity_Saturation_Threshold = 180
ratio_of_saturated_pixels_authorized = 0.1

input_location = r'C:\Users\lilia\Documents\0 ECOLE\5 CV STAGES EMPLOIS\4 Emplois\ALTEN\Projet_intermission_1_Essaim_Drones\Scripts_Python'
output_location_red = r'C:\Users\lilia\Documents\0 ECOLE\5 CV STAGES EMPLOIS\4 Emplois\ALTEN\Projet_intermission_1_Essaim_Drones\Frames_Dark_Current_Noise_Quantification\Frames_low_Red_identification'
output_location_green = r'C:\Users\lilia\Documents\0 ECOLE\5 CV STAGES EMPLOIS\4 Emplois\ALTEN\Projet_intermission_1_Essaim_Drones\Frames_Dark_Current_Noise_Quantification\Frames_low_Green_identification'
output_location_blue = r'C:\Users\lilia\Documents\0 ECOLE\5 CV STAGES EMPLOIS\4 Emplois\ALTEN\Projet_intermission_1_Essaim_Drones\Frames_Dark_Current_Noise_Quantification\Frames_low_Blue_identification'

file_red = open('R_list.txt', 'r')
file_green = open('G_list.txt', 'r')
file_blue = open('B_list.txt', 'r')
file_frames = open('Lframes_list.txt', 'r')

L_red = file_red.readlines()
L_green = file_green.readlines()
L_blue = file_blue.readlines()

L_red = L_red[0].replace('[','').replace(']','').split(",")
L_green = L_green[0].replace('[','').replace(']','').split(",")
L_blue = L_blue[0].replace('[','').replace(']','').split(",")

L_red = np.array(L_red)
L_green = np.array(L_green)
L_blue = np.array(L_blue)

L_red = L_red.astype(np.float64)
L_green = L_green.astype(np.float64)
L_blue = L_blue.astype(np.float64)

List_of_data_intervals_to_process_in_L = []

#charger la vidéo dans la variable cap
cap = cv2.VideoCapture('video_fpv_crazyflie.mp4')
ret, frame = cap.read()
height, width = frame.shape[0:2]
N_pixels = height*width
print("Height=",height," and Width=",width)
N_frames_total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
print("N_frames = ",N_frames_total)
video_time_length = 104 # length of the video in seconds

delta_t = 50 # paramètre décrivant la durée en temps de chaque séquence sur laquelle sera pré-sélectionnée UNE image selon les critères décrit plus haut.
delta_epsilon_t = 0.1 # permet de ne pas avoir de chevauchement entre al dernière image de l'intervalle précédent et la première image de l'intervalle suivant
N_interval = 1+(video_time_length+delta_epsilon_t) / (delta_t + delta_epsilon_t)
N_interval = int(N_interval)
print("N_interval =",N_interval)
for i in range(N_interval):
    t_start = i*(delta_t+delta_epsilon_t)
    t_end = t_start + delta_t
    frame_start = int(t_start*N_frames_total/video_time_length)
    frame_end = int(t_end*N_frames_total/video_time_length)
    List_of_data_intervals_to_process_in_L.append([frame_start,frame_end])

print(List_of_data_intervals_to_process_in_L)

L_index_red     = find_mins_values(L_red.tolist() , List_of_data_intervals_to_process_in_L,input_location , output_location_red,Intensity_Threshold)
#L_index_green = find_mins_values(L_green.tolist() , List_of_data_intervals_to_process_in_L,input_location , output_location_green)
#L_index_blue    = find_mins_values(L_blue.tolist() , List_of_data_intervals_to_process_in_L,input_location , output_location_blue)
