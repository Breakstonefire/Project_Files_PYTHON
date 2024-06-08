# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 02:13:44 2024

@author: lilia
"""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # CLEANING / CLEARING / REMOVING / DELETING EVERYTHING BEFORE EXECUTION # #

g = globals() # Inititalizing g with globals()
d = dir() # Initializing d with dir() - this stores a list of all the variables in the program
for obj in d : # Checking user-defined variables in the directory
  if not obj.startswith('__') : # Checking for built-in variables/functions
    del globals()[obj] # Deleting the said obj, since a user-defined function

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing Libraries # # # # # # # # # # # # # # # # # # # # # #

import os
import csv
import time
import scipy
import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt

from scipy        import signal
from numpy.fft    import fft
from mpl_toolkits import mplot3d
from scipy        import integrate

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # REMOVING PLOTS BEFORE EXECUTION # # # # # # # # # # # # # # # # #
plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Functions # # # # # # # # # # # # # # # # # # # # # # # # # # #

def convert_data_from_str_to_float(x) :
  
  BOOL_Exposant = False
  
  for i in range(len(x)) :
    
    if x[i] == 'E' and BOOL_Exposant == False :
      BOOL_Exposant = True
    
    elif x[i] == 'E' and BOOL_Exposant == True :
      print('\nERROR :\n\tImpossible conversion, two or more exponents ...')
      return None
  
  if BOOL_Exposant == False : return float(x)
  
  elif BOOL_Exposant == True :
    for i in range(len(x)) :
      if x[i] == 'E' or x[i] == 'e' :
        return float(x[0:i]) * pow(10 , float(x[i+1 : len(x)]))

def read_csv(filename , delimiter , rows_predicted , limit_value_to_error_ACCELERATION , STARTING_LINE_FORCED , STOPPING_LINE_FORCED) :
  
  if delimiter == ',' : delimiter_2 = ';'
  else                : delimiter_2 = ','
      
  LTime     = []
  Lacc_x    = []
  Lacc_y    = []
  Lacc_z    = []
  Labs_acc  = []
  
  with open(filename, 'r', newline='') as csvfile :
    
    csv_reader = csv.reader(csvfile, delimiter = delimiter)
    BOOL_first_line = True
    row_number = -1
    previous_percent_shown_reading_csv = 0
    
    for row in csv_reader :  
      row_number += 1
      process_pourcent_reading_csv       = int(max(0 , 100 * (row_number - STARTING_LINE_FORCED)/(STOPPING_LINE_FORCED - STARTING_LINE_FORCED)))
      if int(previous_percent_shown_reading_csv) < int(process_pourcent_reading_csv) :
        print("\nPROCESSING CSV READING : Line {}/{} = {}% d''avancement".format(row_number , STOPPING_LINE_FORCED , process_pourcent_reading_csv))
        previous_percent_shown_reading_csv = int(process_pourcent_reading_csv)
      
      if BOOL_first_line == True :
        
        cpt_rows            = 0
        ind_delim_previous  = -1
        Lindexdelimiters    = [ind_delim_previous]
        Ltitle_tmp          = []
        
        for num_carac in range(0 , len(row[0]) , 1) :
            
          if ((row[0][num_carac] == delimiter) or (row[0][num_carac] == delimiter_2)) :
          
            Lindexdelimiters.append(num_carac)
            title_tmp = row[0][ind_delim_previous + 1 : num_carac]
            Ltitle_tmp.append(title_tmp)
            ind_delim_previous = num_carac
            cpt_rows = cpt_rows + 1
            print("\nTitle = {} ; Index of {} (or {}) = {}".format(title_tmp , delimiter , delimiter_2 , num_carac))
          
          elif (num_carac == len(row[0]) - 1) :
          
            title_tmp = row[0][ind_delim_previous + 1 : len(row[0])]
            Ltitle_tmp.append(title_tmp)
            cpt_rows = cpt_rows + 1
        
        if cpt_rows == rows_predicted :
          
          Title_Time , Title_acc_x , Title_acc_y , Title_acc_z , Title_abs_acc = Ltitle_tmp[:]
        
        else :
          
          print('\nWARNING :\n\tNumber of rows detected = {} =/= {}'.format(cpt_rows , rows_predicted))
          
        BOOL_first_line = False
    
      elif STARTING_LINE_FORCED - 1 <= row_number and row_number <= STOPPING_LINE_FORCED - 1 :
          
        cpt_rows            = 0
        cpt_convert_success = 0
        cpt_convert_error   = 0
        ind_delim_previous  = -1
        Lindexdelimiters    = [ind_delim_previous]
        LValues_str_tmp     = []
        
        for num_carac in range(0 , len(row[0]) , 1) :
          
          if ((row[0][num_carac] == delimiter) or (row[0][num_carac] == delimiter_2)) :
            
            Lindexdelimiters.append(num_carac)
            Values_str_tmp = row[0][ind_delim_previous + 1 : num_carac]
            LValues_str_tmp.append(Values_str_tmp)
            ind_delim_previous = num_carac
            cpt_rows = cpt_rows + 1
          
          elif (num_carac == len(row[0]) - 1) :
          
            Values_str_tmp = row[0][ind_delim_previous + 1 : len(row[0])]
            LValues_str_tmp.append(Values_str_tmp)
            cpt_rows = cpt_rows + 1
        
        if cpt_rows == rows_predicted :
          
          Time , acc_x , acc_y , acc_z , abs_acc = LValues_str_tmp[:]
          
          try :
            Time = convert_data_from_str_to_float(Time)
            LTime.append(Time)
            cpt_convert_success += 1
          except :
            print('\nWARNING :\n\tImpossible conversion from string to float of Time n°{} : {}'.format(row_number , Time))
            LTime.append(limit_value_to_error_ACCELERATION)
            cpt_convert_error += 1
          
          try :
            acc_x = convert_data_from_str_to_float(acc_x)
            Lacc_x.append(acc_x)
            cpt_convert_success += 1
          except :
            print('\nWARNING :\n\tImpossible conversion from string to float of acc_x n°{} : {}'.format(row_number , acc_x))
            Lacc_x.append(limit_value_to_error_ACCELERATION)
            cpt_convert_error += 1
          
          try :
            acc_y = convert_data_from_str_to_float(acc_y)
            Lacc_y.append(acc_y)
            cpt_convert_success += 1
          except :
            print('\nWARNING :\n\tImpossible conversion from string to float of acc_y n°{} : {}'.format(row_number , acc_y))
            Lacc_y.append(limit_value_to_error_ACCELERATION)
            cpt_convert_error += 1
          
          try :
            acc_z = convert_data_from_str_to_float(acc_z)
            Lacc_z.append(acc_z)
            cpt_convert_success += 1
          except :
            print('\nWARNING :\n\tImpossible conversion from string to float of acc_z n°{} : {}'.format(row_number , acc_z))
            Lacc_z.append(limit_value_to_error_ACCELERATION)
            cpt_convert_error += 1
          
          try :
            abs_acc = convert_data_from_str_to_float(abs_acc)
            Labs_acc.append(abs_acc)
            cpt_convert_success += 1
          except :
            print('\nWARNING :\n\tImpossible conversion from string to float of abs_acc n°{} : {}'.format(row_number , abs_acc))
            Labs_acc.append(limit_value_to_error_ACCELERATION)
            cpt_convert_error += 1
            
        else :
          print('\nWARNING :\n\tNumber of rows detected = {} =/= {}'.format(cpt_rows , rows_predicted))
        
        if cpt_convert_error > 0 :
          print('\nWARNING : (Reading CSV file line N°{})'.format(row_number))
          print('\n\tConversion errors      = {}'  .format(cpt_convert_error))
          print('\n\tSuccessful conversions = {}'.format(cpt_convert_success))
  
  print("\nCSV READING DONE ! Line {} = 100% d''avancement".format(row_number))
  
  return LTime , Lacc_x , Lacc_y , Lacc_z , Labs_acc , Title_Time , Title_acc_x , Title_acc_y , Title_acc_z , Title_abs_acc

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Main () # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#############
## PARAMETERS
print("\n\n## SETTING MAIN PARAMETERS")
Name_of_CSV_File  = 'Raw Data.csv'
# Name_of_CSV_File  = 'Raw Data Test.csv'
# FolderName        = 'Acclration sans g 2024-04-10 19-47-43'
# FolderName        = 'Acclration sans g 2024-04-10 19-47-43'
# FolderName        = 'Acclration sans g 2024-03-20 13-13-00'
FolderName          = 'COURSE MARDI 30 AVRIL 2024 CERGY ETANGS CLEAN'
# FolderName          = 'EXTRAIT_DATA_TEST'
filename          = os.path.join(r'C:\Users\lilia\Desktop\Course à pied CERGY' , FolderName , Name_of_CSV_File)
delimiter         = ','
rows_predicted    = 5
# limit_value_to_error_TIME         = 100
limit_value_to_error_ACCELERATION = 1000 # limit value (lower and upper in absolute value) that probably could indicate an error in the measure
# limit_value_to_error_ACCELERATION   = np.nan

# Paramètres de sélection de la plage de temps sur laquelle on prélève les données à traiter (rognage du début / fin des données)
dt_time_step_estimated       = 0.004781 # estimated mean time step in seconds
STARTING_DELAY_SEC_estimated = 0 # estimated time before the interesting data outcomes from the experiment
STOPPING_DELAY_SEC_estimated = 3700 # estimated time after which there is no more interesting data from the experiment
STARTING_LINE_FORCED         = int(STARTING_DELAY_SEC_estimated / dt_time_step_estimated) # Deduced parameter from the upper ones
STOPPING_LINE_FORCED         = int(STOPPING_DELAY_SEC_estimated / dt_time_step_estimated) # Deduced parameter from the upper ones

## Paramètres généraux de la FFT
# FFT Sampling number - nombre de points d'échantillonnage
Ne = pow(10,6) # nombre de points d'échantillonnage

# Mode choisi pour gérer la plage de fréquences à afficher des graphes FFT
# MODE_FFT_Frequences = "auto" # pas de limitation imposée (par défaut on tronque automatiquement la moitié droite du graphe de FFT car symétrie verticale à freq_max/2 théoriquement
MODE_FFT_Frequences = "Min_Max" # Du Min au Max paramétrés par l'utilisateur

# Paramètres d'affichage des fréquences min et max des graphes de FFT
freq_min_FFT_to_show = 0
freq_max_FFT_to_show = 1.5

# Mode choisi pour limiter l'affichage des trop grandes amplitudes des graphes FFT
MODE_FFT_Amplitudes = "auto" # pas de limitation imposée
# NE FONCTIONNE PAS ENCORE # MODE_FFT_Amplitudes = "Tout sauf la plus grande" # De la plus petite à la 2eme plus grande valeur
# MODE_FFT_Amplitudes = "Min_Max" # Du Min au Max paramétrés par l'utilisateur

# Paramètres d'affichage des amplitudes min et max des graphes de FFT
AmpMinParameter_FFT_to_show = -2
AmpMaxParameter_FFT_to_show = 10000

## Initial conditions : Position
x0    = 0
y0    = 0
z0    = 0
OM_0  = np.sqrt(x0**2 + y0**2 + z0**2)

# Initial conditions : Velocity
vx_0  = 0
vy_0  = 0
vz_0  = 0
V_0   = np.sqrt(vx_0**2 + vy_0**2 + vz_0**2)

# Color plots               
color_fig1_1 = [1 , 0 , 0]
color_fig1_2 = [0 , 1 , 0]
color_fig1_3 = [0 , 0 , 1]
color_fig1_4 = [0 , 0 , 0]

# Text sizes
SMALL_SIZE  = 4
MEDIUM_SIZE = 6
BIGGER_SIZE = 8
DefaultLineWidth = 0.2

# Number of bins in the histogram
number_of_bins             = 400
critere_range_histogram    = 'auto' #'std' # 'PourcentageOfValuesHistogram'
nsigma                     = 3 # number of std substracted/added that caracterize the x-axis range
pourcentageValuesHistogram = 95 # value in %

# Paramètres d'interpolation des mesures partiellement récupérées (Temps absent mais accélérations trouvées par ex. ...)
BOOL_Interpolate_Data = True

# Paramètres concernant le lissage / moyennage des données
BOOL_Smooth_Data                        = True # False
chosen_window_size_for_smoothing_filter = 14 # window size used for filtering - Nombre de points pris par moyennage
chosen_order_for_smoothing_filter       = 1 # order of fitted polynomial
# Une période d'oscillation de course semble être d'environ 0.7 secondes soit ~140 points de mesure par période
# => pour conserver un minimum de 10 points par période après lissage on peut donc moyenner 14 points par 14 au maximum
# => FIXONS LE MAX DE POINTS DE LA FENETRE DE LISSAGE A 13 (par sécurité)
# => LARGEUR FENETRE LISSAGE ~= 13 * 0.004781 (dt) = 0.062 secondes = 62 ms
# => TimePeriod_Smooth_ms = 60 # [ms] Période prise comme largeur temporelle de la fenêtre de lissage en [ms]

# Plots parameters
plt.rc('font'   , size      = SMALL_SIZE)       # controls default text sizes
plt.rc('axes'   , titlesize = SMALL_SIZE)       # fontsize of the axes title
plt.rc('axes'   , labelsize = SMALL_SIZE)       # fontsize of the x and y labels
plt.rc('xtick'  , labelsize = SMALL_SIZE)       # fontsize of the tick labels
plt.rc('ytick'  , labelsize = SMALL_SIZE)       # fontsize of the tick labels
plt.rc('legend' , fontsize  = SMALL_SIZE)       # legend fontsize
plt.rc('figure' , titlesize = BIGGER_SIZE)      # fontsize of the figure title
plt.rc('lines'  , linewidth = DefaultLineWidth)

####################
## Main computations
print("\n## STARTING MAIN COMPUTATIONS")
# Saving the starting time
T_START_PROGRAM = time.time()

## Lecture du fichier CSV
print("\n\n # Reading CSV file")
LTime , Lacc_x , Lacc_y , Lacc_z , Labs_acc , Title_Time , Title_acc_x , Title_acc_y , Title_acc_z , Title_abs_acc \
  = read_csv(filename , delimiter , rows_predicted , limit_value_to_error_ACCELERATION , STARTING_LINE_FORCED , STOPPING_LINE_FORCED)
LTime    = np.array(LTime)
Lacc_x   = np.array(Lacc_x)
Lacc_y   = np.array(Lacc_y)
Lacc_z   = np.array(Lacc_z)
Labs_acc = np.array(Labs_acc)

###################################
## PARTIE INTERPOLATION DES DONNEES
if BOOL_Interpolate_Data == False :
  print("\n\n # No data interpolation\n")
else :
  print("\n\n # Intepolating data")
  Lindex_to_interp_TIME = []
  Lindex_to_interp_accX    = []
  Lindex_to_interp_accY    = []
  Lindex_to_interp_accZ    = []
  Lindex_to_interp_abs_acc = []
  for i in range(0,len(LTime),1) :
    if LTime[i]    == np.nan : Lindex_to_interp_TIME.append(i)
    if Lacc_x[i]   == np.nan : Lindex_to_interp_accX.append(i)
    if Lacc_y[i]   == np.nan : Lindex_to_interp_accY.append(i)
    if Lacc_z[i]   == np.nan : Lindex_to_interp_accZ.append(i)
    if Labs_acc[i] == np.nan : Lindex_to_interp_abs_acc.append(i)
  
  # INTERPOLATING TIME array
  print("\n # # Intepolating Time data")
  Linterpolated_TIME = np.interp(Lindex_to_interp_TIME , np.arange(0 , len(LTime) , 1) , LTime)
  cpt_ind_interp_TIME = 0
  for i in range(0 , len(LTime) , 1) :
    if LTime[i] == np.nan :
      LTime[i] = Linterpolated_TIME[cpt_ind_interp_TIME]
      cpt_ind_interp_TIME += 1
  print("\n{}/{}  values interpolated in TIME array".format(cpt_ind_interp_TIME , len(Lindex_to_interp_TIME)))
  
  # INTERPOLATING X acceleration array
  print("\n # # Intepolating acc X data")
  Linterpolated_accX = np.interp(Lindex_to_interp_accX , LTime , Lacc_x)
  cpt_ind_interp_accX = 0
  for i in range(0 , len(Lacc_x) , 1) :
    if Lacc_x[i] == np.nan : 
      Lacc_x[i] = Linterpolated_accX[cpt_ind_interp_accX]
      cpt_ind_interp_accX += 1
  print("\n{}/{}  values interpolated in accX array".format(cpt_ind_interp_accX , len(Lindex_to_interp_accX)))
  
  # INTERPOLATING Y acceleration array
  print("\n # # Intepolating acc Y data")
  Linterpolated_accY = np.interp(Lindex_to_interp_accY , LTime , Lacc_y)
  cpt_ind_interp_accY = 0
  for i in range(0 , len(Lacc_y) , 1) :
    if Lacc_y[i] == np.nan :
      Lacc_y[i] = Linterpolated_accY[cpt_ind_interp_accY]
      cpt_ind_interp_accY += 1
  print("\n{}/{}  values interpolated in accY array".format(cpt_ind_interp_accY , len(Lindex_to_interp_accY)))
  
  # INTERPOLATING Z acceleration array
  print("\n # # Intepolating acc Z data")
  Linterpolated_accZ = np.interp(Lindex_to_interp_accZ , LTime , Lacc_z)
  cpt_ind_interp_accZ = 0
  for i in range(0 , len(Lacc_z) , 1) :
    if Lacc_z[i] == np.nan : 
      Lacc_z[i] = Linterpolated_accZ[cpt_ind_interp_accZ]
      cpt_ind_interp_accZ += 1
  print("\n{}/{}  values interpolated in accZ array".format(cpt_ind_interp_accZ , len(Lindex_to_interp_accZ)))
  
  # INTERPOLATING ABS acceleration array
  print("\n # # Intepolating abs acc data")
  Linterpolated_abs_acc = np.interp(Lindex_to_interp_abs_acc , LTime , Labs_acc)
  cpt_ind_interp_abs_acc = 0
  for i in range(0 , len(Labs_acc) , 1) :
    if Labs_acc[i] == np.nan : 
      Labs_acc[i] = Linterpolated_abs_acc[cpt_ind_interp_abs_acc]
      cpt_ind_interp_abs_acc += 1
  print("\n{}/{} values interpolated in abs_acc array".format(cpt_ind_interp_abs_acc , len(Lindex_to_interp_abs_acc)))

## Filtrage des valeurs non-converties en valeurs 'NaN'
print("\n\n # Filtering acceleration data (|Values| > {} -> np.nan values)".format(limit_value_to_error_ACCELERATION))
# LTime   [abs(LTime)    > abs(limit_value_to_error_TIME)] = np.nan # NO FILTERING FOR TIME LIST
Lacc_x  [abs(Lacc_x)   > abs(limit_value_to_error_ACCELERATION)] = np.nan
Lacc_y  [abs(Lacc_y)   > abs(limit_value_to_error_ACCELERATION)] = np.nan
Lacc_z  [abs(Lacc_z)   > abs(limit_value_to_error_ACCELERATION)] = np.nan
Labs_acc[abs(Labs_acc) > abs(limit_value_to_error_ACCELERATION)] = np.nan

print("\nNAN Values in LTime    = {} = {}%"  .format(np.isnan(LTime).sum()    , 100 * np.isnan(LTime   ).sum()/LTime.size   ))
print(  "NAN Values in Lacc_x   = {} = {}%"  .format(np.isnan(Lacc_x).sum()   , 100 * np.isnan(Lacc_x  ).sum()/Lacc_x.size  ))
print(  "NAN Values in Lacc_y   = {} = {}%"  .format(np.isnan(Lacc_y).sum()   , 100 * np.isnan(Lacc_y  ).sum()/Lacc_y.size  ))
print(  "NAN Values in Lacc_z   = {} = {}%"  .format(np.isnan(Lacc_z).sum()   , 100 * np.isnan(Lacc_z  ).sum()/Lacc_z.size  ))
print(  "NAN Values in Labs_acc = {} = {}%\n".format(np.isnan(Labs_acc).sum() , 100 * np.isnan(Labs_acc).sum()/Labs_acc.size))

## Histogramme de répartition des écarts temporels entre 2 mesures
print("\n\n # Ploting temporal steps histogram")
LTimeStep     = [LTime[i+1] - LTime[i] for i in range(0,len(LTime)-1,1)]
minLTimeStep  = min(LTimeStep)
maxLTimeStep  = max(LTimeStep)
dt            = np.mean(LTimeStep)
std_dt        = np.std(LTimeStep)

# Label de l'histogramme et range choisie des écarts temporels
if critere_range_histogram == 'std' :
  Xmin_histogram  = dt - nsigma*std_dt
  Xmax_histogram  = dt + nsigma*std_dt
# elif 'PourcentageOfValuesHistogram' :
#   Xmin_histogram  = dt - nsigma*std_dt
#   Xmax_histogram  = dt + nsigma*std_dt
#   LTimeStep_sort  = np.sort(LTimeStep)
  
#   previous_gap = abs(max(LTimeStep_sort) - dt)
#   cpt_identical_TimeStep_Gaps = 1
#   for i in range(0 , len(LTimeStep_sort) , 1) :
#     if abs(LTimeStep_sort[i] - dt) < previous_gap :
#       ind_dt = i
#       previous_step = abs(LTimeStep_sort[i] - dt)
#       cpt_identical_TimeStep_Gaps = 1
#     elif abs(LTimeStep_sort[i] - dt) == previous_gap :
#       cpt_identical_TimeStep_Gaps += 1
#     else :
#       cpt_identical_TimeStep_Gaps = 1
#   if cpt_identical_TimeStep_Gaps > 1 :
#     ind_dt
#   for i in range(0 , int(len(LTimeStep_sort)/2) , 1) :
#     if 100 * len(LTimeStep_sort[])/len(LTimeStep) > pourcentageValuesHistogram :
#       Xmin_histogram  = dt - nsigma*std_dt
#       Xmax_histogram  = dt + nsigma*std_dt
#       break
else :
  Xmin_histogram  = LTimeStep[0]
  Xmax_histogram  = LTimeStep[-1]

histogram_label = "dt_mean = {}ms\nstd = {}us".format(round(1000*dt,2) , round(1000000*std_dt,2))

#############################
## PARTIE LISSAGE DES DONNEES
if BOOL_Smooth_Data == True :
  print('\n\n # Smoothing data')
  print("\n # # Smoothing Time data")
  LTime    = signal.savgol_filter(LTime    , chosen_window_size_for_smoothing_filter , chosen_order_for_smoothing_filter) # window size used for filtering ; order of fitted polynomial
  print("\n # # Smoothing acc X data")
  Lacc_x   = signal.savgol_filter(Lacc_x   , chosen_window_size_for_smoothing_filter , chosen_order_for_smoothing_filter) # window size used for filtering ; order of fitted polynomial
  print("\n # # Smoothing acc Y data")
  Lacc_y   = signal.savgol_filter(Lacc_y   , chosen_window_size_for_smoothing_filter , chosen_order_for_smoothing_filter) # window size used for filtering ; order of fitted polynomial
  print("\n # # Smoothing acc Z data")
  Lacc_z   = signal.savgol_filter(Lacc_z   , chosen_window_size_for_smoothing_filter , chosen_order_for_smoothing_filter) # window size used for filtering ; order of fitted polynomial
  print("\n # # Smoothing abs acc data")
  Labs_acc = signal.savgol_filter(Labs_acc , chosen_window_size_for_smoothing_filter , chosen_order_for_smoothing_filter) # window size used for filtering ; order of fitted polynomial

#####################
## PARTIE INTEGRATION
# Intégration des accélérations et calcul des vitesses
# => Cumulative sum of trapezoids for integration, probably the best for a series of points :
vx = integrate.cumtrapz(Lacc_x  , LTime, initial = vx_0)
vy = integrate.cumtrapz(Lacc_y  , LTime, initial = vy_0)
vz = integrate.cumtrapz(Lacc_z  , LTime, initial = vz_0)
V  = integrate.cumtrapz(Labs_acc, LTime, initial = V_0)

# Intégration des vitesses et calcul des positions
# => Cumulative sum of trapezoids for integration, probably the best for a series of points :
X  = integrate.cumtrapz(vx  , LTime, initial = x0)
Y  = integrate.cumtrapz(vy  , LTime, initial = y0)
Z  = integrate.cumtrapz(vz  , LTime, initial = z0)
OM = integrate.cumtrapz(V   , LTime, initial = OM_0)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Plotting the results  # # # # # # # # # # # # # # # # # # # # #
print('\n\n # Plotting the results')

#######
## TIME
print("\n # # Plotting Time data")
fig1 = plt.figure(1)

# Graphe 1
plt.subplot(221)
plt.hist(LTimeStep, range = (Xmin_histogram,Xmax_histogram), bins = number_of_bins , color = color_fig1_1 , edgecolor = 'black' , label = histogram_label)
plt.xlabel('écart temporel inter-mesures')
plt.ylabel("Occurrences\n{}".format(histogram_label))
plt.title('Histogramme de répartition des écarts temporels entre 2 mesures')

# Graphe 2
## FAST FOURIER TRANSFORM - TIME
FFT_LTimeStep         = fft(LTimeStep , Ne)
Abs_FFT_LTimeStep     = np.abs(FFT_LTimeStep)
ind_Amp_max           = Abs_FFT_LTimeStep.argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
Amp_max               = Abs_FFT_LTimeStep[ind_Amp_max]
ind_AmpMax_without_f0 = Abs_FFT_LTimeStep[1 : len(Abs_FFT_LTimeStep)].argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
AmpMax_without_f0     = Abs_FFT_LTimeStep[ind_AmpMax_without_f0]
AmpMean_FFT_to_show   = np.mean(Abs_FFT_LTimeStep)
N                     = len(Abs_FFT_LTimeStep)
T                     = N/Ne
freq                  = np.arange(N)/T # crée la liste des fréquences auxquelles on va identifier un sinus
print("Pas fréquentiel FFT = {}mHz ; Nombre d'échantillons paramétré 'Ne' = {}".format(int((freq[1]-freq[0])*1000) , Ne))
plt.subplot(222)
plt.stem(freq , Abs_FFT_LTimeStep , 'b', markerfmt = " " , basefmt = "-b")
if   MODE_FFT_Frequences == "auto"    : plt.xlim(np.min(freq) , np.max(freq)/2) # on divise par 2 pour ne regarder que la première moitié du spectre (car symétrie verticale)
elif MODE_FFT_Frequences == "Min_Max" : plt.xlim(freq_min_FFT_to_show , freq_max_FFT_to_show)
if   MODE_FFT_Amplitudes == "Min_Max" : plt.ylim(AmpMinParameter_FFT_to_show , AmpMaxParameter_FFT_to_show)
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |FFT_LTimeStep(freq)|')

# Graphe 3
plt.subplot(223)
# plt.title("Titre")
plt.plot( np.linspace(1,len(LTime),len(LTime)) , LTime , color = color_fig1_1, label = 'Time(-)')
plt.xlabel("Arbitrary units")
plt.ylabel(Title_Time)
plt.legend()
plt.grid(True)

plt.show()

################
## ACCELERATIONS

#################
## ACCELERATION X
print("\n # # Plotting X acceleration data")
fig2 = plt.figure(2)

# Graphe 1
plt.subplot(211)
# plt.title("Titre")
plt.plot(LTime , Lacc_x , color = color_fig1_1, label = 'Ax(t)')
plt.xlabel(Title_Time)
plt.ylabel(Title_acc_x)
plt.legend()
plt.grid(True)

# Graphe 2
## FAST FOURIER TRANSFORM - X ACCELERATION
FFT_LaccX             = fft(Lacc_x , Ne)
Abs_FFT_LaccX         = np.abs(FFT_LaccX)
ind_Amp_max           = Abs_FFT_LaccX.argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
Amp_max               = Abs_FFT_LaccX[ind_Amp_max]
ind_AmpMax_without_f0 = Abs_FFT_LaccX[1 : len(Abs_FFT_LaccX)].argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
AmpMax_without_f0     = Abs_FFT_LaccX[ind_AmpMax_without_f0]
AmpMean_FFT_to_show    = np.mean(Abs_FFT_LaccX)
N                     = len(Abs_FFT_LaccX)
 
T                     = N/Ne
freq                  = np.arange(N)/T # crée la liste des fréquences auxquelles on va identifier un sinus
print("Pas fréquentiel FFT = {}mHz ; Nombre d'échantillons paramétré 'Ne' = {}".format(int((freq[1]-freq[0])*1000) , Ne))
plt.subplot(212)
plt.stem(freq , Abs_FFT_LaccX , 'b', markerfmt = " " , basefmt = "-b")
if   MODE_FFT_Frequences == "auto"    : plt.xlim(np.min(freq) , np.max(freq)/2) # on divise par 2 pour ne regarder que la première moitié du spectre (car symétrie verticale)
elif MODE_FFT_Frequences == "Min_Max" : plt.xlim(freq_min_FFT_to_show , freq_max_FFT_to_show)
if   MODE_FFT_Amplitudes == "Min_Max" : plt.ylim(AmpMinParameter_FFT_to_show , AmpMaxParameter_FFT_to_show)
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |FFT_LaccX(freq)|')

plt.show()

#################
## ACCELERATION Y
print("\n # # Plotting Y acceleration data")
fig3 = plt.figure(3)

# Graphe 1
plt.subplot(211)
# plt.title("Titre")
plt.plot(LTime , Lacc_y , color = color_fig1_2, label = 'Ay(t)')
plt.xlabel(Title_Time)
plt.ylabel(Title_acc_y)
plt.legend()
plt.grid(True)

# Graphe 2
## FAST FOURIER TRANSFORM - Y ACCELERATION
FFT_LaccY             = fft(Lacc_y , Ne)
Abs_FFT_LaccY         = np.abs(FFT_LaccY)
ind_Amp_max           = Abs_FFT_LaccY.argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
Amp_max               = Abs_FFT_LaccY[ind_Amp_max]
ind_AmpMax_without_f0 = Abs_FFT_LaccY[1 : len(Abs_FFT_LaccY)].argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
AmpMax_without_f0     = Abs_FFT_LaccY[ind_AmpMax_without_f0]
AmpMean_FFT_to_show    = np.mean(Abs_FFT_LaccY)
N                     = len(Abs_FFT_LaccY)
 
T                     = N/Ne
freq                  = np.arange(N)/T # crée la liste des fréquences auxquelles on va identifier un sinus
print("Pas fréquentiel FFT = {}mHz ; Nombre d'échantillons paramétré 'Ne' = {}".format(int((freq[1]-freq[0])*1000) , Ne))
plt.subplot(212)
plt.stem(freq , Abs_FFT_LaccY , 'b', markerfmt = " " , basefmt = "-b")
if   MODE_FFT_Frequences == "auto"    : plt.xlim(np.min(freq) , np.max(freq)/2) # on divise par 2 pour ne regarder que la première moitié du spectre (car symétrie verticale)
elif MODE_FFT_Frequences == "Min_Max" : plt.xlim(freq_min_FFT_to_show , freq_max_FFT_to_show)
if   MODE_FFT_Amplitudes == "Min_Max" : plt.ylim(AmpMinParameter_FFT_to_show , AmpMaxParameter_FFT_to_show)
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |FFT_LaccY(freq)|')

plt.show()

#################
## ACCELERATION Z
print("\n # # Plotting Z acceleration data")
fig4 = plt.figure(4)

# Graphe 1
plt.subplot(211)
# plt.title("Titre")
plt.plot(LTime , Lacc_z , color = color_fig1_3, label = 'Az(t)')
plt.xlabel(Title_Time)
plt.ylabel(Title_acc_z)
plt.legend()
plt.grid(True)

# Graphe 2
## FAST FOURIER TRANSFORM - Z ACCELERATION
FFT_LaccZ             = fft(Lacc_z , Ne)
Abs_FFT_LaccZ         = np.abs(FFT_LaccZ)
ind_Amp_max           = Abs_FFT_LaccZ.argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
Amp_max               = Abs_FFT_LaccZ[ind_Amp_max]
ind_AmpMax_without_f0 = Abs_FFT_LaccZ[1 : len(Abs_FFT_LaccZ)].argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
AmpMax_without_f0     = Abs_FFT_LaccZ[ind_AmpMax_without_f0]
AmpMean_FFT_to_show    = np.mean(Abs_FFT_LaccZ)
N                     = len(Abs_FFT_LaccZ)
 
T                     = N/Ne
freq                  = np.arange(N)/T # crée la liste des fréquences auxquelles on va identifier un sinus
print("Pas fréquentiel FFT = {}mHz ; Nombre d'échantillons paramétré 'Ne' = {}".format(int((freq[1]-freq[0])*1000) , Ne))
plt.subplot(212)
plt.stem(freq , Abs_FFT_LaccZ , 'b', markerfmt = " " , basefmt = "-b")
if   MODE_FFT_Frequences == "auto"    : plt.xlim(np.min(freq) , np.max(freq)/2) # on divise par 2 pour ne regarder que la première moitié du spectre (car symétrie verticale)
elif MODE_FFT_Frequences == "Min_Max" : plt.xlim(freq_min_FFT_to_show , freq_max_FFT_to_show)
if   MODE_FFT_Amplitudes == "Min_Max" : plt.ylim(AmpMinParameter_FFT_to_show , AmpMaxParameter_FFT_to_show)
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |FFT_LaccZ(freq)|')

plt.show()

#######################
## ACCELERATION ABSOLUE
print("\n # # Plotting absolute acceleration data")
fig5 = plt.figure(5)

# Graphe 1
plt.subplot(211)
# plt.title("Titre")
plt.plot(LTime , Labs_acc , color = color_fig1_4 , label = '||A||(t)')
plt.xlabel(Title_Time)
plt.ylabel(Title_abs_acc)
plt.legend()
plt.grid(True)

# Graphe 2
## FAST FOURIER TRANSFORM - ABSOLUTE ACCELERATION
FFT_Labs_acc          = fft(Labs_acc , Ne)
Abs_FFT_Labs_acc      = np.abs(FFT_LaccZ)
ind_Amp_max           = Abs_FFT_Labs_acc.argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
Amp_max               = Abs_FFT_Labs_acc[ind_Amp_max]
ind_AmpMax_without_f0 = Abs_FFT_Labs_acc[1 : len(Abs_FFT_Labs_acc)].argsort()[-1] # récupère le dernier indice (correspond par défaut au plus grand nombre) du tableau d'indices trié
AmpMax_without_f0     = Abs_FFT_Labs_acc[ind_AmpMax_without_f0]
AmpMean_FFT_to_show    = np.mean(Abs_FFT_Labs_acc)
N                     = len(Abs_FFT_Labs_acc)
 
T                     = N/Ne
freq                  = np.arange(N)/T # crée la liste des fréquences auxquelles on va identifier un sinus
print("Pas fréquentiel FFT = {}mHz ; Nombre d'échantillons paramétré 'Ne' = {}".format(int((freq[1]-freq[0])*1000) , Ne))
plt.subplot(212)
plt.stem(freq , Abs_FFT_Labs_acc , 'b', markerfmt = " " , basefmt = "-b")
if   MODE_FFT_Frequences == "auto"    : plt.xlim(np.min(freq) , np.max(freq)/2) # on divise par 2 pour ne regarder que la première moitié du spectre (car symétrie verticale)
elif MODE_FFT_Frequences == "Min_Max" : plt.xlim(freq_min_FFT_to_show , freq_max_FFT_to_show)
if   MODE_FFT_Amplitudes == "Min_Max" : plt.ylim(AmpMinParameter_FFT_to_show , AmpMaxParameter_FFT_to_show)
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |FFT_Labs_acc(freq)|')

plt.show()

######################
## VITESSES
print("\n # # Plotting velocity data")
fig6 = plt.figure(6)

#graphe 1
plt.subplot(221)
# plt.title("Titre")
plt.plot(LTime , vx , color = color_fig1_1, label = 'Vx(t)')
plt.xlabel(Title_Time)
plt.ylabel("Velocity X axis")
plt.legend()
plt.grid(True)

#graphe 2
plt.subplot(222)
# plt.title("Titre")
plt.plot(LTime , vy , color = color_fig1_2, label = 'Vy(t)')
plt.xlabel(Title_Time)
plt.ylabel("Velocity Y axis")
plt.legend()
plt.grid(True)

#graphe 3
plt.subplot(223)
# plt.title("Titre")
plt.plot(LTime , vz , color = color_fig1_3, label = 'Vz(t)')
plt.xlabel(Title_Time)
plt.ylabel("Velocity Z axis")
plt.legend()
plt.grid(True)

#graphe 4
plt.subplot(224)
# plt.title("Titre")
plt.plot(LTime , V , color = color_fig1_4 , label = 'V(t)')
plt.xlabel(Title_Time)
plt.ylabel("Velocity (absolute)")
plt.legend()
plt.grid(True)

plt.show()

############
## POSITIONS
print("\n # # Plotting position data")
fig7 = plt.figure(7)

#graphe 1
plt.subplot(221)
# plt.title("Titre")
plt.plot(LTime , X , color = color_fig1_1, label = 'X(t)')
plt.xlabel(Title_Time)
plt.ylabel("Position X axis")
plt.legend()
plt.grid(True)

#graphe 2
plt.subplot(222)
# plt.title("Titre")
plt.plot(LTime , Y , color = color_fig1_2, label = 'Y(t)')
plt.xlabel(Title_Time)
plt.ylabel("Position Y axis")
plt.legend()
plt.grid(True)

#graphe 3
plt.subplot(223)
# plt.title("Titre")
plt.plot(LTime , Z , color = color_fig1_3, label = 'Z(t)')
plt.xlabel(Title_Time)
plt.ylabel("Position Z axis")
plt.legend()
plt.grid(True)

#graphe 4
plt.subplot(224)
# plt.title("Titre")
plt.plot(LTime , OM , color = color_fig1_4 , label = 'OM(t)')
plt.xlabel(Title_Time)
plt.ylabel("Position (absolute)")
plt.legend()
plt.grid(True)

plt.show()

#################
## TRAJECTOIRE 3D
print("\n # # Plotting trajectory data")
fig8 = plt.figure(8)
ax8 = plt.axes(projection = '3d')

if True == True :
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
          print("\nWARNING :\n\tWhen computing the colors for plots: doublon at index {}".format(i))
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

# Data for a three-dimensional line
ax8.plot3D(X, Y, Z, 'gray')

# Data for three-dimensional scattered points
color_map_to_use = (mpl.colors.ListedColormap(Lcolors).with_extremes(over=Lcolors[0], under=Lcolors[len(Lcolors)-1]))
ax8.scatter3D(X , Y , Z , c = Z , cmap = color_map_to_use);
plt.show()

#################################
##  AFFICHAGE DE TOUS LES GRAPHES
plt.show()
# plt.savefig('fig5.png')

########################
# Showing the time spent
T_END_PROGRAM = time.time()
print("\nTotal time for computations = {}s\n".format(round(T_END_PROGRAM - T_START_PROGRAM , 2)))