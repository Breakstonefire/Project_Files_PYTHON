# -*- coding: utf-8 -*-
"""
Created on Sun May 19 09:10:03 2024
@author: lilia
"""

# Lilian LAGARRIGUE
# This code is meant to read all the ZetaPotential data (txt files) placed in a parametrized directory and plot the mean curve of one variable in the data.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # CLEANING / CLEARING / REMOVING / DELETING EVERYTHING BEFORE EXECUTION # #

g = globals() # Inititalizing g with globals()
d = dir() # Initializing d with dir() - this stores a list of all the variables in the program
for obj in d : # Checking user-defined variables in the directory # Checking for built-in variables/functions
  if not obj.startswith('__') : del globals()[obj] # Deleting the said obj, since a user-defined function

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing Libraries # # # # # # # # # # # # # # # # # # # # # #
from   datetime import datetime
import math
import matplotlib.pyplot as plt
import numpy             as np
import os
import pandas as pd
import time

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # REMOVING PLOTS BEFORE EXECUTION # # # # # # # # # # # # # # # # #
plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Functions # # # # # # # # # # # # # # # # # # # # # # # # # # #

def convert_data_from_str_to_float(x) :
  
  BOOL_Exposant = False
  for i in range(len(x)) :
    if    x[i] == 'E' and BOOL_Exposant == False  : BOOL_Exposant = True
    elif  x[i] == 'E' and BOOL_Exposant == True   :
      print('\nERROR :\n\tImpossible conversion, two or more exponents ...')
      return None
  if    BOOL_Exposant == False  : return float(x)
  elif  BOOL_Exposant == True   :
    for i in range(len(x)) :
      if x[i] == 'E' or x[i] == 'e' : return float(x[0:i]) * pow(10 , float(x[i+1 : len(x)]))

def ReturnRGBColorCode(color_string = 'red') :
  if    color_string == 'white'   or color_string == 'blanc'  : return [1   , 1   , 1]
  elif  color_string == 'gray'    or color_string == 'gris'   : return [0.5 , 0.5 , 0.5]
  elif  color_string == 'black'   or color_string == 'noir'   : return [0   , 0   , 0]
  elif  color_string == 'red'     or color_string == 'rouge'  : return [1   , 0   , 0]
  elif  color_string == 'orange'                              : return [1   , 0.5 , 0]
  elif  color_string == 'yellow'  or color_string == 'jaune'  : return [1   , 1   , 0]
  elif  color_string == 'apple'   or color_string == 'pomme'  : return [0.5 , 1   , 0]
  elif  color_string == 'green'   or color_string == 'vert'   : return [0   , 1   , 0]
  elif  color_string == 'cyan'                                : return [0   , 1   , 1]
  elif  color_string == 'blue'    or color_string == 'bleu'   : return [0   , 0   , 1]
  elif  color_string == 'indigo'                              : return [0.5 , 0   , 1]
  elif  color_string == 'purple'  or color_string == 'violet' : return [1   , 0   , 1]
  elif  color_string == 'magenta'                             : return [1   , 0   , 0.5]
  else                                                        :
    print("WARNING in return 'ReturnRGBColorCode' function :\n {} is not recognized, Black is returned by default as [0 , 0 , 0] RGB color code.".format(color_string))
    return [0 , 0 , 0]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Main () # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Saving the starting time
TIME_START_PROGRAM = time.time()

#############
## PARAMETERS
print("\n\n## SETTING MAIN PARAMETERS")

# Getting the current date and time
DATE_AND_TIME = datetime.now().strftime("%Y_%m_%d_%H_%M_%S") # Format Year_Month_Day_Hour_Minute_Second

# Getting the current directory
CurrentDirectory = os.getcwd()

# Setting the path leading to all data files
path_to_all_files = CurrentDirectory
# path_to_all_files = r'C:\Users\lilia\Documents\11_Liou_Lili_LeRoux_LaGags\Projet Carac Zeta Potential' # another path whether the user changes data location or not

# Changing to the path leading to all files (the current directory by default)
os.chdir(path_to_all_files)

# Defining the expected file format
file_format = 'txt'

# Defining the 'background data' keyword
Background_KEYWORD = 'Background'

# Defining the expected number of files (without the background file)
Number_of_Expected_Files_without_Background = 3

# Defining whether to mute all the very detailed messages in the command window or not
MUTE_MESSAGE_DETAILED = True # True # False

# Defining the reading pourcentage step between two messages in the command window
PourcentageStep = 10

# Defining the accepted caracters for 'pure data' to identify the header and the real data numbers
AcceptedCaracters_for_PureData = ",.+-E 0123456789\t\n\r" # The space and comma and dot and E letter might be present in the lines containing pure data

# Defining the HEADER delimiter caracter to slice string caracters and correcty extract HEADER variables and their names
HeaderDelimiter = ':' # Double-dot caracter

# Defining the PURE DATA delimiter caracters' list to slice string caracters and correcty extract PURE DATA numbers
PureDataDelimiter1 = ' '  # Space caracter
PureDataDelimiter2 = '\t' # Tabulation
PureDataDelimiter3 = '\n' # Back to the next line
PureDataDelimiter4 = '\r' # Here it is a space caracter
All_PureData_Delimiters_STRING = PureDataDelimiter1 + PureDataDelimiter2 + PureDataDelimiter3 + PureDataDelimiter4

# Defining the number of columns expected in the 'Pure Data' lines <=> number of variables in the 'Pure Data'
columns_predicted = 5

# Defining the error value to put in lists when errors occur
error_value_default = np.nan # Nan instead of a random value to keep the whole data as it really is

## Plots parameters
# FIGURE parameters
DefaultFigureWidth    = 5
DefaultFigureHeight   = 3
DefaultDPIResolution  = 600 # 800 is enough for detailed plots on a PC screen
DefaultFaceColor      = ReturnRGBColorCode('white') # Spyder Black Background color is [25/100 , 35/100 , 45/100]
DefaultEdgeColor      = [0.1 , 0.1 , 0.1] # Spyder Black Background color is [25/100 , 35/100 , 45/100]
DefaultLayout         = 'constrained' # 'constrained', 'compressed', 'tight'

# LINE / MARKER parameters
DefaultLineWidth  = 0.3
DefaultMarkerSize = 0.3
DefaultLineStyle  = '-'

# GRID parameters
DefaultGridValue      = True # True # False
DefaultGridColor      = ReturnRGBColorCode('gray')
DefaultGridLineStyle  = '--'
DefaultGridLineWidth  = 0.4
DefaultGridOpacity    = 0.25

# TEXT parameters
SMALL_SIZE  = 3
MEDIUM_SIZE = 5
BIGGER_SIZE = 7

# Setting the names of variables
Size_label            = "Size [nm]"
Number_label          = "Number [-]"
Concentration_label   = "Concentration [cm-3]"
Volume_label          = "Volume [nm3]"
Area_label            = "Area [nm²]"
Background_label      = "Background"

# Setting the range of X AXIS to show - setting the xmin and xmax to show in every plot
MODE_XminXmax_to_show = "relative"  # "relative" defines Xmin = 0 and Xmax as a given percentage of the absolute Xmax of each original data
                                    # "absolute" defines the Xmin and Xmax values to consider for every plots
if MODE_XminXmax_to_show == "relative" : # -> Parameters for a realtive Xmin Xmax setting
  Xmin_for_plots  = 0
  Percentage      = 10 # Xmax is defined as a given percentage of the absolute Xmax of each original data
elif MODE_XminXmax_to_show == "absolute" : # -> Parameters for an absolute Xmin Xmax setting
  Xmin_for_plots = 0
  Xmax_for_plots = 800

# Setting the defaults curve transparencies, from 0 to 1 (1 is plenty, 0 is fully transparent)
Main_Curve_Transparency       = 1
Std_Curve_Transparency        = 0.2
Background_Curve_Transparency = 0.1

# Setting EXCEL parameters in order to export data to excel format
ExcelName                 = "EXCEL_ZetaPotential"
Excel_Sheet_Default_Name  = "Excel_Sheet_1_from_python"
Excel_Cells_to_Freeze     = "A2"

#################################################################################
#### THE FOLLOWING BLOCK CAN BE COPY-PASTE IN EVERY CODE TO PLOT K colored curves
# K = value corresponding to the number of 'y(x)' curves to plot
# Theoretically the number of triplet-color-values in LColorsREDtoBLACK is equal to 5*(2^8-1) = 1275 # RED TO BLACK
LColorsREDtoBLACK , index = [[]] , 0
for i       in range(0  ,1276 ,1) : LColorsREDtoBLACK.append([0,0,0])
for green   in range(0  ,256,1)   : LColorsREDtoBLACK[index] , index = [255 , green ,0]     , index + 1 # Color is RED    # INCREASING THE VALUE of GREEN
len_col = len(LColorsREDtoBLACK)
for red     in range(254,-1 ,-1)  : LColorsREDtoBLACK[index] , index = [red , 255   ,0]     , index + 1 # Color is YELLOW # DECREASING THE VALUE of RED
len_col = len(LColorsREDtoBLACK)
for blue    in range(1  ,256,1)   : LColorsREDtoBLACK[index] , index = [0   , 255   ,blue]  , index + 1 # Color is GREEN  # INCREASING THE VALUE of BLUE
len_col = len(LColorsREDtoBLACK)
for green   in range(254,-1 ,-1)  : LColorsREDtoBLACK[index] , index = [0   , green ,255]   , index + 1 # Color is CYAN   # DECREASING THE VALUE of GREEN
len_col = len(LColorsREDtoBLACK)
for blue    in range(254,0  ,-1)  : LColorsREDtoBLACK[index] , index = [0   , 0     ,blue]  , index + 1 # Color is BLUE   # DECREASING THE VALUE of BLACK
# if LColorsREDtoBLACK[len(LColorsREDtoBLACK)-1] == LColorsREDtoBLACK[len(LColorsREDtoBLACK)-2] : del LColorsREDtoBLACK[len(LColorsREDtoBLACK)-1] # Color is BLACK  # DELETING THE DOUBLE VALUES IN THE END
len_col = len(LColorsREDtoBLACK)
for i       in range(len(LColorsREDtoBLACK)) : # Then we normalize the values of LColorsREDtoBLACK from [0:255] to [0:1]
  for j     in range(len(LColorsREDtoBLACK[0])) : LColorsREDtoBLACK[i][j] /= 255
LColorsBLACKtoRED = list(LColorsREDtoBLACK)
LColorsREDtoBLUE  = LColorsREDtoBLACK[0 : 4*(2**8-1) + 1]
LColorsBLACKtoRED.reverse()
#### END OF THE GENERIC BLOCK ####
##################################

# Setting the colors for every curves
DefaultColor_Size             = ReturnRGBColorCode("yellow")
DefaultColor_SizeSTD          = ReturnRGBColorCode("orange")
DefaultColor_Number           = ReturnRGBColorCode("magenta")
DefaultColor_NumberSTD        = ReturnRGBColorCode("red")
DefaultColor_Concentration    = ReturnRGBColorCode("purple")
DefaultColor_ConcentrationSTD = ReturnRGBColorCode("magenta")
DefaultColor_Volume           = ReturnRGBColorCode("blue")
DefaultColor_VolumeSTD        = ReturnRGBColorCode("cyan")
DefaultColor_Area             = ReturnRGBColorCode("cyan")
DefaultColor_AreaSTD          = ReturnRGBColorCode("blue")
DefaultColor_Background       = ReturnRGBColorCode("gray")

#####################################################################################
## INTERMEDIARY / NECESSARY / AUTOMATIC COMPUTATIONS BEFORE THE MAIN COMPUTATION PART
# Setting the different plot parameters before all plots
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
print("\n## STARTING MAIN COMPUTATIONS\n")
# Saving the starting time
T_START_PROGRAM = time.time()

# Searching all files having the right format (here it is '.txt')
List_of_Files_and_Directories   = os.listdir(path_to_all_files)
List_of_Files_with_Right_Format = [files for files in List_of_Files_and_Directories if os.path.isfile(os.path.join(path_to_all_files,files)) and files.endswith(file_format)]
Number_of_RightFormatFiles      = len(List_of_Files_with_Right_Format)

# Printing all the files-to-the-right-format found in the given path
print("\nListing all files\n")
for i in range(0 , Number_of_RightFormatFiles , 1) :
  FileName = List_of_Files_with_Right_Format[i]
  print("\t'{}' filename n°{}/{} : '{}'".format(file_format , i+1 , Number_of_RightFormatFiles , FileName))

# Opening all the files one after another
print("\nOpening all files\n")
BOOL_BackgroundFileFound      = False
count_number_of_files_treated = 0 # This variable counts the number of files treated without the background file
List_ALL_LINES_OF_ALL_FILES   = [[''] for u in range(Number_of_Expected_Files_without_Background + 1)] # Initializing the list that will contain every lines of every files AND THE BACKGROUND AT THE FIRST INDEX OF THE LIST

for i in range(0 , Number_of_RightFormatFiles , 1) :
  
  FileName = List_of_Files_with_Right_Format[i]
  print("\n\tOpening file n°{}/{} : '{}'".format(i+1 , Number_of_RightFormatFiles , FileName))
  TXT_FILE_tmp                  = open(                   os.path.join(path_to_all_files , FileName) , 'r')
  List_of_Lines_of_Current_File = TXT_FILE_tmp.readlines()
  
  # Checking whether it is a background file or not
  FileNameUpper = FileName.upper() # To be non-case sensitive to high letters
  if BOOL_BackgroundFileFound == False and FileNameUpper.find(Background_KEYWORD.upper()) > -1 :
    
    print("\t\t-> This file is taken as BACKGROUND FILE")
    BOOL_BackgroundFileFound == True
    List_ALL_LINES_OF_ALL_FILES[0] = List_of_Lines_of_Current_File
    
  elif count_number_of_files_treated < Number_of_Expected_Files_without_Background :
    
    print("\t\t-> This file is NOT taken as BACKGROUND FILE")
    count_number_of_files_treated += 1
    List_ALL_LINES_OF_ALL_FILES[count_number_of_files_treated] = List_of_Lines_of_Current_File

# INITIALIZING THE LISTS THAT WILL CONTAIN EVERY ARRAY OF DATA TAKEN IN EVERY FILES
Number_of_Measurements_MAX        = max([len(List_ALL_LINES_OF_ALL_FILES[u])  for u in range(0 , len(List_ALL_LINES_OF_ALL_FILES) , 1)])
List_of_ALL_Header_Variable_Names = [np.ones(Number_of_Measurements_MAX)      for u in range(Number_of_RightFormatFiles)] # Lists for HEADER DATA ONLY
List_of_ALL_Header_Variable_Value = [np.ones(Number_of_Measurements_MAX)      for u in range(Number_of_RightFormatFiles)] # Lists for HEADER DATA ONLY
List_of_ALL_Size                  = [np.ones(Number_of_Measurements_MAX)      for u in range(Number_of_RightFormatFiles)] # Lists for PURE DATA ONLY
List_of_ALL_Number                = [np.ones(Number_of_Measurements_MAX)      for u in range(Number_of_RightFormatFiles)] # Lists for PURE DATA ONLY
List_of_ALL_Concentration         = [np.ones(Number_of_Measurements_MAX)      for u in range(Number_of_RightFormatFiles)] # Lists for PURE DATA ONLY
List_of_ALL_Volume                = [np.ones(Number_of_Measurements_MAX)      for u in range(Number_of_RightFormatFiles)] # Lists for PURE DATA ONLY
List_of_ALL_Area                  = [np.ones(Number_of_Measurements_MAX)      for u in range(Number_of_RightFormatFiles)] # Lists for PURE DATA ONLY

# Reading the lines and skipping the header of each and every file
for num_file_tmp in range(0 , Number_of_RightFormatFiles , 1) :
  
  FileName            = List_of_Files_with_Right_Format[num_file_tmp]
  print("\n\tReading file n°{}/{} : '{}'".format(num_file_tmp+1 , Number_of_RightFormatFiles , FileName))
  CurrentFile         = open(FileName , 'r' , newline = '').readlines()
  Number_of_lines_MAX = len(CurrentFile)
  print("\n\t\t{} lines in this file\n".format(Number_of_lines_MAX))
  
  # Defining a list of Lines for the Header and one for the 'Pure Data'
  List_of_HeaderLineNumbers   = []
  List_of_PureDataLineNumbers = []
  
  # Reading the file line by line to identify wich lines are part of the header and which lines are part of 'Pure Data'
  PreviousReadingPourcentage = 0
  
  for num_line in range(0 , Number_of_lines_MAX , 1) :
    
    current_line = CurrentFile[num_line]
    
    # Showing reading pourcentage
    ReadingPourcentage = int(max(0 , 100 * num_line / Number_of_lines_MAX))
    if int(PreviousReadingPourcentage) + PourcentageStep <= int(ReadingPourcentage) :
      print("\t\tFILE READING : Line {}/{} = {}%".format(num_line , Number_of_lines_MAX , ReadingPourcentage))
      PreviousReadingPourcentage = int(ReadingPourcentage)
    
    # Reading every caracter to detect the header boundaries first
    BOOLEAN_LineWithPureData        = True # One starts reading the line considering (by default) it will only contain 'Pure Data' (no header etc ...)
    BOOL_First_Time_Showing_Message = True # Boolean to indicate whether the following message has already been shown or not
    
    for num_carac in range(0 , len(current_line) , 1) :
      caracter = current_line[num_carac]
      AcceptedCaracters_for_PureData_UPPER = AcceptedCaracters_for_PureData.upper()
      if AcceptedCaracters_for_PureData_UPPER.find(caracter.upper()) < 0 :
        BOOLEAN_LineWithPureData = False
        if BOOL_First_Time_Showing_Message == True :
          BOOL_First_Time_Showing_Message = False
          if MUTE_MESSAGE_DETAILED == False : print("\t\t\tLINE {} : #{}th Caracter '{}' unaccepted for 'Pure Data' -> HEADER LINE ".format(num_line , num_carac , caracter))
          
    # Affecting the line number to the header or not
    if BOOLEAN_LineWithPureData == True : List_of_PureDataLineNumbers.append(num_line)
    else                                : List_of_HeaderLineNumbers  .append(num_line)
    
    # Showing which type of line it is
    if MUTE_MESSAGE_DETAILED == False :
      if BOOLEAN_LineWithPureData == True : print("\t\t\tLINE {} -> PURE DATA LINE".format(num_line))
      else                                : print("\t\t\tLINE {} -> HEADER DATA LINE".format(num_line))
  
  # Showing the number of lines detected for each line type
  Number_of_lines_HEADER_tmp = len(List_of_HeaderLineNumbers)
  Number_of_lines_DATA_tmp   = len(List_of_PureDataLineNumbers)
  print("\n\t\t{} lines for Header".format(Number_of_lines_HEADER_tmp))
  print("\t\t{} lines for Data".format(Number_of_lines_DATA_tmp))
  
  ## DATA EXTRACTION PART - HEADER + PURE DATA
  print("\n\t\tEXTRACTING DATA FROM HEADER + PURE DATA")
  
  # HEADER DATA EXTRACTION
  # Initializing the lists (to np.nan values) that will contain the converted string chains to real numbers
  LHeaderVariableName_tmp   = ['' for u in range(Number_of_lines_HEADER_tmp)]
  LHeaderVariableValue_tmp  = ['' for u in range(Number_of_lines_HEADER_tmp)]
  
  # Extracting and converting data from 'Header' lines
  print("\n\t\tExtracting and converting data from 'Header' lines")
  for num_line in range(0 , Number_of_lines_HEADER_tmp , 1) :
    current_line_number       = List_of_HeaderLineNumbers[num_line]
    current_line              = CurrentFile[current_line_number]
    ind_previous_delimiter    = -1 # Initializing the previous delimiter position in the string chain to -1
    BOOL_HeaderDelimiterFound = False
    
    # Reading every caracter to identify the different string chains and the delimiters between them
    for num_carac in range(0 , len(current_line) , 1) :
      
      caracter = current_line[num_carac].upper()
      if BOOL_HeaderDelimiterFound == False : # IF ANY HEADER DELIMITER HASN'T BEEN FOUND YET
        
        if caracter == HeaderDelimiter : # IF A HEADER DELIMITER HAS BEEN FOUND
          
          if ind_previous_delimiter + 1 < num_carac : # IF THERE IS AT LEAST ONE CARACTER BETWEEN THE PREVIOUS DELIMITER AND THE CURRENT ONE
            
            BOOL_HeaderDelimiterFound = True
            
            # Adding the header variable name first
            VariableName_str_tmp = current_line[ind_previous_delimiter + 1 : num_carac] # Extracting the string chain
            LHeaderVariableName_tmp[num_line] = VariableName_str_tmp # Listing every string chain
            
            # Adding the header variable value after
            Variables_str_tmp = current_line[num_carac + 1 : len(current_line)] # Extracting the string chain
            LHeaderVariableValue_tmp[num_line] = Variables_str_tmp # Listing every string chain
            
          ind_previous_delimiter = num_carac # Refreshing the value of the previous delimiter
          
        elif (num_carac == len(current_line) - 1) : # IF THE END OF THE LINE HAS BEEN REACHED
        
          # Adding the header variable name first
          VariableName_str_tmp = current_line[ind_previous_delimiter + 1 : num_carac] # Extracting the string chain
          LHeaderVariableName_tmp[num_line] = VariableName_str_tmp # Listing every string chain
          
          # Adding the header variable value after
          Variables_str_tmp = current_line[num_carac + 1 : len(current_line)] # Extracting the string chain
          LHeaderVariableValue_tmp[num_line] = Variables_str_tmp # Listing every string chain
  
  # Listing every array of header data in the big lists that will contain every array of header data of every file
  List_of_ALL_Header_Variable_Names[num_file_tmp] = LHeaderVariableName_tmp
  List_of_ALL_Header_Variable_Value[num_file_tmp] = LHeaderVariableValue_tmp
          
  # PURE DATA EXTRACTION  
  # Initializing the lists (to np.nan values) that will contain the converted string chains to real numbers
  LSize_tmp           = np.nan * np.ones(Number_of_lines_DATA_tmp)
  LNumber_tmp         = np.nan * np.ones(Number_of_lines_DATA_tmp)
  LConcentration_tmp  = np.nan * np.ones(Number_of_lines_DATA_tmp)
  LVolume_tmp         = np.nan * np.ones(Number_of_lines_DATA_tmp)
  LArea_tmp           = np.nan * np.ones(Number_of_lines_DATA_tmp)
  
  # Extracting and converting data from 'Pure Data' lines
  print("\n\t\tExtracting and converting data from 'Pure Data' lines")
  for num_line in range(0 , Number_of_lines_DATA_tmp , 1) :
    
    current_line_number = List_of_PureDataLineNumbers[num_line]
    current_line        = CurrentFile[current_line_number]
    
    cpt_columns             = 0 # Initializing the column counter
    cpt_convert_success     = 0 # Initializing the successful conversion counter
    cpt_convert_error       = 0 # Initializing the unsuccessful conversion counter
    ind_previous_delimiter  = -1 # Initializing the previous delimiter position in the string chain to -1
    Lindexdelimiters        = [ind_previous_delimiter]  # Listing the first position of the very first delimiter (-1 index position is out of the string chain but allows to identify the next ones)
    LValues_str_tmp         = []  # Initializing the empty list of string chains to convert (the numbers one is looking to extract)
    
    # Reading every caracter to identify the different string chains and the delimiters between them
    for num_carac in range(0 , len(current_line) , 1) :
      
      caracter = current_line[num_carac].upper()
      if All_PureData_Delimiters_STRING.find(caracter) > -1 : # IF A DELIMITER HAS BEEN FOUND
        
        if ind_previous_delimiter + 1 < num_carac : # IF THERE IS AT LEAST ONE CARACTER BETWEEN THE PREVIOUS DELIMITER AND THE CURRENT ONE
        
          Values_str_tmp         = current_line[ind_previous_delimiter + 1 : num_carac] # Extracting the string chain corresponding to the value to convert in real number
          LValues_str_tmp.append(Values_str_tmp)    # Listing every string chain to convert later
          cpt_columns            = cpt_columns + 1  # Increasing the number of columns of values identified
          
        Lindexdelimiters.append(num_carac)        # Listing the position of every delimiter
        ind_previous_delimiter = num_carac        # Refreshing the value of the previous delimiter
        
      elif (num_carac == len(current_line) - 1) : # IF THE END OF THE LINE HAS BEEN REACHED
      
        Values_str_tmp = current_line[ind_previous_delimiter + 1 : len(current_line)] # Extracting the string chain corresponding to the value to convert in real number
        LValues_str_tmp.append(Values_str_tmp)  # Listing every string chain to convert later
        cpt_columns    = cpt_columns + 1        # Increasing the number of columns of values identified
    
    if cpt_columns == columns_predicted : # IF THE NUMBER OF COLUMNS (so the number of values to convert) IS EXPECTED
      
      # Putting all expected strings in variables
      Size , Number , Concentration , Volume , Area = LValues_str_tmp[:]
      
      # TRYING TO CONVERT EVERY STRING CHAIN TO REAL NUMBERS
      try :
        Size                = convert_data_from_str_to_float(Size)
        LSize_tmp[num_line] = Size
        cpt_convert_success += 1
      except :
        print('\nWARNING :\n\tImpossible conversion from string to float of Size n°{} : {}'.format(num_line , Size))
        cpt_convert_error += 1
      
      try :
        Number                = convert_data_from_str_to_float(Number)
        LNumber_tmp[num_line] = Number
        cpt_convert_success += 1
      except :
        print('\nWARNING :\n\tImpossible conversion from string to float of Number n°{} : {}'.format(num_line , Number))
        cpt_convert_error += 1
      
      try :
        Concentration             = convert_data_from_str_to_float(Concentration)
        LConcentration_tmp[num_line]  = Concentration
        cpt_convert_success       += 1
      except :
        print('\nWARNING :\n\tImpossible conversion from string to float of Concentration n°{} : {}'.format(num_line , Concentration))
        cpt_convert_error += 1
      
      try :
        Volume              = convert_data_from_str_to_float(Volume)
        LVolume_tmp[num_line]   = Volume
        cpt_convert_success += 1
      except :
        print('\nWARNING :\n\tImpossible conversion from string to float of Volume n°{} : {}'.format(num_line , Volume))
        cpt_convert_error += 1
      
      try :
        Area                = convert_data_from_str_to_float(Area)
        LArea_tmp[num_line]     = Area
        cpt_convert_success += 1
      except :
        print('\nWARNING :\n\tImpossible conversion from string to float of Area n°{} : {}'.format(num_line , Area))
        cpt_convert_error += 1
        
    else : # MORE OR LESS COLUMNS WERE FOUND => THERE MIGHT BE A PROBLEM IN THE CODE OR IN THE DATA
      print('\nWARNING :\n\tNumber of columns detected = {} =/= {}'.format(cpt_columns , columns_predicted))
    
    # WARNING MESSAGE IF ONE OR MORE ERRORS OCCURED DURING CONVERSIONS
    if cpt_convert_error > 0 :
      print('\nWARNING : (Reading file line N°{})'.format(num_line))
      print('\n\tConversion errors      = {}'     .format(cpt_convert_error))
      print('\n\tSuccessful conversions = {}'     .format(cpt_convert_success))
      
  # Listing every array of data in the big lists that will contain every array of data of every file
  List_of_ALL_Size         [num_file_tmp] = LSize_tmp
  List_of_ALL_Number       [num_file_tmp] = LNumber_tmp
  List_of_ALL_Concentration[num_file_tmp] = LConcentration_tmp
  List_of_ALL_Volume       [num_file_tmp] = LVolume_tmp
  List_of_ALL_Area         [num_file_tmp] = LArea_tmp

####################################################################
## CHECKING HEADER DATA VARIABLE NAMES TO BE THE SAMES IN EVERY FILE
print("\nChecking header data variable names to be the sames in every file")
print("\n\t# Checking the number of lines in every headers to be the same")
BOOL_Identical_Variable_Number_HEADER = True
expected_variable_number_header = len(List_of_ALL_Header_Variable_Names[0]) # The first header found is taken as the reference by default
for i in range(1 , len(List_of_ALL_Header_Variable_Names) , 1) :
  current_number_of_measurements = len(List_of_ALL_Header_Variable_Names[i])
  if current_number_of_measurements != expected_variable_number_header :
    BOOL_Identical_Variable_Number_HEADER = False
    print("\t\t{}th list hasn't the same number of variable as expected : {} (list {}) / {} (expect.)".format(i , current_number_of_measurements , i , expected_variable_number_header))

BOOL_Identical_Variable_Name_HEADER_tmp = True
if BOOL_Identical_Variable_Number_HEADER == True :
  print("\n\t# Checking every variable names for every line in the header to be indentical from a file to another")
  BOOL_Identical_Variable_Name_HEADER_tmp = True
  
  for num_variable_name in range(0 , len(List_of_ALL_Header_Variable_Names[0]) , 1) :
    expected_variable_name_header_tmp = List_of_ALL_Header_Variable_Names[0][num_variable_name] # The first header found is taken as the reference by default
    
    for num_file_tmp in range(1 , Number_of_RightFormatFiles , 1) :
      variable_name_header_tmp = List_of_ALL_Header_Variable_Names[num_file_tmp][num_variable_name] # The first header found is taken as the reference by default
      
      if variable_name_header_tmp != expected_variable_name_header_tmp :
        BOOL_Identical_Variable_Name_HEADER_tmp = False
        print("\t\t{}th variable name isn't the same as expected : (var. name {} // expected in first file) => '{}' / '{}'".format(num_variable_name , num_variable_name , variable_name_header_tmp , expected_variable_name_header_tmp))

else : BOOL_Identical_Variable_Name_HEADER_tmp = False

if BOOL_Identical_Variable_Number_HEADER == True and BOOL_Identical_Variable_Name_HEADER_tmp == True :
  print("\t\tThe header variable names are identicals, the first header variable names are then taken as the reference")
  List_Header_Variable_Names_for_EXCEL    = List_of_ALL_Header_Variable_Names[0]
  List_Header_Variable_Strings_for_EXCEL  = List_of_ALL_Header_Variable_Value

else :
  print("\t\tWARNING :")
  print("\t\t\tThe headers are different, the first header variable names are taken as the reference no matter the others")
  List_Header_Variable_Names_for_EXCEL    = List_of_ALL_Header_Variable_Names[0]
  List_Header_Variable_Strings_for_EXCEL  = List_of_ALL_Header_Variable_Value
  
# Deleting the eventual last line of the HEADER because it contains the main names of variables represented in the graphs without more information in the header
print("\t\tDeleting last line of the HEADER because it contains the main names of variables represented in the graphs without more information in the header")
List_Header_Variable_Names_for_EXCEL    = List_Header_Variable_Names_for_EXCEL  [0 : len(List_Header_Variable_Names_for_EXCEL)    - 1]
for num_file_tmp in range(0 , Number_of_RightFormatFiles , 1) : List_Header_Variable_Strings_for_EXCEL[num_file_tmp]  = List_Header_Variable_Strings_for_EXCEL[num_file_tmp][0 : len(List_Header_Variable_Strings_for_EXCEL[num_file_tmp])  - 1]

######################################################################
## CHECKING THE NUMBER OF PURE DATA LINES TO BE THE SAME IN EVERY FILE
print("\nChecking the amount of data to be the same amount in every file")

print("\n\t# Checking for SIZE DATA")
BOOL_Correct_number_of_measurements_SIZE  = True
expected_number_of_measurements_SIZE      = len(List_of_ALL_Size[0]) # The first list found is taken as the reference by default
for i in range(1 , len(List_of_ALL_Size) , 1) :
  current_number_of_measurements = len(List_of_ALL_Size[i])
  if current_number_of_measurements != expected_number_of_measurements_SIZE :
    BOOL_Correct_number_of_measurements_SIZE = False
    print("\t\t{}th list hasn't the same number of measurements as expected : {} (list {}) / {} (expect.)".format(i , current_number_of_measurements , i , expected_number_of_measurements_SIZE))

print("\n\t# Checking for NUMBER DATA")
BOOL_Correct_number_of_measurements_NUMBER  = True
expected_number_of_measurements_NUMBER      = len(List_of_ALL_Number[0]) # The first list found is taken as the reference by default
for i in range(1 , len(List_of_ALL_Number) , 1) :
  current_number_of_measurements = len(List_of_ALL_Number[i])
  if current_number_of_measurements != expected_number_of_measurements_NUMBER :
    BOOL_Correct_number_of_measurements_NUMBER = False
    print("\t\t{}th list hasn't the same number of measurements as expected : {} (list {}) / {} (expect.)".format(i , current_number_of_measurements , i , expected_number_of_measurements_NUMBER))

print("\n\t# Checking for CONCENTRATION DATA")
BOOL_Correct_number_of_measurements_CONCENTRATION = True
expected_number_of_measurements_CONCENTRATION = len(List_of_ALL_Concentration[0]) # The first list found is taken as the reference by default
for i in range(1 , len(List_of_ALL_Concentration) , 1) :
  current_number_of_measurements = len(List_of_ALL_Concentration[i])
  if current_number_of_measurements != expected_number_of_measurements_CONCENTRATION :
    BOOL_Correct_number_of_measurements_CONCENTRATION = False
    print("\t\t{}th list hasn't the same number of measurements as expected : {} (list {}) / {} (expect.)".format(i , current_number_of_measurements , i , expected_number_of_measurements_CONCENTRATION))

print("\n\t# Checking for VOLUME DATA")
BOOL_Correct_number_of_measurements_VOLUME = True
expected_number_of_measurements_VOLUME = len(List_of_ALL_Volume[0]) # The first list found is taken as the reference by default
for i in range(1 , len(List_of_ALL_Volume) , 1) :
  current_number_of_measurements = len(List_of_ALL_Volume[i])
  if current_number_of_measurements != expected_number_of_measurements_VOLUME :
    BOOL_Correct_number_of_measurements_VOLUME = False
    print("\t\t{}th list hasn't the same number of measurements as expected : {} (list {}) / {} (expect.)".format(i , current_number_of_measurements , i , expected_number_of_measurements_VOLUME))

print("\n\t# Checking for AREA DATA")
BOOL_Correct_number_of_measurements_AREA = True
expected_number_of_measurements_AREA = len(List_of_ALL_Area[0]) # The first list found is taken as the reference by default
for i in range(1 , len(List_of_ALL_Area) , 1) :
  current_number_of_measurements = len(List_of_ALL_Area[i])
  if current_number_of_measurements != expected_number_of_measurements_AREA :
    BOOL_Correct_number_of_measurements_AREA = False
    print("\t\t{}th list hasn't the same number of measurements as expected : {} (list {}) / {} (expect.)".format(i , current_number_of_measurements , i , expected_number_of_measurements_AREA))

###########################################################################################
## STOPPING THE PROGRAM WHETHER ONE OR ANOTHER FILE CONTAINS A DIFFERENT NUMBER OF ELEMENTS
if    BOOL_Correct_number_of_measurements_SIZE          == False : exit()
elif  BOOL_Correct_number_of_measurements_NUMBER        == False : exit()
elif  BOOL_Correct_number_of_measurements_CONCENTRATION == False : exit()
elif  BOOL_Correct_number_of_measurements_VOLUME        == False : exit()
elif  BOOL_Correct_number_of_measurements_AREA          == False : exit()

######################################################################
## SPLITTING ALL LISTS IN TWO DIFFERENT PARTS - Linear and Logarithmic
# The lists are all containing a first dataset which is composed of 'linear measurements' and a second one composed of 'logarithmic measurements'
# => Splitting in two and saving away the 'linear measurements' for now and the 'logarithmic measurements' for later

print("\nSplitting data in two different parts\n")
print("\n\tIdentifying the Linear and Logarithmic parts in these files")
print("\nChecking the amount of linear data to be the same in every file")

BOOL_IDENTICAL_LINE_TO_SPLIT_EVERY_DATASET = True
for i in range(Number_of_RightFormatFiles) :
  for j in range(expected_number_of_measurements_SIZE) :
    if List_of_ALL_Size[i][j] == -1 :
      Expected_splitting_dataset_line_SIZE = j
      break
  for j in range(expected_number_of_measurements_NUMBER) :
    if List_of_ALL_Number[i][j] == -1 :
      Expected_splitting_dataset_line_NUMBER = j
      break
  for j in range(expected_number_of_measurements_CONCENTRATION) :
    if List_of_ALL_Concentration[i][j] == -1 :
      Expected_splitting_dataset_line_CONCENTRATION = j
      break
  for j in range(expected_number_of_measurements_VOLUME) :
    if List_of_ALL_Volume[i][j] == -1 :
      Expected_splitting_dataset_line_VOLUME = j
      break
  # for j in range(expected_number_of_measurements_AREA) :
  #   if List_of_ALL_Area[i][j] == -1 :
  #     Expected_splitting_dataset_line_AREA = j
  #     break
  Expected_splitting_dataset_line_AREA = Expected_splitting_dataset_line_VOLUME # EXCEPTIONNALLY : for the AREA, a zero is expected at this specific splitting line but zeros are present elsewhere in the AREA dataset so one ignores to try to find out the splitting line thes same way one did for the other variables
  if Expected_splitting_dataset_line_SIZE               != Expected_splitting_dataset_line_NUMBER \
    or Expected_splitting_dataset_line_NUMBER           != Expected_splitting_dataset_line_CONCENTRATION \
      or Expected_splitting_dataset_line_CONCENTRATION  != Expected_splitting_dataset_line_VOLUME \
        or Expected_splitting_dataset_line_VOLUME       != Expected_splitting_dataset_line_AREA :
          BOOL_IDENTICAL_LINE_TO_SPLIT_EVERY_DATASET = False
          print("\n\tERROR :")
          print("\t\tDifferent lines observed to split every datasets in two parts :")
          print("\t\tLines {} , {} , {} , {} , {} detected (Size / Number / Concentration / Volume / Area)".format(Expected_splitting_dataset_line_SIZE , Expected_splitting_dataset_line_NUMBER , Expected_splitting_dataset_line_CONCENTRATION , Expected_splitting_dataset_line_VOLUME , Expected_splitting_dataset_line_AREA))
          print("\t\t=> IMPOSSIBLE TO SPLIT IDENTICALLY => STOPPING PROGRAM")
          break

# STOPPING THE PROGRAM IF NOT POSSIBLE TO SPLIT DATA INTO TWO PARTS
if BOOL_IDENTICAL_LINE_TO_SPLIT_EVERY_DATASET == False : exit()
print("\n\t... Identical amount of linear data checked successfully")
print("\n\t... Splitting every datasets in two parts")
List_of_ALL_Size_LIN          = [[List_of_ALL_Size[i][j]          for j in range(0                                                  , Expected_splitting_dataset_line_SIZE          , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Size_LOG          = [[List_of_ALL_Size[i][j]          for j in range(Expected_splitting_dataset_line_SIZE + 1           , len(List_of_ALL_Size[i])                      , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Number_LIN        = [[List_of_ALL_Number[i][j]        for j in range(0                                                  , Expected_splitting_dataset_line_NUMBER        , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Number_LOG        = [[List_of_ALL_Number[i][j]        for j in range(Expected_splitting_dataset_line_NUMBER + 1         , len(List_of_ALL_Number[i])                    , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Concentration_LIN = [[List_of_ALL_Concentration[i][j] for j in range(0                                                  , Expected_splitting_dataset_line_CONCENTRATION , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Concentration_LOG = [[List_of_ALL_Concentration[i][j] for j in range(Expected_splitting_dataset_line_CONCENTRATION + 1  , len(List_of_ALL_Concentration[i])             , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Volume_LIN        = [[List_of_ALL_Volume[i][j]        for j in range(0                                                  , Expected_splitting_dataset_line_VOLUME        , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Volume_LOG        = [[List_of_ALL_Volume[i][j]        for j in range(Expected_splitting_dataset_line_VOLUME + 1         , len(List_of_ALL_Volume[i])                    , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Area_LIN          = [[List_of_ALL_Area[i][j]          for j in range(0                                                  , Expected_splitting_dataset_line_AREA          , 1)] for i in range(Number_of_RightFormatFiles)]
List_of_ALL_Area_LOG          = [[List_of_ALL_Area[i][j]          for j in range(Expected_splitting_dataset_line_AREA + 1           , len(List_of_ALL_Area[i])                      , 1)] for i in range(Number_of_RightFormatFiles)]

###################################################
## COMPUTING DATA - SUBSTRACTING BACKGROUND TO DATA
print("\nComputing data")
print("\n\tSubstracting backgroung to data")
if BOOL_Correct_number_of_measurements_SIZE == True :
  LMean_Size                  = np.array([sum(     List_of_ALL_Size_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , Expected_splitting_dataset_line_SIZE)])
  LMean_Size_LOG              = np.array([sum(     List_of_ALL_Size_LOG[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , expected_number_of_measurements_SIZE - Expected_splitting_dataset_line_SIZE - 1)])
  LStd_Size                   = np.array([np.std([ List_of_ALL_Size_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1)])                                 for j in range(0 , Expected_splitting_dataset_line_SIZE)])
  LBackground_on_Size         = np.array(List_of_ALL_Size_LIN[:][0])
  LMean_Size_minus_Background = LMean_Size - LBackground_on_Size
else : print("\n\tERROR :\n\t\'SIZE' lists have not the same number of measurements => IMPOSSIBLE TO COMPUTE THE MEAN VALUES")

if BOOL_Correct_number_of_measurements_NUMBER == True :
  LMean_Number                  = np.array([sum(     List_of_ALL_Number_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , Expected_splitting_dataset_line_NUMBER)])
  LMean_Number_LOG              = np.array([sum(     List_of_ALL_Number_LOG[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , expected_number_of_measurements_NUMBER - Expected_splitting_dataset_line_NUMBER - 1)])
  LStd_Number                   = np.array([np.std([ List_of_ALL_Number_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1)])                                 for j in range(0 , Expected_splitting_dataset_line_NUMBER)])
  LBackground_on_Number         = np.array(          List_of_ALL_Number_LIN[:][0])
  LMean_Number_minus_Background = LMean_Number - LBackground_on_Number
else : print("\n\tERROR :\n\t\t'NUMBER' lists have not the same number of measurements => IMPOSSIBLE TO COMPUTE THE MEAN VALUES")

if BOOL_Correct_number_of_measurements_CONCENTRATION == True :
  LMean_Concentration                   = np.array([sum(     List_of_ALL_Concentration_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1)  for j in range(0 , Expected_splitting_dataset_line_CONCENTRATION)])
  LMean_Concentration_LOG               = np.array([sum(     List_of_ALL_Concentration_LOG[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1)  for j in range(0 , expected_number_of_measurements_CONCENTRATION - Expected_splitting_dataset_line_CONCENTRATION - 1)])
  LStd_Concentration                    = np.array([np.std([ List_of_ALL_Concentration_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1)])                                  for j in range(0 , Expected_splitting_dataset_line_CONCENTRATION)])
  LBackground_on_Concentration          = np.array(          List_of_ALL_Concentration_LIN[:][0])
  LMean_Concentration_minus_Background  = LMean_Concentration - LBackground_on_Concentration
else : print("\n\tERROR :\n\t\t'CONCENTRATION' lists have not the same number of measurements => IMPOSSIBLE TO COMPUTE THE MEAN VALUES")

if BOOL_Correct_number_of_measurements_VOLUME == True :
  LMean_Volume                  = np.array([sum(     List_of_ALL_Volume_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , Expected_splitting_dataset_line_VOLUME)])
  LMean_Volume_LOG              = np.array([sum(     List_of_ALL_Volume_LOG[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , expected_number_of_measurements_VOLUME - Expected_splitting_dataset_line_VOLUME - 1)])
  LStd_Volume                   = np.array([np.std([ List_of_ALL_Volume_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1)])                                 for j in range(0 , Expected_splitting_dataset_line_VOLUME)])
  LBackground_on_Volume         = np.array(          List_of_ALL_Volume_LIN[:][0])
  LMean_Volume_minus_Background = LMean_Volume - LBackground_on_Volume
else : print("\n\tERROR :\n\t\t'VOLUME' lists have not the same number of measurements => IMPOSSIBLE TO COMPUTE THE MEAN VALUES")

if BOOL_Correct_number_of_measurements_AREA == True :
  LMean_Area                  = np.array([sum(     List_of_ALL_Area_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , Expected_splitting_dataset_line_AREA)])
  LMean_Area_LOG              = np.array([sum(     List_of_ALL_Area_LOG[i][j] for i in range(1 , Number_of_RightFormatFiles , 1))/(Number_of_RightFormatFiles - 1) for j in range(0 , expected_number_of_measurements_AREA - Expected_splitting_dataset_line_AREA - 1)])
  LStd_Area                   = np.array([np.std([ List_of_ALL_Area_LIN[i][j] for i in range(1 , Number_of_RightFormatFiles , 1)])                                 for j in range(0 , Expected_splitting_dataset_line_AREA)])
  LBackground_on_Area         = np.array(          List_of_ALL_Area_LIN[:][0])
  LMean_Area_minus_Background = LMean_Area - LBackground_on_Area
else : print("\n\tERROR :\n\t\t'AREA' lists have not the same number of measurements => IMPOSSIBLE TO COMPUTE THE MEAN VALUES")

# Computing the correlation coefficient between Size and an arbitrary unit axis (shows the linearity of Size Data)
R2_Size = np.power(np.corrcoef(np.cumsum(np.ones(len(LMean_Size))) , LMean_Size) , 2) # correlation between the two arrays (symetric correlation corr(X,Y) or corr(Y,X))
round_to_N_digit = 3
R2_Size = math.floor(pow(10,round_to_N_digit) * R2_Size[0][1]) * pow(10,-round_to_N_digit)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Plotting the results  # # # # # # # # # # # # # # # # # # # # #
print('\nPlotting the results in figures, graphs and curves ...')

###########
## FIGURE 1
fig1 = plt.figure(1 , figsize = (DefaultFigureWidth , DefaultFigureHeight) , dpi = DefaultDPIResolution , facecolor = DefaultFaceColor , edgecolor = DefaultEdgeColor , layout = DefaultLayout)
plt.subplot(231)
plt.title("Particle {} over (-) arbitrary unit".format(Size_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Size)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.scatter(X_Axis_Arbitrary_Units , LMean_Size          , label = Size_label + '; r² = {}'.format(R2_Size)  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Size , alpha = Main_Curve_Transparency)
plt.scatter(X_Axis_Arbitrary_Units , LBackground_on_Size , label = Background_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Background , alpha = Background_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Size_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(232)
plt.title("Particle {} over (-) arbitrary unit".format(Number_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Number)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.scatter(X_Axis_Arbitrary_Units , LMean_Number          , label = Number_label      , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Number , alpha = Main_Curve_Transparency)
plt.scatter(X_Axis_Arbitrary_Units , LBackground_on_Number , label = Background_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Background , alpha = Background_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Number_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(233)
plt.title("Particle {} over (-) arbitrary unit".format(Concentration_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Concentration)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.scatter(X_Axis_Arbitrary_Units , LMean_Concentration           , label = Concentration_label , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Concentration , alpha = Main_Curve_Transparency)
plt.scatter(X_Axis_Arbitrary_Units , LBackground_on_Concentration  , label = Background_label    , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Background , alpha = Background_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Concentration_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(234)
plt.title("Particle {} over (-) arbitrary unit".format(Volume_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Volume)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.scatter(X_Axis_Arbitrary_Units , LMean_Volume          , label = Volume_label      , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Volume , alpha = Main_Curve_Transparency)
plt.scatter(X_Axis_Arbitrary_Units , LBackground_on_Volume , label = Background_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Background , alpha = Background_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Volume_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(235)
plt.title("Particle {} over (-) arbitrary unit".format(Area_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Area)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.scatter(X_Axis_Arbitrary_Units , LMean_Area          , label = Area_label        , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Area , alpha = Main_Curve_Transparency)
plt.scatter(X_Axis_Arbitrary_Units , LBackground_on_Area , label = Background_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Background , alpha = Background_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Area_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)

plt.show()

###########
## FIGURE 2
fig2 = plt.figure(2 , figsize = (DefaultFigureWidth , DefaultFigureHeight) , dpi = DefaultDPIResolution , facecolor = DefaultFaceColor , edgecolor = DefaultEdgeColor , layout = DefaultLayout)
plt.subplot(231)
plt.title("Particle {} over (-) arbitrary unit".format(Size_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Size)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.plot(         X_Axis_Arbitrary_Units  , LMean_Size                                      , label = Size_label                  , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Size       , alpha = Main_Curve_Transparency)
plt.fill_between( X_Axis_Arbitrary_Units  , LMean_Size + LStd_Size , LMean_Size - LStd_Size , label = 'Envelope of ' + Size_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_SizeSTD  , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Size_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(232)
plt.title("Particle {} over (-) arbitrary unit".format(Number_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Number)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.plot(         X_Axis_Arbitrary_Units  , LMean_Number                                             , label = Number_label                  , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Number     , alpha = Main_Curve_Transparency)
plt.fill_between( X_Axis_Arbitrary_Units  , LMean_Number + LStd_Number , LMean_Number - LStd_Number  , label = 'Envelope of ' + Number_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_NumberSTD  , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Number_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(233)
plt.title("Particle {} over (-) arbitrary unit".format(Concentration_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Concentration)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.plot(         X_Axis_Arbitrary_Units  , LMean_Concentration                                                                 , label = Concentration_label                  , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Concentration     , alpha = Main_Curve_Transparency)
plt.fill_between( X_Axis_Arbitrary_Units  , LMean_Concentration + LStd_Concentration , LMean_Concentration - LStd_Concentration , label = 'Envelope of ' + Concentration_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_ConcentrationSTD  , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Concentration_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(234)
plt.title("Particle {} over (-) arbitrary unit".format(Volume_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Volume)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.plot(         X_Axis_Arbitrary_Units  , LMean_Volume                                             , label = Volume_label                  , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Volume     , alpha = Main_Curve_Transparency)
plt.fill_between( X_Axis_Arbitrary_Units  , LMean_Volume + LStd_Volume , LMean_Volume - LStd_Volume  , label = 'Envelope of ' + Volume_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_VolumeSTD  , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Volume_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(235)
plt.title("Particle {} over (-) arbitrary unit".format(Area_label))
X_Axis_Arbitrary_Units = np.cumsum(np.ones(len(LMean_Area)))
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(X_Axis_Arbitrary_Units) / 100
plt.plot(         X_Axis_Arbitrary_Units  , LMean_Area                                     , label = Area_label                  , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Area     , alpha = Main_Curve_Transparency)
plt.fill_between( X_Axis_Arbitrary_Units  , LMean_Area + LStd_Area , LStd_Area - LStd_Area , label = 'Envelope of ' + Area_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_AreaSTD  , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel("(-) arbitrary unit")
plt.ylabel(Area_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)

plt.show()

###########
## FIGURE 3
if MODE_XminXmax_to_show == "relative" : Xmax_for_plots = Percentage * np.max(LMean_Size) / 100 # Setting the Xmax for every following plots

fig3 = plt.figure(3 , figsize = (DefaultFigureWidth , DefaultFigureHeight) , dpi = DefaultDPIResolution , facecolor = DefaultFaceColor , edgecolor = DefaultEdgeColor , layout = DefaultLayout)
plt.subplot(221)
plt.title("Particle {} over {}".format(Number_label , Size_label))
plt.scatter( LMean_Size , LMean_Number_minus_Background               , label = Number_label                            , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Number  , alpha = Main_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Number_minus_Background + LStd_Number , label = Number_label + ' + Standard Deviation'  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_NumberSTD , alpha = Std_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Number_minus_Background - LStd_Number , label = Number_label + ' - Standard Deviation'  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_NumberSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Number_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(222)
plt.title("Particle {} over {}".format(Concentration_label , Size_label))
plt.scatter( LMean_Size , LMean_Concentration_minus_Background                      , label = Concentration_label                           , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Concentration    , alpha = Main_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Concentration_minus_Background + LStd_Concentration , label = Concentration_label + ' + Standard Deviation' , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_ConcentrationSTD , alpha = Std_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Concentration_minus_Background - LStd_Concentration , label = Concentration_label + ' - Standard Deviation' , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_ConcentrationSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Concentration_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(223)
plt.title("Particle {} over {}".format(Volume_label , Size_label))
plt.scatter( LMean_Size , LMean_Volume_minus_Background                , label = Volume_label                           , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Volume    , alpha = Main_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Volume_minus_Background + LStd_Volume  , label = Volume_label + ' + Standard Deviation' , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_VolumeSTD , alpha = Std_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Volume_minus_Background - LStd_Volume  , label = Volume_label + ' - Standard Deviation' , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_VolumeSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Volume_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(224)
plt.title("Particle {} over {}".format(Area_label , Size_label))
plt.scatter( LMean_Size , LMean_Area_minus_Background                , label = Area_label                           , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_Area    , alpha = Main_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Area_minus_Background + LStd_Volume  , label = Area_label + ' + Standard Deviation' , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_AreaSTD , alpha = Std_Curve_Transparency)
plt.scatter( LMean_Size , LMean_Area_minus_Background - LStd_Volume  , label = Area_label + ' - Standard Deviation' , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_AreaSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Area_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)

plt.show()

###########
## FIGURE 4
fig4 = plt.figure(4 , figsize = (DefaultFigureWidth , DefaultFigureHeight) , dpi = DefaultDPIResolution , facecolor = DefaultFaceColor , edgecolor = DefaultEdgeColor , layout = DefaultLayout)
plt.subplot(221)
plt.title("Particle {} Standard Deviation over {}".format(Number_label , Size_label))
plt.fill_between( LMean_Size  , + LStd_Number , - LStd_Number , label = 'Envelope of ' + Number_label             , ls = DefaultLineStyle , lw = DefaultLineWidth                         , color = DefaultColor_Number , alpha = Std_Curve_Transparency)
plt.scatter(      LMean_Size  , + LStd_Number                 , label = '+ Standard Deviation of ' + Number_label , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_NumberSTD , alpha = Std_Curve_Transparency)
plt.scatter(      LMean_Size  , - LStd_Number                 , label = '- Standard Deviation of ' + Number_label , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_NumberSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Number_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(222)
plt.title("Particle {} Standard Deviation over {}".format(Concentration_label , Size_label))
plt.fill_between( LMean_Size , + LStd_Concentration   , - LStd_Concentration                                      , label = 'Envelope of ' + Concentration_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Concentration , alpha = Std_Curve_Transparency)
plt.scatter(      LMean_Size , + LStd_Concentration   , label = '+ Standard Deviation of ' + Concentration_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_ConcentrationSTD    , alpha = Main_Curve_Transparency)
plt.scatter(      LMean_Size , - LStd_Concentration   , label = '- Standard Deviation of ' + Concentration_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_ConcentrationSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Concentration_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(223)
plt.title("Particle {} Standard Deviation over {}".format(Volume_label , Size_label))
plt.fill_between( LMean_Size , + LStd_Volume   , - LStd_Volume                                      , label = 'Envelope of ' + Volume_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Volume , alpha = Std_Curve_Transparency)
plt.scatter(      LMean_Size , + LStd_Volume   , label = '+ Standard Deviation of ' + Volume_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_VolumeSTD    , alpha = Main_Curve_Transparency)
plt.scatter(      LMean_Size , - LStd_Volume   , label = '- Standard Deviation of ' + Volume_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_VolumeSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Volume_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)
plt.subplot(224)
plt.title("Particle {} Standard Deviation over {}".format(Area_label , Size_label))
plt.fill_between( LMean_Size , + LStd_Area , - LStd_Area                                        , label = 'Envelope of ' + Area_label , ls = DefaultLineStyle , lw = DefaultLineWidth , color = DefaultColor_Area , alpha = Std_Curve_Transparency)
plt.scatter(      LMean_Size , + LStd_Area   , label = '+ Standard Deviation of ' + Area_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_AreaSTD    , alpha = Main_Curve_Transparency)
plt.scatter(      LMean_Size , - LStd_Area   , label = '- Standard Deviation of ' + Area_label  , ls = DefaultLineStyle , lw = DefaultLineWidth , s = DefaultMarkerSize , color = DefaultColor_AreaSTD , alpha = Std_Curve_Transparency)
plt.xlim(Xmin_for_plots , Xmax_for_plots)
plt.xlabel(Size_label)
plt.ylabel(Area_label)
plt.legend()
plt.grid(visible = DefaultGridValue , which = 'both' , alpha = DefaultGridOpacity , color = DefaultGridColor , linestyle = DefaultGridLineStyle , linewidth = DefaultGridLineWidth)

plt.show()

###################################
## SHOWING ALL GRAPHS PLOTS FIGURES
plt.show()

############################
## SAVING DATA TO EXCEL FILE
print("\nSaving data to Excel file")

# Writing the complete path and Excel file name to save it later
ExcelName_COMPLETE = ExcelName + '_' + DATE_AND_TIME + '.xlsx'
ExcelPath_and_Name = os.path.join(path_to_all_files , ExcelName_COMPLETE)

# Create a dictionary containing your data arrays.
DataDictionnary_for_ExcelFile = {
  'HEADER Variable Name'                                : List_Header_Variable_Names_for_EXCEL ,
  'HEADER 1 : Variable Information'                     : List_Header_Variable_Strings_for_EXCEL[0] ,
  'HEADER 2 : Variable Information'                     : List_Header_Variable_Strings_for_EXCEL[1] ,
  'HEADER 3 : Variable Information'                     : List_Header_Variable_Strings_for_EXCEL[2] ,
  'HEADER 4 : Variable Information'                     : List_Header_Variable_Strings_for_EXCEL[3] ,
  'MEAN '       + Size_label                            : LMean_Size ,
  'MEAN '       + Size_label          + ' - Background' : LMean_Size_minus_Background , 
  'Backgroung ' + Size_label                            : LBackground_on_Size , 
  'STD '        + Size_label                            : LStd_Size , 
  'LOG '        + Size_label                            : LMean_Size_LOG , 
  'MEAN '       + Number_label                          : LMean_Number , 
  'MEAN '       + Number_label        + ' - Background' : LMean_Number_minus_Background , 
  'Backgroung ' + Number_label                          : LBackground_on_Number , 
  'STD '        + Number_label                          : LStd_Number , 
  'LOG '        + Number_label                          : LMean_Number_LOG , 
  'MEAN '       + Concentration_label                   : LMean_Concentration , 
  'MEAN '       + Concentration_label + ' - Background' : LMean_Concentration_minus_Background , 
  'Backgroung ' + Concentration_label                   : LBackground_on_Concentration , 
  'STD '        + Concentration_label                   : LStd_Concentration , 
  'LOG '        + Concentration_label                   : LMean_Concentration_LOG , 
  'MEAN '       + Volume_label                          : LMean_Volume , 
  'MEAN '       + Volume_label        + ' - Background' : LMean_Volume_minus_Background ,
  'Backgroung ' + Volume_label                          : LBackground_on_Volume , 
  'STD '        + Volume_label                          : LStd_Volume , 
  'LOG '        + Volume_label                          : LMean_Volume_LOG , 
  'MEAN '       + Area_label                            : LMean_Area , 
  'MEAN '       + Area_label          + ' - Background' : LMean_Area_minus_Background , 
  'Backgroung ' + Area_label                            : LBackground_on_Area , 
  'STD '        + Area_label                            : LStd_Area , 
  'LOG '        + Area_label                            : LMean_Area_LOG }

# Try to directly convert the dictionnary into a pandas dataframe
try     : DataFrame_for_ExcelFile = pd.DataFrame(DataDictionnary_for_ExcelFile)
except  : # Try to convert the dictionnary transposed to avoid "ValueError: All arrays must be of the same length" ERROR. Then transpose again to get back to the initial shape
  DataFrame_for_ExcelFile = pd.DataFrame.from_dict(DataDictionnary_for_ExcelFile , orient = 'index')
  DataFrame_for_ExcelFile = DataFrame_for_ExcelFile.transpose()

# Save the DataFrame to an Excel file using the to_excel method + freezing the first row
with pd.ExcelWriter(ExcelPath_and_Name , engine = 'openpyxl') as writer :
  DataFrame_for_ExcelFile.to_excel(writer , index = False , sheet_name = Excel_Sheet_Default_Name)
  writer.sheets[Excel_Sheet_Default_Name].freeze_panes = Excel_Cells_to_Freeze

print("\n\tExcel file saved : {}".format(ExcelPath_and_Name))

#################
## END OF PROGRAM
TIME_END_PROGRAM        = time.time()
ELAPSED_TIME_PROGRAM    = round(TIME_END_PROGRAM - TIME_START_PROGRAM , 2)
ELAPSED_TIME_PROGRAM_DD = int( ELAPSED_TIME_PROGRAM / (3600*24))
ELAPSED_TIME_PROGRAM_HH = int((ELAPSED_TIME_PROGRAM - ELAPSED_TIME_PROGRAM_DD * 3600 * 24) / 3600)
ELAPSED_TIME_PROGRAM_MN = int((ELAPSED_TIME_PROGRAM - ELAPSED_TIME_PROGRAM_DD * 3600 * 24 - ELAPSED_TIME_PROGRAM_HH * 3600) / 60)
ELAPSED_TIME_PROGRAM_SS = int( ELAPSED_TIME_PROGRAM - ELAPSED_TIME_PROGRAM_DD * 3600 * 24 - ELAPSED_TIME_PROGRAM_HH * 3600 - ELAPSED_TIME_PROGRAM_MN * 60)

if    ELAPSED_TIME_PROGRAM < 60         : print("Elapsed time : {}s."             .format(ELAPSED_TIME_PROGRAM_SS))
elif  ELAPSED_TIME_PROGRAM < 3600       : print("Elapsed time : {}mn {}s."        .format(ELAPSED_TIME_PROGRAM_MN , ELAPSED_TIME_PROGRAM_SS))
elif  ELAPSED_TIME_PROGRAM < 3600 * 24  : print("Elapsed time : {}h {}mn {}s."    .format(ELAPSED_TIME_PROGRAM_HH , ELAPSED_TIME_PROGRAM_MN , ELAPSED_TIME_PROGRAM_SS))
else                                    : print("Elapsed time : {}d {}h {}mn {}s.".format(ELAPSED_TIME_PROGRAM_DD , ELAPSED_TIME_PROGRAM_HH , ELAPSED_TIME_PROGRAM_MN , ELAPSED_TIME_PROGRAM_SS))