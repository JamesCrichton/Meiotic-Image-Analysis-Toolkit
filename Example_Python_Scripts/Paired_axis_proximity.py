# -*- coding: utf-8 -*-
"""
Created on Sat May 20 17:18:38 2023

@author: jcc234
"""

# Measuring the proximity of paired axes
# Data:
# Users must have first traced axes in their images, renaming the traces with a number shared by both axes in a pair, and a letter unique to each e.g. "1a" and "1b" for two paired chromosomes. Isolated traces can be included, and information regarding positions along the axis will be calculated, however distances to a partner will not.

# Trace files from the SNT ImageJ plugin must be saved in the .traces format.

# Downloads:
# Go to https://github.com/huangziwei/TracePy and clone this repository. You can simply select "Code", "Download ZIP", and after downloading, extract these files

# Specify paths to package locations and image data to be analysed:
# Set the TracePy_path and Meio_toolkit_path paths below to the folders where you saved the TracePy and meiosis_toolkit packages.

# Set the img_library path below to the folder that contains the subfolders of image data to be analysed (trace files in .traces format). The script will cycle through the subfolders in the img_library location, all subfolders containing image metadata will be analysed in one batch.

import os
import sys
##############################################################################################
############### PATHS TO BE MANUALLY EDITED ##################################################
##############################################################################################

#Automatically setting paths to repositories downloaded to the desktop, and image directories. 
desktop = os.path.expanduser("~/Desktop")
Meio_toolkit_path=os.path.join(desktop,"Meiotic-Image-Analysis-Toolkit-main")
TracePy_path= os.path.join(desktop,"TracePy-master")
img_library=os.path.join(Meio_toolkit_path,"SIM_sample_data")

# #Alternatively paths to repositories and image directories can be set manually. This is required if for example a OneDrive Desktop is being used
# #Note: on Windows systems use two backslashes or a forward slash as a file separator when setting paths e.g.
# TracePy_path= "C:/Users/jcc234/OneDrive - University of Exeter/Desktop/TracePy-master"
# Meio_toolkit_path="C:/Users/jcc234/OneDrive - University of Exeter/Desktop/Meiotic-Image-Analysis-Toolkit-main"
# img_library="C:/Users/jcc234/OneDrive - University of Exeter/Meiotic-Image-Analysis-Toolkit-main/SIM_sample_data"

##############################################################################################
##############################################################################################
##############################################################################################

#check if these automatic paths are valid
if not os.path.exists(Meio_toolkit_path):
    sys.exit("Error! Meio_toolkit_path not valid. Define manually")
   
if not os.path.exists(TracePy_path):
    sys.exit("Error! TracePy_path not valid. Define manually")
    
if not os.path.exists(img_library):
    sys.exit("Error! img_library path not valid. Define manually")

#Install required modules
current_wd=os.getcwd()
os.chdir(Meio_toolkit_path)# Navigate to the directorty to access the requirements file. Not working well with longer paths
!pip install -r requirements.txt
os.chdir(current_wd)


# Run paired axis proximity quantification
# This will output results as "Axis_Trace_Measurements.csv" files in the image metadata subfolders.

## Add package addresses to the environment
import sys

## Add package addresses to the environment
if not TracePy_path in sys.path:
    sys.path.append(TracePy_path)

if not Meio_toolkit_path in sys.path:
    sys.path.append(Meio_toolkit_path)

## Import packages
import meiosis_toolkit
import fnmatch

## Identify image data and process
Folders=next(os.walk(img_library))[1]

for Folder in Folders:
    Folder_Path=os.path.join(img_library,Folder)
    trace_found=0
    foci_data_found=0

    #Only proceed if .traces file present for an image in its metadata folder
    Traces_found=any('traces' in file for file in os.listdir(Folder_Path))
            
    if Traces_found:
        for file in os.listdir(Folder_Path):
            if fnmatch.fnmatch(file,'*.traces'):
                trace_found=1
                trace_path=os.path.join(Folder_Path,file)
                trace_df=meiosis_toolkit.Axis_Proximity(trace_path)
                                      
                trace_df["Image"]=Folder#Add the codename to entries in this file. Will help later if concatenating 
                trace_df_path=os.path.join(Folder_Path,"Axis_Trace_Measurements.csv")
                trace_df.to_csv(trace_df_path)       
        else:
            print("ERROR! No .traces file found for "+file+". Is the data in the correct location? Require a collection of metadata folders named in the format \"imagename_Output\" and containing .traces files, within the location you defined:" +img_library)
                
if len(Folders)==0:
    print("ERROR! NO FOLDERS DETECTED. Is the data in the correct location? Require a collection of metadata folders named in the format \"imagename_Output\" and containing .traces files, within the location you defined:" +img_library)