# -*- coding: utf-8 -*-
"""
Created on Fri May 19 10:24:28 2023

@author: jcc234
"""


# Calculating focus positioning on traced chromosome axes
# Data:
# Axial foci should be segmented and measured using the "Image Segmentation" ImageJ macro. This macro outputs foci measurements that are needed for this script, the default name of the foci measurements file is: Watershed_Foci_measurements_to_Axes_Mask.csv

# Users must also first trace chromosome axes using the SNT plugin in ImageJ, to produce a ".traces" file needed for this script.

# Downloads:
# Go to https://github.com/huangziwei/TracePy and clone this repository. You can simply select "Code", "Download ZIP", and after downloading, extract these files

# Specify package and data locations:
# Set the TracePy_path and Meio_toolkit_path paths below to the folders where you saved the TracePy and meiosis_toolkit packages.

# Set the img_library path below to the folder that contains the subfolders of image data to be analysed (trace files in .traces format). The script will cycle through the subfolders in the img_library location, all subfolders containing image metadata will be analysed in one batch.

# The Foci_data variable below should be set to match the name of the foci measurements output file generated by the "Image Segmentation" ImageJ macro. The default output name for this file is "Watershed_Foci_measurements_to_Axes_Mask.csv".

import os
import sys

#Automatically setting paths to repositories downloaded to the desktop, and image directories. 
desktop = os.path.expanduser("~/Desktop")
Meio_toolkit_path=os.path.join(desktop,"Meiotic-Image-Analysis-Toolkit-main")
TracePy_path= os.path.join(desktop,"TracePy-master")
img_library=os.path.join(Meio_toolkit_path,"Widefield_sample_data")

Foci_data="Watershed_Foci_measurements_to_Axes_Mask.csv"

# #Alternatively paths to repositories and image directories can be set manually. This is required if for example a OneDrive Desktop is being used
# #Note: on Windows systems use two backslashes or a forward slash as a file separator when setting paths e.g. "C:\Users\me\Path to folder"

# TracePy_path= "C:/Users/jcc234/OneDrive - University of Exeter/Desktop/TracePy-master"
# Meio_toolkit_path="C:/Users/jcc234/OneDrive - University of Exeter/Desktop/Meiotic-Image-Analysis-Toolkit-main"
# img_library="C:/Users/jcc234/OneDrive - University of Exeter/Meiotic-Image-Analysis-Toolkit-main/Widefield_sample_data"
# Foci_data="Watershed_Foci_measurements_to_Axes_Mask.csv"

#check if these automatic paths are valid
if not os.path.exists(Meio_toolkit_path):
    sys.exit("Error! Meio_toolkit_path not valid. Define manually")
   
if not os.path.exists(TracePy_path):
    sys.exit("Error! TracePy_path not valid. Define manually")
    
if not os.path.exists(img_library):
    sys.exit("Error! img_library path not valid. Define manually")
    
# Run position measurement script:
# This will output results as "Watershed_Foci_measurements_to_Axes_Mask_and_positioning" files and modify the "Axis_Trace_Measurements.csv" files in the image metadata subfolders.

## Add package addresses to the environment
sys.path.append(TracePy_path)
sys.path.append(Meio_toolkit_path)

## Import packages
import meiosis_toolkit
import fnmatch

## Identify image data and process
Folders=next(os.walk(img_library))[1]

for Folder in Folders:
    Folder_Path=img_library+"/"+Folder+"/"
    trace_found=0
    foci_data_found=0

    #Find the .traces file and foci measurement data for an image in its metadata folder
    #Only proceed if both are found 
    for file in os.listdir(Folder_Path):
        if fnmatch.fnmatch(file,'*.traces'):
            trace_found=1
            trace_path=Folder_Path+file
        if fnmatch.fnmatch(file,Foci_data):
            foci_data_found=1
            foci_path=os.path.join(Folder_Path+Foci_data)
        if (foci_data_found & trace_found):
            trace_df, Focus_Info=meiosis_toolkit.Focus_Position_on_Trace(trace_path, foci_path)
                                  
            trace_df["Image"]=Folder#Add the codename to entries in this file. Will help later if concatenating 
            trace_data_path=os.path.join(Folder_Path+"Axis_Trace_Measurements.csv")
            trace_df.to_csv(trace_data_path)
            name=Foci_data.replace(".csv","")
            focus_info_path=os.path.join(Folder_Path, (name+"_and_positioning.csv"))
            Focus_Info.to_csv(focus_info_path)
 