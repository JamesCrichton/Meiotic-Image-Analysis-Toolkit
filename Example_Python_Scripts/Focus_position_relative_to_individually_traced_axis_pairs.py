# -*- coding: utf-8 -*-
"""
Created on Sat May 20 17:43:34 2023

@author: jcc234
"""

# Focus position relative to individually traced homologue axes
# This program calculates the position of segmented foci relative to individually traced homologus axes (typically asynapsed and/or super-resolution images).

# Data:
# Images must have been segmented using the "Image Segmentation" ImageJ macro to generate a foci labelmap (default filename "Watershed_Foci.tif").

# Users will also first need to calculate the proximity of paired axes by running the "Paired_axis_proximity.py" python script to produce an "Axis_Trace_Measurements.csv" file.

# Traces generated from the SNT ImageJ plugin in .traces format should be present in the metadata subfolders, and the original image file should also be present in the parent folder containing the metadata subfolders.

# Specify package and data locations:
# Set the Meio_toolkit_path path below to the folder where you saved the meiosis_toolkit package.

# Set the img_library path below to the folder that contains the subfolders of image data to be analysed (trace files in .traces format). The script will cycle through the subfolders in the img_library location, all subfolders containing image metadata will be analysed in one batch.

# This script needs a foci labelmap .tif outputted by the "Image Segmentation" Image J macro (named "Watershed_Foci.tif" by default). In the variables below, specify the name for this file (Foci_Labelmap_name="RAD51.tif"), its channel in the original image (Foci_Channel=2), and the name to append to the output filenames (e.g. Foci_name="RAD51").

import os
import sys
##############################################################################################
############### PATHS AND VARIABLES TO BE MANUALLY EDITED ####################################
##############################################################################################

#Automatically setting paths to repositories downloaded to the desktop, and image directories. 
desktop = os.path.expanduser("~/Desktop")
Meio_toolkit_path=os.path.join(desktop,"Meiotic-Image-Analysis-Toolkit-main")
img_library=os.path.join(Meio_toolkit_path,"SIM_sample_data")

Foci_Labelmap_name="RAD51.tif"
Foci_Channel=2
Foci_name="RAD51"

# #Alternatively paths to repositories and image directories can be set manually. This is required if for example a OneDrive Desktop is being used
# #Note: on Windows systems use two backslashes or a forward slash as a file separator when setting paths e.g:
# Meio_toolkit_path="C:/Users/jcc234/OneDrive - University of Exeter/Desktop/Meiotic-Image-Analysis-Toolkit-main"
# img_library="C:/Users/jcc234/OneDrive - University of Exeter/Meiotic-Image-Analysis-Toolkit-main/SIM_sample_data"

# Specify foci filters
# Foci included in analysis are filtered by area measured in pixels. Inspect your images to determine reasonable minimum and maximum limits for foci area, and set the Foci_area_min and Foci_area_max variables below.
Foci_area_min=5
Foci_area_max=100

##############################################################################################
##############################################################################################
##############################################################################################

#check if these automatic paths are valid
if not os.path.exists(Meio_toolkit_path):
    sys.exit("Error! Meio_toolkit_path not valid. Define manually")
   
if not os.path.exists(img_library):
    sys.exit("Error! img_library path not valid. Define manually")


# Run position measurement script:
# The focus positions will be saved in the image metadata subfolder with the suffix "homologous_axis_measurements.csv"

#Install required modules
current_wd=os.getcwd()
os.chdir(Meio_toolkit_path)# Navigate to the directorty to access the requirements file. Not working well with longer paths
!pip install -r requirements.txt 
os.chdir(current_wd)
 
## Add package addresses to the environment

sys.path.append(Meio_toolkit_path)

## Import packages
import pandas as pd
import meiosis_toolkit
from skimage import io
from timeit import default_timer


## Identify image data and process
Folders=next(os.walk(img_library))[1]

for Folder in Folders:
    start = default_timer()
    Folder_Path=os.path.join(img_library,Folder)
    
    #Only process if the necessary files are present in an image metadata folder
    Foci_found=any(Foci_Labelmap_name in file for file in os.listdir(Folder_Path))
    Traces_found=any("Axis_Trace_Measurements.csv" in file for file in os.listdir(Folder_Path))
    if (Foci_found & Traces_found):
        Foci_path=os.path.join(Folder_Path,Foci_Labelmap_name)
        Traces_path=os.path.join(Folder_Path,"Axis_Trace_Measurements.csv")
        #Open files
        Foci_Labelmap=io.imread(Foci_path)
        Trace_df=pd.read_csv(Traces_path)
        
        Foci_Channel_Corrected=Foci_Channel-1
        img=Folder.replace("_Output", "")
        Image_Path=os.path.join(img_library,img+".tif")
        Foci_Intensity_Img=io.imread(Image_Path)[:,:,Foci_Channel_Corrected]
        
        Foci_Measurements=meiosis_toolkit.Resolved_homologues_focus_position(Trace_df, Foci_Labelmap, Foci_Intensity_Img, Foci_area_min, Foci_area_max)
        
        #Save output
        output_path=os.path.join(Folder_Path,Foci_name+"_homologous_axis_measurements.csv")
        Foci_Measurements.to_csv(output_path)
        
        print(Folder+" processed")
    
    else:
        print("ERROR! No .traces or focus labelmap file found in "+Folder+". Is the data in the correct location? Require a collection of metadata folders named in the format \"imagename_Output\" and containing .traces and watershed label map files you've named "+Foci_name+", within the location you defined:" +img_library)
            
if len(Folders)==0:
    print("ERROR! NO FOLDERS DETECTED. Is the data in the correct location? Require a collection of metadata folders named in the format \"imagename_Output\" and containing .traces and watershed label map files you've named "+Foci_name+", within the location you defined:" +img_library)