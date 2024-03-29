# -*- coding: utf-8 -*-
"""
Created on Sat May 20 15:54:06 2023

@author: jcc234
"""

# Axial shuffling of foci
# To determine whether different types of foci are genuinely associated with one another it can be informative to artificially shuffle segmented foci around the axial space and compare the observed and randomly-distributed levels of colocalisation/proximity.

# Data:
# Images must have been segmented using the "Image Segmentation" ImageJ macro to generate nuclear masks ("Nuclear_Mask.tif"), axial masks ("Axial_Mask.tif") and foci labelmaps (.tif files). These default file names can be altered below, if others have been used. Axial shuffling requires two focus labelmaps, one (Label_1) remains static and the other (Label_2) is shuffled around the axis mask. Distances from each shuffled Label_2 focus to the closest reference Label 1 focus are measured.

# Specify package and data locations:
# Set the Meio_toolkit_path path below to the folder where you saved the meiosis_toolkit package.

# Set the img_library path below to the folder that contains the subfolders of image data to be analysed (trace files in .traces format). The script will cycle through the subfolders in the img_library location, all subfolders containing image metadata will be analysed in one batch.

# This script needs two different foci labelmaps outputted by the "Image Segmentation" Image J macro. In the variables below, specify the name for the Label_1 (reference) and Label_2 (to be shuffled) foci labelmap files (.tif file format, e.g. Label_1="RPA.tif"), the channel in the original image that each of Label_1 and Label_2 correspond to (e.g. Label_1_Channel=2), and the name for Label_1 and Label_2 to append to the output filenames (e.g. Label_1_name="RAD51").

# The axes_mask variable below should be set to match the names of the axes mask output file generated by the "Image Segmentation" ImageJ macro. The default output names for this file is "Axes_Mask.tif".

import os
import sys 

##############################################################################################
############### PATHS AND VARIABLES FOR MANUAL EDITING #######################################
##############################################################################################


#Automatically setting paths to repositories downloaded to the desktop, and image directories. 
desktop = os.path.expanduser("~/Desktop")
Meio_toolkit_path=os.path.join(desktop,"Meiotic-Image-Analysis-Toolkit-main")
img_library=os.path.join(Meio_toolkit_path,"SIM_sample_data")

# #Alternatively paths to repositories and image directories can be set manually. This is required if for example a OneDrive Desktop is being used
# #Note: on Windows systems use two backslashes or a forward slash as a file separator when setting paths e.g.
# Meio_toolkit_path="C:/Users/jcc234/OneDrive - University of Exeter/Desktop/Meiotic-Image-Analysis-Toolkit-main"
# img_library="C:/Users/jcc234/OneDrive - University of Exeter/Desktop/Meiotic-Image-Analysis-Toolkit-main/SIM_sample_data"

#check if these automatic paths are valid
if not os.path.exists(Meio_toolkit_path):
    sys.exit("Error! Meio_toolkit_path not valid. Define manually")
      
if not os.path.exists(img_library):
    sys.exit("Error! img_library path not valid. Define manually")


Label_1="RAD51.tif"
Label_1_Channel=2
Label_1_name="RAD51"
Label_2="SYCE2.tif"
Label_2_Channel=4
Label_2_name="SYCE2"
axes_mask="Axes_Mask.tif"

# Specify foci filters and shuffling parameters
# Foci included in analysis are filtered by area measured in pixels. Inspect your images to determine reasonable minimum and maximum limits for foci area, and set the area_min and area_max variables for Label_1 and Label_2 foci below.

# The number of times that the script will shuffle the foci data in each image can be set using the shuffles variable below

Label_1_area_min=5
Label_1_area_max=100
Label_2_area_min=5
Label_2_area_max=100
shuffles=10

##############################################################################################
##############################################################################################
##############################################################################################

#Install required modules
current_wd=os.getcwd()
os.chdir(Meio_toolkit_path)# Navigate to the directorty to access the requirements file. Not working well with longer paths
!pip install -r requirements.txt 
os.chdir(current_wd)
 

# Run axial foci shuffling:
# This script will output multiple .csv and .tif files to the image metadata subfolders

## Add package addresses to the environment
if not Meio_toolkit_path in sys.path:
    sys.path.append(Meio_toolkit_path)

## Import packages
import meiosis_toolkit
from skimage import io
from timeit import default_timer

#Correcting the channel numbering for Python
Label_1_Channel_Corrected=Label_1_Channel-1
Label_2_Channel_Corrected=Label_2_Channel-1

## Identify image data and process
Folders=next(os.walk(img_library))[1]

for Folder in Folders:
    start = default_timer()
    Folder_Path=os.path.join(img_library,Folder)
    
    #Only process if the necessary files are present in an image metadata folder
    nucleus_found=any("Nuclear_Mask.tif" in file for file in os.listdir(Folder_Path))
    axes_found=any(axes_mask in file for file in os.listdir(Folder_Path))
    Label1_found=any(Label_1 in file for file in os.listdir(Folder_Path))
    Label2_found=any(Label_2 in file for file in os.listdir(Folder_Path))
    if (nucleus_found & axes_found & Label1_found & Label2_found):
        nuclear_mask_path=os.path.join(Folder_Path,"Nuclear_Mask.tif")
        axis_mask_path=os.path.join(Folder_Path,axes_mask)
        Label_1_path=os.path.join(Folder_Path,Label_1)
        Label_2_path=os.path.join(Folder_Path,Label_2)
        #Open label and mask images
        Label1=io.imread(Label_1_path)
        Label2=io.imread(Label_2_path)
        Mask=io.imread(axis_mask_path)
        #Open intensity images
        img=Folder.replace("_Output", "")
        Image_Path=os.path.join(img_library,img+".tif")
        img_file=io.imread(Image_Path)
        Label1_intensity=img_file[:,:,Label_1_Channel_Corrected]
        Label2_intensity=img_file[:,:,Label_2_Channel_Corrected]
        
        Shuffled_centroid_distances, Original_centroid_distances, Original_weighted_centroid_distances, Selected_Label_1, Selected_Label_2,New_Label_2=meiosis_toolkit.Axial_Shuffle(Mask, Label1, Label2, Label1_intensity, Label2_intensity, Label_1_area_min, Label_1_area_max, Label_2_area_min, Label_2_area_max, shuffles)

        #Save output
        shuffled_centroid_path=os.path.join(Folder_Path,"Axis_shuffled_"+Label_2_name+"_to_"+Label_1_name+".csv")
        Shuffled_centroid_distances.to_csv(shuffled_centroid_path)#Size selected focus distances to axis
        
        original_entroid_path=os.path.join(Folder_Path,"Original_centroid_"+Label_2_name+"_to_"+Label_1_name+".csv")
        Original_centroid_distances.to_csv(original_entroid_path)#Size selected focus distances to axis
        
        original_weighted_centroid_path=os.path.join(Folder_Path,"Original_weighted_centroid_"+Label_2_name+"_to_"+Label_1_name+".csv")
        Original_weighted_centroid_distances.to_csv(original_weighted_centroid_path)#Size selected focus distances to axis
        
        size_selected_Label1_path=os.path.join(Folder_Path,"Size-selected_"+Label_1_name+"foci_labelmap.tif")
        io.imsave(size_selected_Label1_path, Selected_Label_1, check_contrast =False)#Focus labelmap image of the size-selected foci included in processing
        
        size_selected_Label2_path=os.path.join(Folder_Path,"Size-selected_"+Label_2_name+"foci_labelmap.tif")
        io.imsave(size_selected_Label2_path, Selected_Label_2, check_contrast =False)#Focus labelmap image of the size-selected foci included in processing
        axis_shuffled_path=os.path.join(Folder_Path,"Axis-Shuffled_"+Label_2_name+".tif")
        io.imsave(axis_shuffled_path, New_Label_2, check_contrast =False)#Example image of the shuffled foci in nuclear space

        duration = default_timer() - start

        print(Folder+" processed ("+ str(duration) +" seconds).")
        
    else:
        print("ERROR! Files not detected. Is the data in the correct location? Require a collection of metadata folders named in the format \"imagename_Output\" within the location you defined:" +img_library)
        
if len(Folders)==0:
    print("ERROR! NO FOLDERS DETECTED. Is the data in the correct location? Require a collection of metadata folders named in the format \"imagename_Output\" within the location you defined:" +img_library)

