import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import scipy.stats
import scipy.ndimage
import pandas as pd
from skimage import data, img_as_float,color, img_as_ubyte, color, io, measure
from skimage.measure import label, regionprops
import copy
import random
import scipy.spatial 
from scipy.stats import ks_2samp,fisher_exact
from math import log10, floor
import statistics

#Function to measure proximity foci to chromosome axes
def Centroid_Proxomity_and_Overlap(Mask,Label,size_min,size_max): #Label1 is the also add in size restictions

    Size=Label.shape#Make a new array to fill, which has same shape as other images
    Mask_Map=np.zeros((Size))
    SignalMask=Mask>0 #outputs T or F 
    Mask_Map[SignalMask]=1#Make binary image in case it started as a labelmap

    #Get full size Label1 map which overlap Mask
    LabelMap_with_Mask=copy.copy(Label)
    LabelMap_with_Mask[~SignalMask]=0 #makes all labels outside of Maks into "0"
    LabelIDsonMask=np.unique(LabelMap_with_Mask)#all label IDs overlapping with axis 
    LabelIDsonMask=LabelIDsonMask[LabelIDsonMask>0]#greater than 0, as this is background

    #Need to select foci with labels listed in LabelIDsonMask
    Label_Mask_Full=np.zeros((Size),dtype=int)
    for label in LabelIDsonMask:
        Focus_Locus=Label==label
        Label_Mask_Full[Focus_Locus]=label

    #Select foci with areas within size range specified
    Label_Mask_Full_selected=np.zeros((Size),dtype=int)
    Starting_Label_props=regionprops(Label_Mask_Full)#get properties of each label
    Selected_Areas_array_label=[]
    Selected_IDs_array_label=[]
    for reg in Starting_Label_props:
        if reg.area>=size_min:
            if reg.area<=size_max:
                Selected_Areas_array_label.append(reg.area)
                Selected_IDs_array_label.append(reg.label)
    for label in Selected_IDs_array_label:
        Focus_Locus=Label==label
        Label_Mask_Full_selected[Focus_Locus]=label

    Label_Mask_Full=Label_Mask_Full_selected

    #Caluclate overlapping Label2 with Label1 from this Mask-associated proportion of Label2     
    Labs_Total=len(np.unique(Label_Mask_Full)[np.unique(Label_Mask_Full)>0])#unique values will include zeros so remove
    OriginalOverlap=copy.copy(Label_Mask_Full)
    Label1_Sites=Label1>0 #outputs T or F 
    OriginalOverlap[~Label1_Sites]=0
    Original_Coloc=len(np.unique(OriginalOverlap)[np.unique(OriginalOverlap)>0])
    Original_NoColoc=Labs_Total-Original_Coloc

    #2. Select one pixel (currently centroid, but could be a random selection) in this label to become a point overlapping with the SCP3 mask 

    Label2_props=regionprops(Label2_Mask_Full, intensity_image=Intensity2)#get properties of each label
    Label1_props=regionprops(Label1, intensity_image=Intensity1)#get properties of each label
    ##need to extract the measurements for each label from this 
    lab_array_2 = []
    cen_array_2 = []
    weighted_cen_array_2=[]
    for reg in Label2_props:
        cen_array_2.append(reg.centroid)
        weighted_cen_array_2.append(reg.weighted_centroid)
        lab_array_2.append(reg.label)

    lab_array_1 = []
    cen_array_1 = []
    weighted_cen_array_1=[]
    for reg in Label1_props:
        cen_array_1.append(reg.centroid)
        weighted_cen_array_1.append(reg.weighted_centroid)
        lab_array_1.append(reg.label)

    # Calculate closest neighbour distances between Label2 and Label1 foci
    Label1_Tree=scipy.spatial.KDTree(weighted_cen_array_1)#Make a KD tree to order the points
    KDTreeLabel2toLabel1=pd.DataFrame(columns=["Query_Location","Distance","Neighbour_Position","Neighbour_Value"])
    for Focus in weighted_cen_array_2:
        d,i=Label1_Tree.query(Focus,k=1)#provides distance from "test" to closest neighbour, and the index of this in the array
        Neighbour=weighted_cen_array_1[i]
        Neighbour_Value=lab_array_1[i]
        KD_temp_df=pd.DataFrame([[Focus,d,Neighbour,Neighbour_Value]],columns=["Query_Location","Distance","Neighbour_Position","Neighbour_Value"])
        KDTreeLabel2toLabel1=pd.concat([KDTreeLabel2toLabel1,KD_temp_df])

    return KDTreeLabel2toLabel1, Original_Coloc, Original_NoColoc,Label2_Mask_Full# self explanatory. Label2_Mask_Full is size and mask-selected foci

#In the shuffling QC need to check if new focus overlaps old other new ones. If they don't fall completely within the image frame this prevents identifying overlap, so here I trim the focus to fit the frame
def Trim_Location(Coordinates,X_limit,Y_limit):
    In=Coordinates[0]<=X_limit
    Trimmed=Coordinates[0][In],Coordinates[1][In]
    In=Trimmed[0]>=0
    Trimmed=Trimmed[0][In],Trimmed[1][In]
    In=Trimmed[1]<=Y_limit
    Trimmed=Trimmed[0][In],Trimmed[1][In]
    In=Trimmed[1]>=0
    Trimmed=Trimmed[0][In],Trimmed[1][In]

    return Trimmed


def Axial_Shuffle(Mask,Label1,Label2,Intensity1,Intensity2,Label1_size_min,Label1_size_max,Label2_size_min,Label2_size_max, Shuffles):#Mask (can be labelmap), labels to be compared to, labels to be used for shuffling, intensity images for either set of labels, max and min areas for focus selection for each set of foci, number of shuffles
     
    RandomOverlaps=[]
    RandomNonOverlaps=[]
    Label2ShuffledRecordTotal=pd.DataFrame(columns=["Coords","Label","Focus_Number","Shuffle_Number"])
    Size=Label2.shape
    Label1Sites=Label1>0 #outputs T or F 

    Mask_Map=np.zeros((Size))
    SignalMask=Mask>0 #outputs T or F 
    Mask_Map[SignalMask]=1#Make binary image in case it started as a labelmap

    Mask_x,Mask_y=np.where(Mask_Map==1)#splits the x and y coordinates
    Mask_xy=list(zip(Mask_x,Mask_y))#makes xy coords into a tuple for each pixel
    Mask_xy=np.array(Mask_xy)#convert to a numpy array for maths

    #Get full size Label2 map which overlaps Mask
    Label2Map_with_Mask=copy.copy(Label2)
    Label2Map_with_Mask[~SignalMask]=0 #makes all labels outside of Masks into "0"
    Label2IDsonMask=np.unique(Label2Map_with_Mask)#all label IDs overlapping with axis 
    Label2IDsonMask=Label2IDsonMask[Label2IDsonMask>0]#greater than 0, as this is background

    #Select foci with areas within size range specified
    #Need to select foci with labels listed in Label2IDsonMask  
    #Label2 size selection
    Label2_Mask_Full=np.zeros((Size),dtype=int)
    for label in Label2IDsonMask:
            Focus_Locus=Label2==label
            Label2_Mask_Full[Focus_Locus]=label

    #Select foci with areas within size range specified
    Label2_Mask_Full_selected=np.zeros((Size),dtype=int)
    Starting_Label2_props=regionprops(Label2_Mask_Full)#get properties of each label
    Selected_Areas_array_label2=[]
    Selected_IDs_array_label2=[]
    for reg in Starting_Label2_props:
        if reg.area>=Label2_size_min:
            if reg.area<=Label2_size_max:
                Selected_Areas_array_label2.append(reg.area)
                Selected_IDs_array_label2.append(reg.label)
    for label in Selected_IDs_array_label2:
        Focus_Locus=Label2==label
        Label2_Mask_Full_selected[Focus_Locus]=label

    Label2_Mask_Full=Label2_Mask_Full_selected
    Labs_Total=len(np.unique(Label2_Mask_Full)[np.unique(Label2_Mask_Full)>0])#unique values will include zeros so remove
    
    #Label1 size selection
    Label1_Mask_Full_selected=np.zeros((Size),dtype=int)
    Starting_Label1_props=region_props(Label1)
    Selected_Areas_array_label1=[]
    Selected_IDs_array_label1=[]
    for reg in Starting_Label1_props:
        if reg.area>=Label1_size_min:
            if reg.area<=Label1_size_max:
                Selected_Areas_array_label1.append(reg.area)
                Selected_IDs_array_label1.append(reg.label)
    for label in Selected_IDs_array_label1:
        Focus_Locus=Label1==label
        Label1_Mask_Full_selected[Focus_Locus]=label

    Label1_Mask_Full=Label1_Mask_Full_selected
    
    Label2_props=regionprops(Label2_Mask_Full, intensity_image=Intensity2)#get properties of each label
    Label1_props=regionprops(Label1_Mask_Full, intensity_image=Intensity1)#get properties of each label
    
      
    #need to extract the measurements for each label from this 
    lab_array_2 = []
    cen_array_2 = []
    weighted_cen_array_2=[]
    for reg in Label2_props:
        cen_array_2.append(reg.centroid)
        weighted_cen_array_2.append(reg.weighted_centroid)
        lab_array_2.append(reg.label)

    lab_array_1 = []
    cen_array_1 = []
    weighted_cen_array_1=[]
    for reg in Label1_props:
        cen_array_1.append(reg.centroid)
        weighted_cen_array_1.append(reg.weighted_centroid)
        lab_array_1.append(reg.label)

    Label1_Tree=scipy.spatial.KDTree(weighted_cen_array_1)#Make a KD tree to order Label1 points
    New_coord_details_total=[]
    KDTreeShuffledLabel2toLabel1=pd.DataFrame()
    for x in range(Shuffles):
        #Make a new array to fill, which has same shape as other images
        New_Label2=np.zeros((Size),dtype=int)
        #New_Label2_centres=np.zeros((Size),dtype=int)
        New_coord_details=[]
        Label2ShuffledRecordIndividual=pd.DataFrame(columns=["Coords","Label","Focus_Number","Shuffle_Number"])

        for number in range(len(Selected_IDs_array_label2)):#going from 0 
            lab=Selected_IDs_array_label2[number]
            coords=np.where(Label2_Mask_Full==lab)#coordinates of label
            x_cen,y_cen =cen_array_2[number]#unpack the tuple of the coordinates
            x_cen=int(x_cen)#make the values into integers from floats
            y_cen=int(y_cen)

            #calculate relative positions of the labels coordinates around the centroid
            rel_coords=(coords[0]-x_cen,coords[1]-y_cen)    
            #make a random selection from mask to use as the new centroid for a given focus
            New_Centroid=random.choice(Mask_xy)
            New_Focus_Coords=(rel_coords[0]+New_Centroid[0],rel_coords[1]+New_Centroid[1])#moves all the coordinates of the shape, so the centroid is the newly selected pixel from within the SCP3 domain

            limit_x, limit_y=Size
            limit_x=limit_x-1#Size starts from 1, but coordinates start from 0. Need to adjust
            limit_y=limit_y-1
            NewMax_x=max(New_Focus_Coords[0])
            NewMax_y=max(New_Focus_Coords[1])
            NewMin_x=min(New_Focus_Coords[0])
            NewMin_y=min(New_Focus_Coords[1])
            #sum(New_Label2[New_Focus_Coords])!=0
            condition2=limit_x<NewMax_x
            condition3=limit_y<NewMax_y
            condition4=0<NewMin_x
            condition5=0<NewMin_y
            TrimmedCoords=Trim_Location(New_Focus_Coords,limit_x,limit_y)

            while condition2 or condition3 or condition4 or condition5 or sum(New_Label2[TrimmedCoords])!=0:
                New_Centroid=random.choice(Mask_xy)
                New_Focus_Coords=(rel_coords[0]+New_Centroid[0],rel_coords[1]+New_Centroid[1])#moves all the coordinates of the shape, so the centroid is the newly selected pixel from within the SCP3 domain
                NewMax_x=max(New_Focus_Coords[0])
                NewMax_y=max(New_Focus_Coords[1])
                NewMin_x=min(New_Focus_Coords[0])
                NewMin_y=min(New_Focus_Coords[1])
                TrimmedCoords=Trim_Location(New_Focus_Coords,limit_x,limit_y)
                condition2=limit_x<NewMax_x
                condition3=limit_y<NewMax_y
                condition4=0>NewMin_x
                condition5=0>NewMin_y

    
            New_Label2[New_Focus_Coords]=lab#give these pixels the label value
            #Calculate new weighted centroid 
            x_weighted_cen, y_weighted_cen=weighted_cen_array_2[number]#also get the weighted centroid coordinates to mark this point for future measurements relative to the centroid by shape
            rel_weighted_cen_x=int(x_weighted_cen-x_cen)
            rel_weighted_cen_y=int(y_weighted_cen-y_cen)
            New_weighted_cen=(rel_weighted_cen_x+New_Centroid[0],rel_weighted_cen_y+New_Centroid[1])
            #New_Label2_centres[New_weighted_cen]=lab#adds one pixel with the label number to the new location

            df1=pd.DataFrame([[New_weighted_cen,lab,number,x]],columns=["Coords","Label","Focus_Number","Shuffle_Number"])#new df with details of coords, label, focus numer and shuffle number
            Label2ShuffledRecordTotal=pd.concat([Label2ShuffledRecordTotal,df1])
            Label2ShuffledRecordIndividual=pd.concat([Label2ShuffledRecordIndividual,df1])

        #Select Label2 which overlap Label1
        ShuffledOverlap=copy.copy(New_Label2)
        ShuffledOverlap[~Label1Sites]=0
        Shuffled_Coloc=len(np.unique(ShuffledOverlap)[0<np.unique(ShuffledOverlap)])#exclude zeros, this is unique but background
        Shuffled_NoColoc=Labs_Total-Shuffled_Coloc

        RandomOverlaps.append(Shuffled_Coloc)
        RandomNonOverlaps.append(Shuffled_NoColoc)
        New_coord_details_total.append(New_coord_details)

        #Caluclate the nearest neighbours for all the shuffled Label2 to Label1
    for Focus in Label2ShuffledRecordTotal["Coords"]:
        d,i=Label1_Tree.query(Focus,k=1)#provides distance from "test" to closest neighbour, and the index of this in the array
        Neighbour=weighted_cen_array_1[i]
        Neighbour_Value=lab_array_1[i]
        KD_temp_df=pd.DataFrame([[x,Focus,d,Neighbour,Neighbour_Value]],columns=["Shuffle","Query_Location","Distance","Neighbour_Position","Neighbour_Value"])
        KDTreeShuffledLabel2toLabel1=pd.concat([KDTreeShuffledLabel2toLabel1,KD_temp_df])        

        
    return KDTreeShuffledLabel2toLabel1, RandomOverlaps, RandomNonOverlaps,New_Label2



##Measure overlap between foci and mask
def Focus_Centroid_To_Object_Proxomity_and_Overlap(Mask,Label,Intensity,size_min,size_max): #includes size/shape restictions for foci

    Size=Label.shape#Make a new array to fill, which has same shape as other images

    #Select foci with areas within size range specified
    Label_Size_Selected=np.zeros((Size),dtype=int)
    Label_props=regionprops(Label,Intensity)#get properties of each label
    Selected_Areas_array=[]
    Selected_IDs_array=[]
    Selected_WeightedCentroid_array=[]
    for reg in Label_props:
        if reg.area>=size_min:
            if reg.area<=size_max:
                Selected_Areas_array.append(reg.area)
                Selected_IDs_array.append(reg.label)
                Selected_WeightedCentroid_array.append(reg.weighted_centroid)
    for ID in Selected_IDs_array:
        Focus_Locus=Label==ID
        Label_Size_Selected[Focus_Locus]=ID

    SignalMask=Mask>0 #outputs T or F 

    #Which labels from "Label" overlap with "Mask"?
    LabelMap_with_Mask=copy.copy(Label_Size_Selected)
    LabelMap_with_Mask[~SignalMask]=0 #makes all labels outside of Masks into "0"
    LabelIDsonMask=np.unique(LabelMap_with_Mask)#all label IDs overlapping with axis 
    LabelIDsonMask=LabelIDsonMask[LabelIDsonMask>0]#greater than 0, as this is background
    Overlapping_Labels=len(np.unique(LabelIDsonMask)[np.unique(LabelIDsonMask)>0])
    Total_Labels=len(np.unique(Label)[np.unique(Label)>0])
    NonOverlapping_Labels=Total_Labels-Overlapping_Labels

    #Make a KD tree for all Mask coordinates
    Mask_x,Mask_y=np.where(Mask>0)#splits the x and y coordinates
    Mask_xy=list(zip(Mask_x,Mask_y))#makes xy coords into a tuple for each pixel
    Mask_xy=np.array(Mask_xy)#convert to a numpy array for maths
    Mask_KD=scipy.spatial.KDTree(Mask_xy)

    # Calculate closest neighbour distances between Label and Mask
    Label_Distance_to_Mask=pd.DataFrame()
    #First calculate distances from weighted centroid of focus to Mask
    for order in range(len(Selected_WeightedCentroid_array)):
        Focus=Selected_WeightedCentroid_array[order]
        LabelID=Selected_IDs_array[order]
        d_cen,i_cen=Mask_KD.query(Focus,k=1)#provides distance from Focus to closest neighbour, and the index of this in the array
        #Second calculate distance from edge of the whole focus to mask
        Lab_x,Lab_y=np.where(Label_Size_Selected==LabelID)#select all label coordinates
        Distances=[]
        for pixel in range(len(Lab_x)):
            d_pixel,i_pixel=Mask_KD.query((Lab_x[pixel],Lab_y[pixel]),k=1)#provides distance from Focus to closest neighbour, and the index of this in the array
            Distances=np.append(Distances,d_pixel)
        Edge_Distance=min(Distances)#shortest distance between a pixel and the mask reflects the closest edge
        temp_df=pd.DataFrame([[LabelID,Focus,d_cen,Edge_Distance]],columns=["LabelID","Centroid_Coords","Centroid_to_Mask","Edge_to_Mask"])
        Label_Distance_to_Mask=pd.concat([Label_Distance_to_Mask,temp_df])


    return Label_Distance_to_Mask, Overlapping_Labels, NonOverlapping_Labels,Label_Size_Selected# self explanatory. Label_Size_Selected is the size-selected foci



#Distances from foci to mask and shuffle within nucleus with same measurements
def Whole_Nuclear_Shuffle(Nucleus_path,Axis_path,Focus_path,size_min,size_max, Shuffles): #includes size/shape restictions for foci. Nuclear area is >0 signal
    
    Nucleus=io.imread(Nucleus_path)
    Mask=io.imread(Axis_path)
    Label=io.imread(Focus_path)
    
    Label=Label.astype(int) # Make Label image values integer

    #######Analysis for the original foci
    Size=Label.shape#Make a new array to fill, which has same shape as other images

    #Select foci with areas within size range specified
    Label_Size_Selected=np.zeros((Size),dtype=int)
    Label_props=regionprops(Label, Label) #NB Intensity image can be used here to calculate weighted centroid coordinates
    Selected_Areas_array=[]
    Selected_IDs_array=[]
#    Selected_WeightedCentroid_array=[]
    Selected_Centroid_array=[]
    for reg in Label_props:
        if reg.area>=size_min:
            if reg.area<=size_max:
                Selected_Areas_array.append(reg.area)
                Selected_IDs_array.append(reg.label)
                Selected_Centroid_array.append(reg.centroid)
    for ID in Selected_IDs_array:
        Focus_Locus=Label==ID
        Label_Size_Selected[Focus_Locus]=ID
    SignalMask=Mask>0 #outputs T or F 

    #Which labels from "Label" overlap with "Mask"?
    LabelMap_with_Mask=copy.copy(Label_Size_Selected)
    LabelMap_with_Mask[~SignalMask]=0 #makes all labels outside of Masks into "0"
    LabelIDsonMask=np.unique(LabelMap_with_Mask)#all label IDs overlapping with axis 
    LabelIDsonMask=LabelIDsonMask[LabelIDsonMask>0]#greater than 0, as this is background
    Overlapping_Labels=len(np.unique(LabelIDsonMask)[np.unique(LabelIDsonMask)>0])
    Total_Labels=len(np.unique(Label)[np.unique(Label)>0])
    NonOverlapping_Labels=Total_Labels-Overlapping_Labels

    #Make a KD tree for all Mask coordinates
    Mask_x,Mask_y=np.where(Mask>0)#splits the x and y coordinates
    Mask_xy=list(zip(Mask_x,Mask_y))#makes xy coords into a tuple for each pixel
    Mask_xy=np.array(Mask_xy)#convert to a numpy array for maths
    Mask_KD=scipy.spatial.KDTree(Mask_xy)

    # Calculate closest neighbour distances between Label and Mask
    Label_Distance_to_Mask=pd.DataFrame()
    #First calculate distances from weighted centroid of focus to Mask
    for order in range(len(Selected_Centroid_array)):
        Focus=Selected_Centroid_array[order]
        LabelID=Selected_IDs_array[order]
        d_cen,i_cen=Mask_KD.query(Focus,k=1)#provides distance from Focus to closest neighbour, and the index of this in the array
        #Second calculate distance from edge of the whole focus to mask
        Lab_x,Lab_y=np.where(Label_Size_Selected==LabelID)#select all label coordinates
        Distances=[]
        for pixel in range(len(Lab_x)):
            d_pixel,i_pixel=Mask_KD.query((Lab_x[pixel],Lab_y[pixel]),k=1)#provides distance from Focus to closest neighbour, and the index of this in the array
            Distances=np.append(Distances,d_pixel)
        Edge_Distance=min(Distances)#shortest distance between a pixel and the mask reflects the closest edge
        temp_df=pd.DataFrame([[LabelID,Focus,d_cen,Edge_Distance]],columns=["LabelID","Centroid_Coords","Centroid_to_Mask","Edge_to_Mask"])
        Label_Distance_to_Mask=pd.concat([Label_Distance_to_Mask,temp_df])

    ########Shuffling and analysis for shuffled foci

    #Shuffle foci within nucleus       
    Nuclear_x,Nuclear_y=np.where(Nucleus>0)
    Nuclear_xy=list(zip(Nuclear_x,Nuclear_y))#makes xy coords into a tuple for each pixel
    Nuclear_Signal=Nucleus[Nuclear_xy[0]]#example of the signal contained in the nuclear mask for later reference
    NewLabel_Distance_to_Mask=pd.DataFrame()
    ShuffleOverlap=[]
    ShuffleNonOverlap=[]
    LabelShuffledRecordTotal=pd.DataFrame()
    for x in range(Shuffles):
        #Make a new array to fill, which has same shape as other images
        New_Label=np.zeros((Size),dtype=int)
        LabelShuffledRecordIndividual=pd.DataFrame()
        shuffle_number=x+1# for df as counter starts from 0
        for number in range(len(Selected_IDs_array)):#going from 0 
            lab=Selected_IDs_array[number]
            coords=np.where(Label_Size_Selected==lab)#coordinates of label
            x_cen,y_cen =Selected_Centroid_array[number]#unpack the tuple of the coordinates
            x_cen=int(x_cen)#make the values into integers from floats
            y_cen=int(y_cen)

            #calculate relative positions of the labels coordinates around the centroid
            rel_coords=(coords[0]-x_cen,coords[1]-y_cen)    
            #make a random selection from mask to use as the new centroid for a given focus
            New_Centroid=random.choice(Nuclear_xy)
            New_Focus_Coords=(rel_coords[0]+New_Centroid[0],rel_coords[1]+New_Centroid[1])#moves all the coordinates of the shape, so the centroid is the newly selected pixel from within the SCP3 domain

            limit_x, limit_y=Size
            limit_x=limit_x-1#Size starts from 1, but coordinates start from 0. Need to adjust
            limit_y=limit_y-1
            NewMax_x=max(New_Focus_Coords[0])
            NewMax_y=max(New_Focus_Coords[1])
            NewMin_x=min(New_Focus_Coords[0])
            NewMin_y=min(New_Focus_Coords[1])
            #Check new coords are within theimage frame
            condition2=limit_x<NewMax_x
            condition3=limit_y<NewMax_y
            condition4=0<NewMin_x
            condition5=0<NewMin_y
            #May need to trim the coordinates if they fall outside the frame to check for overlap and avoid encountering an error
            TrimmedCoords=Trim_Location(New_Focus_Coords,limit_x,limit_y)

            while condition2 or condition3 or condition4 or condition5 or sum(New_Label[TrimmedCoords])!=0 or statistics.mean(Nucleus[TrimmedCoords])!=Nuclear_Signal:#final two conditions avoid new foci overlapping with eachother, and avoid edges falling outside the bounds of Nucleus
                New_Centroid=random.choice(Nuclear_xy)
                New_Focus_Coords=(rel_coords[0]+New_Centroid[0],rel_coords[1]+New_Centroid[1])#moves all the coordinates of the shape, so the centroid is the newly selected pixel from within the SCP3 domain
                NewMax_x=max(New_Focus_Coords[0])
                NewMax_y=max(New_Focus_Coords[1])
                NewMin_x=min(New_Focus_Coords[0])
                NewMin_y=min(New_Focus_Coords[1])
                TrimmedCoords=Trim_Location(New_Focus_Coords,limit_x,limit_y)
                condition2=limit_x<NewMax_x
                condition3=limit_y<NewMax_y
                condition4=0>NewMin_x
                condition5=0>NewMin_y    

            New_Label[New_Focus_Coords]=lab#give these pixels the label value
            #Calculate new weighted centroid 
            x_cen, y_cen=Selected_Centroid_array[number]#also get the weighted centroid coordinates to mark this point for future measurements relative to the centroid by shape
            rel_cen_x=int(x_cen-x_cen)
            rel_cen_y=int(y_cen-y_cen)
            New_cen=(rel_cen_x+New_Centroid[0],rel_cen_y+New_Centroid[1])

            df1=pd.DataFrame([[New_cen,lab,number,shuffle_number]],columns=["Coords","Label","Focus_Number","Shuffle_Number"])#new df with details of coords, label, focus numer and shuffle number
            LabelShuffledRecordTotal=pd.concat([LabelShuffledRecordTotal,df1])
            LabelShuffledRecordIndividual=pd.concat([LabelShuffledRecordIndividual,df1])
        
        #Which labels from "New_Label" overlap with "Mask"?
        LabelMap_with_Mask=copy.copy(New_Label)
        LabelMap_with_Mask[~SignalMask]=0 #makes all labels outside of Masks into "0"
        LabelIDsonMask=np.unique(LabelMap_with_Mask)#all label IDs overlapping with axis 
        LabelIDsonMask=LabelIDsonMask[LabelIDsonMask>0]#greater than 0, as this is background
        Overlapping_ShuffleLabels=len(np.unique(LabelIDsonMask)[np.unique(LabelIDsonMask)>0])
        Total_Labels=len(np.unique(Label)[np.unique(Label)>0])
        NonOverlapping_ShuffleLabels=Total_Labels-Overlapping_ShuffleLabels
        ShuffleOverlap=np.append(ShuffleOverlap,Overlapping_ShuffleLabels)
        ShuffleNonOverlap=np.append(ShuffleNonOverlap,NonOverlapping_ShuffleLabels)

        # Calculate closest neighbour distances between Label and Mask
        #First calculate distances from centroid of focus to Mask
        NewLabel_props=regionprops(New_Label, New_Label)#get properties of each label
        NewLabel_Centroid_array=[]
        NewLabel_Label_array=[]
        for reg in NewLabel_props:
            NewLabel_Label_array.append(reg.label)
            NewLabel_Centroid_array.append(reg.centroid)

        for order in range(len(NewLabel_Centroid_array)):
            Focus=NewLabel_Centroid_array[order]
            LabelID=NewLabel_Label_array[order]
            d_cen,i_cen=Mask_KD.query(Focus,k=1)#provides distance from Focus to closest neighbour, and the index of this in the array
            #Second calculate distance from edge of the whole focus to mask
            Lab_x,Lab_y=np.where(New_Label==LabelID)#select all label coordinates
            Distances=[]
            for pixel in range(len(Lab_x)):
                d_pixel,i_pixel=Mask_KD.query((Lab_x[pixel],Lab_y[pixel]),k=1)#provides distance from Focus to closest neighbour, and the index of this in the array
                Distances=np.append(Distances,d_pixel)
            Edge_Distance=min(Distances)#shortest distance between a pixel and the mask reflects the closest edge
            temp_df=pd.DataFrame([[LabelID,Focus,d_cen,Edge_Distance,x]],columns=["LabelID","Centroid_Coords","Centroid_to_Mask","Edge_to_Mask","Shuffle"])
            NewLabel_Distance_to_Mask=pd.concat([NewLabel_Distance_to_Mask,temp_df])
    
    #Stats
    real_from_Centroid=Label_Distance_to_Mask["Centroid_to_Mask"]
    real_from_Edge=Label_Distance_to_Mask["Edge_to_Mask"]
    fake_from_Centroid=NewLabel_Distance_to_Mask["Centroid_to_Mask"]
    fake_from_Edge=NewLabel_Distance_to_Mask["Edge_to_Mask"]
    ks_stat_Centroid,ks_p_Centroid=ks_2samp(real_from_Centroid,fake_from_Centroid) 
    ks_stat_Edge,ks_p_Edge=ks_2samp(real_from_Edge,fake_from_Edge) 
    fisher_oddsratio,fisher_p=fisher_exact([[Overlapping_Labels,NonOverlapping_Labels],[sum(ShuffleOverlap),sum(ShuffleNonOverlap)]])

    return Label_Distance_to_Mask, NewLabel_Distance_to_Mask, Label_Size_Selected, New_Label  # Outputs: observed distance from foci to the axial mask, shuffled distances to the axial mask, labelmap image of the size-filtere foci included in the analysis, an example labelmap image of one of the rounds of positional shuffling