import numpy as np
import seaborn as sns; sns.set()
import scipy.stats
import scipy.ndimage
import pandas as pd
import copy
import random
import scipy.spatial 
import statistics
import seaborn as sns
sns.set_style("white")
import scipy
import re 
from skimage import io
from skimage.measure import regionprops
from math import sqrt


#During the shuffling QC, need to check if new focus overlaps old other new ones. If they don't fall completely within the image frame this prevents identifying overlap. This function trims the focus to fit the frame, ensuring the shuffling functions continue to run. NB shuffled focus positions will not be included if they extend beyond the image frame. 
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
    

    return Label_Distance_to_Mask, NewLabel_Distance_to_Mask, Label_Size_Selected, New_Label  # Outputs: observed distance from foci to the axial mask, shuffled distances to the axial mask, labelmap image of the size-filtere foci included in the analysis, an example labelmap image of one of the rounds of positional shuffling

def Axial_Shuffle(Mask,Label1,Label2,Intensity1,Intensity2,Label1_size_min,Label1_size_max,Label2_size_min,Label2_size_max, Shuffles):#Mask (can be labelmap), labels to be compared to, labels to be used for shuffling, intensity images for either set of labels, max and min areas for focus selection for each set of foci, number of shuffles
     
    Label2ShuffledRecordTotal=pd.DataFrame(columns=["Coords","Label","Focus_Number","Shuffle_Number"])
    Size=Label2.shape

    Mask_Map=np.zeros((Size))
    SignalMask=Mask>0 
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
    #Label2 size selection
    Label2_Mask_Full=np.zeros((Size),dtype=int)#First select foci overlapping the axis mask
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
    
    #Label1 size selection
    Label1_Mask_Full_selected=np.zeros((Size),dtype=int)
    Starting_Label1_props=regionprops(Label1)
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
    
      
    #Extract the measurements for each label 
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
    KDTreeShuffledLabel2toLabel1=pd.DataFrame()
    for x in range(Shuffles):
        #Make a new array to fill, which has same shape as other images
        New_Label2=np.zeros((Size),dtype=int)
        
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
            New_Focus_Coords=(rel_coords[0]+New_Centroid[0],rel_coords[1]+New_Centroid[1])#moves all the coordinates of the shape, so the centroid is the newly selected pixel from within the axial domain

            limit_x, limit_y=Size
            limit_x=limit_x-1#Size starts from 1, but coordinates start from 0. Need to adjust
            limit_y=limit_y-1
            NewMax_x=max(New_Focus_Coords[0])
            NewMax_y=max(New_Focus_Coords[1])
            NewMin_x=min(New_Focus_Coords[0])
            NewMin_y=min(New_Focus_Coords[1])
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
            
            
    #Caluclate the nearest neighbours for all the shuffled Label2 to Label1
    for index, row in Label2ShuffledRecordTotal.iterrows():       
        Focus=row["Coords"]
        Label=row["Label"]
        Shuffle_Number=row["Shuffle_Number"]+1
        d,i=Label1_Tree.query(Focus,k=1)#provides distance from "test" to closest neighbour, and the index of this in the array
        Neighbour=weighted_cen_array_1[i]
        Neighbour_Value=lab_array_1[i]
        KD_temp_df=pd.DataFrame([[Label,Focus,d,Neighbour,Neighbour_Value, Shuffle_Number]],columns=["Focus_Label", "Label_Location","Distance","Neighbour_Location","Neighbour_Label", "Shuffle_Number"])
        KDTreeShuffledLabel2toLabel1=pd.concat([KDTreeShuffledLabel2toLabel1,KD_temp_df])        

    #Caluclate the nearest neighbours for all the original (unshuffled, but included in the shuffling) Label2 to Label1 foci
    original_cen_distances=pd.DataFrame()
    for row in range(len(cen_array_2)):
        Focus=cen_array_2[row]
        Label=Selected_IDs_array_label2[row]
        d,i=Label1_Tree.query(Focus,k=1)#provides distance from "test" to closest neighbour, and the index of this in the array
        Neighbour=weighted_cen_array_1[i]
        Neighbour_Value=lab_array_1[i]
        original_cen_distances_temp=pd.DataFrame([[Label,Focus,d,Neighbour,Neighbour_Value]],columns=["Focus_Label", "Label_Location","Distance","Neighbour_Location","Neighbour_Label"])
        original_cen_distances=pd.concat([original_cen_distances,original_cen_distances_temp])        
    
    #Caluclate the nearest neighbours for all the original (unshuffled, but included in the shuffling) Label2 to Label1 foci. NB this time using Label2 centroids weighted by signal instensity
    original_weighted_cen_distances=pd.DataFrame()
    for row in range(len(weighted_cen_array_2)):
        Focus=weighted_cen_array_2[row]
        Label=Selected_IDs_array_label2[row]
        d,i=Label1_Tree.query(Focus,k=1)#provides distance from "test" to closest neighbour, and the index of this in the array
        Neighbour=weighted_cen_array_1[i]
        Neighbour_Value=lab_array_1[i]
        original_weighted_cen_distances_temp=pd.DataFrame([[Label,Focus,d,Neighbour,Neighbour_Value]],columns=["Focus_Label", "Label_Location","Distance","Neighbour_Location","Neighbour_Label"])
        original_weighted_cen_distances=pd.concat([original_weighted_cen_distances,original_weighted_cen_distances_temp])        

        

    return  KDTreeShuffledLabel2toLabel1, original_cen_distances, original_weighted_cen_distances, Label1_Mask_Full, Label2_Mask_Full, New_Label2 # shuffled distances between foci in Label2 and Label1; size_selected foci from Label1, size selected foci from Label2, example of axially shuffed labelmap
 
def Focus_Position_on_Trace(trace_path, label_data):
    import TracePy as tp
    import imp
    imp.reload(tp)
 
   
    df, meta = tp.read_trace(trace_path)

    ########## Adding data to the df:
    #1. First add the pixel number (total, individual and percentage) to each trace path
    all_trace_names=df["name"].unique().tolist()
    all_trace_ids=df["id"].unique().tolist()
    all_paths=tp.get_all_paths(df)

    for ID in range(0,len(all_trace_names)):#loop through each path ID
        single_trace_name=all_trace_names[ID]
        single_trace_id=int(all_trace_ids[ID])
        single_trace_data=df.loc[df["name"]==single_trace_name] #select data for single trace
        single_trace_index=df.index[df.name==single_trace_name]#get index positions in df for each row

        starting_pixel=single_trace_index[0]#for px position calculations
        final_pixel=single_trace_index[-1]
        total_pixels=len(single_trace_index)
        trace = all_paths[single_trace_id]#for euclidian distance calculations
        length,dims=trace.shape
        full_length = np.sum(tp.get_distance(trace[:][1:], trace[:][:-1]))

        #Adding pixel number info
        for point_row in single_trace_index:
            df.at[point_row,"Total_Trace_Pixels"]=total_pixels
            pixel_position=point_row-starting_pixel+1
            df.at[point_row,"Trace_Pixel_Number"]=pixel_position
            df.at[point_row,"Trace_Pixel_Percentage"]=pixel_position/total_pixels*100
        #Adding distance measurements
            distance_on_path = np.sum(tp.get_distance(trace[:pixel_position][1:], trace[:pixel_position][:-1]))
            percentage_distance=distance_on_path/full_length*100
            df.at[point_row,"Distance_euclidian"]=distance_on_path
            df.at[point_row,"Percentage_distance_euclidian"]=percentage_distance
            df.at[point_row,"Full_axis_length_euclidian"]=full_length
            
 
   #Calculate a KD tree for the coordinates of all the traces combined
    #2. Get unique coordinates covered by the traces
    all_coords=[]
    for trace in all_paths:
        trace_coords=all_paths[trace]
        for individual_coord in trace_coords:
            all_coords.append(individual_coord)
    all_coords_tree=scipy.spatial.KDTree(all_coords)# make kd tree

    #3. Calculate labelled focus position relative to the traced axes         
    Label_Info=pd.read_csv(label_data)
        
    for lab in range(0,len(Label_Info)):
        cen_x=Label_Info.at[lab,'centroid_x']
        cen_y=Label_Info.at[lab,'centroid_y']
        
        #Distance from axis trace
        cen_3d=cen_x,cen_y,0 #need to make 3dimensional coordinates to compare to the data from the axis trace
        d,i=all_coords_tree.query(cen_3d,k=1)
        Label_Info.at[lab,"Centroid_Trace_Distance"]=d

        closest_axis_pixel=all_coords[i]
        x,y,z=closest_axis_pixel
        Label_Info.at[lab,"Close_Trace_coord_x"]=x
        Label_Info.at[lab,"Close_Trace_coord_y"]=y
        axis_pixel_info=df.loc[(df["x"]==str(x))&(df["y"]==str(y))]#values in the df are strings, so match in this format
        axis_pixel_index=df.index[(df["x"]==str(x))&(df["y"]==str(y))]#values in the df are strings, so match in this format
        Label_Info.at[lab,"Percentage_Distance_Along_Axis_Euclidean"]=axis_pixel_info.at[axis_pixel_index[0],"Percentage_distance_euclidian"]
        Label_Info.at[lab,"Euclidian_Distance_Along_Trace"]=axis_pixel_info.at[axis_pixel_index[0],"Distance_euclidian"]
        Label_Info.at[lab,"Axis_Full_Length_Euclidian"]=axis_pixel_info.at[axis_pixel_index[0],"Full_axis_length_euclidian"]
        Label_Info.at[lab,"Total_Trace_Pixels"]=axis_pixel_info.at[axis_pixel_index[0],"Total_Trace_Pixels"]
        Label_Info.at[lab,"Px_Number_Along_Axis_Trace"]=axis_pixel_info.at[axis_pixel_index[0],"Trace_Pixel_Number"]
        Label_Info.at[lab,"Trace_Name"]=axis_pixel_info.at[axis_pixel_index[0],"name"]

    return df,Label_Info

def Resolved_homologues_focus_position(Trace_df, Foci_Labelmap, Foci_Intensity_Img, Min_Area, Max_Area):

    

    #Calculate a KD tree for the coordinates of all the traces combined
    #Get unique coordinates covered by the traces
    all_coords=[]
    for index, row in Trace_df.iterrows():
        individual_coord=row["x"],row["y"],row["z"]
        all_coords.append(individual_coord)

    all_coords_tree=scipy.spatial.KDTree(all_coords)

    #Extract relevant label info to the Foci_Data dataframe
    Foci_props=regionprops(Foci_Labelmap,Foci_Intensity_Img)#get properties of each label
    Foci_data=pd.DataFrame()#collect necessary foci data
        
    for lab in range(0,len(Foci_props)):
        reg=Foci_props[lab]
        Foci_data.at[lab,"ID"]=reg.label
        Foci_data.at[lab,"Area"]=reg.area
        cen_x,cen_y=reg.centroid
        Foci_data.at[lab,"Centroid_x"]=cen_x#coordinates have to be unpacked as they're tuples
        Foci_data.at[lab,"Centroid_y"]=cen_y
        w_cen_x,w_cen_y=reg.weighted_centroid
        Foci_data.at[lab,"Weighted_Centroid_x"]=w_cen_x
        Foci_data.at[lab,"Weighted_Centroid_y"]=w_cen_y
        Foci_data.at[lab,"Mean_Intensity"]=reg.mean_intensity
        #Distance from axis trace
        weighted_cen_3d=w_cen_x,w_cen_y,0 #need to make 3dimensional coordinates to compare to the data from the axis trace
        d,i=all_coords_tree.query(weighted_cen_3d,k=1)
        Foci_data.at[lab,"Weighted_Centroid_Axis_Distance"]=d
        
        closest_axis_pixel=all_coords[i]
        x,y,z=closest_axis_pixel
        Foci_data.at[lab,"Close_Axis_coord_x"]=x
        Foci_data.at[lab,"Close_Axis_coord_y"]=y
       
        #Pull out the row of data relating to this trace pixel
        axis_pixel_info=Trace_df.loc[(Trace_df["x"]==x)&(Trace_df["y"]==y)]#find the info matching the trace pixel cooredinates
        axis_pixel_index=Trace_df.index[(Trace_df["x"]==x)&(Trace_df["y"]==y)]
        
        Close_axis_name=axis_pixel_info.at[axis_pixel_index[0],"name"]
        Foci_data.at[lab,"Close_axis_name"]=Close_axis_name
        
        #Distance from the distal partner in the homololog pair
        #Some traces won't have a partner, so only include those which do
        if (axis_pixel_info.at[axis_pixel_index[0],"Trace Partner"]!="No homolog partner found") :
            n_x=axis_pixel_info.at[axis_pixel_index[0],"Neighbour px x"]#selects only the neighbour of the trace occupying the given pixel at the top of the df
            n_y=axis_pixel_info.at[axis_pixel_index[0],"Neighbour px y"]

            Foci_data.at[lab,"Close_axis_homologue_partner_name"]=axis_pixel_info.at[axis_pixel_index[0],"Trace Partner"]
            Foci_data.at[lab,"Distance_between_proximal_axis_and_partner"]=axis_pixel_info.at[axis_pixel_index[0],"Distance from homolog"]
            Foci_data.at[lab,"Euclidian_percentage_along_close_axis"]=axis_pixel_info.at[axis_pixel_index[0],"Percentage_distance_euclidian"]
            Foci_data.at[lab,"Total_Euclidian_length_close_axis"]=max(Trace_df["Distance_euclidian"].loc[(Trace_df["name"]==Close_axis_name)])
            Foci_data.at[lab,"Euclidian_distace_along_axis"]=(Foci_data.at[lab,"Euclidian_percentage_along_close_axis"]*Foci_data.at[lab,"Total_Euclidian_length_close_axis"])/100        
            
            Partner_axis_name=axis_pixel_info.at[axis_pixel_index[0],"Trace Partner"]
            partner_axis_pixel_info=Trace_df.loc[(Trace_df["x"]==n_x)&(Trace_df["y"]==n_y)&(Trace_df["name"]==Partner_axis_name)]
            partner_axis_pixel_index=Trace_df.index[(Trace_df["x"]==n_x)&(Trace_df["y"]==n_y)&(Trace_df["name"]==Partner_axis_name)]
            
            Foci_data.at[lab,"Euclidian_percentage_along_partner_axis"]=partner_axis_pixel_info.at[partner_axis_pixel_index[0],"Percentage_distance_euclidian"]
            Foci_data.at[lab,"Total_Euclidian_length_partner_axis"]=max(Trace_df["Distance_euclidian"].loc[(Trace_df["name"]==Partner_axis_name)])
            Foci_data.at[lab,"Euclidian_distace_along_partner_axis"]=(Foci_data.at[lab,"Euclidian_percentage_along_partner_axis"]*Foci_data.at[lab,"Total_Euclidian_length_partner_axis"])/100    
            Foci_data.at[lab,"Euclidian_distace_along_partner_axis"]
                       
            #Is the focus between or outside axes?
            #If the distance to the partner axis is greater than the distance from the close axis to its partner, then the focus must be outside the axes, if the distance is shorter then it is between the axes
            Foci_data.at[lab,"Weighted_Centroid_Partner_Axis_Distance"]=sqrt(((w_cen_x-n_x)**2)+((w_cen_y-n_y)**2))
            Foci_data.at[lab,"Focus_Between_Axes"]=Foci_data.at[lab,"Weighted_Centroid_Partner_Axis_Distance"]<Foci_data.at[lab,"Distance_between_proximal_axis_and_partner"]
            
            
    return Foci_data

def Axis_Proximity(trace_path): 
    import TracePy as tp
    import imp
    imp.reload(tp)
        
    df, meta = tp.read_trace(trace_path)
    
    # Swap x and y axis in trace data (and xd/yd too to avoid confusion) as these dimensions switch from ImageJ
    df.rename(columns={"x":"y","y":"x","xd":"yd","yd":"xd"},inplace=True)
    
    
    ########## Adding data to the df:
    
    #1. First add the pixel number (total, individual and percentage) and Euclidean distances to each trace path position
    all_trace_names=df["name"].unique().tolist()
    all_trace_ids=df["id"].unique().tolist()
    all_paths=tp.get_all_paths(df)
    
    for ID in range(0,len(all_trace_names)):#loop through each path ID
        single_trace_name=all_trace_names[ID]
        single_trace_id=int(all_trace_ids[ID])#this position needs to be defined, as if any traces are deleted then their original position will cease to exist too, to numerical gaps can exist in the library
        single_trace_data=df.loc[df["name"]==single_trace_name] #select data for single trace
        single_trace_index=df.index[df.name==single_trace_name]#get index positions in df for each row
    
        starting_pixel=single_trace_index[0]#for px position calculations
        total_pixels=len(single_trace_index)
        trace = all_paths[single_trace_id]#for euclidian distance calculations
        length,dims=trace.shape
        full_length = np.sum(tp.get_distance(trace[:][1:], trace[:][:-1]))
    
        #Adding pixel number info
        for point_row in single_trace_index:
            df.at[point_row,"Total_Trace_Pixels"]=total_pixels
            pixel_position=point_row-starting_pixel+1
            df.at[point_row,"Trace_Pixel_Number"]=pixel_position
            df.at[point_row,"Trace_Pixel_Percentage"]=pixel_position/total_pixels*100
        #Adding distance measurements
            distance_on_path = np.sum(tp.get_distance(trace[:pixel_position][1:], trace[:pixel_position][:-1]))
            percentage_distance=distance_on_path/full_length*100
            df.at[point_row,"Distance_euclidian"]=distance_on_path
            df.at[point_row,"Percentage_distance_euclidian"]=percentage_distance
            #print(ID,point_row,pixel_position,distance_on_path,percentage_distance)
    
    #2. Next add the homolog partner ID of each trace (if present and noted as "number""letter", and only has one partner)
    #use reg expressiongs to break up chromosome names and work out what has a partner
    for ID in all_trace_names:#loop through each path ID
        homolog_partner_confirmed=[]#initiate list
        if re.match("([0-9]+)([a-zA-Z]+)", ID): #if a number then a letter then this trace should have a partner (only one)
            temp = re.compile("([0-9]+)([a-zA-Z]+)") 
            res = temp.match(ID).groups() 
            chromosome_trace,homolog=res  
    
            #Find partner chromosome ID
            for ID_Partner in all_trace_names:
                if re.match("([0-9]+)([a-zA-Z]+)", ID_Partner): #ID number and letter?
                    chromosome_partner,homolog_partner = temp.match(ID_Partner).groups() #split into number and letter
                    if chromosome_partner==chromosome_trace:#must be the same chromosome number
                        if homolog!=homolog_partner:#must be a different letter after this number
                            homolog_partner_confirmed.append(ID_Partner)
    
        else:
            Partner="No homolog partner found"
        if len(homolog_partner_confirmed)==1:
            Partner=homolog_partner_confirmed[0]
        if len(homolog_partner_confirmed)>1:
            Partner= "Error. Multiple homolog partners found"
    
        df.loc[df["name"]==ID,"Trace Partner"]=Partner
    
    #3. Calculate distances between axes points
    for ID in range(0,len(all_trace_ids)):#loop through each path ID
        single_trace_name=all_trace_names[ID]
        single_trace_id=int(all_trace_ids[ID])#this position needs to be defined, as if any traces are deleted then their original position will cease to exist too, to numerical gaps can exist in the library
    
        single_trace_data=df.loc[df["name"]==single_trace_name] #select data for single trace
        single_trace_index=df.index[df.name==single_trace_name]#get index positions in df for each row
    
        trace = all_paths[single_trace_id]#for euclidian distance calculations  
        Partner_ID=single_trace_data["Trace Partner"].unique()[0]
        #print(single_trace_name,Partner_ID)
    
        if re.match("([0-9]+)([a-zA-Z]+)", Partner_ID):
            #Partner_ID_index=all_trace_names.index(Partner_ID)
            partner_trace_data=df.loc[df["name"]==Partner_ID] #select data for single trace
            Partner_ID_index=int(partner_trace_data["id"].unique()[0])
    
            Partner_trace = all_paths[Partner_ID_index]#for euclidian distance calculations
            #print(single_trace_ID, Partner_ID)
            Traced_distance_tree=scipy.spatial.KDTree(Partner_trace)#Make a KD tree to order the positions of each pixel in a given trace
    
            for single_point in range(0,len(single_trace_index)):
                query_coords=trace[single_point]
                row_index=single_trace_index[single_point]
                d,i=Traced_distance_tree.query(query_coords,k=1)#provides distance from "test" to closest neighbour, and the index of this in the array
                x,y,z=Partner_trace[i]
                df.at[row_index,"Distance from homolog"]=d
                df.at[row_index,"Neighbour px x"]=int(x)
                df.at[row_index,"Neighbour px y"]=int(y)
                df.at[row_index,"Neighbour px z"]=int(z)

    return df

