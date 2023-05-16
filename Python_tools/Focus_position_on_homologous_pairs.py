# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 16:16:33 2023

@author: jcc234
"""


def Resolved_homologues_focus_position(Trace_df, Foci_Labelmap, Foci_Intensity_Img, Min_Area, Max_Area):

    import pandas as pd
    import scipy
    from skimage.measure import regionprops
    from math import sqrt
    

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


    
