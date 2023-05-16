# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 21:16:07 2023

@author: jcc234
"""

def Axis_Proximity(trace_path):
    
    import numpy as np
    import seaborn as sns
    sns.set_style("white")
    import re 
    import scipy   
    import imp
    import TracePy as tp
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