import os
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style("white")
import scipy
# os.chdir(TracePy_path)
# import TracePy as tp
import TracePy as tp
import imp
imp.reload(tp)


def Focus_Position_on_Trace(trace_path, label_data):
    
   
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
