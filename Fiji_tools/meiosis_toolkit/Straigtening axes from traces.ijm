""" Generating straightened axis images for each trace """

//Choose a directory and loop through the files
//If a tif image is found and a .traces image with the same name, then these traces will be converted to ROIs and used to straighten the axes
//Data will be saved to new metadata folder


source_dir=getDir("Master Directory");//dir containing images and traces files
list=getFileList(source_dir);

setBatchMode("hide");

for (i=0; i<list.length; i++){
	if (endsWith(list[i],".tif")){
		img=i;
		trace_name=replace(img,".tif",".traces");
		for (j=0; j<list.length; j++){
			if (list[j]==trace_name)){
				//open image
				run("Bio-Formats Importer", "open=["+source_dir+img+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
				rename("img");
				
				//loop through roi files
				roi_list=getFileList(source_dir+folder);
				for (k=0; k<roi_list.length; k++){
					if (endsWith(roi_list[k], "roi")){
						print(roi_list[k]);
						roiManager("reset");	
						roiManager("Open", source_dir+folder+roi_list[k]);
						file_name=replace(roi_list[k],".roi","_"+name+"_straight");
						selectImage("img");roiManager("Select", 0);
						run("Straighten...", "title=["+img+"] line=21 process");saveAs("tiff", source_dir+folder+file_name);
						close(file_name);
						}}
				
}}}
close("img");
}


//run("script:Skeletons_and_ROIs/Convert_Reconstruction_To_ROIs.groovy", "recfile=[C:\\Users\\jcc234\\OneDrive - University of Exeter\\Desktop\\Code for GitHub\\Meiotic-Image-Analysis-Toolkit\\SoRa synapsed pachytene data\\100xB6_260918A_D_MLH1_Scp3_Syce2 - trace.traces] impfile=[C:\\Users\\jcc234\\OneDrive - University of Exeter\\Desktop\\Code for GitHub\\Meiotic-Image-Analysis-Toolkit\\SoRa synapsed pachytene data\\100xB6_260918A_D_MLH1_Scp3_Syce2 - Denoise_ai.nd2_MaxProj.tif-Registered.tif] convertpaths=true convertbranchpoints=false converttips=false");


