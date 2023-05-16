

setBatchMode("False");
run("Set Scale...", "distance=0 known=0 unit=pixel global");
//Prepare the environment
roiManager("reset");//clear ROIs from manager
setOption("BlackBackground", true);//setting the background for binary images. 
run("Colors...", "foreground=white background=black selection=darkgray");
close("Log");//ensure log is clear
run("Clear Results");

//Ensure appropriate image is open
if (nImages!=1){
exit("Before running macro, open a single stack image to analyse");
}

if (nSlices<2){
exit("Stack image must have at least two slices");
}


title = getTitle();
title= replace(title, "\.tif", "");

////First assign the channels "Axis" and "Foci"
Dialog.createNonBlocking("Select channels to segment");

slice_range=newArray();
for (i = 1; i <= nSlices; i++) {
  slice_range=Array.concat(slice_range,d2s(i,0));
}
  slice_range_axes=Array.concat(slice_range,"Pre-existing mask");
  slice_range_foci=Array.concat(slice_range,"Pre-existing labelmap");
  
Dialog.addChoice("Axis channel number:", slice_range_axes);
//Dialog.addString("Pre-segmented mask filename", "Axes_Mask");
Dialog.addString("Axis analysis title", "Axes_Mask");
Dialog.addChoice("Foci channel number:", slice_range_foci);
//Dialog.addString("Pre-segmented mask filename", "Watershed_Foci");
Dialog.addString("Foci analysis title", "Watershed_Foci");
Dialog.show();

Axes_channel  = Dialog.getChoice();
//Axes_filename  = Dialog.getString();
Axes_output_name  = Dialog.getString();
Foci_channel  = Dialog.getChoice();
//Foci_filename  = Dialog.getString();
Foci_output_name  = Dialog.getString();

//Choose directory to save folder of output to
dir1= getDirectory("Choose destination directory");
dir2 = dir1 + title +"_Output";
width=getWidth();
height=getHeight();
run("Select None");run("Duplicate...", "title=Original duplicate");

//Does a metadata folder already exist? If so open the Image_Summary.csv and RoiSet.zip
files=getFileList(dir1);
Image_not_previously_scored=1;//1 == True
potential_metadata_path=title+"_Output/";
for (i = 0; i < files.length; i++) {
    if (files[i]==potential_metadata_path){
    	Image_not_previously_scored=0;
		open(dir2+"/Image_Summary.csv");	
		XY_value=getResultString("XY Analysis", 0);
		XY=(XY_value=="Yes");//convert string to bool
		Manual_Stage = getResultString("Stage", 0);
		roiManager("Open",dir2+"/RoiSet.zip");	
		print("Existing image data detected");
    }}


//1. Input basic image details for reference
if (Image_not_previously_scored){
	//make new folder and save files to it
	File.makeDirectory(dir2);
	
	run("Brightness/Contrast...");run("Channels Tool...");
	
	Dialog.createNonBlocking("Image Detials");
	MeioticStage=newArray("Leptotene", "Zygotene","Pachytene","Diplotene","Metaphase I","Asynapsed Pachytene","Other");
	Dialog.addRadioButtonGroup("Meiotic Stage", MeioticStage, 2,3,"");
	Dialog.addCheckbox("Separate XY Scoring", false);
	Dialog.addString("Additional Comments","");
	Dialog.show;
		
	Manual_Stage=Dialog.getRadioButton();
	Comments=Dialog.getString();
	XY=Dialog.getCheckbox();
	close("Channels");
}

//2. Select ROIs

//2a. Manually select nuclear ROI
if (Image_not_previously_scored){
	selectWindow("Original");
	setTool("freehand");
	run("Colors...", "foreground=white background=black selection=yellow");
	waitForUser("Draw around nucleus");
	roiManager("Add");roiManager("Select", 0);roiManager("Rename", "Nucleus_Mask");run("Select None");
	run("Colors...", "foreground=white background=black selection=darkgray");

	newImage("Nuclear_Mask", "8-bit black", width, height, 1);//make nuclear mask image
	roiManager("Select", 0);roiManager("Fill");run("Select None");saveAs("tiff", dir2 + "/Nuclear_Mask");
	close("Nuclear_Mask");
	}


//2b. Axial ROI
reusing_axis_mask=(Axes_channel=="Pre-existing mask");
if (reusing_axis_mask==0){
	selectWindow("Original");run("Duplicate...", "title=Axes duplicate channels="+Axes_channel);
	selectWindow("Axes");run("Duplicate...", "title="+Axes_output_name);
	run("Brightness/Contrast...");
	selectWindow(Axes_output_name);roiManager("Select", 0);run("Clear Outside");run("Select None");//clear signal outside nucleus
	waitForUser("Adjust Axis Brightness/Contrast", "When ready, click \"OK\". Do not click \"Apply\""); //Need to press OK to continue
	
	///Gaussian Filtering of Axis
	//Place holders to set loop
	sigma=0;
	decision="";
	filtered_img_name="filtered_img_name";
	//Select level of Gaussian filtering for axis segmentation
	while (decision!="Proceed"){
		Dialog.createNonBlocking("Gaussian Filter");
		Dialog.addMessage("Set a Gaussian filter to reduce noise in image.This can aid axis segmentation.\nClick \"OK\" to try thresholding. Adjustment is permitted later if necessary");
		Dialog.addString("Radius (sigma)", sigma);
		Dialog.show;
		
		sigma=Dialog.getString();

		filtered_img_name="Preview Sigma = "+sigma;
		selectWindow(Axes_output_name);
		run("Duplicate...", "title=["+filtered_img_name+"]");
		run("Gaussian Blur...", "sigma="+sigma);

		run("Threshold...");call("ij.plugin.frame.ThresholdAdjuster.setMode", "Red");

		Dialog.createNonBlocking("Thresholding Filtered Image");
		Dialog.addMessage("Adjust threshold limits, making axis territory red. \nDo not \"Apply\", but select \"Proceed\" and click \"OK\" if you're happy with the result.\nAlternatively select \"Reset\" to restart the filtering process.");
		filter_options = newArray("Reset", "Proceed");
		Dialog.addRadioButtonGroup("", filter_options, 2, 1, "Reset");
		Dialog.show;
		
		decision=Dialog.getRadioButton();		
					
		if (decision=="Reset"){
			close(filtered_img_name);
			close("Threshold");
			}
		if(decision=="Proceed"){
			run("Convert to Mask");close("Threshold");
			close(Axes_output_name);//close the axial intensity image
			selectWindow(filtered_img_name);rename(Axes_output_name);//Assign the mask this name as it's used later
			run("Create Selection");roiManager("Add");run("Select None");
			Axis_roi=roiManager("count")-1;
			roiManager("Select", Axis_roi);roiManager("Rename", Axes_output_name);run("Select None");
		}}

	//If XY analysis has been selected then separate data needs to be generated for an axis mask lasking the sex chromosomes
	//Removing XY manually using eraser if requested at start 
	if (XY) {
		selectWindow("Axes");run("Duplicate...","title=Axes_Copy");run("8-bit");
		selectWindow(Axes_output_name);run("Select None");run("Duplicate...","title="+Axes_output_name+"_Copy");
		run("Concatenate...", "  title=[Remove XY] image1=["+Axes_output_name+"_Copy] image2=[Axes_Copy] image3=[-- None --]");
		run("Make Composite", "display=Composite");
		
		Stack.setChannel(1);run("Blue");
		Stack.setChannel(2);run("Grays");
		setSlice(1);run("Channels Tool...");
		
		waitForUser("Check Mask","Erase sex chromosomes from the mask (Blue) using \"Eraser Tool\" in \"Drawing Tools\".\nOriginal axis channel (Grey) is a reference and can be turned on/off using the \"Channels Tool\"\nClick \"OK\" when complete.");//users can remove signal with eraser
		selectWindow("Channels");run("Close");
		Autosome_output_name="Autosome_"+Axes_output_name;
		selectWindow("Remove XY");run("Duplicate...", "title="+Autosome_output_name+" channels=1");run("Grays");
		run("Create Selection");roiManager("Add");run("Select None");//add roi of axes withot XY to manager
		Autosome_roi=roiManager("count")-1;
		roiManager("Select", Autosome_roi);roiManager("Rename", "Autosome_"+Axes_output_name);
		}}
	
else {//if reusing existing roi for axes then define whih one. Create dialogue with one opetion if only autosomal, two if XY are to be removed
	waitForUser("Select Axial Roi in Manager");
	Axes_output_name=Roi.getName;
	Axis_roi=roiManager("index");
	//Make mask for axis using the roi
	newImage(Axes_output_name, "8-bit black", width, height, 1);
	roiManager("Select", Axis_roi);roiManager("Fill");run("Select None");
	if (XY){//Repeat for autosome mask. Default to using the next Roi i manager after the one selected, since this is how they're recorded
		Autosome_roi=Axis_roi+1;
		Autosome_output_name="Autosome_"+Axes_output_name;
		//Make mask for axis using the roi
		newImage(Autosome_output_name, "8-bit black", width, height, 1);
		roiManager("Select", Autosome_roi);roiManager("Fill");run("Select None");		
	}}
	

//Skeletonize
//Skeletonize all axes
if (reusing_axis_mask==0){//only skeletonize if this hasn't been done before
	selectWindow(Axes_output_name);run("Select None");run("Duplicate...","title=Axes_Skeleton");run("Skeletonize");
	
	selectWindow("Axes");run("Duplicate...","title=Axes_Copy");run("8-bit");
	run("Concatenate...", "  title=[Axes Skeleton QC] image1=[Axes_Skeleton] image2=Axes_Copy image3=[-- None --]");
	run("Make Composite", "display=Composite");
	
	Stack.setChannel(1);run("Green");
	Stack.setChannel(2);run("Grays");
	setSlice(1);run("Channels Tool...");selectWindow("Channels");
	waitForUser("Clean up skeleton", "Skeleton connections (Green) can be added or removed if necessary.\n \nUse \"Eraser Tool\" in \"Drawing Tools\"and \"Pencil Tool\" set to 1px.\n\nOriginal axis channel (Grey) is a reference and can be turned on/off using the \"Channels Tool\"\nClick \"OK\" when complete.");//users can remove signal with eraser e.g. from SYCP3 aggregates. Switch between original and skeletonized channels for comparison
	close("Channels");
	
	selectWindow("Axes Skeleton QC");run("Duplicate...", "title=Accepted_Skeleton duplicate channels=1");run("Grays");//add roi of skeleton to manager
	run("Create Selection");roiManager("Add");run("Select None");
	ROI_position=roiManager("count")-1;
	roiManager("Select", ROI_position);roiManager("Rename", Axes_output_name+"_Skeleton");
	
	//Measure full skeleton length
	selectWindow("Accepted_Skeleton");
	run("Analyze Skeleton (2D/3D)", "prune=none show");
	close("Results");selectWindow("Branch information");
	skel_len=Table.getColumn("Euclidean distance");
	Array.getStatistics(skel_len, min, max, mean, stdDev);
	number_of_elements=lengthOf(skel_len);
	total_skel_len=mean * number_of_elements;
	close("Branch information");
	
	//Removing XY automatically using autosome mask if requested at start 
	if (XY) {
	selectWindow("Accepted_Skeleton");run("Select None");run("Duplicate...","title=Autosome_Skeleton");
	roiManager("Select", Autosome_roi);run("Clear Outside");run("Select None");
	run("Create Selection");roiManager("Add");run("Select None");//add roi
	ROI_position=roiManager("count")-1;roiManager("Select", ROI_position);roiManager("Rename", "Autosome_"+Axes_output_name+"_Skeleton");
	
	//Measure autosomal skeleton
	selectWindow("Autosome_Skeleton");
	run("Analyze Skeleton (2D/3D)", "prune=none show");
	close("Results");selectWindow("Branch information");
	skel_len=Table.getColumn("Euclidean distance");
	number_of_elements=lengthOf(skel_len);
	autosome_skel_len=mean * number_of_elements;
	close("Branch information");
	}
	
	//Assess skeleton quality
	selectWindow("Accepted_Skeleton");
	Dialog.createNonBlocking("Skeleton Assessment");
	details=newArray("Good for axis measurement", "Poor - exclude");
	Dialog.addRadioButtonGroup("Skeleton Quality", details, 1, 2, "Good for axis measurement");
	Dialog.addString("Additional Comments","");
	Dialog.show;
	
	QualitySelection=Dialog.getRadioButton();
	SkeletonComments=Dialog.getString();
	}

//2c. Segment foci

//When watershedding masks using maxima, the particles without a maxima will be removed. This function adds them back in so that watershedding will only separate foci and increase counts
function add_watershed_foci_to_original_labelmap(Binary_original_foci, Watershed_labelmap_foci, Output_name){
	selectWindow(Watershed_labelmap_foci);run("LabelMap to ROI Manager (2D)");run("Select None");//Add watershed foci to the ROI manager (individual ROIs)
	//Then work out which labels were lost with the watershedding (didn't have a maxima detected)
	selectWindow(Watershed_labelmap_foci);run("Duplicate...", "title=Focus_Mask_New");
	setThreshold(1.0000, 1000000000000000000000000000000.0000);run("Convert to Mask");
	run("Create Selection");roiManager("Add");New_Foci_Mask_Roi=roiManager("count")-1;
	
	selectWindow(Binary_original_foci);run("Select None");run("Duplicate...", "title=Focus_Mask_Original");
	setThreshold(1.0000, 1000000000000000000000000000000.0000);run("Convert to Mask");
	
	selectWindow("Focus_Mask_Original");roiManager("Select", New_Foci_Mask_Roi);
	run("Clear", "slice");run("Select None");

	//Are there any foci missing in the watershed image from the original image? If not, the image will now be blank
	run("Clear Results");run("Set Measurements...", "min redirect=None decimal=3");run("Measure");//If there are no missed foci, no signal will be detected
	if (getResult("Max", 0)>0){
		selectWindow("Focus_Mask_Original");run("Threshold to label map (2D, 3D)", "threshold=1");//If there are missed foci, make these foci a labelmap
		run("LabelMap to ROI Manager (2D)");//add the foci to the ROI manager
		}
		roiManager("Select", New_Foci_Mask_Roi);roiManager("delete");run("Select None");//remove the roi of the total new foci
		
//	Use the foci in the ROI manager to make a labelmap
//	close(Binary_original_foci);close(Watershed_labelmap_foci);//close old labelmap to make way for new one
	run("ROI Manager to LabelMap(2D)");rename(Output_name);
	}


reusing_foci_labelmap=(Foci_channel=="Pre-existing labelmap");
if (reusing_foci_labelmap==0){

	//Loop enables the detection of zero foci
	set_foci_threshold=true;
	while(set_foci_threshold){
		selectWindow("Original");run("Select All");run("Duplicate...", "title=Foci duplicate channels="+Foci_channel);
		selectWindow("Foci");run("Duplicate...","title=Foci_Mask");
		roiManager("Select", 0);run("Clear Outside");run("Select None");//clear signal outside nucleus
		roiManager("Select", 1);run("Add Selection...");run("Select None");//add an outline overlay of the segmented axis mask for reference
		run("Brightness/Contrast...");waitForUser("Adjust Foci Brightness/Contrast", "When ready, click \"OK\". Do not click \"Apply\""); //Need to press OK to continue
		selectWindow("B&C");run("Close");
		//Threshold
		run("Threshold...");call("ij.plugin.frame.ThresholdAdjuster.setMode", "Red");
		waitForUser("Adjust threshold limits, making foci territory red. \nDo not \"Apply\", but click \"OK\" when complete"); //Need to press OK to continue
		
		run("Convert to Mask");
		selectWindow("Threshold");run("Close");
		selectWindow("Foci_Mask");run("Remove Overlay");
		
		//Have any foci been detected? If not this should be counted as 0
		run("Set Measurements...", "min redirect=None decimal=3");
		run("Measure");no_foci=getResult("Max", 0)==0;//is the max value in the mask 0? i.e. have no foci been detected?
		close("Results");
		
		
		if (no_foci){
			
			Dialog.createNonBlocking("No foci detected");
			Dialog.addMessage("No foci have been detected. Is this correct? If not then repeat threshold setting");
			no_foci_options=newArray("OK, continue", "Repeat threshold setting");
			Dialog.addRadioButtonGroup("", no_foci_options, 1, 2, "Repeat threshold setting");
			Dialog.show();
			
			decision=Dialog.getRadioButton();
			
			close("Foci_Mask");
			
			if(decision=="OK, continue"){
				set_foci_threshold=false;
				}
			
		}
		
		if(!no_foci){	
			set_foci_threshold=false;//change this to allow algorithm to progress
			run("Create Selection");roiManager("Add");run("Select None");
			Foci_roi=roiManager("count")-1;
			roiManager("Select", Foci_roi);
			roiManager("Rename", Foci_output_name+"_Total");			
			}}
	
		//If foci are detected then watershed separate
		if(!no_foci){
		selectWindow("Foci");run("Duplicate...","title=Foci_Overlay");run("Yellow");
		run("Colors...", "foreground=white background=black selection=darkgray");
		roiManager("Select", Axis_roi);run("Add Selection...");run("Select None");//add an outline overlay of the segmented axis mask for reference
		run("Colors...", "foreground=white background=black selection=cyan");
		roiManager("Select", Foci_roi);run("Add Selection...");run("Select None");	
		run("Brightness/Contrast...");
		
		///Focus maxima setting for watershedding
		//Place holders to set loop
		prominence=0;
		decision="";
		overlay_with_maxima="overlay_without_maxima";
		//Select level of Gaussian filtering for axis segmentation
		while (decision!="Accept and continue"){
			Dialog.createNonBlocking("Maxima Watershed Setting");
			Dialog.addMessage("Set maxima prominence to mark individual foci at sites where they are merged together.\nLower values result in greater separation.\nRun \"Preview\" first to see if merged foci (cyan) acquire separate points (magenta).\nAdjustment of Brightness & Contrast may be necessary to view original staining pattern. Do not click \"Apply\".");
			Dialog.addString("Prominence", prominence);
			options = newArray("Reset", "Preview", "Accept and continue");
			Dialog.addRadioButtonGroup("", options, 3, 1, "Reset");
			Dialog.show;
			
			prominence=Dialog.getString();
			decision=Dialog.getRadioButton();
			
			if (decision=="Reset"){
				if (overlay_with_maxima!="overlay_without_maxima"){			
					selectWindow("overlay_with_maxima");run("Remove Overlay");run("Select None");
					run("Colors...", "foreground=white background=black selection=darkgray");
					roiManager("Select", Axis_roi);run("Add Selection...");
					run("Colors...", "foreground=white background=black selection=cyan");
					roiManager("Select", Foci_roi);run("Add Selection...");
				}
				else{
					waitForUser("No maxima setting to reset");
					selectWindow("Foci_Overlay");
				}}
			if(decision=="Preview"){
				if(overlay_with_maxima=="overlay_without_maxima"){
					overlay_with_maxima="Preview Maxima Prominence = "+prominence;
					selectWindow("Foci_Overlay");rename(overlay_with_maxima);			
				}
				selectWindow(overlay_with_maxima);
				overlay_with_maxima="Preview Maxima Prominence = "+prominence;//prominence;#reset title to new prominence
				rename(overlay_with_maxima);
				run("Find Maxima...", "prominence="+prominence+" output=[Single Points]");rename("Maxima");
				roiManager("Select", Foci_roi);run("Clear Outside");run("Select None");//clear maxima outside of foci mask		
				run("Colors...", "foreground=white background=black selection=darkgray");
				run("Create Selection");roiManager("Add");Maxima_roi=roiManager("count")-1;//select the maxima points
				
				selectWindow(overlay_with_maxima);run("Remove Overlay");run("Select None");
				run("Colors...", "foreground=white background=black selection=darkgray");
				roiManager("Select", Axis_roi);run("Add Selection...");
				run("Colors...", "foreground=white background=black selection=cyan");
				roiManager("Select", Foci_roi);run("Add Selection...");
				run("Colors...", "foreground=white background=black selection=magenta");
				roiManager("Select", Maxima_roi);run("Add Selection...");run("Select None");	
				close("Maxima");roiManager("Select", Maxima_roi);roiManager("Delete");run("Select None");	
				
				}
				}
		
		run("Colors...", "foreground=white background=black selection=darkgray");
		
		if (overlay_with_maxima!="overlay_without_maxima"){
			selectWindow(overlay_with_maxima);rename("Foci_Overlay");
		}
		
		selectWindow("Foci_Overlay");
		run("Find Maxima...", "prominence="+prominence+" output=[Single Points]");
		roiManager("Select", Foci_roi);run("Clear Outside");run("Select None");rename("Maxima_Binary");//clear maxima outside of foci mask		
		
		
		//Add maxima to rois
		run("Create Selection");roiManager("Add");run("Select None");
		Maxima_roi=roiManager("count")-1;
		roiManager("Select", Maxima_roi);
		roiManager("Rename", Foci_output_name+"_Maxima")
		
		setBatchMode("hide");//enter batch mode
		
		//Save the ROIs and close. A clear ROI manager is needed to combined the 
		roiManager("Save", dir2 + "/RoiSet.zip");roiManager("reset");//save ROIs and close manager
		
		//Convert to watershed labelmap
		run("Threshold to label map (2D, 3D)", "threshold=1");rename("Foci_Maxima");//Convert to labelmap so maxima can act as seed points for watershedding
		close("Maxima_Binary");
		run("Watershed with seed points (2D, 3D)", "image_to_segment=Foci_Mask image_with_seed_points=Foci_Maxima use_threshold threshold=1");
		rename("Watershed_labelmap_foci");
		//Now combined the watershed foci with any foci lost from the original mask (if a focus lacks a maxima point, they will be lost. This Process add them back)
		add_watershed_foci_to_original_labelmap("Foci_Mask", "Watershed_labelmap_foci", Foci_output_name);
		}}

	else {//if reusing existing roi for foci then define whih one. 
		waitForUser("Select Focal Roi in Manager");
		//		setBatchMode("hide");//enter batch mode
		Foci_output_name=Roi.getName;
		Foci_roi=roiManager("index");run("Select None");
		//Make mask for foci using the roi
		newImage("Foci_Mask", "8-bit black", width, height, 1);
		roiManager("Select", Foci_roi);roiManager("Fill");run("Select None");
		newImage("Foci_Maxima", "8-bit black", width, height, 1);
		Foci_Maxima_roi=Foci_roi+1;
		roiManager("Select", Foci_Maxima_roi);roiManager("Fill");run("Select None");
		roiManager("Save", dir2 + "/RoiSet.zip");roiManager("reset");//save ROIs and close manager
		run("Threshold to label map (2D, 3D)", "threshold=1");rename("Foci_Maxima_Labelmap");//Convert to labelmap so maxima can act as seed points for watershedding
		run("Watershed with seed points (2D, 3D)", "image_to_segment=Foci_Mask image_with_seed_points=Foci_Maxima_Labelmap use_threshold threshold=1");
		rename("Watershed_labelmap_foci");
		add_watershed_foci_to_original_labelmap("Foci_Mask", "Watershed_labelmap_foci", Foci_output_name);
		no_foci=false;
		}


///////////////////Batch-mode entered as no more manual intervention is required

//3. Calculate Measurements
//Features to measure
img_name=newArray();
stage=newArray();
label=newArray();
area=newArray();
circularity=newArray();
centroid_x=newArray();
centroid_y=newArray();
distance_axes_mask=newArray();
distance_autosome_mask=newArray();

//If there are foci detected
if(!no_foci){
	selectWindow(Foci_output_name);
	run("Analyze Regions", "area circularity centroid");//Generate area, circularity and centroid measurements
	Table.rename(Foci_output_name+"-Morphometry", "Results");
	//Populate label, area, and centroid arrays
	for (i = 0; i < nResults(); i++) {
	    img_name=Array.concat(img_name,title);
	    stage=Array.concat(stage,Manual_Stage);
	    label = Array.concat(label,i+1);
		area = Array.concat(area,getResult("Area",i));
		circularity=Array.concat(circularity,getResult("Circularity",i));
		centroid_x = Array.concat(centroid_x,getResult("Centroid.X", i));
		centroid_y = Array.concat(centroid_y,getResult("Centroid.Y", i));
		}
	close("Results");
	
	//Foci distance to axis mask
	selectWindow(Axes_output_name);run("Select None");run("Duplicate...","title=Axes_Distance");run("Invert");run("Distance Map");//make distance map of Axes_Mask to measure focus distances
	run("Intensity Measurements 2D/3D", "input=[Axes_Distance] labels="+Foci_output_name+" min");//minimum distance from focus to axis 
	Table.rename("Axes_Distance-intensity-measurements", "Results");
	
	for (i = 0; i < nResults(); i++) {
		distance_axes_mask = Array.concat(distance_axes_mask,getResult("Min",i));
	}
	close("Results");
	}

//Repeat for autosome mask
if(!no_foci){
	if (XY) {
		selectWindow(Autosome_output_name);run("Select None");run("Duplicate...","title=Autosome_Axes_Distance");run("Invert");run("Distance Map");//make distance map of Autosome_Axes_Mask to measure focus distances
		run("Intensity Measurements 2D/3D", "input=[Autosome_Axes_Distance] labels="+Foci_output_name+" min");//minimum distance from focus to autosome axis 
		Table.rename("Autosome_Axes_Distance-intensity-measurements", "Results");
		
		for (i = 0; i < nResults(); i++) {
			distance_autosome_mask = Array.concat(distance_autosome_mask,getResult("Min",i));
		}
		close("Results");
		}}

//Compile data into table
//If Image_Summary.csv has been made before then add to it
if (Image_not_previously_scored){
	if(!no_foci){
		if (XY) {
			Table.showArrays("Focus_measurements", img_name, stage, label, area, circularity, centroid_x, centroid_y, distance_axes_mask, distance_autosome_mask)}
		else {
			Table.showArrays("Focus_measurements", img_name, stage, label, area, circularity, centroid_x, centroid_y, distance_axes_mask)}
		
		//Summarise focus counts
		selectWindow("Focus_measurements");
		distance_axes_mask=Table.getColumn("distance_axes_mask");
		if (XY){
			distance_autosome_mask=Table.getColumn("distance_autosome_mask");
			}
		
		//Count the number of foci on axes
		foci_on_axes=0;
		if (XY){
			foci_on_autosomes=0;
		}
		foci_off_axes=0;
		selectWindow("Focus_measurements");
		for (i = 0; i < distance_axes_mask.length; i++) {
			focus_to_axes=distance_axes_mask[i];
		    if (focus_to_axes==0){
		    	foci_on_axes+=1;
		    	}
		    	else{
		    		foci_off_axes+=1;
		    	}
		   	if(XY){
		   	focus_to_autosomes=distance_autosome_mask[i];
		    if (focus_to_autosomes==0){
		    	foci_on_autosomes+=1;
		    	}}}}
    	
	//Compile a table of basic information 
	setResult("Stage",0,Manual_Stage);
	setResult("Image_Comments",0,Comments);
	setResult(Axes_output_name+"_Sigma_for_Gaussian_Filter", 0, sigma);					
	setResult(Axes_output_name+"_Skeleton_Quality", 0, QualitySelection);					
	setResult(Axes_output_name+"_Skeleton_Comments",0,SkeletonComments);
	setResult(Axes_output_name+"_Euclidean_Total_Skeleton_Length",0,total_skel_len);
	if(!no_foci){
		setResult("Maxima_Prominence_"+Foci_output_name,0,prominence);
		setResult(Foci_output_name+"_Foci_on_"+Axes_output_name,0,foci_on_axes);
		setResult(Foci_output_name+"_Foci_off_"+Axes_output_name,0,foci_off_axes);
	}
	else{
		setResult(Foci_output_name+"_Foci_on_"+Axes_output_name,0,0);
		setResult(Foci_output_name+"_Foci_off_"+Axes_output_name,0,0);
		}
		
	if (XY)	{
		setResult("XY Analysis",0,"Yes");
		setResult("Autosomal_"+Axes_output_name+"_Euclidean_Skeleton_Length",0,autosome_skel_len);
		if(!no_foci){
			setResult(Foci_output_name+"_Foci_on_autosomal_"+Axes_output_name,0,foci_on_autosomes);
		}
		else{
			setResult(Foci_output_name+"_Foci_on_autosomal_"+Axes_output_name,0,0);
		}
	}
	else {
		setResult("XY Analysis",0,"No");
	}
	updateResults();
}

	//if adding new data to existing Image_Summary.csv file
	else{
		if(!no_foci){	
			if (XY) {
				Table.showArrays("Focus_measurements", img_name, stage, label, area, circularity, centroid_x, centroid_y, distance_axes_mask, distance_autosome_mask)}
			else {
				Table.showArrays("Focus_measurements", img_name, stage, label, area, circularity, centroid_x, centroid_y, distance_axes_mask)}
		
			//Summarise focus counts
			selectWindow("Focus_measurements");
			distance_axes_mask=Table.getColumn("distance_axes_mask");
			if (XY){
				distance_autosome_mask=Table.getColumn("distance_autosome_mask");
				}
			
			//Count the number of foci on axes
			foci_on_axes=0;
			foci_on_autosomes=0;
			foci_off_axes=0;
			selectWindow("Focus_measurements");
			for (i = 0; i < distance_axes_mask.length; i++) {
				focus_to_axes=distance_axes_mask[i];
			    if (focus_to_axes==0){
			    	foci_on_axes+=1;
			    	}
			    	else{
			    		foci_off_axes+=1;
			    	}
				if (XY){		    	
				   	focus_to_autosomes=distance_autosome_mask[i];
				    if (focus_to_autosomes==0){
				    	foci_on_autosomes+=1;
				    	}}}}
		    	else{
		    		foci_on_axes=0;
					foci_on_autosomes=0;
					foci_off_axes=0;
		    	}
				    	
	    IJ.renameResults("Image_Summary.csv","Results"); 	
    	if (reusing_axis_mask==0){
    		setResult(Axes_output_name+"_Sigma_for_Gaussian_Filter", 0, sigma);					
	    	setResult(Axes_output_name+"_Skeleton_Quality", 0, QualitySelection);					
			setResult(Axes_output_name+"_Skeleton_Comments",0,SkeletonComments);
			setResult(Axes_output_name+"_Euclidean_Total_Skeleton_Length",0,total_skel_len);
			if (XY)	{
				setResult("Autosomal_"+Axes_output_name+"_Euclidean_Skeleton_Length",0,autosome_skel_len);
				}}
		

		setResult(Foci_output_name+"_Foci_on_"+Axes_output_name,0,foci_on_axes);
		setResult(Foci_output_name+"_Foci_off_"+Axes_output_name,0,foci_off_axes);
		if (XY)	{
			setResult(Foci_output_name+"_Foci_on_autosomal_"+Axes_output_name,0,foci_on_autosomes);
					}
		if(reusing_foci_labelmap==0){
			if(!no_foci){
			setResult("Maxima_Prominence_"+Foci_output_name,0,prominence);
			}}
	updateResults();
		}



//4. Save data 

//save data and masks/labelmaps and rois
selectWindow("Results");saveAs("Results",dir2 + "/Image_Summary.csv");close("Image_Summary.csv");close("Results");

if(!no_foci){
	selectWindow("Focus_measurements");saveAs("Results",dir2 + "/"+Foci_output_name+"_measurements_to_"+Axes_output_name+".csv");close(Foci_output_name+"_measurements_to_"+Axes_output_name+".csv");
}
//save images, labelmaps and ROIs
if (reusing_axis_mask==0){
	selectWindow(Axes_output_name);saveAs("tiff", dir2 + "/"+Axes_output_name);
	if(XY){
	selectWindow(Autosome_output_name);saveAs("tiff", dir2 + "/"+Autosome_output_name);
	}}
if(reusing_foci_labelmap==0){
	if(!no_foci){
	selectWindow(Foci_output_name);saveAs("tiff", dir2 + "/"+Foci_output_name);
	}}


roiManager("Save", dir2 + "/"+Foci_output_name+"_ROIs.zip");//save ROIs

close("*");close("B&C");close("Channels");//close all images and tools
close("ROI Manager");
//print a summary of focs measurements to log
setBatchMode("show");
print("Image analysis complete");
print("Foci on "+Axes_output_name+":"+foci_on_axes);
if(XY){
print("Foci on "+Axes_output_name+" autosomes:"+foci_on_autosomes);
}
print("Foci off "+Axes_output_name+":"+foci_off_axes);
