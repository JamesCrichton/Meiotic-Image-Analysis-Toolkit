//Macro allows user to select a source directory containing files of interest and a destination directory for renamed files. 
//Files are given a new numerical string as a name, and saved in the destination folder as tiffs.
//The original names and new codenames are saved in an .csv spreadsheet

//Clear existing results
run("Clear Results");     

//Set batch mode, define source and destination directories, get a list of files in the source directory
setBatchMode(true);
dir1 = getDirectory("Source directory");
dir2 = getDirectory("Destination directory");
list = getFileList(dir1);

//Initialise an array to store used codenames
codenames_assigned=newArray();

//Loop thorugh the file list from the source folder
for (i=0; i<list.length; i++){
	//Get the path for the first file, if it is not a directory then process it
	path = dir1+list[i];
	if (! File.isDirectory(path)){
		open(path);
		title = File.nameWithoutExtension;
		//Set the codename as the title to start with
		codename = title;
		//Make a random number as string and use a while loop controlled by the name_OK variable to make sure number isn't assigned more than once
		name_OK=0; 
		while (name_OK == 0){
			//Creates random number between 0-1, make it an integer, then round to 0 decimal places and ovewrite the codename
			a=random; 
			b=(a*10000000);
			codename=d2s(b, 0);
			//Codename has been set so OK to exit while loop, reset name_OK to 0 if codename has already been assigned
			name_OK=1; 
			for (j = 0; j < lengthOf(codenames_assigned); j++) {		
	    		if (codenames_assigned[j]==codename){
	    			name_OK=0;
	    		}
	    	}
		}
		
					
		//Add original name and codename to Results table
		row=nResults;
		setResult("Original_Name", row, title);
		setResult("Codename", row, codename);
		updateResults();
		//Store codename in array of assigned codenames
		codenames_assigned=Array.concat(codenames_assigned,codename);
		//Save image under codename then close image windows without saving
		saveAs("Tiff", dir2 + codename);
		close("*");
	}
}
//Save original names and codenames as Encoding_Record.csv then output number of files renames message to user
saveAs("Results", dir2 + "Encoding_Record.csv"); 
print(nResults+" files renamed");
