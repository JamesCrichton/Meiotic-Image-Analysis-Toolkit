# Meiosis_Image_Analysis_Toolkit
A collection of scripts for the analysis of meiotic prophase I fluorescence images. 

Download the contents of this repository by clicking Code > Download ZIP and extracting, or by cloning using git. 


## Running ImageJ macro scripts
ImageJ script guides users through the segmentation of axial and focal staining patterns, measuring the frequecy of foci, length of axes and generating data which can be further analysed using subsequent python scripts. 

- Running this macro first requires the installation of [Fiji](https://fiji.sc/). 
- Open Fiji and drag and drop the downloaded macro file into the Fiji toolbar to open. The code will appear in the script editor.
- Open your image of interest (this should be a stack image). 
- Click "Run" to start the macro.


## Running Python scripts

Python functions created are packaged into a single file "meiosis_toolkit.py". 

Several scripts also require [TracePy](https://github.com/huangziwei/TracePy/tree/master), which should also be downloaded. 

Example scripts employing the these to analyse batches of images are in the "Example_Python_Scripts" folder. These are setup to use the example image data included in this repository. These scripts reqire that the location of the downloaded repositories ("TracePy" and "Meiosis_Image_Analysis_Toolkit") are specified, so they can be installed and sample images accessed. Code is written with the **Desktop** being the default location to which these repositories are saved. 

To run on your own image data or store the repositories elsewhere, paths to image directories must be changed appropriately. 

Example scripts enable users to:
- Calculate segmented focus position along traced axes ("Focus_positioning_on_traces.py")
- Shuffle foci randomly in nuclear space and calculate distances to another marker e.g. axes ("Nuclear_focus_shuffling.py")
- Shuffle foci randomly in axial space and calculate distances to another merger e.g. another axial focus ("Axis_focus_shuffling.py")
- Calculate the proximity of pairs of axis traces ("Paired_axis_proximity.py")
- Calculate the distance of segmented foci from axes and orientation between/outside of pairs of axes ("Focus_position_relative_to_individually_traced_axis_pairs.py")

