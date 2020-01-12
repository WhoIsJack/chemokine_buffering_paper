
// Created on Thu July 16 15:34:56 2015
//
// @author:    Jonas Hartmann @ Gilmour group @ EMBL Heidelberg
// 
// @descript:  Perform conversion of a microscopy stack or time course to an
//             8bit tif based either on a given range ('camera') or on the min
//             and max value in the stack itself.


// USER INPUT: Loading & saving
input_dir = "E:\\path\\to\\input_directory\\";
output_dir = input_dir;
input_ending = ").czi";
stack_order = "XYCZT";

// USER INPUT: Conversion settings
use_minmax = 0;  // If 1, use the stacks minimum and maximum values, else use the 'camera' values
camera_min = newArray(10000, 10000);
camera_max = newArray(35000, 35000);
reorder = "21";
make_gray = 1;


// Function to convert range(camera_min,camera_max) to 8bit
function make_8bit_camera(input_dir, output_dir, filename, stack_order, camera_min, camera_max, reorder) {
	
	// Open stack using Bio-Formats
	inpath = input_dir + filename;
	run("Bio-Formats Importer", "open='" + inpath + "'autoscale color_mode=Default view=Hyperstack stack_order="+stack_order);
	
	// Get stack dimensions
	getDimensions(w, h, channels, slices, frames);

	// Reset histogram for each channel
	for (ch=1; ch<=channels; ch++) {
		Stack.setChannel(ch);
		setMinAndMax(camera_min[ch-1], camera_max[ch-1]);
	}
		
	// Convert to 8-bit
	run("8-bit");

	// Reorder channels as needed
	run("Arrange Channels...", "new="+reorder);

	// Make gray if needed
	if (make_gray==1){
		if (channels>1){
			Stack.setDisplayMode("grayscale")
		}else{
			run("Grays");
		}
	}
	
	// Save output
	outpath = output_dir + substring(filename, 0, indexOf(filename, ".")) + '_8bit.tif';
	saveAs(outpath);
	close();
}


// Function to convert range(min,max) to 8bit
// Note: min and max are taken from z-projections of the first frame!
function make_8bit_minmax(input_dir, output_dir, filename, stack_order, reorder) {
	
	// Open stack using Bio-Formats
	inpath = input_dir + filename;
	run("Bio-Formats Importer", "open='" + inpath + "'autoscale color_mode=Default view=Hyperstack stack_order="+stack_order);
	
	// Get stack dimensions
	getDimensions(w, h, channels, slices, frames);

	// For each channel
	// Note: This needs more z projects than really necessary but is easier to read...
	for (ch=1; ch<=channels; ch++) {
		
		// Maximum project to get overall maximum value
		run("Z Project...", "projection=[Max Intensity]");
		Stack.setChannel(ch);
		getRawStatistics(area, mean, minWrong, max);
		close();
	
		// Minimum project to get overall minimum value
		run("Z Project...", "projection=[Min Intensity]");
		Stack.setChannel(ch);
		getRawStatistics(area, mean, min, maxWrong);
		close();

		// Set channel
		Stack.setChannel(ch);
	
		// Set min and max to the range
		setMinAndMax(min,max);
	}
		
	// Convert to 8-bit
	run("8-bit");

	// Reorder channels as needed
	run("Arrange Channels...", "new="+reorder);
	
	// Save output
	outpath = output_dir + substring(filename, 0, indexOf(filename, ".")) + '_8bit.tif';
	saveAs(outpath);
	close();
}


// Iterate over files
filelist = getFileList(input_dir);
for (i=0; i<filelist.length; i++) {

	// Check if file has correct ending
	if (endsWith(filelist[i],input_ending)) {
		
		// Make 8bit with camera_max
		if (use_minmax==0) {
        	make_8bit_camera(input_dir, output_dir, filelist[i], stack_order, camera_min, camera_max, reorder);
		}
		
        // Make 8bit with minmax
		else {
			make_8bit_minmax(input_dir, output_dir, filelist[i], stack_order, reorder);
		}
	}
}



