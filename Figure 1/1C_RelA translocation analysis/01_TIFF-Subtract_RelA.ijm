//Converts nd2 files to 16-bit TIFF. Also opens up as a dialog.

// --- INPUT ---

//Folder containing images
#@ File (label = "Input/Output Directory", style = "directory", value = "~/ImageAnalysisWorkflow/01_TIFF-Subtract") dir

//Folder containing dark frame images
#@ File (label = "Dark-Frame Directory", style = "directory", value = "~/ImageAnalysisWorkflow/00_Setup/Darkfields") dir_df

//Parameters
#@ Boolean (label = "Run C=0", value = false) runC0
#@ String (label = "LUT C=0", value = "Magenta") LUT_0
#@ String (label = "Dark-Field Image Name C=0", value = "AVG_20230922_Darkfield_SD-W1_Cy5_300ms.tif") df_0

#@ Boolean (label = "Run C=1", value = false) runC1
#@ String (label = "LUT C=1", value = "Green") LUT_1
#@ String (label = "Dark-Field Image Name C=1", value = "AVG_20230922_Darkfield_SD-W1_GFP_200ms.tif") df_1

#@ Boolean (label = "Run C=2", value = false) runC2
#@ String (label = "LUT C=2", value = "Cyan") LUT_2
#@ String (label = "Dark-Field Image Name C=2", value = "AVG_20230922_Darkfield_SD-W1_BFP_200ms.tif") df_2

// --- ND2 TO TIFF CONVERSION ---
dir = dir + "/";
dir_df = dir_df + "/";

ext ="nd2";
inList = getFileList(dir);
list = getFromFileList(ext, inList);

// Get file list of .nd2 files
print("Below is a list of files to be converted:");
printArray(list); // Implemented below

setBatchMode(true);

//Image conversion loop
for (i=0; i<list.length; i++) 
{
  inFullname = dir + list[i];
  print("Converting", i+1, "of", list.length, list[i]);
  
  //Apply function to convert file type
  splitBioFormatToTif(inFullname, dir); //Block if using TIFF using '//' at the beggining

  print("...done.");
}

print("::: ND2 to TIFF Conversion Finished :::");

//Function to convert file type
function splitBioFormatToTif(inFullname, dir)
{
    run("Bio-Formats Macro Extensions");
    Ext.setId(inFullname);
    Ext.getSizeT(numT);
    Ext.getSizeC(numC);
    Ext.getSizeZ(numZ);

    for (c = 1; c <= numC; c++)
    {
        coptions = newArray(c, c, 1);
        toptions = newArray(1, numT, 1);
        zoptions = newArray(1, numZ, 1);

	    name = File.getName(inFullname); 
        lenIndex = lengthOf(name) - 4; 
        title = substring(name, 0, lenIndex); 
		
        outName = title + " C=" + toString(c-1);
        
        if (numC <= 1)
        {            
            outName = title;
        }
   
  		newdir = dir + title + '/';
  		File.makeDirectory(newdir);
  		
        outFullname = newdir + outName + ".tif";
        id = bfImport(inFullname, coptions, zoptions, toptions);
        selectImage(id);
        prop = getMetadata("Info");
        File.saveString(prop, newdir + "Metadata.txt")
        saveAs("Tiff", outFullname);
        close();
    }
}

function bfImport(path,channels,zs,times)
{

  run("Bio-Formats Macro Extensions");
  Ext.setId(path);
  dimOrder = "";
  Ext.getDimensionOrder(dimOrder);
  Ext.getSizeT(numT);
  Ext.getSizeZ(numZ);
  Ext.getSizeC(numC);
  
  options = "open=[" + path + "] view=[Standard ImageJ] stack_order=" + dimOrder + " virtual specify_range ";
  cOpts = "c_begin=" + channels[0] + " c_end=" + channels[1] + " c_step=" + channels[2];
  zOpts = "z_begin=" + zs[0]       + " z_end=" + zs[1]       + " z_step=" + zs[2];
  tOpts = "t_begin=" + times[0]    + " t_end=" + times[1]    + " t_step=" + times[2];
  options = options + cOpts + " " + zOpts + " " + tOpts;
  
  run("Bio-Formats Importer", options);
  id = getImageID();
  return id;
}

function getFromFileList(ext, fileList)
{
  selectedFileList = newArray(fileList.length);
  ext = toLowerCase(ext);
  j = 0;
  for (i=0; i<fileList.length; i++)
    {
      extHere = toLowerCase(getExtension(fileList[i]));
      if (extHere == ext)
        {
          selectedFileList[j] = fileList[i];
          j++;
        }
    }
  selectedFileList = Array.trim(selectedFileList, j);
  return selectedFileList;
}

// Print array items
function printArray(array)
{ 
  for (i=0; i<array.length; i++)
    print(array[i]);
}

function getExtension(filename)
{
  ext = substring( filename, lastIndexOf(filename, ".") + 1 );
  return ext;
}

//Trim tailing backslashes
function trimDirTail(dir)
{
    return replace(dir, "////+$", "");	\\Windows:: return replace(dir, "\\\\+$", "");
}


// --- REMOVE BACKGROUND ---
   count = 0;
   countFiles(dir);	
   n = 0;
   processFiles(dir);
   print(count+" files processed");
   
   function countFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              countFiles(""+dir+list[i]);
          else
              count++;
      }
  }

   function processFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              processFiles(""+dir+list[i]);
          else {
             showProgress(n++, count);
             path = dir + list[i];
            if(runC0 == true)
             	processFile0(path); //For C=0
             if(runC1 == true)
             	processFile1(path); //For C=1
             if(runC2 == true)
             	processFile2(path); //For C=2
          }
      }
  }
  
  //C=0 function, this channel is RelA staining channel
  function processFile0(path) {
  	parent = File.getParent(path) + "/";
	img = substring(path, lengthOf(parent), lengthOf(path)-4);
       if (endsWith(path, " C=0.tif")) {
		//Open image
		open(path);
		//Extract file name
		img = File.nameWithoutExtension();
		open(dir_df + df_0);
		//REMOVE DARK-FIELD
		print("Removing Dark-Field (C=0) from " + img);
		imageCalculator("Subtract create stack", img + ".tif", df_0);
		//Change LUT
		run(LUT_0);
		//Reset brightness and contrast
		resetMinAndMax();
		//Save Dark-Field Removed (DFRm)
		saveAs ("Tiff", dir + img + "-DFRm");
		//Apply median filter with 151-pixel radius
		print("Calculating Median");
		run("Median...", "radius=151 stack");
		//Save median image
		saveAs ("Tiff", dir + img + "-Median");
		//REMOVE MEDIAN
		//Open dark-field
		open(dir + img + "-DFRm.tif");
		//Open Median
		open(dir + img + "-Median.tif");
		//Subtract median from dark-field
		print("Removing Median (C=0) from " + list[i]);
		imageCalculator("Subtract create stack", img + "-DFRm.tif", img + "-Median.tif");
		//Save median removed
		//saveAs ("Tiff", dir + img + "-MedRm");
		saveAs ("Tiff", dir + img + "-zFinal");
		//Close open images
		run("Close All"); 

      }
  }

  //C=1 function, this channel is GFP channel for cell body segmentation, therefore use DFRm
  function processFile1(path) {
  	parent = File.getParent(path) + "/";
	img = substring(path, lengthOf(parent), lengthOf(path)-4);
       if (endsWith(path, " C=1.tif")) {
		//Open image
		open(path);
		//Extract file name
		img = File.nameWithoutExtension();
		open(dir_df + df_1);
		//REMOVE DARK-FIELD
		print("Removing Dark-Field (C=1) from " + img);
		imageCalculator("Subtract create stack", img + ".tif", df_1);
		//Change LUT
		run(LUT_1);
		//Reset brightness and contrast
		resetMinAndMax();
		//Save Dark-Field Removed (DFRm)
		//saveAs ("Tiff", dir + img + "-DFRm");
		saveAs ("Tiff", dir + img + "-zFinal");
		//Close open images
		run("Close All"); 

      }
  }

  //C=2 function, this channel is nucleus channel
  function processFile2(path) {
  	parent = File.getParent(path) + "/";
	img = substring(path, lengthOf(parent), lengthOf(path)-4);
       if (endsWith(path, " C=2.tif")) {
		//Open image
		open(path);
		//Extract file name
		img = File.nameWithoutExtension();
		open(dir_df + df_2);
		//REMOVE DARK-FIELD
		print("Removing Dark-Field (C=2) from " + img);
		imageCalculator("Subtract create stack", img + ".tif", df_2);
		//Change LUT
		run(LUT_2);
		//Reset brightness and contrast
		resetMinAndMax();
		//Save Dark-Field Removed (DFRm)
		saveAs ("Tiff", dir + img + "-DFRm");
		//Apply median filter with 99-pixel radius
		print("Calculating Median");
		run("Median...", "radius=99 stack");
		//Save median image
		saveAs ("Tiff", dir + img + "-Median");
		//REMOVE MEDIAN
		//Open dark-field
		open(dir + img + "-DFRm.tif");
		//Open Median
		open(dir + img + "-Median.tif");
		//Subtract median from dark-field
		print("Removing Median (C=2) from " + list[i]);
		imageCalculator("Subtract create stack", img + "-DFRm.tif", img + "-Median.tif");
		//Save median removed
		//saveAs ("Tiff", dir + img + "-MedRm");
		saveAs ("Tiff", dir + img + "-zFinal");
		//Close open images
		run("Close All"); 

      }
  }
  
print("Done");
