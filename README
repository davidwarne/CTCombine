CTCombine Version 0.195

# Summary
#---------

  Rotates Dicom data as indicated a desired angle about a desired isocentre. 
  Voxelizes data and appends EPID data to the file.


# Authors
#---------

    Mark Dwyer    ( m2.dwyer@qut.edu.au )
    David Warne   ( david.warne@qut.edu.au )

# Institution
#-------------

    High Preformance Computing and Research Support
    Queensland University of Technology

# Licenses
#----------

 Feel free to copy/modify all or any portion of this as so desired.

 No guarantee is made as to the accuracy of this code.  Probably best
 not to use in a nuclear reactor.

# Package structure :
#--------------------

  The directory Version_0.195/ is organized as follows :

    - EGSPhant.h                 : The single (header) file of the file format storage.
    - DicomReader.h              : Code that reads Dicom data from .dcm files (using DICOMParser), also contains rotation/translation function
    - Main.cpp		         : My dodgy kick off file containing the main() 
    - Makefile                   : How I like to compile things
    - README.txt                 : You be readin' it
    - DICOMParser/               : contains code that deals with extracting data from dicom files
    - EPIDexample.txt            : example of the EPID input data format
    - cthead.inp                 : example of .inp file for voxel conversion

    - documentation/ : Ha Ha Ha haaa (wipes tear from eye) Good one!

# Getting Started
#-----------------

    At the command line, just type "make".  This will work on *nix systems.

    To Merge to .egsphant files:
       
         CTCombine -c file1.egsphant file2.egsphant outputfile.egsphant
    example:
          CTCombine -c ./cthead.egsphant ./EPID.egsphant ./cthead_and_EPID.egsphant
    
   
To Rotate about a the isocentre to align with a given gantry angle:

        CTCombine [-g zeroGantry] [-r desiredGantry] -E EPIDdatafile -inp .inpInputFile -o outputfile

   example:
        CTCombine -g 270 90 -r 300 90 -E ./EPIDexample.txt -inp ./cthead.inp -o ./cthead.egsphant

NOTE: the zero gantry will be theta=270 phi=90 by default. The following example will simply convert from 
DICOM to .egsphant and include the EPID: 
   example:
        CTCombine -E ./EPIDexample.txt -inp ./cthead.inp -o ./cthead.egsphant

Other options (these are all optional):

-d dir            Overrides directory stated in .inp file for DICOM data reading
example:          CTCombine -d ./Data/cthead/ -E ./EPIDexample.txt -inp ./cthead.inp -o ./cthead.egsphant

-i x y z          translates volume to this point to create a new isocentre, any rotation is done about this point
example:          CTCombine -i 0.053 -22.2 0 -g 270 90 -r 300 90 -E ./EPIDexample.txt -inp ./cthead.inp -o ./cthead.egsphant      
    
-e dist           specifies the distance (in cm) that the surface of the EPID is to be from the isocentre
example:          CTCombine -e 59.4 -E ./EPIDexample.txt -inp ./cthead.inp -o ./cthead.egsphant
     or           CTCombine -e 59.4 -c ./cthead.egsphant ./EPID.egsphant ./cthead_and_EPID.egsphant

# File Format Description
#-------------------------

EPID Specification file:

xdimesion,zdimension //layer surface area
numberOfMaterials    // count of materials in both CT and EPID
material1name        //List of material names
material2name
material3name
.
.
.
materialNname       //end of list
numberOfLayers      //count of EPID layers 
materialNumber1 materialdensity1 thickness1 //List of layer information
materialNumber2 materialdensity2 thickness2
.
.
.
materialNumberN materialdensityN thicknessN

- for an example see the EPIDexample.txt file

INP file (you already know this...): 

DataType                                                                  // only DICOM is handled
directory                                                                 // the folder containing the data files
xMin, xMax, yMin, yMax, zMin, zMax                                        // boundaries of selected data region
xVoxeldimension, yVoxeldimension,zVoxeldimension                          //size of Voxel Dimensions in cm
numberOfMaterials                                                         // Count of Materials in CT
numberOfControlPoints1 numberOfControlPoints2 ... numberOfControlPointsN  // control point numbers
materialName1                                                             // material and conversion data
CT11,Density11,CT12,Density12,...,CT1N,Density1N, estepe1                 //list of CT,Density pairs forming the transfer function followed by estepe
materialName2                                                
CT21,Density21,CT22,Density22,...,CT2N,Density2N, estepe2
.
.
.
materialNameN                                                
CTN1,DensityN1,CTN2,DensityN2,...,CTNN,DensityNN, estepeN

# History
#--------

  Features:
  	
  	Version 0.19 - 
  		- Output reporting updated to be cleaner and more readable.
  		- Some issues with remapping mesh from selection region are resolved.
  	Version 0.18 - 
  		- The Volume clipping problem when the volume is rotated or translated to the extremes is now fixed
    Version 0.17 - 
		- Dicom reader now takes into account the directional cosines of the patients orientation
		- Dicom Data Types of float, unsigned char, short, and unsigned short are now handled (though only short has been tested)
    
    Version 0.16 - 
		- density data conversion is now more powerful, the user may define there own transfer function for 
		CT to density (NOTE: as a result of this modification the format of the .inp file has changed, see below)
		- x and z cooridnates are set in increasing order for use with Dosxyz
		- Bug with EPID distance is now fixed
		- problem with converting CT values to densities is now removed by implementing the use of the 
		Dicom RescaleOffset value (useful, the poor phantom now has a brain)
		- the selected isocentre is now always set in the centre of the zx plain (0,100,0) for Dosxyz
    
    Version 0.15 -
		- Y axis now shifted such that y<100 is above the isocentre and y>100 is below

    Version 0.14 -
		- volume extended laterally to fully include the EPID

    Version 0.13 -
		- Input validation has been implemented
		- some minor bug fixes 
     
    From Version 0.12 -
		- Reads DICOM data, and converts to .egsphant using voxel sizes and coordinate limts specified in .inp file
		- Rotates the data to a geometrically equivalent orientaion for a new gantry angle (as an input)
		- Translates the isocentre to 0,0,0 and inverts the y axis
		- Adds EPID data to .egsphant as specified from an EPID spec file (see below)
		- Pads the volume with air to fill in space between EPID and phantom, and along the edges where the EPID extends

    From Version 0.11 -
		- Will read two (2) egsphant files, combine them and write out.  No smarts implemented.

    From Version 0.1 -
		- fixed integer writing to a specified number of characters
		
# Known Issues
--------------
Version 0.195 - 
  		- error revealed with region selection for datasets using directional cosines
  		  other than (1,0,0), (0,1,0)... dodgy hack for signed short data sets
  		  has been implemented, but needs to be fixed properly. 
