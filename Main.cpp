/**
 * Main.cpp (how descriptive)
 *
 * @Authors: Mark Dwyer, David Warne
 * @Contact: m2.dwyer@qut.edu.au, david.warne@qut.edu.au
 * @Last Modified: 16/02/2009
 *
 * Summary:
 *    Makes use of DicomReader and EGSPhant objects to re-orientate and convert data from DICOM to EGSPhant format.
 *    also attaches EPID Data
 */



#include "./EGSPhant.h"
#include "./DicomReader.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>

using namespace EGS;
using namespace DICOM;

#define AIRDENSE 0.001
#define VERSION 0.13

template<class T>
void byteswap(T* t)
{
	std::reverse(reinterpret_cast<char*>(t), reinterpret_cast<char*>(t + 1));
}



void Get_EGSPHANT_Stats(const char *filename)
{
	
	EGSPhant data((char *)filename);
	int result = data.Read();
	result = data.PrintStats();		
}


/* Globals to storge input argument and file read data
 */
char dir[255];                          // DICOM data directory
char INPfile[255];
float xmin,xmax,ymin,ymax,zmin,zmax;    // require region boundaries of CT data
float vx,vy,vz=0.5;                     // voxel size 
float ix,iy,iz = 0;                     // isocentre
bool flag = 0;
float Z_gantry[2]= {270,90};            // direction of beam at zero gantry angle
float theta = 270;                      // rotation angles
float phi = 90;
float EPID_dist=0;                      // distance of isocentre to EPID data

char EPIDfileName[255];                 //filename of EPID data
char outputfileName[255];               //EGSPHANT output filename

float EPIDx,EPIDz;                      //EPID surface dimension
int totalEPIDLayers;                    // number of EPID dimensions
int EPIDnumMedia; 
float** EPIDLayerInfo;                  // EPIDLayerInfo[0] = material number
                                        // EPIDLayerInfo[1] = material density
                                        // EPIDLayerInfo[2] = Layer Thickness
std::string* EPIDMaterials;             // EPID material names

//similar info for the CT scan
int numMedia;
std::string* MediaNames;
char* mediaCats;
float** conversionData;                 //stores the conversion data



/*
 * Combine the CT with the EPID data, is probably redundant now we create EPID data on conversion From Dicom
 *
 * Params: ct_filename - the phantom dataset
 * 	      epid_filename - the epid dataset
 *
 */
void Combine_EGSPHANT(const char *ct_filename, const char *epid_filename, const char *output_filename)
{
	int i, j, k;
	int airPadding=0;

	// Read in the phantom data
	EGSPhant ct_data((char *)ct_filename);
	int result1 = ct_data.Read();

	// Read in the EPID data
	EGSPhant epid_data((char *)epid_filename);
	int result2 = epid_data.Read();

	//number of rows of air required to get the desired EDIP distance
	airPadding = (int)floor(EPID_dist-fabs(ct_data.yBoundaries[ct_data.ySize]-iy));

	// Need to ensure that these data may be combined - check dimensions
	if (ct_data.xSize > epid_data.xSize )
	{
		fprintf(stdout, "Cannot be combined: Check number of voxels in Z dimension.\n");
		exit(1);
	}

	if (ct_data.zSize != epid_data.zSize )
	{
		fprintf(stdout, "Cannot be combined: Check number of voxels in Z dimension.\n");
		exit(1);
	}

	/*
	 *	Create the new EGSPhant file
	 */

	// Initialise the combinatory dataset
	EGSPhant combine((char *)output_filename);
	
	// Set the number of Media	
	combine.numberOfMedia = epid_data.numberOfMedia;

	// Set the names of the Media
	combine.mediaNames = new std::string[epid_data.numberOfMedia];

	for (i = 0; i<epid_data.numberOfMedia; i++)
	{
		combine.mediaNames[i] = epid_data.mediaNames[i];
	}

	// Initialise and read in the ESTEPE
	combine.estepe = new float[combine.numberOfMedia];

	for (i = 0; i<combine.numberOfMedia; i++)
	{
		combine.estepe[i] = epid_data.estepe[i];
	}

	// Set the number of Voxels in each dimension
	combine.xSize = epid_data.xSize;
	combine.ySize = ct_data.ySize + epid_data.ySize + airPadding;
	combine.zSize = epid_data.zSize;

	/*
	*	Populate the Boundaries
	*/

	// X dimension - no changes from epid_datas
	combine.xBoundaries = new float[epid_data.xSize + 1 ];

	for (i = 0; i<( epid_data.xSize + 1 ); i++)
	{
		combine.xBoundaries[i] = epid_data.xBoundaries[i];
	}

	// Y dimension - EPID followed by CT
	combine.yBoundaries = new float[epid_data.ySize +airPadding+ ct_data.ySize+1];

	for (i = 0; i<(ct_data.ySize+1); i++)
	{
		combine.yBoundaries[i] = ct_data.yBoundaries[i];
	}

	for (i=1;i<=airPadding;i++)
	{
		combine.yBoundaries[ct_data.ySize + i] = ct_data.yBoundaries[ct_data.ySize] + i*vy;
	}
		
	for (i = 1; i<(epid_data.ySize+1); i++)
	{
		combine.yBoundaries[ct_data.ySize + airPadding + i] = ct_data.yBoundaries[ct_data.ySize + airPadding]+i*vy;
	}

	// Z dimension - No Change
	combine.zBoundaries = new float[ct_data.zSize+1];

	for (i = 0; i<(ct_data.zSize+1); i++)
	{
		combine.zBoundaries[i] = ct_data.zBoundaries[i];
	}
	
	/*
	*	Create the voxel medium numbers 
	*/
	combine.voxelMedium = new char**[combine.zSize];

	for (k = 0; k<combine.zSize; k++)
	{
		combine.voxelMedium[k] = new char*[combine.ySize];

		for (j = 0; j<combine.ySize; j++)
		{
			combine.voxelMedium[k][j] = new char[combine.xSize];
		}
	}

	// the column in which the CT data starts 
	int dataStart = (int)floor((epid_data.xSize - ct_data.xSize)/2);
	
	// Start Initialising
	for (k = 0; k<combine.zSize; k++)
	{
		for (j = 0; j<ct_data.ySize; j++)
		{
			for (i = 0; i<epid_data.xSize; i++)
			{
				//only add the data within its region
				if (i>=dataStart && i<(ct_data.xSize+dataStart))
				{
					combine.voxelMedium[k][j][i] = ct_data.voxelMedium[k][j][i-dataStart];
				}
				else
				{
					//pad with air
					combine.voxelMedium[k][j][i] = '1';
				}
			}
		}

		//pad the gap between the EPID and CT with air
		for (j=0;j<airPadding;j++)
		{
			for (i=0;i<epid_data.xSize;i++)
			{
				combine.voxelMedium[k][ct_data.ySize +j][i] = '1';
			}
		}
	
		//write the EPID
		for (j = 0; j<epid_data.ySize; j++)
		{
			for (i = 0; i<epid_data.xSize; i++)
			{
				combine.voxelMedium[k][ct_data.ySize + airPadding + j][i] = epid_data.voxelMedium[k][j][i];
			}
		}
	}

	/*
	*	Create the voxel densities 
	*/
	combine.voxelDensity = new float**[combine.zSize];

	for (k = 0; k<combine.zSize; k++)
	{
		combine.voxelDensity[k] = new float*[combine.ySize];

		for (j = 0; j<combine.ySize; j++)
		{
			combine.voxelDensity[k][j] = new float[combine.xSize];
		}
	}

	// Start Initialising
	for (k = 0; k<combine.zSize; k++)
	{
		for (j = 0; j<ct_data.ySize; j++)
		{
			for (i = 0; i<ct_data.xSize; i++)
			{
				if (i>=dataStart && i<(ct_data.xSize+dataStart))
				{
					combine.voxelDensity[k][j][i] = ct_data.voxelDensity[k][j][i-dataStart];
				}
				else
				{
					//pad with air
					combine.voxelDensity[k][j][i] = 0.001;
				}
			}
		}
				
		for (j=0;j<airPadding;j++)
		{
			for (i=0;i<epid_data.xSize;i++)
			{
				combine.voxelDensity[k][ct_data.ySize +j][i] = 0.001;
			}
		}
	
		for (j = 0; j<epid_data.ySize; j++)
		{
			for (i = 0; i<epid_data.xSize; i++)
			{
				combine.voxelDensity[k][ct_data.ySize + airPadding + j][i] = epid_data.voxelDensity[k][j][i];
			}
		}
	}

	// Write Out the data file
	combine.Write();
}


/*
* Converts DICOM data as extracted from a DICOMreader object and converts it to the voxelised
* EGSPhant format, incudes EPID data generation
*
*/
void ConvertDICOMToEGSPhant(DICOMReader* Dicom,EGSPhant* EGS){

	
	float rescale = 0.1;          //DICOM data coords are in mm and EGSPhant is in cm
	int i,j,k;                    //loop counters
	int idx,idy,idz;              //indices used as starting point for interpolation
	float lr,ls,lt,r,s,t;         // weights used for interpolation
	float dx,dy,dz;               //voxel size in terms of DICOM data grid indices
	
	float total,CTavg;            // used to calculate the voxel category
	float lowerBound;
	int nonZero,dataStart;
	
	unsigned w;                   //loop counter that we also do some bit-wise stuff with
	
	int airPadding;               //all these are used for patitioning volumn into data and air fill-in zones
	int EPIDVoxels=0; 
	char temp[255];
	int a,b;
	int startLayer;
	int EPIDtop;

	// get the depth of the EPID in voxels
	for(i=0;i<totalEPIDLayers;i++)
	{
		EPIDVoxels += (int)ceil(EPIDLayerInfo[i][2]/vy);
	}

	//get the number of voxels between the data volume and the EPID surface 
	if (EPID_dist)
	{
		airPadding = (int)floor((EPID_dist - fabs(Dicom->yCoords[Dicom->height-1]*rescale - iy))/vy);
	}
	else
	{
		airPadding =0;
	}

	// get the number of columns befor data starts
	dataStart = (int)floor(((EPIDx - fabs(Dicom->xCoords[0] - Dicom->xCoords[Dicom->width -1])*rescale)/vx)*0.5 );
	
	//calculate the number of voxels
	EGS->xSize = (int)floor(EPIDx/vx); // x dimension of the EPID
	EGS->ySize = (int)floor(((fabs(Dicom->yCoords[0] - Dicom->yCoords[Dicom->height -1])*rescale))/vy)+airPadding+EPIDVoxels;
	EGS->zSize = (int)floor((fabs(Dicom->zCoords[0] - Dicom->zCoords[Dicom->numSlices -1])*rescale)/vz);

	// print out some grid stats
	printf("Original Volume Dimensions: %d %d %d\n",Dicom->width,Dicom->height,Dicom->numSlices);
	printf("Voxel Dimesions %f %f %f\n",vx,vy,vz);
	printf("Voxelised Grid Dimensions with EPID: %d %d %d\n",EGS->xSize,EGS->ySize,EGS->zSize);

	// map the size of voxels in cms to size in DICOM data grid indicies
	dx = (vx*(Dicom->width-1))/(fabs(Dicom->xCoords[0] - Dicom->xCoords[Dicom->width -1])*rescale);
	dy = (vy*(Dicom->height-1))/(fabs(Dicom->yCoords[0] - Dicom->yCoords[Dicom->height -1])*rescale);
	dz = (vy*(Dicom->numSlices-1))/(fabs(Dicom->zCoords[0] - Dicom->zCoords[Dicom->numSlices -1])*rescale);

	//include CT medium summary
	EGS->numberOfMedia = EPIDnumMedia;
	EGS->mediaNames = EPIDMaterials;

	//initialise ESTEP array
	EGS->estepe = new float[EPIDnumMedia];

	//initialise boundary arrays
	EGS->xBoundaries = new float[EGS->xSize];
	EGS->yBoundaries = new float[EGS->ySize];
	EGS->zBoundaries = new float[EGS->zSize];

	//initialise the EGS medium and density grids
	EGS->voxelMedium = new char**[EGS->zSize];
	EGS->voxelDensity = new float**[EGS->zSize];

	for (k=0; k<EGS->zSize; k++)
	{
		EGS->voxelMedium[k] = new char*[EGS->ySize];
		EGS->voxelDensity[k] = new float*[EGS->ySize];

		for (j=0; j<EGS->ySize; j++)
		{	
			EGS->voxelMedium[k][j] = new char[EGS->xSize];
			EGS->voxelDensity[k][j] = new float[EGS->xSize];
		}
	}

	//store the ESTEP
	for (i=0;i<numMedia;i++)
	{
		EGS->estepe[i] = conversionData[i][3];
	}

	for (i=numMedia;i<EPIDnumMedia;i++)
	{
		EGS->estepe[i] = 1;
	}

	//for each voxel
	for(k=0;k<EGS->zSize;k++)
	{
		//the section where the CT data will be written
		for(j=0;j<EGS->ySize-airPadding-EPIDVoxels;j++)
		{
			//pad with air in columes before data
			for(i=0;i<dataStart;i++)
			{
				EGS->voxelMedium[k][j][i] = '1';
				EGS->voxelDensity[k][j][i] = AIRDENSE;
			}//end for i

			//categorise all voxels in data grid and interpolate the densities
			for(i=dataStart;i<EGS->xSize-dataStart;i++)
			{
				total =0; //sum the data values on the voxel vertices
				nonZero = 0; // keeps track of non-air vertices for weighting

				/* the verticies may not be on DICOM grid points, 
				 * so we must interpolate the data on the boundaries
				 */
				for (w=0;w<8;w++)
				{	
					// use the voxel indecies to map to the indices of the surrounding DICOM grid points 
					lr = (float)((i-dataStart+(w & 1)))*dx;
					ls = (float)((j+((w & 2)>>1)))*dy;
					lt = (float)((k+((w & 4)>>2)))*dz;
					idx = (int)floor(lr);
					idy = (int)floor(ls);
					idz = (int)floor(lt);
					r = lr - idx;
					s = ls - idy;
					t = lt - idz;
		
					//if the point is within grid bounds form a weighted sum of the surrounding points 
					//(just summed to the total since we where going to add them anyway)
					if (idx >= 0 && idy >= 0 && idz >= 0 && (idx < Dicom->width-1) && (idy < Dicom->height-1) && (idz < Dicom->numSlices-1))
					{	
						//this section is pretty inefficient, but it does work... 
						//Note-to-self: optimise this (or at least make less crap)
						if(Dicom->data_s[idz][idy][idx]>0)
							total += (1-r)*(1-s)*(1-t)*((float)Dicom->data_s[idz][idy][idx]);
						if(Dicom->data_s[idz][idy][idx+1]>0)
							total += (r)*(1-s)*(1-t)*((float)Dicom->data_s[idz][idy][idx+1]);
						if(Dicom->data_s[idz][idy+1][idx]>0)
							total += (1-r)*(s)*(1-t)*((float)Dicom->data_s[idz][idy+1][idx]);
						if(Dicom->data_s[idz+1][idy][idx]>0)
							total += (1-r)*(1-s)*(t)*((float)Dicom->data_s[idz+1][idy][idx]);
						if(Dicom->data_s[idz][idy+1][idx+1])
							total += (r)*(s)*(1-t)*((float)Dicom->data_s[idz][idy+1][idx+1]);
						if(Dicom->data_s[idz+1][idy][idx+1])
							total += (r)*(1-s)*(t)*((float)Dicom->data_s[idz+1][idy][idx+1]);
						if(Dicom->data_s[idz+1][idy+1][idx]>0)
							total += (1-r)*(s)*(t)*((float)Dicom->data_s[idz+1][idy+1][idx]);
						if(Dicom->data_s[idz+1][idy+1][idx+1]>0)
							total += (r)*(s)*(t)*((float)Dicom->data_s[idz+1][idy+1][idx+1]);
					}

					if(total)
					{   //keep a count of non-air boundaries
						nonZero++;
					}
				}//end for w
				
				//for now just assign the voxel the CT number of the mid point as interpolated by the boundaries
				if (nonZero)
				{
					CTavg = total/nonZero;
				}
				else // if its all air, just set to zero
				{
					CTavg =0;
				}

				//categorise Voxel
				lowerBound =0;

				for(w=0;w<numMedia;w++)
				{
					//if its not in this range of this category, then we move on to the nex one
					if( !(CTavg >= lowerBound && CTavg < conversionData[w][0]) && !(CTavg < 0) && !(CTavg > conversionData[numMedia-1][0]))
					{
						lowerBound = conversionData[w][0];
					}
					else //we've found the category and we stop
					{
						break;
					}
				}
				
				//write classification of voxel
				EGS->voxelMedium[k][j][i] = mediaCats[w];
				
				//if we are less than zero than set to lower bound of AIR
				if (CTavg < 0)
				{					
					EGS->voxelDensity[k][j][i] = conversionData[0][1];
				}
				else //otherwise use linear interpolation for the density based on conversion data
				{
					EGS->voxelDensity[k][j][i] = ((conversionData[w][2]-conversionData[w][1])*(CTavg - lowerBound))/(conversionData[w][0]-lowerBound) + conversionData[w][1];
				}
			}//end for i

			//pad with air the columns following the data 
			for(i=EGS->xSize-dataStart; i<EGS->xSize; i++)
			{
				EGS->voxelMedium[k][j][i] = '1';
				EGS->voxelDensity[k][j][i] = AIRDENSE;
			}//end for i
		}//end for j

		// Pad the Rows between the Volume and the and EPID with air
		for(j=EGS->ySize-airPadding-EPIDVoxels; j<EGS->ySize-EPIDVoxels; j++)
		{
			//theses layers are all air
			for(i=0; i<EGS->xSize; i++)
			{
				EGS->voxelMedium[k][j][i] = '1';
				EGS->voxelDensity[k][j][i] = AIRDENSE;
			}
		}

		//now for the EPID

		EPIDtop = EGS->ySize-EPIDVoxels;
		a=0;
		b=0;
		startLayer =0;

		//for each EPID Layer
		for(w=0;w<totalEPIDLayers;w++)
		{
			//we get the start and end of this Layer in voxels
			a = startLayer;
			b = startLayer+(int)ceil(EPIDLayerInfo[w][2]/vy);

			//for the rows in this Layer
			for(j=EPIDtop+a; j<EPIDtop+b; j++)
			{
				// write the Material Category
				for(i=0; i<EGS->xSize; i++)
				{
					sprintf(temp,"%d",(int)EPIDLayerInfo[w][0]);
					EGS->voxelMedium[k][j][i] = temp[0];
					EGS->voxelDensity[k][j][i] = EPIDLayerInfo[w][1];
				}
			}

			//set up for next Layer
			startLayer = b; 
		}
		
	}

	// now just write the boundary points (we can just calculate 'em using the start coords and voxel dimensions)
	for (i=0;i<=EGS->xSize;i++)
	{
		EGS->xBoundaries[i] = (EPIDx*0.5) - i*vx;
	}
	for (j=0;j<=EGS->ySize;j++)
	{
		EGS->yBoundaries[j] = Dicom->yCoords[0]*rescale + j*vy;
	}

	for (k=0;k<=EGS->zSize;k++)
	{
		EGS->zBoundaries[k] = Dicom->zCoords[0]*rescale - k*vz;
	}
}

//Just a simple function that was used for debugging purposes
void ReadDicoms_Test( const char *directory,float deg, int slice)
{
	float axis[3] = {0,0,1};
	
	DICOMReader reader((char *)directory);
	reader.ReadSample();
	reader.ReadAll();

	reader.WriteTestFile();
	reader.Select_Region(-10,10,-100,-90,-10,10);
	reader.PrintAllFiles();

	reader.OutputSliceToFile(slice,"Pre-rotation");
	reader.Rotate3D(deg,axis,ix,iy,iz,flag);
	printf("rotated and angle of %f\n",deg);
	reader.OutputSliceToFile(slice,"Post-rotation");
	

}

// Another Test function used for initial testing of rotation
void TestEGSPhantRotate(const char* file, float deg,const char* outfile)
{	
	float axis[3] = {0,0,1};
	EGSPhant data((char*)file);
	data.Read();
	data.Rotate3D(deg,axis);
	data.Write();
}


// Reader for the EPID info file used to generate EPID Layers in ConvertDICOMtoEGSPhant()
int ReadEPIDspec(const char* filename){

	FILE* file;
	int i;
	char tempName[255];

	// check we can access the file to write
	if(!(file = fopen(filename,"r"))){
		printf("ERROR: could not open EPID data file\n");
		return 1;
	}

	// read the x and z dimensions
	if(!(fscanf(file,"%f %f\n",&EPIDx,&EPIDz)==2))
	{
		printf("ERROR: could not read EPID xz dimensions\n");
		return 1;
	}

	// read number of materials in this file
	if(!(fscanf(file,"%d\n",&EPIDnumMedia)==1))
	{
		printf("ERROR: could not read EPID Material Number\n");
		return 1;
	}

	//initialise materials array
	EPIDMaterials = new std::string[EPIDnumMedia];

	// read in the material names
	for (i=0;i<EPIDnumMedia;i++)
	{
		if(!(fscanf(file,"%s",tempName)==1))
		{
			printf("ERROR: could not read Material %d\n",i);
			return 1;
		}

		EPIDMaterials[i] = tempName;
	}

	// read the number of EPID Layers
	if(!(fscanf(file,"%d\n",&totalEPIDLayers)==1))
	{
		printf("ERROR: could not read number of Layers\n");
		return 1;
	}

	//initailise EPID layer data
	EPIDLayerInfo = new float*[totalEPIDLayers];

	for(i=0;i<totalEPIDLayers;i++)
	{
		EPIDLayerInfo[i] = new float[3];
	}

	// read layer info
	for(i=0;i<totalEPIDLayers;i++)
	{
		if(!(fscanf(file,"%f %f %f",&(EPIDLayerInfo[i][0]),&(EPIDLayerInfo[i][1]),&(EPIDLayerInfo[i][2]))==3))
		{
			printf("ERROR: could not read layer info for layer %d",i);
			return 1;
		}
	}

	return 0;
}

// reads the .inp input file that specifies parameters for conversion to EGSPhant format
int ReadINP(const char* filename)
{
	FILE* file;
	int i;
	char temp[255];
	char s[255];

	// check we can open the file
	if(!(file = fopen(filename,"r")))
	{
		printf("ERROR: could not open .inp data file\n");
		return 1;
	}
	
	//skip the first two lines
	if(!(fscanf(file,"%s\n",temp)==1))
	{
		printf("ERROR: could not read data format\n");
		return 1;
	}

	// check the data format
	if(strcmp(temp,"DICOM"))
	{
		printf("ERROR: Only DICOM data is supported by CTCombine %f\n",VERSION);
		return 1;
	}

	//read DICOM directory
	if(!(fscanf(file,"%s\n",dir)==1))
	{
		printf("ERROR: could not read DICOM directory\n");
		return 1;
	}

	// read data volume bounds
	if(!(fscanf(file,"%f, %f, %f, %f, %f, %f\n",&xmin,&xmax,&ymin,&ymax,&zmin,&zmax)==6))
	{
		printf("ERROR: could not read data boundaries\n");
		return 1;
	}

	// read voxel dimensions
	if(!(fscanf(file,"%f, %f, %f\n",&vx,&vy,&vz)==3))
	{
		printf("ERROR: could not read voxel dimensions\n");
		return 1;
	}

	// read number of media
	if(!(fscanf(file,"%d\n",&numMedia)==1))
	{
		printf("ERROR: could not read number of Materials in CT\n");
		return 1;
	}

	// intialise the Material Name,Category, and conversion data arrays
	MediaNames = new std::string[numMedia];

	mediaCats = new char[numMedia];

	conversionData = new float*[numMedia];
	for (i=0; i<numMedia; i++)
	{
		conversionData[i] = new float[4];
	}

	// now read in the relevant data
	for (i=0;i<numMedia;i++)
	{
		// read the media name 
		if (!(fscanf(file,"%s",temp)==1))
		{
			printf("ERROR: could not read Material name %d\n",i);
			return 1;
		}

		MediaNames[i] = temp;

		//assign a category number
		sprintf(s,"%d",i+1);
		mediaCats[i] = s[0]; 

		// read the conversion data
		if (!(fscanf(file,"%f,%f,%f,%f\n",&(conversionData[i][0]),&(conversionData[i][1]),&(conversionData[i][2]),&(conversionData[i][3]))==4))
		{
			printf("ERROR: could not read CT to density conversion data\n");
			return 1;
		}
	}

	return 0;
}

// just dumps an options menu
int printHelp(){
	printf("CTCombine: version %f\n",VERSION);
	printf("\n");
	printf("The Following switches are used to incicate the data being input.\n");
	printf("all of the following parameters are required:\n");
	printf("\n");
	printf("-E filename              the following parameter is the name of the EPID spec file\n");
	printf("-o outfile               the following parameter is the name of the output.egsphant file\n");
	printf("-inp filename            the following is the file name of the .inp input file\n");
	printf("\n");
	printf("The following are optional:");
	printf("\n");
	printf("-d dir                   the following parameter will override the DICOM directory in the .inp file\n");
	printf("-i x y z                 indicates the desired iso-centre\n");
	printf("-e dist                  indicates the desired distance from the isocentre to the top of the EPID\n");
	printf("-g theta phi             indicates the theta and phi angle of zero gantry angle\n");
	printf("-r theta phi             indicates the theta and phi angle of desired gantry angle\n");
	printf("-c file1 file2 file3     Combines the following two files (file1,file2) into the one output file indicayed by file3");
}

// reads the commandline arguments passed to CTCombine
int ReadArgs(int argc, char* argv[]){

	int i;
	int temp;
	bool EPIDinit=0;
	bool OUTinit=0;
	bool INPinit=0;
	bool overrideDIR=0;


	// check there are some parameters
	if(argc>1)
	{
		//need to put check in to make sure the correct numbe of arguments follow each switch
		for(i=1; i<argc; i++)
		{
			if(!strcmp(argv[i],"-d"))
			{
				sprintf(dir,"%s",argv[++i]);
				temp = i;
				overrideDIR = 1;
			}
			else if (!strcmp(argv[i],"-inp"))//indicates INP input data file
			{
				if(ReadINP(argv[++i]))
				{
					printf("ERROR: Bad read of file %s\n",argv[i]);
					return 1;
				}
				INPinit = 1;
			}
			else if (!strcmp(argv[i],"-i"))//indicates isocentre
			{
				ix = atof(argv[++i]);
				iy = atof(argv[++i]);
				iz = atof(argv[++i]);

				flag = 1;
			}
			else if (!strcmp(argv[i],"-g"))//indicates direction of zero gantry angle
			{
				Z_gantry[0] = atof(argv[++i]);//theta angle
				Z_gantry[1] = atof(argv[++i]);//phi angle
			}
			else if (!strcmp(argv[i],"-r")) //indicates theta and phi rotation angles
			{
				theta = atof(argv[++i]);
				phi = atof(argv[++i]);
			}
			else if (!strcmp(argv[i],"-e")) //indicates distances desired from isocentre to EPID data
			{
				EPID_dist = atof(argv[++i]);
			}
			else if (!strcmp(argv[i],"-c")) //indicates CT number to density convesion data
			{
				Combine_EGSPHANT(argv[i+1],argv[i+2],argv[i+3]);
				return 0;
			}
			else if (!strcmp(argv[i],"-E")) //indicates EPID data file
			{
				if(ReadEPIDspec(argv[++i]))
				{
					printf("ERROR: Bad read of file %s",argv[i]);
					return 1;
				}
				EPIDinit = 1;
			}
			else if (!strcmp(argv[i],"-o")) //indicates outputfile name
			{
				sprintf(outputfileName,"%s",argv[++i]);
				OUTinit = 1;
			}
			else
			{
				printHelp();
				return 1;
			}
		}
	}
	else
	{
		printHelp();
		return 1;
	}

	// check the require arguments have been included
	if (EPIDinit && OUTinit && INPinit)
	{
		// check if the user used the -d switch
		if (overrideDIR)
		{
			sprintf(dir,"%s",argv[temp]);
		}
		return 0;
	}
	else
	{
		printf("ERROR: one of the required files is not specified\n");
		printHelp();
	}
}


int main(int argc, char* argv[])
{
	int i; //for loops
	int scaleFactor = 10; //to convert from cm to mm
	float zaxis[3] = {0,0,1}; // axis of rotation for theta
	float xaxis[3] = {1,0,0}; // axis of rotation for phi
	
	// read commandline arguments
	if(ReadArgs(argc,argv))
	{
		printf("Program cannot complete due to invalid input\n");
		return 1;
	}
	else
	{
		// read the DICOM data
		printf("DICOM Directory: %s\n",dir);
		DICOMReader DICOM(dir); //initialise the DICOM reader
		printf("Reading DICOM data...\n");
		DICOM.ReadSample();
		DICOM.ReadAll();
		printf("Read Complete...\n");

		printf("Rotating for new gantry angle...\n");
		//make roation from zero gantry angle to the new gantry angle 
		DICOM.Rotate3D(-Z_gantry[0] + theta,zaxis,ix*scaleFactor,iy*scaleFactor,iz*scaleFactor,flag);
		DICOM.Rotate3D(-Z_gantry[1] + phi,xaxis,ix*scaleFactor,iy*scaleFactor,iz*scaleFactor,flag);
		printf("Rotation Complete...\n");
		// Select from the rotated data the desired data region
		printf("Setting bounds...\n");
		if(DICOM.Select_Region(xmin*10,xmax*scaleFactor,ymin*scaleFactor,ymax*scaleFactor,zmin*scaleFactor,zmax*scaleFactor))
		{
			printf("Could not complete operation exiting...\n");
			return 1;
		}
		// translate to new isocentre
		
		printf("Converting to .egsphant...\n");
		//  now initailise EGSPhant object to store converted data
		EGSPhant EGSct(outputfileName);
		// convert
		ConvertDICOMToEGSPhant(&DICOM,&EGSct);
		printf("data converted...\n");
		printf("Translating to isocentre...\n");
		EGSct.ShiftToOrigin(ix,iy,iz,flag);
		printf("Flipping y axis...\n");
		EGSct.FlipY();
		printf("Writing EGS file...\n");
		// hmm... what does the next line do? Oh yes, thats right... we write to a file
		EGSct.Write();
	}

	/*
	* The Following is all just test code calls that I dont want to delete
	*
	*/

	//just a test of the .inp reader
	//if (argc ==4){

	//	std::string INPfilename(argv[1]);

	//	if(ReadINP(INPfilename.c_str())){
	//		printf("ERROR: cannot read file...\n");
	//		return 1;
	//	}
	//	else{
	//		theta = atof(argv[2]);
	//		phi = atof(argv[3]);
	//		std::string dirname = "./Data/cthead/";
	//		EGSPhant EGS1("test1.egsphant");
	//		DICOMReader reader((char *)dirname.c_str());
	//		reader.ReadSample();
	//		reader.ReadAll();
	//		printf("Rotating...\n");
	//		
	//		reader.Rotate3D(90 - theta,zaxis);
	//		reader.Rotate3D(270 - phi,xaxis);
	//		
	//		printf("Rotation complete...\n");
	//		printf("converting..\n");
	//		ConvertDICOMToEGSPhant(&reader,&EGS1);
	//		printf("file converted...\n");
	//		EGS1.Write();
	//		/*printf("Rotating...\n");
	//		
	//		printf("Rotation complete...\n");
	//		printf("Converting second file...\n");
	//		ConvertDICOMToEGSPhant(&reader,&EGS2);
	//		printf("Second file Converted");
	//		EGS2.Write();*/

	//		return 0;
	//	}

	//}

	//
	//// Count the command line arguements
	//if (argc < 4 || argc > 4)
	//{
	//	fprintf(stdout, "Wrong number of command line arguements\n");
	//	fprintf(stdout, "Usage: CTCombine phantom.egsphant epid.egsphant someoutputfile.egsphant\n");
	//}
	//else
	//{
	//	std::string ct_filename(argv[1]);
	//	std::string epid_filename(argv[2]);
	//	std::string output_filename(argv[3]);	
	//	Combine_EGSPHANT(ct_filename.c_str(), epid_filename.c_str(), output_filename.c_str());
	//}


	//std::string ct_filename = "./Data/Second/cthead.egsphant";
	//std::string epid_filename = "./Data/Second/EPID_for_cthead.egsphant";	
	//Combine_EGSPHANT(ct_filename.c_str(), epid_filename.c_str(), output_filename.c_str());
	//if(argc == 3){
	//	float deg = atof(argv[1]);
	//	int slc = atoi(argv[2]);

	//	std::string dirname = "./Data/cthead/";
		//std::string testfile = "cthead.egsphant";
		//std::string outfile = "cthead_rotated.egsphant";
		//printf("Testing rotation...\n");
		//TestEGSPhantRotate(testfile.c_str(),deg,outfile.c_str());
	//	ReadDicoms_Test(dirname.c_str(),deg,slc);
	//}	
	return 0;
}



