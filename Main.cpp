/**
 * Main.cpp (how descriptive)
 *
 * @Authors: Mark Dwyer, David Warne
 * @Contact: m2.dwyer@qut.edu.au, david.warne@qut.edu.au
 * @Last Modified: 19/03/2009
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
#define VERSION 0.18
#define HITS_MAX 100

template<class T>
void byteswap(T* t)
{
	std::reverse(reinterpret_cast<char*>(t), reinterpret_cast<char*>(t + 1));
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
int numMedia;                          // total materials in CT
std::string* MediaNames;               // names of said materials
char* mediaCats;                       // category numbers (as ASCII char) of CT materials
float** conversionData;                // stores the conversion data
int* numControlPoints;                 // total control points for each density transfer function
float* estepes;                        // to store the ESTEPEs (whatever they are)
float*** densityTransfers;             // density transfer functions (one per material)


/*
 * COMBINE_EGSPHANT: Combine the CT with the EPID data, is probably redundant now we create 
 *                   EPID data on conversion From Dicom
 *
 * Parameters: 
 *            ct_filename     - the phantom dataset
 * 	          epid_filename   - the epid dataset
 *            output_filename - the name of the output file
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
	airPadding = (int)round(EPID_dist-fabs(ct_data.yBoundaries[ct_data.ySize]-iy));

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
	int xdataStart = (int)round((epid_data.xSize - ct_data.xSize)/2);
	
	// Start Initialising
	for (k = 0; k<combine.zSize; k++)
	{
		for (j = 0; j<ct_data.ySize; j++)
		{
			for (i = 0; i<epid_data.xSize; i++)
			{
				//only add the data within its region
				if (i>=xdataStart && i<(ct_data.xSize+xdataStart))
				{
					combine.voxelMedium[k][j][i] = ct_data.voxelMedium[k][j][i-xdataStart];
				}
				else
				{
					//pad with air
					combine.voxelMedium[k][j][i] = mediaCats[0];
				}
			}
		}

		//pad the gap between the EPID and CT with air
		for (j=0;j<airPadding;j++)
		{
			for (i=0;i<epid_data.xSize;i++)
			{
				combine.voxelMedium[k][ct_data.ySize +j][i] = mediaCats[0];
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
				if (i>=xdataStart && i<(ct_data.xSize+xdataStart))
				{
					combine.voxelDensity[k][j][i] = ct_data.voxelDensity[k][j][i-xdataStart];
				}
				else
				{
					//pad with air
					combine.voxelDensity[k][j][i] = densityTransfers[0][0][1];
				}
			}
		}
				
		for (j=0;j<airPadding;j++)
		{
			for (i=0;i<epid_data.xSize;i++)
			{
				combine.voxelDensity[k][ct_data.ySize +j][i] = densityTransfers[0][0][1];
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
* CONVERTDICOMTOEGSPHANT: Converts DICOM data as extracted from a DICOMreader object and converts it to the 
*                         voxelised EGSPhant format, incudes EPID data generation 
*
* Parameters:
*                DICOM - pointer to DicomReader object
*                EGS   - pointer to EGSPhant object
*
* Pre-Condition: for the conversion to work is that the coordinates have already been
*       translated such that the iso-centre is now the origin
*/
void ConvertDICOMToEGSPhant(DICOMReader* Dicom,EGSPhant* EGS){

	
	float rescale = 0.1;         //DICOM data coords are in mm and EGSPhant is in cm
	int i,j,k;                   //loop counters
	int idx,idy,idz;             //indices used as starting point for interpolation
	float CTavg;                 // used to calculate the voxel category
	float lowerBound;
	unsigned w,c;                //loop counter that we also do some bit-wise stuff with	
	int airPadding;              //all these are used for patitioning volumn into data and air fill-in zones
	int EPIDVoxels=0;
	char temp[255];
	int a,b;
	int startLayer;
	int EPIDtop;
	int xVolVoxels;
	int yVolVoxels;
	int zVolVoxels;
	int*** hits;

	// get the depth of the EPID in voxels
	for(i=0;i<totalEPIDLayers;i++)
	{
		EPIDVoxels += (int)ceil(EPIDLayerInfo[i][2]/vy);
	}

	//get the number of voxels between the data volume and the EPID surface 
	airPadding = (int)round((EPID_dist - Dicom->ylims[1]*rescale)/vy);
	
	if (airPadding < 0)// negatives are un-physical
	{
	    airPadding =0;
	}
	
	xVolVoxels = (int)round(fabs(Dicom->xlims[0] - Dicom->xlims[1])*rescale)/vx;
	yVolVoxels = (int)round(fabs(Dicom->ylims[0] - Dicom->ylims[1])*rescale)/vy;
	zVolVoxels = (int)round(fabs(Dicom->zlims[0] - Dicom->zlims[1])*rescale)/vz;
	
	//calculate the number of voxels
	EGS->xSize = (int)round(EPIDx/vx); // x dimension of the EPID
	EGS->ySize = (int)round(((fabs(Dicom->ylims[0] - Dicom->ylims[1])*rescale))/vy)+airPadding+EPIDVoxels;
	EGS->zSize = (int)round(EPIDz/vz); // z dimesion of EPID  
 	// print out some grid stats
	printf("Original Volume Dimensions: %d %d %d\n",Dicom->width,Dicom->height,Dicom->numSlices);
	printf("Voxel Dimesions %f %f %f\n",vx,vy,vz);
	printf("Voxelised Grid Dimensions with EPID: %d %d %d\n",EGS->xSize,EGS->ySize,EGS->zSize);

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
	hits = new int**[EGS->zSize];

	for (k=0; k<EGS->zSize; k++)
	{
		EGS->voxelMedium[k] = new char*[EGS->ySize];
		EGS->voxelDensity[k] = new float*[EGS->ySize];
		hits[k] = new int*[EGS->ySize];

		for (j=0; j<EGS->ySize; j++)
		{	
			EGS->voxelMedium[k][j] = new char[EGS->xSize];
			EGS->voxelDensity[k][j] = new float[EGS->xSize];
			hits[k][j] = new int[EGS->xSize];
			for (i=0;i<EGS->xSize;i++)
			{
			    EGS->voxelDensity[k][j][i] = densityTransfers[0][0][1];
			    EGS->voxelMedium[k][j][i] = mediaCats[0];
			    hits[k][j][i] = 0;
			}
		}
	}

	//store the ESTEP
	for (i=0;i<numMedia;i++)
	{
		EGS->estepe[i] = estepes[i];
	}

	for (i=numMedia;i<EPIDnumMedia;i++)
	{
		EGS->estepe[i] = 1;
	}
	
	if (Dicom->imageDataType ==2) //short
	{
	    // now map the Dicom Data to the Dicom volume
	    for (k=0;k<Dicom->numSlices;k++)
	    {
	        for(j=0;j<Dicom->height;j++)
	        {
	            for(i=0;i<Dicom->width;i++)
	            {
	               
	                // map to a voxel
	                idx = (int)round(((EPIDx*0.5 + Dicom->xCoords[k][j][i]*rescale)/EPIDx)*EGS->xSize);
	                idy = (int)round((Dicom->yCoords[k][j][i]-Dicom->ylims[0])/(Dicom->ylims[1]-Dicom->ylims[0])*(yVolVoxels));
	                idz = (int)round(((EPIDz*0.5 - Dicom->zCoords[k][j][i]*rescale)/EPIDz)*EGS->zSize);;
	                
	                if (idx >=0 && idx < EGS->xSize && idz >= 0 && idz < EGS->zSize)
	                {
	                    // add the CT value to density array to save storage
	                    EGS->voxelDensity[idz][idy][idx] +=(float) Dicom->data_s[k][j][i];
	                   
	                    //increment the char in the Hits array (might be better to use medium array)
	                    if(hits[idz][idy][idx] < HITS_MAX)
	                    {
	                        hits[idz][idy][idx]++;
	                    }
	                }
	               
	            }
	        }
	    }
	}
	else if (Dicom->imageDataType == 0)//float
	{
	    // now map the Dicom Data to the Dicom volume
	    for (k=0;k<Dicom->numSlices;k++)
	    {
	        for(j=0;j<Dicom->height;j++)
	        {
	            for(i=0;i<Dicom->width;i++)
	            {
	                // map to a voxel
	                idx = (int)round(((EPIDx*0.5 + Dicom->xCoords[k][j][i]*rescale)/EPIDx)*EGS->xSize);
	                idy = (int)round((Dicom->yCoords[k][j][i]-Dicom->ylims[0])/(Dicom->ylims[1]-Dicom->ylims[0])*(yVolVoxels));
	                idz = (int)round(((EPIDz*0.5 - Dicom->zCoords[k][j][i]*rescale)/EPIDz)*EGS->zSize);;
	                
	                if (idx >=0 && idx < EGS->xSize && idz >= 0 && idz < EGS->zSize)
	                {
	                    // add the CT value to density array to save storage
	                    EGS->voxelDensity[idz][idy][idx] +=(float) Dicom->data_f[k][j][i];
	                   
	                    //increment the char in the Hits array (might be better to use medium array)
	                    if(hits[idz][idy][idx] < HITS_MAX)
	                    {
	                        hits[idz][idy][idx]++;
	                    }
	                }
	            }
	        }
	    }
	}
	else if(Dicom->imageDataType == 1) // unsigned char
	{
	    // now map the Dicom Data to the Dicom volume
	    for (k=0;k<Dicom->numSlices;k++)
	    {
	        for(j=0;j<Dicom->height;j++)
	        {
	            for(i=0;i<Dicom->width;i++)
	            {   
	                 // map to a voxel
	                idx = (int)round(((EPIDx*0.5 + Dicom->xCoords[k][j][i]*rescale)/EPIDx)*EGS->xSize);
	                idy = (int)round((Dicom->yCoords[k][j][i]-Dicom->ylims[0])/(Dicom->ylims[1]-Dicom->ylims[0])*(yVolVoxels));
	                idz = (int)round(((EPIDz*0.5 - Dicom->zCoords[k][j][i]*rescale)/EPIDz)*EGS->zSize);;
	                
	                if (idx >=0 && idx < EGS->xSize && idz >= 0 && idz < EGS->zSize)
	                {
	                    // add the CT value to density array to save storage
	                    EGS->voxelDensity[idz][idy][idx] +=(float) Dicom->data_uc[k][j][i];
	                   
	                    //increment the char in the Hits array (might be better to use medium array)
	                    if(hits[idz][idy][idx] < HITS_MAX)
	                    {
	                        hits[idz][idy][idx]++;
	                    }
	                }
	            }
	        }
	    }
	}
	else if (Dicom->imageDataType == 3) // unsigned short
	{
	    // now map the Dicom Data to the Dicom volume
	    for (k=0;k<Dicom->numSlices;k++)
	    {
	        for(j=0;j<Dicom->height;j++)
	        {
	            for(i=0;i<Dicom->width;i++)
	            {
	                // map to a voxel
	                idx = (int)round(((EPIDx*0.5 + Dicom->xCoords[k][j][i]*rescale)/EPIDx)*EGS->xSize);
	                idy = (int)round((Dicom->yCoords[k][j][i]-Dicom->ylims[0])/(Dicom->ylims[1]-Dicom->ylims[0])*(yVolVoxels));
	                idz = (int)round(((EPIDz*0.5 - Dicom->zCoords[k][j][i]*rescale)/EPIDz)*EGS->zSize);;
	                
	                if (idx >=0 && idx < EGS->xSize && idz >= 0 && idz < EGS->zSize)
	                {
	                    // add the CT value to density array to save storage
	                    EGS->voxelDensity[idz][idy][idx] +=(float) Dicom->data_us[k][j][i];
	                   
	                    //increment the char in the Hits array (might be better to use medium array)
	                    if(hits[idz][idy][idx] < HITS_MAX)
	                    {
	                        hits[idz][idy][idx]++;
	                    }
	                }
	            }
	        }
	    }
	}
	
	//use info in the medium and desity arrrays to find actual data
	for (k=0;k<EGS->zSize;k++)
	{
	    for(j=0;j<EGS->ySize-EPIDVoxels;j++)
	    {
	        for(i=0;i<EGS->xSize;i++)
	        {
	            if(hits[k][j][i])
	            {
	                CTavg = EGS->voxelDensity[k][j][i]/hits[k][j][i];
	                //printf("%f\n",CTavg);
	            }
	            else
	            {
	                CTavg = densityTransfers[0][0][0];
	            }
	            
	            if(CTavg > densityTransfers[0][0][0])
	            {
	                //find which category of material this CTavg is within
				    for(w=0;w<numMedia-1;w++)
				    {
					    if(CTavg >= densityTransfers[w][0][0] && CTavg <= densityTransfers[w][numControlPoints[w]-1][0])
					    {
						    break;					    }
				    }
				   
				    // classified
				    EGS->voxelMedium[k][j][i] = mediaCats[w];
				    lowerBound =0;
				    // now find which interval of the transfer function for this material the CTavg in within
				    for(c=0;c<numControlPoints[w]-2;c++)
				    {
					    if( CTavg >= densityTransfers[w][c][0] && CTavg <= densityTransfers[w][c+1][0])
					    {
						    break;
					    }
				    }

				    //interpolate using linear interpolation
				    EGS->voxelDensity[k][j][i] = ((densityTransfers[w][c+1][1]-densityTransfers[w][c][1])*(CTavg - densityTransfers[w][c][0]))/(densityTransfers[w][c+1][0]-densityTransfers[w][c][0]) + densityTransfers[w][c][1];
				}
				else
				{
				    EGS->voxelDensity[k][j][i] = densityTransfers[0][0][1];
				    EGS->voxelMedium[k][j][i] = mediaCats[0];
				}
				
	        }
	    }
	  
	    
	    // now for the EPID
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
		EGS->xBoundaries[i] = -(EPIDx*0.5) + i*vx;
	}
	for (j=0;j<=EGS->ySize;j++)
	{
		EGS->yBoundaries[j] = Dicom->ylims[0]*rescale + j*vy; 
	}

	for (k=0;k<=EGS->zSize;k++)
	{
		EGS->zBoundaries[k] = -(EPIDz*0.5) + k*vz;
	}

	// now remove the cushion
	EGS->RemoveCushion(mediaCats[1],mediaCats[0],mediaCats[1]);
	
}

/* READEPIDSPEC: Reader for the EPID info file used to generate EPID Layers in ConvertDICOMtoEGSPhant()
 *
 * Parameters: 
 *             filename - name of the EPID file
 */
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

/* READINP: Reads the .inp input file that specifies parameters for conversion to EGSPhant format
 * 
 * Parameters: 
 *            filename - name of the .inp file
 */
int ReadINP(const char* filename)
{
	FILE* file;
	int i,j;
	float tmp=0;
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

	//to make the order of the bounds irrelevent, if they aint in the right order, swap 'em 
	if (xmin > xmax)
	{
		tmp = xmin;
		xmin = xmax;
		xmax = tmp;
	}
	if (ymin > ymax)
	{
		tmp = ymin;
		ymin = ymax;
		ymax = tmp;
	}
	if (zmin > zmax)
	{
		tmp = zmin;
		zmin = zmax;
		zmax = tmp;
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

	// read number of control points list
	numControlPoints = new int[numMedia];

	for (i=0;i<numMedia;i++)
	{
		if(!(fscanf(file,"%d ",&(numControlPoints[i]))==1))
		{
			printf("ERROR: could not read control point values\n");
			return 1;
		}
		printf("control points for %d: %d\n",i+1,numControlPoints[i]);
	}
	fscanf(file,"\n",temp);

	// intialise the Material Name,Category, and conversion data arrays
	MediaNames = new std::string[numMedia];

	mediaCats = new char[numMedia];
	estepes = new float[numMedia]
;
	densityTransfers = new float**[numMedia];
	for (i=0; i<numMedia; i++)
	{
		densityTransfers[i] = new float*[numControlPoints[i]];
		for(j=0; j< numControlPoints[i];j++)
		{
			densityTransfers[i][j] = new float[2];
		}
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
		printf("Media name: %s\n",MediaNames[i].c_str());

		//assign a category number
		sprintf(s,"%d",i+1);
		mediaCats[i] = s[0]; 

		for (j=0; j<numControlPoints[i];j++)
		{
			if (!(fscanf(file,"%f,%f,",&(densityTransfers[i][j][0]),&(densityTransfers[i][j][1]))==2))
			{
				printf("ERROR: Could not read transfer function for material %d\n",i+1);
				return 1;
			}
			printf("%f %f\n",densityTransfers[i][j][0],densityTransfers[i][j][1]);
		}
		//read the estepe
		if(!(fscanf(file,"%f\n",&(estepes[i]))==1))
		{
			printf("ERROR: Could not read ESTEPE for material %d\n",i);
			return 1;
		}
		printf("estep: %f\n",estepes[i]);
	}

	return 0;
}

/* PRINTHELP: just dumps an options menu
 */
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

/* READARGS: Reads the commandline arguments passed to CTCombine
 */
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
				if(ReadINP(argv[++i]))// testing new reader
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
		printf("\nDICOM Directory: %s\n",dir);
		DICOMReader DICOM(dir); //initialise the DICOM reader
		printf("Reading DICOM data...\n");
		DICOM.ReadSample();
		DICOM.ReadAll();
		
		printf("Read Complete...\n");
		//return 0;

		// if no isocentre specified then get the mid point
		if(!flag) 
		{
		    printf("\nNo User Suppied Isocentre...");
		    printf("Calculating...\n");
			DICOM.GetMidPoint(&ix,&iy,&iz);
			ix /= scaleFactor;
			iy /= scaleFactor;
			iz /= scaleFactor;
			printf("Calculated Isocentre is (%f %f %f)\n", ix,iy,iz);
		}
		else
		{
		    printf("User Supplied Isocentre is (%f %f %f).\n",ix,iy,iz);
		}
		// Select from the rotated data the desired data region
		printf("\nSetting bounds...\n");
		printf("Oringinal bounds were (%f, %f),(%f, %f),(%f, %f)\n",DICOM.xlims[0],DICOM.xlims[1],DICOM.ylims[0],DICOM.ylims[1],DICOM.zlims[0],DICOM.zlims[1]);
		
		//Case 1: bounds within EPID Dimension
		if(xmax - xmin <= EPIDx && zmax - zmin <= EPIDz)
		{
		    if(DICOM.Select_Region(xmin*scaleFactor,xmax*scaleFactor,ymin*scaleFactor,ymax*scaleFactor,zmin*scaleFactor,zmax*scaleFactor))
		    {
			    printf("Could not complete operation exiting...\n");
			    return 1;
		    }
		    printf("Selected region Was (%f, %f), (%f, %f), (%f, %f)\n",xmin,xmax,ymin,ymax,zmin,zmax);
		}
		else //Case 2: we need to crop the data volume
		{
		    printf("Selected region exceeds EPID dimensions... cropping...\n");
		    
		    while(xmax - xmin > EPIDx)
		    {
		        xmax--;
		        xmin++;
		    }
		    
		    while(zmax - zmin > EPIDz)
		    {
		        zmax--;
		        zmax++;
		    }
		    
		     if(DICOM.Select_Region(xmin*scaleFactor,xmax*scaleFactor,ymin*scaleFactor,ymax*scaleFactor,zmin*scaleFactor,zmax*scaleFactor))
		    {
			    printf("Could not complete operation exiting...\n");
			    return 1;
		    }
		    printf("Selected region Was (%f, %f), (%f, %f), (%f, %f)\n",xmin,xmax,ymin,ymax,zmin,zmax);
		}
			
		//make rotation from zero gantry angle to the new gantry angle 
		printf("\nRotating To new granty angle...\n");
		DICOM.Rotate3D(-Z_gantry[0] + theta,zaxis,ix*scaleFactor,iy*scaleFactor,iz*scaleFactor);
		DICOM.Rotate3D(-Z_gantry[1] + phi,xaxis,ix*scaleFactor,iy*scaleFactor,iz*scaleFactor);
		printf("Rotation Complete...\n");
		
		printf("\nTranslating to isocentre (%f, %f, %f)...\n",ix,iy,iz);
		DICOM.Translate(ix*scaleFactor,iy*scaleFactor,iz*scaleFactor);
		
		printf("\nConverting to .egsphant...\n");
		//  now initailise EGSPhant object to store converted data
		EGSPhant EGSct(outputfileName);
		// convert
		ConvertDICOMToEGSPhant(&DICOM,&EGSct);
		printf("data converted...\n");

		printf("Re-defining coordinates for work with Dosxyz...\n");
		EGSct.ShiftToOrigin(0,-100,0);
		
		printf("\nWriting EGS file...\n");
		// hmm... what does the next line do? Oh yes, thats right... we write to a file
		EGSct.Write();
	}

	return 0;
}

