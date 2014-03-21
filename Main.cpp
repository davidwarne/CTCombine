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


// argument storage

char dir[255]; //DICOM data directory
char INPfile[255];
float xmin,xmax,ymin,ymax,zmin,zmax; //require region boundaries of CT data
float vx,vy,vz=0.5; //voxel size 
float ix,iy,iz=0; // isocentre
float Z_gantry[2]= {270,90}; // direction of beam at zero gantry angle
float theta = 270; // rotation angles
float phi = 90;
float EPID_dist=0; // distance of isocentre to EPID data

char EPIDfileName[255]; //filename of EPID data
char outputfileName[255]; //EGSPHANT output filename

float EPIDx,EPIDz; //EPID surface dimension
int totalEPIDLayers; // number of EPID dimensions
int EPIDnumMedia; 
float** EPIDLayerInfo; // EPIDLayerInfo[0] = material number
                       // EPIDLayerInfo[1] = material density
                       // EPIDLayerInfo[2] = Layer Thickness
std::string* EPIDMaterials; // EPID material names

//similar info for the CT scan
int numMedia;
std::string* MediaNames;
char* mediaCats;
float** conversionData; //stores the conversion data



/*
 * Combine the CT with the EPID data, is probably redundant now we create EPID data on conversion From Dicom
 *
 * Params: ct_filename - the phantom dataset
 * 	   epid_filename - the epid dataset
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
	printf("%f \n",EPID_dist-fabs(ct_data.yBoundaries[ct_data.ySize]-iy));
	airPadding = (int)floor(EPID_dist-fabs(ct_data.yBoundaries[ct_data.ySize]-iy));

	// Need to ensure that these data may be combined - check dimensions
	if (ct_data.xSize > epid_data.xSize )
	{
		exit(1);
	}

	if (ct_data.zSize != epid_data.zSize )
	{
		fprintf(stdout, "Cannot be combined: Check number of voxels in Z dimension.\n");
		exit(1);
	}

	// Initialise the combinatory dataset
	EGSPhant combine((char *)output_filename);
	
	/*
	 *	Create the new EGSPhant file
	 */
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

	/////
	//	Populate the Boundaries
	/////

	// X dimension - no changes from epid_datas
	combine.xBoundaries = new float[epid_data.xSize + 1 ];
	for (i = 0; i<( epid_data.xSize + 1 ); i++)
		combine.xBoundaries[i] = epid_data.xBoundaries[i];
	// Y dimension - EPID followed by CT
	combine.yBoundaries = new float[epid_data.ySize +airPadding+ ct_data.ySize+1];
	for (i = 0; i<(ct_data.ySize+1); i++)
		combine.yBoundaries[i] = ct_data.yBoundaries[i];
	for (i=1;i<=airPadding;i++)
		combine.yBoundaries[ct_data.ySize + i] = ct_data.yBoundaries[ct_data.ySize] + i*vy;
	for (i = 1; i<(epid_data.ySize+1); i++)
		combine.yBoundaries[ct_data.ySize + airPadding + i] = ct_data.yBoundaries[ct_data.ySize + airPadding]+i*vy;
	// Z dimension - No Change
	combine.zBoundaries = new float[ct_data.zSize+1];
	for (i = 0; i<(ct_data.zSize+1); i++)
		combine.zBoundaries[i] = ct_data.zBoundaries[i];
	
	/////
	//	Create the voxel medium numbers 
	/////
	combine.voxelMedium = new char**[combine.zSize]; 
	for (k = 0; k<combine.zSize; k++)
	{
		combine.voxelMedium[k] = new char*[combine.ySize];
		for (j = 0; j<combine.ySize; j++)
			combine.voxelMedium[k][j] = new char[combine.xSize];
	}

	// the column in which the CT data starts 
	int dataStart = (int)floor((epid_data.xSize - ct_data.xSize)/2);
	
	// Start Initialising
	for (k = 0; k<combine.zSize; k++)
	{
		for (j = 0; j<ct_data.ySize; j++){
			for (i = 0; i<epid_data.xSize; i++){
				//only add the data within its region
				if (i>=dataStart && i<(ct_data.xSize+dataStart)){
				combine.voxelMedium[k][j][i] = ct_data.voxelMedium[k][j][i-dataStart];
				}
				else{
					//pad with air
					combine.voxelMedium[k][j][i] = '1';
				}
			}
		}

		//pad the gap between the EPID and CT with air
		for (j=0;j<airPadding;j++)
			for (i=0;i<epid_data.xSize;i++)
				combine.voxelMedium[k][ct_data.ySize +j][i] = '1';
	
		//write the EPID
		for (j = 0; j<epid_data.ySize; j++)
			for (i = 0; i<epid_data.xSize; i++)
				combine.voxelMedium[k][ct_data.ySize + airPadding + j][i] = epid_data.voxelMedium[k][j][i];
	}

	/////
	//	Create the voxel densities 
	/////
	combine.voxelDensity = new float**[combine.zSize]; 
	for (k = 0; k<combine.zSize; k++)
	{
		combine.voxelDensity[k] = new float*[combine.ySize];
		for (j = 0; j<combine.ySize; j++)
			combine.voxelDensity[k][j] = new float[combine.xSize];
	}
	// Start Initialising
	for (k = 0; k<combine.zSize; k++)
	{
		for (j = 0; j<ct_data.ySize; j++)
			for (i = 0; i<ct_data.xSize; i++){
				if (i>=dataStart && i<(ct_data.xSize+dataStart)){
				combine.voxelDensity[k][j][i] = ct_data.voxelDensity[k][j][i-dataStart];
				}
				else{
					//pad with air
					combine.voxelDensity[k][j][i] = 0.001;
				}
			}
				
		for (j=0;j<airPadding;j++)
			for (i=0;i<epid_data.xSize;i++)
				combine.voxelDensity[k][ct_data.ySize +j][i] = 0.001;
	
		for (j = 0; j<epid_data.ySize; j++)
			for (i = 0; i<epid_data.xSize; i++)
				combine.voxelDensity[k][ct_data.ySize + airPadding + j][i] = epid_data.voxelDensity[k][j][i];
	}

	// Write Out the data file
	
	combine.Write();
}


/*
* Converts DICOM data as extracted from a DICOMreader object and converts it to the voxelised
* EGSPhant format
*
*/
void ConvertDICOMToEGSPhant(DICOMReader* Dicom,EGSPhant* EGS){

	//DICOM data coords are in mm and EGSPhant is in cm
	float rescale = 0.1;
	float boundaryData[8]; //temp for voxel value calculation
	FILE* ConvertData;
	int i,j,k;
	int idx,idy,idz;
	float lr,ls,lt,r,s,t,dx,dy,dz,total,CTavg,EPIDdepth;
	unsigned w;
	float lowerBound;
	int nonZero,dataStart,airPadding;
	int EPIDVoxels=0;
	char temp[255];
	int a,b;
	int startLayer;
	int EPIDtop;

	for(i=0;i<totalEPIDLayers;i++){
		EPIDVoxels += (int)ceil(EPIDLayerInfo[i][2]/vy);
	}

	
	airPadding = (int)floor((EPID_dist - fabs(Dicom->yCoords[Dicom->height-1]*rescale - iy))/vy);
	dataStart = (int)floor(((EPIDx - fabs(Dicom->xCoords[0] - Dicom->xCoords[Dicom->width -1])*rescale)/vx)*0.5 );
	//printf("%d",dataStart);
	//calculate the number of voxels
	//EGS->xSize = (int)floor((fabs(Dicom->xCoords[0] - Dicom->xCoords[Dicom->width -1])*rescale)/vx);
	EGS->xSize = (int)floor(EPIDx/vx);
	EGS->ySize = (int)floor(((fabs(Dicom->yCoords[0] - Dicom->yCoords[Dicom->height -1])*rescale))/vy)+airPadding+EPIDVoxels;
	EGS->zSize = (int)floor((fabs(Dicom->zCoords[0] - Dicom->zCoords[Dicom->numSlices -1])*rescale)/vz);

	printf("Original Volume Dimensions: %d %d %d\n",Dicom->width,Dicom->height,Dicom->numSlices);
	printf("Voxel Dimesions %f %f %f\n",vx,vy,vz);
	printf("Voxelised Grid Dimensions with EPID: %d %d %d\n",EGS->xSize,EGS->ySize,EGS->zSize);
	

	//TODO: continue to include the Creation and addition of EPID data
	dx = (vx*(Dicom->width-1))/(fabs(Dicom->xCoords[0] - Dicom->xCoords[Dicom->width -1])*rescale);
	dy = (vy*(Dicom->height-1))/(fabs(Dicom->yCoords[0] - Dicom->yCoords[Dicom->height -1])*rescale);
	dz = (vy*(Dicom->numSlices-1))/(fabs(Dicom->zCoords[0] - Dicom->zCoords[Dicom->numSlices -1])*rescale);
	//printf("%f %f %f\n",dx,dy,dz);

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
	for (k=0;k<EGS->zSize;k++){
		EGS->voxelMedium[k] = new char*[EGS->ySize];
		EGS->voxelDensity[k] = new float*[EGS->ySize];
		for (j=0;j<EGS->ySize;j++){
			//printf("k,j = %d %d\n",k,j);
			EGS->voxelMedium[k][j] = new char[EGS->xSize];
			EGS->voxelDensity[k][j] = new float[EGS->xSize];
		}
	}

	//store the ESTEP
	for (i=0;i<numMedia;i++){
		EGS->estepe[i] = conversionData[i][3];
	}
	for (i=numMedia;i<EPIDnumMedia;i++){
		EGS->estepe[i] = 1;
	}

	


	//printf("do we get to here?");
	//for each voxel
	for(k=0;k<EGS->zSize;k++){

		//the section where the CT data will be written
		for(j=0;j<EGS->ySize-airPadding-EPIDVoxels;j++){

			//pad with air
			for(i=0;i<dataStart;i++){
				EGS->voxelMedium[k][j][i] = '1';
				EGS->voxelDensity[k][j][i] = AIRDENSE;
			}
			//interpolate the actual data
			for(i=dataStart;i<EGS->xSize-dataStart;i++){

				//printf("%d\n",i);
				total =0;
				nonZero = 0;
				//interpolate the data on the boundaries 
				for (w=0;w<8;w++){
							
					lr = (float)((i-dataStart+(w & 1)))*dx;
					ls = (float)((j+((w & 2)>>1)))*dy;
					lt = (float)((k+((w & 4)>>2)))*dz;
					idx = (int)floor(lr);
					idy = (int)floor(ls);
					idz = (int)floor(lt);
					r = lr - idx;
					s = ls - idy;
					t = lt - idz;
					//printf("%f %f %f\n",ls,lr,lt);
					//printf("%d %d %d, %d\n",idx ,idy,idz, Dicom->data_s[idz][idy][idx]);
					//printf("%f %f %f\n",r,s,t);
		
					if (idx >= 0 && idy >= 0 && idz >= 0 && (idx < Dicom->width-1) && (idy < Dicom->height-1) && (idz < Dicom->numSlices-1)){
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

						//if(Dicom->data_s[k][j][i] >=0)printf("%d\n",Dicom->data_s[idz][idy][idx]);
						//printf("%f, (%d %d %d) (%d %d %d) (%f %f %f)\n",total,i,j,k,idx,idy,idz,r,s,t);
						
					}
					else{//just set it to zero
						total +=0;
						//printf("%f, (%d %d %d) (%d %d %d) (%f %f %f)\n",total,i,j,k,idx,idy,idz,r,s,t);
					}// we have to deal with the upper boundary cases better

					if(total)
						nonZero++;
				}
				
				//for now just assign the voxel the CT number of the mid point as interpolated by the boundaries
				if (nonZero)
					CTavg = total/nonZero;
				else
					CTavg =0;

				

				//if(CTavg >=0){printf("%f\n",CTavg);}
				
				//categorise Voxel
				lowerBound =0;
				//printf("do we get here? %f, %f\n",conversionData[0][0],CTavg);
				for(w=0;w<numMedia;w++){
					if( !(CTavg >= lowerBound && CTavg < conversionData[w][0]) && !(CTavg < 0) && !(CTavg > conversionData[numMedia-1][0])){
						lowerBound = conversionData[w][0];
					}
					else{
						break;
					}
				}
				
				//write classification of voxel
				EGS->voxelMedium[k][j][i] = mediaCats[w];
				//if (!(EGS->voxelMedium[k][j][i] ^ '4')){printf("%c\n",EGS->voxelMedium[k][j][i]);}
				
				if (CTavg < 0){//if we are less than zero than set to lower bound of AIR
					EGS->voxelDensity[k][j][i] = conversionData[0][1];
				}
				else{//otherwise use linear interpolation for the density based on conversion data
					EGS->voxelDensity[k][j][i] = ((conversionData[w][2]-conversionData[w][1])*(CTavg - lowerBound))/(conversionData[w][0]-lowerBound) + conversionData[w][1];
					//if (!(EGS->voxelMedium[k][j][i] ^ '3')){printf("%f\n",EGS->voxelDensity[k][j][i]);}
				}
			}
			//pad with air 
			for(i=EGS->xSize-dataStart;i<EGS->xSize;i++){
				EGS->voxelMedium[k][j][i] = '1';
				EGS->voxelDensity[k][j][i] = AIRDENSE;
			}
		}

		for(j=EGS->ySize-airPadding-EPIDVoxels;j<EGS->ySize-EPIDVoxels;j++){
			//theses layers are all air
			for(i=0;i<EGS->xSize;i++){
				EGS->voxelMedium[k][j][i] = '1';
				EGS->voxelDensity[k][j][i] = AIRDENSE;
			}
		}

		//now for the EPID

	
		EPIDtop = EGS->ySize-EPIDVoxels;
		a=0;
		b=0;
		startLayer =0;
		//get the start of each layer
		for(w=0;w<totalEPIDLayers;w++){

			a = startLayer;
			b = startLayer+(int)ceil(EPIDLayerInfo[w][2]/vy);
			for(j=EPIDtop+a;j<EPIDtop+b;j++){
				for(i=0;i<EGS->xSize;i++){
					sprintf(temp,"%d",(int)EPIDLayerInfo[w][0]);
					EGS->voxelMedium[k][j][i] = temp[0];
					EGS->voxelDensity[k][j][i] = EPIDLayerInfo[w][1];
				}
			}
			startLayer = b; 
			
		}
		
	}
	//printf("do we get here?");

	for (i=0;i<=EGS->xSize;i++){
				EGS->xBoundaries[i] = (EPIDx*0.5) - i*vx;
	}
	for (j=0;j<=EGS->ySize;j++){
			EGS->yBoundaries[j] = -Dicom->yCoords[0]*rescale - j*vy;
	}

	for (k=0;k<=EGS->zSize;k++){
		EGS->zBoundaries[k] = Dicom->zCoords[0]*rescale + k*vz;
	}

	//calculate voxel coordinate boundaries
	

}


void ReadDicoms_Test( const char *directory,float deg, int slice)
{
	//printf("do we enter?");
	DICOMReader reader((char *)directory);
	
	reader.ReadSample();
	float axis[3] = {0,0,1};
	reader.ReadAll();
	reader.WriteTestFile();
	//reader.Select_Region(-10,10,-100,-90,-10,10);
	//reader.PrintAllFiles();
	//for (int i = 0; i<reader.numSlices; i++)
	//reader.OutputSliceToFile(slice,"Pre-rotation");
	reader.Rotate3D(deg,axis,ix,iy,iz);
	//printf("rotated and angle of %f\n",deg);
	//reader.OutputSliceToFile(slice,"Post-rotation");
	

}

void TestEGSPhantRotate(const char* file, float deg,const char* outfile){
	
	float axis[3] = {0,0,1};
	EGSPhant data((char*)file);
	data.Read();
	//data.RotateAboutZaxis(deg);
	data.Rotate3D(deg,axis);
	//data.Write();
	

}

int ReadEPIDspec(const char* filename){

	FILE* file;
	int i;
	

	char tempName[255];
	if(!(file = fopen(filename,"r"))){
		printf("ERROR: could not open EPID data file\n");
		return 1;
	}

	fscanf(file,"%f %f\n",&EPIDx,&EPIDz);
	fscanf(file,"%d\n",&EPIDnumMedia);
	EPIDMaterials = new std::string[EPIDnumMedia];
	for (i=0;i<EPIDnumMedia;i++){
		fscanf(file,"%s",tempName);
		EPIDMaterials[i] = tempName;
		printf("%s\n",EPIDMaterials[i].c_str());
	}


	fscanf(file,"%d\n",&totalEPIDLayers);

	//initailise EPID layer data
	EPIDLayerInfo = new float*[totalEPIDLayers];
	for(i=0;i<totalEPIDLayers;i++)
		EPIDLayerInfo[i] = new float[3];

	for(i=0;i<totalEPIDLayers;i++){
		fscanf(file,"%f %f %f",&(EPIDLayerInfo[i][0]),&(EPIDLayerInfo[i][1]),&(EPIDLayerInfo[i][2]));
		//printf("%f %f %f\n",EPIDLayerInfo[i][0],EPIDLayerInfo[i][1],EPIDLayerInfo[i][2]);
	}


}

int ReadINP(const char* filename){

	FILE* file;
	int i;
	
	if(!(file = fopen(filename,"r"))){
		printf("ERROR: could not open CT to Density convert data file\n");
		return 1;
	}
	char temp[255];
	char s[255];
	
	//skip the first two lines
	fscanf(file,"%s\n",temp);
	fscanf(file,"%s\n",temp);
	if(!(fscanf(file,"%f, %f, %f, %f, %f, %f\n",&xmin,&xmax,&ymin,&ymax,&zmin,&zmax)==6))
		return 1;
	if(!(fscanf(file,"%f, %f, %f\n",&vx,&vy,&vz)==3))
		return 1;
	if(!(fscanf(file,"%d\n",&numMedia)==1))
		return 1;
	

	MediaNames = new std::string[numMedia];
	mediaCats = new char[numMedia];

	conversionData = new float*[numMedia];

	for (i=0;i<numMedia;i++)
		conversionData[i] = new float[4];

	for (i=0;i<numMedia;i++){
		if (!(fscanf(file,"%s",temp)==1))
			return 1;
		MediaNames[i] = temp;
		sprintf(s,"%d",i+1);
		mediaCats[i] = s[0];
		if (!(fscanf(file,"%f,%f,%f,%f\n",&(conversionData[i][0]),&(conversionData[i][1]),&(conversionData[i][2]),&(conversionData[i][3]))==4))
			return 1;
	}

	return 0;
}

int printHelp(){
	printf("CTCombine: version 0.12\n");
	printf("\n");
	printf("The Following switches are used to incicate the data being input.\n");
	printf("all of the following parameters are required:\n");
	printf("\n");
	printf("-d dir                   the following parameter is the directory if the DICOM data files\n");
	printf("-E filename              the following parameter is the name of the EPID spec file\n");
	printf("-o outfile               the following parameter is the name of the output.egsphant file\n");
	printf("-inp filename            the following is the file name of the .inp input file\n");
	printf("\n");
	printf("The following are optional:");
	printf("\n");
	printf("-i x y z                 indicates the desired iso-centre\n");
	printf("-e dist                  indicates the desired distance from the isocentre to the top of the EPID\n");
	printf("-g theta phi             indicates the theta and phi angle of zero gantry angle\n");
	printf("-r theta phi             indicates the theta and phi angle of desired gantry angle\n");
	printf("-c file1 file2 file3     Combines the following two files (file1,file2) into the one output file indicayed by file3");
}

int ReadArgs(int argc, char* argv[]){

	int i;

	if(argc>1){
	//need to put check in to make sure the correct numbe of arguments follow each switch
	for(i=1;i<argc;i++){

		if(!strcmp(argv[i],"-d")){//indicates INP input data file

			sprintf(dir,"%s",argv[++i]);
			//printf("read Directory %s\n",dir);
		}
		else if (!strcmp(argv[i],"-inp")){
			ReadINP(argv[++i]);
		}
		else if (!strcmp(argv[i],"-i")){//indicates isocentre

			ix = atof(argv[++i]);
			iy = atof(argv[++i]);
			iz = atof(argv[++i]);
			//printf("read isocentre: (%f %f %f)\n",ix,iy,iz);
		}
		else if (!strcmp(argv[i],"-g")){//indicates direction of zero gantry angle
			Z_gantry[0] = atof(argv[++i]);//theta angle
			Z_gantry[1] = atof(argv[++i]);//phi angle
			//printf("read zero gantry angle: theta = %f, phi = %f\n",Z_gantry[0],Z_gantry[1]);
		}
		else if (!strcmp(argv[i],"-r")){//indicates theta and phi rotation angles
			theta = atof(argv[++i]);
			phi = atof(argv[++i]);
			//printf("read desired angle: theta: %f, phi %f\n",theta,phi);
		}
		else if (!strcmp(argv[i],"-e")){//indicates distances desired from isocentre to EPID data
			EPID_dist = atof(argv[++i]);

			//printf("read disired distance from EPID surface to isocentre %f cm",EPID_dist);
		}
		else if (!strcmp(argv[i],"-c")){//indicates CT number to density convesion data
			
			Combine_EGSPHANT(argv[i+1],argv[i+2],argv[i+3]);
			return 0;
		}
		else if (!strcmp(argv[i],"-E")){//indicates EPID data file
			
			ReadEPIDspec(argv[++i]);
			
			//printf("read EPID data file %s\n",EPIDfileName);
		}
		else if (!strcmp(argv[i],"-o")){//indicates outputfile name
			
			sprintf(outputfileName,"%s",argv[++i]);
			//printf("read output file name: %s\n",outputfileName);
		}
		else{
			printHelp();
			return 1;
		
		}
	}
	}
	else{
		printHelp();
		return 1;
	}
	return 0;
}




int main(int argc, char* argv[])
{
	int i;
	float xaxis[3] = {1,0,0};
	
	float zaxis[3] = {0,0,1};

	if(ReadArgs(argc,argv)){
		printf("Program cannot complete due to invalid input\n");
	}
	else
	{
		printf("dir: %s",dir);
		DICOMReader DICOM(dir);
		printf("Reading DICOM data...\n");
		DICOM.ReadSample();
		DICOM.ReadAll();
		printf("Read Complete...\n");

		printf("Rotating for new gantry angle...\n");
		DICOM.Rotate3D(-Z_gantry[0] + theta,zaxis,ix*10,iy*10,iz*10);
		DICOM.Rotate3D(-Z_gantry[1] + phi,xaxis,ix*10,iy*10,iz*10);
		DICOM.Select_Region(xmin*10,xmax*10,ymin*10,ymax*10,zmin*10,zmax*10);
		DICOM.Translate(ix*10,iy*10,iz*10);
		printf("Rotation Complete...\n");
		printf("Converting to .egsphant...\n");
		EGSPhant EGSct(outputfileName);
		ConvertDICOMToEGSPhant(&DICOM,&EGSct);
		printf("data converted...\n");
		printf("Writing EGS file...\n");
		EGSct.Write();
		//printf("Combining EGS files...\n");
		//Combine_EGSPHANT(EGSct.input_filename,EPIDfileName,outputfileName);
		//printf("Data Combined...\n");



	}

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



