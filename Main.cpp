#include "./EGSPhant.h"

#include <stdio.h>
#include <string>

using namespace EGS;


void Get_EGSPHANT_Stats(const char *filename)
{
	
	EGSPhant data((char *)filename);
	int result = data.Read();
	result = data.PrintStats();		
	
}


/*
 * Combine the CT with the EPID data
 *
 * Params: ct_filename - the phantom dataset
 * 	   epid_filename - the epid dataset
 *
 */
void Combine_EGSPHANT(const char *ct_filename, const char *epid_filename, const char *output_filename)
{
	int i, j, k;

	
	

	// Read in the phantom data
	EGSPhant ct_data((char *)ct_filename);
	int result1 = ct_data.Read();

	// Read in the EPID data
	EGSPhant epid_data((char *)epid_filename);
	int result2 = epid_data.Read();
	
	// TODO: Need to ensure that these data may be combined - check dimensions
	if (ct_data.xSize != epid_data.xSize )
	{
		fprintf(stdout, "Cannot be combined: Check number of voxels in X dimension.\n");
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
	combine.ySize = ct_data.ySize + epid_data.ySize;
	combine.zSize = epid_data.zSize;

	/////
	//	Populate the Boundaries
	/////

	// X dimension - no changes from ct_datas
	combine.xBoundaries = new float[ ct_data.xSize + 1 ];
	for (i = 0; i<( ct_data.xSize + 1 ); i++)
		combine.xBoundaries[i] = ct_data.xBoundaries[i];
	// Y dimension - EPID followed by CT
	combine.yBoundaries = new float[epid_data.ySize + ct_data.ySize + 1];
	for (i = 0; i<(epid_data.ySize+1); i++)
		combine.yBoundaries[i] = epid_data.yBoundaries[i];
	for (i = 1; i<(ct_data.ySize+1); i++)
		combine.yBoundaries[epid_data.ySize + i] = ct_data.yBoundaries[i];
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
	// Start Initialising
	for (k = 0; k<combine.zSize; k++)
	{
		for (j = 0; j<ct_data.ySize; j++)
			for (i = 0; i<ct_data.xSize; i++)
				combine.voxelMedium[k][j][i] = ct_data.voxelMedium[k][j][i];
		for (j = 0; j<epid_data.ySize; j++)
			for (i = 0; i<epid_data.xSize; i++)
				combine.voxelMedium[k][ct_data.ySize + j][i] = epid_data.voxelMedium[k][j][i];
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
			for (i = 0; i<ct_data.xSize; i++)
				combine.voxelDensity[k][j][i] = ct_data.voxelDensity[k][j][i];
		for (j = 0; j<epid_data.ySize; j++)
			for (i = 0; i<epid_data.xSize; i++)
				combine.voxelDensity[k][ct_data.ySize + j][i] = epid_data.voxelDensity[k][j][i];
	}

	// Write Out the data file
	combine.Write();
}


int main(int argc, char* argv[])
{

	
	// Count the command line arguements
	if (argc < 4 || argc > 4)
	{
		fprintf(stdout, "Wrong number of command line arguements\n");
		fprintf(stdout, "Usage: CTCombine phantom.egsphant epid.egsphant someoutputfile.egsphant\n");
	}
	else
	{
		std::string ct_filename(argv[1]);
		std::string epid_filename(argv[2]);
		std::string output_filename(argv[3]);	
		Combine_EGSPHANT(ct_filename.c_str(), epid_filename.c_str(), output_filename.c_str());
	}
	/*
	std::string output_filename = "output.egsphant";	
	Get_EGSPHANT_Stats(ct_filename.c_str());

	std::string ct_filename = "./Data/Second/cthead.egsphant";
	std::string epid_filename = "./Data/Second/EPID_for_cthead.egsphant";	
	Combine_EGSPHANT(ct_filename.c_str(), epid_filename.c_str(), output_filename.c_str());
	*/

	return 0;
}

