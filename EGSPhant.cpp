/**
 * EGSPhant.cpp (Inplementation)
 *
 * @Author Mark Dwyer
 * @Contact m2.dwyer@qut.edu.au
 * @Created 07/10/2008
 *
 * Most of the crap here has been ripped from:
 * http://www.irs.inms.nrc.ca/BEAM/user_manuals/pirs794/node99.html
 *
 */


#include "EGSPhant.h"
namespace EGS
{
	/****************************
	 *		Constructors
	 ***************************/

	/**
 	 * Default Constructor.
  	 */
	EGSPhant::EGSPhant(void)
	{
		// Nothing to see here - move on.
	};
	
	/**
 	 * Main Constructor.
  	 */
	EGSPhant::EGSPhant(char *filename)
	{
		input_filename = filename;
	};

	/*
	 * Prepare the datastructures and populate from input_filename.
	 *
	 * There could be a future error in the read of the "medium number"
	 * It is hardcoded to read in 1024 char values and will be an error
  	 * if you try and read more in the future - ie. micro_ct
	 * To fix, just change the char from 1024 to the desired amount
	 * See "TODO" for appropriate line of code (LOC)
	 */
	int EGSPhant::Read(void)
	{
		FILE *file;
		int success;
		int i, j, k;		
		if((file=fopen(input_filename,"r"))==NULL)
			fprintf(stdout, "EGSPhant::PrepareRead - Error opening file\n");

		// Read in the number of media
		success = fscanf(file, " %d\n", &numberOfMedia);
		
		// Initialise the names of media and read
		mediaNames = new std::string[numberOfMedia];
		char name[255];
		for (i = 0; i<numberOfMedia; i++)
		{
			success = fscanf(file, "%s", name);
			mediaNames[i] = name;
		}

		// Initialise and read in the ESTEPE
		estepe = new float[numberOfMedia];
		for (i = 0; i<numberOfMedia; i++)
		{
			success = fscanf(file, "%f", &estepe[i]);
		}

		// Read in the number of voxels in each dimension
		success = fscanf(file, "%d", &xSize);
		success = fscanf(file, "%d", &ySize);
		success = fscanf(file, "%d", &zSize);

		// Initialise and populate xBoundaries
		xBoundaries = new float[xSize+1];
		for (i = 0; i<(xSize+1); i++)
			success = fscanf(file, "%f", &xBoundaries[i]);
		// Initialise and populate yBoundaries
		yBoundaries = new float[ySize+1];
		for (i = 0; i<(ySize+1); i++)
			success = fscanf(file, "%f", &yBoundaries[i]);

		// Initialise and populate zBoundaries
		zBoundaries = new float[zSize+1];
		for (i = 0; i<(zSize+1); i++)
			success = fscanf(file, "%f", &zBoundaries[i]);

		// Initialise and start reading in the voxel medium numbers
		// TODO
		char voxel[1024];
		voxelMedium = new char**[zSize]; 
		for (k = 0; k<zSize; k++)
		{
			voxelMedium[k] = new char*[ySize];
			for (j = 0; j<ySize; j++)
			{
				voxelMedium[k][j] = new char[xSize];
			}
		}
		// Start Reading
		for (k = 0; k<zSize; k++)
		{
			for (j = 0; j<ySize; j++)
			{
				success = fscanf(file, "%s", voxel);		
				for (i = 0; i<xSize; i++)
				{	
					voxelMedium[k][j][i] = voxel[i];
				}
			}
		}

		// Initialise and start reading in the densities for each voxel
		float value;
		voxelDensity = new float**[zSize]; 
		for (k = 0; k<zSize; k++)
		{
			voxelDensity[k] = new float*[ySize];
			for (j = 0; j<ySize; j++)
				voxelDensity[k][j] = new float[xSize];
		}
		// Start Reading
		for (k = 0; k<zSize; k++)
			for (j = 0; j<ySize; j++)		
				for (i = 0; i<xSize; i++)		
					success = fscanf(file, "%f", &voxelDensity[k][j][i]);

		fclose(file);
		fprintf(stdout, "%s successfully read.\n", input_filename);
	};


	/*
	 * Prepare the datastructures and populate to input_filename.
	 *
	 */
	int EGSPhant::Write(void)
	{
		FILE *file;
		int success;
		int i, j, k;		
		if((file=fopen(input_filename,"w"))==NULL)
			fprintf(stdout, "EGSPhant::PrepareRead - Error opening file\n");

		// Write the number of media
		fprintf(file, " %d\n", numberOfMedia);
		
		// Initialise the names of media and read
		for (i = 0; i<numberOfMedia; i++)
		{
			fprintf(file, "%s\n", mediaNames[i].c_str());
		}

		// Write the ESTEPE
		for (i = 0; i<numberOfMedia; i++)
		{
			fprintf(file, "%f ", estepe[i]);
		}
		fprintf(file, "\n");
		// Write the number of voxels in each dimension
		//fprintf(file, "%d ", xSize);
		WriteInteger(file, xSize, 5);
		//fprintf(file, "%d ", ySize);
		WriteInteger(file, ySize, 5);
		//fprintf(file, "%d ", zSize);
		WriteInteger(file, zSize, 5);
		fprintf(file, "\n");
		// Write xBoundaries
		for (i = 0; i<(xSize+1); i++)
		{
			fprintf(file, "%f ", xBoundaries[i]);
			if (i%5 == 0 && i!=0)
				fprintf(file, "\n");
		}
		fprintf(file, "\n");
		for (i = 0; i<(ySize+1); i++)
		{
			fprintf(file, "%f ", yBoundaries[i]);
			if (i%5 == 0 && i!=0)
				fprintf(file, "\n");
		}
		fprintf(file, "\n");
		for (i = 0; i<(zSize+1); i++)
		{
			fprintf(file, "%f ", zBoundaries[i]);
			if (i%5 == 0 && i!=0)
				fprintf(file, "\n");
		}
		fprintf(file, "\n");
		// Write the voxel medium numbers
		for (k = 0; k<zSize; k++)
		{
			for (j = 0; j<ySize; j++)
			{	
				for (i = 0; i<xSize; i++)
				{	
					fprintf(file, "%c",voxelMedium[k][j][i]);
					//fprintf(file," ");
				}
				fprintf(file, "\n");
			}
			fprintf(file, "\n");
		}

		// Write the densities for each voxel
		for (k = 0; k<zSize; k++)
		{
			for (j = 0; j<ySize; j++)		
			{
				for (i = 0; i<xSize; i++)
				{		
					fprintf(file, "%f ", voxelDensity[k][j][i]);
					if (i%5 == 0 && i!=0)
						fprintf(file, "\n");
				}
				fprintf(file, "\n");
			}
			fprintf(file, "\n");
		}
		fclose(file);
	};

	/*
	 *  Write an Integer to file that has a character length
	 *
	 *
	 */
	int EGSPhant::WriteInteger(FILE *file, int num, int length)
	{
		char numAsChar[length+1];
		if (num < 10)
			sprintf(numAsChar, "    %d", num);
		else if (num < 100)
			sprintf(numAsChar, "   %d", num);
		else if (num < 1000)
			sprintf(numAsChar, "  %d", num);
		else if (num < 10000)
			sprintf(numAsChar, " %d", num);

		fprintf(file, "%s",numAsChar);
	};

	/*
	* shifts the grid from by -(x,y,z)
	*/
	int EGSPhant::ShiftToOrigin(float x,float y,float z){

		int i;

		//point is to be our new origin
		for(i=0;i<=xSize;i++)
			xBoundaries[i] -= x;
		for(i=0;i<=ySize;i++)
			yBoundaries[i] -= y;
		for(i=0;i<=zSize;i++)
			zBoundaries[i] -= z;
		
	};

	/* Removes the cushion using a very simple (and possibly primitive) algorithm
	 * basically recognises that the cushion is simply a contiguous body of non-air voxels
	 * out side the skin, only a prototype, wont be accurate for very small voxel sizes (were the skin
	 * could be more than 1 voxel thick),  
	 * 
	 */
	int EGSPhant::RemoveCushion(char cushMedium,char airMedium,char skinMedium)
	{
		int i,j,k; // ah yes, our friendly loop pointers again
		
		char boundary[8];
		short t;
		bool result=true;
		// for each voxel in the volume
		for (k=1;k<zSize-1;k++)
		{
			for (j=1;j<ySize-1;j++)
			{
				for (i=1; i<xSize-1; i++)
				{
					result = true;
					
					if(voxelMedium[k][j][i] == cushMedium)
					{
						// get the surrounding voxels (the this slice)
						boundary[0] = voxelMedium[k][j-1][i-1];
						boundary[1] = voxelMedium[k][j-1][i];
						boundary[2] = voxelMedium[k][j-1][i+1];
						boundary[3] = voxelMedium[k][j][i-1];
						boundary[4] = voxelMedium[k][j][i+1];
						boundary[5] = voxelMedium[k][j+1][i-1];
						boundary[6] = voxelMedium[k][j+1][i];
						boundary[7] = voxelMedium[k][j+1][i+1];
					
						t=0;
						// check if they are all air of cushion
						while (result && t<8)
						{
							if(boundary[t] == airMedium || boundary[t] == cushMedium || boundary[t] == skinMedium)
							{
								result = true;
							}
							else
							{
								result = false;
							}

							t++;
						}
						// if this point is surrounded by cushion and air then we clear it 
						if (result)
						{
							voxelMedium[k][j][i] = airMedium;
							// is always air
							voxelDensity[k][j][i] = voxelDensity[0][0][0];
						}
						//printf("medium: %c\n",voxelMedium[k][j][i]);
					}
				}
			}
		}
		//printf("%d %d %d",i,j,k);

	};


	int EGSPhant::PrintStats(void)
	{
		fprintf(stdout,"\n------------------\n");
		fprintf(stdout,"ESGPHANT Information:\n");
		fprintf(stdout,"------------------\n");
		fprintf(stdout,"File Name:\t\t %s\n",input_filename);
		fprintf(stdout,"Number of Slices:\t %d\n",zSize);
		fprintf(stdout,"Slice Resolution:\t %d x %d\n",xSize,ySize);
		
		fprintf(stdout, "Number of Media:\t %d\n", numberOfMedia);
		fprintf(stdout, "Media Names:\n");
		for (int i = 0; i<numberOfMedia; i++)
		{
			fprintf(stdout, "\t\t\t%d. %s\n",i+1, mediaNames[i].c_str());
		}
	};
}
