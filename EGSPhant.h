/**
 * EGSPhant.h (Definition)
 *
 * @Author Mark Dwyer
 * @Contact m2.dwyer@qut.edu.au
 * @Created 07/10/2008
 *
 * Most of the crap here has been ripped from:
 * http://www.irs.inms.nrc.ca/BEAM/user_manuals/pirs794/node99.html
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#ifndef EGSPHANT_H
#define EGSPHANT_H

template <class T> const T& min ( const T& a, const T& b ) {
  return (a<b)?a:b;    
}

template <class T> const T& max ( const T& a, const T& b ) {
  return (a>b)?a:b;     
}

const float PI = 3.14159265358979323846;



namespace EGS
{
	/**
	 * EGSPhant
	 *
	 * Pretty standard implementation
	 */
	
	class EGSPhant
	{
		public:
			// Constructor
			EGSPhant(void);
			EGSPhant(char *filename);

			int Read(void);
			int Write(void);
			int ShiftToOrigin(float x,float y,float z);
			int WriteInteger(FILE *file, int num, int length);
			int RemoveCushion(char cushMedium,char airMedium,char skinMedium);
			int PrintStats(void);

		public:
			// Variables
			char *input_filename;
			
			// The number of media in the phantom
			int numberOfMedia;
			// The names of the media
			std::string *mediaNames;
			// The ESTEPE for each medium (now a dummy point)
			float *estepe;
			// The number of voxels in the X, Y and Z directions
			int xSize, ySize, zSize;
			// A list of the voxel boundaries in the X direction
			float *xBoundaries;
			// A list of the voxel boundaries in the Y direction
			float *yBoundaries;
			// A list of the voxel boundaries in the Z direction
			float *zBoundaries;
			// For each Z slice, an XY matrix containing the medium number in each voxel
			char ***voxelMedium;
			// For each slice in the Z direction, an XY matrix containing the densities in each voxel
			float ***voxelDensity;
			
	};

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
		fprintf(stdout, "%s successfully written!\n", input_filename);
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
		fprintf(stdout, "Number of Media = %d\n", numberOfMedia);
		for (int i = 0; i<numberOfMedia; i++)
		{
			fprintf(stdout, "%s\n", mediaNames[i].c_str());
		}
		
		for (int i = 0; i<numberOfMedia; i++)
		{
			fprintf(stdout, "%f ", estepe[i]);
		}
		fprintf(stdout, "\n");
		fprintf(stdout, "%d %d %d\n", xSize, ySize, zSize);

	};

} //namespace EGS

#endif //MYMATH_MATRIX_H
