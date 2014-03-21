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
			int WriteRotatedData(const char* outfilename,float*** voxelMediumRotated,float*** voxelDensityRotated); // just an overload for quicker output
			int ShiftToOrigin(float x,float y,float z,int flag);
			void FlipY();
			int WriteInteger(FILE *file, int num, int length);
			int RotateAboutZaxis(float theta);
			int RemoveCushion(char cushMedium,char airMedium,char skinMedium);
			int Rotate3D(float theta, float* axis);
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

	/* Overload: write to specified file
	 * Prepare the datastructures and populate to input_filename.
	 *
	 */
	int EGSPhant::WriteRotatedData(const char* outfilename,float*** voxelMediumRotated,float*** voxelDensityRotated)
	{
		FILE *file;
		int success;
		int i, j, k;		
		if((file=fopen((char*)outfilename,"w"))==NULL)
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
					fprintf(file, "%f ", voxelMediumRotated[k][j][i]);
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
					fprintf(file, "%f", voxelDensityRotated[k][j][i]);
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
	 * rotates the volume by angle theta in the xy (yet to be optimised)
	 */
	int EGSPhant::RotateAboutZaxis(float theta)
	{
		int i,j,k; // the usual loop counters
		float rads;
		float sintheta;
		float costheta;
		float* xycoords;
		float** xyoldcoords;
		float ox,oy;
		//float** old_coords;
		float*** VoxelMediumRotated;
		float*** voxelDensityRotated;
		float xstepsize,ystepsize;
		float xRight,xLeft,yUp,yLow;
		float xmax,xmin,ymax,ymin;
		char data[1]; //wont be needed in final version
		int indx,indy;
		float sum;
		float r,s,x,y,xNew,yNew;
		
		//make sure angle is in [0,360)
		while (theta > 360)
			theta -=360;
		while (theta < 0)
			theta+=360;

		//we now need to define a new 2d grid (again will be updated to 3d)
			VoxelMediumRotated = new float**[zSize];
			voxelDensityRotated = new float**[zSize];
		    for (k = 0; k<zSize; k++)
		    {
			    VoxelMediumRotated[k] = new float*[ySize];
				voxelDensityRotated[k] = new float*[ySize];
			    for (j = 0; j<ySize; j++)
			    {
				    VoxelMediumRotated[k][j] = new float[xSize];
					voxelDensityRotated[k][j] = new float[xSize];
			    }
		    }



		//idenitify if angle is a special case (that is 90, 180,270, or 360)
		if (!(theta - round(theta))&& !(((int)theta)%90)){//All too easy...

			if(theta == 90){//flip about yaxis and then relfect about diagonal
				printf("90 degrees...\n");
				for(k=0;k<zSize;k++){
				    for (i=0;i<xSize;i++){
					    for(j=0;j<ySize;j++){
						    //printf("(i,j)= (%d,%d)\n",i,j);
						    if((ySize-j-1) < xSize && i <ySize){
							    VoxelMediumRotated[k][j][i] = voxelMedium[k][i][ySize-j-1];
								voxelDensityRotated[k][j][i] = voxelDensity[k][i][ySize-j-1];
						    }
					    }
				    }
				}
			}
			else if (theta == 180){
				printf("180 degrees...\n");
				for(k=0;k<zSize;k++){
				    for (i=0;i<xSize;i++){
					    for(j=0;j<ySize;j++){
						    //printf("(i,j)= (%d,%d)\n",i,j);
						    VoxelMediumRotated[k][j][i] = voxelMedium[k][ySize-j-1][xSize-i-1];
							voxelDensityRotated[k][j][i] = voxelDensity[k][ySize-j-1][xSize-i-1];
					    }
				    }
				}
			}
			else if (theta == 270){//flip about xaxis and reflect about the diagonal
				printf("270 degrees...\n");
				for(k=0;k<zSize;k++){
				    for (i=0;i<xSize;i++){
					    for(j=0;j<ySize;j++){
						    //printf("(i,j)= (%d,%d)\n",i,j);
						    if((xSize-i-1) < ySize && j < xSize){
							    VoxelMediumRotated[k][j][i] = voxelMedium[k][xSize-i-1][j];
								voxelDensityRotated[k][j][i] = voxelDensity[k][xSize-i-1][j];
						    }
					    }
					}
				}
			}
			else //oh hey it must be that theta = n*360 where n is a member of the natural numbers...
			{
				//in simple terms we tell the user that the rotation request is silly (but we say it in a nice way)
				printf("Angles of integer multiples of 360 degrees will result in unchanged image, no work need be done\n");
			}
		}
		else //now we do this the hard way
		{
			//convert degrees to radians in [0, 2*PI)
			rads = ((theta*PI)/180);

			sintheta = sin(rads);
			costheta = cos(rads);

			/*2d rotation matix is
			 *
			 *  cos(theta) -sin(theta)
			 *  sin(theta)  cos(theta)
			 */

			//initialise the xy coords storage
			xycoords = new float[2];   

			xyoldcoords = new float*[(xSize)*(ySize)];
			for(i=0; i<(xSize)*(ySize);i++)
				xyoldcoords[i] = new float[2];

			 
		
			k=0;
			for (i=0;i < xSize;i++){
				for (j=0;j<ySize;j++){
					//printf("do we get here? (i,j)= (%d,%d)\n",i,j);

					//data considered to be at the geometric centre (will not be needed in final version, just for testing)
					xyoldcoords[k][0] = (xBoundaries[xSize-i] + xBoundaries[xSize-i-1])*0.5;
					xyoldcoords[k][1] = (yBoundaries[ySize-j] + yBoundaries[ySize-j-1])*0.5;
					k++;
				}
			}

			// origin point to rotate about
			ox = (xyoldcoords[0][0]+xyoldcoords[xSize*ySize-1][0])*0.5;
			oy = (xyoldcoords[0][1]+xyoldcoords[xSize*ySize-1][1])*0.5;

			xstepsize = (xyoldcoords[0][0] - xyoldcoords[xSize*ySize-1][0])/(xSize-1);
			ystepsize = (xyoldcoords[0][1] - xyoldcoords[xSize*ySize-1][1])/(ySize-1);

			for(k=0;k<zSize;k++){
			    for (i=0;i<xSize;i++){
				    for (j=0;j<ySize;j++){
						VoxelMediumRotated[k][j][i] = 1;
						voxelDensityRotated[k][j][i] = 0;
					    //printf("Whats going on here? (i,j) = (%d,%d)\n",i,j);
					    //grab the points coordinates 
						x = xyoldcoords[i*ySize+j][0];
						y = xyoldcoords[i*ySize+j][1];
					
						//rotate by -theta
						xNew = costheta*(x-ox) +sintheta*(y-oy)+ox;
						yNew = (-sintheta)*(x-ox) + costheta*(y-oy) +oy;
					
					    if (xNew < xyoldcoords[0][0] && yNew < xyoldcoords[0][1] && yNew > xyoldcoords[xSize*ySize-1][1] && xNew > xyoldcoords[xSize*ySize-1][0]){
					        
							//get the surrounding data values from the old grid
						    indx = (int) floor((xyoldcoords[0][0] - xNew)/xstepsize);
						    indy = (int) floor((xyoldcoords[0][1] - yNew)/ystepsize);
						

						    r = xyoldcoords[indx*ySize+indy][0] - xNew;
						    s = xyoldcoords[indx*ySize+indy][1] - yNew;
						    sum =0;
						
							//intepolate the data
						    //printf("indx = %d indy = %d, x = %f y = %f, steps = (%f %f), coords = %f %f,s = %f, r = %f origin = %f %f\n",indx,indy,xNew,yNew,xstepsize,ystepsize,x,y,s,r,ox,oy);
						    data[0] = voxelMedium[k][indy][indx];
						    sum += (1-s)*(1-r)*atof(data);
						    //printf("does it get past here?\n" );
						    data[0] = voxelMedium[k][indy][indx+1];
						    sum += (1-s)*r*atof(data);
						    data[0] = voxelMedium[k][indy+1][indx];
						    sum += s*(1-r)*atof(data);
						    data[0] = voxelMedium[k][indy+1][indx+1];
						    sum += s*r*atof(data);
						
							// assign to the point the interpolated data
						    VoxelMediumRotated[k][j][i] = sum;

							sum =0;
						
							//intepolate the data
						    //printf("indx = %d indy = %d, x = %f y = %f, steps = (%f %f), coords = %f %f,s = %f, r = %f origin = %f %f\n",indx,indy,xNew,yNew,xstepsize,ystepsize,x,y,s,r,ox,oy);
						    sum += (1-s)*(1-r)*voxelDensity[k][indy][indx];
						    //printf("does it get past here?\n" );
						    sum += (1-s)*r*voxelDensity[k][indy][indx+1];
						    sum += s*(1-r)*voxelDensity[k][indy+1][indx];			
						    sum += s*r*voxelDensity[k][indy+1][indx+1];
						
							// assign to the point the interpolated data
						    voxelDensityRotated[k][j][i] = sum;
					    }
					    else{
						    //its out of bounds so who really cares?
					    }
				    }
			    }
			}
		}

		//test output
		std::string outputfile = "test_rotate.txt";
		WriteRotatedData(outputfile.c_str(),VoxelMediumRotated,voxelDensityRotated);
		
		return 0;

	};

	


	/*
	 * rotates volume by theta about an arbitrary vector axis wrt to the mid-point of the grid
	 */
	int EGSPhant::Rotate3D(float theta, float* axis){

		int i,j,k;
		int indx,indy,indz;
		float s,r,t;
		float ox,oy,oz;
		float x,y,z;
		float xr,yr,zr;
		float nx,ny,nz;
		float sintheta,costheta;
		float*** voxelMediumRotated;
		float*** voxelDensityRotated;
		float xStep,yStep,zStep;
		float* xPoints;
		float* yPoints;
		float* zPoints;
		float error,mag;
		char data[1];
		float sum;

		error = 0.00001;

		
		//make sure angle is in [0,360)
		while (theta > 360)
			theta -=360;
		while (theta < 0)
			theta+=360;

		float rads = (theta*PI)/180;

		//store the cosine and sine of this angle
		sintheta = sin(rads);
		costheta = cos(rads);

		//check that axis is a unit vector (if not then normalise)
		if ( fabs(mag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2])-1) > error){

			axis[0] /= mag;
			axis[1] /= mag;
			axis[2] /= mag;
		}

		//we now need to define a new 3d grid of same dimension;		
		voxelMediumRotated = new float**[zSize];		
		voxelDensityRotated = new float**[zSize];
		for (k = 0; k<zSize; k++)
		{
			voxelMediumRotated[k] = new float*[ySize];
			voxelDensityRotated[k] = new float*[ySize];
			for (j = 0; j<ySize; j++)
			{
				voxelMediumRotated[k][j] = new float[xSize];
				voxelDensityRotated[k][j] = new float[xSize];
			}
		}

		//for now we ignore special cases -- though they should be added

		//initialise point arrays
		xPoints = new float[xSize];
		yPoints = new float[ySize];
		zPoints = new float[zSize];



		//gets the midpoints of the voxels which we use for the calculation
		for (i=0;i<xSize;i++)
			xPoints[i] = (xBoundaries[xSize-i] + xBoundaries[xSize-i-1])*0.5;
		for (i=0;i<ySize;i++)
			yPoints[i] = (yBoundaries[ySize-i] + yBoundaries[ySize-i-1])*0.5;
		for (i=0;i<zSize;i++)
			zPoints[i] = (zBoundaries[zSize-i] + zBoundaries[zSize-i-1])*0.5;

		
		//get the mid point (an origin shift)
		ox = (xPoints[0]+xPoints[xSize-1])*0.5;
		oy = (yPoints[0]+yPoints[ySize-1])*0.5;
		oz = (zPoints[0]+zPoints[zSize-1])*0.5;

		//get stepsizes
		xStep = (xPoints[0]-xPoints[xSize-1])/(xSize-1);
		yStep = (yPoints[0]-yPoints[ySize-1])/(ySize-1);
		zStep = (zPoints[0]-zPoints[zSize-1])/(zSize-1);
		printf("delta = %f %f %f\n",xStep,yStep,zStep);
		
		//for each point in the new grid
		for (k=0;k<zSize;k++){
			for (j=0;j<ySize;j++){
				for (i=0;i<xSize;i++){
					
					voxelMediumRotated[k][j][i] = 1;
					voxelDensityRotated[k][j][i] = 0;

					//get the corresponding coordinates of the grid node
					x = xPoints[i];
					y = yPoints[j];
					z = zPoints[k];
					nx = axis[0];
					ny = axis[1];
					nz = axis[2];

					//apply general 3D rotation matrix to this point (can be improved but this is the standard from)
					xr = (1+(1-costheta)*(nx*nx-1))*(x-ox)+(-nz*sintheta+(1-costheta)*nx*ny)*(y-oy)+(ny*sintheta+(1-costheta)*nx*nz)*(z-oz)+ox;
					yr = (nz*sintheta+(1-costheta)*nx*ny)*(x-ox)+(1+(1-costheta)*(ny*ny-1))*(y-oy)+(-nx*sintheta+(1-costheta)*ny*nz)*(z-oz)+oy;
					zr = (-ny*sintheta+(1-costheta)*nx*nz)*(x-ox)+(nx*sintheta+(1-costheta)*ny*nz)*(y-oy)+(1+(1-costheta)*(nz*nz-1))*(z-oz)+oz;

					//check the rotated point is in the grid
					if(xr < xPoints[0] && xr > xPoints[xSize-1] && yr < yPoints[0] && yr > yPoints[ySize-1] && zr < zPoints[0] && zr > zPoints[zSize-1]){

						//now get the indices of the nearest grid point in the upper right voxel corner
						indx = (int)floor((xPoints[0] - xr)/xStep);
						indy = (int)floor((yPoints[0] - yr)/yStep);
						indz = (int)floor((zPoints[0] - zr)/zStep);

						printf("do we get here? i,j,k = %d %d %d, indxyz = %d %d %d\n ",i,j,k,indx,indy,indz);
						//use trilinear interpolation to get the data value for this node
						s = xPoints[indx] - xr;
						r = yPoints[indy] - yr;
						t = zPoints[indz] - zr;

						//only using this data[] to allow atof() to work, wont be needed in final version 
						data[0] = voxelMedium[indz][indy][indx];
						sum += (1-t)*(1-r)*(1-s)*atof(data);
						data[0] = voxelMedium[indz][indy][indx+1];
						sum += (1-t)*(1-r)*(s)*atof(data);
						data[0] = voxelMedium[indz+1][indy][indx];
						sum += (t)*(1-r)*(1-s)*atof(data);
						data[0] = voxelMedium[indz][indy+1][indx];
						sum += (1-t)*(r)*(1-s)*atof(data);
						data[0] = voxelMedium[indz][indy+1][indx+1];
						sum += (1-t)*(r)*(s)*atof(data);
						data[0] = voxelMedium[indz+1][indy][indx+1];
						sum += (t)*(1-r)*(s)*atof(data);
						data[0] = voxelMedium[indz+1][indy+1][indx];
						sum += (t)*(r)*(1-s)*atof(data);
						data[0] = voxelMedium[indz+1][indy+1][indx+1];
						sum += (t)*(r)*(s)*atof(data);

						// assign to the point the interpolated data
						voxelMediumRotated[k][j][i] = sum;

					    sum =0;

						sum += (1-t)*(1-r)*(1-s)*voxelDensity[indz][indy][indx];
						sum += (1-t)*(1-r)*(s)*voxelDensity[indz][indy][indx+1];
						sum += (t)*(1-r)*(1-s)*voxelDensity[indz+1][indy][indx];
						sum += (1-t)*(r)*(1-s)*voxelDensity[indz][indy+1][indx];
						sum += (1-t)*(r)*(s)*voxelDensity[indz][indy+1][indx+1];
						sum += (t)*(1-r)*(s)*voxelDensity[indz+1][indy][indx+1];
						sum += (t)*(r)*(1-s)*voxelDensity[indz+1][indy+1][indx];
						sum += (t)*(r)*(s)*voxelDensity[indz+1][indy+1][indx+1];

						voxelDensityRotated[k][j][i] = sum;
					}
				}
			}
		}

		std::string outputfile = "test_rotate.txt";
		WriteRotatedData(outputfile.c_str(),voxelMediumRotated,voxelDensityRotated);

	};

	/*
	* shifts the grid from by -(x,y,z)
	*/
	int EGSPhant::ShiftToOrigin(float x,float y,float z,int flag){

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

	void EGSPhant::FlipY(){

		int i;

		for (i=0;i<=ySize;i++)
		{
			yBoundaries[i] = -yBoundaries[i];
		}
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

		/*
		for (int j = 0; j<ySize; j++)
		{
			for (int i = 0; i<xSize; i++)
			{
				fprintf(stdout, "%c", voxelMedium[0][j][i]);
			}
			fprintf(stdout, "\n");
		}
		*/		
	};

} //namespace EGS

#endif //MYMATH_MATRIX_H
