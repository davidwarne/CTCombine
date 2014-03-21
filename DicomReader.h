/**
 * DicomReader.h (Definition)
 *
 * @Author Mark Dwyer
 * @Contact m2.dwyer@qut.edu.au
 * @Created 24/10/2008
 *
 * Makes use of the Open Source DicomParser from
 * sourceforge.net
 *
 */
#include "./DICOMParser/DICOMParser.h"
#include "./DICOMParser/DICOMFile.h"
#include "./DICOMParser/DICOMAppHelper.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <dirent.h>
#include <math.h>

#ifndef DICOMREADER_H
#define DICOMREADER_H

const float Pi = 3.14159265358979323846;


namespace DICOM
{
	/**
	 * DicomReader
	 *
	 * Pretty standard implementation
	 */
	
	class DICOMReader
	{
		public:
			// Constructor
			DICOMReader( void );
			DICOMReader( char *directory );

			void PrintAllFiles( void );
			void PrintSlice(int slicenum);
			void ReadSample( void );
			void ReadFile( char* filename );
			void ReadAll( void );

			void OutputSliceToFile(int slice,char* outFileName);
			void FastRotate(float angle, int slice);
			void Rotate3D(float theta, float* axis,float ox,float oy,float oz);
			void FlipAxis(float* axis, int length1);
			void Translate(float x, float y, float z);
			void Select_Region(float xmin,float xmax,float ymin, float ymax, float zmin, float zmax);
			void WriteTestFile();

		public:
			char* directory;
			std::vector<std::string> data_files;
			
			int height;
			int width;
			float *spacing;

			// Data type
			// 0 - Float
			// 1 - Unsigned Char
			// 2 - Short
			// 3 - Unsigned Short
			int imageDataType;

			float ***data_f;
			unsigned char ***data_uc;
			unsigned short ***data_us;
			short ***data_s;

			int numSlices;
			float* xCoords; 
		    float* yCoords;
		    float* zCoords;
	private:
			DICOMParser *parser;
			DICOMAppHelper *helper;

	};

	/****************************
	 *		Constructors
	 ***************************/

	/**
 	 * Default Constructor.
  	 */
	DICOMReader::DICOMReader(void)
	{
		// Nothing to see here - move on.
	};

	/**
 	 * Main Constructor.
	 *
	 * Supplied with a directory, a read of the contents of the directory
	 * will occur.  
	 *
  	 */
	DICOMReader::DICOMReader( char *dir )
	{
		directory = dir;
		parser = new DICOMParser();
		helper = new DICOMAppHelper();
	
		DIR* d = opendir( directory );
    	struct dirent* dirp;

    	if (d)
    	{
        	while ( (dirp = readdir(d)) != NULL )
        	{
            	// check that not current or parent directory
				if ( 0 == strcmp( ".", dirp->d_name ) || 0 == strcmp( "..",dirp->d_name ) || 0 == strcmp( ".DS_Store",dirp->d_name )) 
					continue;

				std::string file(dirp->d_name);
				std::string::size_type loc = file.find(".dcm");
				if (loc!=std::string::npos)
					data_files.push_back(dirp->d_name);
        	}
    	}
    	closedir(d);
		// Sort the directory listing
		std::sort( data_files.begin(), data_files.end() );
    
    	std::vector<std::string>::iterator new_end_pos;
    	new_end_pos = std::unique( data_files.begin(), data_files.end() );

		numSlices = data_files.size();
		fprintf(stdout, "%d number of slices detected.\n", numSlices);
	};

	/**
 	 * Read a sample Dicom file.
	 *
	 * Reads a sample file to get dimensions and allocate memory  
	 *
  	 */
	void DICOMReader::ReadSample( void )
	{
		int i, j, k;

		helper->Clear();
    	parser->ClearAllDICOMTagCallbacks();
		
		char filename[255];
		sprintf(filename, "%s%s", directory, (char *)data_files[0].c_str());
		//fprintf(stdout, "%s\n", filename);
    	parser->OpenFile(filename);
    	helper->RegisterCallbacks(parser);

    	    	
		helper->RegisterPixelDataCallback(parser);
    	parser->ReadHeader();

		void* imgData = NULL;
    	DICOMParser::VRTypes dataType;
    	unsigned long imageDataLength;
		helper->GetImageData(imgData, dataType, imageDataLength);

		height = helper->GetHeight();
		width = helper->GetWidth();
		spacing = helper->GetPixelSpacing();
		fprintf(stdout, "Width = %d\n", width);
		fprintf(stdout, "Height = %d\n", height);

		/* Data is loaded, now to figure out what type it is */
		int bit_depth = helper->GetBitsAllocated();
  		int num_comp = helper->GetNumberOfComponents();
  		bool isFloat = helper->RescaledImageDataIsFloat();
  		bool sgn = helper->RescaledImageDataIsSigned();

		//initialise arrays that store coordinate data
		zCoords = new float[numSlices];
		xCoords = new float[width];
		yCoords = new float[height];


		if (isFloat)
    	{
			// Data is of type float
			fprintf(stdout, "DICOM data of type: float\n");
    		imageDataType = 0;
			data_f = new float**[numSlices]; 
			for (k = 0; k<numSlices; k++)
			{
				data_f[k] = new float*[height];
				for (j = 0; j<height; j++)
				{
					data_f[k][j] = new float[width];
				}
			}
    	}
  		else if (bit_depth <= 8)
    	{
			// Data is of type unsigned char
    		//fprintf(stdout, "Data: Unsigned char\n");
			fprintf(stdout, "DICOM data of type: unsigned char\n");			
			imageDataType = 1;
			// unsigned char
			data_uc = new unsigned char**[numSlices]; 
			for (k = 0; k<numSlices; k++)
			{
				data_uc[k] = new unsigned char*[height];
				for (j = 0; j<height; j++)
				{
					data_uc[k][j] = new unsigned char[width];
				}
			}
			
    	}
  		else
    	{
    		if (sgn)
    	  	{
				// Data is of type short
    	  		//fprintf(stdout, "Data: Short\n");
				imageDataType = 2;
				fprintf(stdout, "DICOM data of type: short\n");
				data_s = new short**[numSlices]; 
				for (k = 0; k<numSlices; k++)
				{
					data_s[k] = new short*[height];
					for (j = 0; j<height; j++)
					{
						data_s[k][j] = new short[width];
					}
				}
    	  	}
    		else
    	  	{
				// Data is of type unsigned short
   		   		//fprintf(stdout, "Data: Unsigned Short\n");
				imageDataType = 3;
				fprintf(stdout, "DICOM data of type: unsigned short\n");
				data_us = new unsigned short**[numSlices]; 
				for (k = 0; k<numSlices; k++)
				{
					data_us[k] = new unsigned short*[height];
					for (j = 0; j<height; j++)
					{
						data_us[k][j] = new unsigned short[width];
					}
				}
   		   	}
    	}
	};


	/**
 	 * Read the whole set.
	 *
	 *   
	 *
  	 */
	void DICOMReader::ReadAll( void )
	{
		int i, j, k;
		void* imgData = NULL;
		DICOMParser::VRTypes dataType;
		unsigned long imageDataLength;
		
		float* xyzStart;

		short *iData_s;
		unsigned short *iData_us;
		unsigned char *iData_uc;
		float *iData_f;
		char filename[255];
		

		if (imageDataType == 2)
		{
			for (k = 0; k<numSlices; k++)
			{
				helper->Clear();
    			parser->ClearAllDICOMTagCallbacks();
				sprintf(filename, "%s%s", directory, (char *)data_files[k].c_str());
    			parser->OpenFile(filename);
    			helper->RegisterCallbacks(parser);
				helper->RegisterPixelDataCallback(parser);
    			parser->ReadHeader();

				imgData = NULL;
				helper->GetImageData(imgData, dataType, imageDataLength);
				xyzStart = helper->GetImagePositionPatient();
				
				zCoords[k] = xyzStart[2];
				// short
				iData_s = (short *)imgData;
				for (j = 0; j<height; j++)
				{
					yCoords[j] = xyzStart[1] + j*spacing[1];
					for (i = 0; i<width; i++)
					{
						if (iData_s[j*width+i]>0){
							data_s[k][j][i] = iData_s[j*width+i];
						}
						else{
							data_s[k][j][i]=0;
						}
						xCoords[i] = xyzStart[0] + i*spacing[0];

						//printf("x,y,z = %f %f %f\n",xCoords[i],yCoords[j],zCoords[k]);
						
					}
					
				}
			}
			
		}
		else if (imageDataType = 1)
		{
			for (k = 0; k<numSlices; k++)
			{
				helper->Clear();
    			parser->ClearAllDICOMTagCallbacks();
    			sprintf(filename, "%s%s", directory, (char *)data_files[k].c_str());
    			parser->OpenFile(filename);
    			helper->RegisterCallbacks(parser);
				helper->RegisterPixelDataCallback(parser);
    			parser->ReadHeader();

				imgData = NULL;
				helper->GetImageData(imgData, dataType, imageDataLength);
				xyzStart = helper->GetImagePositionPatient();
				zCoords[k] = xyzStart[2];
				//printf("Start xyz = %f %f %f\n",xyzStart[0],xyzStart[1],xyzStart[2]);
				// unsigned char
				iData_uc = (unsigned char *)imgData;
				for (j = 0; j<height; j++)
				{
					yCoords[j] = xyzStart[1] + j*spacing[1];
					for (i = 0; i<width; i++)
					{
						if (iData_uc[j*width+i]>0){
							data_uc[k][j][i] = iData_uc[j*width+i];
						}
						else{
							data_uc[k][j][i]=0;
						}
						xCoords[i] = xyzStart[0] + i*spacing[0];

						//printf("x,y,z = %f %f %f\n",xCoords[i],yCoords[j],zCoords[k]);
						
					}
					
				}
			}
		}
		else if (imageDataType == 3)
		{
			for (k = 0; k<numSlices; k++)
			{
				helper->Clear();
    			parser->ClearAllDICOMTagCallbacks();
    			sprintf(filename, "%s%s", directory, (char *)data_files[k].c_str());
    			parser->OpenFile(filename);
    			helper->RegisterCallbacks(parser);
				helper->RegisterPixelDataCallback(parser);
    			parser->ReadHeader();

				imgData = NULL;
				helper->GetImageData(imgData, dataType, imageDataLength);
				xyzStart = helper->GetImagePositionPatient();

				zCoords[k] = xyzStart[2];
				//printf("Start xyz = %f %f %f\n",xyzStart[0],xyzStart[1],xyzStart[2]);

				iData_us = (unsigned short *)imgData;
				for (j = 0; j<height; j++)
				{
					yCoords[j] = xyzStart[1] + j*spacing[1];
					for (i = 0; i<width; i++)
					{
						if (iData_us[j*width+i]>0){
							data_us[k][j][i] = iData_us[j*width+i];
						}
						else{
							data_us[k][j][i]=0;
						}
						xCoords[i] = xyzStart[0] + i*spacing[0];

						//printf("x,y,z = %f %f %f\n",xCoords[i],yCoords[j],zCoords[k]);
						
					}
					
				}
			}
		}
		else if (imageDataType == 0)
		{
			for (k = 0; k<numSlices; k++)
			{
				helper->Clear();
    			parser->ClearAllDICOMTagCallbacks();
    			sprintf(filename, "%s%s", directory, (char *)data_files[k].c_str());
    			parser->OpenFile(filename);
    			helper->RegisterCallbacks(parser);
				helper->RegisterPixelDataCallback(parser);
    			parser->ReadHeader();

				imgData = NULL;
				helper->GetImageData(imgData, dataType, imageDataLength);
				xyzStart = helper->GetImagePositionPatient();

				zCoords[k] = xyzStart[2];
				//printf("Start xyz = %f %f %f\n",xyzStart[0],xyzStart[1],xyzStart[2]);

				iData_f = (float *)imgData;
				for (j = 0; j<height; j++)
				{
					yCoords[j] = xyzStart[1] + j*spacing[1];
					for (i = 0; i<width; i++)
					{
						if (iData_s[j*width+i]>0){
							data_f[k][j][i] = iData_f[j*width+i];
						}
						else{
							data_f[k][j][i]=0;
						}
						xCoords[i] = xyzStart[0] + i*spacing[0];

						//printf("x,y,z = %f %f %f\n",xCoords[i],yCoords[j],zCoords[k]);
						
					}
					
				}
			}
		}	
		else
		{
			fprintf(stdout, "Data type in Dicom not handled.  Not Short, unsigned short, float or unsigned char.");
		}

	};


	/**
 	 * Read a Dicom file.
	 *
	 * Supplied with a file, a read of the contents of the file
	 * will occur.  
	 *
  	 */
	void DICOMReader::ReadFile( char* filename )
	{
		helper->Clear();
    	parser->ClearAllDICOMTagCallbacks();

    	parser->OpenFile(filename);
    	helper->RegisterCallbacks(parser);

    	    	
		helper->RegisterPixelDataCallback(parser);
    	parser->ReadHeader();

		void* imgData = NULL;
    	DICOMParser::VRTypes dataType;
    	unsigned long imageDataLength;
		helper->GetImageData(imgData, dataType, imageDataLength);
  		
	};

	/**
 	 * Print all files in the directory
	 *
	 * Essentially a tester method
	 *
  	 */
	void DICOMReader::PrintAllFiles( void )
	{
		int i;
		for (i = 0; i<data_files.size(); i++)
		{
			fprintf(stdout, "%s\n", data_files[i].c_str());
		}
	};

	void DICOMReader::OutputSliceToFile(int slice, char* outfileName)
	{
		if (slice < 0 || slice > (numSlices-1))
		{
			fprintf(stdout, "Slice out of bounds.\n");
			return;
		}
		FILE *file;
		
		int i, j, k;
		char filename[255];
		sprintf(filename, "%s_%d.txt",outfileName, slice);		
		if((file=fopen(filename,"w"))==NULL)
			fprintf(stdout, "Error opening file\n");

		if ( imageDataType == 2 )
		{
			fprintf(stdout, "Writing: short\n");
			for (j = 0; j<height; j++)
			{
				for (i = 0; i<width; i++)
				{
					fprintf(file, "%d ", data_s[slice][j][i]);
				}
				fprintf(file, "\n");
			}
		}
		else if ( imageDataType == 1 )
		{
			fprintf(stdout, "Writing: unsigned char\n");
			for (j = 0; j<height; j++)
			{
				for (i = 0; i<width; i++)
				{
					fprintf(file, "%d ", data_uc[slice][j][i]);
				}
				fprintf(file, "\n");
			}
		}
		else if ( imageDataType == 3 )
		{
			fprintf(stdout, "Writing: unsigned short\n");
			for (j = 0; j<height; j++)
			{
				for (i = 0; i<width; i++)
				{
					fprintf(file, "%d ", data_us[slice][j][i]);
				}
				fprintf(file, "\n");
			}
		}
		else if( imageDataType == 0 )
		{
			fprintf(stdout, "Writing: float\n");
			for (j = 0; j<height; j++)
			{
				for (i = 0; i<width; i++)
				{
					fprintf(file, "%f ", data_f[slice][j][i]);
				}
				fprintf(file, "\n");
			}
		}
		else
			fprintf(stdout, "Data type in Dicom not handled.  Not Short, unsigned short, float or unsigned char.");
		fclose(file);
	};

	/*
	 *	FastRotate
	 *
	 *  Rotates a slice
	 */
	void DICOMReader::FastRotate(float angle, int slice)
	{
		// need to ensure that angle is [0, 360)
		while (angle < 0)
			angle+=360.0;
		while (angle > 360.0)
			angle-=360.0;

		short * dest = new short[width*height];
	
		// Here need to insert the 0,90,180,360 special cases.
		const float rad = (float) ((angle*Pi)/180.0);
		const float	ca=(float)cos(rad);
		const float sa=(float)sin(rad);

		const float ux  = (float)(fabs(width*ca));
		const float uy  = (float)(fabs(width*sa));
	  	const float vx  = (float)(fabs(height*sa));
		const float vy  = (float)(fabs(height*ca));
	  	const float w2  = 0.5f*(width-1);
		const float h2  = 0.5f*(height-1);
	  	const float dw2 = 0.5f*(ux+vx);
		const float dh2 = 0.5f*(uy+vy); // dw2, dh2 are the dimentions for rotated image without cropping.
	
		int X,Y; // Locations in the source matrix
		float X_f, Y_f;
		int x,y; // For counters
		int image_offset;
		int slice_offset = slice*width*height;
		float temp_y, temp_x;	
		int numPixels = width*height;
		int temp_value;

		int sum = 0;
		int average;
		int counter = 0;
		int threshold = -750;
		int c_x, c_y, f_x, f_y, r_x, r_y;
		/*
		 *  Really need to check for the multiple of 90 degree cases.
		 */
		if (angle == 180.0)
		{
			// Nothing much needs to be done!
			// Just populate the destination array by transposition 
			for (y = 0; y<height; y++)
			{
				image_offset = y*width;
				// flip in x axis
				for (x = 0; x<width; x++)
				{
					dest[image_offset+x] = data_s[slice][height-1-y][width-1-x];
				}
			}
			
		}
		else if (angle == 360.0 || angle == 0.0)
		{
			// Nothing needs to be done!
			// Just populate the destination array
			for (y = 0; y<height; y++)
			{
				image_offset = y*width;
				for (x = 0; x<width; x++)
				{
					dest[image_offset+x] = data_s[slice][y][x];
				}
			}
			
		}
		else
		{
			// not a multiple of 90 degrees - phew!
			for(y=0;y<height;y++)
			{
				image_offset=y*width;
				temp_y = y;
				for(x=0;x<width;x++)
				{
					temp_x = x;
					X_f = (w2 + (temp_x-w2)*ca + (temp_y-h2)*sa+0.5); // Source X
					Y_f = (h2 - (temp_x-w2)*sa + (temp_y-h2)*ca+0.5); // Source Y
					
					r_x = (int)round(X_f);
					r_y = (int)round(Y_f);

					//Method 2 - better success but inaccurate with single degree rotations					
					sum = 0;
					counter = 0;	
					if (r_x < 1 || r_y < 1 || r_x >= (width-1) || r_y >= (height-1))
					{
						dest[image_offset+x] = -1024;   
					}
					else
					{
						sum+=4*data_s[slice][r_y][r_x];
						sum+=data_s[slice][r_y+1][r_x];
						sum+=data_s[slice][r_y-1][r_x];
						sum+=data_s[slice][r_y][r_x+1];
						sum+=data_s[slice][r_y][r_x-1];
						dest[image_offset+x] = (int)round((float)sum/8.0);
					}
					

					/*  First attempt - mostly a fail.
					if (c_x<0 || c_y<0 || c_x>=width || c_y>=height)
					{

					}
					else
					{
						// both ceilings
						temp_value = data_s[slice][c_y][c_x];
						if (temp_value > threshold)
						{
							sum+=temp_value;
							counter++;
						}
					}
						
					
					if (f_x<0 || f_y<0 || f_x>=width || f_y>=height)
					{
						
					}	
					else
					{
						// both floors
						temp_value = data_s[slice][f_y][f_x];
						if (temp_value > threshold)
						{
							sum+=temp_value;
							counter++;
						}
					}				
					if (c_x<0 || f_y<0 || c_x>=width || f_y>=height)
					{
						
					}
					else
					{
						// ceil X and floor Y
						temp_value = data_s[slice][f_y][c_x];
						if (temp_value > threshold)
						{
							sum+=temp_value;
							counter++;
						}
					}	
					if (f_x<0 || c_y<0 || f_x>=width || c_y>=height)
					{
						
					}
					else
					{
						// floor X and ceil Y
						temp_value = data_s[slice][c_y][f_x];
						if (temp_value > threshold)
						{
							sum+=temp_value;
							counter++;
						}
					}
					if (counter > 0)
					{	
						average = (int)round((float)sum/(float)counter);
						dest[image_offset+x] = (short)average;
					}
					else
						dest[image_offset+x] = -1024;
					*/
					
					X = (int)round(X_f);
					Y = (int)round(Y_f);				
					dest[image_offset+x] = (X<0 || Y<0 || X>=width || Y>=height)?-1024:data_s[slice][Y][X];
				}
			}
		}
	
		
		FILE *file;
		
		int i, j, k;
		char filename[255];
		sprintf(filename, "Rotate_%d.txt", slice);		
		if((file=fopen(filename,"w"))==NULL)
			fprintf(stdout, "Error opening file\n");

		for (j = 0; j<height; j++)
		{
			for (i = 0; i<width; i++)
			{
				fprintf(file, "%d ", dest[j*width+i]);
			}
			fprintf(file, "\n");
			//fprintf(stdout, "%d\n", j);
		}
		fclose(file);
		
	};

	//rotates the volume by theta degrees about arbitrary axis
	//
	//NOTE: should add the ability to perform multilp rotation in the one function
	void DICOMReader::Rotate3D(float theta, float* axis,float ox,float oy,float oz){

		int i,j,k;
		int indx,indy,indz;
		float s,r,t;
		//float ox,oy,oz;
		float x,y,z;
		float xr,yr,zr;
		float nx,ny,nz;
		short*** newData_s;
		unsigned char*** newData_uc;
		unsigned short*** newData_us;
		float*** newData_f;
		float sintheta,costheta;
		float xStep,yStep,zStep;
		float* xPoints;
		float* yPoints;
		float* zPoints;
		float error,mag;
		float sum,rads;

		// need to ensure that angle is [0, 360)
		while (theta < 0)
			theta+=360.0;
		while (theta > 360.0)
			theta-=360.0;

		//special angle cases should be added

		//for now we assume the data is of type short
		newData_s = new short**[numSlices];
		for (k=0;k<numSlices;k++){
			newData_s[k] = new short*[height];
			for(j=0;j<height;j++){
				newData_s[k][j] = new short[width];
			}
		}

		

			rads = (theta*PI)/180;

		
			sintheta = sin(rads);
			costheta = cos(rads);
			//printf("%f %f \n",sintheta,costheta);

			//check that axis is a unit vector (if not then normalise)
			if ( fabs(mag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2])-1) > error){

				axis[0] /= mag;
				axis[1] /= mag;
				axis[2] /= mag;
			}

			nx = axis[0];
			ny = axis[1];
			nz = axis[2];

			//get the midpoint of rotation if none has been specified
			if (!ox && !oy && !oz ){
			  ox = (xCoords[0] + xCoords[width-1])*0.5;
			  oy = (yCoords[0] + yCoords[height-1])*0.5;
			  oz = (zCoords[0] + zCoords[numSlices-1])*0.5;
			}

			xStep = (xCoords[0] - xCoords[width-1])/(width-1);
			yStep = (yCoords[0] - yCoords[height-1])/(height-1);
			zStep = (zCoords[0] - zCoords[numSlices-1])/(numSlices-1);

			if ( imageDataType == 2 )//short
			{
				for (k=0;k<numSlices;k++){
					for (j=0;j<height;j++){
						for (i=0;i<width;i++){

							newData_s[k][j][i] =0;

							x = xCoords[i];
							y = yCoords[j];
							z = zCoords[k];

							//apply general 3D rotation matrix to this point (can be improved but this is the standard from)
							xr = (1+(1-costheta)*(nx*nx-1))*(x-ox)+(-nz*sintheta+(1-costheta)*nx*ny)*(y-oy)+(ny*sintheta+(1-costheta)*nx*nz)*(z-oz)+ox;
							yr = (nz*sintheta+(1-costheta)*nx*ny)*(x-ox)+(1+(1-costheta)*(ny*ny-1))*(y-oy)+(-nx*sintheta+(1-costheta)*ny*nz)*(z-oz)+oy;
							zr = (-ny*sintheta+(1-costheta)*nx*nz)*(x-ox)+(nx*sintheta+(1-costheta)*ny*nz)*(y-oy)+(1+(1-costheta)*(nz*nz-1))*(z-oz)+oz;
							//TODO: add special cases for rotaion about coordinate axes
					
							//now get the indices of the nearest grid point 
								indx = (int)floor((xCoords[0] - xr)/xStep);
								indy = (int)floor((yCoords[0] - yr)/yStep);
								indz = (int)floor((zCoords[0] - zr)/zStep);

							//check the point in within the grid bounds
							if(indx >=0 && indx < (width-1) &&  indy >=0 && indy < (height-1) && indz >=0 && indz < (numSlices-1)){

								

								//printf("do we get here? i,j,k = %d %d %d, indxyz = %d %d %d\n ",i,j,k,indx,indy,indz);
								//use trilinear interpolation to get the data value for this node
								s = xCoords[indx] - xr;
								r = yCoords[indy] - yr;
								t = zCoords[indz] - zr;
								//printf("%f %f %f\n",s,r,t);

								sum =0;

								sum += (1-t)*(1-r)*(1-s)*((float)data_s[indz][indy][indx]);
								sum += (1-t)*(1-r)*(s)*((float)data_s[indz][indy][indx+1]);
								sum += (t)*(1-r)*(1-s)*((float)data_s[indz+1][indy][indx]);
								sum += (1-t)*(r)*(1-s)*((float)data_s[indz][indy+1][indx]);
								sum += (1-t)*(r)*(s)*((float)data_s[indz][indy+1][indx+1]);
								sum += (t)*(1-r)*(s)*((float)data_s[indz+1][indy][indx+1]);
								sum += (t)*(r)*(1-s)*((float)data_s[indz+1][indy+1][indx]);
								sum += (t)*(r)*(s)*((float)data_s[indz+1][indy+1][indx+1]);
								//printf("%f\n",sum);

								newData_s[k][j][i] = (short)round(sum);
								//printf("%d\n",newData_s[k][j][i]);

							}
						}
					}
				}

				data_s = newData_s;
			}

			//printf("(%f,%f), (%f,%f) (%f,%f)\n",xCoords[0],xCoords[width-1],yCoords[0],yCoords[height-1],zCoords[0],zCoords[numSlices-1]);
	};

	void FlipAxis(float* Axis, int Length){

		float temp;
		int i;
		for (i=0;i<Length;i++){
			temp = Axis[i];
			Axis[i] = Axis[Length-1-i];
			Axis[Length-1-i] = temp;
		}

	};

	void DICOMReader::Translate(float x, float y, float z){

		int i;

		for(i=0;i<width;i++)
			xCoords[i] -= x;
		for(i=0;i<height;i++)
			yCoords[i] -= y;
		for(i=0;i<numSlices;i++)
			zCoords[i] -= z;
	};

	void DICOMReader::WriteTestFile(){

		int i,j,k;
		FILE *file;

		if(!(file = fopen("DICOMtestFile.txt","w")))
			printf("Error opening file\n");

		fprintf(file,"%d %d %d\n",width,height,numSlices);
		if ( imageDataType == 2 )
		{
			fprintf(stdout, "Writing: short\n");
			for(k=0;k<numSlices;k++){
				for (j = 0; j<height; j++)
				{
				    for (i = 0; i<width; i++)
				    {
					     fprintf(file, "%d ", data_s[k][j][i]);
				    }
				    fprintf(file, "\n");
			    }
				fprintf(file,"\n");
			}
		}
		else if ( imageDataType == 1 )
		{
			fprintf(stdout, "Writing: unsigned char\n");
			for(k=0;k<numSlices;k++){
				for (j = 0; j<height; j++)
				{
				    for (i = 0; i<width; i++)
				    {
					     fprintf(file, "%d ", data_uc[k][j][i]);
				    }
				    fprintf(file, "\n");
			    }
				fprintf(file,"\n");
			}
		}
		else if ( imageDataType == 3 )
		{
			fprintf(stdout, "Writing: unsigned short\n");
			for(k=0;k<numSlices;k++){
				for (j = 0; j<height; j++)
				{
				    for (i = 0; i<width; i++)
				    {
					     fprintf(file, "%d ", data_us[k][j][i]);
				    }
				    fprintf(file, "\n");
			    }
				fprintf(file,"\n");
			}
		}
		else if( imageDataType == 0 )
		{
			fprintf(stdout, "Writing: float\n");
			for(k=0;k<numSlices;k++){
				for (j = 0; j<height; j++)
				{
				    for (i = 0; i<width; i++)
				    {
					     fprintf(file, "%f ", data_f[k][j][i]);
				    }
				    fprintf(file, "n");
			    }
				fprintf(file,"\n");
			}
		}
		else
			fprintf(stdout, "Data type in Dicom not handled.  Not Short, unsigned short, float or unsigned char.");
		fclose(file);

	};

	// removes all data outside the given bounds
	void DICOMReader::Select_Region(float xmin,float xmax,float ymin, float ymax, float zmin, float zmax){

		short*** temp;
		float xStep,yStep,zStep;
		int xStart,yStart,zStart;
		int xSize,ySize,zSize;
		float* xtemp;
		float* ytemp;
		float* ztemp;
		int k,j,i;

		xStep = (xCoords[0] - xCoords[width-1])/(width-1);
		yStep = (yCoords[0] - yCoords[height-1])/(height-1);
		zStep = (zCoords[0] - zCoords[numSlices-1])/(numSlices-1);
		//printf("%f %f %f\n",xStep,yStep,zStep);
		//printf("%f %f %f\n",xCoords[0],yCoords[0],zCoords[0]);
		//printf("%f %f %f\n",xCoords[width-1],yCoords[height-1],zCoords[numSlices-1]);
		
		//initialise memory for selected region
		xSize = (int)floor((xmin - xmax)/xStep);
		ySize = (int)floor((ymin - ymax)/yStep);
		zSize = (int)floor((zmin - zmax)/zStep);
		//printf("%d %d %d\n",width,height,numSlices);
		//printf("%d %d %d\n",xSize,ySize,zSize);

		temp = new short**[zSize];
		for (k=0;k<zSize;k++){
			temp[k] = new short*[ySize];
			for(j=0;j<ySize;j++){
				temp[k][j] = new short[xSize];
			}
		}

		xtemp = new float[xSize];
		ytemp = new float[ySize];
		ztemp = new float[zSize];

		xStart = (int)floor(fabs((xCoords[0] - xmin)/xStep));
		yStart = (int)floor(fabs((yCoords[0] - ymin)/yStep));
		zStart = (int)floor(fabs((zCoords[0] - zmin)/zStep));
		//printf("%d %d %d\n",xStart,yStart,zStart);

		for (k=0;k<zSize;k++){
			for (j=0;j<ySize;j++){
				for (i=0;i<xSize;i++){
					temp[k][j][i] = data_s[zStart+k][yStart+j][xStart+i];
					xtemp[i] = xCoords[i];
				}
				ytemp[j] = yCoords[j];
			}
			ztemp[k] = zCoords[k];
		}

		width = xSize;
		height = ySize;
		numSlices = zSize;
		printf("%d %d %d\n",width,height,numSlices);


		data_s = temp;
		xCoords = xtemp;
		yCoords = ytemp;
		zCoords = ztemp;
	};

} //namespace DICOM

#endif //DICOMREADER_H
