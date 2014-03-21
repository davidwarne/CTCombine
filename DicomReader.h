/**
 * DicomReader.h (Definition)
 *
 * @Author: Mark Dwyer, David Warne
 * @Contact: m2.dwyer@qut.edu.au, david.warne@qut.edu.au
 * @Created: 24/10/2008
 * @Last Modified: 12/03/2009
 *
 * Makes use of the Open Source DicomParser from
 * sourceforge.net
 *
 * NOTE: currently this file is in the process of re-engineering
 *
 * Summary:
 *   Reads DICOM data and contains functions for rotations, translations, and data region selestion
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
    
            // file read Functions
			
			void ReadSample( void );
			void ReadAll( void );
			
			// file write Functions
            void OutputSliceToFile(int slice,char* outFileName);
            void WriteTestFile();
            
            // print to stdout Functions
            void PrintAllFiles( void );
			void PrintSlice(int slicenum);
            
            //Data manipulation functions
			void Rotate3D(float theta, float* axis,float ox,float oy,float oz);
			void Translate(float x, float y, float z);
			int Select_Region(float xmin,float xmax,float ymin, float ymax, float zmin, float zmax);
			void GetMidPoint(float* x,float* y,float* z);

		public:
			char* directory;
			std::vector<std::string> data_files;
			
			// grid dimensions z*y*x
			int numSlices;
			int height;
			int width;
			
			float xlims[2];
			float ylims[2];
			float zlims[2];
			
			// x,y,z coordinate spacing
			float *spacing;
			float mx,my,mz; // to store the original mid-point
			// Data type
			// 0 - Float
			// 1 - Unsigned Char
			// 2 - Short
			// 3 - Unsigned Short
			int imageDataType;
			
			// data shift value
			float offset; 
			float slope;
			
			// holds directional cosines of the rows (X) and
			// columns (Y)
			float XdirCosine[3];
			float YdirCosine[3];
			
			// to store the data			
			float ***data_f;
			unsigned char ***data_uc;
			unsigned short ***data_us;
			short ***data_s;
			
			float*** xCoords; 
		    float*** yCoords;
		    float*** zCoords;
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
 	 * ReadSample: Read a sample Dicom file.
	 *
	 * Reads a sample file to get dimensions and allocate memory  
	 *
  	 */
	void DICOMReader::ReadSample( void )
	{
		int i, j, k;
        void* imgData = NULL;
        unsigned long imageDataLength;
        std::vector<std::string> v;
		float* Orientation;
		char filename[255];
		float* xyzStart;
		DICOMParser::VRTypes dataType;
				
		helper->Clear();
    	parser->ClearAllDICOMTagCallbacks();
		
		sprintf(filename, "%s%s", directory, (char *)data_files[0].c_str());
		
    	parser->OpenFile(filename);
    	helper->RegisterCallbacks(parser);    	
		helper->RegisterPixelDataCallback(parser);
    	parser->ReadHeader();
 
		offset = helper->GetRescaleOffset();
		slope = helper->GetRescaleSlope();
		printf("Rescale offset and slope: %f %f \n",offset,slope);
	 	
		helper->GetImageData(imgData, dataType, imageDataLength);
		height = helper->GetHeight();
		width = helper->GetWidth();
		spacing = helper->GetPixelSpacing();
		
		fprintf(stdout, "Width = %d\n", width);
		fprintf(stdout, "Height = %d\n", height);
				
		Orientation = helper->GetOrientation();
		
		printf("XY directional Cosines: (%f %f %f), (%f %f %f)\n",Orientation[0],Orientation[1],Orientation[2],Orientation[3],Orientation[4],Orientation[5]);
		XdirCosine[0] = (float)Orientation[0];
		XdirCosine[1] = (float)Orientation[1];
		XdirCosine[2] = (float)Orientation[2];
		YdirCosine[0] = (float)Orientation[3];
		YdirCosine[1] = (float)Orientation[4];
		YdirCosine[2] = (float)Orientation[5];
		//Data is loaded, now to figure out what type it is 
		int bit_depth = helper->GetBitsAllocated();
  		int num_comp = helper->GetNumberOfComponents();
  		bool isFloat = helper->RescaledImageDataIsFloat();
  		bool sgn = helper->RescaledImageDataIsSigned();
  		xyzStart = helper->GetImagePositionPatient();
  		
  	
  		xlims[0] = xyzStart[0];
  		xlims[1] = xyzStart[0];
  		ylims[0] = xyzStart[1];
  		ylims[1] = xyzStart[1];
  		zlims[0] = xyzStart[2];
  		zlims[1] = xyzStart[2];
  		
		//initialise arrays that store coordinate data
		xCoords = new float**[numSlices];
		yCoords = new float**[numSlices];
		zCoords = new float**[numSlices];
		for (k=0; k<numSlices; k++)
		{
		    xCoords[k] = new float*[height];
		    yCoords[k] = new float*[height];
		    zCoords[k] = new float*[height];
		    for (j=0; j<height; j++)
		    {
		        xCoords[k][j] = new float[width];
		        yCoords[k][j] = new float[width];
		        zCoords[k][j] = new float[width];
		    } 
		}

        // allocate the array for data of the given type
		if (isFloat)
    	{
			// Data is of type float
			fprintf(stdout, "DICOM data of type: float\n");
    		imageDataType = 0;
			data_f = new float**[numSlices]; 
			for (k=0;k<numSlices;k++)
			{
			    data_f[k] = new float*[height];
			    for (j=0;j<height;j++)
			    {
			        data_f[k][j] = new float[width];
			    }
			}
    	}
  		else if (bit_depth <= 8)
    	{
			// Data is of type unsigned char
    		fprintf(stdout, "DICOM data of type: unsigned char\n");			
			imageDataType = 1;
			data_uc = new unsigned char**[numSlices]; 
			for (k=0;k<numSlices;k++)
			{
			    data_uc[k] = new unsigned char*[height];
			    for (j=0;j<height;j++)
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
				imageDataType = 2;
				fprintf(stdout, "DICOM data of type: short\n");
				data_s = new short**[numSlices]; 
			    for (k=0;k<numSlices;k++)
			    {
			        data_s[k] = new short*[height];
			        for (j=0;j<height;j++)
			         {
			            data_s[k][j] = new short[width];
			         }
			    }
    	  	}
    		else
    	  	{
				// Data is of type unsigned short
				imageDataType = 3;
				fprintf(stdout, "DICOM data of type: unsigned short\n");
				data_us = new unsigned short**[numSlices]; 
			    for (k=0;k<numSlices;k++)
			    {
			        data_us[k] = new unsigned short*[height];
			        for (j=0;j<height;j++)
			        {
			            data_us[k][j] = new unsigned short[width];
			        }
			    }
   		   	}
    	}
	};


	/**
 	 * ReadAll: Read the whole set.
	 *
	 * NOTE: if ReadSample has not been called this will crash
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
		
		//check what data type was recorded
		if (imageDataType == 2)
		{
		    //type was short
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
				
				//printf("%f %f %f\n",xyzStart[0],xyzStart[1],xyzStart[2]);

				iData_s = (short *)imgData;
				for (j = 0; j<height; j++)
				{
					for (i = 0; i<width; i++)
					{
				        // the data      
						data_s[k][j][i] = ((short)slope)*iData_s[j*width+i] - (short)offset;
                       
                        //the coordinates in room coordinates
						xCoords[k][j][i] = XdirCosine[0]*((float)i)*spacing[0] + YdirCosine[0]*((float)j)*spacing[1] + xyzStart[0];
						yCoords[k][j][i] = XdirCosine[1]*((float)i)*spacing[0] + YdirCosine[1]*((float)j)*spacing[1] + xyzStart[1];
						zCoords[k][j][i] = XdirCosine[2]*((float)i)*spacing[0] + YdirCosine[2]*((float)j)*spacing[1] + xyzStart[2];
					
					    // updata bounds info
						if (xCoords[k][j][i] < xlims[0])
						{
						    xlims[0] = xCoords[k][j][i];
						}
						if(xCoords[k][j][i] > xlims[1])
						{
						    xlims[1] = xCoords[k][j][i];
						}
						
						if (yCoords[k][j][i] < ylims[0])
						{
						    ylims[0] = yCoords[k][j][i];
						}
						if(yCoords[k][j][i] > ylims[1])
						{
						    ylims[1] = yCoords[k][j][i];
						}
						
						if (zCoords[k][j][i] < zlims[0])
						{
						    zlims[0] = zCoords[k][j][i];
						}
						if(zCoords[k][j][i] > zlims[1])
						{
						    zlims[1] = zCoords[k][j][i];
						}
					}
					
				}
			}
			
			//get the midpoint of the grid 	
			mx = xlims[0] + (xlims[1] - xlims[0])*0.5;
			my = ylims[0] + (ylims[1] - ylims[0])*0.5;
			mz = zlims[0] + (zlims[1] - zlims[0])*0.5;
			
		}
		else if(imageDataType == 0)
		{
		     //type was float
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
				
				//printf("%f %f %f\n",xyzStart[0],xyzStart[1],xyzStart[2]);

				iData_f = (float *)imgData;
				for (j = 0; j<height; j++)
				{
					for (i = 0; i<width; i++)
					{
				        // the data
				       
						data_f[k][j][i] = ((float)slope)*iData_f[j*width+i] - (float)offset;
              
                        //the coordinates in room coordinates
						xCoords[k][j][i] = XdirCosine[0]*((float)i)*spacing[0] + YdirCosine[0]*((float)j)*spacing[1] + xyzStart[0];
						yCoords[k][j][i] = XdirCosine[1]*((float)i)*spacing[0] + YdirCosine[1]*((float)j)*spacing[1] + xyzStart[1];
						zCoords[k][j][i] = XdirCosine[2]*((float)i)*spacing[0] + YdirCosine[2]*((float)j)*spacing[1] + xyzStart[2];
						
						// updata bounds info
						if (xCoords[k][j][i] < xlims[0])
						{
						    xlims[0] = xCoords[k][j][i];
						}
						if(xCoords[k][j][i] > xlims[1])
						{
						    xlims[1] = xCoords[k][j][i];
						}
						
						if (yCoords[k][j][i] < ylims[0])
						{
						    ylims[0] = yCoords[k][j][i];
						}
						if(yCoords[k][j][i] > ylims[1])
						{
						    ylims[1] = yCoords[k][j][i];
						}
						
						if (zCoords[k][j][i] < zlims[0])
						{
						    zlims[0] = zCoords[k][j][i];
						}
						if(zCoords[k][j][i] > zlims[1])
						{
						    zlims[1] = zCoords[k][j][i];
						}
					}
					
				}
			}
			
			//get the midpoint of the grid 	
			mx = xlims[0] + (xlims[1] - xlims[0])*0.5;
			my = ylims[0] + (ylims[1] - ylims[0])*0.5;
			mz = zlims[0] + (zlims[1] - zlims[0])*0.5;
		}
		else if (imageDataType == 1)
		{
		     //type was unsigned char
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
				
				

				iData_uc = (unsigned char *)imgData;
				for (j = 0; j<height; j++)
				{
					for (i = 0; i<width; i++)
					{
				        // the data
				       
						data_uc[k][j][i] = ((unsigned char)slope)*iData_uc[j*width+i] - (unsigned char)offset;
                       
                        //the coordinates in room coordinates
						xCoords[k][j][i] = XdirCosine[0]*((float)i)*spacing[0] + YdirCosine[0]*((float)j)*spacing[1] + xyzStart[0];
						yCoords[k][j][i] = XdirCosine[1]*((float)i)*spacing[0] + YdirCosine[1]*((float)j)*spacing[1] + xyzStart[1];
						zCoords[k][j][i] = XdirCosine[2]*((float)i)*spacing[0] + YdirCosine[2]*((float)j)*spacing[1] + xyzStart[2];
						
						// updata bounds info
						if (xCoords[k][j][i] < xlims[0])
						{
						    xlims[0] = xCoords[k][j][i];
						}
						if(xCoords[k][j][i] > xlims[1])
						{
						    xlims[1] = xCoords[k][j][i];
						}
						
						if (yCoords[k][j][i] < ylims[0])
						{
						    ylims[0] = yCoords[k][j][i];
						}
						if(yCoords[k][j][i] > ylims[1])
						{
						    ylims[1] = yCoords[k][j][i];
						}
						
						if (zCoords[k][j][i] < zlims[0])
						{
						    zlims[0] = zCoords[k][j][i];
						}
						if(zCoords[k][j][i] > zlims[1])
						{
						    zlims[1] = zCoords[k][j][i];
						}
					}
					
				}
			}
			
			//get the midpoint of the grid 	
			mx = xlims[0] + (xlims[1] - xlims[0])*0.5;
			my = ylims[0] + (ylims[1] - ylims[0])*0.5;
			mz = zlims[0] + (zlims[1] - zlims[0])*0.5;
		}
		else if(imageDataType == 3)
		{
		     //type was unsigned short
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
				
				

				iData_us = (unsigned short *)imgData;
				for (j = 0; j<height; j++)
				{
					for (i = 0; i<width; i++)
					{
				        // the data
				       
						data_us[k][j][i] = ((unsigned short)slope)*iData_us[j*width+i] - (unsigned short)offset;
                        
                        //the coordinates in room coordinates
						xCoords[k][j][i] = XdirCosine[0]*((float)i)*spacing[0] + YdirCosine[0]*((float)j)*spacing[1] + xyzStart[0];
						yCoords[k][j][i] = XdirCosine[1]*((float)i)*spacing[0] + YdirCosine[1]*((float)j)*spacing[1] + xyzStart[1];
						zCoords[k][j][i] = XdirCosine[2]*((float)i)*spacing[0] + YdirCosine[2]*((float)j)*spacing[1] + xyzStart[2];
						
						// updata bounds info
						if (xCoords[k][j][i] < xlims[0])
						{
						    xlims[0] = xCoords[k][j][i];
						}
						if(xCoords[k][j][i] > xlims[1])
						{
						    xlims[1] = xCoords[k][j][i];
						}
						
						if (yCoords[k][j][i] < ylims[0])
						{
						    ylims[0] = yCoords[k][j][i];
						}
						if(yCoords[k][j][i] > ylims[1])
						{
						    ylims[1] = yCoords[k][j][i];
						}
						
						if (zCoords[k][j][i] < zlims[0])
						{
						    zlims[0] = zCoords[k][j][i];
						}
						if(zCoords[k][j][i] > zlims[1])
						{
						    zlims[1] = zCoords[k][j][i];
						}
					}
					
				}
			}
			
			//get the midpoint of the grid 	
			mx = xlims[0] + (xlims[1] - xlims[0])*0.5;
			my = ylims[0] + (ylims[1] - ylims[0])*0.5;
			mz = zlims[0] + (zlims[1] - zlims[0])*0.5;
		}
		else
		{
			fprintf(stdout, "Data type in Dicom not handled.  Not Short, unsigned short, float or unsigned char.");
		}

	};

	/**
 	 * PrintAllFiles: Print all files in the directory
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

	/**
	 *OutpuSliceToFile: Writes a selected slice to file
	 * 
	 * Used for testing 2D rotation before extending to 3D
	 */
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

	/**
	 * rotates the volume by theta degrees about arbitrary axis, centred at the give point (ox,oy,oz)
	 *
	 * Axis is a vector parallel to the axis of rotation
	 * if flag is true then the user has selected to use another isocentre
	 * 
	 */
	void DICOMReader::Rotate3D(float theta, float* axis,float ox,float oy,float oz)
	{
		int i,j,k;                      //for loops
		float x,y,z;                    // coordinate of original point
		float xr,yr,zr;                 // coordinate of rotated point
		float nx,ny,nz;                 // the nomalised vector components of the axis of rotation
		float tx,ty,tz;                 // temps to hold tho original centre (we need it for bounds checking)
		float sintheta,costheta;        // to store sine and cosine of theta
		float tol;                      // tolerance level
		float mag;                      // vector magnitude
		float sum;                      // for weigthed sum of intepolation
		float rads;                     // theta in radians

		tol = 0.0000001;

		// need to ensure that angle is [0, 360)
		while (theta < 0)
		{
			theta+=360.0;
		}

		while (theta > 360.0)
		{
			theta-=360.0;
		}

		
		//printf("volume dimensions: %d %d %d",width,height,numSlices);
		
		// convert to radians
		rads = (theta*PI)/180;

		// get trig ratios
		sintheta = sin(rads);
		costheta = cos(rads);
		
		//check that vector parallel to the axis is a unit vector (if not then normalise)
		if ( fabs(mag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2])-1) > tol)
		{
			// nomalise
			axis[0] /= mag;
			axis[1] /= mag;
			axis[2] /= mag;
		}

		// to avoid excessive indexing
		nx = axis[0];
		ny = axis[1];
		nz = axis[2];

		//printf("(%f %f) (%f %f) (%f %f)\n",xlims[0],xlims[1],ylims[0],ylims[1],zlims[0],zlims[1]);
		//printf("%f %f %f",mx,my,mz);
		// calculate the distance between grid nodes (assumes regular grid)
		

		// for each grid point
		for (k=0; k<numSlices; k++)
		{
			for (j=0; j<height; j++)
			{
				for (i=0; i<width; i++)
				{
					// store the coordinates at this point
					x = xCoords[k][j][i];
					y = yCoords[k][j][i];
					z = zCoords[k][j][i];

					//apply general 3D rotation matrix, centred about the isocentre 
					xr = (1+(1-costheta)*(nx*nx-1))*(x-ox)+(-nz*sintheta+(1-costheta)*nx*ny)*(y-oy)+(ny*sintheta+(1-costheta)*nx*nz)*(z-oz)+ox;
					yr = (nz*sintheta+(1-costheta)*nx*ny)*(x-ox)+(1+(1-costheta)*(ny*ny-1))*(y-oy)+(-nx*sintheta+(1-costheta)*ny*nz)*(z-oz)+oy;
					zr = (-ny*sintheta+(1-costheta)*nx*nz)*(x-ox)+(nx*sintheta+(1-costheta)*ny*nz)*(y-oy)+(1+(1-costheta)*(nz*nz-1))*(z-oz)+oz;
					//TODO: add special cases for rotaion about coordinate axes
					
					xCoords[k][j][i] = xr;
					yCoords[k][j][i] = yr;
					zCoords[k][j][i] = zr;
						
					// updata bounds info
					if (xr < xlims[0])
					{
					    xlims[0] = xr;
					}
					if(xr > xlims[1])
 	    			{
					    xlims[1] = xr;
					}
					
					if(yr < ylims[0])
					{
					    ylims[0] = yr;
					}
					if(yr > ylims[1])
					{
					    ylims[1] = yr;
					}
					
					if (zr < zlims[0])
					{
					    zlims[0] = zr;
					}
					if(zr > zlims[1])
					{
					    zlims[1] = zr;
					}
				}
			}
		}
		
		
	};

    /**
     * Translate: shifts coordinate system such that x,y,z is the new origin
     */
	void DICOMReader::Translate(float x, float y, float z)
	{

		int i,j,k;

		for(k=0; k<numSlices;k++)
		{
		    for(j=0;j<height;j++)
		    {
		        for(i=0;i<width;i++)
		        {
		            xCoords[k][j][i] -=x;
		            yCoords[k][j][i] -=y;
		            zCoords[k][j][i] -=z;
		        }
		    }
		}
		
		xlims[0] -=x;
		xlims[1] -=x;
		ylims[0] -=y;
		ylims[1] -=y;
		zlims[0] -=z;
		zlims[1] -=z;
	};

    /**
     * WriteTestFile: outputs a file dump of data files
     *
     * NOTE: Only used for testing and debugging
     */
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

	/**Select_Region: removes all data outside the given bounds
	 *
	 * This resizes the data grid
	 */
	int DICOMReader::Select_Region(float xmin,float xmax,float ymin, float ymax, float zmin, float zmax)
	{

		int k,j,i; //for iterating through original volume
		int u,v,w; //for mapping to new volume
		
		// new grid dimension
		int newdimx = (int)round(((xmax-xmin)/(xlims[1]-xlims[0]))*width);
		int newdimy = (int)round(((ymax-ymin)/(ylims[1]-ylims[0]))*height);
		int newdimz = (int)round(((zmax-zmin)/(zlims[1]-zlims[0]))*numSlices);
		
		if (imageDataType==2)// short
		{
		    short*** new_data_s;
		    float*** newxCoords;
		    float*** newyCoords;
		    float*** newzCoords;
		    
		    new_data_s = new short**[newdimz];
		    newxCoords = new float**[newdimz];
		    newyCoords = new float**[newdimz];
		    newzCoords = new float**[newdimz];
		    for (k=0;k<newdimz;k++)
		    {
		        new_data_s[k] = new short*[newdimy];
		        newxCoords[k] = new float*[newdimy];
		        newyCoords[k] = new float*[newdimy];
		        newzCoords[k] = new float*[newdimy];
		        for(j=0;j<newdimy;j++)
		        {
		            new_data_s[k][j] = new short[newdimx];
		            newxCoords[k][j] = new float[newdimx];
		            newyCoords[k][j] = new float[newdimx];
		            newzCoords[k][j] = new float[newdimx];
		        }
		    }
		    
		    
		    
		    u=0;
	    	// this could be smater
		    for (k=0;k<numSlices && u<newdimz;k++)
		    {
		        v=0;
		        for(j=0;j<height && v<newdimy;j++)
		        {   
		            w=0;
		            for(i=0; i<width && w<newdimx;i++)
		            {
		                // if its in the bounds then add it to the new array
		                if(xCoords[k][j][i] < xmax && xCoords[k][j][i] > xmin && 
		                    yCoords[k][j][i] < ymax && yCoords[k][j][i] > ymin &&
		                    zCoords[k][j][i] < zmax && zCoords[k][j][i] > zmin)
		                {
		                    new_data_s[u][v][w] = data_s[k][j][i]; 
		                    newxCoords[u][v][w] = xCoords[k][j][i];
		                    newyCoords[u][v][w] = yCoords[k][j][i];
		                    newzCoords[u][v][w] = zCoords[k][j][i];
		                    w++;
		                }
		            }
		            if (w) // only increment if we found some stuff
		            {
		                v++;
		            }
		        }
		        if(v) // same as above
		        {
		            u++;
		        }
		    }
		    
		    // deallocate the old data
		    for (k=0;k<numSlices;k++)
		    {
		        for (j=0;j<height;j++)
		        {
		            delete xCoords[k][j];
		            delete yCoords[k][j];
		            delete zCoords[k][j];
		            delete data_s[k][j];
		        }
		        delete xCoords[k];
		        delete yCoords[k];
		        delete zCoords[k];
		        delete data_s[k];
		    }
		    delete xCoords;
		    delete yCoords;
		    delete zCoords;
		    delete data_s;
		    
		    // now point to the new data
		    data_s = (short***) new_data_s;
		    xCoords = (float ***) newxCoords;
		    yCoords = (float ***) newyCoords;
		    zCoords = (float ***) newzCoords;
		    
		    // update relevant data
		    numSlices = newdimz;
		    height = newdimy;
		    width = newdimx;
		    
		    xlims[0] = xmin;
		    xlims[1] = xmax;
		    ylims[0] = ymin;
		    ylims[1] = ymax;
		    zlims[0] = zmin;
		    zlims[1] = zmax;
		    
		    //NOTE: mx,my,mz are not updated as they represent the midpoint of the original volume
		}
		else if(imageDataType == 0)// float
		{
		    float*** new_data_f;
		    float*** newxCoords;
		    float*** newyCoords;
		    float*** newzCoords;
		    
		    new_data_f = new float**[newdimz];
		    newxCoords = new float**[newdimz];
		    newyCoords = new float**[newdimz];
		    newzCoords = new float**[newdimz];
		    for (k=0;k<newdimz;k++)
		    {
		        new_data_f[k] = new float*[newdimy];
		        newxCoords[k] = new float*[newdimy];
		        newyCoords[k] = new float*[newdimy];
		        newzCoords[k] = new float*[newdimy];
		        for(j=0;j<newdimy;j++)
		        {
		            new_data_f[k][j] = new float[newdimx];
		            newxCoords[k][j] = new float[newdimx];
		            newyCoords[k][j] = new float[newdimx];
		            newzCoords[k][j] = new float[newdimx];
		        }
		    }
		    
		    
		    
		    u=0;
	    	// this could be smater
		    for (k=0;k<numSlices && u<newdimz;k++)
		    {
		        v=0;
		        for(j=0;j<height && v<newdimy;j++)
		        {   
		            w=0;
		            for(i=0; i<width && w<newdimx;i++)
		            {
		                // if its in the bounds then add it to the new array
		                if(xCoords[k][j][i] < xmax && xCoords[k][j][i] > xmin && 
		                    yCoords[k][j][i] < ymax && yCoords[k][j][i] > ymin &&
		                    zCoords[k][j][i] < zmax && zCoords[k][j][i] > zmin)
		                {
		                    new_data_f[u][v][w] = data_f[k][j][i]; 
		                    newxCoords[u][v][w] = xCoords[k][j][i];
		                    newyCoords[u][v][w] = yCoords[k][j][i];
		                    newzCoords[u][v][w] = zCoords[k][j][i];
		                    w++;
		                }
		            }
		            if (w) // only increment if we found some stuff
		            {
		                v++;
		            }
		        }
		        if(v) // same as above
		        {
		            u++;
		        }
		    }
		    
		    // deallocate the old data
		    for (k=0;k<numSlices;k++)
		    {
		        for (j=0;j<height;j++)
		        {
		            delete xCoords[k][j];
		            delete yCoords[k][j];
		            delete zCoords[k][j];
		            delete data_f[k][j];
		        }
		        delete xCoords[k];
		        delete yCoords[k];
		        delete zCoords[k];
		        delete data_f[k];
		    }
		    delete xCoords;
		    delete yCoords;
		    delete zCoords;
		    delete data_f;
		    
		    // now point to the new data
		    data_f = (float***) new_data_f;
		    xCoords = (float ***) newxCoords;
		    yCoords = (float ***) newyCoords;
		    zCoords = (float ***) newzCoords;
		    
		    // update relevant data
		    numSlices = newdimz;
		    height = newdimy;
		    width = newdimx;
		    
		    xlims[0] = xmin;
		    xlims[1] = xmax;
		    ylims[0] = ymin;
		    ylims[1] = ymax;
		    zlims[0] = zmin;
		    zlims[1] = zmax;
		    
		    //NOTE: mx,my,mz are not updated as they represent the midpoint of the original volume
		}
		else if(imageDataType==1)//unsigned char
		{
		    unsigned char*** new_data_uc;
		    float*** newxCoords;
		    float*** newyCoords;
		    float*** newzCoords;
		    
		    new_data_uc = new unsigned char**[newdimz];
		    newxCoords = new float**[newdimz];
		    newyCoords = new float**[newdimz];
		    newzCoords = new float**[newdimz];
		    for (k=0;k<newdimz;k++)
		    {
		        new_data_uc[k] = new unsigned char*[newdimy];
		        newxCoords[k] = new float*[newdimy];
		        newyCoords[k] = new float*[newdimy];
		        newzCoords[k] = new float*[newdimy];
		        for(j=0;j<newdimy;j++)
		        {
		            new_data_uc[k][j] = new unsigned char[newdimx];
		            newxCoords[k][j] = new float[newdimx];
		            newyCoords[k][j] = new float[newdimx];
		            newzCoords[k][j] = new float[newdimx];
		        }
		    }
		    
		    
		    
		    u=0;
	    	// this could be smater
		    for (k=0;k<numSlices && u<newdimz;k++)
		    {
		        v=0;
		        for(j=0;j<height && v<newdimy;j++)
		        {   
		            w=0;
		            for(i=0; i<width && w<newdimx;i++)
		            {
		                // if its in the bounds then add it to the new array
		                if(xCoords[k][j][i] < xmax && xCoords[k][j][i] > xmin && 
		                    yCoords[k][j][i] < ymax && yCoords[k][j][i] > ymin &&
		                    zCoords[k][j][i] < zmax && zCoords[k][j][i] > zmin)
		                {
		                    new_data_uc[u][v][w] = data_uc[k][j][i]; 
		                    newxCoords[u][v][w] = xCoords[k][j][i];
		                    newyCoords[u][v][w] = yCoords[k][j][i];
		                    newzCoords[u][v][w] = zCoords[k][j][i];
		                    w++;
		                }
		            }
		            if (w) // only increment if we found some stuff
		            {
		                v++;
		            }
		        }
		        if(v) // same as above
		        {
		            u++;
		        }
		    }
		    
		    // deallocate the old data
		    for (k=0;k<numSlices;k++)
		    {
		        for (j=0;j<height;j++)
		        {
		            delete xCoords[k][j];
		            delete yCoords[k][j];
		            delete zCoords[k][j];
		            delete data_uc[k][j];
		        }
		        delete xCoords[k];
		        delete yCoords[k];
		        delete zCoords[k];
		        delete data_uc[k];
		    }
		    delete xCoords;
		    delete yCoords;
		    delete zCoords;
		    delete data_uc;
		    
		    // now point to the new data
		    data_uc = (unsigned char***) new_data_uc;
		    xCoords = (float ***) newxCoords;
		    yCoords = (float ***) newyCoords;
		    zCoords = (float ***) newzCoords;
		    
		    // update relevant data
		    numSlices = newdimz;
		    height = newdimy;
		    width = newdimx;
		    
		    xlims[0] = xmin;
		    xlims[1] = xmax;
		    ylims[0] = ymin;
		    ylims[1] = ymax;
		    zlims[0] = zmin;
		    zlims[1] = zmax;
		    
		    //NOTE: mx,my,mz are not updated as they represent the midpoint of the original volume
		}
		else if(imageDataType==3)//unsigned short
		{
		    unsigned short*** new_data_us;
		    float*** newxCoords;
		    float*** newyCoords;
		    float*** newzCoords;
		    
		    new_data_us = new unsigned short**[newdimz];
		    newxCoords = new float**[newdimz];
		    newyCoords = new float**[newdimz];
		    newzCoords = new float**[newdimz];
		    for (k=0;k<newdimz;k++)
		    {
		        new_data_us[k] = new unsigned short*[newdimy];
		        newxCoords[k] = new float*[newdimy];
		        newyCoords[k] = new float*[newdimy];
		        newzCoords[k] = new float*[newdimy];
		        for(j=0;j<newdimy;j++)
		        {
		            new_data_us[k][j] = new unsigned short[newdimx];
		            newxCoords[k][j] = new float[newdimx];
		            newyCoords[k][j] = new float[newdimx];
		            newzCoords[k][j] = new float[newdimx];
		        }
		    }
		    
		    
		    
		    u=0;
	    	// this could be smater
		    for (k=0;k<numSlices && u<newdimz;k++)
		    {
		        v=0;
		        for(j=0;j<height && v<newdimy;j++)
		        {   
		            w=0;
		            for(i=0; i<width && w<newdimx;i++)
		            {
		                // if its in the bounds then add it to the new array
		                if(xCoords[k][j][i] < xmax && xCoords[k][j][i] > xmin && 
		                    yCoords[k][j][i] < ymax && yCoords[k][j][i] > ymin &&
		                    zCoords[k][j][i] < zmax && zCoords[k][j][i] > zmin)
		                {
		                    new_data_us[u][v][w] = data_us[k][j][i]; 
		                    newxCoords[u][v][w] = xCoords[k][j][i];
		                    newyCoords[u][v][w] = yCoords[k][j][i];
		                    newzCoords[u][v][w] = zCoords[k][j][i];
		                    w++;
		                }
		            }
		            if (w) // only increment if we found some stuff
		            {
		                v++;
		            }
		        }
		        if(v) // same as above
		        {
		            u++;
		        }
		    }
		    
		    // deallocate the old data
		    for (k=0;k<numSlices;k++)
		    {
		        for (j=0;j<height;j++)
		        {
		            delete xCoords[k][j];
		            delete yCoords[k][j];
		            delete zCoords[k][j];
		            delete data_us[k][j];
		        }
		        delete xCoords[k];
		        delete yCoords[k];
		        delete zCoords[k];
		        delete data_us[k];
		    }
		    delete xCoords;
		    delete yCoords;
		    delete zCoords;
		    delete data_us;
		    
		    // now point to the new data
		    data_us = (unsigned short***) new_data_us;
		    xCoords = (float ***) newxCoords;
		    yCoords = (float ***) newyCoords;
		    zCoords = (float ***) newzCoords;
		    
		    // update relevant data
		    numSlices = newdimz;
		    height = newdimy;
		    width = newdimx;
		    
		    xlims[0] = xmin;
		    xlims[1] = xmax;
		    ylims[0] = ymin;
		    ylims[1] = ymax;
		    zlims[0] = zmin;
		    zlims[1] = zmax;
		    
		    //NOTE: mx,my,mz are not updated as they represent the midpoint of the original volume
		}

		return 0;
	};

	// simple... returns the midpoint (assumed Iso-centre unless specified by the user)
	void DICOMReader::GetMidPoint(float* x,float* y,float* z)
	{
		*x = mx;
		*y = my;
		*z = mz;
	};

} //namespace DICOM

#endif //DICOMREADER_H
