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


#ifndef DICOMREADER_H
#define DICOMREADER_H
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
}//namespace DICOM

#endif //DICOMREADER_H
