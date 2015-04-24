/* CTCombine: A Tool for egsphant creation and manipulation.
 * Copyright (C) 2015 David Warne and Mark Dwyer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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



#ifndef EGSPHANT_H
#define EGSPHANT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

template <class T> const T& min ( const T& a, const T& b ) {
  return (a<b)?a:b;    
}

template <class T> const T& max ( const T& a, const T& b ) {
  return (a>b)?a:b;     
}

#define PI  3.14159265358979323846



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

	

} //namespace EGS

#endif //MYMATH_MATRIX_H
