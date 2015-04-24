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
/*file: EGSPhantTodx.c
 * 
 * Author: David Warne
 * Date: 03/02/2009
 *
 * Converts .egsphant format to openDx .dx format
 *
 * Syntax: EGSPhant EGSPhantFileName.egsphant dxOutFileName.dx
 */

#include<stdio.h>
#include<stdlib.h>
#include"EGSPhant.h"

using namespace EGS;



int main(int argc, char** argv){

	EGSPhant EGSPhantom(argv[1]);
	FILE* DXFILE;
	float ox,oy,oz,xstep,ystep,zstep;

	int i,j,k;

	//first check and read the input
	if(argc == 3){

		//if(!EGSPhantom(argv[1])){
			//printf("ERROR: A problem occured when trying initialise the ESG object [error code: (1/sqrt(pi))*int_0^inf e^(-v^2) dv]\n");
			//printf("Please call the IT helpdesk (or a mathematician to help with the error code)\n");
		//}
		//else{

			EGSPhantom.Read();
			if(!(DXFILE = fopen(argv[2],"w"))){
				printf("ERROR: Could not read file %s",argv[2]);
				return 1;
			}

			//first write the grid positions data
			fprintf(DXFILE,"object \"nodes\" class gridpositions\n");
			fprintf(DXFILE,"counts %d %d %d\n",EGSPhantom.xSize,EGSPhantom.ySize,EGSPhantom.zSize);
			ox = (EGSPhantom.xBoundaries[EGSPhantom.xSize] + EGSPhantom.xBoundaries[0])*0.5;
			oy = (EGSPhantom.yBoundaries[EGSPhantom.ySize] + EGSPhantom.yBoundaries[0])*0.5;
			oz =(EGSPhantom.zBoundaries[EGSPhantom.zSize] + EGSPhantom.zBoundaries[0])*0.5;
			fprintf(DXFILE,"origin %f %f %f\n",ox,oy,oz);
			xstep = (EGSPhantom.xBoundaries[EGSPhantom.xSize] - EGSPhantom.xBoundaries[0])/EGSPhantom.xSize;
			ystep = (EGSPhantom.xBoundaries[EGSPhantom.ySize] - EGSPhantom.yBoundaries[0])/EGSPhantom.ySize;
			zstep = (EGSPhantom.xBoundaries[EGSPhantom.zSize] - EGSPhantom.zBoundaries[0])/EGSPhantom.zSize;
			fprintf(DXFILE,"delta %f 0 0\n",xstep);
			fprintf(DXFILE,"delta 0 %f 0\n",ystep);
			fprintf(DXFILE,"delta 0 0 %f\n",zstep);

			//write the grid connection data
			fprintf(DXFILE,"object \"cubes\" class gridconnections\n");
			fprintf(DXFILE,"counts %d %d %d\n",EGSPhantom.xSize,EGSPhantom.ySize,EGSPhantom.zSize);
			fprintf(DXFILE,"attribute \"ref\" string \"positions\"\n");
			fprintf(DXFILE,"attribute \"element type\" string \"cubes\"\n");

			//write the meduim data
			fprintf(DXFILE,"object\"meduims\" class array type byte rank 0 items %d\n",EGSPhantom.xSize*EGSPhantom.ySize*EGSPhantom.zSize);

			for(k=0;k<EGSPhantom.zSize;k++){
				for(j=0;j<EGSPhantom.ySize;j++){
					for(i=0;i<EGSPhantom.xSize;i++){

						fprintf(DXFILE,"%c\n",EGSPhantom.voxelMedium[k][j][i]);
					}
				}
			}
			fprintf(DXFILE,"attribut \"dep\" string \"positions\"\n");

			//write the meduim data
			fprintf(DXFILE,"object\"densities\" class array type byte rank 0 items %d\n",EGSPhantom.xSize*EGSPhantom.ySize*EGSPhantom.zSize);

			for(k=0;k<EGSPhantom.zSize;k++){
				for(j=0;j<EGSPhantom.ySize;j++){
					for(i=0;i<EGSPhantom.xSize;i++){

						fprintf(DXFILE,"%f\n",EGSPhantom.voxelDensity[k][j][i]);
					}
				}
			}
			fprintf(DXFILE,"attribut \"dep\" string \"positions\"\n");

			//now just the footer
			fprintf(DXFILE,"object \"%s\" class field\n");
			fprintf(DXFILE,"component \"positions\" value \"nodes\"\n");
			fprintf(DXFILE,"component \"connections\" value \"cubes\"\n");
			fprintf(DXFILE,"component \"data\" value \"mediums\"\n");
			fprintf(DXFILE,"component \"data\" value \"densities\"\n");

			fprintf(DXFILE,"end\n");

			fclose(DXFILE);

			return 0;
		//}
		
	}
}
