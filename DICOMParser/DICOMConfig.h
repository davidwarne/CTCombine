/*=========================================================================

  Program:   DICOMParser
  Module:    $RCSfile: DICOMConfig.h,v $
  Language:  C++
  Date:      $Date: 2003/10/23 19:13:34 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 2003 Matt Turek
  All rights reserved.
  See Copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __DICOM_CONFIG_H_
#define __DICOM_CONFIG_H_

//
// CMake Hook
//
//#include "DICOMCMakeConfig.h"

//
// BEGIN Toolkit (ITK,VTK, etc) specific
// 
#ifdef vtkDICOMParser_EXPORTS
  #define DICOM_EXPORT_SYMBOLS
#endif

#ifdef ITKDICOMParser_EXPORTS
  #define DICOM_EXPORT_SYMBOLS
#endif

//
// END toolkit (ITK, VTK, etc) specific
//

#define dicom_stl std
#define dicom_stream std

#include <iostream>
#include <fstream>
#include <iomanip>




#define DICOM_EXPORT 


#endif // __DICOM_CONFIG_H_
