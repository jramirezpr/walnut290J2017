/*
 * QdomainBasisFunctions.h
 *
 *  Created on: May 29, 2014
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright 2015 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

#ifndef _QdomainBasisFunctions_
#define _QdomainBasisFunctions_

#include <iostream>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdexcept>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

namespace UCLAQ
{

class QdomainBasisFunctions
{
	public:

	QdomainBasisFunctions()
	{
	verboseFlag = false;
	}

	void setVerboseFlag(bool flagValue = true)
	{
	verboseFlag = flagValue;
	}

	void clearVerboseFlag()
	{
	verboseFlag = false;
	}

	void outputToFile(vector<double>& basisIndexValues, vector<UCLAQ::GridFunction3d>& basisVectors, string basisFileName)
	{
    string outputPath;
	if(!outputPath.empty())
	{
	    outputPath.append("/");
	    outputPath.append(basisFileName);
	    basisFileName = outputPath;
	}

	double eigVal;
	FILE* binaryDataFile;

	if( (binaryDataFile = fopen(basisFileName.c_str(), "wb" )) == NULL )
	{
	string msg ="\nXXX QdomainBasisFunctions::outputToFile(...) XXX\n";
	msg.append("Basis set file  \"" + basisFileName + "\" could not be opened. \n");
    throw std::runtime_error(msg);
	}

    long basisCount = basisVectors.size();

	double doubleBasisCount = basisCount + 1.0e-07;
	fwrite(&doubleBasisCount, sizeof(double), 1, binaryDataFile);

	if(verboseFlag)
    {
	printf("Writing %ld basis vectors to %s \n",basisCount,basisFileName.c_str());
	printf("Basis index values : \n");
	}

	for(long k = 0; k < basisCount; k++ )
	{
	    eigVal = basisIndexValues[k];
	    fwrite(&eigVal,  sizeof(double), 1, binaryDataFile);
	    if(verboseFlag)
	    {
	    printf("%ld : %10.5e \n",k,basisIndexValues[k]);
	    }
	}
	//
	// First basis function  output contains the structure information
	//

	gUtility3D.outputToBinaryDataFile(basisVectors[0], binaryDataFile);

	//
	// Remaining basis function output just contains the values
	//
	for(long k = 1; k < basisCount; k++ )
	{
	gUtility3D.appendValuesToBinaryDataFile(basisVectors[k], binaryDataFile);
	}

	fclose(binaryDataFile);
	}


    int getFileBasisCount(const string& basisFileName)
    {
    FILE* binaryDataFile;

	if( (binaryDataFile = fopen(basisFileName.c_str(), "rb" )) == NULL )
	{
	string msg ="\nXXX QdomainBasisFunctions::inputFromFile(...) XXX\n";
	msg.append("Basis set file  \"" + basisFileName + "\" not found. \n");
    throw std::runtime_error(msg);
	}

	double  doubleBasisCount;
	long    basisCount;

    //
    // Read in sizes and all of the basis index values in the input file
    //

	fread(&doubleBasisCount,  sizeof(double), 1, binaryDataFile);
	basisCount  = (int)(doubleBasisCount);
    return basisCount;
    }

	void inputFromFile(vector<double>& basisIndexValues, vector<UCLAQ::GridFunction3d>& basisVectors,
	const string& basisFileName, long basisInputSize = -1)
	{
	FILE* binaryDataFile;

	if( (binaryDataFile = fopen(basisFileName.c_str(), "rb" )) == NULL )
	{
	string msg ="\nXXX QdomainBasisFunctions::inputFromFile(...) XXX\n";
	msg.append("Basis set file  \"" + basisFileName + "\" not found. \n");
    throw std::runtime_error(msg);
	}

	double  doubleBasisCount;
	long    basisCount;

    //
    // Read in sizes and all of the basis index values in the input file
    //

	fread(&doubleBasisCount,  sizeof(double), 1, binaryDataFile);
	basisCount  = (int)(doubleBasisCount);

	basisIndexValues.clear();
	basisIndexValues.resize(basisCount,0.0);

	for(long k = 0; k < basisCount; k++ )
	{
	  fread(&basisIndexValues[k],  sizeof(double), 1, binaryDataFile);
	}

    //
    // Reset size if requesting fewer than total number of basis elements
    //

    if((basisInputSize > 0)&&(basisInputSize < basisCount))
	{
	basisCount = basisInputSize;
	}
	basisIndexValues.resize(basisCount);

	if(verboseFlag)
    {
	    printf("Reading in %ld basis vectors from %s \n",basisCount,basisFileName.c_str());
	    printf("Basis index values : \n");
	    for(long k = 0; k < basisCount; k++ )
	    {
	    printf("%ld : %10.5e \n",k,basisIndexValues[k]);
	    }
	}

    basisVectors.clear();
	basisVectors.resize(basisCount);
	gUtility3D.inputFromBinaryDataFile(basisVectors[0], binaryDataFile);

	//
	// Initialize the remaining basis functions with the first, and then
	// read in their values
	//
	for(long k = 1; k < basisCount; k++ )
	{
		   cout << k << endl;
	       basisVectors[k].initialize(basisVectors[0]);
	       gUtility3D.inputValuesFromBinaryDataFile(basisVectors[k],binaryDataFile);
	}
	fclose(binaryDataFile);
	}

	void outputSquaredToVTK(vector<UCLAQ::GridFunction3d>& basisVectors, string basisFileName, double xScalingFactor = -1.0)
	{
	string baseName;
	int lastindex = basisFileName.find_last_of(".");
	if(lastindex != (int)string::npos)
	{
    baseName.assign(basisFileName.substr(0, lastindex));
    }
    else
    {
    baseName.assign(basisFileName);
    }
    ostringstream basisfile;
    ostringstream basisName;

    long basisCount = basisVectors.size();
    GridFunction3d eigVecTemp(basisVectors[0]);

	if(verboseFlag)
    {
	printf("Writing %ld basis vectors in VTK format : \n",basisCount);
	}
	for(long k = 0; k < basisCount; k++ )
	{
       basisName.str("");
       basisfile.str("");
       if(k <= 9)
       {
       basisName << baseName<< "_2_0"  << k << ends;
       basisfile << baseName <<"_2_0"  << k << ".vtk" << ends;
       }
       else
       {
       basisName << baseName<< "_2_"  << k << ends;
       basisfile << baseName <<"_2_"  << k << ".vtk" << ends;
       }

       eigVecTemp   =  basisVectors[k];
       eigVecTemp.squareValues();

   	   gUtility3D.outputDataToVTKfile(eigVecTemp,basisfile.str().c_str(),basisName.str().c_str(),xScalingFactor);
   	   if(verboseFlag)
   	   {
       printf("%s \n",(basisfile.str()).c_str());
       }
	}
	}

	void outputToVTK(vector<UCLAQ::GridFunction3d>& basisVectors, string basisFileName,double xScalingFactor=-1)
	{
	string baseName;
	int lastindex = basisFileName.find_last_of(".");
	if(lastindex != (int)string::npos)
	{
    baseName.assign(basisFileName.substr(0, lastindex));
    }
    else
    {
    baseName.assign(basisFileName);
    }
    ostringstream basisfile;
    ostringstream basisName;
    long basisCount = basisVectors.size();

    if(verboseFlag)
    {
	printf("Writing %ld basis vectors in VTK format : \n",basisCount);
	}

	for(long k = 0; k < basisCount; k++ )
	{
       basisName.str("");
       basisfile.str("");
       if(k <= 9)
       {
       basisName << baseName<< "_0"  << k << ends;
       basisfile << baseName <<"_0"  << k << ".vtk" << ends;
       }
       else
       {
       basisName << baseName<< "_"  << k << ends;
       basisfile << baseName <<"_"  << k << ".vtk" << ends;
       }

   	   gUtility3D.outputDataToVTKfile(basisVectors[k],basisfile.str().c_str(),basisName.str().c_str(),xScalingFactor);
   	   if(verboseFlag)
   	   {
       printf("%s \n",(basisfile.str()).c_str());
       }
	}
	}

    bool                        verboseFlag;
	UCLAQ::GridFunction3dUtility gUtility3D;
};


} // UCLAQ  namesapce

#endif /* _QdomainBasisFunctions_ */
