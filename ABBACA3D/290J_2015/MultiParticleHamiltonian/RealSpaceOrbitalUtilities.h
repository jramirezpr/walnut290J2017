/*
#############################################################################
#
# Copyright 2015-16 Chris Anderson
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


#ifndef __RealSpaceOrbitalUtilities__
#define __RealSpaceOrbitalUtilities__

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
using namespace std;

#include "GridFunction3d.h"
#include "GridFunction3dUtility.h"
#include "DoubleArray1D.h"
#include "DoubleArray2D.h" 
#include "QdomainBasisFunctions.h"

#include "StringUtilities.h"

namespace UCLAQ
{

class RealSpaceOrbitalUtilities 
{

public: 	

RealSpaceOrbitalUtilities()
{
basisVectorsPtr = 0;
}
	
void inputDensityMatrix(string densityMatrixFileName,
long& basisCount, DoubleArray2D& densityMatrix)
{

	ifstream densityMatrixFile(densityMatrixFileName.c_str());

	if(!densityMatrixFile)
	{
	printf("##########   Error    ##########################\n"); 
	printf("Density Matrix File Not Found\n\n"); 
	printf("Input File Specified: \n\n%s\n\n",densityMatrixFileName.c_str());
	printf("################################################\n");
	exit(1);
	}
	//
	// 
	// Read lines of input file as strings to extract basis count.
	// Stop reading when "Density Values" is found as a leading token.
	//
	int readStop = 0;
	
	std::string inputLine;
	std::string delim(":");
	std::string token;
	std::vector<string> inputTokens;
	
	while((densityMatrixFile)&&(readStop == 0))
	{
	inputLine.clear();
	inputTokens.clear();
	getline(densityMatrixFile,inputLine);
	Tokenize(inputLine,inputTokens,delim);
	token = trim(inputTokens[0]);
	if(token == string("Density Values")) {readStop = 1;}
	else
	{
		if(inputTokens.size() > 1)
		{
		if(token == string("Orbital Count"))
		{
		token = trim(inputTokens[1]);
		istringstream iss (token.c_str(),istringstream::in);
		iss >> basisCount;
		}
		}
	}
	}
	//
	// Now read in the density matrix data 
	//
	densityMatrix.initialize(basisCount,basisCount);
	
	long i;
	long j;
    for(i = 0;  i < basisCount; i++)
    {
		for(j = 0; j < basisCount; j++)
		{
    	densityMatrixFile >> densityMatrix(i,j);
    	//printf("%20.15e ",densityMatrix(i,j));
    	}
    	printf("\n");
    }
	
	densityMatrixFile.close();
}
	
void initializeOrbitalBasis(string orbitalBasisFileName,long orbitalBasisSize = -1)
{
    vector<double>  basisIndexValues;
    QdomainBasisFunctions gbBasisFunctions;
    gbBasisFunctions.inputFromFile(basisIndexValues,internalBasisVectors,orbitalBasisFileName,orbitalBasisSize);
    basisVectorsPtr = &internalBasisVectors;
}

void createDensity(DoubleArray2D& densityMatrix, bool density3Dflag,
string densityMatrixFileName, bool logDensityFlag, bool binaryDensityFlag)
{
    vector <GridFunction3d>& basisVectors = *basisVectorsPtr;
	GridFunction3d gridFunTemp(basisVectors[0]);
	gridFunTemp.setToValue(0.0);

	GridFunction3d densityTemp(basisVectors[0]);
	densityTemp.setToValue(0.0);
	
	gridFunTemp  = basisVectors[0];
	gridFunTemp *= basisVectors[0];

	long i; long j;
	
	long basisCount = densityMatrix.getIndex1Size();
	
	for(i = 0; i < basisCount; i++)
	{
		for(j = 0; j < basisCount; j++)
		{
		gridFunTemp  = basisVectors[i];
		gridFunTemp *= basisVectors[j];
		gridFunTemp *= densityMatrix(i,j);
		densityTemp += gridFunTemp;
	}}
	
	printf("\n#######################################################\n\n");
    printf("Density integral            :  %-15.5e\n", densityTemp.integral());
  

    const char* outputFile = "orbitalDensity.vtk";
	if(density3Dflag == 1)
	{
	gUtility.outputDataToVTKfile(outputFile,"OrbitalDensity",densityTemp);
    printf("Density written to          :  %s\n",outputFile);
	}

    FILE* binaryDataFile;
    string binaryDataFileName = "orbitalDensity.dat";

	if(binaryDensityFlag == true)
	{
    if( (binaryDataFile = fopen(binaryDataFileName.c_str(), "wb" )) == NULL )
    {
          printf( "The file %s could not be  opened\n",binaryDataFileName.c_str());
          exit(1);
    }

    gUtility.outputToBinaryDataFile(densityTemp,binaryDataFileName.c_str());
    fclose(binaryDataFile);
    printf("Binary density written to   :  %s\n",binaryDataFileName.c_str());
	}


    outputFile = "orbitalDensityLog.vtk";
    double        *dataPtr;
    long    totalDataCount;

	if(logDensityFlag == true)
	{
    totalDataCount = densityTemp.Values.getIndex1Size()*densityTemp.Values.getIndex2Size()*densityTemp.Values.getIndex3Size();
    dataPtr        = densityTemp.Values.getDataPointer();
    for(long j = 0; j < totalDataCount; j++)
    {
    	if(dataPtr[j] < 1.0e-06) dataPtr[j] = 1.0e-06;
    	dataPtr[j] = log10(dataPtr[j]);
    }
	gUtility.outputDataToVTKfile(outputFile,"LogOrbitalDensity",densityTemp);
    printf("Log Density written to      :  %s\n",outputFile);
	}

}

    vector <GridFunction3d>   internalBasisVectors;
    vector <GridFunction3d>*       basisVectorsPtr;
    
	GridFunction3d globalGridFunTemp;
	GridFunction3d  localGridFunTemp;
	
	GridFunction3dUtility   gUtility;

};
} // UCLAQ namespace

#endif  // __RealSpaceOrbitalUtilities__
