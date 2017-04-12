/*
#############################################################################
#
# Copyright 2015-2016 Chris Anderson
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
//
//  Removed nuclear interpolation for evaluation of nuclear potential contribution
//  and code for incremental construction
//
//
//  Thurs. Dec. 17, 2015
//
//  author: Chris Anderson
//  (C) UCLA

#ifdef _USE_ASSERTS_
#include <cassert>
#endif

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <sys/param.h>
#include <unistd.h>
using namespace std;


#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

#include "MultiParticleHamiltonian/RealSpatialOrbitalIntegrals.h"

#include "Timing/ClockIt.h"

#ifndef _QdomainOrbitalIntegrals_
#define _QdomainOrbitalIntegrals_

#define  _DEFAULT_INV_LAPLACIAN_EXTENSION_FACTOR_ 2.0
#define  _DEFAULT_OUTPUT_FILENAME_      "orbitalIntegrals.dat"
#define  _DEFAULT_laplaceOpCoefficient 1.0/(4.0*3.14159265358979323846)
//
// Member function required of InverseLaplaceOp for instantiation
//
// InverseLaplaceOp()
// void initialize();
// void initialize(double laplaceOpCoefficient, long xPanel, double xMin, double xMax,
//                                              long yPanel, double yMin, double yMax,
//                                              long zPanel, double zMin, double zMax, double invLapExtFactor);
// void applyInverseOp(const UCLAQ::GridFunction3d& V);
//
//

namespace UCLAQ
{

template < class InverseLaplaceOp, class SingleParticleOp  > class QdomainOrbitalIntegrals
{

public :

QdomainOrbitalIntegrals()
{
    initialize();
}

void initialize()
{
    doubleOrbitalIntegralData.clear();
    doubleOrbitalIntegralCount = 0;

    singleOrbitalIntegralData.clear();
    singleOrbitalIntegralCount = 0;

    singleOrbitalIntegralIndices.clear();
    doubleOrbitalIntegralIndices.clear();

    singleOrbitalTime = 0.0;
	doubleOrbitalTime = 0.0;
	startupTime       = 0.0;
	totalTime         = 0.0;
    timingConstantsFlag    = false;
    binaryOutputFlag       = false;
    incrementalFlag        = false;

	xMin   = 0.0; xMax   = 0.0;  xPanel = 0;
	yMin   = 0.0; yMax   = 0.0;  yPanel = 0;
	zMin   = 0.0; zMax   = 0.0;  zPanel = 0;

	basisVectorsPtr       = NULL;

    invLaplaceOp3d.initialize();
	invLapExtFactor       = _DEFAULT_INV_LAPLACIAN_EXTENSION_FACTOR_;
	outputFileName        = _DEFAULT_OUTPUT_FILENAME_;
	laplaceOpCoefficient  = _DEFAULT_laplaceOpCoefficient;

	verboseFlag     = false;
}


~QdomainOrbitalIntegrals()
{}

void setVerbose(bool flagValue = true)
{verboseFlag = flagValue;}

void clearVerbose()
{verboseFlag = false;}

//
//  Coefficient of Laplace operator so that inverse Laplacian
//  with infinite boundary conditions yields electron-electron
//  interaction potential.
//
//  The default value = 1/(4*pi) is the coefficient associated with
//  the Schroedinger operator in atomic units.
//
void setLaplaceOpCoefficient(double value)
{
	laplaceOpCoefficient = value;
}


void setInvLaplaceExtFactor(double invLapExtFactor)
{
this->invLapExtFactor = invLapExtFactor;
}

//
// The domain characteristics are captured from either the external potential or
// one of the input orbital basis vectors. Consistency of the specification
// of the potential and orbital basis vectors is checked.
//

void attachBasisVectors(vector<UCLAQ::GridFunction3d>& externalBasisVectors)
{
	basisVectorsPtr = &externalBasisVectors;
	vector<UCLAQ::GridFunction3d>& basisVectors = *basisVectorsPtr;

    try{captureDomainSpecificationAndInitialize(basisVectors[0]);} catch (const runtime_error& e){throw e;}

    //
    // Create default basis index values
    //
    basisIndexValues.resize(basisVectors.size());
    for(long k = 0; k < basisIndexValues.size(); k++)
    {
    basisIndexValues[k] = (double)k;
    }
}

void attachSingleParticleOperator(SingleParticleOp* sOperatorPtr)
{
	singleParticleSop = sOperatorPtr;
}


void createOrbitalIntegrals(RealSpatialOrbitalIntegrals& realSpatialOrbitalIntegrals, long orbitalCount = -1)
{

    if(orbitalCount <= 0)
    {
    orbFunBasisCount =  basisVectorsPtr->size();
    }
    else if (orbFunBasisCount > basisVectorsPtr->size())
    {
    orbFunBasisCount =  basisVectorsPtr->size();
    }

    globalClockIt.start();
    clockIt.start();
    orbitalProgressValue = 0;

    createSingleOrbitalIntegrals();

    clockIt.stop();
    singleOrbitalTime = clockIt.getMilliSecElapsedTime();

    clockIt.start();
    orbitalProgressValue = 0;
    createDoubleOrbitalIntegrals();


    clockIt.stop();
    doubleOrbitalTime = clockIt.getMilliSecElapsedTime();

    // Pack orbital integrals into a RealSpatialOrbitalIntegrals instance

    realSpatialOrbitalIntegrals.initialize(orbFunBasisCount,&basisIndexValues[0],
    singleOrbitalIntegralCount, &singleOrbitalIntegralData[0], doubleOrbitalIntegralCount, &doubleOrbitalIntegralData[0]);

    globalClockIt.stop();

    if(verboseFlag)
    {
    outputTimingData();
    }
}



void outputOrbitalIntegrals(RealSpatialOrbitalIntegrals& realSpatialOrbitalIntegrals, string& orbitalIntegralsFile)
{
	outputFileName.assign(orbitalIntegralsFile);
    if((not binaryOutputFlag)&&(not timingConstantsFlag))
    {
        realSpatialOrbitalIntegrals.outputOrbitalIntegralData(outputFileName);
    }
    else if(not timingConstantsFlag)
    {
       outputOrbitalIntegralDataBinary();
    }
}


void captureDomainSpecificationAndInitialize(UCLAQ::GridFunction3d& gF)
{
	xMin = gF.getXmin();
	xMax = gF.getXmax();
	yMin = gF.getYmin();
	yMax = gF.getYmax();
	zMin = gF.getZmin();
	zMax = gF.getZmax();
	xPanel = gF.getXpanelCount();
	yPanel = gF.getYpanelCount();
	zPanel = gF.getZpanelCount();

	//  Initialize required grid functions

	oFunI.initialize(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);
	oFunJ.initialize(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);

	oFunIxJ.initialize(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);
	oFunMxN.initialize(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);

	// Initialize infinite domain Poisson solver

    invLaplaceOp3d.initialize(laplaceOpCoefficient, xPanel, xMin, xMax,yPanel, yMin, yMax,zPanel, zMin,  zMax, invLapExtFactor);
}


void createSingleOrbitalIntegrals()
{
    long i; long j;

    vector<UCLAQ::GridFunction3d>& basisVectors = *basisVectorsPtr;

    if(not incrementalFlag)
    {
    singleOrbitalIntegralCount     = (orbFunBasisCount)*(orbFunBasisCount+1)/2;
    singleOrbitalIntegralData.resize(singleOrbitalIntegralCount);
    }
    else
    {
    singleOrbitalIntegralData.resize(singleIntegralCountIncremental);
    singleOrbitalIntegralIndices.resize(singleIntegralCountIncremental);
    }


    long dataIndex  = 0;
    long ijIndex    = 0;
    long ijCount    = 0;

    //
    // When constructing integrals start indices counting at 1
    // and offset by -1 when accessing basis function array.
    //
    if(verboseFlag)
    {
    	if(incrementalFlag != 1)
    	{
    		printf("Computing %12ld one-electron integrals \n",singleOrbitalIntegralCount);
    	}
    	else
    	{
    		printf("Computing %12ld one-electron integrals \n",singleIntegralCountIncremental);
    	}
    }

    for(i = 1; i <= orbFunBasisCount; i++)
    {
    oFunI    = basisVectors[i-1];
    singleParticleSop->apply(oFunI);

    for(j = 1; j <= i; j++)
    {
    ijIndex = ((i-1)*i)/2 + j;
    oFunJ   = basisVectors[j-1];

    if(incrementalFlag != 1)
    {
    	singleOrbitalIntegralData[dataIndex]= oFunJ.scaledDot(oFunI);
    }
    else if((ijIndex >= ijIndexStart)&&(ijIndex <= ijIndexEnd))
    {
    	singleOrbitalIntegralData[ijCount]    = oFunJ.scaledDot(oFunI);
    	singleOrbitalIntegralIndices[ijCount] = dataIndex;
        ijCount++;
    }
    dataIndex++;
    if(verboseFlag)
    {
    if(not incrementalFlag)
    {
    	progressMessage(dataIndex,singleOrbitalIntegralCount);
    }
    else
    {
    	progressMessage(ijCount,singleIntegralCountIncremental);
    }}

    }}
}

//
// This method computes the single orbital integral values and then overwrites the
// corresponding values in the input realSpatialOrbitalIntegrals.
//
// The modified potential that necessitates an invocation to this member function
// is specified in the single particle Schroedinger operator which is referenced
// by this class. If the Schroedinger operator instance used to create the updated single orbital integrals
// differs from that used for the original invocation of createOrbitalIntegrals, then the new
// operator must be attached to this class using a call to the member function
//
// void attachSingleParticleOperator(QdomainSingleParticleSopBase* sOperatorPtr)
//
// before calling this update member function.
//
// !!! This method assumes that there has been a previous invocation of createOrbitalIntegrals
//      Error checking is NOT done to validate this assumption.
//
//  Assumptions:
//
// (*) The reference to the single particle Schroedinger operator has been set to a valid instance.
// (*) References to the orbital basis has been set to a valid instance.
// (*) The number of orbitals has been specified.
// (*) The nuclear interpolant has been initialized (if needed)
//
//

void updateSingleOrbitalIntegrals(RealSpatialOrbitalIntegrals& realSpatialOrbitalIntegrals)
{
    string errMessage;

    if(singleOrbitalIntegralCount != realSpatialOrbitalIntegrals.singleOrbitalIntegralCount)
    {
    errMessage  = "\nXXX QdomainOrbitalIntegrals XXX\n";
	errMessage += "Error in member function updateSingleOrbitalIntegrals(...)\n\n";
	errMessage += "Allocation for the updated single orbital integral data in \n";
	errMessage += "input realSpatialOrbitalIntegrals instance is \n";
	errMessage += "insufficient for updated values.\n";
	throw std::runtime_error(errMessage);
	}

    globalClockIt.start();
    clockIt.start();
    orbitalProgressValue = 0;

    createSingleOrbitalIntegrals();

    clockIt.stop();
    singleOrbitalTime = clockIt.getMilliSecElapsedTime();

    // Re-pack single orbital integral data into realSpatialOrbitalIntegrals;

    for(long i = 0; i < singleOrbitalIntegralCount; i++)
    {realSpatialOrbitalIntegrals.singleOrbitalIntegralData[i] = singleOrbitalIntegralData[i];}

    if(verboseFlag)
    {
    outputTimingData();
    }
}

void createDoubleOrbitalIntegrals()
{
	vector<UCLAQ::GridFunction3d>& basisVectors = *basisVectorsPtr;

    long i; long j;
    long m; long n;

    long p                   = orbFunBasisCount;
    long ijIndex    = 0;
    long dCount     = 0;


    doubleOrbitalIntegralCount = (p*(p*p*p + 2*p*p + 3*p +2))/8;

    if(not incrementalFlag)
    {
    if(verboseFlag)
    {
    printf("Computing %12ld two-electron integrals \n",doubleOrbitalIntegralCount);
    }
    doubleOrbitalIntegralData.resize(doubleOrbitalIntegralCount);
    for(i = 0; i < doubleOrbitalIntegralCount; i++)
    {doubleOrbitalIntegralData[i] = 0.0;}
    }
    else
    {
    if(verboseFlag)
    {
    printf("Computing %12ld two-electron integrals \n",doubleIntegralCountIncremental);
    }
    doubleOrbitalIntegralData.resize(doubleIntegralCountIncremental);
    doubleOrbitalIntegralIndices.resize(doubleIntegralCountIncremental);
    for(i = 0; i < doubleIntegralCountIncremental; i++)
    {doubleOrbitalIntegralData[i] = 0.0;}
    }


    long dataIndex  = 0;

    //
    // When constructing integrals start indices counting at 1
    // and offset by -1 when accessing basis function array.
    //

    for(i = 1; i <= orbFunBasisCount; i++)
    {
    for(j = 1; j <= i; j++)
    {
    ijIndex = ((i-1)*i)/2 + j;

    if(not incrementalFlag)
    {
    oFunIxJ  = basisVectors[i-1];
    oFunIxJ *= basisVectors[j-1];
    //
    // Evaluate the convolution of the Green's function with oFunI x oFunJ
    //
    invLaplaceOp3d.applyInverseOp(oFunIxJ);
    }
    else if((ijIndex >= ijIndexStart)&&(ijIndex <= ijIndexEnd))
    {
    oFunIxJ  = basisVectors[i-1];
    oFunIxJ *= basisVectors[j-1];
    //
    // Evaluate the convolution of the Green's function with oFunI x oFunJ
    //
    invLaplaceOp3d.applyInverseOp(oFunIxJ);
    }

    // m < i

    for(m = 1; m  < i;  m++)
    {
        for(n = 1; n <= m;  n++)
        {

        if(incrementalFlag != 1)
        {
        oFunMxN  = basisVectors[m-1];
        oFunMxN *= basisVectors[n-1];

        // - sign included because the Green's function is -Delta^(-1)

        doubleOrbitalIntegralData[dataIndex] = -oFunMxN.scaledDot(oFunIxJ);
        }
        else if((ijIndex >= ijIndexStart)&&(ijIndex <= ijIndexEnd))
        {

        oFunMxN  = basisVectors[m-1];
        oFunMxN *= basisVectors[n-1];

        // - sign included because the Green's function is -Delta^(-1)

         doubleOrbitalIntegralData[dCount]   = -oFunMxN.scaledDot(oFunIxJ);
         doubleOrbitalIntegralIndices[dCount] = dataIndex;
    	 dCount++;
        }

        dataIndex++;
        if(verboseFlag)
        {
        	if(incrementalFlag != 1)
        	{progressMessage(dataIndex,doubleOrbitalIntegralCount);}
        	else
        	{progressMessage(dCount,doubleIntegralCountIncremental);}
        }

        }
    }

    m = i;
    for(n = 1; n <= j;  n++)
    {
    	if(incrementalFlag != 1)
        {
        oFunMxN  = basisVectors[m-1];
        oFunMxN *= basisVectors[n-1];

        // - sign included because the Green's function is -Delta^(-1)

        doubleOrbitalIntegralData[dataIndex] = -oFunMxN.scaledDot(oFunIxJ);
        }
        else if((ijIndex >= ijIndexStart)&&(ijIndex <= ijIndexEnd))
        {
        oFunMxN  = basisVectors[m-1];
        oFunMxN *= basisVectors[n-1];

        // - sign included because the Green's function is -Delta^(-1)

        doubleOrbitalIntegralData[dCount]    = -oFunMxN.scaledDot(oFunIxJ);
        doubleOrbitalIntegralIndices[dCount] = dataIndex;
        dCount++;
        }

        dataIndex++;
        if(verboseFlag)
        {
        if(incrementalFlag != 1)
        {progressMessage(dataIndex,doubleOrbitalIntegralCount);}
        else
        {progressMessage(dCount,doubleIntegralCountIncremental);}
        }
    }

    }}
}

void outputOrbitalIntegralDataBinary()
{
    ostringstream s;
    FILE *binaryDataFile;

    if(incrementalFlag != 1)
    {
    if( (binaryDataFile = fopen(outputFileName.c_str(), "wb" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",outputFileName.c_str());
      exit(1);
    }

    int  orbFunBasisCountInt           = orbFunBasisCount;
    int  singleOrbitalIntegralCountInt = singleOrbitalIntegralCount;
    int  doubleOrbitalIntegralCountInt = doubleOrbitalIntegralCount;

    fwrite(&orbFunBasisCountInt,             sizeof(int), 1, binaryDataFile);
    fwrite(&singleOrbitalIntegralCountInt,   sizeof(int), 1, binaryDataFile);
    fwrite(&doubleOrbitalIntegralCountInt,   sizeof(int), 1, binaryDataFile);

    //fwrite(&orbFunBasisCount,             sizeof(long), 1, binaryDataFile);
    //fwrite(&singleOrbitalIntegralCount,   sizeof(long), 1, binaryDataFile);
    //fwrite(&doubleOrbitalIntegralCount,   sizeof(long), 1, binaryDataFile);


    fwrite(&basisIndexValues[0],sizeof(double),orbFunBasisCount,binaryDataFile);
    fwrite(&singleOrbitalIntegralData[0],     sizeof(double),singleOrbitalIntegralCount,binaryDataFile);
    fwrite(&doubleOrbitalIntegralData[0],     sizeof(double),doubleOrbitalIntegralCount,binaryDataFile);


    fclose(binaryDataFile);


    printf("\n#######################################################\n\n");
    printf("Orbital Integral Data File (binary format): %s   \n",outputFileName.c_str());
    printf("\n#######################################################\n");
    return;
    }

    //
    // Create file names for incremental results
    //
    string binFilePrefix(outputFileName);
    string singleIndicesName   = binFilePrefix + "SingleIndices.bin";
    string doubleIndicesName   = binFilePrefix + "DoubleIndices.bin";
    string singleIntegralsName = binFilePrefix + "SingleIntegrals.bin";
    string doubleIntegralsName = binFilePrefix + "DoubleIntegrals.bin";
    //
    // write out the number of indices output and then the indices
    //
    if( (binaryDataFile = fopen(singleIndicesName.c_str(), "wb" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",singleIndicesName.c_str());
      exit(1);
    }
    fwrite(&singleIntegralCountIncremental,sizeof(int), 1, binaryDataFile);
    fwrite(&singleOrbitalIntegralIndices[0], sizeof(int),singleIntegralCountIncremental, binaryDataFile);
    fclose(binaryDataFile);

    //
    // write out the number of integrals output and then the integrals
    //
    if( (binaryDataFile = fopen(singleIntegralsName.c_str(), "wb" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",singleIntegralsName.c_str());
      exit(1);
    }
    fwrite(&singleIntegralCountIncremental,sizeof(int), 1, binaryDataFile);
    fwrite(&singleOrbitalIntegralData[0], sizeof(double),singleIntegralCountIncremental, binaryDataFile);
    fclose(binaryDataFile);

    //
    // write out the number of indices output and then the indices
    //
    if( (binaryDataFile = fopen(doubleIndicesName.c_str(), "wb" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",doubleIndicesName.c_str());
      exit(1);
    }
    fwrite(&doubleIntegralCountIncremental,sizeof(int), 1, binaryDataFile);
    fwrite(&doubleOrbitalIntegralIndices[0], sizeof(int),doubleIntegralCountIncremental, binaryDataFile);
    fclose(binaryDataFile);
    //
    // write out the number of integrals output and then the integrals
    //
    if( (binaryDataFile = fopen(doubleIntegralsName.c_str(), "wb" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",doubleIntegralsName.c_str());
      exit(1);
    }
    fwrite(&doubleIntegralCountIncremental,sizeof(int), 1, binaryDataFile);
    fwrite(&doubleOrbitalIntegralData[0], sizeof(double),doubleIntegralCountIncremental, binaryDataFile);
    fclose(binaryDataFile);


    printf("\n#######################################################\n\n");
    printf("Orbital Integral Data Files Output (binary format)  \n");
    printf("=== %s   \n",singleIndicesName.c_str());
    printf("=== %s   \n",singleIntegralsName.c_str());
    printf("=== %s   \n",doubleIndicesName.c_str());
    printf("=== %s   \n",doubleIntegralsName.c_str());


    char* fullPath;
    printf("\n");
    fullPath =  realpath(singleIndicesName.c_str(),NULL);
    printf("singleIndicesPath : %s \n",fullPath);
    free(fullPath);
    fullPath =  realpath(singleIntegralsName.c_str(),NULL);
    printf("singleIntegralsPath : %s \n",fullPath);
    free(fullPath);
    fullPath =  realpath(doubleIndicesName.c_str(),NULL);
    printf("doubleIndicesPath : %s \n",fullPath);
    free(fullPath);
    fullPath =  realpath(doubleIntegralsName.c_str(),NULL);
    printf("doubleIntegralsPath : %s \n",fullPath);
    free(fullPath);
    printf("\n#######################################################\n");
}

void progressMessage(long dataIndex,long totalIndexSize)
{
    long percentDone = long ((double(dataIndex)/double(totalIndexSize))*10.0);
    if(percentDone > orbitalProgressValue)
    {
    printf("Completed %5ld %% of the required integrals \n",percentDone*10);
    orbitalProgressValue = percentDone;
    }
}

void outputTimingData()
{
	totalTime = globalClockIt.getMilliSecElapsedTime();
    printf("Clock Time Taken   (minutes)         : %15.6f  \n",totalTime/(60000.0));
    printf("Startup         Time Taken (seconds) : %15.8e  \n",startupTime/(1000.0));
    printf("Single Integral Time Taken (seconds) : %15.8e  \n",singleOrbitalTime/(1000.0));
    printf("Double Integral Time Taken (seconds) : %15.8e  \n",doubleOrbitalTime/(1000.0));
    printf("OrbitalIntegrals_Time_M              : %15.6f  \n",totalTime/(60000.0));
}

//
// This routine times the dot product component of the double orbital integral
// calculation.
//
double timeDots(long repetitions)
{
	vector<UCLAQ::GridFunction3d>& basisVectors = *basisVectorsPtr;

    long i; long j;
    long m; long n;

    UCLAQ::GridFunction3d oFunIxJ(basisVectors[0]);
    UCLAQ::GridFunction3d oFunMxN(basisVectors[0]);

    oFunIxJ  = basisVectors[0];
    oFunMxN  = basisVectors[1];
    double dotVal;

    clockIt.start();
    long k;
    for(k = 1; k <= repetitions; k++)
    {
    for(i = 1; i <= orbFunBasisCount; i++)
    {
    for(j = 1; j <= i; j++)
    {
    	for(m = 1; m  < i;  m++)
    	{
        for(n = 1; n <= m;  n++)
        {
        dotVal = oFunMxN.scaledDot(oFunIxJ);
        }}

    	m = i;
    	for(n = 1; n <= j;  n++)
    	{
        dotVal = oFunMxN.scaledDot(oFunIxJ);
        }
    }}

    }
    clockIt.stop();
	double dotTime = (clockIt.getMilliSecElapsedTime())/((double)(repetitions));
    dotTime        =  dotTime/double(doubleOrbitalIntegralCount);
	return dotTime;
}
//
// This routine times the grid function product component of the
// double orbital integral calculation.
//
double timeProducts(long repetitions)
{
    vector<UCLAQ::GridFunction3d>& basisVectors = *basisVectorsPtr;

    long i; long j;
    long m; long n;

    UCLAQ::GridFunction3d oFunIxJ(basisVectors[0]);
    UCLAQ::GridFunction3d oFunMxN(basisVectors[0]);

    oFunIxJ  = basisVectors[0];
    oFunMxN  = basisVectors[1];
    double dotVal;

    clockIt.start();
    long k;
    for(k = 1; k <= repetitions; k++)
    {
    for(i = 1; i <= orbFunBasisCount; i++)
    {
    for(j = 1; j <= i; j++)
    {
    	oFunIxJ  = basisVectors[i-1];
        oFunIxJ *= basisVectors[j-1];
    	for(m = 1; m  < i;  m++)
    	{
        for(n = 1; n <= m;  n++)
        {
        oFunMxN  = basisVectors[m-1];
        oFunMxN *= basisVectors[n-1];
        }}

    	m = i;
    	for(n = 1; n <= j;  n++)
    	{
        oFunMxN  = basisVectors[m-1];
        oFunMxN *= basisVectors[n-1];
        }
    }}

    }
    clockIt.stop();
    double prodTime = (clockIt.getMilliSecElapsedTime())/((double)(repetitions));
    prodTime        = prodTime/double(doubleOrbitalIntegralCount + singleOrbitalIntegralCount);
	return prodTime;
}

   // void run();

    //void readParameterLists();
    //void readParameterLists(const char* parameterInputFile);
    //void readParameterLists(XML_ParameterListArray& paramsInput);

    //void getExtenralPotential();
    //void getOrbitalBasis();

    void setTimingConstantsFlag(bool flagValue = true)
    {timingConstantsFlag = flagValue;}

    void clearTimingConstantsFlag()
    {timingConstantsFlag = false;}
//
//  Input/output parameters
//
    long        orbFunBasisCount;
    string   orbFunBasisFileName;
    string        outputFileName;
    bool        binaryOutputFlag;
    string  programInputFileName;
//
//  Parameters required to set up the 3D grid
//
	double xMin; double xMax;   long xPanel;
	double yMin; double yMax;   long yPanel;
	double zMin; double zMax;   long zPanel;
//
//  Grid functions for computing single and double orbital integrals
//
	UCLAQ::GridFunction3d oFunI;
	UCLAQ::GridFunction3d oFunJ;

	UCLAQ::GridFunction3d oFunIxJ;
    UCLAQ::GridFunction3d oFunMxN;

//
//  The use of a pointer for the orbital basis vectors
//  facilitates the use of this class to create orbitals when
//  basis functions are externally managed.
//
    vector<double>                    basisIndexValues;
    vector< UCLAQ::GridFunction3d >*    basisVectorsPtr;

//  Pointer to single particle Schroedinger operator

    SingleParticleOp* singleParticleSop;
//
//  Inverse Laplace operator and requisite parameter
//
    double          invLapExtFactor;
    InverseLaplaceOp invLaplaceOp3d;

//
//  Coefficient of Laplace operator so that inverse Laplacian
//  with infinite boundary conditions yields electron-electron
//  interaction potential.
//
//  The default value = 1/(4*pi) is the coefficient associated with
//  the Schroedinger operator in atomic units.
//
    double     laplaceOpCoefficient;
//
//  Parameters for parallel construction
//
    bool incrementalFlag;
    long iStart;
    long iEnd;
    long jStart;
    long jEnd;
    long ijIndexStart;
    long ijIndexEnd;
    long doubleIntegralCountIncremental;
    long singleIntegralCountIncremental;
    vector<int> singleOrbitalIntegralIndices;
    vector<int> doubleOrbitalIntegralIndices;

//
// Orbital data 
//
    long          singleOrbitalIntegralCount;
    vector<double> singleOrbitalIntegralData;

    long          doubleOrbitalIntegralCount;
    vector<double> doubleOrbitalIntegralData;

    long                orbitalProgressValue;
    
//
// Timing data and routines
//
    
    ClockIt             clockIt;
    ClockIt       globalClockIt;
    bool    timingConstantsFlag;
	double   singleOrbitalTime;
	double   doubleOrbitalTime;
	double         startupTime;
	double           totalTime;

	bool verboseFlag;
};

} // UCLAQ namespace

#undef  _DEFAULT_INV_LAPLACIAN_EXTENSION_FACTOR_
#undef  _DEFAULT_LAPLACE_OP_ORDER_
#undef  _DEFAULT_OUTPUT_FILENAME_
#undef  _DEFAULT_laplaceOpCoefficient
#endif



 
