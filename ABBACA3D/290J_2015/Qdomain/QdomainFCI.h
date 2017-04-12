/*
 * QdomainFCI.h
 *
 *  Created on: Jun 30, 2014
 *      Author: anderson
 */
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

#ifndef  _QdomainFCI_
#define  _QdomainFCI_

#include <sstream>
#include <iostream>
#include <cstdio>
#include <stdexcept>

#include "MultiParticleHamiltonian/RealSpatialOrbitalIntegrals.h"
#include "MultiParticleHamiltonian/MultiParticleHamiltonianMatrix.h"
#include "RayleighChebyshev/RayleighChebyshev.h"
#include "RandOpNd/RandOp1d.h"

#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"

#include "Timing/ClockIt.h"

#define _DEFAULT_DENSITY_FILE_PREFIX_ "FCI"

#define  _DEFAULT_eigCountMax            5
#define  _DEFAULT_subspaceTol            5.0e-06
#define  _DEFAULT_subspaceIncrementSize  10
#define  _DEFAULT_bufferSize             5

namespace UCLAQ
{


class QdomainFCI
{
public:

	QdomainFCI()
    {
		initialize();
    }

    void initialize()
    {
    electronCount        = 0;
    spatialOrbBasisCount = 0;
    binaryIntegralsFlag  = false;
    iterativeEigFlag     = 0;
    dropTolerance        = 0.0;

    nuclearEnergyShift      = 0.0;
    outputEigenvectorFlag   = false;
    verboseEigensystemFlag  = false;
    energyDecompositionFlag = true;

    rOrbitalsPtr             = 0;
    verboseFlag              = false;
    verboseEigFlag           = false;


    eigCountMax           = _DEFAULT_eigCountMax;
    subspaceTol           = _DEFAULT_subspaceTol;
	subspaceIncrementSize = _DEFAULT_subspaceIncrementSize;
	bufferSize            = _DEFAULT_bufferSize;
    }

    void setVerbose(bool flagValue = true)
    { verboseFlag = flagValue;}

    void clearVerbose()
    { verboseFlag =  false;}

    void setVerboseEigensystemComputation(bool flagValue = true)
    { verboseEigFlag  = flagValue;}

    void clearVerboseEigensystemComputation()
    { verboseEigFlag  =  false;}


    void setEigCountMax(long eigCountMax)
    {
    this->eigCountMax = eigCountMax;
    }

    void setSubspaceTol(double subspaceTol)
    {
    this->subspaceTol = subspaceTol;
    }

    void setSubspaceIncrementSize(long subspaceIncrementSize)
    {
    this->subspaceIncrementSize = subspaceIncrementSize;
    }

    void setBufferSize(long bufferSize)
    {
    this->bufferSize = bufferSize;
    }


    void run(long electronCount, long spatialOrbitalBasisCount, double nuclearEnergyShift,
    RealSpatialOrbitalIntegrals& rOrbitalIntegrals, vector < string >& spinBlockIdentifiers,
    vector <vector <double> >& spinBlockEnergies, vector < vector <UCLAQ::DoubleVector2d> >& spinBlockDensities)
    {

    if(verboseFlag)
    {
    printf("\n");
    printf("#############\n");
    printf("############# QdomainFCI \n");
    printf("#############\n");
    printf("\n");
    }

    string errMessage;
    std::ostringstream errStream;

    // Localize input parameters

    this->electronCount        = electronCount;
    this->spatialOrbBasisCount = spatialOrbitalBasisCount;
    this->nuclearEnergyShift   = nuclearEnergyShift;
    this->rOrbitalsPtr         = &rOrbitalIntegrals;

    // Start the clock

    globalClockIt.start();

    RealSpatialOrbitalIntegrals& rOrbitals = *rOrbitalsPtr;

    long basisCount = rOrbitals.getOrbitalBasisCount();

    if(spatialOrbBasisCount == -1)
    {
    spatialOrbBasisCount = basisCount;
    }

    if(verboseFlag)
    {
    	printf("Electron count                         : %ld \n",electronCount);
    	printf("Number of spatial orbitals in data set : %ld \n", basisCount);
    	printf("--- Orbital Indices ---  \n");

    	for(long i = 0; i < basisCount; i++)
    	{
    		printf("%-10.5g \n",rOrbitals.getOrbitalEnergy(i));
    	}
    	printf("\n");
    }
//
//  Number of electrons
//
    long N  = electronCount;
//
//  Create an array of Spin Orbitals
//
    long MspatialOrbitals = spatialOrbBasisCount;

    if(verboseFlag)
    {
    printf("Number of spatial orbitals used in computation : %ld \n",spatialOrbBasisCount);
    }

    if(basisCount < spatialOrbBasisCount)
    {
    errMessage  = "\nXXX QdomainFCI XXX\n";
	errMessage +="Insufficient number of basis functions in data set \n\n";
	errStream << spatialOrbBasisCount;
	errMessage +="Requested basis size: " + errStream.str() + "\n";
	errStream.str("");
	errStream << basisCount;
    errMessage +="Available basis size: " + errStream.str() + "\n";
	throw std::runtime_error(errMessage);
    }
//
//  Input spatial orbitals and energies
//
    double* spatialOrbitalEnergies = new double[MspatialOrbitals];

    for(long i = 0; i < MspatialOrbitals; i++)
    {
    spatialOrbitalEnergies[i] = float(i+1);
    }
//
//  Create an array of SpinOrbitals from the spatialOrbitals
//
    long M = 2*MspatialOrbitals;

    SpinOrbital* spinOrbitals = new SpinOrbital[M];

    long     spatialOrbitalIndex;
    int                     spin;
    double                energy;

    spin = 1;
    for(long i= 0; i < MspatialOrbitals; i++)
    {
    spatialOrbitalIndex = i;
    energy              = spatialOrbitalEnergies[i];
    spin                =  1;
    spinOrbitals[2*i].initialize(spatialOrbitalIndex,spin,energy);

    spin                = -1;
    spinOrbitals[2*i+1].initialize(spatialOrbitalIndex,spin,energy);
    }
//
//  Create a Slater determinant basis from the list of spin orbitals
//
    SlaterDeterminantBasis sdBasis;
    sdBasis.initialize(N,M,spinOrbitals);
    long sdSize = sdBasis.getBasisCount();

    if(verboseFlag)
    {
    printf("Slater Determinant Basis Size: %ld \n",sdSize);
    }
//
//  Construct the Hamiltonian matrix
//
    MultiParticleHamiltonianMatrix hMatrix(sdBasis); // attach Slater determinant basis
    ElectronOperators eOp(rOrbitals);                // attach single and double electron
                                                     // operators that use rOrbital info.

    hMatrix.clearVerboseFlag();
    if(verboseFlag)
    {
    hMatrix.setVerboseFlag();
    }

    long k; long m; long n;

    if(verboseFlag)
    {
    for(k = 0; k < hMatrix.spinBlockCount; k++)
    {
    printf("Spin block: %ld   Spin Block Size:  %ld \n",k,sdBasis.getSpinBlockSize(k));
    }
    }

    long blockStart = 0;
    long blockEnd   = hMatrix.spinBlockCount / 2 + hMatrix.spinBlockCount % 2 - 1;

    long    eigCount;

    UCLAQ::DoubleVector1d Vstar;
    UCLAQ::DoubleVector2d     A;

    UCLAQ::DoubleVector1d   energyKandP;
    UCLAQ::DoubleVector1d energyCoulomb;

    UCLAQ::DoubleVector2d densityMatrix;
    densityMatrix.initialize(spatialOrbBasisCount, spatialOrbBasisCount);
    vector < UCLAQ::DoubleVector2d > densityMatrices;


    vector<double>                      eigValues;
    vector <UCLAQ::DoubleVector1d>     eigVectors;

    UCLAQ::RandOp1d                  randOp;
    UCLAQ::DoubleVector1d              vTmp;

    RayleighChebyshev <UCLAQ::DoubleVector1d, SparseSpinBlockMatrix<UCLAQ::DoubleVector1d> , UCLAQ::RandOp1d>  rayleighChebyshev;
    long rcSubspaceSize;
    long  rcBufferSize;

    if(verboseEigFlag)
    {
    rayleighChebyshev.setVerboseFlag();
    }

    SparseSpinBlockMatrix<UCLAQ::DoubleVector1d> sparseSpinBlockMatrix;

	int opType;
	string        spinBlockIdentifier;
	ostringstream eOutput;

    for(k = blockStart; k <= blockEnd;  k++)
    {
    // Null spin block

    if(sdBasis.getSpinBlockSize(k) == 0)
    {
    densityMatrices.clear();
    eigValues.clear();
    spinBlockEnergies.push_back(eigValues);
    spinBlockDensities.push_back(densityMatrices);
    spinBlockIdentifiers.push_back("NULL");
    }

    // Non-null spin block

    if(sdBasis.getSpinBlockSize(k) > 0)
    {

    //
    // construct spin block matrix
    //
    hMatrix.initializeSparseSpinBlockMatrix(k,eOp,dropTolerance,sparseSpinBlockMatrix);
    m = sdBasis.getSpinBlockSize(k);

    if(m == 1)
    {
    eigCount = 1;
    vTmp.initialize(1);
    eigValues.resize(1);
    eigValues[0] =  sparseSpinBlockMatrix(0,0);
    eigVectors.resize(1,vTmp);
    eigVectors[0].setToValue(1.0);
    }
    else
    {
    rcSubspaceSize = subspaceIncrementSize;
    rcBufferSize   = bufferSize;
    eigCount = eigCountMax;

    if(eigCount > m)
    {
    	eigCount = m;
    	rcSubspaceSize = m;
    	rcBufferSize   = 0;
    }
    if(m < subspaceIncrementSize + bufferSize )
    {
        rcSubspaceSize = m;
    	rcBufferSize   = 0;
    }

    vTmp.initialize(m);

    if(verboseFlag)
    {
    printf("Iterative Eigensystem Construction Start \n");
    rayleighChebyshev.setEigDiagnosticsFlag();
    }
    if(verboseEigensystemFlag)
    {
	rayleighChebyshev.setVerboseFlag();
    }

    eigValues.clear();
    eigVectors.clear();

    clockIt.start();

    cout << " Eig count XXXX "  << eigCount << endl;

	rayleighChebyshev.getMinEigenSystem(eigCount, subspaceTol, rcSubspaceSize, rcBufferSize, vTmp,sparseSpinBlockMatrix, randOp, eigValues, eigVectors);

	clockIt.stop();
	if(verboseFlag)
	{
    printf("Iterative Eigensystem Construction End : Time (minutes)  : %-10.4f \n",(clockIt.getMilliSecElapsedTime())*(1.0/60000.0));
    }
    }

    if(verboseFlag)
    {
    printf("\n##### Block %2ld  ###### \n",k);
    }

    spinBlockIdentifiers.push_back(getSpinBlockIdentifier(electronCount, k, hMatrix));

    spinBlockEnergies.push_back(eigValues);
    for(long i = 0; i < eigCount; i++)
    {
    spinBlockEnergies.back()[i] += nuclearEnergyShift;
    }


    if(verboseFlag)
    {
    for(long i = 0; i < eigCount; i++)
    {
    printf("%10.8e  \n",eigValues[i] + nuclearEnergyShift);
    }
    printf("\n\n");

	printf("\n\nTotal Energy \n");
	for(long i = 0; i < eigCount; i++)
	{
		eOutput.str("");
		spinBlockIdentifier =  getSpinBlockIdentifier(electronCount, k, hMatrix);
		eOutput << spinBlockIdentifier << "_Energy_" << i + 1 << " : ";
		printf("%s %20.15e \n",eOutput.str().c_str(),eigValues[i] + nuclearEnergyShift);
	}
	}

    //
    // Capture spin block densities
    //

    densityMatrices.clear();
    densityMatrices.resize(eigCount,densityMatrix);

    for(long i = 0; i < eigCount; i++)
	{
	densityMatrix.setToValue(0.0);
    hMatrix.createDensityMatrix(k, eigVectors[i], densityMatrix);
    densityMatrices[i] = densityMatrix;
    }

    spinBlockDensities.push_back(densityMatrices);


    /*
	outputDensityMatrices(densityFilePrefix.c_str(),spatialOrbBasisCount,electronCount,k,eigValues,
	eigVectors,hMatrix);


	if(outputEigenvectorFlag)
	{
		outputEigenSystem(densityFilePrefix.c_str(),spatialOrbBasisCount,electronCount,k,sdBasis, spatialOrbitalEnergies,
		eigValues, eigVectors);
	}
	*/


	//
	// Create and output the energy decomposition
	//

	if(energyDecompositionFlag)
	{
		energyKandP.initialize(eigCount);
		energyCoulomb.initialize(eigCount);

		m = sdBasis.getSpinBlockSize(k);
		n = m;
		A.initialize(m,n);
		Vstar.initialize(m);

	    opType = MultiParticleHamiltonianMatrix::ONE_PTCLE_OP;
	    hMatrix.clearVerboseFlag();


		hMatrix.initializeSpinBlockMatrices(eOp,opType);


		for(long i = 0; i < m; i++)
		{
		for(long j = 0; j < n; j++)
		{
			A(i,j) = hMatrix.spinBlockMatrix[k](i,j);
		}}

		for(long j = 0; j < eigCount; j++)
		{
		for(long i = 0; i < m; i++)
		{
			Vstar(i) = eigVectors[j](i);
		}
		energyKandP(j) = Vstar.dot(AtimesV(A,Vstar));
		}

    	opType = MultiParticleHamiltonianMatrix::TWO_PTCLE_OP;
		hMatrix.initializeSpinBlockMatrices(eOp,opType);

		for(long i = 0; i < m; i++)
		{
		for(long j = 0; j < n; j++)
		{
			A(i,j) = hMatrix.spinBlockMatrix[k](i,j);
		}}

		for(long j = 0; j < eigCount; j++)
		{
		for(long i = 0; i < m; i++)
		{
			Vstar(i) = eigVectors[j](i);
		}
		energyCoulomb(j) = Vstar.dot(AtimesV(A,Vstar));
		}

		if(verboseFlag)
		{
		printf("\n\nFCI Energy Decomposition\n\n");
		printf("Kinetic and Potential (a.u.)  :  E-E   (a.u.)        : Total (a.u.) \n");
		for(long i = 0; i < eigCount; i++)
		{
		printf("%10.8e               :  %10.8e      : %10.8e \n",
			energyKandP(i) + nuclearEnergyShift ,energyCoulomb(i),energyKandP(i) + energyCoulomb(i) + nuclearEnergyShift);
		}
		printf("\n");


		printf("\n\nKinetic + Potential Energy \n");
		for(long i = 0; i < eigCount; i++)
		{
		eOutput.str("");
		spinBlockIdentifier =  getSpinBlockIdentifier(electronCount, k, hMatrix);
		eOutput << spinBlockIdentifier << "_KP_Energy_" << i +1  << " : ";
		printf("%s %20.15e \n",eOutput.str().c_str(),energyKandP(i) + nuclearEnergyShift);
		}

		printf("\n\nElectron-Electron Energy \n");
		for(long i = 0; i < eigCount; i++)
		{
		eOutput.str("");
		spinBlockIdentifier =  getSpinBlockIdentifier(electronCount, k, hMatrix);
		eOutput << spinBlockIdentifier << "_EE_Energy_" << i +1 << " : ";
		printf("%s %20.15e \n",eOutput.str().c_str(),energyCoulomb(i) );
		}

		printf("\nXXXXXXXXXXXXXXXXXXX \n");
	    hMatrix.setVerboseFlag();
	    }
	}

    }} // end of spin block loop

    globalClockIt.stop();

    if(verboseFlag)
    {
	printf("FCI_Clock_Time_Taken_M         : %-10.6f \n",(globalClockIt.getMilliSecElapsedTime())*(1.0/60000.0));
	}
//
//  Clean up
//

    if(spatialOrbitalEnergies != 0) delete [] spatialOrbitalEnergies;
    if(spinOrbitals != 0)           delete [] spinOrbitals;

    }


string getSpinBlockIdentifier(long electronCount, int spinBlockIndex,
MultiParticleHamiltonianMatrix& hMatrix)
{
    string spinBlockIdentifier;

    //
    // Construct identifier for spin block consisting of U's and D's for up and down
    // spin
    //

    SlaterDeterminantReference sI;
    long IblockIndex;
    long IspinBlockGlobalIndex;

    long i = 0;
    IspinBlockGlobalIndex = hMatrix.sdBasisPtr->getSpinBlockGlobalIndex(spinBlockIndex,i);
    IblockIndex = hMatrix.sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex];
    sI          = hMatrix.sdBasisPtr->operator()(IblockIndex);


    int upSpinCount   = 0;
    int downSpinCount = 0;

    for(long j = 0; j < electronCount; j++)
    {
    if(sI.getOrbitalSpin(j) == 1) {upSpinCount++; }
    else                          {downSpinCount++;}
    }
    for(long j = 0; j < upSpinCount; j++)
    {spinBlockIdentifier.append("U");}

    for(long j = 0; j < downSpinCount; j++)
    {spinBlockIdentifier.append("D");}

    return spinBlockIdentifier;
}

UCLAQ::DoubleVector1d AtimesV(UCLAQ::DoubleVector2d& A, UCLAQ::DoubleVector1d& X)
{
	UCLAQ::DoubleVector1d Y(X);
	double sum = 0.0;
	long iMin = 0;
	long iMax = A.getIndex1Size();
	long jMin = 0;
	long jMax = A.getIndex2Size();

	long i; long j;
	for(i = iMin; i < iMax; i++)
	{
		sum = 0.0;
		for(j = jMin; j <  jMax; j++)
		{
			sum += A(i,j)*X(j);
		}
		Y(i) = sum;
	}
	return Y;

}
void outputDensityMatrices(const char* filePrefix,long spatialOrbitalCount, long electronCount,
int spinBlockIndex,vector <double>& eigValues,vector <UCLAQ::DoubleVector1d>& eigVectors,
MultiParticleHamiltonianMatrix& hMatrix)
{
//
// This routine outputs the density matrices associated with an FCI
// computation.
//
// Existing data files are overwritten.
//
// The output file name for the eigenvalues is of the form
//
// FilePrefix_SbX_denMatrixY_M_N.dat
//
// This file contains the spin associated with each spin block
// and the eigenvalues for each spin block.
//
// The orbitalIntegrals.dat file name is also recorded.
//
// The output file name for each eigenvector is
// of the form FilePrefix_SbX_eigvecY_M_N.dat
//
// X = index of spin block
//
// Y = index of eigenvector (starting with 1), and
// printing out the lowest eigenvectors of each
// spin block.
//
// M = number of spatial (not spin) orbitals
//
// N = number of particles
//
// For output use a fixed format ASCII file (so I don't
// have to worry about difference in integer storage
// formats).
//
    FILE* dataFile;

    struct tm when;
    time_t now; time(&now);
    when = *localtime( &now );

    UCLAQ::DoubleVector2d densityMatrix;
    densityMatrix.initialize(spatialOrbitalCount,spatialOrbitalCount);
    densityMatrix.setToValue(0.0);

	ostringstream s;
    string densityMatrixFileName;
    string spinBlockIdentifier;
    printf("\nDensity Matrices Output To Files: \n\n");

    long i; long j; long k;
    for(k = 1; k <= (long)eigVectors.size(); k++)
    {
    densityMatrix.initialize(spatialOrbitalCount,spatialOrbitalCount);
    densityMatrix.setToValue(0.0);
    hMatrix.createDensityMatrix(spinBlockIndex, eigVectors[k-1], densityMatrix);

    densityMatrixFileName.clear();
    densityMatrixFileName.assign(filePrefix);
    densityMatrixFileName.append("_Sb");
    s.str("");
    s << spinBlockIndex;
    densityMatrixFileName.append((s.str()).c_str());
    densityMatrixFileName.append("_denMatrix");
    s.str("");
    s << k;
    densityMatrixFileName.append((s.str()).c_str());
    s.str("");
    s << "_" << spatialOrbitalCount << "_" << electronCount;
    densityMatrixFileName.append((s.str()).c_str());
    densityMatrixFileName.append(".dat");



    //
    // Construct identifier for spin block consisting of U's and D's for up and down
    // spin
    //

    SlaterDeterminantReference sI;
    long IblockIndex;
    long IspinBlockGlobalIndex;

    i = 0;
    IspinBlockGlobalIndex = hMatrix.sdBasisPtr->getSpinBlockGlobalIndex(spinBlockIndex,i);
    IblockIndex = hMatrix.sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex];
    sI          = hMatrix.sdBasisPtr->operator()(IblockIndex);

    spinBlockIdentifier.clear();

    int upSpinCount   = 0;
    int downSpinCount = 0;

    for(j = 0; j < electronCount; j++)
    {
    if(sI.getOrbitalSpin(j) == 1) {upSpinCount++; }
    else                          {downSpinCount++;}
    }
    for(j = 0; j < upSpinCount; j++)
    {spinBlockIdentifier.append("U");}

    for(j = 0; j < downSpinCount; j++)
    {spinBlockIdentifier.append("D");}


    printf("State %ld: %s (%s) \n",k-1,densityMatrixFileName.c_str(),spinBlockIdentifier.c_str());
   	//
   	// Open file and output data
   	//
    if( (dataFile = fopen(densityMatrixFileName.c_str(), "w+" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",densityMatrixFileName.c_str());
      return;
    }
    fprintf(dataFile,"Spin Block     : %s  \n",spinBlockIdentifier.c_str());
    fprintf(dataFile,"Energy         : %20.15e \n",eigValues[k-1]);
    fprintf(dataFile,"Electron Count : %ld \n",electronCount);
    fprintf(dataFile,"Orbital Count  : %ld \n",spatialOrbitalCount);
    fprintf(dataFile,"FCI Run Date     : %s", asctime( &when ) );
    fprintf(dataFile,"FCI Code Version : %s \n", "XXXX");
    fprintf(dataFile,"Density Values   :   \n");

    for(i = 0;  i < spatialOrbitalCount; i++)
    {
		for(j = 0; j < spatialOrbitalCount; j++)
		{
    	fprintf(dataFile,"%20.15e ",densityMatrix(i,j));
    	}
    	fprintf(dataFile,"\n");
    }

    fclose(dataFile);
    }
    printf("\n");
}

void outputEigenSystem(const char* filePrefix,long Morbital, long Nparticle,
int spinBlockIndex, SlaterDeterminantBasis& sdBasis,double* spatialOrbitalEnergies,
vector<double>& eigValues,vector <UCLAQ::DoubleVector1d>& eigVectors)
{
//
// This routine outputs the eigenvectors and eigenvalues
// associated with the FCI computation. The eigenvectors
// are output using the standard ordering of the Slater
// determinant basis as implemented by the SlaterDeterminantBasis
// class.
//
// Existing data files are overwritten.
//
// The output file name for the eigenvalues is of the form
//
// FilePrefix_eigval_M_N.dat
//
// This file contains the spin associated with each spin block
// and the eigenvalues for each spin block.
//
// The orbitalIntegrals.dat file name is also recorded.
//
// The output file name for each eigenvector is
// of the form FilePrefix_SbX_eigvecY_M_N.dat
//
// X = index of spin block
//
// Y = index of eigenvector (starting with 1), and
// printing out the lowest eigenvectors of each
// spin block.
//
// M = number of spatial (not spin) orbitals
//
// N = number of particles
//
// For output use a fixed format ASCII file (so I don't
// have to worry about difference in integer storage
// formats).
//
// Print out size of the vector
// Print out data entries
//
//
    FILE* dataFile;
    SlaterDeterminantReference R;

	ostringstream s;
    string eigenvalueFileName;
    string eigenvectorFileName;

    eigenvalueFileName.clear();
    eigenvalueFileName.assign(filePrefix);
    eigenvalueFileName.append("_eigval_");
    s.str("");
    s << Morbital << "_" << Nparticle;
    eigenvalueFileName.append((s.str()).c_str());
    eigenvalueFileName.append(".dat");

    long basisSize = eigVectors[0].getSize();
    long i;  long k;
    for(k = 1; k <= (long)eigVectors.size(); k++)
    {
    eigenvectorFileName.clear();
    eigenvectorFileName.assign(filePrefix);
    eigenvectorFileName.append("_Sb");
    s.str("");
    s << spinBlockIndex;
    eigenvectorFileName.append((s.str()).c_str());
    eigenvectorFileName.append("_eigvec");
    s.str("");
    s << k;
    eigenvectorFileName.append((s.str()).c_str());
    s.str("");
    s << "_" << Morbital << "_" << Nparticle;
    eigenvectorFileName.append((s.str()).c_str());
    eigenvectorFileName.append(".dat");
    cout << eigenvectorFileName << endl;


    if( (dataFile = fopen(eigenvectorFileName.c_str(), "w+" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",eigenvectorFileName.c_str());
      return;
    }

    fprintf(dataFile,"%ld \n",basisSize);
    for(i = 0;  i < basisSize; i++)
    {
    	/*
    	basisIndex = sdBasis.getSpinBlockGlobalIndex(spinBlockIndex,i);
    	stdIndex   = sdBasis.spinBlockToStd[basisIndex];
    	R = sdBasis(stdIndex);
    	energySum = 0.0;
    	for(j = 0; j < Nparticle; j++)
    	{
    	spatialIndex = R.getSpatialOrbitalIndex(j);
    	printf("%ld ",spatialIndex);
    	energySum += spatialOrbitalEnergies[spatialIndex];
        }
    	printf(" %20.15e  %20.15e \n",energySum,eigVectors[k](i));

    	*/
    	fprintf(dataFile,"%20.15e \n",eigVectors[k-1](i));
    }
    fclose(dataFile);
    }

}

    void outputA(UCLAQ::DoubleVector2d& A)
    {
    long i; long j;

    long m = A.getIndex1Size();
    long n = A.getIndex2Size();

    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
        printf("%6.3e ",A(i,j));
        }
        printf("\n");
    }

    }


    long             eigCountMax;
    double           subspaceTol;
	long    subspaceIncrementSize;
	long               bufferSize;


    bool     verboseFlag;
    bool  verboseEigFlag;

    ClockIt clockIt;
    ClockIt globalClockIt;

    RealSpatialOrbitalIntegrals* rOrbitalsPtr;

    long        electronCount;
    long spatialOrbBasisCount;
    bool  binaryIntegralsFlag;
    bool     iterativeEigFlag;
    double     dropTolerance;


    double nuclearEnergyShift;
    bool outputEigenvectorFlag;
    bool verboseEigensystemFlag;
    bool energyDecompositionFlag;

    string     densityFilePrefix;

};
} // UCLAQ namespace

#undef _DEFAULT_eigCountMax
#undef _DEFAULT_subspaceTol
#undef _DEFAULT_subspaceIncrementSize
#undef _DEFAULT_bufferSize

#endif /*  _QdomainFCI_*/
