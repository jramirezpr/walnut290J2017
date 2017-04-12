//
// MultiparticleHamiltonainMatrix.h
//
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
#ifndef __MultiParticleHamiltonianMatrix__
#define __MultiParticleHamiltonianMatrix__

#include <iostream>
#include <cmath>
using namespace std;

#include "SpinOrbital.h"
#include "SlaterDeterminantReference.h"
#include "SlaterDeterminantBasis.h"
#include "SpinOrbitalUtilities.h"


#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"
#include "LongArrayNd/UCLAQ_LongArray1d.h"
#include "LongArrayNd/UCLAQ_LongArray2d.h"

#include "ElectronOperators.h"
#include "SparseSpinBlockMatrix.h"

namespace UCLAQ
{

class MultiParticleHamiltonianMatrix
{
    public: 

	MultiParticleHamiltonianMatrix();

    MultiParticleHamiltonianMatrix(SlaterDeterminantBasis& sdB);
    ~MultiParticleHamiltonianMatrix();
//
//  Dense matrix representation
//
//  OpType =  ONE_PTCLE_OP or TWO_PTCLE_OP or ONE_AND_TWO_PTCLE_OP

    void initializeSpinBlockMatrices(ElectronOperators& eOp, int opType);
    
    void initializeSpinBlockMatrices_FullLoops();
//
//  Sparse matrix representation
//
//  OpType =  ONE_PTCLE_OP or TWO_PTCLE_OP or ONE_AND_TWO_PTCLE_OP

    void initializeSparseSpinBlockMatrices(ElectronOperators& eOp,
    double dropTolerance);
    
    void initializeSparseSpinBlockMatrix(int spinBlockIndex, 
	ElectronOperators& eOp, double dropTolerance, SparseSpinBlockMatrix<UCLAQ::DoubleVector1d>& sparseSpinBlockMatrix);
//
//  Matrix free representation
//
    void setMatrixFreeSpinBlock(long i) 
    {matrixFreeSpinBlock = i;}

    void setElectronOp(ElectronOperators& electronOp)
    {this->eOpPtr = &electronOp;}

    void applyForwardOp(UCLAQ::DoubleVector1d& V);


	void createDensityMatrix(int spinBlockIndex, UCLAQ::DoubleVector1d& eigVector,
	UCLAQ::DoubleVector2d& densityMatrix);
//
    void outputMatrix();
    
    void outputSpinBlockMatrices(const char* outputFilePrefix, 
    const char* outputFormatString,const char* inputIntegralDataFile);

    void outputSpinBlockMatricesKformat(const char* outputFilePrefix, 
    const char* outputFormatString,const char* inputIntegralDataFile);
 
    void destroyData();

    void setVerboseFlag(bool flagValue = true)
    {verboseFlag = flagValue;}

    void clearVerboseFlag()
    {verboseFlag = false;}

    enum {ONE_PTCLE_OP, TWO_PTCLE_OP, ONE_AND_TWO_PTCLE_OP};

    SlaterDeterminantBasis*        sdBasisPtr;
    SpinOrbitalUtilities spinOrbitalUtilities;

    long* orbitalIndexList;
    long* complimentaryOrbitalIndexList;

    long sdSize;
    long localSpinBlockSize;
    long spinBlockCount;
//
//  Dense matrix representation data
//
    UCLAQ::DoubleVector2d*                 spinBlockMatrix;
//
//  Sparse matrix representation data
//
    SparseSpinBlockMatrix<UCLAQ::DoubleVector1d>* sparseSpinBlockMatrix;
//
//  Matrix Free Implementation data
//
    long matrixFreeSpinBlock;     // sets spin block index for matrix free operator
    UCLAQ::DoubleVector1d vTmp;                // temporary array for matrix free operator
    ElectronOperators*    eOpPtr; // pointer to ElectronOperator integrals
//
//  Utility routines
//
	int compareOrbitalIndices(long iIndex,long jIndex, long nParticle,
	UCLAQ::LongArray2d& orbitalIndices);

	bool verboseFlag;

};

} // UCLAQ namespace

#endif
 
