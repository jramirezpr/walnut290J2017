#include <iostream>
#include <sstream>
#include <string>
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

#include <cmath>
#include <stdexcept>
using namespace std;

#include "Timing/ClockIt.h"
#include "MultiParticleHamiltonianMatrix.h"

namespace UCLAQ
{


MultiParticleHamiltonianMatrix::MultiParticleHamiltonianMatrix()
{
	spinBlockMatrix               = 0;
	orbitalIndexList              = 0;
    complimentaryOrbitalIndexList = 0;
    sparseSpinBlockMatrix         = 0;
    sdBasisPtr                    = 0;

    matrixFreeSpinBlock           = 0;
    spinBlockCount                = 0;
    sdSize                        = 0;
    localSpinBlockSize            = 0;
    eOpPtr                        = 0;
    verboseFlag                   = false;
}


void MultiParticleHamiltonianMatrix::createDensityMatrix(int spinBlockIndex, 
UCLAQ::DoubleVector1d& eigVector, UCLAQ::DoubleVector2d& densityMatrix)
{
	long i; long j; long k; long m; 
    long a; long b; 
    

	long rowSize;
    rowSize = sdBasisPtr->getSpinBlockSize(spinBlockIndex);
	
	SlaterDeterminantReference sI;
    SlaterDeterminantReference sJ;

    long IblockIndex;

    long N = sdBasisPtr->getElectronCount();
    long M = sdBasisPtr->getSpinOrbitalCount();
    
    double aI;
    double aJ;

    long  orbitalTmp;

    int aSpin; int bSpin;

    long  basisIndex; 
    long IspinBlockLocalIndex;
    long JspinBlockLocalIndex;

    long IspatialOrbIndex;
    long JspatialOrbIndex;
    double abSign;
 

//
// Although the following is a contiguous loop over all elements of the basis, it
// is effectively a sequential loop over the variables comprising each spin block.   
//
// Here i is the global spin block index.
//
// The original(std) basis index associated with this global spin block index is
// obtained using the call sdBasisPtr->spinBlockToStd[i].
//  
// The spin block and the local spin block index associated with this index is
// obtained using the call sdBasisPtr->getSpinBlockLocalIndex(i,IspinBlock,IspinBlockLocalIndex)
//
//
    
	for(IspinBlockLocalIndex = 0; IspinBlockLocalIndex < rowSize; IspinBlockLocalIndex++)
	{
//
//  Diagonal element 
//
    i            = sdBasisPtr->getSpinBlockGlobalIndex(spinBlockIndex,IspinBlockLocalIndex);
    IblockIndex  = sdBasisPtr->spinBlockToStd[i];
    sI.initialize(IblockIndex, sdBasisPtr);
    
    aI           = eigVector(IspinBlockLocalIndex);
    
	for(m = 0; m < N; m++)
	{
        IspatialOrbIndex  = sI.getSpatialOrbitalIndex(m);
        densityMatrix(IspatialOrbIndex,IspatialOrbIndex) += aI*aI;
	}
//
//  Off diagonal elements
//

// Create list of orbitals

    for(k = 0; k < N; k++)
    {
    orbitalIndexList[k] = sI.getSpinOrbitalIndex(k);
    }

    // Create complimentary list of orbitals

    spinOrbitalUtilities.createComplementarySet(orbitalIndexList,
    complimentaryOrbitalIndexList, N,M); 

    //
    // Loop over all slater determinants that differ from sI by exactly one
    // orbital. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each orbital index, obtaining a new slater determinant by replacing that
    // orbital with one of the elements in the complementary list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //
    // When replacing the orbital, one only needs to consider those cases where the
    // inserted orbital has the same spin; otherwise the total spin will be 
    // different and the inner product of the spins will lead to a zero result. 
    // 
    for(m = 0; m < N; m++)
    {
    orbitalTmp        = orbitalIndexList[m];
    aSpin             = sdBasisPtr->spinOrbitals[orbitalTmp].getSpin();
	a                 = m;
	IspatialOrbIndex  = sI.getSpatialOrbitalIndex(a);
	        
    for(j = 0; j < (M-N); j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    bSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();

    if(aSpin == bSpin)
    {
        basisIndex  = sdBasisPtr->getBasisIndexFromSpinOrbitals(orbitalIndexList);
        sJ.initialize(basisIndex, sdBasisPtr);
        JspinBlockLocalIndex  = sdBasisPtr->stdToSpinBlockLocalIndex[basisIndex];

        aJ  = eigVector(JspinBlockLocalIndex);

        // find orbital indices a of sI and b of sJ of the differing orbital 

        spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],b);
		JspatialOrbIndex  = sJ.getSpatialOrbitalIndex(b);
		
		abSign = 1.0;               
        if(((a+b)%2) != 0 ) abSign = -1.0; 
        densityMatrix(IspatialOrbIndex,JspatialOrbIndex) += abSign*aI*aJ;
    }
    }
    orbitalIndexList[m] = orbitalTmp;
    }

    } // loop over rows in spin block
    
}

void MultiParticleHamiltonianMatrix::initializeSparseSpinBlockMatrix(int spinBlockIndex, 
ElectronOperators& eOp, double dropTolerance,SparseSpinBlockMatrix<UCLAQ::DoubleVector1d>& sparseSpinBlockMatrix)
{
    long i; long j; long k; long m; long n;
    long a; long b; long c; long d;


    spinBlockCount        = sdBasisPtr->getElectronCount() + 1;
    
    string errMessage;
    if(spinBlockIndex >= spinBlockCount) 
    {
    errMessage  = "\nXXX MultiParticleHamiltonianMatrix XXX\n";
    errMessage += "Sparse Matrix Construction Error \n";
    errMessage += "Index of spin block out of range \n";
    throw std::runtime_error(errMessage);
    }
    
    long rowSize;
    rowSize = sdBasisPtr->getSpinBlockSize(spinBlockIndex);
    sparseSpinBlockMatrix.initialize(rowSize);
 
    SlaterDeterminantReference sI;
    SlaterDeterminantReference sJ;

    long IblockIndex;

    long N = sdBasisPtr->getElectronCount();
    long M = sdBasisPtr->getSpinOrbitalCount();

    long  orbitalTmp;
    long  orbitalTmp1;
    long  orbitalTmp2;

    int aSpin; int bSpin;
    int cSpin; int dSpin;

    long  basisIndex; 
    long IspinBlockLocalIndex;
    long JspinBlockLocalIndex;


    if(verboseFlag)
    {
    printf("\nSparse Matrix Construction Start \n");
    }

    ClockIt clockTime;
    ClockIt countTime;

    clockTime.start();
    countTime.start();

    //
    // Execute the construction loop twice. 
    //
    // loopCount = 0: Count the non-zero elements to get sizes.
    //
    // loopCount = 1: Evaluate and pack the coefficients.
    //
    //
    int loopCount;
    double elementValue;
    for(loopCount = 0; loopCount <= 1; loopCount++)
    {
    
    //
    // if loopCount == 1 then create coefficient arrays
    //
    if(loopCount == 1)
    {
        sparseSpinBlockMatrix.createCoefficientArrays();
    }

    //
    // Although the following is a contiguous loop over all elements of the basis, it
    // is effectively a sequential loop over the variables comprising each spin block.   
    //
    // Here i is the global spin block index.
    //
    // The original(std) basis index associated with this global spin block index is
    // obtained using the call sdBasisPtr->spinBlockToStd[i].
    //  
    // The spin block and the local spin block index assocated with this index is 
    // obtained using the call sdBasisPtr->getSpinBlockLocalIndex(i,IspinBlock,IspinBlockLocalIndex)
    //
    //

    for(IspinBlockLocalIndex = 0; IspinBlockLocalIndex < rowSize; IspinBlockLocalIndex++)
    {
    //
    //  Diagonal element 
    //
    i            = sdBasisPtr->getSpinBlockGlobalIndex(spinBlockIndex,IspinBlockLocalIndex);
    IblockIndex  = sdBasisPtr->spinBlockToStd[i];
    //sI           = sdBasisPtr->operator()(IblockIndex);
    sI.initialize(IblockIndex, sdBasisPtr);

    if(loopCount == 0)
    {
    sparseSpinBlockMatrix.incrementRowCount(IspinBlockLocalIndex);
    }
    else
    {
    elementValue = eOp.oneElectronDiagOp(sI) +  eOp.twoElectronDiagOp(sI);
    sparseSpinBlockMatrix.packCoefficientData
                                     (IspinBlockLocalIndex,IspinBlockLocalIndex,elementValue);
    }
  
    //
    //  Off diagonal elements
    //

    // Create list of orbitals

    for(k = 0; k < N; k++)
    {
    orbitalIndexList[k] = sI.getSpinOrbitalIndex(k);
    }

    // Create complimentary list of orbitals

    spinOrbitalUtilities.createComplementarySet(orbitalIndexList,
    complimentaryOrbitalIndexList, N,M); 

    //
    // Loop over all slater determinants that differ from sI by exactly one
    // orbital. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each orbital index, obtaining a new slater determinant by replacing that
    // orbital with one of the elements in the complementary list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //
    // When replacing the orbital, one only needs to consider those cases where the
    // inserted orbital has the same spin; otherwise the total spin will be 
    // different and the inner product of the spins will lead to a zero result. 
     // 


    for(m = 0; m < N; m++)
    {
    orbitalTmp = orbitalIndexList[m];
    aSpin      = sdBasisPtr->spinOrbitals[orbitalTmp].getSpin();

    for(j = 0; j < (M-N); j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    bSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();

    if(aSpin == bSpin)
    {
        //basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
        basisIndex          = sdBasisPtr->getBasisIndexFromSpinOrbitals(orbitalIndexList);
        
        //sJ                  = sdBasisPtr->operator()(basisIndex);
        sJ.initialize(basisIndex, sdBasisPtr);
         
        //JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
        //sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);
        JspinBlockLocalIndex =  sdBasisPtr->stdToSpinBlockLocalIndex[basisIndex];

        // find orbital indices a of sI and b of sJ of the differing orbital 

        a = m;
        spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],b);

        if(loopCount == 0)
        {
            elementValue = eOp.oneElectronOffDiagOp(sI, a, sJ, b) 
                         + eOp.twoElectronOffDiagOp(sI, a, sJ, b);
            if(fabs(elementValue) > dropTolerance)
            {
            sparseSpinBlockMatrix.incrementRowCount(IspinBlockLocalIndex);
            }
        }
        else
        {
            elementValue = eOp.oneElectronOffDiagOp(sI, a, sJ, b) 
                         + eOp.twoElectronOffDiagOp(sI, a, sJ, b);
            if(fabs(elementValue) > dropTolerance)
            {
             sparseSpinBlockMatrix.packCoefficientData
                                     (IspinBlockLocalIndex,JspinBlockLocalIndex,elementValue);
            }
        }
     }

    }
    orbitalIndexList[m] = orbitalTmp;
    }

    

    //
    // Loop over all slater determinants that differ from sI by exactly two
    // orbitals. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each pair of orbital indicies, obtaining a new slater determinant 
    // by replacing that pair with a pair of elements from the complementary 
    // list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //


    for(m = 0; m < N-1; m++)
    {
    orbitalTmp1 = orbitalIndexList[m];
    aSpin      = sdBasisPtr->spinOrbitals[orbitalTmp1].getSpin();

    for(n = m+1; n < N; n++)
    {
    orbitalTmp2 = orbitalIndexList[n];
    bSpin      = sdBasisPtr->spinOrbitals[orbitalTmp2].getSpin();


    for(j = 0; j < (M-N)-1; j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    cSpin      = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();
    for(k = j+1; k < (M-N); k++)
    {
    orbitalIndexList[n] = complimentaryOrbitalIndexList[k];
    dSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[n]].getSpin();

    if((aSpin + bSpin)  == (cSpin + dSpin)) // insures we are in the same spin block
    {
    //basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
    basisIndex          = sdBasisPtr->getBasisIndexFromSpinOrbitals(orbitalIndexList);
    
    //sJ                  = sdBasisPtr->operator()(basisIndex);
    sJ.initialize(basisIndex, sdBasisPtr);
    
    // 
    // find indices of orbitals that differ between sI and sJ 
    //

    // Indices for sI (by construction)

    a = m;  b = n;

    // Find indices for sJ 

    spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],
                                                orbitalIndexList[n],c,d);

    
    //JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
    //sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);
    JspinBlockLocalIndex =  sdBasisPtr->stdToSpinBlockLocalIndex[basisIndex];

    if(loopCount == 0)
    {
        elementValue = eOp.twoElectronOffDiagOp(sI, a, b, sJ, c, d);
        if(fabs(elementValue) > dropTolerance)
        {
           sparseSpinBlockMatrix.incrementRowCount(IspinBlockLocalIndex);
        }
    }
    else
    {
        elementValue = eOp.twoElectronOffDiagOp(sI, a, b, sJ, c, d);
        if(fabs(elementValue) > dropTolerance)
        {
        sparseSpinBlockMatrix.packCoefficientData
                                     (IspinBlockLocalIndex,JspinBlockLocalIndex,elementValue);
        }
    }

    }

    }}

    orbitalIndexList[m] = orbitalTmp1;
    orbitalIndexList[n] = orbitalTmp2;
    }}

    }

	if(loopCount == 0)
	{countTime.stop();}
 
    } // end of loop count

    sparseSpinBlockMatrix.sortRows();

    if(verboseFlag)
    {
	printf("Sparse Matrix Count Time (minutes)         : %-10.6f \n",
	(countTime.getMilliSecElapsedTime())*(1.0/60000.0));
    
    clockTime.stop();
	printf("Sparse Matrix Construction Time (minutes)  : %-10.6f \n",
    (clockTime.getMilliSecElapsedTime())*(1.0/60000.0));
    }

    //
    // now sort the indices of each row and print out statistics

    double rowSizeD;
    rowSizeD = (double)sdBasisPtr->getSpinBlockSize(spinBlockIndex);

    if(verboseFlag)
    {
    printf("Non-zero Spin Block Elements = %10ld out of %10.f   Percentage = %10.2f \n",
    sparseSpinBlockMatrix.getTotalDataCount(),rowSizeD*rowSizeD,
        100.0*double(sparseSpinBlockMatrix.getTotalDataCount())/double(rowSizeD*rowSizeD));
    }
}



MultiParticleHamiltonianMatrix::MultiParticleHamiltonianMatrix(SlaterDeterminantBasis& sdB) 
    : vTmp()
{
	sdBasisPtr          = &sdB;
    sdSize              = sdBasisPtr->getBasisCount();
    spinBlockCount      = sdBasisPtr->getElectronCount() + 1;
    spinBlockMatrix     = 0;

    matrixFreeSpinBlock = 0;

    long N = sdBasisPtr->getElectronCount();
    long M = sdBasisPtr->getSpinOrbitalCount();

    spinOrbitalUtilities.initialize(M,N);

    orbitalIndexList              = new long[N];
    complimentaryOrbitalIndexList = new long[M-N];

    sparseSpinBlockMatrix = 0;
    verboseFlag           = false;
    localSpinBlockSize    = 0;
    eOpPtr                = 0;

}


MultiParticleHamiltonianMatrix::~MultiParticleHamiltonianMatrix()
{destroyData();}


void MultiParticleHamiltonianMatrix::destroyData()
{
    if(spinBlockMatrix != 0) 
    {delete [] spinBlockMatrix;}

    if(orbitalIndexList != 0) 
    {delete [] orbitalIndexList;}

    if(complimentaryOrbitalIndexList != 0) 
    {delete [] complimentaryOrbitalIndexList;}

    if(sparseSpinBlockMatrix != 0)
    {delete [] sparseSpinBlockMatrix;}


    spinBlockMatrix               = 0;
	orbitalIndexList              = 0;
    complimentaryOrbitalIndexList = 0;
    sparseSpinBlockMatrix         = 0;
}


void MultiParticleHamiltonianMatrix::initializeSpinBlockMatrices_FullLoops()
{
    long i; long j; long k;

    SlaterDeterminantReference sI;
    SlaterDeterminantReference sJ;


    long IblockIndex;
    long JblockIndex;

    long IspinBlockGlobalIndex;
    long JspinBlockGlobalIndex;

    long *IorbitalList;
    long *JorbitalList;

    long diffCount;
    long electronCount = sdBasisPtr->operator()(0).getElectronCount();
    long nonZeroElementCount;
    double rowSizeD;
    

    ClockIt countTime;

    if(verboseFlag)
    {
    printf("\nAlternate Matrix Construction Start \n");
    }

    countTime.start();
    
    for(k = 0; k < spinBlockCount; k++)
    {
    nonZeroElementCount = 0;
    localSpinBlockSize  = sdBasisPtr->getSpinBlockSize(k);

    for(i = 0; i < localSpinBlockSize; i++)              
    {
    IspinBlockGlobalIndex = sdBasisPtr->getSpinBlockGlobalIndex(k,i);// global spin block index
    IblockIndex  = sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex]; // std. basis index
    sI           = sdBasisPtr->operator()(IblockIndex);                          // ith element of basis
    IorbitalList = sI.getSpinOrbitalIndexList();

    for(j = 0; j < localSpinBlockSize; j++)             
    {
    JspinBlockGlobalIndex = sdBasisPtr->getSpinBlockGlobalIndex(k,j);// global spin block index
    JblockIndex = sdBasisPtr->spinBlockToStd[JspinBlockGlobalIndex]; // std. basis index
    sJ = sdBasisPtr->operator()(JblockIndex);                                   // jth element of basis
    JorbitalList = sJ.getSpinOrbitalIndexList();
    //
    // 
    diffCount = spinOrbitalUtilities.getOrbitalDifferenceCount(IorbitalList,
        JorbitalList, electronCount);
    if(diffCount <= 2)
    {
    nonZeroElementCount += 1;
    //spinBlockMatrix[k](i,j) += SingleParticleOperator(sI,sJ);
    //spinBlockMatrix[k](i,j) += DoubleParticleOperator(sI,sJ);
    }
    }}
    
	rowSizeD = (double)localSpinBlockSize;
	if(verboseFlag)
	{
	printf("Non-zero Spin Block Elements = %10ld out of %10.f   Percentage = %10.2f \n",
	nonZeroElementCount,rowSizeD*rowSizeD,
	100.0*double(nonZeroElementCount)/double(rowSizeD*rowSizeD));
	}
    
    }
    countTime.stop();

    if(verboseFlag)
    {
    printf("Alternate Matrix Count Time (minutes)         : %-10.6f \n",
	(countTime.getMilliSecElapsedTime())*(1.0/60000.0));
	}
}


// OpType =  ONE_PTCLE_OP or TWO_PTCLE_OP or ONE_AND_TWO_PTCLE_OP

void MultiParticleHamiltonianMatrix::initializeSparseSpinBlockMatrices
                                  (ElectronOperators& eOp, double dropTolerance)
{
    long i; long j; long k; long m; long n;
    long a; long b; long c; long d;

    //
    // Initialize array of matrices for sparse spin blocks
    //

    spinBlockCount        = sdBasisPtr->getElectronCount() + 1;
    sparseSpinBlockMatrix = new SparseSpinBlockMatrix<UCLAQ::DoubleVector1d>[spinBlockCount];
    long rowSize;

    for(k = 0; k < spinBlockCount; k++)
    {
        rowSize = sdBasisPtr->getSpinBlockSize(k);
        sparseSpinBlockMatrix[k].initialize(rowSize);
    }

    SlaterDeterminantReference sI;
    SlaterDeterminantReference sJ;

    long IblockIndex;
    long JblockIndex;

    long N = sdBasisPtr->getElectronCount();
    long M = sdBasisPtr->getSpinOrbitalCount();

    long  orbitalTmp;
    long  orbitalTmp1;
    long  orbitalTmp2;

    int aSpin; int bSpin;
    int cSpin; int dSpin;

    long  basisIndex; 

    long IspinBlock;
    long IspinBlockLocalIndex;

    long JspinBlock;
    long JspinBlockLocalIndex;

    ClockIt clockTime;
    ClockIt countTime;
    
    clockTime.start();
    countTime.start();

    if(verboseFlag)
    {
    printf("\nSparse Matrix Construction Start \n");
    }

    //
    // Execute the construction loop twice. 
    //
    // loopCount = 0: Count the non-zero elements to get sizes. 
    //
    // loopCount = 1: Evaluate and pack the coefficients.
    //
    //
    int loopCount;
    double elementValue;
    for(loopCount = 0; loopCount <= 1; loopCount++)
    {
    
    //
    // if loopCount == 1 then create coefficient arrays
    //
    if(loopCount == 1)
    {
        for(k = 0; k < spinBlockCount; k++)
        {
        sparseSpinBlockMatrix[k].createCoefficientArrays();
        }
    }

    //
    // Although the following is a contiguous loop over all elements of the basis, it
    // is effectively a sequential loop over the variables comprising each spin block.   
    //
    // Here i is the global spin block index.
    //
    // The original(std) basis index associated with this global spin block index is
    // obtained using the call sdBasisPtr->spinBlockToStd[i].
    //  
    // The spin block and the local spin block index assocated with this index is 
    // obtained using the call sdBasisPtr->getSpinBlockLocalIndex(i,IspinBlock,IspinBlockLocalIndex)
    //
    //
    for(i = 0; i < sdSize; i++)
    {
    //
    //  Diagonal element 
    //
    IblockIndex  = sdBasisPtr->spinBlockToStd[i];
    sI           = sdBasisPtr->operator()(IblockIndex);
    sdBasisPtr->getSpinBlockLocalIndex(i,IspinBlock,IspinBlockLocalIndex);

    if(loopCount == 0)
    {
    sparseSpinBlockMatrix[IspinBlock].incrementRowCount(IspinBlockLocalIndex);
    }
    else
    {
    elementValue = eOp.oneElectronDiagOp(sI) +  eOp.twoElectronDiagOp(sI);
    sparseSpinBlockMatrix[IspinBlock].packCoefficientData
                                     (IspinBlockLocalIndex,IspinBlockLocalIndex,elementValue);
    }
  
    //
    //  Off diagonal elements
    //

    // Create list of orbitals

    for(k = 0; k < N; k++)
    {
    orbitalIndexList[k] = sI.getSpinOrbitalIndex(k);
    }

    // Create complimentary list of orbitals

    spinOrbitalUtilities.createComplementarySet(orbitalIndexList,
    complimentaryOrbitalIndexList, N,M); 

    //
    // Loop over all slater determinants that differ from sI by exactly one
    // orbital. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each orbital index, obtaining a new slater determinant by replacing that
    // orbital with one of the elements in the complementary list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //
    // When replacing the orbital, one only needs to consider those cases where the
    // inserted orbital has the same spin; otherwise the total spin will be 
    // different and the inner product of the spins will lead to a zero result. 
     // 


    for(m = 0; m < N; m++)
    {
    orbitalTmp = orbitalIndexList[m];
    aSpin      = sdBasisPtr->spinOrbitals[orbitalTmp].getSpin();

    for(j = 0; j < (M-N); j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    bSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();

    if(aSpin == bSpin)
    {
        basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
        sJ                  = sdBasisPtr->operator()(basisIndex);
        JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
        sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);

        // find orbital indices a of sI and b of sJ of the differing orbital 

        a = m;
        spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],b);

        if(loopCount == 0)
        {
            elementValue = eOp.oneElectronOffDiagOp(sI, a, sJ, b) 
                         + eOp.twoElectronOffDiagOp(sI, a, sJ, b);
            if(fabs(elementValue) > dropTolerance)
            {
            sparseSpinBlockMatrix[IspinBlock].incrementRowCount(IspinBlockLocalIndex);
            }
        }
        else
        {
            elementValue = eOp.oneElectronOffDiagOp(sI, a, sJ, b) 
                         + eOp.twoElectronOffDiagOp(sI, a, sJ, b);
            if(fabs(elementValue) > dropTolerance)
            {
             sparseSpinBlockMatrix[IspinBlock].packCoefficientData
                                     (IspinBlockLocalIndex,JspinBlockLocalIndex,elementValue);
            }
        }
     }

    }
    orbitalIndexList[m] = orbitalTmp;
    }

    

    //
    // Loop over all slater determinants that differ from sI by exactly two
    // orbitals. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each pair of orbital indicies, obtaining a new slater determinant 
    // by replacing that pair with a pair of elements from the complementary 
    // list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //


    for(m = 0; m < N-1; m++)
    {
    orbitalTmp1 = orbitalIndexList[m];
    aSpin      = sdBasisPtr->spinOrbitals[orbitalTmp1].getSpin();

    for(n = m+1; n < N; n++)
    {
    orbitalTmp2 = orbitalIndexList[n];
    bSpin      = sdBasisPtr->spinOrbitals[orbitalTmp2].getSpin();


    for(j = 0; j < (M-N)-1; j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    cSpin      = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();
    for(k = j+1; k < (M-N); k++)
    {
    orbitalIndexList[n] = complimentaryOrbitalIndexList[k];
    dSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[n]].getSpin();

    if((aSpin + bSpin)  == (cSpin + dSpin)) // insures we are in the same spin block
    {
    basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
    sJ                  = sdBasisPtr->operator()(basisIndex);
    
    // 
    // find indices of orbitals that differ between sI and sJ 
    //

    // Indices for sI (by construction)

    a = m;  b = n;

    // Find indices for sJ 

    spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],
                                                orbitalIndexList[n],c,d);

    
    JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
    sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);
   

    if(loopCount == 0)
    {
        elementValue = eOp.twoElectronOffDiagOp(sI, a, b, sJ, c, d);
        if(fabs(elementValue) > dropTolerance)
        {
           sparseSpinBlockMatrix[IspinBlock].incrementRowCount(IspinBlockLocalIndex);
        }
    }
    else
    {
        elementValue = eOp.twoElectronOffDiagOp(sI, a, b, sJ, c, d);
        if(fabs(elementValue) > dropTolerance)
        {
        sparseSpinBlockMatrix[IspinBlock].packCoefficientData
                                     (IspinBlockLocalIndex,JspinBlockLocalIndex,elementValue);
        }
    }

    }

    }}

    orbitalIndexList[m] = orbitalTmp1;
    orbitalIndexList[n] = orbitalTmp2;
    }}

    }

	if(loopCount == 0)
	{countTime.stop();}
 
    } // end of loop count

   for(k = 0; k < spinBlockCount; k++)
   {
    sparseSpinBlockMatrix[k].sortRows();
   }


    clockTime.stop();
    
    if(verboseFlag)
    {
	printf("Sparse Matrix Count Time (minutes)         : %-10.6f \n",
	(countTime.getMilliSecElapsedTime())*(1.0/60000.0));
    
	printf("Sparse Matrix Construction Time (minutes)  : %-10.6f \n",
    (clockTime.getMilliSecElapsedTime())*(1.0/60000.0));
    }

    //
    // now sort the indices of each row 
    // and print out statistics 

    double rowSizeD;

    if(verboseFlag)
    {
    for(k = 0; k < spinBlockCount; k++)
    {
        rowSizeD = (double)sdBasisPtr->getSpinBlockSize(k);

        printf("Non-zero Spin Block Elements = %10ld out of %10.f   Percentage = %10.2f \n",
        sparseSpinBlockMatrix[k].getTotalDataCount(),rowSizeD*rowSizeD,
        100.0*double(sparseSpinBlockMatrix[k].getTotalDataCount())/double(rowSizeD*rowSizeD));
    }
    }
}


void MultiParticleHamiltonianMatrix::applyForwardOp(UCLAQ::DoubleVector1d& V)
{

    long i; long j; long k; long m; long n;
    long a; long b; long c; long d;

    long spinBlockSize = sdBasisPtr->getSpinBlockSize(matrixFreeSpinBlock);

    vTmp.initialize(spinBlockSize);
     
    SlaterDeterminantReference sI;
    SlaterDeterminantReference sJ;

    long IblockIndex = matrixFreeSpinBlock;
    long JblockIndex;

    long JspinBlock;
    long JspinBlockLocalIndex;


    long N = sdBasisPtr->getElectronCount();
    long M = sdBasisPtr->getSpinOrbitalCount();

    long  orbitalTmp;
    long  orbitalTmp1;
    long  orbitalTmp2;

    int aSpin; int bSpin;
    int cSpin; int dSpin;

    long  basisIndex; 

    long iStdIndex;
    long iSpinBlockGlobal;

    // loop over each row of the local spin block 

    for(i = 0; i < spinBlockSize; i++)
    {
    vTmp(i) = 0.0;

    // determine ith basis element 

    iSpinBlockGlobal = sdBasisPtr->getSpinBlockGlobalIndex(IblockIndex,i);
    iStdIndex        = sdBasisPtr->spinBlockToStd[iSpinBlockGlobal];
    sI               = sdBasisPtr->operator()(iStdIndex);

    //
    //  Diagonal accumulation 
    //
    vTmp(i) += V(i)*(eOpPtr->oneElectronDiagOp(sI) + eOpPtr->twoElectronDiagOp(sI));

    // Create list of orbitals

    for(k = 0; k < N; k++)
    {
    orbitalIndexList[k] = sI.getSpinOrbitalIndex(k);
    }

    // Create complimentary list of orbitals

    spinOrbitalUtilities.createComplementarySet(orbitalIndexList,
    complimentaryOrbitalIndexList, N, M); 

    //
    // Loop over all slater determinants that differ from sI by exactly one
    // orbital. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each orbital index, obtaining a new slater determinant by replacing that
    // orbital with one of the elements in the complementary list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //
    // When replacing the orbital, one only needs to consider those cases where the
    // inserted orbital has the same spin; otherwise the total spin will be 
    // different and the inner product of the spins will lead to a zero result. 
    // 

    for(m = 0; m < N; m++)
    {
    orbitalTmp = orbitalIndexList[m];
    aSpin      = sdBasisPtr->spinOrbitals[orbitalTmp].getSpin();

    for(j = 0; j < (M-N); j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    bSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();

    if(aSpin == bSpin) // this insures the same spin block 
    {
        basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
        sJ                  = sdBasisPtr->operator()(basisIndex);
        JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
        sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);

        // find orbital indices a of sI and b of sJ of the differing orbital 

        a = m;
        spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],b);

        vTmp(i) += V(JspinBlockLocalIndex)*(  eOpPtr->oneElectronOffDiagOp(sI, a, sJ, b) 
                                            + eOpPtr->twoElectronOffDiagOp(sI, a, sJ, b));
    }}

    orbitalIndexList[m] = orbitalTmp;
    }

    //
    // Loop over all slater determinants that differ from sI by exactly two
    // orbitals. 
    //

    
    for(m = 0; m < N-1; m++)
    {
    orbitalTmp1 = orbitalIndexList[m];
    aSpin      = sdBasisPtr->spinOrbitals[orbitalTmp1].getSpin();

    for(n = m+1; n < N; n++)
    {
    orbitalTmp2 = orbitalIndexList[n];
    bSpin      = sdBasisPtr->spinOrbitals[orbitalTmp2].getSpin();


    for(j = 0; j < (M-N)-1; j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    cSpin      = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();
    for(k = j+1; k < (M-N); k++)
    {
    orbitalIndexList[n] = complimentaryOrbitalIndexList[k];
    dSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[n]].getSpin();

    if((aSpin + bSpin)  == (cSpin + dSpin)) // insures we are in the same spin block 
    {
    basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
    sJ                  = sdBasisPtr->operator()(basisIndex);
    
    // 
    // find indices of orbitals that differ between sI and sJ 
    //

    // Indices for sI (by construction)

    a = m;  b = n;

    // Find indices for sJ 

    spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],
                                                orbitalIndexList[n],c,d);

    
    JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
    sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);
    
   
    vTmp(i) += V(JspinBlockLocalIndex)*eOpPtr->twoElectronOffDiagOp(sI, a, b, sJ, c, d);

    }

    }}

    orbitalIndexList[m] = orbitalTmp1;
    orbitalIndexList[n] = orbitalTmp2;
    }}
   

    } // loop over spin block rows

    V = vTmp;
}


// OpType =  ONE_PTCLE_OP or TWO_PTCLE_OP or ONE_AND_TWO_PTCLE_OP

void MultiParticleHamiltonianMatrix::initializeSpinBlockMatrices(ElectronOperators& eOp, int opType)
{

    long i; long j; long k; long m; long n;
    long a; long b; long c; long d;

    //
    // Allocate array of matrices for spin blocks
    //

    spinBlockCount  = sdBasisPtr->getElectronCount() + 1;

    if(spinBlockMatrix != 0) {delete [] spinBlockMatrix;}
    spinBlockMatrix = new UCLAQ::DoubleVector2d[spinBlockCount];

    double* spinBlockDataPtr;

    for(k = 0; k < spinBlockCount; k++)
    {
        localSpinBlockSize = sdBasisPtr->getSpinBlockSize(k);
        spinBlockMatrix[k].initialize(localSpinBlockSize,localSpinBlockSize);
        spinBlockDataPtr = spinBlockMatrix[k].getDataPointer();

        for(i = 0; i < localSpinBlockSize*localSpinBlockSize; i++)
        {
        spinBlockDataPtr[i] = 0;
        }
    }

    long elementCount = 0;

    SlaterDeterminantReference sI;
    SlaterDeterminantReference sJ;

    long IblockIndex;

    long N = sdBasisPtr->getElectronCount();
    long M = sdBasisPtr->getSpinOrbitalCount();

    long  orbitalTmp;
    long  orbitalTmp1;
    long  orbitalTmp2;

    int aSpin; int bSpin;
    int cSpin; int dSpin;

    long  basisIndex; 

    long IspinBlock;
    long IspinBlockLocalIndex;
    long JspinBlockLocalIndex;

    ClockIt clockTime;
    clockTime.start();


    if(verboseFlag)
    {
    printf("\nDiagonal Matrix Construction Start \n");
    clockTime.start();
    }

    for(i = 0; i < sdSize; i++)
    {
    //
    //  Diagonal element 
    //
    IblockIndex  = sdBasisPtr->spinBlockToStd[i];
    //sI           = sdBasisPtr->operator()(IblockIndex);
    sI.initialize(IblockIndex, sdBasisPtr);

    sdBasisPtr->getSpinBlockLocalIndex(i,IspinBlock,IspinBlockLocalIndex);
    

    elementCount++;

    if(     opType == ONE_PTCLE_OP)
    {
    spinBlockMatrix[IspinBlock](IspinBlockLocalIndex,IspinBlockLocalIndex)
    += eOp.oneElectronDiagOp(sI);
    }
    else if(opType == TWO_PTCLE_OP)
    {
    spinBlockMatrix[IspinBlock](IspinBlockLocalIndex,IspinBlockLocalIndex)
    += eOp.twoElectronDiagOp(sI);
    }
    else // ONE_AND_TWO_PTCLE_OP
    {
    spinBlockMatrix[IspinBlock](IspinBlockLocalIndex,IspinBlockLocalIndex)
    += eOp.oneElectronDiagOp(sI) +  eOp.twoElectronDiagOp(sI);
    }
    }

    if(verboseFlag)
    {
    clockTime.stop();
    printf("Diagonal Matrix Construction End : Time (minutes)  : %-10.6f \n",
    (clockTime.getMilliSecElapsedTime())*(1.0/60000.0));
    }

    
    if(verboseFlag)
    {
    printf("\nOne Electron Matrix Construction Start \n");
    clockTime.start();;
    }

    //
    //  Off diagonal elements
    //
    for(i = 0; i < sdSize; i++)
    {
    IblockIndex  = sdBasisPtr->spinBlockToStd[i];
    //sI           = sdBasisPtr->operator()(IblockIndex);
    sI.initialize(IblockIndex, sdBasisPtr);
    sdBasisPtr->getSpinBlockLocalIndex(i,IspinBlock,IspinBlockLocalIndex);
   
    // Create the list of orbitals in sI

    for(k = 0; k < N; k++)
    {
    orbitalIndexList[k] = sI.getSpinOrbitalIndex(k);
    }

    // Create complimentary list of orbitals

    spinOrbitalUtilities.createComplementarySet(orbitalIndexList,
    complimentaryOrbitalIndexList, N,M); 

    //
    // Loop over all slater determinants that differ from sI by exactly one
    // orbital. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each orbital index, obtaining a new slater determinant by replacing that
    // orbital with one of the elements in the complementary list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //
    // When replacing the orbital, one only needs to consider those cases where the
    // inserted orbital has the same spin; otherwise the total spin will be 
    // different and the inner product of the spins will lead to a zero result. 
     // 


    for(m = 0; m < N; m++)
    {
    orbitalTmp = orbitalIndexList[m];
    aSpin      = sdBasisPtr->spinOrbitals[orbitalTmp].getSpin();

    for(j = 0; j < (M-N); j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    bSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();

    if(aSpin == bSpin)
    {
        elementCount++;
        //basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
        basisIndex          = sdBasisPtr->getBasisIndexFromSpinOrbitals(orbitalIndexList);
        //sJ                  = sdBasisPtr->operator()(basisIndex);
        sJ.initialize(basisIndex, sdBasisPtr);
        //JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
        //sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);
        JspinBlockLocalIndex =  sdBasisPtr->stdToSpinBlockLocalIndex[basisIndex];

        // find orbital indices a of sI and b of sJ of the differing orbital 

        a = m;
        spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],b);

        if(     opType == ONE_PTCLE_OP)
        {
        spinBlockMatrix[IspinBlock](IspinBlockLocalIndex,JspinBlockLocalIndex) 
        += eOp.oneElectronOffDiagOp(sI, a, sJ, b);
        }
        else if (opType == TWO_PTCLE_OP)
        {
        spinBlockMatrix[IspinBlock](IspinBlockLocalIndex,JspinBlockLocalIndex) 
        += eOp.twoElectronOffDiagOp(sI, a, sJ, b);
        }
        else // ONE_AND_TWO_PTCLE_OP
        {
        spinBlockMatrix[IspinBlock](IspinBlockLocalIndex,JspinBlockLocalIndex) 
        += eOp.oneElectronOffDiagOp(sI, a, sJ, b) 
         + eOp.twoElectronOffDiagOp(sI, a, sJ, b);
        }
    }}

    orbitalIndexList[m] = orbitalTmp;
    }

    }

    if(verboseFlag)
    {
    clockTime.stop();
    printf("One Electron Matrix Construction End : Time (minutes)  : %-10.6f \n",
    (clockTime.getMilliSecElapsedTime())*(1.0/60000.0));
    }

    //printf("\n        Returning early  \n");
    //printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    //printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    //return;
    
    if(verboseFlag)
    {
    printf("\nTwo Electron Matrix Construction Start \n");
    clockTime.start();
    }


    for(i = 0; i < sdSize; i++)
    {
    IblockIndex  = sdBasisPtr->spinBlockToStd[i];
    sI           = sdBasisPtr->operator()(IblockIndex);
    sdBasisPtr->getSpinBlockLocalIndex(i,IspinBlock,IspinBlockLocalIndex);

    //
    //  Off diagonal elements
    //

    // Create list of orbitals in sI

    for(k = 0; k < N; k++)
    {
    orbitalIndexList[k] = sI.getSpinOrbitalIndex(k);
    }

    // Create complimentary list of orbitals

    spinOrbitalUtilities.createComplementarySet(orbitalIndexList,
    complimentaryOrbitalIndexList, N,M); 

    //
    // Loop over all slater determinants that differ from sI by exactly two
    // orbitals. 
    //
    // This loop consists of looping over the orbital indices in sI, and, 
    // for each pair of orbital indicies, obtaining a new slater determinant 
    // by replacing that pair with a pair of elements from the complementary 
    // list of orbitals
    //
    // What makes this work is the ability to find the Slater determinant basis set
    // index given a list of orbital indices. 
    //


    for(m = 0; m < N-1; m++)
    {
    orbitalTmp1 = orbitalIndexList[m];
    aSpin       = sdBasisPtr->spinOrbitals[orbitalTmp1].getSpin();

    for(n = m+1; n < N; n++)
    {
    orbitalTmp2 = orbitalIndexList[n];
    bSpin      = sdBasisPtr->spinOrbitals[orbitalTmp2].getSpin();


    for(j = 0; j < (M-N)-1; j++)
    {
    orbitalIndexList[m] = complimentaryOrbitalIndexList[j];
    cSpin      = sdBasisPtr->spinOrbitals[orbitalIndexList[m]].getSpin();
    for(k = j+1; k < (M-N); k++)
    {
    orbitalIndexList[n] = complimentaryOrbitalIndexList[k];
    dSpin               = sdBasisPtr->spinOrbitals[orbitalIndexList[n]].getSpin();

    if((aSpin + bSpin)  == (cSpin + dSpin)) // insures we are in the same spin block
    {
    elementCount++;
    //basisIndex          = spinOrbitalUtilities.getBasisIndex(orbitalIndexList,N,M);
    basisIndex          = sdBasisPtr->getBasisIndexFromSpinOrbitals(orbitalIndexList);
    //sJ                  = sdBasisPtr->operator()(basisIndex);
    sJ.initialize(basisIndex, sdBasisPtr);
    
    // 
    // find indices of orbitals that differ between sI and sJ 
    //

    // Indices for sI (by construction)

    a = m;  b = n;

    // Find indices for sJ 

    spinOrbitalUtilities.getOrbitalIndices(sJ,orbitalIndexList[m],
                                                orbitalIndexList[n],c,d);

    
    //JblockIndex         = sdBasisPtr->stdToSpinBlock[basisIndex];
    //sdBasisPtr->getSpinBlockLocalIndex(JblockIndex,JspinBlock,JspinBlockLocalIndex);
    JspinBlockLocalIndex =  sdBasisPtr->stdToSpinBlockLocalIndex[basisIndex];
    
    if (opType != ONE_PTCLE_OP)
    {
    spinBlockMatrix[IspinBlock](IspinBlockLocalIndex,JspinBlockLocalIndex) 
    += eOp.twoElectronOffDiagOp(sI, a, b, sJ, c, d);
    }
    }

    }}

    orbitalIndexList[m] = orbitalTmp1; // restore original SI index list values
    orbitalIndexList[n] = orbitalTmp2;
    }}

    }


    if(verboseFlag)
    {
    clockTime.stop();
    printf("Two Electron Matrix Construction End : Time (minutes)  : %-10.6f \n",
    (clockTime.getMilliSecElapsedTime())*(1.0/60000.0));

    printf("\nTotal Elements Computed %ld out of %ld \n",elementCount,sdSize*sdSize);
    printf("Sparsity Percentage : %10.5f \n\n",double(elementCount)/double(sdSize*sdSize));
    }
}


void MultiParticleHamiltonianMatrix::outputMatrix()
{
    long i; long j; long k; 

    long N = sdBasisPtr->getElectronCount();

    // long M = sdBasisPtr->getSpinOrbitalCount();

    SlaterDeterminantReference sI;
    long IblockIndex;

    long IspinBlockGlobalIndex;

    for(k = 0; k < spinBlockCount; k++)
    {
    localSpinBlockSize = sdBasisPtr->getSpinBlockSize(k);

    for(i = 0; i < localSpinBlockSize; i++)              
    {
    IspinBlockGlobalIndex = sdBasisPtr->getSpinBlockGlobalIndex(k,i);// global spin block index
    IblockIndex = sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex]; // std. basis index
    sI          = sdBasisPtr->operator()(IblockIndex);                          // ith element of basis

    for(j = 0; j < N; j++)
    {
    if(sI.getOrbitalSpin(j) == 1)
    {
    printf("%2ld+", sI.getSpatialOrbitalIndex(j)+1);
    }
    else
    {
    printf("%2ld-", sI.getSpatialOrbitalIndex(j)+1);
    }
    }

    printf(" : ");
    for(j = 0; j < localSpinBlockSize; j++)
    {
    printf("%6.3e ",spinBlockMatrix[k](i,j));
    }

    printf(" \n");
    }}
}

void MultiParticleHamiltonianMatrix::outputSpinBlockMatrices
(const char* outputFilePrefix, const char* outputFormatString,const char* inputIntegralDataFile)
{
    ostringstream fileName;  
    ostringstream indexFileName;  
    ostringstream outputFormat;
    
    FILE *MatrixOutputFile;
    FILE *MatrixIndexOutputFile;
    
    long i; long j; long k; 

    long N = sdBasisPtr->getElectronCount();

    // long M = sdBasisPtr->getSpinOrbitalCount();

    SlaterDeterminantReference sI;
    long IblockIndex;
    long IspinBlockGlobalIndex;

    cout << endl;
    cout << "Spin Block Matrix Output Files: " << endl << endl;
    
    for(k = 0; k < spinBlockCount; k++)
    {
    localSpinBlockSize = sdBasisPtr->getSpinBlockSize(k);
    //
    // Compose file name based on spin
    //
    i = 0;
    IspinBlockGlobalIndex = sdBasisPtr->getSpinBlockGlobalIndex(k,i);// global spin block index
    IblockIndex = sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex]; // std. basis index
    sI          = sdBasisPtr->operator()(IblockIndex);                          // ith element of basis
    
    fileName.str("");
    indexFileName.str("");
    fileName      << outputFilePrefix << "_";
    indexFileName << outputFilePrefix << "Index_";
    for(j = 0; j < N; j++)
    {
    if(sI.getOrbitalSpin(j) == 1) {fileName << "U"; indexFileName << "U";}
    else                          {fileName << "D"; indexFileName << "D";}
    }

    fileName << ".dat";
    indexFileName << ".dat";

    cout << (fileName.str()).c_str() << "  " << (indexFileName.str()).c_str() << endl;
 
    MatrixOutputFile      = fopen((fileName.str()).c_str(), "w");
    MatrixIndexOutputFile = fopen((indexFileName.str()).c_str(), "w");
   
    localSpinBlockSize = sdBasisPtr->getSpinBlockSize(k);
    
    fprintf(MatrixIndexOutputFile,"################################################################################\n");
    fprintf(MatrixIndexOutputFile,"Indexing for spin block matrix file : %s \n",(fileName.str()).c_str());
    fprintf(MatrixIndexOutputFile,"Orbital integrals data file         : %s \n",inputIntegralDataFile);
    fprintf(MatrixIndexOutputFile,"################################################################################\n\n");
    fprintf(MatrixIndexOutputFile,"ROW : Spin Orbitals \n");
    fprintf(MatrixIndexOutputFile,"=================== \n");
    
    for(i = 0; i < localSpinBlockSize; i++)              
    {
    IspinBlockGlobalIndex = sdBasisPtr->getSpinBlockGlobalIndex(k,i);// global spin block index
    IblockIndex = sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex]; // std. basis index
    sI          = sdBasisPtr->operator()(IblockIndex);                          // ith element of basis

    fprintf(MatrixIndexOutputFile,"%3ld : ",i+1);
    for(j = 0; j < N; j++)
    {
    if(sI.getOrbitalSpin(j) == 1)
    {fprintf(MatrixIndexOutputFile,"%4ld+", sI.getSpatialOrbitalIndex(j)+1);}
    else
    {fprintf(MatrixIndexOutputFile,"%4ld-", sI.getSpatialOrbitalIndex(j)+1);}
    }
    fprintf(MatrixIndexOutputFile,"\n");
    
    outputFormat.str("");
    outputFormat << outputFormatString << " ";
    for(j = 0; j < localSpinBlockSize; j++)
    {
    fprintf(MatrixOutputFile,(outputFormat.str()).c_str(),spinBlockMatrix[k](i,j));
    }
    fprintf(MatrixOutputFile,"\n");
    }
    fclose(MatrixOutputFile);
    fclose(MatrixIndexOutputFile);
    }	
    
    cout << endl;
}
//
// Routine to output spin block matrices in A.K's desired format : 
// dictionary ordering with upSpin orbitals first.
// 

void MultiParticleHamiltonianMatrix::outputSpinBlockMatricesKformat
(const char* outputFilePrefix, const char* outputFormatString,const char* inputIntegralDataFile)
{
    UCLAQ::LongArray2d    indexAndSpin;   // Local storage for given basis element spatial indices and spin
    UCLAQ::LongArray2d    orbitalIndices; // Integer value associated with a given basis element (for sorting)
    UCLAQ::LongArray1d    indexMap;       // ith element stores the index of the original basis element
                                            // mapped to the ith element of the new basis.
    
    
    UCLAQ::DoubleVector1d swapSign;

    ostringstream fileName;  
    ostringstream indexFileName;  
    ostringstream outputFormat;
    
    FILE *MatrixOutputFile;
    FILE *MatrixIndexOutputFile;
    
    long i; long j; long k; 

    long N = sdBasisPtr->getElectronCount();

    // long M = sdBasisPtr->getSpinOrbitalCount();

    SlaterDeterminantReference sI;
    long IblockIndex;
    long IspinBlockGlobalIndex;

    cout << endl;
    cout << "Spin Block Matrix Output Files: " << endl << endl;

    int upSpinCount;
    int downSpinCount;
    int       swapFlag;
    double swapSignTmp;
    int       indexTmp;
    int spinTmp;
   
    
    for(k = 0; k < spinBlockCount; k++)
    {
    localSpinBlockSize = sdBasisPtr->getSpinBlockSize(k);

    //
    // Compose file name based on spin
    //
    i = 0;
    IspinBlockGlobalIndex = sdBasisPtr->getSpinBlockGlobalIndex(k,i);// global spin block index
    IblockIndex = sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex]; // std. basis index
    sI          = sdBasisPtr->operator()(IblockIndex);                          // ith element of basis
    
    fileName.str("");
    indexFileName.str("");
    fileName      << outputFilePrefix << "_";
    indexFileName << outputFilePrefix << "Index_";

    upSpinCount   = 0;
    downSpinCount = 0; 

    for(j = 0; j < N; j++)
    {
    if(sI.getOrbitalSpin(j) == 1) {upSpinCount++; }
    else                          {downSpinCount++;}
    }
    for(j = 0; j < upSpinCount; j++)
    {fileName << "U"; indexFileName << "U";}

    for(j = 0; j < downSpinCount; j++)
    {fileName << "D"; indexFileName << "D";}

    fileName << ".dat";
    indexFileName << ".dat";

    cout << (fileName.str()).c_str() << "  " << (indexFileName.str()).c_str() << endl;
 
    MatrixOutputFile      = fopen((fileName.str()).c_str(), "w");
    MatrixIndexOutputFile = fopen((indexFileName.str()).c_str(), "w");
    //
    // Print out indexing header 
    //
    localSpinBlockSize = sdBasisPtr->getSpinBlockSize(k);
    
    fprintf(MatrixIndexOutputFile,"################################################################################\n");
    fprintf(MatrixIndexOutputFile,"Indexing for spin block matrix file : %s \n",(fileName.str()).c_str());
    fprintf(MatrixIndexOutputFile,"Orbital integrals data file         : %s \n",inputIntegralDataFile);
    fprintf(MatrixIndexOutputFile,"################################################################################\n\n");
    fprintf(MatrixIndexOutputFile,"ROW : Spin Orbitals \n");
    fprintf(MatrixIndexOutputFile,"=================== \n");
    //
    // Output index list
    
    indexAndSpin.initialize(N,2);
    orbitalIndices.initialize(localSpinBlockSize,N);
    indexMap.initialize(localSpinBlockSize);
    swapSign.initialize(localSpinBlockSize);

    for(i = 0; i < localSpinBlockSize; i++)              
    {
    IspinBlockGlobalIndex = sdBasisPtr->getSpinBlockGlobalIndex(k,i);// global spin block index
    IblockIndex = sdBasisPtr->spinBlockToStd[IspinBlockGlobalIndex]; // std. basis index
    sI          = sdBasisPtr->operator()(IblockIndex);                          // ith element of basis
    swapSign(i) = 1.0;
    
    for(j = 0; j < N; j++)
    {
    if(sI.getOrbitalSpin(j) == 1)
    {indexAndSpin(j,0) = sI.getSpatialOrbitalIndex(j)+1; indexAndSpin(j,1) =  1;}
    else
    {indexAndSpin(j,0) = sI.getSpatialOrbitalIndex(j)+1; indexAndSpin(j,1) = -1;}
    }
    //
    // Rearrange the order of the basis elements, so those with upSpin come first
    // Every pairwise swap introduces a sign change in the basis element, so keep track of
    // this sign change.
    // 
    swapFlag    = 1;
    swapSignTmp = 1.0;
    while(swapFlag == 1)
    {
      swapFlag = 0;
      for(j = 0; j < N-1; j++)
      {
        if(indexAndSpin(j,1) < indexAndSpin(j+1,1))
        {
        swapFlag     = 1;
        swapSignTmp *= -1.0;
        indexTmp     = indexAndSpin(j,0);
        spinTmp      = indexAndSpin(j,1);
        indexAndSpin(j,0)   = indexAndSpin(j+1,0);
        indexAndSpin(j,1)   = indexAndSpin(j+1,1);
        indexAndSpin(j+1,0) = indexTmp;
        indexAndSpin(j+1,1) = spinTmp;
        }
      }
    }
    //
    // Rearrange the order of the basis elements so that the orbital index increases
    // from left to right for those of the same spin. 
    //
    // upSpin
    //
    swapFlag    = 1;
    while(swapFlag == 1)
    {
      swapFlag = 0;
      for(j = 0; j < upSpinCount-1; j++)
      {
        if(indexAndSpin(j,0) > indexAndSpin(j+1,0))
        {
        swapFlag     = 1;
        swapSignTmp *= -1.0;
        indexTmp     = indexAndSpin(j,0);
        spinTmp      = indexAndSpin(j,1);
        indexAndSpin(j,0)   = indexAndSpin(j+1,0);
        indexAndSpin(j,1)   = indexAndSpin(j+1,1);
        indexAndSpin(j+1,0) = indexTmp;
        indexAndSpin(j+1,1) = spinTmp;
        }
      }
    }
    //
    // downSpin
    //
    swapFlag    = 1;
    while(swapFlag == 1)
    {
      swapFlag = 0;
      for(j = upSpinCount; j < N-1; j++)
      {
        if(indexAndSpin(j,0) > indexAndSpin(j+1,0))
        {
        swapFlag     = 1;
        swapSignTmp *= -1.0;
        indexTmp     = indexAndSpin(j,0);
        spinTmp      = indexAndSpin(j,1);
        indexAndSpin(j,0)   = indexAndSpin(j+1,0);
        indexAndSpin(j,1)   = indexAndSpin(j+1,1);
        indexAndSpin(j+1,0) = indexTmp;
        indexAndSpin(j+1,1) = spinTmp;
        }
      }
    }
    //
    // Capture swap sign and the orbital indices associated with the new orbital ordering.
    //
    swapSign(i) = swapSignTmp;
    for(j = 0; j < N; j++)
    {
    orbitalIndices(i,j) = indexAndSpin(j,0);
    }
    //
    // Initialize index mapping
    //
    indexMap(i) = i;
    }
    //
    // Sort the basis by their spatial orbital index, using dictionary order 
    // Use bubble sort because it's simple and this routine will not be 
    // used for very large size systems. 
    
    swapFlag    = 1;
    while(swapFlag == 1)
    {
    	swapFlag = 0;
    	for(i = 0; i < localSpinBlockSize-1; i++)              
    	{
        	if(compareOrbitalIndices(indexMap(i),indexMap(i+1),N,orbitalIndices) > 0)
        	{	
        		swapFlag = 1;
        		indexTmp     = indexMap(i);
        		indexMap(i)  = indexMap(i+1);
        		indexMap(i+1) = indexTmp;
        	}
    	}
    }
    // 
    // Output the orbital coefficients 
    //
    for(i = 0; i < localSpinBlockSize; i++)              
    {
    fprintf(MatrixIndexOutputFile,"%3ld : ",i+1);
    for(j = 0; j < N; j++)
    {
    if(j < upSpinCount)
    {fprintf(MatrixIndexOutputFile,"%4ld+", orbitalIndices(indexMap(i),j));}
    else
    {fprintf(MatrixIndexOutputFile,"%4ld-", orbitalIndices(indexMap(i),j));}
    }
    
    fprintf(MatrixIndexOutputFile,"\n");
    }	
    
    fclose(MatrixIndexOutputFile);
    //
    // Output matrix elements
    //
    
    for(i = 0; i < localSpinBlockSize; i++)              
    {
    outputFormat.str("");
    outputFormat << outputFormatString << " ";
    for(j = 0; j < localSpinBlockSize; j++)
    {
    fprintf(MatrixOutputFile,(outputFormat.str()).c_str(),swapSign(indexMap(i))*swapSign(indexMap(j))*spinBlockMatrix[k](indexMap(i),indexMap(j)));
    }
    fprintf(MatrixOutputFile,"\n");
    }
    
    fclose(MatrixOutputFile);
    }	

    cout << endl;
}
//
// This routine compares the ith and jth set of spatial orbital indices and determines if 
// the ith is "greater" or "less" than the jth set with respect to a 
// dictionary ordering. This routine is used for sorting the orbitals to create output 
// where the index of the the spin block entries is via a dictionary ordering of the 
// spatial orbitals. 
//
// Returns 1  if ith > jth 
// Returns  0 if ith = jth
// Returns -1 if ith < jth
//
int MultiParticleHamiltonianMatrix::compareOrbitalIndices(long i,long j,long nParticle,
UCLAQ::LongArray2d& orbitalIndices)
{
	long k;
	for(k = 0; k < nParticle; k++)
	{
		if(orbitalIndices(i,k) > orbitalIndices(j,k)) {return  1;}
	    if(orbitalIndices(i,k) < orbitalIndices(j,k)) {return -1;}
	}
	return 0;
}

} // UCLAQ namespace

 
