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

//
// In this basis the elements are ordered with spin orbital indices
// in increasing order.
//
//


#include "SlaterDeterminantReference.h"
#include "SlaterDetermiantListConstructor.h"
#include "SpinOrbital.h"

#ifndef __SlaterDeterminantBasis__
#define __SlaterDeterminantBasis__

namespace UCLAQ
{

class SlaterDeterminantBasis
{
    public :

    SlaterDeterminantBasis()
    {
    this->sdBasisListData      = 0;
    this->sdBasisSpinData      = 0;
    this->sdBasisSpinBlockSize = 0;

    this->spinBlockToStd       = 0;
    this->stdToSpinBlock       = 0;
    this->spinBlockLocalIndex      = 0;
    this->stdToSpinBlockLocalIndex = 0;

    this->spinOrbitals         = 0;

    this->electronCount    = 0;
    this->spinOrbitalCount = 0; 
    this->basisCount       = 0;
    
    this->indexList        = 0;
	this->dataIndex        = 0;
	this->spinOrbitalIndices = 0;
    }

    SlaterDeterminantBasis(long electionCount, long spinOrbitalCount,
    SpinOrbital* spinOrbitalArray)
    {
    this->sdBasisListData      = 0;
    this->sdBasisSpinData      = 0;
    this->sdBasisSpinBlockSize = 0;

    this->spinBlockToStd       = 0;
    this->stdToSpinBlock       = 0;
    this->spinBlockLocalIndex  = 0;
    this->stdToSpinBlockLocalIndex = 0;

    this->spinOrbitals         = 0;

    this->electronCount    = 0;
    this->spinOrbitalCount = 0; 
    this->basisCount       = 0;
    
    this->indexList        = 0;
	this->dataIndex        = 0;
	this->spinOrbitalIndices = 0;
	
    initialize(electionCount, spinOrbitalCount,spinOrbitalArray);
    }


    ~SlaterDeterminantBasis()
    {destroyData();}


    void initialize()
    {
    destroyData();
    this->electronCount    = 0;
    this->spinOrbitalCount = 0; 
    this->basisCount       = 0;
    }

    void initialize(long electionCount, long spinOrbitalCount,
    SpinOrbital* spinOrbitalArray)
    {
    destroyData();

    long i;
    
    this->electronCount    = electionCount;
    this->spinOrbitalCount = spinOrbitalCount; 

    //
    // Capture spinOrbitalArray
    //

    spinOrbitals = new SpinOrbital[spinOrbitalCount];
    for(i = 0; i < spinOrbitalCount; i++)
    {
    spinOrbitals[i].initialize(spinOrbitalArray[i]);
    }

    //
    // Compute the number of Slater determininat basis elements
    //

    SlaterDetermiantListConstructor createSD;

    basisCount 
    = createSD.countDeterminants(spinOrbitalCount,electionCount);

    //
    // Allocate array for basis information 
    // (consist of integers identifying spin orbitals
    // in the individual basis element). 
    // 

    sdBasisListData = new long[basisCount*electronCount];

    //SlaterDetermiantListConstructor::sdListData = sdBasisListData;
    //SlaterDetermiantListConstructor::M          = spinOrbitalCount;
    //SlaterDetermiantListConstructor::N          = electionCount;
    //SlaterDetermiantListConstructor::count      = 0;

    createSD.initialize(spinOrbitalCount,electionCount,0,sdBasisListData);
    createSD.createSlaterDerminantList();

    //
    // Initialize total spin data 
    //

    sdBasisSpinData = new int[basisCount];

    long            k; 
    int     totalSpin;
    long* orbitalList; 

    for(i = 0; i < basisCount; i++)
    {
        totalSpin   = 0;
        orbitalList = &sdBasisListData[i*electronCount];

        for(k = 0; k < electronCount; k++)
        {
        totalSpin +=spinOrbitals[orbitalList[k]].spin;
        }
    
        sdBasisSpinData[i] = totalSpin; 
    }

    //
    // Determine Spin Block sizes
    //

    sdBasisSpinBlockSize = new long[electronCount+1];

    for(i = 0; i < electronCount+1; i++)
    {sdBasisSpinBlockSize[i] = 0;}


    for(i = 0; i < basisCount; i++)
    {
    k = (sdBasisSpinData[i] + electronCount)/2;
    sdBasisSpinBlockSize[k] += 1;
    }

    //
    // create mapping to from an index blocking based on
    // total spin  to the "standard" index. This indexing is
    // the global spin block indexing. 
    //
  
    stdToSpinBlock = new long[basisCount];
    spinBlockToStd = new long[basisCount];

    long* spinBlockIndices = new long[electronCount+1];

    spinBlockIndices[0] = 0;

    for(i = 1; i < electronCount+1; i++)
    {
      spinBlockIndices[i] = sdBasisSpinBlockSize[i-1] + spinBlockIndices[i-1];
    }

    for(i = 0; i < basisCount; i++)
    {
    k = (sdBasisSpinData[i] + electronCount)/2;
    stdToSpinBlock[i]                   = spinBlockIndices[k];
    spinBlockToStd[spinBlockIndices[k]] = i;
    spinBlockIndices[k] += 1;
    }

    delete [] spinBlockIndices;
    
    //
    // Create a mapping from global spin block indices to local spin block index
    //
    
    spinBlockLocalIndex      = new long[basisCount];
    stdToSpinBlockLocalIndex = new long[basisCount];
    long spinBlockTmp;
    long blockIndex;

    for(i = 0; i < basisCount; i++)
    {
    getSpinBlockLocalIndex(i,spinBlockTmp,spinBlockLocalIndex[i]);
    }
    
    for(i = 0; i < basisCount; i++)
    {
    blockIndex = stdToSpinBlock[i];
    stdToSpinBlockLocalIndex[i] =  spinBlockLocalIndex[blockIndex];
    }
    //
    // Determine size of arrays for spin orbital list -> slater determinant basis mapping


	long*  levelSize;
	long*  indexOffset;
	long*  allocationOffset;
	
	long level;
	long indexMin;
	long indexMax;
	
	long sizeD = 0;
	long sizeP = 0;
	
	long cumulativeSum;
	
	//
	// Return if electron count is equal to 1, as the orbital index -> Slater determinant index is
	// the identity.
	//

	if(electronCount == 1) return;
	
    levelSize        = new long[electronCount-1];
	indexOffset      = new long[electronCount-1];
	allocationOffset = new long[electronCount-1];
	
	
    indexMin =  0;
	indexMax =  (spinOrbitalCount-electronCount);
	sizeD    = 0;
	sizeP    = (indexMax - indexMin) + 1;
	level = 1;
	
	for(i = 0; i < electronCount-1; i++)
	{
    levelSize[i]        = 0;
    indexOffset[i]      = 0;
    allocationOffset[i] = 0;
	}
	levelSize[0] = sizeP;
	getInverseMappingSizes(indexMin, indexMax, level, sizeD, sizeP, levelSize,electronCount);

	indexList = new  long[sizeP];
	dataIndex = new  long[sizeD];
	cumulativeSum = 0;
	for(i = 0; i < electronCount-1; i++)
	{
	indexOffset[i]      =  cumulativeSum;
	allocationOffset[i] =  cumulativeSum;
	cumulativeSum      +=  levelSize[i];
	}
	
	indexMin = 0;
	indexMax = (spinOrbitalCount-electronCount);
	sizeD    = 0;
	level    = 1;
	levelSize[0] = (indexMax - indexMin) + 1;
	setInverseMappingIndices(indexMin,indexMax,level,indexList, dataIndex,sizeD, indexOffset,
	allocationOffset,electronCount);
	
 
	delete []       levelSize;
	delete []      indexOffset;
	delete [] allocationOffset;
	
	//
	// Pack indices of the Slater determinants into dataIndex 
	//
	long kOffset;
	for(i = 0; i < basisCount; i++)
    {
        orbitalList = &sdBasisListData[i*electronCount];
        k = 0; 
        kOffset = indexList[orbitalList[k]];
        for(k = 1; k < electronCount-1; k++)
        {
        kOffset = indexList[kOffset +(orbitalList[k] -orbitalList[k-1])-1];
        }
        k = electronCount-1;
		dataIndex[kOffset +(orbitalList[k] -orbitalList[k-1])-1] = i;
    }

    //
    // Allocate data for inverse mapping function
    //
    spinOrbitalIndices = new long[electronCount];
    
    
}

long getBasisIndexFromSpinOrbitals(long* orbitalList)
{
	if(electronCount == 1) return orbitalList[0];
	long k;
	long kOffset;
    	
	for(k = 0; k < electronCount; k++)
    {
    spinOrbitalIndices[k] = orbitalList[k];
    }

    shellSort(spinOrbitalIndices, electronCount);
    
	k = 0;
	kOffset = indexList[spinOrbitalIndices[k]];
	for(k = 1; k < electronCount-1; k++)
	{
			kOffset = indexList[kOffset +(spinOrbitalIndices[k] - spinOrbitalIndices[k-1])-1];
	}
	k = electronCount-1;
	return dataIndex[kOffset +(spinOrbitalIndices[k] - spinOrbitalIndices[k-1])-1];
}	

void shellSort(long* vals,long n)
{
	long i,j,inc;

	long   valTmp;
	inc=1;
	do 
    {
    inc *= 3;
    inc++;
	} while (inc <= n);
	do 
    {
		inc /= 3;
		for (i=inc+1;i<=n;i++) 
        {
			valTmp=vals[i-1];
			j=i;
			while (vals[(j-inc)-1] > valTmp) {
				vals[j-1]  = vals[(j-inc)-1];
				j -= inc;
				if (j <= inc) break;
			}
			vals[j-1]  = valTmp;
		}
	} while (inc > 1);
}

SlaterDeterminantReference operator()(long i)
{
    SlaterDeterminantReference R(i, this);
    return R;
}


long  getBasisCount()         {return basisCount;}
long  getElectronCount()      {return electronCount;}
long  getSpinOrbitalCount()   {return spinOrbitalCount;}

//
// Spin block utilities 
//
// Returns the size of the ith spin block, where the index 
// i ranges from 0 to electronCount 
//
// (e.g. there are electronCount+1 spin blocks). 
// 

long  getSpinBlockSize(long i){return sdBasisSpinBlockSize[i];}

//
// Given a globals spin block index, this returns the local spin block
// index.
//

void  getSpinBlockLocalIndex(long spinBlockGlobalIndex, long& spinBlock,
long& spinBlockIndexLocal)
{
    spinBlock              = 0;
    long spinBlockSizeSum  = sdBasisSpinBlockSize[spinBlock];
    long spinDiff          = spinBlockGlobalIndex - spinBlockSizeSum;

    while(spinDiff >= 0)
    {
    spinBlock++;
    spinBlockSizeSum += sdBasisSpinBlockSize[spinBlock];
    spinDiff          = spinBlockGlobalIndex - spinBlockSizeSum;
    }

    spinBlockSizeSum   -= sdBasisSpinBlockSize[spinBlock];
    spinBlockIndexLocal = spinBlockGlobalIndex - spinBlockSizeSum;
}


//
// Given a local spin block index, this returns the global spin block
// index.
//

long getSpinBlockGlobalIndex(long spinBlock,long spinBlockIndexLocal)
{
    long k;

    long globalSpinBlockIndex = 0;

    for(k = 0; k <= spinBlock-1; k++)
    {
    globalSpinBlockIndex += sdBasisSpinBlockSize[k];
    }

    globalSpinBlockIndex += spinBlockIndexLocal;
    
    return globalSpinBlockIndex;
}

    
void destroyData()
{
    if(sdBasisListData != 0)      delete [] sdBasisListData; 
    if(sdBasisSpinData != 0)      delete [] sdBasisSpinData; 
    if(sdBasisSpinBlockSize != 0) delete [] sdBasisSpinBlockSize; 

    if(spinBlockToStd != 0)       delete [] spinBlockToStd;
    if(stdToSpinBlock != 0)       delete [] stdToSpinBlock;
    
    if(spinBlockLocalIndex != 0)      delete [] spinBlockLocalIndex;
    if(stdToSpinBlockLocalIndex != 0) delete [] stdToSpinBlockLocalIndex;

    if(spinOrbitals    != 0)      delete [] spinOrbitals;

    sdBasisListData       = 0;
    sdBasisSpinData       = 0;
    sdBasisSpinBlockSize  = 0;
    spinBlockToStd        = 0;
    stdToSpinBlock        = 0;
    spinBlockLocalIndex   = 0;
    stdToSpinBlockLocalIndex = 0;

    spinOrbitals          = 0;
    
    if(indexList          != 0) delete [] indexList;
    if(dataIndex          != 0) delete [] dataIndex;
    if(spinOrbitalIndices != 0) delete [] spinOrbitalIndices;
    indexList          = 0;
	dataIndex          = 0;
	spinOrbitalIndices = 0;
}
//
//  Utility routines for constructing the mapping from a list of orbital indices to
//  the slater determinant basis index.
//
void getInverseMappingSizes(long indexMin, long indexMax, long level, long &sizeD, long &sizeP, 
long* levelSizes, long N)
{
	long j;
	
	if(level == (N-1))
	{
    for(j = indexMin; j <= indexMax; j++)
	{
	    sizeD += (indexMax-j) + 1;
	}
	return;
	}
	
	for(j = indexMin; j <= indexMax; j++)
	{
	    sizeP             += (indexMax-j) + 1;
	    levelSizes[level] += (indexMax-j) + 1;
	}
	
	for(j = indexMin; j <= indexMax; j++)
	{
	getInverseMappingSizes(j+1, indexMax+1, level+1, sizeD, sizeP, levelSizes,N);
	}
}

void setInverseMappingIndices(long indexMin, long indexMax, long level,long*indexList, long* dataIndex,
long &sizeD, long* indexOffset,long* allocationOffset,long N)
{
	long j;
	
	if(level == (N-1))
	{
    for(j = indexMin; j <= indexMax; j++)
	{
		indexList[indexOffset[level-1]] = sizeD;
		indexOffset[level-1]++;
	    sizeD += (indexMax-j) + 1;
	}
	return;
	}
	

	for(j = indexMin; j <= indexMax; j++)
	{
	    indexList[indexOffset[level-1]] = allocationOffset[level];
	    indexOffset[level-1]++;
	    allocationOffset[level]  += (indexMax-j) + 1;
	}
	
	for(j = indexMin; j <= indexMax; j++)
	{
	setInverseMappingIndices(j+1, indexMax+1, level+1,indexList,dataIndex,sizeD, indexOffset,allocationOffset,N);
	}
}


    long*      sdBasisListData;
    int*       sdBasisSpinData;

    long* sdBasisSpinBlockSize;

    long* spinBlockToStd;
    long* stdToSpinBlock;

    SpinOrbital* spinOrbitals;

    long        electronCount;
    long     spinOrbitalCount;
    long           basisCount;
    
    long* spinBlockLocalIndex; 
    long* stdToSpinBlockLocalIndex;
    
    //
    // Data for list of orbital indices -> slater determinant index mapping
    //
    long*  indexList;
	long*  dataIndex;
	long*  spinOrbitalIndices;
	
};
} // UCLAQ namespace

#endif
 
