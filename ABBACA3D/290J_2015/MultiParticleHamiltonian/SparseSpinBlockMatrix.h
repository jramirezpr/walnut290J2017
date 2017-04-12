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

#ifndef __SparseSpinBlockMatrix__
#define __SparseSpinBlockMatrix__
//
// Indexing starts at 0
//
// To create a sparse matrix representation using this class 
// requires the following steps:
//
// Step #1 Call initialize(...) with the matrix size
//
// Step #2 Apply the loops that determine the matrix elements, and
//         call incrementRowCount(i) when a non-zero element at (i,j) 
//         occurs.
//
// Step #3 Call createCoefficientArrays() to allocate storage for the
//         non-zero matrix elements.
//
// Step #4 Apply the loops that determine the matrix elements, but this
//         time call packCoefficientData(i,j,coeffValue) to pack each
//         matrix element in the array.
//
// Step #5 If desired call sortRowIndicies to sort the row indices so
//         that the indices are in increasing order.
//
//
/* 
   This class is templated with respect to a vector class Vtype.  
   ---------
   Vtype is a Vector class with (minimally) the following member functions:

   Vtype()                            (null constructor)

   initialize(const Vtype&)           (copy initializer)
   operator =                         (duplicate assignemnt)

*/

namespace UCLAQ
{


template <class Vtype> class SparseSpinBlockMatrix
{

public :
SparseSpinBlockMatrix()
{
    rowCount   = 0;

    rowSizes       = 0;
    coeffData      = 0;
    coeffIndex     = 0;
    rowOffset      = 0;
    rowFilledSizes = 0;
}

~SparseSpinBlockMatrix()
{
	destroyData();
}

SparseSpinBlockMatrix(SparseSpinBlockMatrix& S)
{
	rowCount       = 0;
    rowSizes       = 0;
    coeffData      = 0;
    coeffIndex     = 0;
    rowOffset      = 0;
    rowFilledSizes = 0;

    initialize(S);
}

void initialize(SparseSpinBlockMatrix& S)
{
	destroyData();
	rowCount             = S.rowCount;
	rowSizes             = new long[rowCount];
    rowFilledSizes       = new long[rowCount];
    rowOffset            = new long[rowCount];
    
    for(long i = 0; i < rowCount; i++)
    {
    rowSizes[i]       = S.rowSizes[i];
    rowFilledSizes[i] = S.rowFilledSizes[i];
    rowOffset[i]      = S.rowOffset[i];
    }
    
    totalDataCount = 0;
    for(long i = 0; i < rowCount; i++)
    {
    totalDataCount += rowSizes[i];
    }

    coeffData  = new double[totalDataCount];
    coeffIndex = new long[totalDataCount];

    for(long i = 0; i < totalDataCount; i++)
    {
    coeffData[i]  = S.coeffData[i];
    coeffIndex[i] = S.coeffIndex[i];
    }
}


void initialize(long rowCount)
{
    destroyData();
    
    this->rowCount       = rowCount;
    rowSizes             = new long[rowCount];
    rowFilledSizes       = new long[rowCount];
    rowOffset            = new long[rowCount];
    long i; 

    for(i = 0; i < rowCount; i++)
    {
    rowSizes[i]       = 0;
    rowFilledSizes[i] = 0;
    rowOffset[i]      = 0;
    }
}

void incrementRowCount(long rowIndex)
{
    rowSizes[rowIndex]++;
}

//
// Assuming that the row sizes have been set, this routine
// allocates memory for the coefficient arrays
//
void createCoefficientArrays()
{
    long i; 

    totalDataCount = 0;
    for(i = 0; i < rowCount; i++)
    {
    totalDataCount += rowSizes[i];
    }
//
//  Clean up any existing coefficient arrays 
    if(coeffData  != 0) delete [] coeffData;
    if(coeffIndex != 0) delete [] coeffIndex;
//
//  Allocate data for coefficient arrays 
//
    coeffData  = new double[totalDataCount];
    coeffIndex = new long[totalDataCount];
    
    for(i = 0; i < totalDataCount; i++)
    {
    coeffData[i]   = 0.0;
    coeffIndex[i]  = 0; 
    }
//
//  Set up row offsets
//
    rowOffset[0] = 0;
    for(i = 1; i < rowCount; i++)
    {
    rowOffset[i] = rowOffset[i-1] + rowSizes[i-1];
    }
    
}

void packCoefficientData(long rowIndex, long colIndex, double coeffValue)
{

#ifdef _DEBUG
    if(rowFilledSizes[rowIndex] > rowSizes[rowIndex]) 
    {
    cerr << "Error packing coefficients of SparseSpinBlockMatrix"   << endl;
    cerr << "Number of coefficients exceeds data allocation"        << endl;
    cerr << "rowIndex = " << rowIndex << " colIndex = " << colIndex << endl;
    exit(1);
    }
#endif 

    long index = rowOffset[rowIndex] + rowFilledSizes[rowIndex];
    coeffData[index]  = coeffValue;
    coeffIndex[index] = colIndex;
    rowFilledSizes[rowIndex]++;
}

double  operator()(long i, long j)
{
    for(long k = 0; k < rowFilledSizes[i]; k++)
    {
        if(j == coeffIndex[rowOffset[i] + k])
        {
        return coeffData[rowOffset[i] + k];
        }
    }
    return 0.0;
}

void apply(Vtype& V)
{
    vTmp.initialize(V); // initialize and copy V to vTmp
    long iBase;
    long i; long j;
    for(i = 0; i < rowCount; i++)
    {
        V(i) = 0.0;
        iBase   = rowOffset[i];

        for(j = 0; j < rowSizes[i]; j++)
        {
            V(i) += vTmp(coeffIndex[iBase + j])*coeffData[iBase + j];
        }
    }
}


void sortWithData(long* indexVals, double* coeffVals, long n)
{
	long i,j,inc;
//
//  Shell sort
//
	long       iTmp;
    double     cTmp;
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
			iTmp  =indexVals[i-1];
            cTmp  =coeffVals[i-1];
			j=i;
			while (indexVals[(j-inc)-1] > iTmp) 
            {
			    indexVals[j-1]      = indexVals[(j-inc)-1];
                coeffVals[j-1]      = coeffVals[(j-inc)-1];
				j -= inc;
				if (j <= inc) break;
			}
			indexVals[j-1]  =   iTmp;
            coeffVals[j-1]  =   cTmp;
		}
	} while (inc > 1);
}

//
// This routine sorts the elements of each row so that the
// indices associated with the coefficients are in increasing 
// order -- hoping that this reduces the number of page accesses.
//

void sortRows()
{
    long i; 
    double* coeffDataPtr;
    long  * coeffIndexPtr;
    
    for(i = 0; i < rowCount; i++)
    {
    coeffDataPtr  = &coeffData[rowOffset[i]];
    coeffIndexPtr = &coeffIndex[rowOffset[i]];
    sortWithData(coeffIndexPtr, coeffDataPtr, rowSizes[i]);
    }
}


//
// Returns the number of non-zero elements in the matrix
//
long getTotalDataCount()
{return totalDataCount;}

void destroyData()
{
    if(coeffData        != 0) delete [] coeffData;
    if(coeffIndex       != 0) delete [] coeffIndex;
    if(rowSizes         != 0) delete [] rowSizes;
    if(rowOffset        != 0) delete [] rowOffset;
    if(rowFilledSizes   != 0) delete [] rowFilledSizes;
    
    coeffData      = 0;
    coeffIndex     = 0;
    rowSizes       = 0;
    rowOffset      = 0;
    rowFilledSizes = 0;
}

long    rowCount;       // number of rows in the matrix
long*   rowSizes;       // the number of non-zero elements in each row
long*   rowFilledSizes; // the number of filled elements in each row 
long*   rowOffset;      // the base index of each row in the coeff arrays
long    totalDataCount; // the total number of non-zero element in the matrix

double* coeffData;      // array for coefficient values
long*   coeffIndex;     // column indices of the coefficients

Vtype vTmp;             // temporary for applyForwardOp;

};

} // UCLAQ namespace
#endif 
