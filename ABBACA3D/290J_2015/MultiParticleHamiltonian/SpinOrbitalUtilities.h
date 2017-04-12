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

#ifndef __SpinOrbitalUtilities__
#define __SpinOrbitalUtilities__

#include <cstdlib>
using namespace std;

#include "SlaterDeterminantReference.h"

#define MAX_ELECTRON_COUNT 256

//
// This class impements a table for values of
// m choose n with m < M_max, n < N_max and
// m >= n.
//
// The operator()(m,n) returns m choose n.
//
//
// The table is contsructed using the recurrance
//
// (m choose n) + (m choose n+1) = (m+1 choose n+1)
//
//
namespace UCLAQ
{

class MchooseNlookup
{
    public : 

    MchooseNlookup()
    {
    	tableData = 0;    
    	this->M_max = 0;
        this->N_max = 0;
    }
    

    MchooseNlookup(long M_max, long N_max)
    {
    this->M_max = M_max;
    this->N_max = N_max;
    tableData = new long[(N_max*(1-N_max))/2 + M_max*N_max];

    long i; long j; 

    for(i = 1; i <= M_max; i++)
    {
    	setTableValue(i,1,i);
    }
    
    for(j = 2; j <= N_max;j++)
    {
    	i = j;
    	setTableValue(i,j,1);
        for(i = j+1; i <= M_max; i++)
        {
        setTableValue(i,j,this->operator()(i-1,j) + this->operator()(i-1,j-1));
        if(this->operator()(i,j) < 0) 
        {
        cout << "Overflow Error in MchooseNlookup construction";
        exit(1); 
        }
        }
    }
    }
    
    void initialize(long M_max, long N_max)
    {
    this->M_max = M_max;
    this->N_max = N_max;
    
    if(tableData != 0) delete [] tableData;
    tableData = new long[(N_max*(1-N_max))/2 + M_max*N_max];

    long i; long j; 

    for(i = 1; i <= M_max; i++)
    {
    	setTableValue(i,1, i);
    }
     
    for(j = 2; j <= N_max;j++)
    {
    	i = j;
    	setTableValue(i,j,1);
        for(i = j+1; i <= M_max; i++)
        {
        setTableValue(i,j,this->operator()(i-1,j) + this->operator()(i-1,j-1));
        if(this->operator()(i,j) < 0) 
        {
        cout << "Overflow Error in MchooseNlookup construction";
        exit(1); 
        }
        }
    } 
    }

    ~MchooseNlookup()
    {delete [] tableData;}

    void setTableValue(long m, long n, long val)
    {
       tableData[((n-1)*(2*M_max - n + 2))/2 + (m-n)] = val;
    }
    
    long operator()(long m, long n)
    {
    if(n == 0) return 1;
#ifdef _DEBUG
	if((m < n)||(m > M_max)||(n > N_max)||(n <= 0)||(m <= 0))
	{
	cerr <<" Error in MchooseNlookup " << endl;
	cerr <<" Initial range  : M_max  = " << M_max << ", N_Max = " << N_max << endl;
	cerr <<" Input arguments:      m = " << m << ",         n = " << n << endl;
	exit(1);
	} 
	return tableData[((n-1)*(2*M_max - n + 2))/2 + (m-n)];
#else
    return tableData[((n-1)*(2*M_max - n + 2))/2 + (m-n)];
#endif

    }


    long M_max;
    long N_max;

    long* tableData;

};


class SpinOrbitalUtilities
{

public :

    SpinOrbitalUtilities()
    {
    orbitalCount_max  = 0;
    electronCount_max = 0;
    };

    SpinOrbitalUtilities(long orbitalCount_max,long electronCount_max) :
    mChooseNlookup(orbitalCount_max,electronCount_max)
    {
    this->orbitalCount_max  = orbitalCount_max;
    this->electronCount_max = electronCount_max;
    }

    void initialize(long orbitalCount_max,long electronCount_max) 
    {
    this->orbitalCount_max  = orbitalCount_max;
    this->electronCount_max = electronCount_max;
    mChooseNlookup.initialize(orbitalCount_max,electronCount_max);
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

//
// This routine returns the number of orbitals that differ between
// the two SlaterDeterminantReferences 
//
// return 0 => spinOrb1 and spinOrb2 have same orbitals. 
//        1 => spinOrb1 and spinOrb2 differ in one   orbital.
//        2 => spinOrb1 and spinOrb2 differ in two   orbitals.
//        3 => spinOrb1 and spinOrb2 differ in three orbitals.
//        etc.
//
long getOrbitalDifferenceCount(SlaterDeterminantReference& sI,
SlaterDeterminantReference& sJ)
{
    long n         = sI.getElectronCount();
    long *spinOrb1 = sI.getSpinOrbitalIndexList();
    long *spinOrb2 = sJ.getSpinOrbitalIndexList();
    return getOrbitalDifferenceCount(spinOrb1, spinOrb2, n);
}


//
// This routine returns the number of orbitals that differ between
// the two lists spinOrb1 and spinOrb2 of size n.
//
// return 0 => spinOrb1 and spinOrb2 have same orbitals. 
//        1 => spinOrb1 and spinOrb2 differ in one   orbital.
//        2 => spinOrb1 and spinOrb2 differ in two   orbitals.
//        3 => spinOrb1 and spinOrb2 differ in three orbitals.
//        etc.
//
long getOrbitalDifferenceCount(long* spinOrb1, long* spinOrb2, long n)
{
    long* spinConcat = new long[2*n];
    long k;

    for(k = 0; k < n; k++)
    {
    spinConcat[k]     = spinOrb1[k];
    spinConcat[n + k] = spinOrb2[k];
    }

    shellSort(spinConcat,2*n);

    long pairCount = 0;
    k              = 0;
    while(k < 2*n -1)
    {
    if(spinConcat[k] == spinConcat[k+1])
    { k+=2; pairCount++;}
    else
    { k+=1;}
    }
    delete [] spinConcat;
    return n-pairCount;
}

void createComplementarySet(long* spinOrb,
long* complementarySpinOrb, long n, long m)
{
    long* spinOrbSpinList = new long[n+m];

    long k;

    for(k = 0; k < m; k++)
    {
    spinOrbSpinList[k] = k;
    }

    for(k = 0; k < n; k++)
    {
    spinOrbSpinList[k+m] = spinOrb[k];
    } 

    shellSort(spinOrbSpinList,n+m);

    k              = 0;
    long diffIndex = 0;
    while(k < (m+n) -1)
    {
        if(spinOrbSpinList[k] == spinOrbSpinList[k+1])
        {k+=2;}
        else
        {
        complementarySpinOrb[diffIndex] = spinOrbSpinList[k];
        diffIndex++;
        k+=1; 
        }
    }
    if(k == (m+n) - 1)
    {
    complementarySpinOrb[diffIndex] = spinOrbSpinList[k];
    }
    delete [] spinOrbSpinList;
}

//
//  Given a slater determinant basis whose spinOrbital indices 
//  are linearly ordered in increasing order,
//  this routine returns the index of the basis set element 
//  with that collection of spinOrbital indices. 
//
long getBasisIndex(long* spinOrbitalIndexList, 
long electronCount, long spinOrbitalCount)
 {
    long spinOrbitalIndices[MAX_ELECTRON_COUNT];

    long j; 
    long k;

    long Mval;
    long Nval;

    long index = 0;
    long N     = electronCount;
    long M     = spinOrbitalCount;

    for(j = 0; j < N; j++)
    {
    spinOrbitalIndices[j] = spinOrbitalIndexList[j];
    }

    shellSort(spinOrbitalIndices, N);

    for(k = 1; k <= spinOrbitalIndices[0]; k++)
    {
    Nval = N-1;Mval = M-k;
    //index += MchooseN(Mval,Nval);
    index += mChooseNlookup(Mval,Nval);
    }

    for(j = 0; j < N; j++)
    {
        for(k = spinOrbitalIndices[j]+2; k <= spinOrbitalIndices[j+1]; k++)
        {
        Nval = N-(j+2);
        Mval = M-k;
        //index += MchooseN(Mval,Nval);
        index += mChooseNlookup(Mval,Nval);
        }
    }

    return index;
}

// Given two slater determinant orbitals that differ by a single
// orbital, this routine evalues the index of the orbitals that differ.
//
// On return : a = index in sI of the differing orbital
//             b = index in sJ of the differing orbital
//
// 
void getSingleDifferenceIndices(SlaterDeterminantReference& sI,
SlaterDeterminantReference& sJ,long& a, long& b)
{

    long diffTemp[MAX_ELECTRON_COUNT];

    long n         = sI.getElectronCount();
    long *spinOrb1 = sI.getSpinOrbitalIndexList();
    long *spinOrb2 = sJ.getSpinOrbitalIndexList();

    long i;
    for(i = 0; i < n; i++)
    {
    diffTemp[i] = spinOrb1[i] - spinOrb2[i];
    }


    long iLower =   0;
    while(diffTemp[iLower] == 0){iLower++;}

    long iUpper =   n-1;
    while(diffTemp[iUpper] == 0){iUpper--;}

    if(iLower == iUpper) 
    {
    a = iUpper; 
    b = iUpper; 
    return;
    }

    if(spinOrb1[iLower] == spinOrb2[iLower+1])
    {
    b = iLower; a = iUpper; return;
    }
    else
    {
    a = iLower; b = iUpper; return;
    }

    a = -1;
    b = -1;
    return;
}

//
// Given a SlaterDeterminant, this routine determines the indices of 
// the aOrbital within it. 
//
void getOrbitalIndices(SlaterDeterminantReference& sI,
long aOrbital,long& aOrbitalIndex)
{
    int  found;
    long index;

    found = 0;
    index = 0;

    long *spinOrb = sI.getSpinOrbitalIndexList();

    while(found == 0)
    {
    if(spinOrb[index] == aOrbital) {found = 1;}
    else                           {index++;}
    }
    aOrbitalIndex = index;
}


//
// Given a SlaterDeterminant, this routine determines the indices of 
// the aOrbital and bOrbital within it. 
//
void getOrbitalIndices(SlaterDeterminantReference& sI,
long aOrbital,long bOrbital,long& aOrbitalIndex, long& bOrbitalIndex)
{
    int  found;
    long index;

    found = 0;
    index = 0;

    long *spinOrb = sI.getSpinOrbitalIndexList();

    while(found == 0)
    {
    if(spinOrb[index] == aOrbital) {found = 1;}
    else                           {index++;}
    }
    aOrbitalIndex = index;

    found = 0;
    index = 0;
    while(found == 0)
    {
    if(spinOrb[index] == bOrbital) {found = 1;}
    else                           {index++;}
    }
    bOrbitalIndex = index;
}



// Given two slater determinant orbitals that differ by a two
// orbitals, this routine evalues the index of the orbitals that differ.
//
// On return : (a,b) = indices in sI of the differing orbitals
//             (c,d) = indices in sJ of the differing orbitals
//
// This is currently implemented the "dumb" O(n^2) way, getting it
// to O(n) is someting to do later. 
//
void getDoubleDifferenceIndices(SlaterDeterminantReference& sI,
SlaterDeterminantReference& sJ,long& a, long& b,long& c, long& d)
{
    a = -1; b = -1; c = -1; d = -1;

    long n         = sI.getElectronCount();
    long *spinOrb1 = sI.getSpinOrbitalIndexList();
    long *spinOrb2 = sJ.getSpinOrbitalIndexList();

    long found = 0; 

    long i; 
    long j;

    int isIn;

    i = 0;
    while((found < 2)&&(i < n))
    {
    isIn  = 0;

    for(j = 0; j < n; j++){if(spinOrb1[i] == spinOrb2[j]){isIn = 1;}}

    if(isIn == 0)
    {
    found++; if(found == 1) {a = i;}; if(found == 2) {b = i;}
    }

    i++;
    }

    found = 0;
    i     = 0;
    while((found < 2)&&(i < n))
    {
    isIn  = 0;

    for(j = 0; j < n; j++){if(spinOrb2[i] == spinOrb1[j]){isIn = 1;}}

    if(isIn == 0)
    {
    found++; if(found == 1) {c = i;}; if(found == 2) {d = i;};
    }

    i++;
    }
}

    MchooseNlookup mChooseNlookup;
    long         orbitalCount_max; 
    long        electronCount_max;
};

} // namespace
#endif 
