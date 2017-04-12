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

#ifndef __SlaterDetermiantListConstructor__
#define __SlaterDetermiantListConstructor__

#include <iostream>
#include <cstdlib>
using namespace std;

#include "SpinOrbitalUtilities.h"

namespace UCLAQ
{

class SlaterDetermiantListConstructor
{
public:

	SlaterDetermiantListConstructor()
	{
	M          = 0;
	N          = 0;
	count      = 0;
	sdListData = 0;
	}

	void initialize(long M, long N, long count, long* sdListData)
	{
	this->M     = M;
	this->N     = N;
	this->count = count;
	this->sdListData = sdListData;
	}

	void fillOrbitals(long electronIndex, long orbitalStartIndex,long baseCount)
	{
    long i; long k;

    if(electronIndex == N) 
    {
    count++; return;
    }

    baseCount = count;
    for(i = orbitalStartIndex; i <= M - (N-electronIndex); i++)
    {
     for(k = 0; k < electronIndex; k++)
     {
     sdListData[count*N + k] = sdListData[baseCount*N +k];
     }
     sdListData[count*N + electronIndex]    = i;
     fillOrbitals(electronIndex + 1, i+1,baseCount);
    }
	}

	void createSlaterDerminantList()
	{
    fillOrbitals(0,0,0);
	}

	long countDeterminants(long Morb, long Nelec)
	{
	MchooseNlookup M(Morb,Nelec);
	return M(Morb,Nelec);
	}
//
// The following routine needs to be made robust with respect 
// to sizes of M and N and the possiblity of integer overflow. 
/*
static long MchooseN(long M, long& N)
{
    if(M < N)
    {
    cerr << " M < N in MchooseN routine " << endl; 
    exit(1); 
    }
    if(M == N) return 1;
    long val = (M*MchooseN(M-1,N))/(M-N);

    if(val <=0)
    {
    cerr << " Integer Overflow in MchooseN routine " << endl; 
    exit(1); 
    }

    return val;
};
*/
	long M;
	long N;
	long count;
	long*  sdListData; // pointer to external data
	                   // (data not managed by this class)
};

} // UCLAQ namespace


#endif 
