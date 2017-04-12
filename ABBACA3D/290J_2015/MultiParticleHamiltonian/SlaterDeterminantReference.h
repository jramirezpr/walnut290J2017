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

#ifndef __SlaterDeterminantBasis__
namespace UCLAQ
{
class SlaterDeterminantBasis;
}
#endif 

#ifndef __SlaterDeterminantReference__
#define __SlaterDeterminantReference__

#include <iostream>
using namespace std;

namespace UCLAQ
{

class SlaterDeterminantReference
{
    public :

    SlaterDeterminantReference()
    {
    this->sdIndex      = -1;
    this->sdBasis      =  0;
    }

    SlaterDeterminantReference(long sdIndex, SlaterDeterminantBasis* sdBasis)
    {
    this->sdIndex      = sdIndex;
    this->sdBasis      = sdBasis;
    }
   
    SlaterDeterminantReference(const SlaterDeterminantReference& sdR)
    {
    this->sdIndex      = sdR.sdIndex;
    this->sdBasis      = sdR.sdBasis;
    }

    void initialize()
    {
    this->sdIndex      = -1;
    this->sdBasis      =  0;
    }

    void initialize(long sdIndex, SlaterDeterminantBasis* sdBasis)
    {
    this->sdIndex      = sdIndex;
    this->sdBasis      = sdBasis;
    }

    void initialize(const SlaterDeterminantReference& sdR)
    {
    this->sdIndex      = sdR.sdIndex;
    this->sdBasis      = sdR.sdBasis;
    }

    long  getElectronCount();

    long  getBasisCount();

    long  getSpinOrbitalCount();

    long* getSpinOrbitalIndexList();

    int   getTotalSpin();

    // returns the index of the spin orbital at the ith location

    long  getSpinOrbitalIndex(long i);

    // returns the index of the spatial orbital associated with the
    // spin orbtial at the ith location

    long  getSpatialOrbitalIndex(long i);

    // returns the spin of the spin orbital at the ith location

    int getOrbitalSpin(long i);

//
//  Output 
//
    friend ostream& operator <<(ostream& out_stream, const SlaterDeterminantReference& R);

    SlaterDeterminantBasis* sdBasis;
    long                    sdIndex;

};

} // UCLAQ namespace

#endif 
