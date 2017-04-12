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

#include "SlaterDeterminantBasis.h"
#include "SlaterDeterminantReference.h"

namespace UCLAQ
{


long SlaterDeterminantReference::getSpinOrbitalIndex(long i)
{
    return sdBasis->sdBasisListData[sdBasis->electronCount*sdIndex + i];
}


long  SlaterDeterminantReference::getSpatialOrbitalIndex(long i)
{
    long spinOrbIndex =
    sdBasis->sdBasisListData[sdBasis->electronCount*sdIndex + i];
    return sdBasis->spinOrbitals[spinOrbIndex].getSpatialOrbitalIndex();
}

int  SlaterDeterminantReference::getOrbitalSpin(long i)
{
    long spinOrbIndex =
    sdBasis->sdBasisListData[sdBasis->electronCount*sdIndex + i];
    return sdBasis->spinOrbitals[spinOrbIndex].getSpin();

}

int  SlaterDeterminantReference::getTotalSpin()
{
    return sdBasis->sdBasisSpinData[sdIndex];
}


long SlaterDeterminantReference::getElectronCount()
{
    return sdBasis->getElectronCount();
}

long SlaterDeterminantReference::getBasisCount()
{
    return sdBasis->getBasisCount();
}

long SlaterDeterminantReference::getSpinOrbitalCount()
{
    return sdBasis->getSpinOrbitalCount();
}

long* SlaterDeterminantReference::getSpinOrbitalIndexList()
{
    return &sdBasis->sdBasisListData[sdIndex*(sdBasis->electronCount)];
}

ostream& operator <<(ostream& out_stream, const SlaterDeterminantReference& R)
{
    long nElectron = R.sdBasis->getElectronCount();
    long i;
    for(i = 0; i < nElectron; i++)
    {
    out_stream << R.sdBasis->sdBasisListData[R.sdIndex*nElectron + i] << " ";
    }
    if(R.sdBasis->sdBasisSpinData[R.sdIndex] < 0)
    {
    out_stream << " Total Spin : " << R.sdBasis->sdBasisSpinData[R.sdIndex] << "/2";
    }
    else
    {
    out_stream << " Total Spin :  " << R.sdBasis->sdBasisSpinData[R.sdIndex] << "/2";
    }
    return out_stream;
} 

} // UCLAQ namespace
 
