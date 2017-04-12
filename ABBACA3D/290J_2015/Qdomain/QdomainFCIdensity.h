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

#ifndef __QdomainFCIdensity__
#define __QdomainFCIdensity__

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"

namespace UCLAQ
{


class QdomainFCIdensity
{

public: 	

QdomainFCIdensity()
{
initialize();
}

QdomainFCIdensity(int electronCount, vector < string >& spinBlockIdentifiers,
		          vector <vector <double> >& spinBlockEnergies,
		          vector < vector <UCLAQ::DoubleVector2d> >& spinBlockDensities,
				  vector <UCLAQ::GridFunction3d>&  externalBasisVectors)
{
	initialize(electronCount,spinBlockIdentifiers,spinBlockEnergies,spinBlockDensities, externalBasisVectors);
}

void initialize()
{
	electronCount         = 0;
	basisVectorsPtr       = 0;
	spinBlockDensitiesPtr = 0;
	spinBlockEnergies.clear();
	spinBlockIdentifiers.clear();
}


void initialize(int electronCount, vector < string >& spinBlockIdentifiers,
		        vector <vector <double> >& spinBlockEnergies,
		        vector < vector <UCLAQ::DoubleVector2d> >& spinBlockDensities,
				vector<UCLAQ::GridFunction3d>&  externalBasisVectors)

{
	this->electronCount         = electronCount;
	this->spinBlockIdentifiers  = spinBlockIdentifiers;
	this->spinBlockEnergies     = spinBlockEnergies;
    this->spinBlockDensitiesPtr = &spinBlockDensities;
	this->basisVectorsPtr       = &externalBasisVectors;
}




string getSpinBlockString(int spinBlockIndex)
{
    return spinBlockIdentifiers[spinBlockIndex];
}

void getDensityMatrix(int waveFunctionIndex,
string& spinBlockDesignator, UCLAQ::DoubleVector2d& densityMatrix)
{
    try{waveFunctionIndexBoundsCheck(waveFunctionIndex,spinBlockDesignator);}
    catch(const runtime_error&e ){throw e;}

    int spinBlockIndex = getSpinBlockIndex(spinBlockDesignator);
    densityMatrix      = (*spinBlockDensitiesPtr)[spinBlockIndex][waveFunctionIndex];
}

double getEnergy(int waveFunctionIndex, string spinBlockDesignator)
{
    try{waveFunctionIndexBoundsCheck(waveFunctionIndex,spinBlockDesignator);}
    catch(const runtime_error&e ){throw e;}

    int spinBlockIndex = getSpinBlockIndex(spinBlockDesignator);
    return spinBlockEnergies[spinBlockIndex][waveFunctionIndex];
}

void getFCIdensity(int waveFunctionIndex,string spinBlockDesignator, UCLAQ::GridFunction3d& FCIdensity)
{
    try{waveFunctionIndexBoundsCheck(waveFunctionIndex,spinBlockDesignator);}
    catch(const runtime_error&e ){throw e;}

    int spinBlockIndex = getSpinBlockIndex(spinBlockDesignator);
    createDensity((*spinBlockDensitiesPtr)[spinBlockIndex][waveFunctionIndex],FCIdensity);
}

int getSpinBlockIndex(string& spinBlockDesignator)
{
    string spinBlockString = spinBlockDesignator;

    std::transform(spinBlockString.begin(),spinBlockString.end(),spinBlockString.begin(), ::toupper);
    int   upCount = 0;
    int downCount = 0;
    for(int i = 0; i < (int)spinBlockString.size(); i++	)
    {
    if(spinBlockString[i] == 'U')   upCount++;
    if(spinBlockString[i] == 'D') downCount++;
    }

    if(upCount > downCount) return electronCount - upCount;
    return electronCount - downCount;
    }

void waveFunctionIndexBoundsCheck(int waveFunctionIndex, string& spinBlockDesignator)
{
    int spinBlockIndex = getSpinBlockIndex(spinBlockDesignator);

    string errMessage;
    std::ostringstream s;

    if((long)spinBlockDesignator.size() != electronCount)
    {
    s << electronCount;
    errMessage  = "\nXXX QdomainExec Error XXX\n";
    errMessage += "Unacceptable spin configuration specification.\n";
    errMessage += "Incorrect length of spin configuration specification.\n";
    errMessage += "Offending spin configuration specification : " + spinBlockDesignator + "\n";
    errMessage += "Correct specification requires exactly " + s.str() + " designators.\n";
    throw runtime_error(errMessage);
    }

    if((long)spinBlockDensitiesPtr->size() <= spinBlockIndex)
    {
    errMessage  = "\nXXX QdomainExec Error XXX\n";
    errMessage += "Requested spin configuration wave function not found.\n";
    errMessage += "Offending spin configuration specification :" + spinBlockDesignator + "\n";
    throw runtime_error(errMessage);
    }


    if((long)(*spinBlockDensitiesPtr)[spinBlockIndex].size() <= waveFunctionIndex)
    {
    errMessage  = "\nXXX QdomainExec XXX\n";
    errMessage += "Requested wave function index out of bounds.\n";
    s << (*spinBlockDensitiesPtr)[spinBlockIndex].size()-1;
    errMessage += "Acceptable wave function index range      : 0 to " + s.str() + "\n";
    s.str("");
    s << waveFunctionIndex;
    errMessage += "Offending wave function index             : " + s.str() + "\n";
    throw runtime_error(errMessage);
    }
}


protected:

void createDensity(UCLAQ::DoubleVector2d& densityMatrix,UCLAQ::GridFunction3d& density)
{
    vector <UCLAQ::GridFunction3d>& basisVectors = *basisVectorsPtr;

	gridFunTemp.initialize(basisVectors[0]);
	gridFunTemp.setToValue(0.0);

    density.initialize(basisVectors[0]);
	density.setToValue(0.0);

	gridFunTemp  = basisVectors[0];
	gridFunTemp *= basisVectors[0];

	long i; long j;

	long basisCount = densityMatrix.getIndex1Size();

	for(i = 0; i < basisCount; i++)
	{
		for(j = 0; j < basisCount; j++)
		{
		gridFunTemp  = basisVectors[i];
		gridFunTemp *= basisVectors[j];
		gridFunTemp *= densityMatrix(i,j);
		density += gridFunTemp;
	}}
}

    // Data obtained from FCI

    int                                         electronCount;
    vector < string >                           spinBlockIdentifiers;
    vector <vector <double> >                   spinBlockEnergies;
    vector < vector <UCLAQ::DoubleVector2d> >*  spinBlockDensitiesPtr;
    vector <UCLAQ::GridFunction3d>*             basisVectorsPtr;
    
	UCLAQ::GridFunction3d                gridFunTemp;
	UCLAQ::GridFunction3dUtility            gUtility;

};

} // UCLAQ namespace
#endif  // __QdomainFCIdensity__
