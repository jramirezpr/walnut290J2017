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

#ifndef __ElectronOperators__
#define __ElectronOperators__

#include "RealSpatialOrbitalIntegrals.h"
#include "SlaterDeterminantReference.h"

namespace UCLAQ
{

class ElectronOperators
{
public :

    ElectronOperators(RealSpatialOrbitalIntegrals& rOrb):
    rOrbitalIntregrals(rOrb)
    {};

    double oneElectronDiagOp(SlaterDeterminantReference& sI)
    {
        long N = sI.getElectronCount();

        long m;

        long spatialOrbIndex;

        double val = 0.0;
        for(m = 0; m < N; m++)
        {
        spatialOrbIndex   = sI.getSpatialOrbitalIndex(m);
        val += rOrbitalIntregrals(spatialOrbIndex, spatialOrbIndex);
        }
        return val;
    }


    double twoElectronDiagOp(SlaterDeterminantReference& sI)
    {
        long N = sI.getElectronCount();

        long m; 
        long n;

        long spatialOrbIndex_m;
        long spatialOrbIndex_n;

        double val = 0.0;

        int spin_m;
        int spin_n;

        for(m = 0; m < N; m++)
        {
        spatialOrbIndex_m   = sI.getSpatialOrbitalIndex(m);
        spin_m              = sI.getOrbitalSpin(m);
        for(n = 0; n < N; n++)
        {
        if(n != m)
        {
        spatialOrbIndex_n   = sI.getSpatialOrbitalIndex(n);
        spin_n              = sI.getOrbitalSpin(n);
        
        val += rOrbitalIntregrals(spatialOrbIndex_m,spatialOrbIndex_n,
                                  spatialOrbIndex_m,spatialOrbIndex_n);

        if((spin_m + spin_n) != 0)
        {
        val -= rOrbitalIntregrals(spatialOrbIndex_m,spatialOrbIndex_n,
                                  spatialOrbIndex_n,spatialOrbIndex_m);
        }
        
        }}}

        val *= 0.5;
        return val;
    }

//
//  Off diagonal operator when sI and sJ differ by exactly one spin orbital.
//  a is the index of the differing orbital in sI and b is the index
//  of the differing orbital in sJ.
//
    double oneElectronOffDiagOp(SlaterDeterminantReference& sI,
           long a, SlaterDeterminantReference& sJ, long b)
    {
    double abSign;

    int  spin_a;
    int  spin_b;

    long spatialOrbIndex_a;
    long spatialOrbIndex_b;

    spin_a = sI.getOrbitalSpin(a);
    spin_b = sJ.getOrbitalSpin(b);

    if(spin_a != spin_b) return 0.0;

    abSign = 1.0;
    if(((a+b)%2) != 0 ) abSign = -1.0; // Using the fact that 
                                       // (-1)^[(a+1) + (b+1)] = (-1)^[a+b]

    spatialOrbIndex_a = sI.getSpatialOrbitalIndex(a);
    spatialOrbIndex_b = sJ.getSpatialOrbitalIndex(b);

    return abSign*rOrbitalIntregrals(spatialOrbIndex_a, 
           spatialOrbIndex_b);
    }
//
//  Off diagonal operator when sI and sJ differ by exactly one spin orbital.
//  a is the index of the differing orbital in sI and b is the index
//  of the differing orbital in sJ.
//
    double twoElectronOffDiagOp(SlaterDeterminantReference& sI,
           long a, SlaterDeterminantReference& sJ, long b)
{
    long N = sI.getElectronCount();

    long n;

    int  spin_a;
    int  spin_b;
    int  spin_n;

    long spatialOrbIndex_a;
    long spatialOrbIndex_b;
    long spatialOrbIndex_n;

    spin_a = sI.getOrbitalSpin(a);
    spin_b = sJ.getOrbitalSpin(b);

    spatialOrbIndex_a = sI.getSpatialOrbitalIndex(a);
    spatialOrbIndex_b = sJ.getSpatialOrbitalIndex(b);

    double val = 0.0;

    //
    // First term 
    // 
    if(spin_a == spin_b)
    {
        for(n = 0; n < N; n++)
        {
        if(n != a)
        {
        spatialOrbIndex_n = sI.getSpatialOrbitalIndex(n);
        val += rOrbitalIntregrals(spatialOrbIndex_a,spatialOrbIndex_n,
                                  spatialOrbIndex_b,spatialOrbIndex_n);
        }
        }
    }
    //
    // Second term 
    //  
    for(n = 0; n < N; n++)
    {
        if(n != a)
        {
        spin_n            = sI.getOrbitalSpin(n);
        spatialOrbIndex_n = sI.getSpatialOrbitalIndex(n);
        if((spin_a == spin_n)&&(spin_b == spin_n))
        {
        val -= rOrbitalIntregrals(spatialOrbIndex_a,spatialOrbIndex_n,
                                  spatialOrbIndex_n,spatialOrbIndex_b);
        }
        }
    }

    double abSign = 1.0;               // Using the fact that
    if(((a+b)%2) != 0 ) abSign = -1.0; // (-1)^[(a+1) + (b+1)] = (-1)^[a+b]

    return abSign*val; 
}
//
//  Off diagonal operator when sI and sJ differ by exactly two spin orbitals.
//  (a,b) are the pair of indices of the differing orbitals in sI and 
//  (c,d) are the pair of indices of the differning orbtials in sJ. 
//
    double twoElectronOffDiagOp(SlaterDeterminantReference& sI,
    long a, long b, SlaterDeterminantReference& sJ, long c, long d)
    {
    double abcdSign = 1.0;               
    if(((a+b+c+d)%2) != 0 ) abcdSign = -1.0; 


    int  spin_a;
    int  spin_b;
    int  spin_c;
    int  spin_d;

    long spatialOrbIndex_a;
    long spatialOrbIndex_b;
    long spatialOrbIndex_c;
    long spatialOrbIndex_d;

    spin_a = sI.getOrbitalSpin(a);
    spin_b = sI.getOrbitalSpin(b);
    spin_c = sJ.getOrbitalSpin(c);
    spin_d = sJ.getOrbitalSpin(d);

    spatialOrbIndex_a = sI.getSpatialOrbitalIndex(a);
    spatialOrbIndex_b = sI.getSpatialOrbitalIndex(b);

    spatialOrbIndex_c = sJ.getSpatialOrbitalIndex(c);
    spatialOrbIndex_d = sJ.getSpatialOrbitalIndex(d);

    double val = 0.0;

    if((spin_a == spin_c)&&(spin_b == spin_d))
    {
    val += abcdSign*rOrbitalIntregrals(spatialOrbIndex_a,spatialOrbIndex_b,
                                  spatialOrbIndex_c,spatialOrbIndex_d);
    }

    if((spin_a == spin_d)&&(spin_b == spin_c))
    {
    val -= abcdSign*rOrbitalIntregrals(spatialOrbIndex_a,spatialOrbIndex_b,
                              spatialOrbIndex_d,spatialOrbIndex_c);
    }

    return val;
    }

    RealSpatialOrbitalIntegrals& rOrbitalIntregrals;
};

} // UCLAQ namespace
#endif

