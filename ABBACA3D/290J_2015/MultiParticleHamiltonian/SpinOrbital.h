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
#ifndef __SpinOrbital__
#define __SpinOrbital__

namespace UCLAQ
{


class SpinOrbital
{
    public :

    SpinOrbital()
    {
    this->spin                = 0;
    this->energy              = 0.0; 
    this->spatialOrbitalIndex = -1;
    }

    SpinOrbital(long spatialOrbitalIndex, int spin)
    {
    this->spin                = spin;
    this->energy              = 0.0; 
    this->spatialOrbitalIndex = spatialOrbitalIndex;
    }

    SpinOrbital(long spatialOrbitalIndex, int spin, double energy)
    {
    this->spin                = spin;
    this->energy              = energy; 
    this->spatialOrbitalIndex = spatialOrbitalIndex;
    }

    SpinOrbital(const SpinOrbital& spO)
    {
    this->spin                = spO.spin;
    this->energy              = spO.energy; 
    this->spatialOrbitalIndex = spO.spatialOrbitalIndex;
    }

    void initialize()
    {
    this->spin                = 0;
    this->energy              = 0.0; 
    this->spatialOrbitalIndex = -1;
    }

    void initialize(long spatialOrbitalIndex, int spin)
    {
    this->spin                = spin;
    this->energy              = 0.0; 
    this->spatialOrbitalIndex = spatialOrbitalIndex;
    }

    void initialize(long spatialOrbitalIndex, int spin, double energy)
    {
    this->spin                = spin;
    this->energy              = energy; 
    this->spatialOrbitalIndex = spatialOrbitalIndex;
    }

    void initialize(const SpinOrbital& spO)
    {
    this->spin                = spO.spin;
    this->energy              = spO.energy; 
    this->spatialOrbitalIndex = spO.spatialOrbitalIndex;
    }

    int getSpin()
    {return this->spin;}

    long getSpatialOrbitalIndex()
    {return this->spatialOrbitalIndex;}

    void operator=(const SpinOrbital& spO)
    {
    this->spin                = spO.spin;
    this->energy              = spO.energy; 
    this->spatialOrbitalIndex = spO.spatialOrbitalIndex;
    }

    int                 spin;
    double            energy; 
    long spatialOrbitalIndex;
};

} // UCLAQ namespace
#endif 
