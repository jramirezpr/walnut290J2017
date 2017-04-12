/*
 * UCLAQ_GridFunction3d.h
 *
 *  Created on: Jun 27, 2015
 *      Author: anderson
 *
 *
 *
 * Decisions: Extending DoubleVector3D so that move semantics can be incorporated into
 * the underlying vector operations on grid values without having to duplicate all
 * member functions utilizing move semantics.
 *
 * Providing a Values member function so that data values can be accessed as if this
 * grid function implementation contained an instance of DoubleVector1D Values.
 *
 * Revised: Nov. 26, 2015 
 *          Jan. 26, 2016
 */

/*
#############################################################################
#
# Copyright  2015 Chris Anderson
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

#ifndef _UCLAQ_GridFunction3d_
#define _UCLAQ_GridFunction3d_


#include "UCLAQ_GridFunction1d.h"
#include "UCLAQ_GridFunction2d.h"

#include <functional>
#include <iostream>

#include "../DoubleVectorNd/UCLAQ_DoubleVector3d.h"
using namespace std;

namespace UCLAQ
{
class GridFunction3d : public DoubleVector3d
{

public :

GridFunction3d() : DoubleVector3d()
{
    this->xPanels = 0;
    this->yPanels = 0;
    this->zPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;
    this->hz      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
    this->zMin = 0.0;
    this->zMax = 1.0;
}

GridFunction3d(const GridFunction3d& G) : DoubleVector3d(G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;
}


GridFunction3d(DoubleVector3d&& G) : DoubleVector3d((DoubleVector3d&&)G)
{
    this->xPanels = 0;
    this->yPanels = 0;
    this->zPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;
    this->hz      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
    this->zMin = 0.0;
    this->zMax = 1.0;
}


GridFunction3d(long xPanels, double hx, long yPanels, double hy, long zPanels, double hz)
: DoubleVector3d(xPanels+1,yPanels+1,zPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->hx     = hx;
    this->hy     = hy;
    this->hz     = hz;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;
    this->zMin = -(zPanels*hz)/2.0;
    this->zMax =  (zPanels*hz)/2.0;

	this->setToValue(0.0);
}


GridFunction3d(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax,
               long zPanels, double zMin, double zMax) : DoubleVector3d(xPanels+1,yPanels+1,zPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;
    this->zMin = zMin;
    this->zMax = zMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);

	this->setToValue(0.0);
}


GridFunction3d(GridFunction3d&& G) : DoubleVector3d((DoubleVector3d&&)G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;
}

virtual ~GridFunction3d(){}

GridFunction3d& operator=(const GridFunction3d& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;
}

// Propagate values

DoubleVector3d::operator=(G);

return *this;
}

GridFunction3d& operator=(GridFunction3d&& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;
}

// Propagate values

DoubleVector3d::operator=((DoubleVector3d&&)G);
return *this;
}

GridFunction3d& operator=(DoubleVector3d&& G)
{
//
// Since underlying operations that result in DoubleVector3d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector3d::operator=((DoubleVector3d&&)G);
return *this;
}


GridFunction3d& operator=(DoubleVector3d& G)
{
//
// Since underlying operations that result in DoubleVector3d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector3d::operator=((DoubleVector3d&)G);
return *this;
}

//####  Incremental operators  ###

void operator*=(double alpha)
{UCLAQ::DoubleVector3d::operator*=(alpha);}

void operator/=(double alpha)
{UCLAQ::DoubleVector3d::operator/=(alpha);}


void operator+=(const GridFunction3d& G)
{UCLAQ::DoubleVector3d::operator+=(G);}

void operator-=(const GridFunction3d& G)
{UCLAQ::DoubleVector3d::operator-=(G);}

void operator*=(const GridFunction3d& G)
{UCLAQ::DoubleVector3d::operator*=(G);}

void operator/=(const GridFunction3d& G)
{UCLAQ::DoubleVector3d::operator/=(G);}


void operator+=(const DoubleVector3d& G)
{UCLAQ::DoubleVector3d::operator+=(G);}

void operator-=(const DoubleVector3d& G)
{UCLAQ::DoubleVector3d::operator-=(G);}

void operator*=(const DoubleVector3d& G)
{UCLAQ::DoubleVector3d::operator*=(G);}

void operator/=(const DoubleVector3d& G)
{UCLAQ::DoubleVector3d::operator/=(G);}


void operator*=(const  std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) *= F(xPos,yPos,zPos);
    }}}
}

void operator/=(const std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) /= F(xPos,yPos,zPos);
    }}}
}


void operator+=(const  std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) += F(xPos,yPos,zPos);
    }}}
}

void operator-=(const  std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) -= F(xPos,yPos,zPos);
    }}}
}

// ######################################################################

// Initialization

void initialize()
{
    DoubleVector3d::initialize();
    this->xPanels = 0;
    this->yPanels = 0;
    this->zPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;
    this->hz      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
    this->zMin = 0.0;
    this->zMax = 1.0;
}

void initialize(const GridFunction3d& G)
{
	DoubleVector3d::initialize(G);
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;
}

void initialize(long xPanels, double hx, long yPanels, double hy, long zPanels, double hz)
{
    DoubleVector3d::initialize(xPanels+1,yPanels+1,zPanels+1);
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->hx     = hx;
    this->hy     = hy;
    this->hz     = hz;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;
    this->zMin = -(zPanels*hz)/2.0;
    this->zMax =  (zPanels*hz)/2.0;

	this->setToValue(0.0);
}


void initialize(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax,
               long zPanels, double zMin, double zMax)
{
    DoubleVector3d::initialize(xPanels+1,yPanels+1,zPanels+1);
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;
    this->zMin = zMin;
    this->zMax = zMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);

	this->setToValue(0.0);
}

GridFunction3d* newDuplicate() const
{
    GridFunction3d* Mptr = new GridFunction3d(*this);
    return Mptr;
}

bool isNull()
{
if(dataPtr == nullptr) return true;
return false;
}

void createProductFunction(const GridFunction1d& funX, const GridFunction1d& funY, const GridFunction1d& funZ)
{
	this->xPanels = funX.getXpanelCount();
    this->yPanels = funY.getXpanelCount();
    this->zPanels = funZ.getXpanelCount();

    this->xMin = funX.getXmin();
    this->xMax = funX.getXmax();
    this->yMin = funY.getXmin();
    this->yMax = funY.getXmax();
    this->zMin = funZ.getXmin();
    this->zMax = funZ.getXmax();

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);

	DoubleVector3d::initialize(xPanels+1,yPanels+1,zPanels+1);

	long i; long j; long k;

	double fX; double fY; double fZ;

	for(i = 0; i <= xPanels; i++)
	{
	fX = funX.Values(i);
	for(j = 0; j <= yPanels; j++)
	{
	fY = funY.Values(j);
	for(k = 0; k <= zPanels; k++)
	{
	fZ = funZ.Values(k);

	Values(i,j,k) = fX*fY*fZ;
	}}}
}

void createProductFunction(const GridFunction2d& funXY,  const GridFunction1d& funZ)
{
	this->xPanels = funXY.getXpanelCount();
    this->yPanels = funXY.getYpanelCount();
    this->zPanels = funZ.getXpanelCount();

    this->xMin = funXY.getXmin();
    this->xMax = funXY.getXmax();

    this->yMin = funXY.getYmin();
    this->yMax = funXY.getYmax();

    this->zMin = funZ.getXmin();
    this->zMax = funZ.getXmax();

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);

	DoubleVector3d::initialize(xPanels+1,yPanels+1,zPanels+1);

	long i; long j; long k;

	double fXY; double fZ;

	for(i = 0; i <= xPanels; i++)
	{
    for(j = 0; j <= yPanels; j++)
	{
	fXY = funXY.Values(i,j);
	for(k = 0; k <= zPanels; k++)
	{
	fZ = funZ.Values(k);

	Values(i,j,k) = fXY*fZ;
	}}}
}

void createProductFunction(const GridFunction1d& funX,const  GridFunction2d& funYZ)
{
	this->xPanels = funX.getXpanelCount();
    this->yPanels = funYZ.getXpanelCount();
    this->zPanels = funYZ.getYpanelCount();

    this->xMin = funX.getXmin();
    this->xMax = funX.getXmax();

    this->yMin = funYZ.getXmin();
    this->yMax = funYZ.getXmax();

    this->zMin = funYZ.getYmin();
    this->zMax = funYZ.getYmax();

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);

	DoubleVector3d::initialize(xPanels+1,yPanels+1,zPanels+1);

	long i; long j; long k;
	double fX; double fYZ;

	for(j = 0; j <= yPanels; j++)
	{
	for(k = 0; k <= zPanels; k++)
	{
	fYZ = funYZ.Values(j,k);
	for(i = 0; i <= xPanels; i++)
	{
	fX = funX.Values(i);
	Values(i,j,k) = fX*fYZ;
	}}}
}


void specify(std::function<double(double,double,double)> F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) = F(xPos,yPos,zPos);
    }}}
}

void squareValues()
{
    transformValues([](double x){return x*x;});
}

void zeroNegativePart()
{
    transformValues([](double x){if(x < 0.0){return 0.0;} return x;});
}

//
// Dot product of all values scaled with mesh widths
//           (Trapezoidal Method)
//

double scaledDot(const GridFunction3d& V) const
{
    double dotVal = 0.0;

    long i; long j; long k;
    double hxhyhz = hx*hy*hz;

//  Interior points weighted by 1

    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz;
    }}}
//
//  Interior face points (weighted by 1/2)
//
    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    j = yPanels;
    for(i = 1; i <  xPanels; i++)
    {
    for(k = 1; k <  zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    //
    k = 0;
    for(i = 1; i <  xPanels; i++)
    {
    for(j = 1; j <  yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}
//
//  Interior corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = 0;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }


    i = 0;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = 0;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }
//
    j = 0;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    j = yPanels;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    j = 0;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    j = yPanels;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }
//
//  Corner points (weighted by 1/8)
//
    dotVal += Values(0,0,0)*V.Values(0,0,0)*hxhyhz*.125;
    dotVal += Values(xPanels,0,0)*V.Values(xPanels,0,0)*hxhyhz*.125;

    dotVal += Values(0,yPanels,0)*V.Values(0,yPanels,0)*hxhyhz*.125;
    dotVal += Values(xPanels,yPanels,0)*V.Values(xPanels,yPanels,0)*hxhyhz*.125;

    dotVal += Values(0,0,zPanels)*V.Values(0,0,zPanels)*hxhyhz*.125;
    dotVal += Values(xPanels,0,zPanels)*V.Values(xPanels,0,zPanels)*hxhyhz*.125;

    dotVal += Values(0,yPanels,zPanels)*V.Values(0,yPanels,zPanels)*hxhyhz*.125;
    dotVal += Values(xPanels,yPanels,zPanels)*V.Values(xPanels,yPanels,zPanels)*hxhyhz*.125;

    return dotVal;
}



//  Trapezoidal method

double integralTrapezoidal(std::function<double(double)> F) const
{
    double intVal = 0.0;
    long i; long j; long k;

    double hxhyhz = hx*hy*hz;
//
//  Interior points weighted by 1
//
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz;
    }}}
//
//  Interior face points (weighted by 1/2)
//
    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    j = yPanels;
    for(i = 1; i < xPanels; i++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    //
    k = 0;
    for(i = 1; i <  xPanels; i++)
    {
    for(j = 1; j <  yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}
//
//  Interior corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = 0;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    //

    i = 0;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = 0;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }
//
    j = 0;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    j = yPanels;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    j = 0;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    j = yPanels;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }
//
//  Corner points (weighted by 1/8)
//
    intVal += F(Values(0,0,0))*hxhyhz*.125;
    intVal += F(Values(xPanels,0,0))*hxhyhz*.125;

    intVal += F(Values(0,yPanels,0))*hxhyhz*.125;
    intVal += F(Values(xPanels,yPanels,0))*hxhyhz*.125;

    intVal += F(Values(0,0,zPanels))*hxhyhz*.125;
    intVal += F(Values(xPanels,0,zPanels))*hxhyhz*.125;

    intVal += F(Values(0,yPanels,zPanels))*hxhyhz*.125;
    intVal += F(Values(xPanels,yPanels,zPanels))*hxhyhz*.125;

    return intVal;
}

double norm1() const
{
  return integralTrapezoidal([](double x){return abs(x);});
}

// Trapezoidal Method Integral Approximation

double integralTrapezoidal() const
{
    return integralTrapezoidal([](double x){return x;});
}

// 2-Norm^2 computed using the trapezoidal rule approximation

double norm2() const
{
    return sqrt(abs(integralTrapezoidal([](double x){return x*x;})));
}

// 2-Norm^2 computed using the trapezoidal rule approximation

double norm2squared() const
{
    return integralTrapezoidal([](double x){return x*x;});
}

//
// Riemann Sum
//

double  integralRiemann() const
{
//  All points weighted by hxhy

    double intVal = 0.0;
    double hxhy = hx*hy;

    long dCount  = (xPanels+1)*(yPanels+1)*(zPanels+1);

    for(long k =0; k < dCount; k++)
    {
    intVal += dataPtr[k];
    }

    intVal *= hxhy;

    return intVal;
}

//  Average computed using Riemann sum approximation of the integral

double getRiemannAverage() const
{
    double avgValue;
    avgValue = integralRiemann();
    return avgValue/((xMax-xMin)*(yMax-yMin)*(zMax-zMin));
}

//  Average computed using Trapezoidal approximation of the integral

double getTrapezoidalAverage() const
{
    double avgValue;
    avgValue = integralTrapezoidal();
    return avgValue/((xMax-xMin)*(yMax-yMin)*(zMax-zMin));
}


double min() const
{
//  Compute function minimum

    long dCount      = (xPanels+1)*(yPanels+1)*(zPanels+1);
    double minValue  = dataPtr[0];
    for(long k = 1; k < dCount; k++)
    {
    minValue = ( minValue < dataPtr[k] ) ? minValue : dataPtr[k];
    }
    return minValue;

}

double max() const
{
//  Compute function maximum

    long dCount  = (xPanels+1)*(yPanels+1)*(zPanels+1);

    double maxValue  = dataPtr[0];
    for(long k = 1; k < dCount; k++)
    {
    maxValue = ( maxValue > dataPtr[k] ) ? maxValue : dataPtr[k];
    }
    return maxValue;
}

GridFunction2d getConstantZslice(long zIndex) const //(x-y function)
{
    GridFunction2d R(xPanels,xMin,xMax,yPanels,yMin,yMax);
    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    R.Values(i,j) = Values(i,j,zIndex);
    }}
    return R;
}

GridFunction2d getConstantYslice(long yIndex) const //(x-z function)
{
	GridFunction2d R(xPanels,xMin,xMax,zPanels,zMin,zMax);
    for(long i = 0; i <= xPanels; i++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    R.Values(i,k) = Values(i,yIndex,k);
    }}
    return R;
}

GridFunction2d getConstantXslice(long xIndex) const //(y-z function)
{
	GridFunction2d R(yPanels,yMin,yMax,zPanels,zMin,zMax);
    for(long j = 0; j <= yPanels; j++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    R.Values(j,k) = Values(xIndex,j,k);
    }}
    return R;
}

GridFunction1d getConstantYZslice(long yIndex, long zIndex) const  // (x function)
{
	GridFunction1d R(xPanels,xMin,xMax);
	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = Values(i,yIndex,zIndex);
	}
	return R;
}

GridFunction1d getConstantXZslice(long xIndex, long zIndex) const  //( y function)
{
	GridFunction1d R(yPanels,yMin,yMax);
	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = Values(xIndex,j,zIndex);
	}
	return R;
}

GridFunction1d getConstantXYslice(long xIndex, long yIndex) const  //( z function)
{
	GridFunction1d R(zPanels,zMin,zMax);
	for(long k = 0; k <= zPanels; k++)
	{
	R.Values(k) = Values(xIndex,yIndex,k);
	}
	return R;
}


virtual int isXperiodic() const {return 0;}
virtual int isYperiodic() const {return 0;}
virtual int isZperiodic() const {return 0;}

double getHx()   const  {return            hx;}
double getXmin() const  {return          xMin;}
double getXmax() const  {return          xMax;}
long   getXpanelCount() const {return xPanels;}

double getHy()   const  {return            hy;}
double getYmin() const  {return          yMin;}
double getYmax() const  {return          yMax;}
long   getYpanelCount() const {return yPanels;}

double getHz()   const  {return            hz;}
double getZmin() const  {return          zMin;}
double getZmax() const  {return          zMax;}
long   getZpanelCount() const {return zPanels;}


DoubleVector3d getValues() const
{
    return DoubleVector3d(*this);
}

DoubleVector3d* getValuesPointer()
{
    return (DoubleVector3d*)this;
}

const DoubleVector3d* getValuesPointer() const
{
    return (DoubleVector3d*)this;
}

/// Set boundary values to specified value

void setBoundaryValues(double value)
{
    long i; long j; long k;

	i = 0;
	for(j = 0; j <= yPanels; j++)
	{
	for(k = 0; k <= zPanels; k++)
	{
		Values(i,j,k) = value;
	}}

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
    for(k = 0; k <= zPanels; k++)
    {
     Values(i,j,k) = value;
    }}


    j = 0;
    for(i = 0; i <= xPanels; i++)
    {
    for(k = 0; k <= zPanels; k++)
    {
    Values(i,j,k) = value;
    }}

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(k = 0; k <= zPanels; k++)
    {
    Values(i,j,k) = value;
    }}

//

    k = 0;
    for(i = 0; i <=  xPanels; i++)
    {
    for(j = 0; j <=  yPanels; j++)
    {
    Values(i,j,k) = value;
    }}

    k = zPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(j = 0; j <= yPanels; j++)
    {
    Values(i,j,k) = value;
    }}
}

void enforcePeriodicity()
{
    long i; long j; long k;

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
    for(k = 0; k <= zPanels; k++)
    {
     Values(i,j,k) = Values(0,j,k);
    }}

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(k = 0; k <= zPanels; k++)
    {
    Values(i,j,k) = Values(i,0,k);
    }}

    k = zPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(j = 0; j <= yPanels; j++)
    {
    Values(i,j,k) = Values(i,j,0);
    }}
}

//  Returns true if the input grid function is structurally identical,
//  e.g. a coincident computational domain and mesh size
//

bool isCoincident(const GridFunction3d& V)
{
	double domainDiffTol    = 1.0e-13;

    long discretizationDiff = abs(xPanels - V.xPanels)
                            + abs(yPanels - V.yPanels)
                            + abs(zPanels - V.zPanels);

    double domainDiff = (abs(xMin - V.xMin) + abs(xMax - V.xMax))/abs(xMax-xMin)
                      + (abs(yMin - V.yMin) + abs(yMax - V.yMax))/abs(yMax-yMin)
                      + (abs(zMin - V.zMin) + abs(zMax - V.zMax))/abs(zMax-zMin);

    if((discretizationDiff > 0) || (domainDiff > domainDiffTol)) return false;
    return true;
}

//  Grid Geometry

    double xMin; double xMax;   // Computational Region is [xMin,xMax]x
    double yMin; double yMax;   //                         [yMin,yMax]x
    double zMin; double zMax;   //                         [zMin,zMax]

    long   xPanels;             // Number of panels
    long   yPanels;
    long   zPanels;

    double hx;                  // mesh width
    double hy;
    double hz;

//###################################################################
//          Values Access (Alternative to operator())
//###################################################################

#ifdef _DEBUG
    double&  Values(long i1, long i2, long i3)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    assert(boundsCheck(i3, 0, index3Size-1,3));
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };

    const double&  Values(long i1, long i2, long i3) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    assert(boundsCheck(i3, 0, index3Size-1,3));
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };
#else
    /*!
    Returns a reference to the element with index (i1,i2,i3) - indexing
    starting at (0,0,0).
    */
    inline double&  Values(long i1, long i2, long i3)
    {
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };

    /*!
    Returns a reference to the element with index (i1,i2,i3) - indexing
    starting at (0,0,0).
    */
    inline const double&  Values(long i1, long i2, long i3) const
    {
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };
#endif




//###################################################################
//                      Bounds Checking
//###################################################################

#ifdef _DEBUG
        bool domainCheck() const
        {
        cerr << "XXX UCLAQ::GridFunction3d Error  XXX" << endl;
        cerr << "Left side of assignment must be a non-null GridFunction3d instance." << endl;
        cerr << endl;
        return false;
        }
#else
        bool domainCheck() const {return true;}
#endif


};

}


#endif /* UCLAQ_GRIDFUNCTION_UCLAQ_GRIDFUNCTION1D_H_ */
