/*
 * UCLAQ_GridFunction2d.h
 *
 *  Created on: Jun 27, 2015
 *      Author: anderson
 *
 *
 *
 * Decisions: Extending DoubleVector2D so that move semantics can be incorporated into
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

#ifndef _UCLAQ_GridFunction2d_
#define _UCLAQ_GridFunction2d_

#include <functional>
#include <iostream>
using namespace std;

#include "../GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "../DoubleVectorNd/UCLAQ_DoubleVector2d.h"

namespace UCLAQ
{
class GridFunction2d : public DoubleVector2d
{

public :

GridFunction2d() : DoubleVector2d()
{
    this->xPanels = 0;
    this->yPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
}

GridFunction2d(const GridFunction2d& G) : DoubleVector2d(G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
}


GridFunction2d(DoubleVector2d&& G) : DoubleVector2d((DoubleVector2d&&)G)
{
    this->xPanels = 0;
    this->yPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
}


GridFunction2d(long xPanels, double hx, long yPanels, double hy)
: DoubleVector2d(xPanels+1,yPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->hx     = hx;
    this->hy     = hy;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;
	this->setToValue(0.0);
}


GridFunction2d(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax)
: DoubleVector2d(xPanels+1,yPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
	this->setToValue(0.0);
}


GridFunction2d(GridFunction2d&& G) : DoubleVector2d((DoubleVector2d&&)G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
}

virtual ~GridFunction2d(){}

// Initialization

void initialize()
{
    DoubleVector2d::initialize();
    this->xPanels = 0;
    this->yPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
}

void initialize(const GridFunction2d& G)
{
	DoubleVector2d::initialize(G);
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
}

void initialize(long xPanels, double hx, long yPanels, double hy)
{
	DoubleVector2d::initialize(xPanels+1,yPanels+1);

    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->hx     = hx;
    this->hy     = hy;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;
	this->setToValue(0.0);
}

void initialize(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax)
{
    DoubleVector2d::initialize(xPanels+1,yPanels+1);

    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
	this->setToValue(0.0);
}


GridFunction2d* newDuplicate() const
{
    GridFunction2d* Mptr = new GridFunction2d(*this);
    return Mptr;
}

bool isNull()
{
if(dataPtr == nullptr) return true;
return false;
}

double getHx()   const  {return hx;}
double getXmin() const  {return   xMin;}
double getXmax() const  {return   xMax;}
long   getXpanelCount() const {return xPanels;}

double getHy()   const  {return  hy;}
double getYmin() const  {return yMin;}
double getYmax() const  {return yMax;}
long   getYpanelCount() const {return yPanels;}

virtual bool isXperiodic() const {return false;}
virtual bool isYperiodic() const {return false;}


DoubleVector2d getValues() const
{
    return DoubleVector2d(*this);
}

DoubleVector2d* getValuesPointer()
{
    return (DoubleVector2d*)this;
}

const DoubleVector2d* getValuesPointer() const
{
    return (DoubleVector2d*)this;
}



//###################################################################
//                  Values Access
//
//       An alternative to using operator()(...)
//###################################################################

#ifdef _DEBUG
    double& Values(long i1, long i2)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +  i2 + i1*index2Size);
    };

    const double&  Values(long i1, long i2) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +   i2  + i1*index2Size);
    };
#else
    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
    */
    inline double&  Values(long i1, long i2)
    {
    return *(dataPtr +  i2 + i1*index2Size);
    };

    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
     */
    inline const double&  Values(long i1, long i2) const
    {
    return *(dataPtr +  i2  + i1*index2Size);
    };
#endif



GridFunction2d& operator=(const GridFunction2d& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
}

// Propagate values

DoubleVector2d::operator=(G);

return *this;
}

GridFunction2d& operator=(GridFunction2d&& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
}

// Propagate values

DoubleVector2d::operator=((DoubleVector2d&&)G);
return *this;
}

GridFunction2d& operator=(DoubleVector2d& G)
{
//
// Since underlying operations that result in DoubleVector2d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector2d::operator=((DoubleVector2d&)G);
return *this;
}

GridFunction2d& operator=(DoubleVector2d&& G)
{
//
// Since underlying operations that result in DoubleVector2d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector2d::operator=((DoubleVector2d&&)G);
return *this;
}

//####  Incremental operators  ###


void operator*=(double alpha)
{UCLAQ::DoubleVector2d::operator*=(alpha);}

void operator/=(double alpha)
{UCLAQ::DoubleVector2d::operator/=(alpha);}


void operator+=(const GridFunction2d& G)
{UCLAQ::DoubleVector2d::operator+=(G);}

void operator-=(const GridFunction2d& G)
{UCLAQ::DoubleVector2d::operator-=(G);}

void operator*=(const GridFunction2d& G)
{UCLAQ::DoubleVector2d::operator*=(G);}

void operator/=(const GridFunction2d& G)
{UCLAQ::DoubleVector2d::operator/=(G);}



void operator+=(const DoubleVector2d& G)
{UCLAQ::DoubleVector2d::operator+=(G);}

void operator-=(const DoubleVector2d& G)
{UCLAQ::DoubleVector2d::operator-=(G);}

void operator*=(const DoubleVector2d& G)
{UCLAQ::DoubleVector2d::operator*=(G);}

void operator/=(const DoubleVector2d& G)
{UCLAQ::DoubleVector2d::operator/=(G);}


void operator*=(const  std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) *= F(xPos,yPos);
    }}
}

void operator/=(const std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) /= F(xPos,yPos);
    }}
}


void operator+=(const  std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) += F(xPos,yPos);
    }}
}


void operator-=(const  std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) -= F(xPos,yPos);
    }}
}

//############################################################################

void createProductFunction(const GridFunction1d& funX, const GridFunction1d& funY)
{
	this->xPanels = funX.getXpanelCount();
    this->yPanels = funY.getXpanelCount();


    this->xMin = funX.getXmin();
    this->xMax = funX.getXmax();
    this->yMin = funY.getXmin();
    this->yMax = funY.getXmax();

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);

	DoubleVector2d::initialize(xPanels+1,yPanels+1);

	long i; long j;

	double fX; double fY;

	for(i = 0; i <= xPanels; i++)
	{
	fX = funX.Values(i);
	for(j = 0; j <= yPanels; j++)
	{
	fY = funY.Values(j);
	Values(i,j) = fX*fY;
	}}
}


void specify(const std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) = F(xPos,yPos);
    }}
}

void squareValues()
{
    transformValues([](double x){return x*x;});
}

void zeroNegativePart()
{
    transformValues([](double x){if(x < 0.0){return 0.0;} return x;});
}


// Dot product of all values scaled with mesh widths
// Trapezoidal method approximation of the integral.

double scaledDot(const GridFunction2d& V) const
{
    double dotVal = 0.0;

    long i; long j;
    double hxhy = hx*hy;

//  Interior points weighted by 1
//
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy;
    }}
//
//  Interior face points (weighted by 1/2)
//
    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

    j = yPanels;
    for(i = 1; i <  xPanels; i++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

//
//  Corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;


    i = 0;
    j = yPanels;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;

    i = xPanels;
    j = 0;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;

    i = xPanels;
    j = yPanels;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;

    return dotVal;
}

// Trapezoidal method approximation of the integral of F(Values(i,j))

double integralTrapezoidal(std::function<double(double)> F) const
{
    double intVal = 0.0;

    long i; long j;
    double hxhy = hx*hy;

//  Interior points weighted by 1

    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j))*hxhy;
    }}

//  Interior face points (weighted by 1/2)

    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

    j = yPanels;
    for(i = 1; i <  xPanels; i++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

//
//  Corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    intVal += F(Values(i,j))*hxhy*0.25;

    i = 0;
    j = yPanels;
    intVal += F(Values(i,j))*hxhy*0.25;

    i = xPanels;
    j = 0;
    intVal += F(Values(i,j))*hxhy*0.25;

    i = xPanels;
    j = yPanels;
    intVal += F(Values(i,j))*hxhy*0.25;

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

    for(long k =0; k < (xPanels+1)*(yPanels+1); k++)
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
    return avgValue/((xMax-xMin)*(yMax-yMin));
}

//  Average computed using Trapezoidal approximation of the integral

double getTrapezoidalAverage() const
{
    double avgValue;
    avgValue = integralTrapezoidal();
    return avgValue/((xMax-xMin)*(yMax-yMin));
}

double min() const
{
//  Compute function minimum

    double minValue  = dataPtr[0];
    for(long k = 1; k < (xPanels+1)*(yPanels+1); k++)
    {
    minValue = ( minValue < dataPtr[k] ) ? minValue : dataPtr[k];
    }
    return minValue;
}

double max() const
{
//  Compute function maximum

    double maxValue  = dataPtr[0];
    for(long k = 1; k < (xPanels+1)*(yPanels+1); k++)
    {
    maxValue = ( maxValue > dataPtr[k] ) ? maxValue : dataPtr[k];
    }
    return maxValue;
}

GridFunction1d getConstantYslice(long yIndex) const  // (x function)
{
	GridFunction1d R(xPanels,xMin,xMax);
	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = Values(i,yIndex);
	}
	return R;
}

GridFunction1d getConstantXslice(long xIndex) const  //( y function)
{
	GridFunction1d R(yPanels,yMin,yMax);
	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = Values(xIndex,j);
	}
	return R;
}

void setBoundaryValues(double value)
{
    long i; long j;

	i = 0;
	for(j = 0; j <= yPanels; j++)
	{
		Values(i,j) = value;
	}

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
     Values(i,j) = value;
    }

    j = 0;
    for(i = 0; i <= xPanels; i++)
    {
    Values(i,j) = value;
    }

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    Values(i,j) = value;
    }
}

void enforcePeriodicity()
{
    long i; long j;

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
     Values(i,j) = Values(0,j);
    }

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    Values(i,j) = Values(i,0);
    }
}

//  Returns true if the input grid function is structurally identical,
//  e.g. a coincident computational domain and mesh size
//

bool isCoincident(const GridFunction2d& V)
{
	double domainDiffTol    = 1.0e-13;

    long discretizationDiff = abs(xPanels - V.xPanels)
                            + abs(yPanels - V.yPanels);

    double domainDiff = (abs(xMin - V.xMin) + abs(xMax - V.xMax))/abs(xMax-xMin)
                      + (abs(yMin - V.yMin) + abs(yMax - V.yMax))/abs(yMax-yMin);

    if((discretizationDiff > 0) || (domainDiff > domainDiffTol)) return false;
    return true;
}

//  Grid Geometry

double xMin; double xMax;   // Computational Region is [xMin,xMax]x
double yMin; double yMax;   //                         [yMin,yMax]

long   xPanels;             // Number of panels
long   yPanels;


double hx;                  // mesh width
double hy;


//###################################################################
//                     Error Checking
//###################################################################

protected:

#ifdef _DEBUG
        bool domainCheck() const
        {
        cerr << "XXX UCLAQ::GridFunction2d Error  XXX" << endl;
        cerr << "Left side of assignment must be a non-null GridFunction2d instance." << endl;
        cerr << endl;
        return false;
        }
#else
        bool domainCheck() const {return true;}
#endif

};

}


#endif /* UCLAQ_GRIDFUNCTION_UCLAQ_GRIDFUNCTION1D_H_ */
