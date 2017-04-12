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

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include <cstdlib>
using namespace std;

#ifndef _QdomainQdomainKineticOp_
#define _QdomainQdomainKineticOp_

//
// This class computes the discrete Poisson operator for functions vanishing
// outside the computation domain. For higher order difference approximations,
// the function values outside the computational domain that are required for the
// difference stencil are assumed to be identically zero.
//
// H(V) = beta*div(grad(V)) (constant coefficient case)
//
//     --or--
//
// H(V) = div(a*grad(V))
//
//    --or--
//
// H(V) = [D_x(beta_x a D_x] V + [D_y(beta_y a D_y] + [D_z(beta_z a D_z]
//
//
// Constant coefficient operator orders : 2,4,6,8
// Variable coefficient operator orders : 2,4,6
//
//
// In the variable coefficient case, the values of the coefficients are
// extended using constant value extension normal to the domain boundary.
//

namespace UCLAQ
{


class QdomainKineticOp
{
public:

	QdomainKineticOp()
	{initialize();}

	QdomainKineticOp(const QdomainKineticOp& Op)
	{initialize(Op);}

	QdomainKineticOp(int order, double beta)
	{
	initialize(order, beta);
	}

    QdomainKineticOp(int order, double beta_x,double beta_y, double beta_z)
	{
	initialize(order, beta_x, beta_y, beta_z);
	}

	QdomainKineticOp(int order, const UCLAQ::GridFunction3d& a)
	{
	initialize(order,a);
	}

    QdomainKineticOp(int order, const UCLAQ::GridFunction3d& a, double beta_x, double beta_y, double beta_z)
	{
	initialize(order,a,beta_x,beta_y,beta_z);
	}


    void initialize()
    {
	this->order  = 2;
	this->xPanel = 0;    this->yPanel = 0;    this->zPanel = 0;
	this->xPanelExt = 0; this->yPanelExt = 0; this->zPanelExt = 0;
	this->valuesTemp.initialize();
	this->aValues.initialize();
	this->beta_x =      0.0;
	this->beta_y =      0.0;
	this->beta_z =      0.0;
	this->variableCoeffFlag = false;
    }


    void initialize(const QdomainKineticOp& Op)
    {
    this->order  = Op.order;
	this->xPanel = Op.xPanel;       this->yPanel = Op.yPanel;       this->zPanel = Op.zPanel;
	this->xPanelExt = Op.xPanelExt; this->yPanelExt = Op.yPanelExt; this->zPanelExt = Op.zPanelExt;
	this->beta_x =  Op.beta_x;
	this->beta_y =  Op.beta_y;
	this->beta_z =  Op.beta_z;
	this->variableCoeffFlag = Op.variableCoeffFlag;

    if((Op.xPanel == 0)||(Op.yPanel==0)||(Op.zPanel==0))   // Appropriately initialize a null instance
    {
    return;
    }

	this->valuesTemp.initialize(Op.valuesTemp);
	this->aValues.initialize(Op.aValues);
    }

    void initialize(int order, double beta)
    {
    initialize();
	this->variableCoeffFlag = false;

    this->order  = order;
    if(this->order > 8)  {this->order = 8;}
	this->beta_x   = beta;
	this->beta_y   = beta;
	this->beta_z   = beta;
    }

    void initialize(int order, double beta_x, double beta_y, double beta_z)
    {
    initialize();
	this->variableCoeffFlag = false;

    this->order  = order;
    if(this->order > 8)  {this->order = 8;}
	this->beta_x   = beta_x;
	this->beta_y   = beta_y;
	this->beta_z   = beta_z;
    }


    void initialize(int order, const UCLAQ::GridFunction3d& a, double beta_x, double beta_y, double beta_z)
    {
    initialize(order,a);
    this->beta_x   = beta_x;
	this->beta_y   = beta_y;
	this->beta_z   = beta_z;
    }

    void initialize(int order, const UCLAQ::GridFunction3d& a)
    {
    initialize();


    xPanel = a.getXpanelCount();
    yPanel = a.getYpanelCount();
    zPanel = a.getZpanelCount();


	this->beta_x   = 1.0;
	this->beta_y   = 1.0;
	this->beta_z   = 1.0;

    this->order = order;

   	this->variableCoeffFlag = true;
    if(this->order > 6)  {this->order = 6;}

	setDimensions(xPanel,yPanel, zPanel);

	long i; long j; long k;

	long offset;

	if(order == 2) offset = 1;
	if(order == 4) offset = 2;
	if(order == 6) offset = 4;


	for(i = 0; i <= xPanel; i++)
	{
	for(j = 0; j <= yPanel; j++)
	{
	for(k = 0; k <= zPanel; k++)
	{
	aValues(i+offset,j+offset,k+offset) = a(i,j,k);
	}}}

	//
	// Constant extension normal to the boundary for exterior
	// coefficient values
	//

	// X-coordinate extension

	long indexSize = aValues.getIndex1Size();

    for(i = 0; i < offset; i++)
    {
    	for(j = 0; j <= yPanel; j++)
    	{
    	for(k = 0; k <= zPanel; k++)
    	{
    	aValues(i,j+offset,k+offset) = a(0,j,k);
    	}}
    }

    for(i = xPanel+1; i <  indexSize; i++)
    {
    	for(j = 0; j <= yPanel; j++)
    	{
    	for(k = 0; k <= zPanel; k++)
    	{
    		aValues(i,j+offset,k+offset) = a(xPanel,j,k);
    	}}
    }

    // Y-coordinate extension

    indexSize = aValues.getIndex2Size();

	for(j = 0; j <  offset; j++)
	{
		for(i = 0; i <= xPanel; i++)
		{
		for(k = 0; k <= zPanel; k++)
		{
			aValues(i+offset,j,k+offset) = a(i,0,k);
		}}
	}

	for(j = yPanel+1; j <  indexSize; j++)
	{
		for(i = 0; i <= xPanel; i++)
		{
		for(k = 0; k <= zPanel; k++)
		{
			aValues(i+offset,j,k+offset) = a(i,yPanel,k);
		}}
	}

	// Z-coordinate extension

    indexSize = aValues.getIndex3Size();

	for(k = 0; k < offset; k++)
	{
	for(i = 0; i <= xPanel; i++)
	{
	for(j = 0; j <= yPanel; j++)
	{
	aValues(i+offset,j+offset,k) = a(i,j,0);
	}}}

	for(k = zPanel+1; k < indexSize; k++)
	{
	for(i = 0; i <= xPanel; i++)
	{
	for(j = 0; j <= yPanel; j++)
	{
	aValues(i+offset,j+offset,k) = a(i,j,zPanel);
	}}}
    }

	void apply(UCLAQ::GridFunction3d& V)
	{
	setDimensions(V.getXpanelCount(),V.getYpanelCount(),V.getZpanelCount());
	//
	// Copy over data to the interior of the temp array
	//
	long i; long j; long k;

    long offset;

	if(order == 2) offset = 1;
	if(order == 4) offset = 2;
	if(order == 6) offset = 4;

	if(not variableCoeffFlag)
	{
	if(order == 8) offset = 4;
	}

	for(i = 0; i <= xPanel; i++)
	{
	for(j = 0; j <= yPanel; j++)
	{
	for(k = 0; k <= zPanel; k++)
	{
	valuesTemp(i+offset,j+offset,k+offset) = V(i,j,k);
	}}}

	if(not variableCoeffFlag)
	{
		switch(order)
		{
		case 2: applyForwardOp2thOrder(V); break;
		case 4: applyForwardOp4thOrder(V); break;
		case 6: applyForwardOp6thOrder(V); break;
		case 8: applyForwardOp8thOrder(V); break;
		}
	}
	else
	{
		switch(order)
		{
		case 2: applyForwardOp2thOrderVarCoeff(V); break;
		case 4: applyForwardOp4thOrderVarCoeff(V); break;
		case 6: applyForwardOp6thOrderVarCoeff(V); break;
		}
	}

	}


	void applyForwardOp2thOrderVarCoeff(UCLAQ::GridFunction3d& V)
	{
	double aP;
	double a0;
	double aM;

    long i; long j; long k;
    long iOff; long jOff; long kOff;

    double hx = V.getHx();
    double hy = V.getHy();
    double hz = V.getHz();

    double ddX; double ddY; double ddZ;

    for(i = 1, iOff = 0; iOff <= xPanel; i++,iOff++)
    {
    for(j = 1, jOff = 0; jOff <= yPanel; j++,jOff++)
    {
    for(k = 1, kOff = 0; kOff <= zPanel; k++,kOff++)
    {
    aM  = (aValues(i-1,j,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i+1,j,k))/2.0;
	a0  = aM + aP;
    ddX = (aP*valuesTemp(i+1,j,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i-1,j,k))/(hx*hx);

    aM  = (aValues(i,j-1,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j+1,k))/2.0;
	a0  = aM + aP;
    ddY = (aP*valuesTemp(i,j+1,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j-1,k))/(hy*hy);

    aM  = (aValues(i,j,k-1)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j,k+1))/2.0;
	a0  = aM + aP;
    ddZ = (aP*valuesTemp(i,j,k+1)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j,k-1))/(hz*hz);

    V(iOff,jOff,kOff) = beta_x*ddX + beta_y*ddY + beta_z*ddZ;
    }}}


	}
	void applyForwardOp4thOrderVarCoeff(UCLAQ::GridFunction3d& V)
	{
	double aP;
	double a0;
	double aM;

    long i; long j; long k;
    long iOff; long jOff; long kOff;

	double ddX; double ddY; double ddZ;
	double ddX2; double ddY2; double ddZ2;

    double hx = V.getHx();
    double hy = V.getHy();
    double hz = V.getHz();

    for(i = 2, iOff = 0; iOff <= xPanel; i++,iOff++)
    {
    for(j = 2, jOff = 0; jOff <= yPanel; j++,jOff++)
    {
    for(k = 2, kOff = 0; kOff <= zPanel; k++,kOff++)
    {
    aM  = (aValues(i-1,j,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i+1,j,k))/2.0;
	a0  = aM + aP;
    ddX = (aP*valuesTemp(i+1,j,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i-1,j,k))/(hx*hx);

    aM  = (aValues(i,j-1,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j+1,k))/2.0;
	a0  = aM + aP;
    ddY = (aP*valuesTemp(i,j+1,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j-1,k))/(hy*hy);

    aM  = (aValues(i,j,k-1)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j,k+1))/2.0;
	a0  = aM + aP;
    ddZ = (aP*valuesTemp(i,j,k+1)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j,k-1))/(hz*hz);

    //
    // 2*h stencil
    //

    aM  = (aValues(i-2,j,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i+2,j,k))/2.0;
	a0  = aM + aP;
    ddX2 = (aP*valuesTemp(i+2,j,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i-2,j,k))/(4.0*hx*hx);

    aM  = (aValues(i,j-2,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j+2,k))/2.0;
	a0  = aM + aP;
    ddY2 = (aP*valuesTemp(i,j+2,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j-2,k))/(4.0*hy*hy);

    aM  = (aValues(i,j,k-2)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j,k+2))/2.0;
	a0  = aM + aP;
    ddZ2 = (aP*valuesTemp(i,j,k+2)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j,k-2))/(4.0*hz*hz);


    V(iOff,jOff,kOff) = beta_x*((4.0/3.0)*ddX - (1.0/3.0)*ddX2)
    		                 + beta_y*((4.0/3.0)*ddY - (1.0/3.0)*ddY2)
    		                 + beta_z*((4.0/3.0)*ddZ - (1.0/3.0)*ddZ2);
    }}}

	}

	void applyForwardOp6thOrderVarCoeff(UCLAQ::GridFunction3d& V)
	{
	double aP;
	double a0;
	double aM;

    long i; long j; long k;
    long iOff; long jOff; long kOff;

	double ddX;  double ddY;  double ddZ;
	double ddX2; double ddY2; double ddZ2;
	double ddX4; double ddY4; double ddZ4;

	double hx = V.getHx();
    double hy = V.getHy();
    double hz = V.getHz();


    for(i = 4, iOff = 0; iOff <= xPanel; i++,iOff++)
    {
    for(j = 4, jOff = 0; jOff <= yPanel; j++,jOff++)
    {
    for(k = 4, kOff = 0; kOff <= zPanel; k++,kOff++)
    {
    aM  = (aValues(i-1,j,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i+1,j,k))/2.0;
	a0  = aM + aP;
    ddX = (aP*valuesTemp(i+1,j,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i-1,j,k))/(hx*hx);

    aM  = (aValues(i,j-1,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j+1,k))/2.0;
	a0  = aM + aP;
    ddY = (aP*valuesTemp(i,j+1,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j-1,k))/(hy*hy);

    aM  = (aValues(i,j,k-1)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j,k+1))/2.0;
	a0  = aM + aP;
    ddZ = (aP*valuesTemp(i,j,k+1)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j,k-1))/(hz*hz);

    //
    // 2*h stencil
    //

    aM  = (aValues(i-2,j,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i+2,j,k))/2.0;
	a0  = aM + aP;
    ddX2 = (aP*valuesTemp(i+2,j,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i-2,j,k))/(4.0*hx*hx);

    aM  = (aValues(i,j-2,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j+2,k))/2.0;
	a0  = aM + aP;
    ddY2 = (aP*valuesTemp(i,j+2,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j-2,k))/(4.0*hy*hy);

    aM  = (aValues(i,j,k-2)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j,k+2))/2.0;
	a0  = aM + aP;
    ddZ2 = (aP*valuesTemp(i,j,k+2)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j,k-2))/(4.0*hz*hz);

    //
    // 4*h stencil
    //

    aM  = (aValues(i-4,j,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i+4,j,k))/2.0;
	a0  = aM + aP;
    ddX4 = (aP*valuesTemp(i+4,j,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i-4,j,k))/(16.0*hx*hx);

    aM  = (aValues(i,j-4,k)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j+4,k))/2.0;
	a0  = aM + aP;
    ddY4 = (aP*valuesTemp(i,j+4,k)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j-4,k))/(16.0*hy*hy);

    aM  = (aValues(i,j,k-4)+ aValues(i,j,k))/2.0;
	aP  = (aValues(i,j,k)  + aValues(i,j,k+4))/2.0;
	a0  = aM + aP;
    ddZ4 = (aP*valuesTemp(i,j,k+4)-a0*valuesTemp(i,j,k) + aM*valuesTemp(i,j,k-4))/(16.0*hz*hz);

    V(iOff,jOff,kOff) =

    beta_x*((16.0/15.0)*((4.0/3.0)*ddX - (1.0/3.0)*ddX2) - (1.0/15.0)*((4.0/3.0)*ddX2 - (1.0/3.0)*ddX4))
    +
    beta_y*((16.0/15.0)*((4.0/3.0)*ddY - (1.0/3.0)*ddY2) - (1.0/15.0)*((4.0/3.0)*ddY2 - (1.0/3.0)*ddY4))
    +
    beta_z*((16.0/15.0)*((4.0/3.0)*ddZ - (1.0/3.0)*ddZ2) - (1.0/15.0)*((4.0/3.0)*ddZ2 - (1.0/3.0)*ddZ4));
    }}}

	}
	//
	// These routines assume that the interior values of valuesTemp have been
	// initialized with the values of V
	//
	// For the 6th order operator, the offset used is larger than that actually
	// required by the stencil, but the larger offset is used to be consistent
	// with the offset required by the 6th order variable coefficient operator.
	//

	void applyForwardOp2thOrder(UCLAQ::GridFunction3d& V)
	{
    long i; long j; long k;
    long iOff; long jOff; long kOff;

    double hx = V.getHx();
    double hy = V.getHy();
    double hz = V.getHz();

	long offset = 1;

    for(i = offset, iOff = 0; iOff <= xPanel; i++,iOff++)
    {
    for(j = offset, jOff = 0; jOff <= yPanel; j++,jOff++)
    {
    for(k = offset, kOff = 0; kOff <= zPanel; k++,kOff++)
    {
    V(iOff,jOff,kOff) =
    beta_x*((valuesTemp(i+1,j,k)-2.0*valuesTemp(i,j,k) + valuesTemp(i-1,j,k))/(hx*hx))
    +
    beta_y*((valuesTemp(i,j+1,k)-2.0*valuesTemp(i,j,k) + valuesTemp(i,j-1,k))/(hy*hy))
    +
    beta_z*((valuesTemp(i,j,k+1)-2.0*valuesTemp(i,j,k) + valuesTemp(i,j,k-1))/(hz*hz));
    }}}
	}

	void applyForwardOp4thOrder(UCLAQ::GridFunction3d& V)
	{
    long i; long j; long k;
    long iOff; long jOff; long kOff;

    double hx = V.getHx();
    double hy = V.getHy();
    double hz = V.getHz();

	long offset = 2;

    for(i = offset, iOff = 0; iOff <= xPanel; i++,iOff++)
    {
    for(j = offset, jOff = 0; jOff <= yPanel; j++,jOff++)
    {
    for(k = offset, kOff = 0; kOff <= zPanel; k++,kOff++)
    {
	V(iOff,jOff,kOff) =
    beta_x*((-valuesTemp(i+2,j,k)+ 16.0*valuesTemp(i+1,j,k)-30.0*valuesTemp(i,j,k) + 16.0*valuesTemp(i-1,j,k) - valuesTemp(i-2,j,k))/(12.0*hx*hx))
    +
	beta_y*((-valuesTemp(i,j+2,k)+ 16.0*valuesTemp(i,j+1,k)-30.0*valuesTemp(i,j,k) + 16.0*valuesTemp(i,j-1,k) - valuesTemp(i,j-2,k))/(12.0*hy*hy))
    +
    beta_z*((-valuesTemp(i,j,k+2)+ 16.0*valuesTemp(i,j,k+1)-30.0*valuesTemp(i,j,k) + 16.0*valuesTemp(i,j,k-1) - valuesTemp(i,j,k-2))/(12.0*hz*hz));
    }}}
	}

	void applyForwardOp6thOrder(UCLAQ::GridFunction3d& V)
	{
    long i; long j; long k;
    long iOff; long jOff; long kOff;

    double hx = V.getHx();
    double hy = V.getHy();
    double hz = V.getHz();

	long offset = 4;

    for(i = offset, iOff = 0; iOff <= xPanel; i++,iOff++)
    {
    for(j = offset, jOff = 0; jOff <= yPanel; j++,jOff++)
    {
    for(k = offset, kOff = 0; kOff <= zPanel; k++,kOff++)
    {
    V(iOff,jOff,kOff) =
    beta_x*((2.0*valuesTemp(i+3,j,k) -27.0*valuesTemp(i+2,j,k)+ 270.0*valuesTemp(i+1,j,k)-490.0*valuesTemp(i,j,k) + 270.0*valuesTemp(i-1,j,k) - 27.0*valuesTemp(i-2,j,k) + 2.0*valuesTemp(i-3,j,k))/(180.0*hx*hx))
    +
	beta_y*((2.0*valuesTemp(i,j+3,k) -27.0*valuesTemp(i,j+2,k)+ 270.0*valuesTemp(i,j+1,k)-490.0*valuesTemp(i,j,k) + 270.0*valuesTemp(i,j-1,k) - 27.0*valuesTemp(i,j-2,k) + 2.0*valuesTemp(i,j-3,k))/(180.0*hy*hy))
	+
	beta_x*((2.0*valuesTemp(i,j,k+3) -27.0*valuesTemp(i,j,k+2)+ 270.0*valuesTemp(i,j,k+1)-490.0*valuesTemp(i,j,k) + 270.0*valuesTemp(i,j,k-1) - 27.0*valuesTemp(i,j,k-2) + 2.0*valuesTemp(i,j,k-3))/(180.0*hz*hz));
    }}}
	}


	void applyForwardOp8thOrder(UCLAQ::GridFunction3d& V)
	{
	long i; long j; long k;
	long iOff; long jOff; long kOff;

	double hx = V.getHx();
	double hy = V.getHy();
	double hz = V.getHz();

	for(i = 4, iOff = 0; iOff <= xPanel; i++,iOff++)
	{
	for(j = 4, jOff = 0; jOff <= yPanel; j++,jOff++)
	{
	for(k = 4, kOff = 0; kOff <= zPanel; k++,kOff++)
	{
	V(iOff,jOff,kOff) =
	beta_x*((-9.0*valuesTemp(i+4,j,k) + 128.0*valuesTemp(i+3,j,k) -1008.0*valuesTemp(i+2,j,k)+ 8064.0*valuesTemp(i+1,j,k)-14350.0*valuesTemp(i,j,k)	+  8064.0*valuesTemp(i-1,j,k) - 1008.0*valuesTemp(i-2,j,k) + 128.0*valuesTemp(i-3,j,k) - 9.0*valuesTemp(i-4,j,k))/(5040.0*hx*hx))
	+
	beta_y*((-9.0*valuesTemp(i,j+4,k) + 128.0*valuesTemp(i,j+3,k) -1008.0*valuesTemp(i,j+2,k)+ 8064.0*valuesTemp(i,j+1,k)-14350.0*valuesTemp(i,j,k)	+  8064.0*valuesTemp(i,j-1,k) - 1008.0*valuesTemp(i,j-2,k) + 128.0*valuesTemp(i,j-3,k) - 9.0*valuesTemp(i,j-4,k))/(5040.0*hy*hy))
	+
	beta_z*((-9.0*valuesTemp(i,j,k+4) + 128.0*valuesTemp(i,j,k+3) -1008.0*valuesTemp(i,j,k+2)+ 8064.0*valuesTemp(i,j,k+1)-14350.0*valuesTemp(i,j,k)	+  8064.0*valuesTemp(i,j,k-1) - 1008.0*valuesTemp(i,j,k-2) + 128.0*valuesTemp(i,j,k-3) - 9.0*valuesTemp(i,j,k-4))/(5040.0*hz*hz));
	}}}
	}


	// order == 2 ==> offset = 1 => panelCount + 2
	// order == 4 ==> offset = 2 => panelCount + 4
	// order == 6 ==> offset = 4 => panelCount + 8

	void setDimensions(long xPanelInput, long yPanelInput, long zPanelInput)
	{
	if((xPanelInput == 0)||(yPanelInput == 0)||(zPanelInput == 0))
	{
    return;
	}
	//
    // If dimensions have changed, then re-allocated the temporary array
    //
	long extSize;

	if(order == 2) {extSize  = 2;}
	if(order == 4) {extSize  = 4;}
    if(order == 6) {extSize  = 8;}
    if(order == 8) {extSize  = 8;} // Only if constant coefficient

    if((xPanelExt != xPanelInput + extSize)||
       (yPanelExt != yPanelInput + extSize)||
       (zPanelExt != zPanelInput + extSize))
    {
	    xPanel = xPanelInput; xPanelExt = xPanel + extSize;
	    yPanel = yPanelInput; yPanelExt = yPanel + extSize;
	    zPanel = zPanelInput; zPanelExt = zPanel + extSize;
	    valuesTemp.initialize(xPanelExt+1, yPanelExt+1, zPanelExt+1);
	    valuesTemp.setToValue(0.0);

	    if(this->variableCoeffFlag == true)
	    {
	    aValues.initialize(xPanelExt+1, yPanelExt+1, zPanelExt+1);
	    aValues.setToValue(0.0);
	    }
	    else
	    {
	    aValues.initialize();
	    }
    }


	}
//
//  Utility routine for unit testing.
//
	void forwardOpUnitTestUtility(UCLAQ::GridFunction3d& Vext, UCLAQ::GridFunction3d& V)
	{
	setDimensions(V.getXpanelCount(),V.getYpanelCount(),V.getZpanelCount());
	//
	// Copy over Vext values to values Temp, and then overwrite with V at interior
	// points to set up appropriate boundary values
	//
	valuesTemp = Vext;

	long offset;
	if(order == 2) offset = 1;
	if(order == 4) offset = 2;
	if(order == 6) offset = 4;

	long i; long j; long k;

	for(i = 0; i <= xPanel; i++)
	{
	for(j = 0; j <= yPanel; j++)
	{
	for(k = 0; k <= zPanel; k++)
	{
		valuesTemp(i+offset,j+offset,k+offset) = V(i,j,k);
	}}}

    if(not variableCoeffFlag)
    {
	switch(order)
	{
		case 2: applyForwardOp2thOrder(V); break;
		case 4: applyForwardOp4thOrder(V); break;
		case 6: applyForwardOp6thOrder(V); break;
	}
    }
    else
    {
    switch(order)
    {
	case 2: applyForwardOp2thOrderVarCoeff(V); break;
	case 4: applyForwardOp4thOrderVarCoeff(V); break;
	case 6: applyForwardOp6thOrderVarCoeff(V); break;
    }
    }
	}

	UCLAQ::DoubleVector3d valuesTemp;
    int order;
	long xPanel;    long yPanel;    long zPanel;
	long xPanelExt; long yPanelExt; long zPanelExt;

	double 			        beta_x;
	double                  beta_y;
	double                  beta_z;
	UCLAQ::DoubleVector3d  aValues;
	bool variableCoeffFlag;

};

} // namespace
#endif
