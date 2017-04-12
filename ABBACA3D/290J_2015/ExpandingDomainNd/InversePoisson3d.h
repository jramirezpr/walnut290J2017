
//###########################################################
//                InversePoisson3d
//###########################################################
//
// This class uses the Fourier X Z procedure to compute
// the evaluate the infinite domain solution to
//
// [ laplaceCoeff (d^2/dx^2) + laplaceCoeff (d^2/dy^2) + laplaceCoeff (d^2/dz^2) +   screenCoeff ]  u
//
// in the region [xMin,xMax] x [yMin,yMax] x [zMin, zMax]
//
// The product of laplaceCoeff and screenCoeff must be less than or equal to 0.
//
// Chris Anderson, July 12, 2015
// (C) UCLA
//
//
/*
#############################################################################
#
# Copyright 2015 Chris Anderson
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
#include "PoissonNpole3d.h"
#include "ScreenedNpole3d.h"
#include "InversePoisson1d.h"
#include "InversePoisson2d.h"

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_2d.h"
#include "FFTW3_InterfaceNd/UCLAQ_FFT_Nvalues.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _CLOCKIT_
#include "ClockIt.h"
#endif

#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>
using namespace std;

#ifndef _InversePoisson3d_
#define _InversePoisson3d_

#define _DEFAULT_EXTENSION_FACTOR_     2.0
#define _DEFAULT_MAX_NPOLE_ORDER_      2
#define _DEFAULT_DIFFERENTIABILITY_     6
#define _DEFAULT_SCREEN_BOUND_         10.0

class InversePoisson3d
{
    public:

	InversePoisson3d()
	{
	initialize();
	}

	InversePoisson3d(double laplaceCoeff, double screenCoeff, long xPanel, double xMin, double xMax,
	long yPanel, double yMin, double yMax, long zPanel, double zMin, double zMax, double extensionFactor =-1.0)
	{
	this->nPoleDiffOrder  =  _DEFAULT_DIFFERENTIABILITY_;
    this->maxNpoleOrder   =  _DEFAULT_MAX_NPOLE_ORDER_;

	initialize(laplaceCoeff, screenCoeff,xPanel, xMin, xMax,yPanel,yMin,yMax,zPanel,zMin,zMax,extensionFactor);
	}

    InversePoisson3d(double laplaceCoeff, long xPanel, double xMin, double xMax,
	long yPanel, double yMin, double yMax, long zPanel, double zMin, double zMax, double extensionFactor =-1.0)
	{
    this->nPoleDiffOrder  =  _DEFAULT_DIFFERENTIABILITY_;
    this->maxNpoleOrder   =  _DEFAULT_MAX_NPOLE_ORDER_;

	initialize(laplaceCoeff,0.0,xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax,extensionFactor);
	}

    InversePoisson3d (const InversePoisson3d& H)
    {
    initialize(H);
    }

    void setMaxNpoleOrder(int maxNpoleOrder)
    {
    this->maxNpoleOrder = maxNpoleOrder;
	if(abs(screenCoeff) < 1.0e-10)
	{
			nPole.initialize(xCent,yCent,zCent,rBar,maxNpoleOrder,laplaceCoeff);
	}
	else if(abs(screenCoeff) < screenCoeffBound)
	{
		    SnPole.initialize(xCent,yCent,zCent,rBar,maxNpoleOrder,laplaceCoeff,screenCoeff);
	}

    }

    void setNpoleDiffOrder(long nPoleDiffOrder)
    {
    	this->nPoleDiffOrder  = nPoleDiffOrder;
    	if(abs(screenCoeff) < 1.0e-10)
    	{
			nPole.setDifferentiability(nPoleDiffOrder);
		}
		else if(abs(screenCoeff) < screenCoeffBound)
		{
			SnPole.setSourceDifferentiability(nPoleDiffOrder);
		}
    }


    void setScreenCoeffBound(double screenCoeffBound)
    {
    this->screenCoeffBound = screenCoeffBound;
    }


	void initialize()
	{

	// Set default values of parameters that are unlikely to
	// need to be changed. However, to change the values for a particular
	// instance, the mutators are called after the constructor, but
	// before the non-trivial initialize(*) member function.

    this->nPoleDiffOrder  =  _DEFAULT_DIFFERENTIABILITY_;
    this->maxNpoleOrder   =  _DEFAULT_MAX_NPOLE_ORDER_;


    this->nx    = 0;   this->ny   = 0;   this->nz   = 0;
    this->xMin  = 0.0; this->yMin = 0.0; this->zMin = 0.0;
    this->xMax  = 0.0; this->yMax = 0.0; this->zMax = 0.0;

    this->laplaceCoeff     = 0.0;
	this->screenCoeff      = 0.0;
    this->screenCoeffBound = _DEFAULT_SCREEN_BOUND_;


	this->xCent = 0.0;
	this->yCent = 0.0;
	this->zCent = 0.0;
	this->rBar  = 0.0;

	this->extensionFactor =  _DEFAULT_EXTENSION_FACTOR_;


	#ifdef _OPENMP
    int ompSize = DFT1dArray.size();
    for(long k = 0; k < ompSize; k++)
    {
    	DFT1dArray[k].initialize();
    	inRealArray1D[k].initialize();
    	inImagArray1D[k].initialize();
    	outRealArray1D[k].initialize();
    	outImagArray1D[k].initialize();
    }

    ompSize = invPoissonOp2dArray.size();
    for(long k = 0; k < ompSize; k++)
    {
        vTransformArray2D[k].initialize();
        invPoissonOp2dArray[k].initialize();
    }
    realData1Dx2D.initialize();
    imagData1Dx2D.initialize();
	#else
    DFT_1d.initialize();
    inReal1Dx.initialize();
    inImag1Dx.initialize();
    outReal1Dx.initialize();
    outImag1Dx.initialize();

    realData1Dx2D.initialize();
    imagData1Dx2D.initialize();
	#endif

    nPole.initialize();
    SnPole.initialize();
	nPoleSource.initialize();
   	nPolePotential.initialize();
	}

    void initialize(double laplaceCoeff, long xPanel, double xMin, double xMax,
	long yPanel, double yMin, double yMax, long zPanel, double zMin, double zMax, double extFactor = -1.0)
	{
	initialize(laplaceCoeff,0.0,xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax,extFactor);
	}

	void initialize(double laplaceCoeff, double screenCoeff, long xPanel, double xMin, double xMax,
	long yPanel, double yMin, double yMax, long zPanel, double zMin, double zMax, double extFactor = -1.0)
	{
	//
	//  A fatal error if the signs of screenCoeff and laplaceCoeff are the same.
	//
	if(abs(screenCoeff) > 1.0e-10)
	{
	if(laplaceCoeff*screenCoeff > 0)
	{
    printf("XXXX Error :  InversePoisson3d  XXXX \n");
    printf("Coefficients of the Helmholtz operator are of the same sign. \n\n");
    printf("XXXX Execution Halted XXXX\n");
    exit(0);
    }}

    this->nx    = xPanel;   this->ny   = yPanel; this->nz   = zPanel;
    this->xMin  =  xMin;    this->yMin = yMin;   this->zMin = zMin;
    this->xMax  =  xMax;    this->yMax = yMax;   this->zMax = zMax;

	this->laplaceCoeff     =  laplaceCoeff;
	this->screenCoeff      =  screenCoeff;
	this->screenCoeffBound = _DEFAULT_SCREEN_BOUND_;

	if(extFactor  < 0)  {extensionFactor = _DEFAULT_EXTENSION_FACTOR_;}
	else                {extensionFactor = extFactor;}

	//  Determine extended domain and panel counts

    double hx = (xMax-xMin)/(double)(nx);
    double hy = (yMax-yMin)/(double)(ny);
    double hz = (zMax-zMin)/(double)(nz);

    UCLAQ::FFT_Nvalues fftNvalues;

	// Extend the domain in the x direction only

    extSizeX = (xMax-xMin)*extensionFactor;
    extSizeY = (yMax-yMin);
    extNx    = extSizeX/hx;
    extNy    = ny;
    if(extNx == nx) extNx += 2;

    extNx      = fftNvalues.getFFT_N(extNx);
    extXoffset = (extNx-nx)/2;


    // Diagnostic output for debugging
    /*
    double xMinExt = xMin - extXoffset*hx;
    double xMaxExt = xMinExt + extNx*hx;
    double yMinExt = yMin;
    double yMaxExt = yMinExt + extNy*hy;

    printf("[ nx, ny, nz ] : [ %ld , %ld, %ld  ] \n",nx,ny,nz);
    printf("[ extNx, extNy, extNz ] : [ %ld , %ld, %ld ] \n",extNx, extNy, nz);
    printf("[ Xmin, Xmax ] = [ %10.5f , %10.5f ] \n",xMin,xMax);
    printf("[ Ymin, Ymax ] = [ %10.5f , %10.5f ] \n",yMin,yMax);
    printf("[ Zmin, Zmax ] = [ %10.5f , %10.5f ] \n",zMin,zMax);
    printf("[ extXmin, extXmax ] = [ %10.5f , %10.5f ] \n",xMinExt,xMaxExt);
    printf("[ extYmin, extYmax ] = [ %10.5f , %10.5f ] \n",yMinExt,yMaxExt);
    printf("[ extYmin, extYmax ] = [ %10.5f , %10.5f ] \n",zMin,zMax);
    */

    #ifdef _OPENMP
	int threadCount = omp_get_max_threads();

	inRealArray1D.resize(threadCount);
	inImagArray1D.resize(threadCount);
	outRealArray1D.resize(threadCount);
	outImagArray1D.resize(threadCount);

    DFT1dArray.resize(threadCount);
    invPoissonOp2dArray.resize(threadCount);
    vTransformArray2D.resize(threadCount);

    for(long k = 0; k < threadCount; k++)
    {
	inRealArray1D[k].initialize(extNx);
	inImagArray1D[k].initialize(extNx);
	outRealArray1D[k].initialize(extNx);
	outImagArray1D[k].initialize(extNx);
	DFT1dArray[k].initialize(extNx);

	vTransformArray2D[k].initialize(ny,yMin,yMax,nz,zMin,zMax);
	vTransformArray2D[k].setToValue(0.0);
	invPoissonOp2dArray[k].initialize(laplaceCoeff,0.0,ny,yMin,yMax,nz,zMin,zMax,extensionFactor);
    }
    realData1Dx2D.initialize((long)(extNx/2) + 1 ,ny+1,nz+1);
    imagData1Dx2D.initialize((long)(extNx/2) + 1 ,ny+1,nz+1);
    #else
    realData1Dx2D.initialize((long)(extNx/2) + 1 ,ny+1,nz+1);
    imagData1Dx2D.initialize((long)(extNx/2) + 1 ,ny+1,nz+1);

    DFT_1d.initialize(extNx);

	inReal1Dx.initialize(extNx);
    inImag1Dx.initialize(extNx);

	outReal1Dx.initialize(extNx);
    outImag1Dx.initialize(extNx);

    vTransform2D.initialize(ny,yMin,yMax,nz,zMin,zMax);
    vTransform2D.setToValue(0.0);
    inversePoissonOp2d.initialize(laplaceCoeff,0.0,ny,yMin,yMax,nz,zMin,zMax,extensionFactor);
    #endif

    // Initialize N-pole data

    xCent      = xMin + (xMax-xMin)/2.0;
    yCent      = yMin + (yMax-yMin)/2.0;
    zCent      = zMin + (zMax-zMin)/2.0;

    long   xCentIndex = ((xCent-xMin) + hx/2.0)/hx;
    long   yCentIndex = ((yCent-yMin) + hy/2.0)/hy;
    long   zCentIndex = ((zCent-zMin) + hz/2.0)/hz;

    xCent = xMin + xCentIndex*hx;
    yCent = yMin + yCentIndex*hy;
    zCent = zMin + zCentIndex*hz;

    // Determine the maximal radius for the central monopole
    // This is chosen to be the radius from the (xCent,yCent)
    // so that the outer boundary is at least two mesh panels
    // away from any boundary

    rBar = abs(xMin-xCent) - 2.0*hx;
    rBar         = (rBar > (abs(xMax-xCent) - 2.0*hx)) ? (abs(xMax-xCent) - 2.0*hx) : rBar;
    rBar         = (rBar > (abs(yMin-yCent) - 2.0*hy)) ? (abs(yMin-yCent) - 2.0*hy) : rBar;
    rBar         = (rBar > (abs(yMax-yCent) - 2.0*hy)) ? (abs(yMax-yCent) - 2.0*hy) : rBar;
    rBar         = (rBar > (abs(zMin-zCent) - 2.0*hz)) ? (abs(zMin-zCent) - 2.0*hz) : rBar;
    rBar         = (rBar > (abs(zMax-zCent) - 2.0*hz)) ? (abs(zMax-zCent) - 2.0*hz) : rBar;

    // Diagnostic output for debugging
    /*
    printf("[ nx, ny, nz] : [ %ld , %ld , %ld ] \n",nx,ny,nz);
    printf("[ Xmin, Xmax ] = [ %10.5f , %10.5f ] \n",xMin,xMax);
    printf("[ Ymin, Ymax ] = [ %10.5f , %10.5f ] \n",yMin,yMax);
    printf("[ Zmin, Zmax ] = [ %10.5f , %10.5f ] \n",zMin,zMax);

    printf("[ xC , yC , zC ] = [ %10.5f , %10.5f , %10.5f ]\n",xCent,yCent,zCent);
    printf("rBar = %10.5f \n",rBar);
    */

    if(abs(screenCoeff) < 1.0e-10)
    {
	nPole.initialize(xCent,yCent,zCent,rBar,maxNpoleOrder,laplaceCoeff);
	nPole.setDifferentiability(nPoleDiffOrder);
	}
	else if(abs(screenCoeff) < screenCoeffBound)
	{
    SnPole.initialize(xCent,yCent,zCent,rBar,maxNpoleOrder,laplaceCoeff,screenCoeff);
	SnPole.setSourceDifferentiability(nPoleDiffOrder);
	}

	nPoleSource.initialize(nx,xMin,xMax,ny,yMin,yMax,nz,zMin,zMax);
   	nPolePotential.initialize(nx,xMin,xMax,ny,yMin,yMax,nz,zMin,zMax);
	}

	void initialize(const InversePoisson3d& H)
	{
	if(this->nx == 0) initialize();

	this->nPoleDiffOrder  = H.nPoleDiffOrder;
	this->maxNpoleOrder   = H.maxNpoleOrder;

	initialize(H.laplaceCoeff, H.screenCoeff, H.nx, H.xMin, H.xMax, H.ny, H.yMin, H.yMax, H.nz, H.zMin, H.zMax, H.extensionFactor);
	this->screenCoeffBound = H.screenCoeffBound;
	}


	void setCoefficients(double laplaceCoeff, double screenCoeff)
	{
	this->laplaceCoeff  = laplaceCoeff;
	this->screenCoeff   = screenCoeff;

	if(abs(screenCoeff) < 1.0e-10)
    {
	nPole.initialize(xCent,yCent,zCent,rBar,maxNpoleOrder,laplaceCoeff);
	nPole.setDifferentiability(nPoleDiffOrder);
	}
	else if(abs(screenCoeff) < screenCoeffBound)
	{
    SnPole.initialize(xCent,yCent,zCent,rBar,maxNpoleOrder,laplaceCoeff,screenCoeff);
	SnPole.setSourceDifferentiability(nPoleDiffOrder);
	}
	}

    void applyInverseOp(UCLAQ::GridFunction3d& V)
    {
	double hx = (xMax-xMin)/(double)(nx);


#ifdef _CLOCKIT_
	correctionTime = 0.0;
	fftTime        = 0.0;
	copyTime       = 0.0;
#endif


    // Diagnostic output for debugging
    /*
    vector<double>                            moments;
	nPole.getMoments3d(xCent, yCent, zCent, V, moments);

	cout << "Matching Moments " << endl;
	for(long i = 0; i < moments.size(); i++)
	{
	cout << i << " : " << moments[i] << endl;
	}
	cout << endl << endl;
    */

	// Create moment matching function


    #ifdef _CLOCKIT_
	clockIt.start();
	#endif

    if(abs(screenCoeff) < 1.0e-10)
    {
	nPole.createMomentMatchedNpole(V);
    nPole.evaluateSource(nPoleSource);
    nPole.evaluatePotential(nPolePotential);

	V -= nPoleSource;
	}
	else if(abs(screenCoeff) < screenCoeffBound)
    {
    SnPole.createMomentMatchedNpole(V);
    SnPole.evaluateSource(nPoleSource);
    SnPole.evaluatePotential(nPolePotential);

	V -= nPoleSource;
    }

    #ifdef _CLOCKIT_
	clockIt.stop();
	correctionTime = clockIt.getMilliSecElapsedTime();
	clockIt.start();
	#endif


    // Diagnostic output for debugging
    /*

	cout << " Modified V norm " << V.norm2() << " " << V.normInf() << endl;
    nPole.getMoments3d(xCent, yCent, zCent, V, moments);
	cout << endl;
	for(long i = 0; i < moments.size(); i++)
	{
	cout << i << " : " << moments[i] << endl;
	}
	cout << endl << endl;
	*/


    int threadIndex;
    long i; long j; long p;

    double dx = 1.0/double(extNx);
    long k1;
    long k1Index;
    long kConjIndex;
    double opTransformFactor;
    double pi = .3141592653589793e+01;

	#ifdef _OPENMP
#pragma omp parallel for  \
		private(i,j,p,threadIndex)\
		schedule(static,1)
    for(j = 0; j <= ny; j++)
    {
    	threadIndex = omp_get_thread_num();
    	for(p = 0; p <= nz; p++)
		{
   		inRealArray1D[threadIndex].setToValue(0.0);
   		inImagArray1D[threadIndex].setToValue(0.0);

		for(i = 0; i <= nx; i++)
    	{
    	inRealArray1D[threadIndex](i + extXoffset) = V.Values(i,j,p);
    	}

		DFT1dArray[threadIndex].fftw1d_forward(inRealArray1D[threadIndex],inImagArray1D[threadIndex],
		                                     outRealArray1D[threadIndex],outImagArray1D[threadIndex]);

		for(i = 0; i < (long)(extNx/2) +1 ; i++)
		{
	    realData1Dx2D(i,j,p) = outRealArray1D[threadIndex](i);
        imagData1Dx2D(i,j,p) = outImagArray1D[threadIndex](i);
    	}
    	}
    }
#else
	for(long j = 0; j <= ny; j++)
    {
    for(long p = 0; p <= nz; p++)
	{
	inImag1Dx.setToValue(0.0);
	inReal1Dx.setToValue(0.0);

	for(long i = 0; i <= nx; i++)
    {
    inReal1Dx(i + extXoffset) = V.Values(i,j,p);
    }

    DFT_1d.fftw1d_forward(inReal1Dx, inImag1Dx, outReal1Dx, outImag1Dx);

	for(long i = 0; i < (long)(extNx/2) + 1; i++)
	{
	realData1Dx2D(i,j,p) = outReal1Dx(i);
    imagData1Dx2D(i,j,p) = outImag1Dx(i);
	}
    }}
#endif


#ifdef _OPENMP
#pragma omp parallel for  \
	private(j,p,k1,k1Index,opTransformFactor,threadIndex)\
	schedule(static,1)
    for(k1 = -(extNx/2); k1 <= 0; k1++)
    {
		threadIndex = omp_get_thread_num();
		k1Index           =  k1 + (extNx/2);
    	opTransformFactor = -laplaceCoeff*(((2.0*pi*k1*dx)*(2.0*pi*k1*dx))/(hx*hx)) + screenCoeff;
        invPoissonOp2dArray[threadIndex].setCoefficients(laplaceCoeff,opTransformFactor);

    	for(j = 0; j <= ny; j++)
    	{
    	for(p = 0; p <= nz; p++)
    	{
    	vTransformArray2D[threadIndex].Values(j,p) = realData1Dx2D(k1Index,j,p);
   		}}

        invPoissonOp2dArray[threadIndex].applyInverseOp(vTransformArray2D[threadIndex]);

    	for(j = 0; j <= ny; j++)
    	{
    	for(p = 0; p <= nz; p++)
    	{
        realData1Dx2D(k1Index,j,p) = vTransformArray2D[threadIndex].Values(j,p);
    	}}

    	for(j = 0; j <= ny; j++)
    	{
    	for(p = 0; p <= nz; p++)
    	{
        vTransformArray2D[threadIndex].Values(j,p) = imagData1Dx2D(k1Index,j,p);
   		}}

        invPoissonOp2dArray[threadIndex].applyInverseOp(vTransformArray2D[threadIndex]);

    	for(j = 0; j <= ny; j++)
    	{
    	for(p = 0; p <= nz; p++)
    	{
        imagData1Dx2D(k1Index,j,p) = vTransformArray2D[threadIndex].Values(j,p);
    	}}
 }
 #else

     // for(k1 = -(extNx/2); k1 <= (extNx-1)/2; k1++)

    for(k1 = -(extNx/2); k1 <= 0; k1++)
    {
    k1Index           =  k1 + (extNx/2);
    opTransformFactor = -laplaceCoeff*(((2.0*pi*k1*dx)*(2.0*pi*k1*dx))/(hx*hx)) + screenCoeff;
    inversePoissonOp2d.setCoefficients(laplaceCoeff,opTransformFactor);

    // Real component

    for(long j = 0; j <= ny; j++)
    {
    for(long p = 0; p <= nz; p++)
    {
    vTransform2D.Values(j,p) = realData1Dx2D(k1Index,j,p);
    }}

    inversePoissonOp2d.applyInverseOp(vTransform2D);

    for(long j = 0; j <= ny; j++)
    {
    for(long p = 0; p <= nz; p++)
    {
    realData1Dx2D(k1Index,j,p) = vTransform2D.Values(j,p);
    }}

    // Imag component

    for(long j = 0; j <= ny; j++)
    {
    for(long p = 0; p <= nz; p++)
    {
    vTransform2D.Values(j,p)     = imagData1Dx2D(k1Index,j,p);
    }}

    inversePoissonOp2d.applyInverseOp(vTransform2D);

    for(long j = 0; j <= ny; j++)
    {
    for(long p = 0; p <= nz; p++)
    {
    imagData1Dx2D(k1Index,j,p) = vTransform2D.Values(j,p);
    }}
}
#endif


#ifdef _OPENMP

// Transform back

#pragma omp parallel for  \
		private(i,j,p,k1,k1Index,kConjIndex,threadIndex)\
		schedule(static,1)
    for(j = 0; j <= ny; j++)
    {
    	// Each thread

    	threadIndex = omp_get_thread_num();
    	for(p = 0; p <= nz; p++)
		{
    	for(k1 = -(extNx/2); k1 <= 0; k1++)
   		{
   	    k1Index           =  k1 + (extNx/2);
   	    inRealArray1D[threadIndex](k1Index) = realData1Dx2D(k1Index,j,p);
        inImagArray1D[threadIndex](k1Index) = imagData1Dx2D(k1Index,j,p);
        }


        for(k1 = 1; k1 <= (extNx-1)/2; k1++)
    	{
    	k1Index            =  k1 + (extNx/2);
    	kConjIndex         = -k1 + (extNx/2);
 		inRealArray1D[threadIndex](k1Index) =  realData1Dx2D(kConjIndex,j,p);
   	 	inImagArray1D[threadIndex](k1Index) = -imagData1Dx2D(kConjIndex,j,p);
    	}

 		DFT1dArray[threadIndex].fftw1d_inverse(inRealArray1D[threadIndex],inImagArray1D[threadIndex],
		                                      outRealArray1D[threadIndex],outImagArray1D[threadIndex]);
		for(i = 0; i <= nx; i++)
    	{
    	V.Values(i,j,p) = outRealArray1D[threadIndex](i + extXoffset);
   		}
    	}
    }
#else

// Inverse transform

	for(long j = 0; j <= ny; j++)
	{
	for(long p = 0; p <= nz; p++)
	{
    for(k1 = -(extNx/2); k1 <= 0; k1++)
    {
    k1Index           =  k1 + (extNx/2);
 	inReal1Dx(k1Index) = realData1Dx2D(k1Index,j,p);
    inImag1Dx(k1Index) = imagData1Dx2D(k1Index,j,p);
    }

    for(k1 = 1; k1 <= (extNx-1)/2; k1++)
    {
    k1Index            =  k1 + (extNx/2);
    kConjIndex         = -k1 + (extNx/2);
 	inReal1Dx(k1Index) =  realData1Dx2D(kConjIndex,j,p);
    inImag1Dx(k1Index) = -imagData1Dx2D(kConjIndex,j,p);
    }

    DFT_1d.fftw1d_inverse(inReal1Dx, inImag1Dx, outReal1Dx, outImag1Dx);

    for(long i = 0; i <= nx; i++)
	{
	V.Values(i,j,p) = outReal1Dx(i + extXoffset);
	}

	}

	};
#endif



    #ifdef _CLOCKIT_
	clockIt.stop();
	fftTime = clockIt.getMilliSecElapsedTime();
	clockIt.start();
	#endif

//  Add in nPole correction


    if(fabs(screenCoeff) < 1.0e-10)
    {
    V += nPolePotential;
    }
    else if(abs(screenCoeff) < screenCoeffBound)
    {
    V += nPolePotential;
    }


    #ifdef _CLOCKIT_
	clockIt.stop();
	copyTime = clockIt.getMilliSecElapsedTime();
	#endif

	}

    // Diagnostic utility routines for disabling/enabling the
    // use of a corrected 2D screened potential.
    //
	// The default implementation uses the corrected 2D screened potential
	// with a current default coefficient value of 10.0
	//

	void disable2DscreenedPotentialCorrection()
	{
	#ifdef _OPENMP
    for(long i = 0; i < (long)invPoissonOp2dArray.size(); i++)
    {
        invPoissonOp2dArray[i].setScreenCoeffBound(0.0);
    }
	#else
	inversePoissonOp2d.setScreenCoeffBound(0.0);
	#endif
	}

    void enable2DscreenedPotentialCorrection(double coeffBound = 10.0)
	{
	#ifdef _OPENMP
    for(long i = 0; i < (long)invPoissonOp2dArray.size(); i++)
    {
        invPoissonOp2dArray[i].setScreenCoeffBound(coeffBound);
    }
	#else
	inversePoissonOp2d.setScreenCoeffBound(coeffBound);
	#endif
	}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    double     extensionFactor;
	long       extXoffset;  long   extYoffset;
    double     extSizeX;    double extSizeY;
    long       extNx;       long    extNy;


	double     screenCoeff;
	double    laplaceCoeff;

	long        nx;     long ny;   long  nz;
	double    xMin; double yMin; double zMin;
	double    xMax; double yMax; double zMax;

    // 1D Fourier Transform Persistent data

    UCLAQ::DoubleVector3d         realData1Dx2D;
    UCLAQ::DoubleVector3d         imagData1Dx2D;

#ifdef _OPENMP
    // 1D Fourier Transform

   vector<UCLAQ::fftw3_1d>          DFT1dArray;
   vector<UCLAQ::DoubleVector1d>  inRealArray1D;
   vector<UCLAQ::DoubleVector1d>  inImagArray1D;
   vector<UCLAQ::DoubleVector1d> outRealArray1D;
   vector<UCLAQ::DoubleVector1d> outImagArray1D;

    // 2D Components

   vector<UCLAQ::GridFunction2d>         vTransformArray2D;
   vector<InversePoisson2d>            invPoissonOp2dArray;
#else

    // 1D Fourier Transform

    UCLAQ::fftw3_1d                   DFT_1d;
    UCLAQ::DoubleVector1d           inReal1Dx;
    UCLAQ::DoubleVector1d           inImag1Dx;
    UCLAQ::DoubleVector1d          outReal1Dx;
    UCLAQ::DoubleVector1d          outImag1Dx;

     // 2D Components

    UCLAQ::GridFunction2d       vTransform2D;
    InversePoisson2d      inversePoissonOp2d;

#endif

    // N-pole and Screened N-pole data

    PoissonNpole3d                nPole;
    ScreenedNpole3d              SnPole;
   	UCLAQ::GridFunction3d   nPoleSource;
   	UCLAQ::GridFunction3d nPolePotential;

   	double xCent; double yCent; double zCent; double rBar;

   	int        maxNpoleOrder;
   	int       nPoleDiffOrder;
   	double  screenCoeffBound;


   	#ifdef _CLOCKIT_
	ClockIt clockIt;
	double correctionTime;
	double fftTime;
	double copyTime;
	#endif

};

#undef _DEFAULT_EXTENSION_FACTOR_
#undef _DEFAULT_MAX_NPOLE_ORDER_
#undef _DEFAULT_DIFFERENTIABILITY_
#undef _DEFAULT_SCREEN_BOUND_
#endif

