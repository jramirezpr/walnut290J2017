#include "fftw3.h"
#include "../DoubleVectorNd/UCLAQ_DoubleVector2d.h"
#include "../GridFunctionNd/UCLAQ_GridFunction2d.h"

#ifndef _UCLAQ_fftw3_2d_
#define _UCLAQ_fftw3_2d_
//
//####################################################################
//          fftw3_2d.h : Version Nov. 11, 2015
//####################################################################
//
/**<p>

    This class provides forward and inverse transforms of two dimensional
    data stored in UCLAQ::GridFunction2d instances or data stored in
    UCLAQ::DoubleVector2d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical transforms of
    periodic functions defined over intervals of size LX x LY.

    The transforms are implemented using the dft_2d(...) routines from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/). The
    FFTW3 headers and library files must be available to create any executable
    that uses this class.

    This transform interface is designed be used with the data values of
    discrete  periodic functions associated with a uniform grid discretization
    consisting of nx grid panels x ny grid panels spanning a rectangular
    domain of size LX x LY.

    When DoubleVector2d's are used for input and/or output their
    size is nx x ny; the number of independent data values associated
    with the representation of periodic functions of such discretizations.

    For computational efficiency nx and ny should be chosen to be a
    product of small primes less than or equal to 13  (if possible).

    ===============================================================

    What's computed assuming the input data is specified as a
    GridFunction2d instance initialized with the parameters
    (nx, xMin, xMax, ny, yMin, yMax):

    If (kx,ky) is a wave number vector with kx, and ky in the ranges

    -(nx/2) <= kx <=  (nx-1)/2  and -(ny/2) <= ky <=  (ny-1)/2

    (the range calculated with integer division), then two GridFunction2d
    instances whose values are the real and imaginary parts of

    d  = sqrt(1/LX)*exp(2*pi*I*kx*(x-xMin)/LX)
        *sqrt(1/LY)*exp(2*pi*I*ky*(y-yMin)/LY)

    evaluated the points

    x = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny

    get mapped by the forward transform to the two DoubleVector2d instances
    of size nx x ny, fhat_real and fhat_imag, where 

    fhat_real(kx + (nx/2), ky + (ny/2)) = 1
    fhat_imag(kx + (nx/2), ky + (ny/2)) = 1

    and all other values are zero.

    Similarly the inverse transform maps the DoubleVector2d instances fhat_real and fhat_imag
    that have zero value except for a value of 1 in the (kx + (nx/2), ky + (ny/2))
    entry to two GridFunction1d instances whose values consist of the real and imaginary components
    of

    d  = sqrt(1/LX)*exp(2*pi*I*kx*(x-xMin)/LX)
        *sqrt(1/LY)*exp(2*pi*I*ky*(y-yMin)/LY)

    evaluated the points

    x = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny

    After computing the transform of the data, one often wants to
    work with the discrete Fourier coefficients. Typically this is done by
    looping over the wave numbers, and then obtaining the coefficient by
    computing the correct offset. For example, one uses code constructs
    as indicated by the fragment below.

    //
    // Loop over wave numbers
    //

    long kxIndex; long kyIndex;

    for(kx = -(nx/2); kx <= (nx-1)/2; kx++)
    {
    for(ky = -(ny/2); ky <= (ny-1)/2; ky++)
    {
    kxIndex = kx + (nx/2);             // (kxIndex, kyIndex) = index of the transform
    kyIndex = ky + (ny/2);             //  coefficient for the (kx,ky)'th  mode.

    rk = realCoeff(kxIndex,kyIndex);   // real coefficient of the (kx,ky)'th mode
    ik = imagCoeff(kxIndex,kyIndex);   // imag coefficient of the (kx,ky)'th mode

             *
             *
             *
    }}

*/
/*
#############################################################################
#
# Copyright 2014-2015 Chris Anderson
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


namespace UCLAQ
{


class fftw3_2d
{
public:

fftw3_2d()
{
	forwardplan = 0;
	inverseplan = 0;

    in  = 0;
    out = 0;

    nx  = 0;
    ny  = 0;

    LX = 1.0;
    LY = 1.0;

#ifdef _FFTW_OPENMP
     fftw_init_threads();
#endif
}

fftw3_2d(long nx, long ny, double LX = 1.0, double LY = 1.0)
{
    this->nx    = nx;
    this->ny    = ny;
    this->LX    = LX;
    this->LY    = LY;

    in          = 0;
    out         = 0;
	forwardplan = 0;
	inverseplan = 0;

#ifdef _FFTW_OPENMP
     fftw_init_threads();
#endif

	initialize(nx,ny);
}

virtual ~fftw3_2d()
{
    if(forwardplan != 0) fftw_destroy_plan(forwardplan);
    if(inverseplan != 0) fftw_destroy_plan(inverseplan); 
    
    if(in  != 0) fftw_free(in); 
    if(out != 0) fftw_free(out);

#ifdef _FFTW_OPENMP
    fftw_cleanup_threads();
#endif
}

#ifdef _FFTW_OPENMP
void initialize(long nx, long ny, int nthreads, double LX = 1.0, double LY = 1.0)
{
    fftw_plan_with_nthreads(nthreads);
    initialize(nx,ny,LX,LY);
}
#endif

void initialize()
{
    if(forwardplan != 0) fftw_destroy_plan(forwardplan);
    if(inverseplan != 0) fftw_destroy_plan(inverseplan); 
    
    if(in  != 0) fftw_free(in); 
    if(out != 0) fftw_free(out);

 	forwardplan = 0;
	inverseplan = 0;

    in  = 0;
    out = 0;

    nx  = 0;
    ny  = 0;

    LX = 1.0;
    LY = 1.0;
}

void initialize(long nx, long ny, double LX = 1.0, double LY = 1.0)
{
    if((this->nx != nx)||(this->ny != ny))
    {
    this-> nx = nx;
    this-> ny = ny;
    if(forwardplan != 0) fftw_destroy_plan(forwardplan);
    if(inverseplan != 0) fftw_destroy_plan(inverseplan); 

    if(in  != 0) fftw_free(in); 
    if(out != 0) fftw_free(out);
    
    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*ny);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*ny);
    

    forwardplan = fftw_plan_dft_2d(nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    inverseplan = fftw_plan_dft_2d(nx, ny, in, out, FFTW_BACKWARD,FFTW_ESTIMATE);
    
    //
    // Note: In the FFTW_BACKWARD plan the "out" argument is transformed 
    //       into the "in" argument when the plan is executed. 
    //
    }

    this->LX = LX;
    this->LY = LY;
}

// fftw2d_forward argument sizes:
//
// DoubleVector2d: Initialized of size nx by ny


void fftw2d_forward(UCLAQ::DoubleVector2d& inReal, UCLAQ::DoubleVector2d& inImag,
UCLAQ::DoubleVector2d& outReal, UCLAQ::DoubleVector2d& outImag)
{
	long i,j,k;

	if((nx != inReal.getIndex1Size())
    || (ny != inReal.getIndex2Size()))
    { 
    initialize(inReal.getIndex1Size(),
               inReal.getIndex2Size());
    }

	//reorder input
	
	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
		k = j + ny*i;

		in[k][0] = inReal(i,j);
		in[k][1] = inImag(i,j);
	}}

    fftw_execute(forwardplan);
	
    double scalingfactor =sqrt(LX*LY)/(nx*ny);

	for(k=0; k < nx*ny; k++)
	{
	
		i = (k/ny + (nx/2) )%nx ;//NB: uses integer division of odd numbers
		j = (k%ny + (ny/2) )%ny ;//NB: uses integer division of odd numbers

		outReal(i,j) = out[k][0]*scalingfactor;
		outImag(i,j) = out[k][1]*scalingfactor;
	}
}


// fftw2d_inverse argument sizes:
//
// DoubleVector2d: Initialized of size nx by ny
//

void fftw2d_inverse(UCLAQ::DoubleVector2d& inReal, UCLAQ::DoubleVector2d& inImag,
UCLAQ::DoubleVector2d& outReal, UCLAQ::DoubleVector2d& outImag)
{
	long i,j,k;
	long newindex1, newindex2;

	if((nx != inReal.getIndex1Size())
    || (ny != inReal.getIndex2Size()))
    { 
    initialize(inReal.getIndex1Size(),
               inReal.getIndex2Size());
    }

	//reorder input
	
	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
		newindex1 = (i - (nx/2) + nx )%nx ;//NB: uses integer division of odd numbers
		newindex2 = (j - (ny/2) + ny )%ny ;//NB: uses integer division of odd numbers

		k = newindex2 + ny*newindex1;

		in[k][0] = inReal(i, j);
		in[k][1] = inImag(i, j);
	}}

	//transform

	fftw_execute(inverseplan);

	double scalingfactor =1.0/sqrt(LX*LY);

	for(k=0; k < nx*ny; k++)
	{
		i = k/ny;
		j = k%ny;
		outReal(i,j) = out[k][0]*scalingfactor;
		outImag(i,j) = out[k][1]*scalingfactor;
	}
}


// fftw2d_forward argument sizes:
//
// GridFunction2d: Initialized with nx by ny panels
// DoubleVector2d: Initialized of size nx by ny


void fftw2d_forward(UCLAQ::GridFunction2d& inReal, UCLAQ::GridFunction2d& inImag,
UCLAQ::DoubleVector2d& outReal, UCLAQ::DoubleVector2d& outImag)
{
	long i,j,k;

	this->LX = inReal.getXmax() - inReal.getXmin();
 	this->LY = inReal.getYmax() - inReal.getYmin();

	if((nx != inReal.getXpanelCount())
    || (ny != inReal.getYpanelCount()))
    {
    initialize(inReal.getXpanelCount(),
               inReal.getYpanelCount());
    }

	// Capture and re-order input. This extraction
	// ignores the periodically defined values at
	// the outer edges of the computational domain.

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
		k = j + ny*i;//for changing from column to row ordering

		in[k][0] = inReal(i,j);
		in[k][1] = inImag(i,j);
	}}

    fftw_execute(forwardplan);

    double scalingfactor = sqrt(LX*LY)/(nx*ny);

	for(k=0; k < nx*ny; k++)
	{
		i = (k/ny + (nx/2) )%nx ;//NB: uses integer division of odd numbers
		j = (k%ny + (ny/2) )%ny ;//NB: uses integer division of odd numbers

		outReal( i,j) = out[k][0]*scalingfactor;
		outImag( i,j) = out[k][1]*scalingfactor;
	}
}


// fftw2d_inverse argument sizes:
//
// DoubleVector2d: Initialized of size nx by ny
// GridFunction2d: Initialized with nx by ny panels
//

void fftw2d_inverse(UCLAQ::DoubleVector2d& inReal, UCLAQ::DoubleVector2d& inImag,
UCLAQ::GridFunction2d& outReal, UCLAQ::GridFunction2d& outImag)
{
	long i,j,k;
	long newindex1, newindex2;

    this->LX = outReal.getXmax() - outReal.getXmin();
 	this->LY = outReal.getYmax() - outReal.getYmin();

	if((nx != inReal.getIndex1Size())
    || (ny != inReal.getIndex2Size()))
    {
    initialize(inReal.getIndex1Size(),
               inReal.getIndex2Size());
    }

	//reorder input

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
		newindex1 = (i - (nx/2) + nx )%nx ;//NB: uses integer division of odd numbers
		newindex2 = (j - (ny/2) + ny )%ny ;//NB: uses integer division of odd numbers

		k = newindex2 + ny*newindex1;

		in[k][0] = inReal(i, j);
		in[k][1] = inImag(i, j);
	}}

	//transform

	fftw_execute(inverseplan);

	double scalingfactor = 1.0/sqrt(LX*LY);

	for(k=0; k < nx*ny; k++)
	{
		i = k/ny;
		j = k%ny;
		outReal(i,j) = out[k][0]*scalingfactor;
		outImag(i,j) = out[k][1]*scalingfactor;
	}

	// Enforce periodicity of output GridFunction2d

	outReal.enforcePeriodicity();
	outImag.enforcePeriodicity();
}


private:

    long nx;   long ny;
    double LX; double LY;

	fftw_plan forwardplan;
	fftw_plan inverseplan;

    fftw_complex*  in;
    fftw_complex* out;


};
}
#endif
