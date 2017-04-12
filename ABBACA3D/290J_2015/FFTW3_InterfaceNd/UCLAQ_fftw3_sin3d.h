#include "fftw3.h"

#include "../DoubleVectorNd/UCLAQ_DoubleVector3d.h"    // Array class header
#include "../GridFunctionNd/UCLAQ_GridFunction3d.h"    // Grid function class header

#ifdef _FFTW_OPENMP
#include <omp.h>
#endif

#ifndef _UCLAQ_fftw3_sin3d_
#define _UCLAQ_fftw3_sin3d_

//
//####################################################################
//          UCLAQ_fftw3_sin3d.h : Version Nov. 12, 2015
//####################################################################
//
/**<p>
    This class provides forward and inverse sin transforms of three dimensional
    data stored in UCLAQ::GridFunction3d instances or data stored in
    UCLAQ::DoubleVector3d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical sin transforms of
    functions defined over intervals of size LX x LY x LZ.

    The transforms are implemented using the r2r_3d(...) routines from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/). The
    FFTW3 headers and library files must be available to create any executable
    that uses this class.

    This transform interface is designed be used with the data values of
    discrete functions with homogeneous boundary values associated with
    a uniform grid discretization consisting of nx grid panels by ny grid panels
    by nx grid panels spanning a domain of size LX x LY x LZ.

    When DoubleVector3d's are used for input and/or output their size
    is (nx-1) x (ny-1) x (nz-1); the number of interior data values associated with
    such discretizations.

    For consistency with other FFTW_InterfaceNd classes, in the
    construction or initialization of class instances, the number
    of panels nx, ny and nz of a discretization is specified.

    For computational efficiency nx, ny and nz should be chosen to be a
    product of small primes less than or equal to 13  (if possible).
    
    ===============================================================

    What's computed assuming the input data is specified as a
    GridFunction3d instance initialized with the parameters
    (nx, xMin, xMax, ny, yMin, yMax, nz, zMin, zMax):

    If (kx,ky,kz) is a wave number vector with kx, ky, kz in the ranges

    1 <= kx <=  (nx-1) , 1 <= ky <= (ny-1) and 1 <= kz <= (nz-1)
    
    then a GridFunction3d instance whose values consist of

    d  = sqrt(2/LX)*sin(kx*pi*(x_i-xMin)/LX)
        *sqrt(2/LY)*sin(ky*pi*(y_j-yMin)/LY)
        *sqrt(2/LZ)*sin(ky*pi*(z_k-zMin)/LZ)

    evaluated the points

    x_i = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y_j = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny
    z_k = zMin + k*hz, k = 0,1,...,nz-1, hz = LZ/nz

    gets mapped by the forward transform to a DoubleVector3d of size
    (nx-1)x*(ny-1)x(nz-1), f_hat, where

    f_hat(kx-1, ky-1, kz-1) = 1

    and all other values are zero.

    Similarly the inverse transform maps the DoubleVector3d of size
    (nx-1)x*(ny-1)x(nz-1), f_hat, that is zero value except for a value of 1 in
    (kx-1,ky-1,kz-1) place to the GridFunction3d instance with values
    corresponding to

    d  = sqrt(2/LX)*sin(kx*pi*(x_i-xMin)/LX)
        *sqrt(2/LY)*sin(ky*pi*(y_j-yMin)/LY)
        *sqrt(2/LZ)*sin(ky*pi*(z_k-zMin)/LZ)

    evaluated the points

    x_i = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y_j = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny
    z_k = zMin + k*hz, k = 0,1,...,nz-1, hz = LZ/nz

    After computing the transform of the data, one often wants to
    work with the discrete Fourier coefficients. Typically this is done by
    looping over the wave numbers, and then obtaining the coefficient by
    computing the correct offset. For example, one uses code constructs
    as indicated by the fragment below.

    //
    // Loop over sin wave numbers
    //

    for(long kx = 1; kx <= nx-1; kx++)
    {
    for(long ky = 1; ky <= ny-1; ky++)
    {
    for(long kz = 1; kz <= nz-1; kz++)
    {
    	sCoeff  = f_hat(kx-1,ky-1,kz-1);  // (kx,ky,kz)'th  sin coefficient is
    	        *                         // the (kx-1,ky-1,kz-1) entry of the transform
    	        *
    	        *
    }}}
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
//####################################################################
// Chris Anderson (C) UCLA                             Nov. 12, 2015
//####################################################################
//
namespace UCLAQ
{
class fftw3_sin3d
{
public:

fftw3_sin3d()
{
	plan = 0;
    in   = 0;
    out  = 0;

    nx  = 0;
    ny  = 0;
    nz  = 0;

    LX = 1.0;
    LY = 1.0;
    LZ = 1.0;

    nSampleX = 0;
    nSampleY = 0;
    nSampleZ = 0;

    #ifdef _FFTW_OPENMP
     fftw_init_threads();
    #endif
}

fftw3_sin3d(long nx, long ny, long nz, double LX = 1.0, double LY = 1.0, double LZ = 1.0)
{
    this->nx    = nx;
    this->ny    = ny;
    this->nz    = nz;

    this->LX = LX;
    this->LY = LY;
    this->LZ = LZ;

    nSampleX   = nx-1;
    nSampleY   = ny-1;
    nSampleZ   = nz-1;

    in          = 0;
    out         = 0;
	plan        = 0;

	#ifdef _FFTW_OPENMP
     fftw_init_threads();
    #endif


	initialize(nx,ny,nz);
}

virtual ~fftw3_sin3d()
{
    if(plan != 0) fftw_destroy_plan(plan);

    if(in  != 0) fftw_free(in);
    if(out != 0) fftw_free(out);

    #ifdef _FFTW_OPENMP
    fftw_cleanup_threads();
    #endif
}

void initialize()
{
    if(plan != 0) fftw_destroy_plan(plan);
    if(in   != 0) fftw_free(in);
    if(out  != 0) fftw_free(out);

 	plan = 0;

    in  = 0;
    out = 0;

    nx  = 0;
    ny  = 0;
    nz  = 0;

    LX = 1.0;
    LY = 1.0;
    LZ = 1.0;

    nSampleX = 0;
    nSampleY = 0;
    nSampleZ = 0;
}

void initialize(long nx, long ny, long nz, double LX = 1.0, double LY = 1.0, double LZ = 1.0)
{
    if((this->nx != nx)||(this->ny != ny)||(this->nz != nz))
    {
    this->nx   = nx;
    this->ny   = ny;
    this->nz   = nz;
    nSampleX   = nx-1;
    nSampleY   = ny-1;
    nSampleZ   = nz-1;

    if(plan != 0) fftw_destroy_plan(plan);

    if(in  != 0) fftw_free(in);
    if(out != 0) fftw_free(out);

    in  = (double*) fftw_malloc(sizeof(double) * nSampleX*nSampleY*nSampleZ);
    out = (double*) fftw_malloc(sizeof(double) * nSampleX*nSampleY*nSampleZ);

    plan = fftw_plan_r2r_3d(nSampleX, nSampleY, nSampleZ, in, out,FFTW_RODFT00,
           FFTW_RODFT00, FFTW_RODFT00,FFTW_ESTIMATE);

    // Efficiency: From the documentation for RODFT00 transform, nSampleX+1, nSampleY+1, nSampleZ+1
    // should be a product of small primes ==> nx, ny and nz should be product of small primes.
    }

    this->LX = LX;
    this->LY = LY;
    this->LZ = LZ;
}

#ifdef _FFTW_OPENMP
void initialize(long nx, long ny, long nz, int nthreads, double LX = 1.0, double LY = 1.0, double LZ = 1.0)
{
    fftw_plan_with_nthreads(nthreads);
    initialize(nx,ny,nz);
}
#endif


// fftw3d_sin_forward argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//
// This operator ignores the perimeter values of the input GridFunction3d

void fftw3d_sin_forward(UCLAQ::GridFunction3d& inF, UCLAQ::DoubleVector3d& outF)
{
	this->LX = inF.getXmax() - inF.getXmin();
 	this->LY = inF.getYmax() - inF.getYmin();
 	this->LZ = inF.getZmax() - inF.getZmin();

	long i,j,k;

	if((nx != inF.getXpanelCount()) || (ny != inF.getYpanelCount())|| (nz != inF.getZpanelCount()))
    {
    initialize(inF.getXpanelCount(),inF.getYpanelCount(),inF.getZpanelCount());
    }

	for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		in[k + nSampleZ*(j+ i*nSampleY)] = inF(i+1,j+1,k+1);
	}}}

    fftw_execute(plan);

    // Capture transform values and scale appropriately

    double scalingfactor =  sqrt(LX*LY*LZ)/(2.0*sqrt(2.0)*((double)(nx)*(double)(ny)*(double)(nz)));

    for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		outF(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
	}}}
}

// fftw3d_sin_inverse argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_inverse(UCLAQ::DoubleVector3d& inF, UCLAQ::GridFunction3d& outF)
{
	long i,j,k;

	this->LX = outF.getXmax() - outF.getXmin();
 	this->LY = outF.getYmax() - outF.getYmin();
 	this->LZ = outF.getZmax() - outF.getZmin();

	if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1)|| (nz != inF.getIndex3Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1,inF.getIndex3Size()+1);
    }

	for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		in[k + nSampleZ*(j+ i*nSampleY)] = inF(i,j,k);
	}}}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(2.0*LX*LY*LZ));

    for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		outF(i+1,j+1,k+1) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
	}}}

    outF.setBoundaryValues(0.0);
}


    
// fftw3d_sin_forward argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_forward(UCLAQ::DoubleVector3d& F)
{
	long i,j,k;

	if((nx != F.getIndex1Size()+1) || (ny != F.getIndex2Size()+1)|| (nz != F.getIndex3Size()+1))
    {
    initialize(F.getIndex1Size()+1,F.getIndex2Size()-1,F.getIndex3Size()+1);
    }

	for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		in[k + nSampleZ*(j+ i*nSampleY)] = F(i,j,k);
	}}}

    fftw_execute(plan);

    // Capture transform values and scale appropriately

    double scalingfactor = sqrt(LX*LY*LZ)/(2.0*sqrt(2.0)*((double)(nx)*(double)(ny)*(double)(nz)));

    for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		F(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
	}}}
}

// fftw3d_sin_forward argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_forward(UCLAQ::DoubleVector3d& inF, UCLAQ::DoubleVector3d& outF)
{
	long i,j,k;

	if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1)|| (nz != inF.getIndex3Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1,inF.getIndex3Size()+1);
    }

	for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		in[k + nSampleZ*(j+ i*nSampleY)] = inF(i,j,k);
	}}}

    fftw_execute(plan);

    // Capture transform values and scale appropriately
    
    double scalingfactor =  sqrt(LX*LY*LZ)/(2.0*sqrt(2.0)*((double)(nx)*(double)(ny)*(double)(nz)));

    for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		outF(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
	}}}
}
	
// fftw3d_sin_inverse argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_inverse(UCLAQ::DoubleVector3d& inF, UCLAQ::DoubleVector3d& outF)
{
	long i,j,k;

	if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1)|| (nz != inF.getIndex3Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1,inF.getIndex3Size()+1);
    }

	for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		in[k + nSampleZ*(j+ i*nSampleY)] = inF(i,j,k);
	}}}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(2.0*LX*LY*LZ));

    for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		outF(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
	}}}
}

// fftw3d_sin_inverse argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//
	
void fftw3d_sin_inverse(UCLAQ::DoubleVector3d&  F)
{
	long i,j,k;

	if((nx != F.getIndex1Size()+1) || (ny != F.getIndex2Size()+1)|| (nz != F.getIndex3Size()+1))
    {
    initialize(F.getIndex1Size()+1,F.getIndex2Size()+1,F.getIndex3Size()+1);
    }

	for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		in[k + nSampleZ*(j+ i*nSampleY)] = F(i,j,k);
	}}}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(2.0*LX*LY*LZ));

    for(i=0; i < nSampleX; i++)
	{
	for(j=0; j < nSampleY; j++)
	{
	for(k=0; k < nSampleZ; k++)
	{
		F(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
	}}}
}


private:


    long nx;   long ny;     long nz;
    double LX; double LY; double LZ;
    
    long nSampleX;
    long nSampleY;
    long nSampleZ;

	fftw_plan plan;

    double*  in;
    double* out;
};
}
#endif

/*

int fftw_init_threads(void);

This function, which need only be called once, performs any one-time initialization required
to use threads on your system. It returns zero if there was some error (which should not happen under normal circumstances)
and a non-zero value otherwise.

Third, before creating a plan that you want to parallelize, you should call:
void fftw_plan_with_nthreads(int nthreads);

The nthreads argument indicates the number of threads you want
FFTW to use (or actually, the maximum number). All plans subsequently created with any planner
routine will use that many threads. You can call fftw_plan_with_nthreads, create some plans,
call fftw_plan_with_nthreads again with a different argument, and create some more
plans for a new number of threads. Plans already created before a call to fftw_plan_with_nthreads
are unaffected. If you pass an nthreads argument of 1 (the default), threads are disabled for subsequent plans.

With OpenMP, to configure FFTW to use all of the currently running OpenMP threads
(set by omp_set_num_threads(nthreads) or by the OMP_NUM_THREADS environment variable), you can do:

fftw_plan_with_nthreads(omp_get_max_threads()). (The ‘omp_’ OpenMP functions are declared via #include <omp.h>.)

Given a plan, you then execute it as usual with fftw_execute(plan), and the execution will
use the number of threads specified when the plan was created. When done, you destroy it as usual with fftw_destroy_plan. As described in Thread safety, plan execution is thread-safe, but plan creation and destruction are not: you should create/destroy plans only from a single thread, but can safely execute multiple plans in parallel.

There is one additional routine: if you want to get rid of all memory and other resources
allocated internally by FFTW, you can call:

void fftw_cleanup_threads(void);

 */
