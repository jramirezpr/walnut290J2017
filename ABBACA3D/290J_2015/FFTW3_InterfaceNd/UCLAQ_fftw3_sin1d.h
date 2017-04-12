#include "fftw3.h"             // fftw3 headers

#include "../DoubleVectorNd/UCLAQ_DoubleVector1d.h"    // Array class header
#include "../GridFunctionNd/UCLAQ_GridFunction1d.h"    // Grid function class header

#ifndef _UCLAQ_fftw3_sin1d_
#define _UCLAQ_fftw3_sin1d_
//
//####################################################################
//                fftw3_sin1d.h : Version Nov 12, 2015
//####################################################################
/**<p>

    This class provides forward and inverse sin transforms of one dimensional
    data stored in UCLAQ::GridFunction1d instances or data stored in
    UCLAQ::DoubleVector1d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical sin transforms of
    functions defined over intervals of size LX.

    The transforms are implemented using the r2r_1d(...) routine from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/). The
    FFTW3 headers and library files must be available to create any executable
    that uses this class.

    This transform interface is designed be used with the data values of
    discrete functions with homogeneous boundary values associated with
    a uniform grid discretization consisting of nx grid panels spanning a
    domain of size LX.

    When DoubleVector1d's are used for input and/or output their size
    is (nx-1); the number of interior data values associated with
    such discretizations.

    For consistency with other FFTW_InterfaceNd classes, in the
    construction or initialization of class instances, the number
    of panels nx of a discretization is specified.

    For computational efficiency nx should be chosen to be a product of small
    primes (if possible).

    What's computed assuming the input data is specified as a
    GridFunction1d instance initialized with the parameters
    (nx, xMin, xMax):

    If kx is a wave number 1 <= kx <= nx-1

    then GridFunction1d instances with values s(i) defined by

    s(i)  = sqrt(2/LX)*sin(kx*pi*(x_i-xMin)/LX)

    with x_i = xMin + i/nx and i = 0,1,...,nx

    are mapped by the forward transform to a DoubleVection1d instance
    of size (nx-1), f_hat,  with a single non-zero entry

    f_hat(kx-1) = 1.0

    All other values of f_hat are identically zero.

    Similarly the inverse transform maps the double vector f_hat
    that is zero except for f_hat(kx-1) = 1.0

    to the GridFunction1d instance whose ith value is given by

    s(i)  = sqrt(2.0/LX)*sin(kx*pi*(x_i-xMin)/LX)

    with x_i = xMin + i/nx and i = 0,1,...,nx


*/
/*
#############################################################################
#
# Copyright 2013-2015 Chris Anderson
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
class fftw3_sin1d
{
public:

fftw3_sin1d()
{
	plan     = 0;
    in       = 0;
    out      = 0;
    nx       = 0;
    LX       = 1.0;
    nSamples = 0;
}

fftw3_sin1d(long nx, double LX = 1.0)
{
    this->nx    = nx;
    this->LX    = LX;
    this->nSamples = nx-1;

    in          = 0;
    out         = 0;
	plan        = 0;

	initialize(nx);
}

~fftw3_sin1d()
{
    if(plan != 0) fftw_destroy_plan(plan);
    if(in  != 0)  fftw_free(in);
    if(out != 0)  fftw_free(out);
    plan = 0;
    in   = 0;
    out  = 0;
}

void initialize()
{
    if(plan != 0) fftw_destroy_plan(plan);
    if(in   != 0)  fftw_free(in);
    if(out  != 0) fftw_free(out);

    plan     = 0;

    in       = 0;
    out      = 0;
    nx       = 0;
    LX       = 1.0;
    nSamples = 0;
}

void initialize(long nx, double LX = 1.0)
{
    if((this->nx != nx))
    {
    this->nx       = nx;
    this->nSamples = nx-1;
    if(plan != 0) fftw_destroy_plan(plan);

    if(in  != 0) fftw_free(in); 
    if(out != 0) fftw_free(out);

	in  = (double*) fftw_malloc(sizeof(double) * nSamples);
    out = (double*) fftw_malloc(sizeof(double) * nSamples);
    plan = fftw_plan_r2r_1d(nSamples, in, out, FFTW_RODFT00, FFTW_ESTIMATE);

    // Efficiency: From the documentation for RODFT00 transform, nSampleX+1
    // should be a product of small primes ==> nx should be product of small primes.
    }

    this->LX = LX;
}


// fftw1d_sin_forward argument sizes:
//
// This operator ignores the perimeter values of the input vector.
//
// DoubleVector1d: outF size nx-1

void fftw1d_sin_forward(UCLAQ::GridFunction1d&  inF,  UCLAQ::DoubleVector1d& outF)
{
    long k;

    // Capture domain size

    this->LX = inF.getXmax() - inF.getXmin();

	if(nx != inF.getXpanelCount())
    {
    initialize(inF.getXpanelCount());
    }

	// copy input ignoring perimeter values

	for(k=0; k < nSamples; k++)
	{
		in[k] = inF(k+1);
	}

    fftw_execute(plan);

    double scalingfactor =  sqrt(LX)/(sqrt(2.0)*((double)(nx)));

	for(k=0; k < nSamples; k++)
	{
		outF(k) = out[k]*scalingfactor;
	}

}

// fftw1d_sin_inverse argument sizes:
//
// DoubleVector1d: inF size nx-1
//

void fftw1d_sin_inverse(UCLAQ::DoubleVector1d&  inF,  UCLAQ::GridFunction1d& outF)
{

	// Capture domain size

    this->LX = outF.getXmax() - outF.getXmin();

	if(nx != inF.getSize()+1)
    {
    initialize(inF.getSize()+1);
    }

	//copy input

	for(long k=0; k < nSamples; k++)
	{
		in[k] = inF(k);
	}

    fftw_execute(plan);

	double scalingfactor = 1.0/sqrt(2.0*LX);

	//reorder AND scale output

	for(long k=0; k < nSamples; k++)
	{
		outF(k+1) = out[k]*scalingfactor;
	}

    // zero perimeter values

	outF.setBoundaryValues(0.0);
}

// fftw1d_sin_forward argument sizes:
//
// DoubleVector1d's are of size nx-1
//
void fftw1d_sin_forward(UCLAQ::DoubleVector1d&  inF,  UCLAQ::DoubleVector1d& outF)
{
    long k;
    double scalingfactor;

	if(nx != inF.getSize()+1)
    {
    initialize(inF.getSize()+1);
    }

	//copy input

	for(k=0; k < nSamples; k++)
	{
		in[k] = inF(k);
	}

    fftw_execute(plan);

    scalingfactor =  sqrt(LX)/(sqrt(2.0)*((double)(nx)));

	for(k=0; k < nSamples; k++)
	{
		outF(k) = out[k]*scalingfactor;
	}
}

// fftw1d_sin_inverse argument sizes:
//
// DoubleVector1d's are of size nx-1
//

void fftw1d_sin_inverse(UCLAQ::DoubleVector1d&  inF,  UCLAQ::DoubleVector1d& outF)
{
	long k;
	double scalingfactor;

	if(nx != inF.getSize()+1)
    { 
    initialize(inF.getSize()+1);
    }

	//copy input
	
	for(k=0; k < nSamples; k++)
	{
		in[k] = inF(k);
	}

    fftw_execute(plan);
	
	scalingfactor = 1.0/sqrt(2.0*LX);

	//reorder AND scale output

	for(k=0; k < nSamples; k++)
	{
		outF(k) = out[k]*scalingfactor;
	}
}

//
// "In place transform"
//
// DoubleVector1d is of size nx-1

void fftw1d_sin_forward(UCLAQ::DoubleVector1d&  F)
{
    long k;

	if(nx != F.getSize()+1)
    { 
    initialize(F.getSize()+1);
    }

	//copy input
	
	for(k=0; k < nSamples; k++)
	{
		in[k] = F(k);
	}

    fftw_execute(plan);

    double scalingfactor = sqrt(LX)/(sqrt(2.0)*((double)(nx)));
	
	for(k=0; k < nSamples; k++)
	{
		F(k) = out[k]*scalingfactor;
	}
}

// DoubleVector1d is of size nx-1

void fftw1d_sin_inverse(UCLAQ::DoubleVector1d&  F)
{
	long k; 

	if(nx != F.getSize()+1)
    { 
    initialize(F.getSize()+1);
    }

	//copy input
	
	for(k=0; k < nSamples; k++)
	{
		in[k] = F(k);
	}

    fftw_execute(plan);
	
	double scalingfactor = 1.0/sqrt(2.0*LX);
	
	//reorder AND scale output
	
	for(k=0; k < nSamples; k++)
	{
		F(k) = out[k]*scalingfactor;
	}
}

private:

    long nx; double LX;

    long nSamples;

	fftw_plan plan;

    double*  in;
    double* out;

};
}
#endif
