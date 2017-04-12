
//
// fftw3_1dTest.cpp : Routine that tests the 1 dimensional FFTW3 routines
//
// Nov. 10, 2015
//

#include <cmath>
#include <cstdio>
using namespace std;

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"


// The following two routines compute the real and imaginary
// part of the k1th discrete Fourier basis function
//
//  exp(2*pi*I*k1*x)
//
//  for x in [0,1]
//

double realPartK1(long k1,double x)
{
    double pi = 3.1415926535897932384;
    return   cos(double(k1)*2.0*pi*x);
}

double imagPartK1(long k1, double x)
{
    double pi = 3.1415926535897932384;
    return   sin(double(k1)*2.0*pi*x);
}



int main(int argc, char* argv[])
{

    long nx   = 10;              // Number of panels in x-direction
    double dx = 1.0/double(nx);  // Mesh size

	long    i;
    double  x;

    long k1;                     // physical space wave number
    long k1Index;                // transform space vector index

	UCLAQ::DoubleVector1d	inReal(nx),    inImag(nx);
	UCLAQ::DoubleVector1d	outReal(nx),  outImag(nx);
    UCLAQ::DoubleVector1d	invReal(nx),  invImag(nx);

    //
    // Declare and initialize an instance of the fftw interface class
    //

	UCLAQ::fftw3_1d DFT;
	DFT.initialize(nx);

    //
    // Loop over wave numbers corresponding to a Fourier basis of
    // physical space.
    // 

    printf("#### Forward Transform Results #### \n");
    double errK1;

    for(k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    {

	for(i=0; i < nx; i++)
	{
        x = double(i)*dx;
		inReal(i) = realPartK1(k1,x);
        inImag(i) = imagPartK1(k1,x);
	}

    //
    // Compute transform
    //

	DFT.fftw1d_forward(inReal, inImag,outReal,outImag);

    k1Index = k1 + (nx/2); // k1Index = array index of the transform value
                           //           of the k1th mode.

    //
    // The transform of k1th basis function should be a 
    // vector whose imaginary component is zero and 
    // whose real component has the value 1 in the
    // k1Index = k1 + (nx/2) location.
    //
    
    // 
    // Subtract what the solution should be and 
    // then take norms of the result for the error 
    //

    outReal(k1Index)= outReal(k1Index) - 1.0;

    errK1  = sqrt(outReal.dot(outReal));
    errK1 += sqrt(outImag.dot(outImag));

    printf("Error in forward transform of %+ld mode %+10.5e \n",k1, errK1);
    }

    printf("\n#### Inverse Transform Results #### \n");
//
//  Check Inverse transform
//
    for(k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    {

    //
    // set up transform of k1th Fourier mode 
    // 

    outReal.setToValue(0.0);
    outImag.setToValue(0.0);

    k1Index = k1 + (nx/2);
    outReal(k1Index) = 1.0;

    //
    // Compute inverse transform
    //

    DFT.fftw1d_inverse(outReal,outImag, invReal, invImag);

    //
    // Check to see if it maps back to the k1th Fourier mode
    //

	for(i=0; i < nx; i++)
	{
        x = double(i)*dx;
		inReal(i) = realPartK1(k1,x);
        inImag(i) = imagPartK1(k1,x);
	}

    // subtract off what it should be and then take norms 

    invReal -= inReal;
    invImag -= inImag;
    
    errK1  = sqrt(invReal.dot(invReal));
    errK1 += sqrt(invImag .dot(invImag ));
    printf("Error in inverse transform of %+ld mode %+10.5e \n",k1, errK1);
    }
 
	return 0;
}
