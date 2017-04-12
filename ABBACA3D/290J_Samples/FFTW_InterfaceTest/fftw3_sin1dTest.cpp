//
// fftw3_sin1dTest.cpp  : Routine that tests the 1D 
//                        FFTW sin transform routines
//
// Nov. 6, 2015
//

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"

#include <cmath>
#include <cstdio>

using namespace std;

int main(int argc, char* argv[])
{
    double pi = 3.1415926535897932384;
    
    long nx   = 10;               // Number of panels in x-direction
    double dx = 1.0/double(nx);  // Mesh size

	long    i;
    double  x;

    long k;                     // physical space wave number

    //
    // Vectors of size nx-1 corresponding to nx-1 interior data values on
    // a grid discretized with nx panels
    //
	UCLAQ::DoubleVector1d	   inF(nx-1);
	UCLAQ::DoubleVector1d	  outF(nx-1);
    UCLAQ::DoubleVector1d	  invF(nx-1);

    // Declare and initialize an instance of the fftw interface class

	UCLAQ::fftw3_sin1d DFT;
	DFT.initialize(nx);

    // Loop over wave numbers corresponding to a Fourier basis of
    // physical space.


    printf("#### Forward Transform Results #### \n");
    double errK;

    long kIndex;

    for(k = 1; k <= nx-1; k++)
    {
    kIndex = k-1;
    //
    // Create kth orthonormal sin function of interior values
    //
	for(i=1; i < nx; i++)
	{
        x = double(i)*dx;
		inF(i-1) = sqrt(2.0)*sin(double(k)*pi*x);
	}

    //
    // Compute in/out transform 
    //

	DFT.fftw1d_sin_forward(inF, outF);
	
    //
    // The transform of kth basis function should be a 
    // vector that has the the value 1 in the kth location
    // where k = 1 to nx-1
    //
    
    outF(kIndex) = outF(kIndex) - 1.0;
    errK         = sqrt(outF.dot(outF));
    
    //
    // Compute and check "in place" transform 
    //

	DFT.fftw1d_sin_forward(inF);
	
    inF(kIndex) = inF(kIndex) - 1.0;
    
    errK   += sqrt(inF.dot(inF));
    
    printf("Error in forward transform of %+ld mode %+10.5e \n",k, errK);
    }

    printf("\n#### Inverse Transform Results #### \n");
//
//  Check Inverse transform
//
    for(k = 1; k <= nx-1; k++)
    {
    kIndex = k-1;
    //
    // set up transform of k1th Fourier mode 
    // 

    outF.setToValue(0.0);
    outF(kIndex) = 1.0;

    //
    // Compute inverse transform
    //

    DFT.fftw1d_sin_inverse(outF,invF);

    //
    // Check to see if it maps back to the kth orthonormal sin
    // mode 
    //

	for(i=1; i <= nx-1; i++)
	{
        x = double(i)*dx;
		inF(i-1) = sqrt(2.0)*sin(double(k)*pi*x);
	}

    // subtract off what it should to evaluate the error

    invF -= inF;
    
    errK  = sqrt(invF.dot(invF));
    
    //
    // Check "in place" inverse transform
    //
    
    DFT.fftw1d_sin_inverse(outF);
    
    outF -= inF;
    
    errK  += sqrt(outF.dot(outF));
    
    printf("Error in inverse transform of %+ld mode %+10.5e \n",k, errK);
    }

	return 0;
}
