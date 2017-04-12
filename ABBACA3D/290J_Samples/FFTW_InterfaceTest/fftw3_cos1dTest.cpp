
//
// fftw3_cos1dTest.cpp : Routine that tests the 1 dimensional FFTW3 cos transform
//                       routines
//
// Nov. 10, 2015
//

#include <cmath>
#include <cstdio>
using namespace std;


#include "FFTW3_InterfaceNd/UCLAQ_fftw3_cos1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"


int main(int argc, char* argv[])
{
    double pi = 3.1415926535897932384;
    
    long nx   = 10;              // Number of panels in x-direction
    double dx = 1.0/double(nx);  // Mesh size

	long    i;
    double  x;

    long k;                     // physical space wave number

	UCLAQ::DoubleVector1d	 inF(nx+1);
	UCLAQ::DoubleVector1d	outF(nx+1);
	
    UCLAQ::DoubleVector1d	 invF(nx+1);

    //
    // Declare and initialize an instance of the fftw interface class
    //

	UCLAQ::fftw3_cos1d DFT;
	DFT.initialize(nx);

    //
    // Loop over wave numbers corresponding to a Fourier basis of
    // physical space.
    // 

    printf("#### Forward Transform Results #### \n");
    double errK;

    double normFactor;

    for(k = 0; k <= nx; k++)
    {

    //
    // Create kth orthonormal cos function
    //
    if((k == 0)||(k == nx)){normFactor = 1.0;}
    else                   {normFactor = sqrt(2.0);}

	for(i=0; i <= nx; i++)
	{
        x = double(i)*dx;
		inF(i) = normFactor*cos(double(k)*pi*x);
	}

    //
    // Compute in/out transform 
    //

	DFT.fftw1d_cos_forward(inF,outF);
	
    //
    // The transform of kth basis function should be a 
    // vector that has the the value 1 in the kth location
    // where k = 0 to nx
    //
    
    outF(k) = outF(k) - 1.0;
    errK   = sqrt(outF.dot(outF));
    
    //
    // Compute and check "in place" transform 
    //

	DFT.fftw1d_cos_forward(inF);
	
    inF(k) = inF(k) - 1.0;
    
    errK   += sqrt(inF.dot(inF));
    
    printf("Error in forward transform of %+ld mode %+10.5e \n",k, errK);
    }

    printf("\n\n");

    printf("\n#### Inverse Transform Results #### \n");
//
//  Check Inverse transform
//
    for(k = 0; k <= nx; k++)
    {

    //
    // set up transform of k1th Fourier mode 
    // 

    outF.setToValue(0.0);
    outF(k) = 1.0;

    //
    // Compute inverse transform
    //

    DFT.fftw1d_cos_inverse(outF,invF);

    //
    // Check to see if it maps back to the kth orthonormal sin
    // mode 
    //
    if((k == 0)||(k == nx)){normFactor = 1.0;}
    else                   {normFactor = sqrt(2.0);}

	for(i=0; i <= nx; i++)
	{
        x = double(i)*dx;
		inF(i) = normFactor*cos(double(k)*pi*x);
	}


    if(k == nx+6)
    {
	for(i=0; i <= nx; i++)
	{
	printf("%10.5e  %10.5e \n",invF(i),inF(i));
	}

    return 0;
    }

    // subtract off what it should to evaluate the error

    invF -= inF;
    
    errK  = sqrt(invF.dot(invF));
    

    //
    // Check "in place" inverse transform
    //
    
    DFT.fftw1d_cos_inverse(outF);
    
    outF -= inF;
    
    errK  += sqrt(outF.dot(outF));
    
    printf("Error in inverse transform of %+ld mode %+10.5e \n",k, errK);
    }

	return 0;
}
