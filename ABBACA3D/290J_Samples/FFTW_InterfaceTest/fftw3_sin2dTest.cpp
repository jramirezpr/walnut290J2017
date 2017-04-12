//
// fftw3_sin2dTest.cpp  : Routine that tests the 2D 
//                        FFTW sin transform routines
//
// Nov. 10, 2015
//


#include <cmath>
#include <cstdio>
using namespace std;

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin2d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"

// (k1,k2)th orthonormal sin function (w.r.t. mesh width scaled inner product).

double sinXsinY(long k1, long k2,double x, double y)
{
    double pi = 3.1415926535897932384;
    return  sqrt(2.0)*sin(double(k1)*pi*x)*sqrt(2.0)*sin(double(k2)*pi*y);
}

int main(int argc, char* argv[])
{
    long nx = 8;                  // number of panels in x-direction
    long ny = 7;                  // number of panels in y-direction

    double dx = 1.0/double(nx);   // reciprocal of panel counts 
    double dy = 1.0/double(ny);   // reciprocal of panel counts

	long i,j;

	UCLAQ::DoubleVector2d	 inF(nx-1,ny-1);
	UCLAQ::DoubleVector2d   outF(nx-1,ny-1);
    UCLAQ::DoubleVector2d	invF(nx-1,ny-1);

    //
    // Declare and initialize an instance of the 2D fftw sin transform 
    // interface class. 
    // 
    //
	UCLAQ::fftw3_sin2d DFT;
	DFT.initialize(nx,ny);
	

    double x; double y;
    long k1;   long k2;       // wave number indices 
    //
    // Loop over wave numbers corresponding to a sin basis of
    // physical space.
    // 
    printf("#### Forward Transform Results #### \n");

    double errK1K2;
    double errK1K2max = 0.0;

    for(k1 = 1; k1 <= nx-1; k1++)
    {
    for(k2 = 1; k2 <= ny-1; k2++)
    {
    //
    // construct (k1,k2)th sin mode
    //

	for(i=1; i <= nx-1; i++)
	{
	for(j=1; j <= ny-1; j++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
		inF(i-1,j-1) = sinXsinY(k1,k2,x,y);
	}}

    //
    // Compute transform
    //

	DFT.fftw2d_sin_forward(inF, outF);
	
    //
    // The transform of (k1th,k2th) basis function should be a 
    // vector whose imaginary component is zero and 
    // whose real component has the value 1 in the
    // (k1, k2) location.
    //
    
    // To compute the error, subtract off what the solution should
    // and then take norms of the result 
    //
    
    outF(k1-1,k2-1) = outF(k1-1,k2-1) - 1.0;
    errK1K2     = sqrt(outF.dot(outF));
    
    errK1K2max = ( errK1K2max < abs(errK1K2) ) ? errK1K2 : errK1K2max;

    //
    // Check in place transform
    //

	DFT.fftw2d_sin_forward(inF);
	
	inF(k1-1,k2-1) = inF(k1-1,k2-1) - 1.0;
    errK1K2   += sqrt(inF.dot(inF));
	
    
    //printf("Error in forward transform of (%+ld , %+ld) mode %+10.5e \n",k1,k2, errK1K2);
    }

    }

    printf("Maximal error in forward transform of modes in range (%+d , %+d) to (%+ld , %+ld) :   %+10.5e \n"
            ,1,1,nx-1,ny-1, errK1K2max);

    errK1K2max = 0.0;
    printf("\n#### Inverse Transform Results #### \n");
//
//  Check Inverse transform
//
    for(k1 = 1; k1 <= nx-1; k1++)
    {
    for(k2 = 1; k2 <= ny-1; k2++)
    {
    //
    // inF is the exact inverse transform 
    //
	for(i=1; i <= nx-1; i++)
	{
	for(j=1; j <= ny-1; j++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
		inF(i-1,j-1) = sinXsinY(k1,k2,x,y);
	}}
	
    //
    // Set up transform value of (k1,k2)th Fourier mode 
    // 
    outF.setToValue(0.0);
    outF(k1-1,k2-1) = 1.0;

    //
    // Compute the inverse transform and then the error 
    //

    DFT.fftw2d_sin_inverse(outF,invF);
    
    invF    -= inF;
    errK1K2  = sqrt( invF.dot( invF));
    
    errK1K2max = ( errK1K2max < abs(errK1K2) ) ? errK1K2 : errK1K2max;

    // Check in place transform
    
    DFT.fftw2d_sin_inverse(outF);
    
    outF     -= inF;
    errK1K2  += sqrt(outF.dot(outF));
    
    //printf("Error in inverse transform of (%+ld , %+ld) mode %+10.5e \n",k1,k2, errK1K2);

    }}

    printf("Maximal error in inverse transform of modes in range (%+d , %+d) to (%+ld , %+ld) :   %+10.5e \n"
            ,1,1,nx-1,ny-1, errK1K2max);

	return 0;
}
