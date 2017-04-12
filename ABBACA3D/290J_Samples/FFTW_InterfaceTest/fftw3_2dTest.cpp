//
// fftw3_2dTest.cpp  : Routine that tests the 2D FFTW routines
//
// Nov. 10, 2015
//

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_2d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"

#include <cmath>
#include <cstdio>
using namespace std;

#ifdef _FFTW_OPENMP
#include <omp.h>
#endif



// The following two routines compute the real and imaginary
// part of the (k1,k2) discrete Fourier basis function
//
//  exp(2*pi*I*k1*x)*exp(2*pi*I*k2*y) 
//
//  for (x,y) in [0,1]X[0,1]
//

double realPartK1K2(long k1, long k2,double x, double y)
{
    double pi =  3.1415926535897932384;
    return   cos(double(k1)*2.0*pi*x)*cos(double(k2)*2.0*pi*y)
           - sin(double(k1)*2.0*pi*x)*sin(double(k2)*2.0*pi*y);

}

double imagPartK1K2(long k1, long k2,double x, double y)
{
    double pi =  3.1415926535897932384;
    return   sin(double(k1)*2.0*pi*x)*cos(double(k2)*2.0*pi*y)
           + cos(double(k1)*2.0*pi*x)*sin(double(k2)*2.0*pi*y);
}

int main(int argc, char* argv[])
{


    #ifdef _FFTW_OPENMP
    string threadCountInput = "-1";
    int threadIndex;
    int threadCount;
    if(not threadCountInput.empty())
    {
    threadCount = atoi(threadCountInput.c_str());
    }
    if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
    if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);

	printf("\n");
    printf("#############\n");
	printf("############# Using OpenMP With %d Threads\n",omp_get_max_threads());
	printf("#############\n");
	printf("\n");
    #endif

    long nx = 8;                  // number of panels in x-direction
    long ny = 8;                  // number of panels in y-direction

    double dx = 1.0/double(nx);   // mesh sizes 
    double dy = 1.0/double(ny);

	long i,j;

	UCLAQ::DoubleVector2d	inReal(nx,ny),    inImag(nx,ny);
	UCLAQ::DoubleVector2d	outReal(nx,ny),  outImag(nx,ny);
    UCLAQ::DoubleVector2d	invReal(nx,ny),  invImag(nx,ny);

    //
    // Declare and initialize an instance of the 2D fftw interface class
    //

	UCLAQ::fftw3_2d DFT;

	#ifdef _FFTW_OPENMP
	DFT.initialize(nx, ny , threadCount);
	#else
	DFT.initialize(nx, ny );
	#endif

	DFT.initialize(nx,ny);
	
    double x; 
    double y;

    long k1;      long k2;       // wave number indices 
    long k1Index; long k2Index;  // transform indices
    //
    // Loop over wave numbers corresponding to a Fourier basis of
    // physical space.
    // 

    printf("#### Forward Transform Results #### \n");
    double errK1K2;

    for(k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    {
    for(k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    {

    //
    // Create real and imaginary parts of (k1,k2)th Fourier mode
    //

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
		inReal(i,j) = realPartK1K2(k1,k2,x,y);
        inImag(i,j) = imagPartK1K2(k1,k2,x,y);
	}}

    //
    // Compute transform
    //

	DFT.fftw2d_forward(inReal, inImag,outReal,outImag);

    k1Index = k1 + (nx/2); // (k1Index,k2Index) = array index of the 
    k2Index = k2 + (ny/2); //  transform value of the (k1,k2)th mode. 

    //
    // The transform of (k1Index,k2Index)  basis function should be a
    // vector whose imaginary component is zero and 
    // whose real component has the value 1 in the
    // (k1Index, k2Index) location.
    //
    
    // To compute the error, subtract off what the solution should
    // and then take norms of the result 
    //

    outReal(k1Index,k2Index) = outReal(k1Index,k2Index)  - 1.0;

    errK1K2  = sqrt(outReal.dot(outReal));
    errK1K2 += sqrt(outImag.dot(outImag));

    printf("Error in forward transform of (%+ld , %+ld) mode %+10.5e \n",k1,k2, errK1K2);
    }

    }

    printf("\n#### Inverse Transform Results #### \n");
//
//  Check Inverse transform
//
    for(k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    {
    for(k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    {

    //
    // Set up transform value of (k1,k2)th Fourier mode 
    // 
    k1Index = k1 + (nx/2); 
    k2Index = k2 + (ny/2); 
    outReal.setToValue(0.0);
    outImag.setToValue(0.0);
    outReal(k1Index,k2Index) = 1.0;

    //
    // Compute the inverse transform
    //

    invReal.setToValue(0.0);
    invImag.setToValue(0.0);

    DFT.fftw2d_inverse(outReal,outImag,invReal, invImag);

    //
    // Check to see if it maps back to the (k1,k2)th Fourier mode
    //

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
		inReal(i,j) = realPartK1K2(k1,k2,x,y);
        inImag(i,j) = imagPartK1K2(k1,k2,x,y);
	}}

    // subtract off what it should be and then take norms 

    invReal -= inReal;
    invImag -= inImag;
    
    errK1K2  = sqrt(invReal.dot(invReal));
    errK1K2 += sqrt(invImag .dot(invImag ));
    printf("Error in inverse transform of (%+ld , %+ld) mode %+10.5e \n",k1,k2, errK1K2);

    }}
	return 0;
}
