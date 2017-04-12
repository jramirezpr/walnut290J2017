//
// fftw3_sin3d.cpp  : Routine that tests the 3D FFTW routines
//
// Nov. 10, 2015
//

#ifdef _FFTW_OPENMP
#include <omp.h>
#endif


#include "FFTW3_InterfaceNd/UCLAQ_fftw3_3d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "Timing/ClockIt.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <cstdlib>

using namespace std;

// The following two routines compute the real and imaginary
// part of the (k1,k2,k3) discrete Fourier basis function
//
//  exp(2*pi*I*k1*x)*exp(2*pi*I*k2*y)*exp(2*pi*I*k3*z)
//
//  for (x,y,z) in [0,1]X[0,1]X[0,1]
/*
(simplify(convert(exp(2*Pi*I*k1*x)*exp(2*Pi*I*k2*y)*exp(2*Pi*I*k3*z),trig)));
cos(2*Pi*k1*x)*cos(2*Pi*k2*y)*cos(2*Pi*k3*z)
-cos(2*Pi*k1*x)*sin(2*Pi*k2*y)*sin(2*Pi*k3*z)
+I*sin(2*Pi*k1*x)*cos(2*Pi*k2*y)*cos(2*Pi*k3*z)
+I*cos(2*Pi*k1*x)*cos(2*Pi*k2*y)*sin(2*Pi*k3*z)
+I*cos(2*Pi*k1*x)*sin(2*Pi*k2*y)*cos(2*Pi*k3*z)
-I*sin(2*Pi*k1*x)*sin(2*Pi*k2*y)*sin(2*Pi*k3*z)
-sin(2*Pi*k1*x)*cos(2*Pi*k2*y)*sin(2*Pi*k3*z)
-sin(2*Pi*k1*x)*sin(2*Pi*k2*y)*cos(2*Pi*k3*z)
*/
double realPartK1K2K3(long k1, long k2,long k3, double x, double y, double z)
{
    double Pi = 4.0*atan(1.0);

return cos(2*Pi*k1*x)*cos(2*Pi*k2*y)*cos(2*Pi*k3*z) -cos(2*Pi*k1*x)*sin(2*Pi*k2*y)*sin(2*Pi*k3*z)
	  -sin(2*Pi*k1*x)*cos(2*Pi*k2*y)*sin(2*Pi*k3*z) -sin(2*Pi*k1*x)*sin(2*Pi*k2*y)*cos(2*Pi*k3*z);
}

double imagPartK1K2K3(long k1, long k2, long k3, double x, double y, double z)
{
    double Pi = 4.0*atan(1.0);
    return   sin(2*Pi*k1*x)*cos(2*Pi*k2*y)*cos(2*Pi*k3*z)+cos(2*Pi*k1*x)*cos(2*Pi*k2*y)*sin(2*Pi*k3*z)
             +cos(2*Pi*k1*x)*sin(2*Pi*k2*y)*cos(2*Pi*k3*z)-sin(2*Pi*k1*x)*sin(2*Pi*k2*y)*sin(2*Pi*k3*z);
}


int main(int argc, char* argv[])
{
    string threadCountInput = "-1";

	#ifdef _FFTW_OPENMP
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


    double errFwd = 0.0;
    double errBck = 0.0;


    long nx = 10;                  // number of panels in x-direction
    long ny = 8;                  // number of panels in y-direction
    long nz = 9;                  // number of panels in z-direction


    double dx = 1.0/double(nx);   // mesh sizes 
    double dy = 1.0/double(ny);
    double dz = 1.0/double(nz);

	long i,j,k;

	UCLAQ::DoubleVector3d	inReal(nx,ny,nz),    inImag(nx,ny,nz);
	UCLAQ::DoubleVector3d	outReal(nx,ny,nz),  outImag(nx,ny,nz);
    UCLAQ::DoubleVector3d	invReal(nx,ny,nz),  invImag(nx,ny,nz);

    //
    // Declare and initialize an instance of the 3D fftw interface class
    //

	UCLAQ::fftw3_3d DFT;

    #ifdef _FFTW_OPENMP
	DFT.initialize(nx, ny , nz,threadCount);
	#else
	DFT.initialize(nx, ny , nz);
	#endif
	
    double x; 
    double y;
    double z;

    long k1;      long k2;      long k3;      // wave number indices
    long k1Index; long k2Index; long k3Index; // transform indices
    //
    // Loop over wave numbers corresponding to a Fourier basis of
    // physical space.
    // 


    double errK1K2K3;

    for(k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    {
    for(k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    {
    for(k3 = -(nz/2); k3 <= (nz-1)/2; k3++)
    {

    // Create real and imaginary parts of (k1,k2,k3)th Fourier mode

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
    for(k=0; k < nz; k++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
        z = double(k)*dz;
		inReal(i,j,k) = realPartK1K2K3(k1,k2,k3,x,y,z);
        inImag(i,j,k) = imagPartK1K2K3(k1,k2,k3,x,y,z);
	}}}


    //
    // Compute transform
    //

	DFT.fftw3d_forward(inReal, inImag,outReal,outImag);

    k1Index = k1 + (nx/2); // (k1Index,k2Index,k3Index) = array index of the
    k2Index = k2 + (ny/2); //  transform value of the (k1,k2,k3)th mode.
    k3Index = k3 + (nz/2); //

    //
    // The transform of (k1,k2,k3)th basis function should be a
    // vector whose imaginary component is zero and 
    // whose real component has the value 1 in the
    // (k1Index, k2Index,k3Index) location.
    //
    
    // To compute the error, subtract off what the solution should
    // and then take norms of the result 
    //

    outReal(k1Index,k2Index,k3Index) = outReal(k1Index,k2Index,k3Index) - 1.0;

    errK1K2K3  = sqrt(outReal.dot(outReal));
    errK1K2K3 += sqrt(outImag.dot(outImag));

    errFwd = (errFwd < errK1K2K3) ? errK1K2K3 : errFwd;
    // printf("Error in forward transform of (%+ld , %+ld, %+ld) mode %+10.5e \n",k1,k2,k3,errK1K2K3);
    }

    }}

//
//  Check Inverse transform
//
    for(k1 = -(nx/2); k1 <= (nx-1)/2; k1++)
    {
    for(k2 = -(ny/2); k2 <= (ny-1)/2; k2++)
    {
    for(k3 = -(nz/2); k3 <= (nz-1)/2; k3++)
    {
    //
    // Set up transform value of (k1,k2,k3)th Fourier mode
    // 
    k1Index = k1 + (nx/2); 
    k2Index = k2 + (ny/2); 
    k3Index = k3 + (nz/2);
    outReal.setToValue(0.0);
    outImag.setToValue(0.0);

    outReal(k1Index,k2Index,k3Index) = 1.0;

    //
    // Compute the inverse transform
    //

    invReal.setToValue(0.0);
    invImag.setToValue(0.0);

    DFT.fftw3d_inverse(outReal,outImag,invReal, invImag);

    //
    // Check to see if it maps back to the (k1,k2,k3)th Fourier mode
    //

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
    for(k=0; k < nz; k++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
        z = double(k)*dz;
		inReal(i,j,k) = realPartK1K2K3(k1,k2,k3,x,y,z);
        inImag(i,j,k) = imagPartK1K2K3(k1,k2,k3,x,y,z);
	}}}

    // subtract off what it should be and then take norms 

    invReal -= inReal;
    invImag -= inImag;
    
    errK1K2K3  = sqrt(invReal.dot(invReal));
    errK1K2K3 += sqrt(invImag .dot(invImag));

    errBck = (errBck < errK1K2K3) ? errK1K2K3 : errBck;
    //printf("Error in inverse transform of (%+ld , %+ld, %+ld) mode %+10.5e \n",k1,k2,k3, errK1K2K3);
    }}}

    printf("Max Error in forward transform of  %+10.5e \n",errFwd);
	printf("Max Error in inverse transform of  %+10.5e \n",errBck);

	//
	// Timing for 3D FFT of 100x100x100
	//

	long kTimes = 20;
	ClockIt clockIt;

    nx = 100;                  // number of panels in x-direction
    ny = 100;                  // number of panels in y-direction
    nz = 100;                  // number of panels in z-direction

    dx = 1.0/double(nx);   // mesh sizes
    dy = 1.0/double(ny);
    dz = 1.0/double(nz);


	inReal.initialize(nx,ny,nz);    inImag.initialize(nx,ny,nz);
	outReal.initialize(nx,ny,nz);  outImag.initialize(nx,ny,nz);
    invReal.initialize(nx,ny,nz);  invImag.initialize(nx,ny,nz);

    // Test Mode

    k1 = (nx-1)/4;
    k2 = (ny-1)/4;
    k3 = (nz-1)/4;

    #ifdef _FFTW_OPENMP
	DFT.initialize(nx, ny , nz,threadCount);
	#else
	DFT.initialize(nx, ny , nz);
	#endif

    //
    // Create real and imaginary parts of Fourier mode
    //

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
    for(k=0; k < nz; k++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
        z = double(k)*dz;
		inReal(i,j,k) = realPartK1K2K3(k1,k2,k3,x,y,z);
        inImag(i,j,k) = imagPartK1K2K3(k1,k2,k3,x,y,z);
	}}}

    //
    // Compute transform
    //

    clockIt.start();
    for(long count = 0; count < kTimes; count++)
    {
	DFT.fftw3d_forward(inReal, inImag,outReal,outImag);
	}
	clockIt.stop();

	printf("Forward_Time_MS : %-10.5f \n",clockIt.getMilliSecElapsedTime()/double(kTimes));

    k1Index = k1 + (nx/2); //  (k1Index,k2Index,k3Index) = array index of the
    k2Index = k2 + (ny/2); //  transform value of the (k1,k2,k3)th mode.
    k3Index = k3 + (nz/2); //

    //
    // The transform of (k1,k2,k3)th basis function should be a
    // vector whose imaginary component is zero and
    // whose real component has the value nx*ny*nz in the
    // (k1Index, k2Index,k3Index) location.
    //

    // To compute the error, subtract off what the solution should
    // and then take norms of the result
    //

    outReal(k1Index,k2Index,k3Index) = outReal(k1Index,k2Index,k3Index) - 1.0;

    errK1K2K3  = sqrt(outReal.dot(outReal));
    errK1K2K3 += sqrt(outImag.dot(outImag));

    printf("Error in forward transform of (%+ld , %+ld, %+ld) mode %+10.5e \n",k1,k2,k3,errK1K2K3);


    k1Index = k1 + (nx/2);
    k2Index = k2 + (ny/2);
    k3Index = k3 + (nz/2);
    outReal.setToValue(0.0);
    outImag.setToValue(0.0);

    outReal(k1Index,k2Index,k3Index) = 1.0;

    //
    // Compute the inverse transform
    //

    invReal.setToValue(0.0);
    invImag.setToValue(0.0);

    clockIt.start();
    for(long count = 0; count < kTimes; count++)
    {
    DFT.fftw3d_inverse(outReal,outImag,invReal, invImag);
    }
	clockIt.stop();
	printf("Inverse_Time_MS : %-10.5f \n",clockIt.getMilliSecElapsedTime()/double(kTimes));

    //
    // Check to see if it maps back to the (k1,k2,k3)th Fourier mode
    //

	for(i=0; i < nx; i++)
	{
	for(j=0; j < ny; j++)
	{
    for(k=0; k < nz; k++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
        z = double(k)*dz;
		inReal(i,j,k) = realPartK1K2K3(k1,k2,k3,x,y,z);
        inImag(i,j,k) = imagPartK1K2K3(k1,k2,k3,x,y,z);
	}}}

    // subtract off what it should be and then take norms

    invReal -= inReal;
    invImag -= inImag;

    errK1K2K3  = sqrt(invReal.dot(invReal));
    errK1K2K3 += sqrt(invImag.dot(invImag));
    printf("Error in inverse transform of (%+ld , %+ld, %+ld) mode %+10.5e \n",k1,k2,k3, errK1K2K3);

	return 0;
}
