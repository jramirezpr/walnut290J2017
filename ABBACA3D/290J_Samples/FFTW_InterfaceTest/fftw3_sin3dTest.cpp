//
// fftw3_sin3dTest.cpp  : Routine that tests the 2D
//                        FFTW sin transform routines
//
// Specify _FFTW_OPENMP=1 on the invocation line of the makefile to enable multi-threaded execution
//
// Nov. 10, 2015
//
#include <cmath>
#include <cstdio>
using namespace std;

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin3d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"

//
// (k1,k2,k3)th orthonormal sin function (w.r.t. mesh width
// scaled inner product). 
//
double sinXsinYsinZ(long k1, long k2, long k3, double x, double y, double z)
{
    double pi = 3.1415926535897932384;
    return  sqrt(2.0)*sin(double(k1)*pi*x)*sqrt(2.0)*sin(double(k2)*pi*y)*sqrt(2.0)*sin(double(k3)*pi*z);
}

int main(int argc, char* argv[])
{
    long nx = 8;                  // number of panels in x-direction
    long ny = 7;                  // number of panels in y-direction
    long nz = 12;                 // number of panels in z-direction

    double dx = 1.0/double(nx);   // reciprocal of panel counts 
    double dy = 1.0/double(ny);   // reciprocal of panel counts
    double dz = 1.0/double(nz);   // reciprocal of panel counts

	long i,j,p;

	UCLAQ::DoubleVector3d	inF(nx-1,ny-1,nz-1);
	UCLAQ::DoubleVector3d   outF(nx-1,ny-1,nz-1);
    UCLAQ::DoubleVector3d	invF(nx-1,ny-1,nz-1);

    //
    // Declare and initialize an instance of the 2D fftw sin transform 
    // interface class. 
    // 
    //
	UCLAQ::fftw3_sin3d DFT;
	DFT.initialize(nx,ny,nz);
	
    double x; double y; double z;
    long k1;   long k2;   long k3;    // wave number indices
    //
    // Loop over wave numbers corresponding to a sin basis of
    // physical space.
    // 
    printf("#### Forward Transform Results #### \n");

    double                errK1K2K3;
    double       errK1K2K3max = 0.0;

    for(k1 = 1; k1 <= nx-1; k1++)
    {
    for(k2 = 1; k2 <= ny-1; k2++)
    {
    for(k3 = 1; k3 <= nz-1; k3++)
    {
    //
    // construct (k1,k2,k3)th sin mode
    //

	for(i=1; i <= nx-1; i++)
	{
	for(j=1; j <= ny-1; j++)
	{
	for(p=1; p <= nz-1; p++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
        z = double(p)*dz;
		inF(i-1,j-1,p-1) = sinXsinYsinZ(k1,k2,k3,x,y,z);
	}}}

    // Compute transform

	DFT.fftw3d_sin_forward(inF, outF);
	
    // The transform of (k1th,k2th,k3th) basis function should be a
    // vector whose imaginary component is zero and 
    // whose real component has the value 1 in the
    // (k1, k2,k3) location.
    
    // To compute the error, subtract off what the solution should be
    // and then take norms of the result 
    //
    
    outF(k1-1,k2-1,k3-1) = outF(k1-1,k2-1,k3-1) - 1.0;
    errK1K2K3      = sqrt(outF.dot(outF));
    
    errK1K2K3max = ( errK1K2K3max < abs(errK1K2K3) ) ? errK1K2K3: errK1K2K3max;

    //
    // Check in place transform
    //

	DFT.fftw3d_sin_forward(inF);
	
	inF(k1-1,k2-1,k3-1) = inF(k1-1,k2-1,k3-1) - 1.0;
    errK1K2K3    += sqrt(inF.dot(inF));
	
    
    //printf("Error in forward transform of (%+ld , %+ld, %+ld) mode error : %+10.5e \n",k1,k2,k3,errK1K2K3);
    }

    }}

    printf("Maximal error in forward transform of modes in range (%+d , %+d, %+d) to (%+ld , %+ld, %+ld) :   %+10.5e \n"
            ,1,1,1, nx-1,ny-1,nz-1, errK1K2K3max);


    printf("\n#### Inverse Transform Results #### \n");

    errK1K2K3max = 0.0;
//
//  Check Inverse transform
//
    for(k1 = 1; k1 <= nx-1; k1++)
    {
    for(k2 = 1; k2 <= ny-1; k2++)
    {
    for(k3 = 1; k3 <= nz-1; k3++)
    {

    // inF is the exact inverse transform 

	for(i=1; i <= nx-1; i++)
	{
	for(j=1; j <= ny-1; j++)
	{
	for(p=1; p <= nz-1; p++)
	{
        x = double(i)*dx;
        y = double(j)*dy;
        z = double(p)*dz;
		inF(i-1,j-1,p-1) = sinXsinYsinZ(k1,k2,k3,x,y,z);
	}}}

    // Set up transform value of (k1,k2)th Fourier mode 

    outF.setToValue(0.0);
    outF(k1-1,k2-1,k3-1) =    1.0;

    // Compute the inverse transform and then the error 

    DFT.fftw3d_sin_inverse(outF,invF);
    
    invF      -= inF;
    errK1K2K3  = sqrt( invF.dot( invF));
    errK1K2K3max = ( errK1K2K3max < abs(errK1K2K3) ) ? errK1K2K3: errK1K2K3max;
    
    // Check in place transform
    
    DFT.fftw3d_sin_inverse(outF);
    
    outF     -= inF;
    errK1K2K3  += sqrt(outF.dot(outF));
    
    //printf("Error in inverse transform of (%+ld , %+ld, %+ld) mode %+10.5e \n",k1,k2,k3, errK1K2K3);

    }}}


    printf("Maximal error in inverse transform of modes in range (%+d , %+d, %+d) to (%+ld , %+ld, %+ld) :   %+10.5e \n"
            ,1,1,1, nx-1,ny-1,nz-1, errK1K2K3max);

	return 0;
}
