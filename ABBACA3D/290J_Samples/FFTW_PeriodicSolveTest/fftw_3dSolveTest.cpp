#include <iostream>
#include <cmath>
#include <functional>
using namespace std;

#ifdef _FFTW_OPENMP
#include <omp.h>
#endif
//
// fftw_3dSolveTest.cpp
//
// A test code that tests the use of FFTW to create a spectral solution of the 3d
// Laplace equation in a periodic domain.
//
// The equation to be solved is 
//
// alpha*(u_xx + u_yy + u_zz) = f
//
// where the exact solution u is the function (1 - x^2 - y^2 - z^2)^q - C  scaled and centered so that
// the support of u(x,y,z) is contained within the computational domain. The constant C
// is chosen so the exact solution has zero average value, and q, the exponent, is
// a value that controls the smoothness of the exact solution is an integer >= 2.
//
// The right hand side f is determined by analytically differentiating u(x,y,z).
//
// This program depends on an installation of the FFTW3 libraries. You may need to modify
// the include paths, library paths, and library to link to an specific installation of
// FFTW3. This program was created an tested on an Ubuntu linux system in which the FFTW3
// libraries were installed using the synaptic package manager.
//
// The test program assumes the directories 290J_2015 and 290J_Samples are organized as follows
//
// Course Directory
// ---------------------------------------------------------------
//      |                |                |            |
// Project_1 Project_2   ....        290J_Samples   290J_2015
//
//
// where the directory 290J_2015 contains the DoubleVectorNd, GridFunctionNd and 
// FFTW3_InterfaceNd source directories.
//
// The command line compilation command for single threaded execution is
//
// g++ fftw_3dSolveTest.cpp -std=c++11 -I../../290J_2015 -lfftw3 -o fftw_3dSolveTest.exe
//
// for multi-threaded execution the command is
//
// g++ fftw_3dSolveTest.cpp -D_FFTW_OPENMP -fopenmp -std=c++11 -I../../290J_2015 -lfftw3_omp -lfftw3 -o fftw_3dSolveTest.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/FFTW_SinSolveTest is
//
// make -f fftw_3dSolveTest.mk release
//
// or, for multi-threaded execution, set _FFTW_OPENMP=1
//
// make -f   fftw_3dSolveTest.mk release  _FFTW_OPENMP=1
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Nov. 11, 2015
//
//
#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"


#include "FFTW3_InterfaceNd/UCLAQ_fftw3_3d.h"
#include "FFTW3_InterfaceNd/UCLAQ_FFT_Nvalues.h"

// define _OUTPUT_PLOTDATA if one wants output in vtk format

#undef  _OUTPUT_PLOTDATA

int main()
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

	long xPanels  = 43;
	long yPanels  = 34;
    long zPanels  = 23;

    long exponent = 5;  // This is the exponent for the exact solution.
                        // The exponent must be >= 2

    // Reset panel count so efficient FFT's are used; the new
    // panel count is the next larger value that's a product
    // of primes < 13

    UCLAQ::FFT_Nvalues fft_Nvalues;

    printf("Original xPanels = %ld ",xPanels);
    xPanels = fft_Nvalues.getFFT_N(xPanels);
    printf("::: New xPanels = %ld \n",xPanels);

    printf("Original yPanels = %ld ",yPanels);
    yPanels = fft_Nvalues.getFFT_N(yPanels);
    printf("::: New yPanels = %ld \n",yPanels);

    printf("Original zPanels = %ld ",zPanels);
    zPanels = fft_Nvalues.getFFT_N(zPanels);
    printf("::: New zPanels = %ld \n",zPanels);

//   Problem set up.
//
//   Use non-unit domain so that this program demonstrates
//   how to scale the transforms of derivatives appropriately
 
	double xMin   = -3.0;
	double xMax   =  2.0;

	double yMin   = -2.0;
	double yMax   =  2.0;

	double zMin   = -1.0;
	double zMax   =  2.0;

	double LX     = (xMax-xMin);
	double LY     = (yMax-yMin);
	double LZ     = (zMax-zMin);

	double pi     =  3.141592653589793238;

	double alpha  = -1.0;  // Coefficient of laplace operator
	
	// Create test vector consisting of the exact eigenvalue of the discrete operator.

    UCLAQ::GridFunction3d       f(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d       u(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d  uExact(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d  fExact(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d  uError(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d  fError(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

    f.setToValue(0.0);
    u.setToValue(0.0);

	// uFunction is the function (1-r^2))^exponent scaled and centered so that
    // the support of u(x,y,z) spans the width of the computational domain.

	double xCent  = (xMin + xMax)/2.0;
	double yCent  = (yMin + yMax)/2.0;
	double zCent  = (zMin + zMax)/2.0;
	double radius      = (xMax-xMin)/2.0;
	radius = ( radius  < (yMax-yMin)/2.0) ? radius: (yMax-yMin)/2.0;
	radius = ( radius  < (zMax-zMin)/2.0) ? radius: (zMax-zMin)/2.0;

    std::function<double(double,double,double)> uFunction = [xCent,yCent,zCent,radius,exponent](double x, double y,double z)
	{
    double r2 = ((x-xCent)*(x-xCent) + (y - yCent)*(y-yCent) + (z - zCent)*(z-zCent))/(radius*radius);
	if(r2 > 1.0) return 0.0;
    return pow(1.0 - r2,exponent);
	};

    // fFunction is the analytic evaluation of alpha*(uFunction_xx + uFunction_yy + uFunction_zz)

    std::function<double(double,double,double)> fFunction = [xCent,yCent,zCent,radius,alpha,exponent](double x,double y,double z)
	{
    double r2 = ((x-xCent)*(x-xCent) + (y - yCent)*(y-yCent) + (z - zCent)*(z-zCent))/(radius*radius);
	if(r2 > 1.0) return 0.0;
    return alpha*(2.0*exponent*pow(1.0 - r2,exponent-2)*((2.0*exponent+1.0)*r2  - 3.0))*(1.0/(radius*radius));
	};

    // Capture the exact right hand side and the exact solution

    fExact.specify(fFunction);
    uExact.specify(uFunction);

    // Initialize the FFTW3 interface routine

    UCLAQ::fftw3_3d DFT; // Discrete transform interface


    #ifdef _FFTW_OPENMP
	DFT.initialize(xPanels,yPanels,zPanels,threadCount);
	#else
	DFT.initialize(xPanels,yPanels,zPanels);
	#endif

    // Checking spectral evaluation of the alpha*(u_xx + u_yy + u_zz)


    UCLAQ::GridFunction3d dataReal(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d dataImag(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

    UCLAQ::DoubleVector3d transDataReal(xPanels,yPanels,zPanels);
    UCLAQ::DoubleVector3d transDataImag(xPanels,yPanels,zPanels);

    // Computing forward transform

    dataReal = uExact;
    dataImag.setToValue(0.0);

    DFT.fftw3d_forward(dataReal,dataImag,transDataReal,transDataImag);

    // Loop over wave numbers (kx,ky,kz) and multiply corresponding coefficients
    // of the transform by the transform of the alpha*(u_xx + u_yy + u_zz ) operator.

    long kxIndex; long kyIndex; long kzIndex;
    double opFactor;


    for(long kx = -(xPanels/2); kx <= (xPanels-1)/2; kx++)
    {
    for(long ky = -(yPanels/2); ky <= (yPanels-1)/2; ky++)
    {
    for(long kz = -(zPanels/2); kz <= (zPanels-1)/2; kz++)
    {

    kxIndex = kx + (xPanels/2);        // (kxIndex, kyIndex,kzIndex) = index of the transform
    kyIndex = ky + (yPanels/2);        //  coefficient for the (kx,ky,kz)'th mode.
    kzIndex = kz + (zPanels/2);        //

    opFactor = -alpha*( (((kx*2.0*pi)/LX)*((kx*2.0*pi)/LX))
    		         +  (((ky*2.0*pi)/LY)*((ky*2.0*pi)/LY))
					 +  (((kz*2.0*pi)/LZ)*((kz*2.0*pi)/LZ)) );

    transDataReal(kxIndex,kyIndex,kzIndex)  *= opFactor;
    transDataImag(kxIndex,kyIndex,kzIndex)  *= opFactor;
    }}}

    // Computing inverse transform

    DFT.fftw3d_inverse(transDataReal,transDataImag,dataReal,dataImag);

    f      = dataReal;

    fError = f-  fExact;

    printf("Error in spectral evaluation of alpha*(u_xx + u_yy  +u_zz)                     : %10.5e \n",fError.norm2());

    // Solving the equation alpha*(u_xx + u_yy + u_zz) = f using a spectral discretization
    //
    // This is a singular problem, and the solution is only defined up to a constant.
    // Since the problem is self-adjoint and the right hand side is constructed by
    // differentiating a periodic function, the right hand side will automatically
    // be orthogonal to constant values, and hence a solution exists.
    // To render the computed solution unique, the spectral solution
    // returned is set to have zero average value -- and the exact solution used
    // to evaluate the error must be modified so that it has zero average as well.
    //

    // Computing forward transform

    dataReal = fExact;
    dataImag.setToValue(0.0);

    DFT.fftw3d_forward(dataReal,dataImag,transDataReal,transDataImag);

    // Loop over wave numbers and multiply corresponding coefficients
    // of the transform by the inverse of transform of the alpha*(u_xx + u_yy + u_zz ) operator.


    for(long kx = -(xPanels/2); kx <= (xPanels-1)/2; kx++)
    {
    for(long ky = -(yPanels/2); ky <= (yPanels-1)/2; ky++)
    {
    for(long kz = -(zPanels/2); kz <= (zPanels-1)/2; kz++)
    {

    kxIndex = kx + (xPanels/2);        // (kxIndex, kyIndex,kzIndex) = index of the transform
    kyIndex = ky + (yPanels/2);        //  coefficient for the (kx,ky,kz)'th mode.
    kzIndex = kz + (zPanels/2);        //


    if((kx == 0)&&(ky == 0)&&(kz ==0)) // OpFactor = 0 when kx=ky=kz=0 to create a solution with zero constant component
    {
    opFactor = 0;
	}
    else
    {
    opFactor = 1.0/(-alpha*( (((kx*2.0*pi)/LX)*((kx*2.0*pi)/LX))
    		         +  (((ky*2.0*pi)/LY)*((ky*2.0*pi)/LY))
					 +  (((kz*2.0*pi)/LZ)*((kz*2.0*pi)/LZ)) ));
    }

    transDataReal(kxIndex,kyIndex,kzIndex)  *= opFactor;
    transDataImag(kxIndex,kyIndex,kzIndex)  *= opFactor;
    }}}

    // Computing inverse transform

    DFT.fftw3d_inverse(transDataReal,transDataImag,dataReal,dataImag);

    u = dataReal;

    // Fix up exact solution to have zero average value. Using
    // the trapezoidal integration because of it's spectral accuracy
    // when applied to periodic functions.

    double uAverage = uExact.getTrapezoidalAverage();
    uExact.addValue(-uAverage);

    // Compare computed with exact

    uError = u - uExact;

    printf("Error in spectral evaluation of the solution to alpha*(u_xx + u_yy + u_zz) = f : %10.5e \n",uError.norm2());

#ifdef _OUTPUT_PLOTDATA
    UCLAQ::GridFunction3dUtility gUtility;
    gUtility.outputDataToVTKfile(uError,"uError.vtk","uError");
    gUtility.outputDataToVTKfile(u,"u.vtk","u");
    gUtility.outputDataToVTKfile(uExact,"uExact.vtk","uExact");
#endif
}

