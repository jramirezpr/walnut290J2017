#include <iostream>
#include <cmath>
#include <functional>
using namespace std;
//
// fftw_2dSolveTest.cpp
//
// A test code that tests the use of FFTW to create a spectral solution of the 2d
// Laplace equation in a periodic domain.
//
// The equation to be solved is 
//
// alpha*(u_xx + u_yy) = f
//
// where the exact solution u is the function (1 - x^2 - y^2)^q - C  scaled and centered so that
// the support of u(x,y) is contained within the computational domain. The constant C
// is chosen so the exact solution has zero average value, and q, the exponent, is
// a value that controls the smoothness of the exact solution is an integer >= 2.
//
// The right hand side f is determined by analytically differentiating u(x,y).
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
// The command line compilation command is
//
// g++ fftw_2dSolveTest.cpp -std=c++11 -I../../290J_2015 -lfftw3 -o fftw_2dSolveTest.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/fftw_PeriodicSolveTest is
//
// make -f fftw_2dSolveTest.mk release
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Nov. 11, 2015
//
//
#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"


#include "FFTW3_InterfaceNd/UCLAQ_fftw3_2d.h"
#include "FFTW3_InterfaceNd/UCLAQ_FFT_Nvalues.h"

// define _OUTPUT_PLOTDATA if one wants output data that can be plotted using GNUplot

#undef  _OUTPUT_PLOTDATA

int main()
{
	long xPanels  = 46;
	long yPanels  = 28;

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

//   Problem set up.
//
//   Use non-unit domain so that this program demonstrates
//   how to scale the transforms of derivatives appropriately
 
	double xMin   = -3.0;
	double xMax   =  2.0;

	double yMin   = -2.0;
	double yMax   =  2.0;

	double LX     = (xMax-xMin);
	double LY     = (yMax-yMin);

	double pi     =  3.141592653589793238;

	double alpha  = -1.0;  // Coefficient of laplace operator
	
	// Create test vector consisting of the exact eigenvalue of the discrete operator.

    UCLAQ::GridFunction2d       f(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d       u(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d  uExact(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d  fExact(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d  uError(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d  fError(xPanels,xMin,xMax,yPanels,yMin,yMax);

    f.setToValue(0.0);
    u.setToValue(0.0);

	// uFunction is the function (1-r^2))^exponent scaled and centered so that
    // the support of u(x,y) spans the width of the computational domain.

	double xCent  = (xMin + xMax)/2.0;
	double yCent  = (yMin + yMax)/2.0;
	double radius = ( (xMax-xMin)  < (yMax-yMin)) ? (xMax-xMin)/2.0 : (yMax-yMin)/2.0;

    std::function<double(double,double)> uFunction = [xCent,yCent,radius,exponent](double x, double y)
	{
    double r2 = ((x-xCent)*(x-xCent) + (y - yCent)*(y-yCent))/(radius*radius);
	if(r2 > 1.0) return 0.0;
    return pow(1.0 - r2,exponent);
	};

    // fFunction is the analytic evaluation of alpha*(uFunction_xx + uFunction_yy)


    std::function<double(double,double)> fFunction = [xCent,yCent,radius,alpha,exponent](double x,double y)
	{
    double r2 = ((x-xCent)*(x-xCent) + (y - yCent)*(y-yCent))/(radius*radius);
	if(r2 > 1.0) return 0.0;

    return alpha*(4.0*exponent*pow(1.0 - r2,exponent-2)*(exponent*r2  - 1.0))*(1.0/(radius*radius));
	};

    // Capture the exact right hand side and the exact solution

    fExact.specify(fFunction);
    uExact.specify(uFunction);



#ifdef _OUTPUT_PLOTDATA
    UCLAQ::GridFunction2dUtility gUtility2d;
    gUtility2d.outputToGNUplot(uExact,"uExact.dat");
    gUtility2d.outputToGNUplot(fExact,"fExact.dat");
#endif

    // Initialize the FFTW3 interface routine

    UCLAQ::fftw3_2d DFT;               // Discrete transform interface
    DFT.initialize(xPanels,yPanels);   // Specify the number of grid panels!

    // Checking spectral evaluation of the alpha*(u_xx + u_yy)

    UCLAQ::GridFunction2d dataReal(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d dataImag(xPanels,xMin,xMax,yPanels,yMin,yMax);

    UCLAQ::DoubleVector2d transDataReal(xPanels,yPanels);
    UCLAQ::DoubleVector2d transDataImag(xPanels,yPanels);

    // Computing forward transform

    dataReal = uExact;
    dataImag.setToValue(0.0);

    DFT.fftw2d_forward(dataReal,dataImag,transDataReal,transDataImag);

    // Loop over wave numbers (kx,ky) and multiply corresponding coefficients
    // of the transform by the transform of the alpha*(u_xx + u_yy) operator.

    long kxIndex; long kyIndex;
    double opFactor;

    for(long kx = -(xPanels/2); kx <= (xPanels-1)/2; kx++)
    {
    for(long ky = -(yPanels/2); ky <= (yPanels-1)/2; ky++)
    {
    kxIndex = kx + (xPanels/2);        // (kxIndex, kyIndex) = index of the transform
    kyIndex = ky + (yPanels/2);        //  coefficient for the (kx,ky)'th mode.

    opFactor = -alpha*( (((kx*2.0*pi)/LX)*((kx*2.0*pi)/LX))
    		         +  (((ky*2.0*pi)/LY)*((ky*2.0*pi)/LY)) );

    transDataReal(kxIndex,kyIndex)  *= opFactor;
    transDataImag(kxIndex,kyIndex)  *= opFactor;
    }}

    // Computing inverse transform

    DFT.fftw2d_inverse(transDataReal,transDataImag,dataReal,dataImag);

    f = dataReal;

    fError = f-  fExact;

    printf("Error in spectral evaluation of alpha*(u_xx + u_yy)                     : %10.5e \n",fError.norm2());

#ifdef _OUTPUT_PLOTDATA
    gUtility2d.outputToGNUplot(f,"f.dat");
    gUtility2d.outputToGNUplot(fExact,"fError.dat");
#endif

    // Solving the equation alpha*(u_xx + u_yy) = f using a spectral discretization
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

    DFT.fftw2d_forward(dataReal,dataImag,transDataReal,transDataImag);

    // Loop over wave numbers (kx,ky) and multiply corresponding coefficients
    // of the transform by the inverse transform of the alpha*(u_xx + u_yy) operator.

    for(long kx = -(xPanels/2); kx <= (xPanels-1)/2; kx++)
    {
    for(long ky = -(yPanels/2); ky <= (yPanels-1)/2; ky++)
    {
    kxIndex = kx + (xPanels/2);        // (kxIndex, kyIndex) = index of the transform
    kyIndex = ky + (yPanels/2);        //  coefficient for the (kx,ky)'th mode.


    if((kx == 0)&&(ky == 0)) // OpFactor = 0 when kx=ky=0 to create a solution with zero constant component
    {
    opFactor = 0;
	}
    else
    {
    opFactor = 1.0/(-alpha*( (((kx*2.0*pi)/LX)*((kx*2.0*pi)/LX))
    		                +(((ky*2.0*pi)/LY)*((ky*2.0*pi)/LY)) ));
    }

    transDataReal(kxIndex,kyIndex)  *= opFactor;
    transDataImag(kxIndex,kyIndex)  *= opFactor;
    }}

    // Computing inverse transform

    DFT.fftw2d_inverse(transDataReal,transDataImag,dataReal,dataImag);

    u = dataReal;

    // Fix up exact solution to have zero average value

    // Fix up exact solution to have zero average value. Using
    // the trapezoidal intergration because of it's spectral accuracy
    // when applied to periodic functions.

    double uAverage = uExact.getTrapezoidalAverage();
    uExact.addValue(-uAverage);

    // Compare computed with exact

    uError = u - uExact;

    printf("Error in spectral evaluation of the solution to alpha*(u_xx + u_yy) = f : %10.5e \n",uError.norm2());


#ifdef _OUTPUT_PLOTDATA
    gUtility2d.outputToGNUplot(u,"u.dat");
    gUtility2d.outputToGNUplot(uExact,"uError.dat");
#endif

}

