#include <iostream>
#include <cmath>
#include <functional>
using namespace std;
//
// fftw_sin2dSolveTest.cpp
//
// A test code that tests the use of FFTW to create a spectral solution of the 2d
// Laplace equation with homogeneous Dirichlet boundary conditions.
//
// The equation to be solved is 
//
// alpha*(u_xx + u_yy) = f
//
// where the exact solution u is the function (1-(x^2 + y^2))^q scaled and centered so that
// the support of u(x,y) is contained within the computational domain. The right hand
// side f is determined by analytically differentiating u(x,y).
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
// g++ fftw_sin2dSolveTest.cpp -std=c++11 -I../../290J_2015 -lfftw3 -o fftw_sin2dSolveTest.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/FFTW_SinSolveTest is
//
// make -f fftw_sin2dSolveTest.mk release
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Nov. 10, 2015
//
//
#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"


#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin2d.h"
#include "FFTW3_InterfaceNd/UCLAQ_FFT_Nvalues.h"

// define _OUTPUT_PLOTDATA if one wants output data that can be plotted using GNUplot

#undef  _OUTPUT_PLOTDATA

int main()
{
	long xPanels  = 43;
	long yPanels  = 34;

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

    UCLAQ::fftw3_sin2d DFT;            // Discrete sin transform interface
    DFT.initialize(xPanels,yPanels);   // Specify the number of grid panels!

    // Checking spectral evaluation of the alpha*(u_xx + u_yy)

    UCLAQ::DoubleVector2d uTransform(xPanels-1,yPanels-1);
    uTransform.setToValue(0.0);

    DFT.fftw2d_sin_forward(uExact,uTransform); // Note: We can use a GridFunction2d as an argument since it's extended
                                               // from a DoubleVector2d

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    	uTransform(kx-1,ky-1) *= -alpha*( (((kx*pi)/LX)*((kx*pi)/LX)) +  (((ky*pi)/LY)*((ky*pi)/LY)) ) ;
    }}

    DFT.fftw2d_sin_inverse(uTransform,f);      // Note: We can use a GridFunction2d as an argument since it's extended
                                               // from a DoubleVector2d

    fError = f-  fExact;

    printf("Error in spectral evaluation of alpha*(u_xx + u_yy)                     : %10.5e \n",fError.norm2());

#ifdef _OUTPUT_PLOTDATA
    gUtility2d.outputToGNUplot(f,"f.dat");
    gUtility2d.outputToGNUplot(fExact,"fError.dat");
#endif

    // Solving the equation alpha*(u_xx + u_yy) = f using a spectral discretization

    UCLAQ::DoubleVector2d fTransform(xPanels-1,yPanels-1);
    fTransform.setToValue(0.0);

    DFT.fftw2d_sin_forward(fExact,fTransform); // Note: We can use a GridFunction2d as an argument since it's extended
                                               // from a DoubleVector2d

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    	fTransform(kx-1,ky-1) *= -1.0/(alpha*( (((kx*pi)/LX)*((kx*pi)/LX)) +  (((ky*pi)/LY)*((ky*pi)/LY))));
    }}

    DFT.fftw2d_sin_inverse(fTransform,u);      // Note: We can use a GridFunction2d as an argument since it's extended
                                               // from a DoubleVector2d

    uError = u - uExact;

    printf("Error in spectral evaluation of the solution to alpha*(u_xx + u_yy) = f : %10.5e \n",uError.norm2());


#ifdef _OUTPUT_PLOTDATA
    gUtility2d.outputToGNUplot(u,"u.dat");
    gUtility2d.outputToGNUplot(uExact,"uError.dat");
#endif
}

