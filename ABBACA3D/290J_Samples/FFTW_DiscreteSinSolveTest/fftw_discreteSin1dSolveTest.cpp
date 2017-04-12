#include <iostream>
#include <cmath>
#include <functional>
using namespace std;
//
// fftw_discreteSin1dSolveTest.cpp
//
// A test code that tests the use of FFTW to create a solution of the 1d
// discrete Laplace equation with homogeneous Dirichlet boundary conditions.
//
// The equation to be solved is 
//
// alpha*Delta_h u   = f
//
// where Delta_h is the three point approximation to the discrete Laplace operator.
//  
// the exact solution u is the function (1-x^2)^q scaled and centered so that
// the support of u(x) spans the width of the computational domain. The right hand
// side f is determined by applying alpha*Delta_h to u.
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
// g++ fftw_discreteSin1dSolveTest.cpp -std=c++11 -I../../290J_2015 -lfftw3 -o fftw_discreteSin1dSolveTest.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/FFTW_DiscreteSinSolveTest is
//
// make -f fftw_discreteSin1dSolveTest.mk release
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Nov. 17, 2015
//
//
#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin1d.h"

// Use a relative path to LaplaceOp1d so that additions to the 
// include path don't have to be set
 

#include "../LaplaceOp1dTest/LaplaceOp1d.h"

// define _OUTPUT_PLOTDATA if one wants output data that can be plotted using GNUplot

#undef _OUTPUT_PLOTDATA

int main()
{
	long xPanels  = 20;
    long exponent = 5;  // This is the exponent for the exact solution and must be >= 2

//   Problem set up.
//
//   Use non-unit domain so that this program demonstrates
//   how to scale the transforms of derivatives appropriately
 
	double xMin = -3.0;
	double xMax =  2.0;

	double LX   = (xMax-xMin);

	double pi   =  3.141592653589793238;

	double alpha  = -1.0;  // Coefficient of Laplace operator


	// Instantiate the Laplace operator

    LaplaceOp1d laplaceOp(alpha);

	// Create test vector consisting of the exact eigenvalue of the discrete operator.

    UCLAQ::GridFunction1d       f(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d       u(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d  uExact(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d  fExact(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d  uError(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d  fError(xPanels,xMin,xMax);

    f.setToValue(0.0);
    u.setToValue(0.0);

	// uFunction is the function (1-x^2)^exponent scaled and centered so that
    // the support of u(x) spans the width of the computational domain.

    std::function<double(double)> uFunction = [xMin,xMax,exponent](double x)
	{
    double xBar  = (xMin+xMax)/2.0;
    double xStar = (x-xBar)/((xMax-xMin)/2.0);
	if(x < xMin) return 0.0;
	if(x > xMax) return 0.0;
    return pow(1.0 - xStar*xStar,exponent);
	};

  
    // Capture the exact solution and create the right hand side 

    uExact.specify(uFunction);

    // Create right hand side by applying alpha*Delta_h

    fExact = uExact;
    laplaceOp.apply(fExact);

    // Set boundary values fExact to 0 since they 
    // will be ignored by the sin transform solution
    // procedure 

    fExact.setBoundaryValues(0.0);


#ifdef _OUTPUT_PLOTDATA
    UCLAQ::GridFunction1dUtility gUtility1d;
    gUtility1d.outputToGNUplot(uExact,"uExact.dat");
    gUtility1d.outputToGNUplot(fExact,"fExact.dat");
#endif

    // Initialize the FFTW3 interface routine

    UCLAQ::fftw3_sin1d DFT;    // Discrete sin transform interface
    DFT.initialize(xPanels);   // Specify the number of grid panels!

    // Checking  evaluation of the alpha*Delta_h u

    UCLAQ::DoubleVector1d uTransform(xPanels-1);
    uTransform.setToValue(0.0);

    DFT.fftw1d_sin_forward(uExact,uTransform); // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    double hx = f.getHx();
    double gamma_k;

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
        gamma_k = alpha*(-4.0/(hx*hx))*sin(kx*pi*(hx/2.0)/LX)*sin(kx*pi*(hx/2.0)/LX);

    	uTransform(kx-1) *= gamma_k;
    }

    DFT.fftw1d_sin_inverse(uTransform,f);      // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    fError = f-  fExact;

    printf("Error in evaluation of alpha*Delta_h u                      : %10.5e \n",fError.norm2());


#ifdef _OUTPUT_PLOTDATA
    gUtility1d.outputToGNUplot(f,"f.dat");
    gUtility1d.outputToGNUplot(fExact,"fError.dat");
#endif


    // Solving the equation alpha*u_xx = f using a spectral discretization

    UCLAQ::DoubleVector1d fTransform(xPanels-1);
    fTransform.setToValue(0.0);

    DFT.fftw1d_sin_forward(fExact,fTransform); // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    	gamma_k = alpha*(-4.0/(hx*hx))*sin(kx*pi*(hx/2.0)/LX)*sin(kx*pi*(hx/2.0)/LX);
    	fTransform(kx-1) *= 1.0/gamma_k;
    }

    DFT.fftw1d_sin_inverse(fTransform,u);      // Note: We can use a GridFunction1d as an argument since it's extended
                                               // from a DoubleVector1d

    uError = u - uExact;

    printf("Error in evaluation  of the solution to alpha*Delta_h u = f : %10.5e \n",uError.norm2());


#ifdef _OUTPUT_PLOTDATA
    gUtility1d.outputToGNUplot(u,"u.dat");
    gUtility1d.outputToGNUplot(uExact,"uError.dat");
#endif


}

