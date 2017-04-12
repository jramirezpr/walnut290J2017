#include <iostream>
#include <cmath>
using namespace std;
//
// LaplaceOp2dTest.cpp
//
// A test code that tests the LaplaceOp2d class which implements a five-point
// discretization of the 2d Laplace operator with homogeneous boundary conditions.
//
// The test program assumes the directories 290J_2015 and 290J_Samples are organized as follows
//
// Course Directory
// ---------------------------------------------------------------
//      |                |                |            |
// Project_1 Project_2   ....        290J_Samples   290J_2015
//
//
// where the directory 290J_2015 contains the DoubleVectorNd and GridFunctionNd
// source directories.
//
// The command line compilation command is
//
// g++ LaplaceOp2dTest.cpp -std=c++11 -I../../290J_2015 -o LaplaceOp2dTest.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/LaplaceOp2dTest is
//
// make -f LaplaceOp2dTest.mk release
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Oct. 28, 2015
//
//
#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "LaplaceOp2d.h"

int main()
{
	long xPanels;
	long yPanels;

	double    kx;
	double    ky;

	xPanels = 50;
	yPanels = 50;
	kx      = 2.0;
	ky      = 2.0;

//   Problem set up. Unit domain.

	double xMin = 0.0;
	double xMax = 1.0;

	double yMin = 0.0;
	double yMax = 1.0;

	double LX   = (xMax-xMin);
	double hx   = (xMax-xMin)/(double)(xPanels);

	double LY   = (yMax-yMin);
	double hy   = (yMax-yMin)/(double)(yPanels);

	double pi   =  3.141592653589793238;

	double alpha  = -1.0;  // Coefficient of discrete Laplace operator
	double lambda =
			alpha*(
			(2.0/(hx*hx))*(cos((kx*pi*hx)/LX) - 1.0)
		  + (2.0/(hy*hy))*(cos((ky*pi*hy)/LY) - 1.0));

	// Create test vector consisting of the exact eigenvalue of the discrete operator.

    UCLAQ::GridFunction2d       f(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d  fExact(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d  fError(xPanels,xMin,xMax,yPanels,yMin,yMax);

	double x; double y;

	// Standard method for specifying the grid values

	for(long i = 0; i <= xPanels; i++)
	{
	for(long j = 0; j <= yPanels; j++)
	{
		x = xMin + i*hx;
		y = yMin + j*hy;
		f(i,j) = sin(kx*pi*(x-xMin)/LX)*sin(ky*pi*(y-yMin)/LY);
	}}


    // Alternate method to specify the grid values using C++11 lambda function
    // capabilities and the GridFunction2D specify member function.

	f.specify(
	          [kx,ky,pi,xMin,yMin,LX,LY]                                // "captured" values used in function
	          (double x, double y)                                     // function prototype
	          {
	            return sin(kx*pi*(x-xMin)/LX)*sin(ky*pi*(y-yMin)/LY);  // function specification
	          }
	);

	// Instantiate operator

	LaplaceOp2d Lop2D(alpha);

    // Verify by checking operator on exact eigenvector of the discrete  operator.

	fExact = lambda*f;

	Lop2D.apply(f);

	// Check results(*)
	//
	// Note: The norm2() member function of the GridFunction2D class is the 2-norm
	// of the grid values scaled by the mesh width appropriately.
	//

	fError = fExact - f;

	cout << endl;
	cout << "XXXX  Laplacian 2D Operator Test Output XXXX "          << endl;

	cout << "X-Panel Count : " << xPanels << endl;
	cout << "X-Wavenumber  : " << kx << endl;

	cout << "Y-Panel Count : " << yPanels << endl;
	cout << "Y-Wavenumber  : " << ky << endl;

	cout << "L_2    Error in operator = " << fError.norm2()  << endl;
	cout << "L_Inf  Error in operator = " << fError.normInf() << endl;
}

