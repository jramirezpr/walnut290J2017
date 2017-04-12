#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;
//
// LaplaceOp2dEigTest.cpp
//
// A test code that demonstrates the use of the RayleighChebyshev procedure
// to determine a specified number of eigenvectors and eigenvalues of
// the five point discretization to the opposite of the Laplace operator with
// homogeneous Dirichlet boundary conditions.
//
// The LaplaceOp2d class used in this test is implemented in
// 290J_Samples/LaplaceOp2dTest/LaplaceOp2d.h
//
// The test program assumes the directories 290J_2015 and 290J_Samples are organized as follows
//
// Course Directory
// ---------------------------------------------------------------
//      |                |                |            |
// Project_1 Project_2   ....        290J_Samples   290J_2015
//
//
// where the directory 290J_2015 contains the source directories :
//
//    DoubleVectorNd
//    GridFunctionNd
//    RandOpNd
//    RayleighChebyshev
//
// The command line compilation command is
//
// g++ LaplaceOp2dEigTest.cpp ../../290J_2015/RandOpNd/randomc/sfmt.cpp -std=c++11 -I../../290J_2015 -o LaplaceOp2dEigTest.exe
//
// To enable OpenMP add -fopenmp to the above command line (or the equivalent flag for your compiler).
//
// Alternately, if one is using MakeScripts
//
// The build command executed from within 290J_Samples/LaplaceOp2dEigTest is
//
// make -f LaplaceOp2dEigTest.mk release
//
// To build without OpenMP, use
//
// make -f LaplaceOp2dEigTest.mk release OpenMP=0
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Oct. 27, 2015
//

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "RayleighChebyshev/RayleighChebyshev.h"
#include "RandOpNd/RandOpDirichlet2d.h"

// Operator used in test

#include "../LaplaceOp2dTest/LaplaceOp2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{

// OpenMP setup

#ifdef _OPENMP
    int threadCount = -1;
    cout << "Enter in number of threads: ";
    cin >> threadCount;

    if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
    if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);

    printf("\n");
    printf("#############\n");
    printf("############# Using OpenMP With %d Threads\n",omp_get_max_threads());
    printf("#############\n");
    printf("\n");
#endif


//   Problem set up.

	double xMin = 0.0;
	double xMax = 1.0;

	double yMin = 0.0;
	double yMax = 1.0;

	long xPanels = 200;
	long yPanels = 200;

	double LX   = (xMax-xMin);
	double hx   = (xMax-xMin)/(double)(xPanels);

	double LY   = (yMax-yMin);
	double hy   = (yMax-yMin)/(double)(yPanels);

	long nx     = xPanels - 1;
	long ny     = yPanels - 1;


	double pi     =  3.141592653589793238;

	double alpha  = -1.0;  // Coefficient of discrete Laplace operator

    // Compute exact eigenvalues

    vector<double> exactEigValues(nx*ny,0.0);

    long index = 0;
	for(long k1 = 1; k1 <= nx; k1++)
	{
	for(long k2 = 1; k2 <= ny; k2++)
	{
	exactEigValues[index] = alpha*(
			(2.0/(hx*hx))*(cos((k1*pi*hx)/LX) - 1.0)
		  + (2.0/(hy*hy))*(cos((k2*pi*hy)/LY) - 1.0));
    index++;
	}}

    // Create sorted list of algebraically smallest to algebraically largest

    sort(exactEigValues.begin(),exactEigValues.end());

    // Instantiate Laplace Operator

    LaplaceOp2d Lop2d(alpha);

    // Instantiate Random Operator

    UCLAQ::RandOpDirichlet2d  randomOp;

    // Allocate arrays for eigenvectors and eigenvalues

    vector <UCLAQ::GridFunction2d>  eigVectors;
    vector <double>                  eigValues;

    // Declare an instance of the Raylegh-Chebyshev eigensystem procedure

    RayleighChebyshev < UCLAQ::GridFunction2d, LaplaceOp2d , UCLAQ::RandOpDirichlet2d > RCeigProcedure;
    //                    |                |                       |
    //                    |                |                       |
    //             vector class    linear operator class     randomize operator class

    RCeigProcedure.setEigDiagnosticsFlag();
    RCeigProcedure.setVerboseFlag();

    UCLAQ::GridFunction2d vTmp(xPanels,xMin,xMax,yPanels,yMin,yMax);     // A temporary vector is required as input. This vector must
                                                                         // be a non-null instance of the vector class

    double dimension           = vTmp.getDimension();
    double subspaceTol         = 2.0e-6;
    long subspaceIncrementSize = 3;
    long bufferSize            = 3;
    long eigCount              = dimension < 10 ? dimension : 10;

    RCeigProcedure.getMinEigenSystem(eigCount, subspaceTol, subspaceIncrementSize, bufferSize, vTmp,
    		Lop2d, randomOp, eigValues, eigVectors);

    printf("\n\nXXXX   RC_OperatorEig_Test Results XXXX\n\n");
    printf("Tolerance Specified : %10.5e\n\n",subspaceTol);

    printf("       Eigenvalue     Error       Relative Error \n");
    for(long  k = 0; k < eigCount; k++ )
    {
    	printf("%-5ld %-10.5e  %10.5e   %-10.5e\n", k+1, eigValues[k], abs(eigValues[k] - exactEigValues[k]),abs(eigValues[k] -exactEigValues[k])/abs(exactEigValues[k]));
    }


	printf("\nXXX Execution Completed XXXX\n");
	return 0;

}

