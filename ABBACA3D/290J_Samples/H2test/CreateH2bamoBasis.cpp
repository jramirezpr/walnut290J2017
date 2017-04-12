#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
using namespace std;
//
// CreateH2bamoBasis.cpp
//
// A test code that demonstrates the use of a finite difference
// procedure to create an orbital basis consisting of the 
// eigenfunctions of the single particle Schroedinger operator 
// with the potential induced by two nucleii located at A and B, 
// points in R^3 (Bare Atomic Molecular Orbitals). 
//
// The basis set consists of the first M eigenfunctions (M specified) 
// of the operator  
//
//  [ -(1/2)*DELTA  + V_epsilon(x) ]
//
// where the Laplace operator DELTA is a high order finite difference approximation
// with homogeneous boundary conditions and V(x) is a smoothed nuclear potential 
// of unit strength with centers at (xA,yA,zA) and (xB,yB,zB). The smoothing
// radius is epsilon. 
//

#include "Timing/ClockIt.h"

// Problem specific classes

// Single particle operator combining a high order 
// approximatio of the Laplace operator with 
// a smoothed nuclear potential

#include "SmoothedH2op.h" 

// Classes for basis set construction as eigenfunctions of the
// single particle operator (Bare Atomic Molecular Orbitals)

#include "RandOpNd/RandOpDirichlet3d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "RayleighChebyshev/RayleighChebyshev.h"

#include "Qdomain/QdomainBasisFunctions.h"


int main()
{
	ClockIt clockIt;

	string basisFileName     = "H2bamoBasis.dat";
	string potentialFilename = "H2nuclearPotential.dat";


	long panelCount               = 50;
	double domainWidth            = 5.0;
	double xDomainWidthExtension  = 0.5;
    double epsilon                = 0.5;

	int threadCount = 8;

	#ifdef _OPENMP
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

	long xPanels = panelCount; long yPanels = panelCount; long zPanels = panelCount;

	double xMin = -domainWidth - xDomainWidthExtension;
	double xMax =  domainWidth + xDomainWidthExtension;

	double yMin = -domainWidth; double yMax = domainWidth;
	double zMin = -domainWidth; double zMax = domainWidth;

	// Location of nuclear charges (using ground state separation)

	vector<double> A = {-0.70055,0.0,0.0};
	vector<double> B = { 0.70055,0.0,0.0};

//##################################################################
//   Initializing smoothed single particle H2 Op
//##################################################################

	SmoothedH2op sOp(A,B, epsilon);

//##################################################################
//   Creating basis elements
//##################################################################

	// Instantiate Random Operator

    UCLAQ::RandOpDirichlet3d  randomOp;

    // Allocate arrays for eigenvectors and eigenvalues

    vector <UCLAQ::GridFunction3d>  eigVectors;
    vector <double>                  eigValues;

    // Declare an instance of the Raylegh-Chebyshev eigensystem procedure

    RayleighChebyshev < UCLAQ::GridFunction3d, SmoothedH2op , UCLAQ::RandOpDirichlet3d > RCeigProcedure;
    //                    |                      |                       |
    //                    |                      |                       |
    //             vector class        linear operator class     randomize operator class

    RCeigProcedure.setEigDiagnosticsFlag();
    RCeigProcedure.setVerboseFlag();

    UCLAQ::GridFunction3d vTmp(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);     // A temporary vector is required as input. This vector must
                                                                                           // be a non-null instance of the vector class

    long basisSize             = 20;
    long dimension             = vTmp.getDimension();
    double subspaceTol         = 2.0e-5;
    long subspaceIncrementSize = basisSize;
    long bufferSize            = 5;
    long eigCount              = dimension < 10 ? dimension : 10;


    printf("Creating Eigensystem : ");
    fflush(stdout);

    clockIt.start();

    RCeigProcedure.getMinEigenSystem(eigCount, subspaceTol, subspaceIncrementSize, bufferSize, vTmp, sOp, randomOp, eigValues, eigVectors);

    clockIt.stop();
	printf("\n\n");
	printf("Done --- Time (sec) : %10.5f \n\n",clockIt.getSecElapsedTime());

    printf("\n\nXXXX   RC_OperatorEig_Test Results XXXX\n\n");
    printf("Tolerance Specified : %10.5e\n\n",subspaceTol);

    printf("       Eigenvalue    \n");
    for(long  k = 0; k < eigCount; k++ )
    {
    	printf("%-5ld %-10.5e   \n", k+1, eigValues[k]);
    }

    double eigNorm = 0.0;

    for(long  k = 0; k < eigCount; k++ )
    {
    	eigNorm = sqrt(abs(eigVectors[k].scaledDot(eigVectors[k])));
        eigVectors[k] /= eigNorm;
    }

//##################################################################
//   Output the orbitals using an instance of QdomainBasisFunction
//##################################################################

	UCLAQ::QdomainBasisFunctions qDomainBasisFunction;

	printf("Writing QdomainBasisFunctions Instance To : %s \n",basisFileName.c_str());

    qDomainBasisFunction.outputToFile(eigValues, eigVectors, basisFileName);


//##################################################################
//   Output the smoothed nuclear potential
//##################################################################

    UCLAQ::GridFunction3dUtility gUtility3d;
    gUtility3d.outputToBinaryDataFile(sOp.potential, potentialFilename);

   // Print out data

    printf("\n");
    printf("H2 Nuclear Potential Output To : %s \n\n",potentialFilename.c_str());
    printf(" Potential Grid Structure \n\n");
    printf(" [xPanels, xMin, xMax] : %ld, %10.5f %10.5f \n",sOp.potential.xPanels,sOp.potential.xMin, sOp.potential.xMax);
    printf(" [yPanels, yMin, yMax] : %ld, %10.5f %10.5f \n",sOp.potential.yPanels,sOp.potential.yMin, sOp.potential.yMax);
    printf(" [zPanels, zMin, zMax] : %ld, %10.5f %10.5f \n",sOp.potential.zPanels,sOp.potential.zMin, sOp.potential.zMax);
    printf("\n");
    printf(" Potential Norm_Inf: %20.15e   Norm_L2: %20.15e  \n\n",sOp.potential.normInf(),sOp.potential.norm2());

	printf("XXXX Execution Complete XXXXX\n");

}

