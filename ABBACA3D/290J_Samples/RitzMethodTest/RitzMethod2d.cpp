#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <algorithm>
using namespace std;
//
// Ritzmethod2d.cpp
//
// A test code that demonstrates the use of a Ritz procedure to
// create approximate eigenfunctions.
//
// The operator whose eigenvalues and eigenvectors are being computed
// is the standard five-point discretization of
//
//  [ alpha*DELTA  + beta*(x^2 + y^2) ]  
//
// where the Laplace operator DELTA is specified with homogeneous boundary conditions.
//
// The basis used is a product of sin functions.
//
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
// g++ RitzMethod2d.cpp -O2 -std=c++11 -I../../290J_2015 -I../LaplaceOp2dTest -o RitzMethod2d.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/RitzmethodTest is
//
// make -f RitzMethod2d.mk release
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Nov. 26, 2015
//

#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "LaplaceOp2d.h"

#include "RandOpNd/RandOp1d.h"
#include "RayleighChebyshev/RayleighChebyshev.h"


//##################################################################
//   ParabolicSchroedingerOp2d:
//   A class whose apply member function implements the operator
//
//   alpha*Delta_h + beta*(x^2 + y^2)
//##################################################################

class ParabolicSchroedingerOp2d
{

public:

	ParabolicSchroedingerOp2d(double alpha, double beta) : laplaceOp2d(alpha)
	{
		this->alpha = alpha;
		this->beta  = beta;

		parabolicPotential = [this](double x, double y) 
	    {return this->beta*(x*x + y*y);};
	}

	void initialize(double alpha, double beta)
	{
		this->alpha = alpha;
		this->beta  = beta;
		laplaceOp2d.initialize(alpha);

	    parabolicPotential = [this](double x, double y)  
	    {return this->beta*(x*x + y*y);};
	}

	void apply(UCLAQ::GridFunction2d& V)
	{
		Vtmp.initialize(V);
		Vtmp *= parabolicPotential;

		laplaceOp2d.apply(V);
		V += Vtmp;
	}

	UCLAQ::GridFunction2d      Vtmp;
	UCLAQ::GridFunction2d potential;
	LaplaceOp2d         laplaceOp2d;  // Laplace operator with homogeneous Dirichlet boundary conditions
	double alpha;                     // coefficient of the laplace operator
	double beta;                      // coefficient of the parabolic potential

    std::function< double(double, double) > parabolicPotential;
};

//##################################################################
//   Class used to define a Matrix operator
//   that can be used with the Rayleigh-Chebyshve procedure
//##################################################################

class MatrixOp
{
	public: MatrixOp(const UCLAQ::DoubleVector2d& M)
	{
		Mptr = &M;
	}
	void apply(UCLAQ::DoubleVector1d& v)
	{
		vTmp.initialize(v);
		for(long i = 0; i < Mptr->getIndex1Size(); i++)
		{
		v(i) = 0.0;
		for(long j = 0; j < Mptr->getIndex2Size(); j++)
		{
			v(i) += Mptr->operator()(i,j)*vTmp(j);
		}}
	}
	const UCLAQ::DoubleVector2d* Mptr;
	UCLAQ::DoubleVector1d  vTmp;
};

int main()
{
	// Maximal sin mode in each direction

	long kxMax = 20;
	long kyMax = 20;

	//   Problem set up.

	long xPanels = 40; long yPanels = 40;

	double xMin = -2.0; double xMax = 2.0;
	double yMin = -2.0; double yMax = 2.0;

	double LX   = (xMax-xMin);
	double LY   = (yMax-yMin);

	double pi   =  3.141592653589793238;

//##################################################################
//               Creating basis set
//##################################################################

	// Use a lambda function to specify the sin basis functions.

	// !!! Note the use of &kx and &ky so that when this function
	// is used in the initialization loop, the values of kx
	// and ky are implicitly passed into the lambda function.

	double kx; double ky;

	std::function< double(double, double) > sinXsinY
    = [&kx,  &ky, xMin, yMin, LX, LY, pi] (double x, double y)
	{
	return (2.0/sqrt(LX*LY))*sin(kx*pi*((x-xMin)/LX))*sin(ky*pi*((y-yMin)/LY));
	};

	// Create a vector of grid functions for the basis functions

	UCLAQ::GridFunction2d  uTemp(xPanels,xMin,xMax,yPanels,yMin,yMax);
	vector <UCLAQ::GridFunction2d > basisFunctions(kxMax*kyMax,uTemp);

	// Initialize values of each element of the basis set

	long basisIndex = 0;
	for(kx = 1; kx <= kxMax; kx++)
	{
	for(ky = 1; ky <= kyMax; ky++)
	{
    basisFunctions[basisIndex].specify(sinXsinY);
    basisIndex++;
	}}

	// Verify orthonormality of the basis functions using scaled dot product.

	double dotProd;
	double dotProdError      = 0.0;
    double dotProdErrorMax   = 0.0;

	for(long i = 0; i < (long)basisFunctions.size(); i++)
	{
    for(long j = 0; j < (long)basisFunctions.size(); j++)
	{
		dotProd = basisFunctions[i].scaledDot(basisFunctions[j]);
		if(i == j) dotProdError = abs(dotProd -1.0);
		else       dotProdError = abs(dotProd);

	    dotProdErrorMax = (dotProdErrorMax < dotProdError) ? dotProdError : dotProdErrorMax;
	}}

	printf(" Basis set orthonomality error : %10.5e \n",dotProdErrorMax);

//##################################################################
//               Initializing operator
//##################################################################

	double alpha = -0.5;
	double beta  =  5.0;

	ParabolicSchroedingerOp2d pSop(alpha, beta);

//##################################################################
//            Construct entries of Ritz matrix
//         Use symmetry to reduce computational work
//##################################################################

	UCLAQ::DoubleVector2d Hbar(basisFunctions.size(),basisFunctions.size());

	for(long i = 0; i < (long)basisFunctions.size(); i++)
	{
		uTemp = basisFunctions[i];
		pSop.apply(uTemp);

		for(long j = i; j < (long)basisFunctions.size(); j++)
		{
		dotProd = basisFunctions[j].scaledDot(uTemp);
		Hbar(i,j) = dotProd;
		Hbar(j,i) = dotProd;
		}
	}

//##################################################################
//   Instantiate  RayleighChebyshev instance and compute
//   smallest 10 eigenvalues of the Ritz matrix
//##################################################################

	long ritzSystemSize = basisFunctions.size();

	// Set up eigensystem solution

	UCLAQ::RandOp1d  randomOp;

	MatrixOp ritzOp(Hbar);

    // Allocate arrays for eigenvectors and eigenvalues

    vector <UCLAQ::DoubleVector1d>  eigVectors;
    vector <double>                  eigValues;

    // Declare an instance of the Raylegh-Chebyshev eigensystem procedure

    RayleighChebyshev < UCLAQ::DoubleVector1d, MatrixOp, UCLAQ::RandOp1d > RCeigProcedure;

    RCeigProcedure.setEigDiagnosticsFlag();
    RCeigProcedure.setVerboseFlag();

    UCLAQ::DoubleVector1d vTmp(ritzSystemSize);     // A temporary vector is required as input. This vector must
                                                    // be a non-null instance of the vector class

    long  dimension            = vTmp.getDimension();
    double subspaceTol         = 2.0e-6;
    long subspaceIncrementSize = 3;
    long bufferSize            = 3;
    long eigCount              = dimension < 10 ? dimension : 10;

    RCeigProcedure.getMinEigenSystem(eigCount, subspaceTol, subspaceIncrementSize, bufferSize, vTmp,
    ritzOp, randomOp, eigValues, eigVectors);

//##################################################################
//   Evaluate the exact eigenvalues and compare to the computed
//   results.
//##################################################################

    printf("\n\nXXXX   RC_OperatorEig_Test Results XXXX\n\n");
    printf("Tolerance Specified : %10.5e\n\n",subspaceTol);

    vector<double> exactEigValues(100,0.0);

    long index = 0;
	for(long k1 = 0; k1 < 10; k1++)
	{
	for(long k2 = 0; k2 < 10; k2++)
	{
	exactEigValues[index] = sqrt(beta*abs(alpha))*(2.0*k1 + 2.0*k2 + 2.0);
    index++;
	}}

    // Create sorted list of algebraically smallest to algebraically largest

    sort(exactEigValues.begin(),exactEigValues.end());

    printf("       Eigenvalue     Exact   Error     Relative Error \n");
    for(long  k = 0; k < eigCount; k++ )
    {
    	printf("%-5ld %-10.5e  %-10.5e  %10.5e   %-10.5e\n", k+1, eigValues[k],  exactEigValues[k], abs(eigValues[k] - exactEigValues[k]),abs(eigValues[k] -exactEigValues[k])/abs(exactEigValues[k]));
    }

}

