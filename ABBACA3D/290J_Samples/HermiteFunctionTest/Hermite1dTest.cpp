#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <functional>
using namespace std;

//
// Hermite1dTest.cpp : a test program for product Hermite functions in 1 dimension
//
// Each Hermite function is parameterized by a real value gamma
// and non-negative integer values.

// For real values of beta and alpha of opposite sign, if
// gamma = ||beta/alpha||^(1/4) then the  Hermite functions
// for any triplet of non-negative integer values (m,n) is an eigenfunction of
//
// L =  alpha*DELTA + beta*(x^2)
//
// Where DELTA is the Laplace operator.
//
// In this test program:
//
// (I) A vector of UCLAQ::GridFunction1d instances is created with grid values given by
// a Hermite function associated with a value of gamma.

// (II) The orthonormality of this collection of grid functions is checked with respect to the
// mesh scaled inner product. As will be observed the orthonormality of the function is not
// especially accurate because of numerical approximation of the inner products and
// the use of a finite domain.
//
// (III) The vector of  grid functions is then orthonormalized using a
// modified Gram-Schmidt process and the orthonormality is checked to verify
// their discrete orthonomality with respect to the scaled dot product.
//

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"

#include "OrthoFunctionNd/UCLAQ_HermiteFunction.h"

void orthoNormalizeMGS1d(vector< UCLAQ::GridFunction1d > & V); // MGS routine (source after main(...)

int main(int argc, char* argv[])
{
	// Grid parameters - purposely choosing a coarse mesh to demonstrate
	// that the analytically orthonormal functions are not discretely
	// orthonormal.

	double xMin    = -5.0;
	double xMax    =  5.0;
    long   xPanelsCoarse  =   30;
    long   xPanelsFine    =  100;

	//
    // (I) Construction of an vector of UCLAQ::GridFunction3d's
    //     whose values are given by  Hermite project functions
    //

    //  Specify coefficients of linear operator that will be used to determine gamma

    double alpha = -0.5;
	double beta  =  5.0;

	double gamma    =  pow(abs(beta/alpha),0.25); // Coefficient in exponent of Hermite function
	double xShift   =  0.2;                       // Shift of origin of the Hermite function
	long indexBound =   11;                       // index of the Hermite functions < indexBound

	UCLAQ::HermiteFunction hermiteFunction;
	hermiteFunction.initialize(gamma,xShift);

    // Create an array UCLAQ::GridFunction1d with values given by
	// the hermite functions with index i = 0 ... indexBound-1

	vector< UCLAQ::GridFunction1d > hermiteGridFunArray(indexBound);

    for(long i = 0; i < indexBound; i++)
    {
    	hermiteGridFunArray[i].initialize(xPanelsCoarse,xMin,xMax);
    	hermiteGridFunArray[i].specify(hermiteFunction.getNthHermiteFunction(i));
    }


	//
    // (II) Checking orthonormality
    //


    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    printf("      1D Coarse Mesh Results \n");
    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n\n");

    printf("Inner products before discrete orthonormalization : \n\n");

    // Verify orthonormality of the basis functions using scaled dot product.

	double dotProd;
	double dotProdError      = 0.0;
    double dotProdErrorMax   = 0.0;

	for(long i = 0; i < (long)hermiteGridFunArray.size(); i++)
	{
    for(long j = 0; j < (long)hermiteGridFunArray.size(); j++)
	{
		dotProd = hermiteGridFunArray[i].scaledDot(hermiteGridFunArray[j]);
		if(i == j) dotProdError = abs(dotProd -1.0);
		else       dotProdError = abs(dotProd);

	    dotProdErrorMax = (dotProdErrorMax < dotProdError) ? dotProdError : dotProdErrorMax;

	    printf("%5.3e ",abs(dotProd));
	}
	printf("\n");
	}

	printf("\n");
	printf("Maximum orthonomality error : %10.5e \n\n",dotProdErrorMax);

	//
    // (III) Orthonormalizing using modified Gram-Schmidt and
	//       re-checking orthonormality.
    //

	orthoNormalizeMGS1d(hermiteGridFunArray);

	printf("Inner products after discrete orthonormalization : \n\n");

    // Verify orthonormality of the basis functions using scaled dot product.

	dotProdError      = 0.0;
    dotProdErrorMax   = 0.0;

	for(long i = 0; i < (long)hermiteGridFunArray.size(); i++)
	{
    for(long j = 0; j < (long)hermiteGridFunArray.size(); j++)
	{
		dotProd = hermiteGridFunArray[i].scaledDot(hermiteGridFunArray[j]);
		if(i == j) dotProdError = abs(dotProd -1.0);
		else       dotProdError = abs(dotProd);

	    dotProdErrorMax = (dotProdErrorMax < dotProdError) ? dotProdError : dotProdErrorMax;

	    printf("%5.3e ",abs(dotProd));
	}
	printf("\n");
	}

	printf("\n");
	printf("Maximum orthonomality error : %10.5e \n\n\n",dotProdErrorMax);

	// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	//        Repeat the test with a fine grid
    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Create an array UCLAQ::GridFunction1d with values given by
	// the hermite functions with index i = 0 ... indexBound-1


    for(long i = 0; i < indexBound; i++)
    {
    	hermiteGridFunArray[i].initialize(xPanelsFine,xMin,xMax);
    	hermiteGridFunArray[i].specify(hermiteFunction.getNthHermiteFunction(i));
    }

    // Check orthonormality

    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    printf("      1D Fine Mesh Results \n");
    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n\n");

    printf("Inner products before discrete orthonormalization : \n\n");

    // Verify orthonormality of the basis functions using scaled dot product.

	dotProdError      = 0.0;
    dotProdErrorMax   = 0.0;

	for(long i = 0; i < (long)hermiteGridFunArray.size(); i++)
	{
    for(long j = 0; j < (long)hermiteGridFunArray.size(); j++)
	{
		dotProd = hermiteGridFunArray[i].scaledDot(hermiteGridFunArray[j]);
		if(i == j) dotProdError = abs(dotProd -1.0);
		else       dotProdError = abs(dotProd);

	    dotProdErrorMax = (dotProdErrorMax < dotProdError) ? dotProdError : dotProdErrorMax;

	    printf("%5.3e ",abs(dotProd));
	}
	printf("\n");
	}

	printf("\n");
	printf("Maximum orthonomality error : %10.5e \n\n",dotProdErrorMax);


	// Orthonormalize using modified Gram-Schmidt

	 orthoNormalizeMGS1d(hermiteGridFunArray);

	printf("Inner products after discrete orthonormalization : \n\n");

    // Verify orthonormality of the basis functions using scaled dot product.

	dotProdError      = 0.0;
    dotProdErrorMax   = 0.0;

	for(long i = 0; i < (long)hermiteGridFunArray.size(); i++)
	{
    for(long j = 0; j < (long)hermiteGridFunArray.size(); j++)
	{
		dotProd = hermiteGridFunArray[i].scaledDot(hermiteGridFunArray[j]);
		if(i == j) dotProdError = abs(dotProd -1.0);
		else       dotProdError = abs(dotProd);

	    dotProdErrorMax = (dotProdErrorMax < dotProdError) ? dotProdError : dotProdErrorMax;

	    printf("%5.3e ",abs(dotProd));
	}
	printf("\n");
	}

	printf("\n");
	printf("Maximum orthonomality error : %10.5e \n\n",dotProdErrorMax);

    printf("XXXX Execution Complete XXXXX\n");

}

//  Orthonormalize the vectors using Modified Gram-Schmidt

void  orthoNormalizeMGS1d(vector< UCLAQ::GridFunction1d > & V)
{
	long subspaceSize = V.size();
	long j;
	long k;
	double rkk;
	double rkj;

    for(k = 0; k < subspaceSize; k++)
    {
        rkk  = sqrt(abs(V[k].scaledDot(V[k])));
        V[k].scal(1.0/rkk);
        for(j = k+1; j < subspaceSize; j++)
        {
            rkj  =   V[j].scaledDot(V[k]);
            V[j].axpy(-rkj,V[k]);
        }
    }
}


