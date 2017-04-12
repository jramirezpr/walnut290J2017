#include <iostream>
#include <string>
#include <cstdio>
#include <functional>
using namespace std;


//
// Hermite2dTest.cpp : a test program for product Hermite functions in 2 dimensions
//
// Each term in a product Hermite function is parameterized by a real value gamma
// and non-negative integer values.

// For real values of beta and alpha of opposite sign, if
// gamma = ||beta/alpha||^(1/4) then the product Hermite functions
// for any triplet of non-negative integer values (m,n) is an eigenfunction of
//
// L =  alpha*DELTA + beta*(x^2 + y^2)
//
// Where DELTA is the Laplace operator.
//
// In this test program:
//
// (I) A vector of UCLAQ::GridFunction2d instances is created with grid values given by
// a Hermite product functions associated with a value of gamma.

// (II) The orthonormality of this collection of grid functions is checked with respect to the
// mesh scaled inner product. As will be observed the orthonormality of the function is not
// especially accurate because of numerical approximation of the inner products and
// the use of a finite domain.
//
// (III) The vector of  grid functions is then orthonormalized using a
// modified Gram-Schmidt process and the orthonormality is checked to verify
// their discrete orthonomality with respect to the scaled dot product.
//

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"

#include "OrthoFunctionNd/UCLAQ_HermiteProductFunction2d.h"

void orthoNormalizeMGS2d(vector< UCLAQ::GridFunction2d > & V); // MGS routine (source after main(...)


#define  _OUTPUT_PLOTS // #define to output 2d data 

int main(int argc, char* argv[])
{
	// Grid parameters - purposely choosing a coarse mesh to demonstrate
	// that the analytically orthonormal functions are not discretely
	// orthonormal.

	double xMin = -5.0;
	double xMax  = 5.0;

	double yMin = -5.0;
	double yMax =  5.0;

    long xPanelsCoarse = 30;
    long yPanelsCoarse = 30;

    long xPanelsFine   = 100;
    long yPanelsFine   = 100;

	//
    // (I) Construction of an vector of UCLAQ::GridFunction3d's
    //     whose values are given by  Hermite project functions
    //

    //  Specify coefficients of linear operator that will be used to determine gamma

	double alpha = -0.5;
	double beta  =  5.0;

	double gamma   =   pow(abs(beta/alpha),0.25); // Coefficient in exponent of Hermite function
	long indexBound =  3;                         // index of x and y Hermite functions in product functions < indexBound


	double gammaX = gamma;   // Using equal exponential coefficient in each direction
	double gammaY = gamma;

	double xShift = -0.1;    // Shift of origin of the Hermite function
	double yShift =  0.2;

	UCLAQ::HermiteProductFunction2d H;
	H.initialize(gammaX,gammaY,xShift,yShift);

	vector < UCLAQ::GridFunction2d > hermiteGridFunArray(indexBound*indexBound);

	long functionIndex = 0;
	for(long n = 0; n < indexBound; n++)
	{
	for(long m = 0; m < indexBound; m++)
	{
	hermiteGridFunArray[functionIndex].initialize(xPanelsCoarse,xMin,xMax,yPanelsCoarse,yMin,yMax);
    hermiteGridFunArray[functionIndex].specify(H.getHermiteProductFunction2d(m,n));
    functionIndex++;
	}}


	//
    // (II) Checking orthonormality
    //

    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    printf("      2D Coarse Mesh Results \n");
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


	orthoNormalizeMGS2d(hermiteGridFunArray);

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

	functionIndex = 0;
	for(long n = 0; n < indexBound; n++)
	{
	for(long m = 0; m < indexBound; m++)
	{
	hermiteGridFunArray[functionIndex].initialize(xPanelsFine,xMin,xMax,yPanelsFine,yMin,yMax);
    hermiteGridFunArray[functionIndex].specify(H.getHermiteProductFunction2d(m,n));
    functionIndex++;
	}}

    // Check orthonormality

    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    printf("      2D Fine Mesh Results \n");
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

	orthoNormalizeMGS2d(hermiteGridFunArray);

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

void  orthoNormalizeMGS2d(vector< UCLAQ::GridFunction2d > & V)
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

