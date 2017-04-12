#include <iostream>
#include <string>
#include <cstdio>
#include <functional>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

#include "OrthoFunctionNd/UCLAQ_HermiteProductFunction3d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//
// Hermite3dTest.cpp : a test program for product Hermite functions in 3 dimensions
//
// Each term in a product Hermite function is parameterized by a real value gamma
// and non-negative integer values.

// For real values of beta and alpha of opposite sign, if
// gamma = ||beta/alpha||^(1/4) then the product Hermite functions
// for any triplet of non-negative integer values (m,n,p) is an eigenfunction of
//
// L =  alpha*DELTA + beta*(x^2 + y^2 + z^2)
//
// Where DELTA is the Laplace operator.
//
// In this test program:
//
// (I) A vector of UCLAQ::GridFunction3d instances is created with grid values given by
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

void  orthoNormalizeMGS3d(vector< UCLAQ::GridFunction3d > & V);

#define  _OUTPUT_PLOTS // #define to output 3d data for VTK

int main(int argc, char* argv[])
{

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


	// Grid parameters - purposely choosing a coarse mesh to demonstrate
	// that the analytically orthonormal functions are not discretely
	// orthonormal.

	double xMin = -5.0;
	double xMax  = 5.0;

	double yMin = -5.0;
	double yMax =  5.0;

	double zMin = -5.0;
	double zMax =  5.0;

	long xPanelsCoarse = 30;
    long yPanelsCoarse = 30;
    long zPanelsCoarse = 30;

    long xPanelsFine = 100;
    long yPanelsFine = 100;
    long zPanelsFine = 100;


    //
    // (I) Construction of an vector of UCLAQ::GridFunction3d instances
    //     whose values are given by a Hermite project functions
    //

    //
    //  Specify coefficients of linear operator that will be used to determine gamma

	double alpha = -0.5;
	double beta  =  5.0;

	double gamma    =   pow(abs(beta/alpha),0.25); // Coefficient in exponent of Hermite function
	long indexBound =   2;                         // index of x, y and z Hermite functions in product functions < indexBound

    double gammaX = gamma;  // Using equal exponential coefficient in each direction
	double gammaY = gamma;
    double gammaZ = gamma;

	double xShift = -0.1;  // Shift of origin of the Hermite function
	double yShift =  0.2;
	double zShift =  0.333;

	UCLAQ::HermiteProductFunction3d H;
	H.initialize(gammaX,gammaY,gammaZ,xShift,yShift,zShift);

	vector < UCLAQ::GridFunction3d > hermiteGridFunArray(indexBound*indexBound*indexBound);

	long functionIndex = 0;
	for(long n = 0; n < indexBound; n++)
	{
	for(long m = 0; m < indexBound; m++)
	{
	for(long p = 0; p < indexBound; p++)
	{
	hermiteGridFunArray[functionIndex].initialize(xPanelsCoarse,xMin,xMax,yPanelsCoarse,yMin,yMax,zPanelsCoarse,zMin,zMax);
    hermiteGridFunArray[functionIndex].specify(H.getHermiteProductFunction3d(m,n,p));
    functionIndex++;
	}}}


	//
    // (II) Checking orthonormality
    //

    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    printf("      3D Coarse Mesh Results \n");
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


	orthoNormalizeMGS3d(hermiteGridFunArray);

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
	for(long p = 0; p < indexBound; p++)
	{
	hermiteGridFunArray[functionIndex].initialize(xPanelsFine,xMin,xMax,yPanelsFine,yMin,yMax,zPanelsFine,zMin,zMax);
    hermiteGridFunArray[functionIndex].specify(H.getHermiteProductFunction3d(m,n,p));
    functionIndex++;
	}}}

    // Check orthonormality

    printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");
    printf("      3D Fine Mesh Results \n");
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

	orthoNormalizeMGS3d(hermiteGridFunArray);

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
//
//  If OpenMP is defined then multi-thread the inner loop
//

void  orthoNormalizeMGS3d(vector< UCLAQ::GridFunction3d > & V)
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

#ifdef _OPENMP
		#pragma omp parallel for \
		private(j,rkj) \
		schedule(static,1)
#endif
        for(j = k+1; j < subspaceSize; j++)
        {
            rkj  =   V[j].scaledDot(V[k]);
            V[j].axpy(-rkj,V[k]);
        }
    }
}
