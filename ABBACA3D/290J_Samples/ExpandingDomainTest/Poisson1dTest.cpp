#include <iostream>
#include <string>
#include <cstdio>
#include <functional>

#include "ExpandingDomainNd/InversePoisson1d.h"
#include "MollifierNd/SmoothPolyMollifier1d.h"
#include "MollifierNd/SmoothPolyPotential1d.h"

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"

using namespace std;

// Use  #define _OUTPUT_PLOTS to output 1d data for gnuplot

#undef  _OUTPUT_PLOTS

int main(int argc, char* argv[])
{

//################################################################
//   Problem specification
//################################################################

//  In one dimension, the convergence of the approximation only depends
//  on the mesh size and the smoothness of the right hand side;
//  not the domain size.

//  Test case parameters typically varied

    double hx = 0.1;                      // Grid mesh size
    int    sourceDifferentiability =   6; // Range 0 to 9

//  Test case static parameters

    double laplaceCoeff = -1.0;           // Coefficient of the Laplace operator

	double xMin = -2.0;                   // Domain size
	double xMax  = 2.0;

	double xBar     = -0.2;               // Parameters for SmoothPolyMollifier
	double sigma    =  0.7;               // used as the right hand side
	double strength =  1.0;

	long xPanels = (xMax-xMin)/hx;
	hx           = (xMax-xMin)/(double)xPanels;

//  Instantiate an exact potential due to a SmooothPolyMollifier source
//  with a specified differentiability

	SmoothPolyPotential1d potentialFunction1d(xBar,sigma, strength, laplaceCoeff);

	potentialFunction1d.setSourceDifferentiability(sourceDifferentiability);

//  Set up grid functions

	UCLAQ::GridFunction1d source1d(xPanels,xMin,xMax);
	source1d.setToValue(0.0);

	UCLAQ::GridFunction1d potential1d(xPanels,xMin,xMax);
	potential1d.setToValue(0.0);

    UCLAQ::GridFunction1d exactPotential1d(xPanels,xMin,xMax);
	exactPotential1d.setToValue(0.0);

	UCLAQ::GridFunction1d error1d(xPanels,xMin,xMax);
	error1d.setToValue(0.0);

// Capture the exact potential and source values

	exactPotential1d.specify(potentialFunction1d.getEvaluationPtr());

	source1d.specify(potentialFunction1d.getSourceEvaluationPtr());

//################################################################
//  Instantiate and initialize infinite domain Poisson solver
//################################################################

    InversePoisson1d inversePoissonOp(laplaceCoeff,xPanels,xMin,xMax);

//################################################################
//  Solve the test case problem and evaluate the error
//################################################################

//  Compute solution

    potential1d = source1d;
    inversePoissonOp.applyInverseOp(potential1d);

//  Evaluate the error and print out the results

    error1d = potential1d - exactPotential1d;

    printf("############################################\n");
    printf("###         InversePoisson1d             ###\n");
    printf("### Infinite domain potential error      ###\n\n");

    printf("HX                      : %-10.5f \n",(xMax-xMin)/(double)xPanels);
    printf("Source differentiabiliy : %d \n\n",sourceDifferentiability);

    printf("Domain Size             : %10.5f\n\n", xMax-xMin);

    printf("Error (Inf)  : %10.5e \n",error1d.normInf());
    printf("Error (L2)   : %10.5e \n",error1d.norm2());
    printf("\n\n\n");

#ifdef _OUTPUT_PLOTS
    UCLAQ::GridFunction1dUtility gUtility1d;
    gUtility1d.outputToGNUplot(potential1d,"computed1d.dat");
    gUtility1d.outputToGNUplot(exactPotential1d,"exact1d.dat");
    gUtility1d.outputToGNUplot(error1d,"error1d.dat");
#endif

//################################################################
//  Solve the test case problem using a mesh size 1/2 the size
//################################################################

// Re-initialize grid functions

    xPanels *= 2;
    hx = (xMax-xMin)/(double)(xPanels);

	source1d.initialize(xPanels,xMin,xMax);
	potential1d.initialize(xPanels,xMin,xMax);
    exactPotential1d.initialize(xPanels,xMin,xMax);
	error1d.initialize(xPanels,xMin,xMax);

// Capture the exact potential and source values

	exactPotential1d.specify(potentialFunction1d.getEvaluationPtr());

	source1d.specify(potentialFunction1d.getSourceEvaluationPtr());

//#####################################################
// Re-initialize infinite domain Poisson solver
//#####################################################

    inversePoissonOp.initialize(laplaceCoeff,xPanels,xMin,xMax);

//################################################################
//  Solve the test case problem and evaluate the error
//################################################################

    potential1d = source1d;
    inversePoissonOp.applyInverseOp(potential1d);

    error1d = potential1d - exactPotential1d;

    printf("############################################\n");
    printf("###         InversePoisson1d             ###\n");
    printf("### Infinite domain potential error      ###\n");
    printf("###   using 1/2 previous mesh size       ###\n\n");

    printf("HX                      : %-10.5f \n",(xMax-xMin)/(double)xPanels);
    printf("Source differentiabiliy : %d \n\n",sourceDifferentiability);

    printf("Domain Size             : %10.5f\n\n", (xMax-xMin));

    printf("Error (Inf)  : %10.5e \n",error1d.normInf());
    printf("Error (L2)   : %10.5e \n",error1d.norm2());
    printf("\n\n\n");

    printf("XXXX Execution Complete XXXXX\n");

}

