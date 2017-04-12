#include <iostream>
#include <string>
#include <cstdio>
#include <functional>

#include "ExpandingDomainNd/InversePoisson2d.h"
#include "MollifierNd/SmoothPolyMollifier2d.h"
#include "MollifierNd/SmoothPolyPotential2d.h"

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"

using namespace std;

#undef  _OUTPUT_PLOTS // #define to output 2d data for gnuplot

int main(int argc, char* argv[])
{
//
// The extension factor is the factor by which the
// x-coordinate domain (1st coordinate) is implicitly
// expanded in the InversePoisson2d instance.
//
// Convergence to the solution is obtained as a double
// limit of the mesh size -> 0 and extension factor -> oo
// Typically one converges solutions w.r.t. the mesh
// size first and then increases the extension factor.
// Specifically, on fixes the extension factor and then
// determines a mesh size so that relative changes
// induced by a reduction in the mesh size are acceptably
// small. One then fixes the mesh size and increases the
// extension factor until relative changes induced by the
// increase are acceptably small.
//
//

    double extensionFactor = 2.0;

	double xMin = -1.0;
	double xMax  = 1.0;

	double yMin = -1.0;
	double yMax =  1.0;

    long xPanels = 100;
    long yPanels = 100;

	double xBar     = -0.1;
	double yBar     =  0.2;
	double sigma    =  0.5;
	double strength =  1.0;

    double laplaceCoeff            =   1.0;
	int    sourceDifferentiability =   8; // Source differentiability range = 0 to 9

	// Instantiate an exact potential with a specified source differentiability

	SmoothPolyPotential2d potentialFunction2d(xBar,yBar, sigma, strength, laplaceCoeff);
	potentialFunction2d.setSourceDifferentiability(sourceDifferentiability);

	UCLAQ::GridFunction2d source2d(xPanels,xMin,xMax,yPanels,yMin,yMax);
	source2d.setToValue(0.0);

	UCLAQ::GridFunction2d potential2d(xPanels,xMin,xMax,yPanels,yMin,yMax);
	potential2d.setToValue(0.0);

    UCLAQ::GridFunction2d exactPotential2d(xPanels,xMin,xMax,yPanels,yMin,yMax);
	exactPotential2d.setToValue(0.0);

	UCLAQ::GridFunction2d error2d(xPanels,xMin,xMax,yPanels,yMin,yMax);
	error2d.setToValue(0.0);

	// Capture the exact potential and source values on GridFunction2d
	// instances

	exactPotential2d.specify(potentialFunction2d.getEvaluationPtr());
	source2d.specify(potentialFunction2d.getSourceEvaluationPtr());

    //#####################################################
    // Initialize InversePoisson2d instance
    //#####################################################

    InversePoisson2d inversePoissonOp(laplaceCoeff,xPanels,xMin,xMax,yPanels,yMin,yMax,extensionFactor);

    potential2d = source2d;

    inversePoissonOp.applyInverseOp(potential2d);
    error2d = potential2d - exactPotential2d;

    printf("############################################\n");
    printf("###         InversePoisson2d             ###\n");
    printf("### Infinite domain potential error      ###\n\n");

    printf("Source differentiabiliy : %d \n\n",sourceDifferentiability);

    printf("Comp. Domain Size       : %8.5f X %-8.5f  \n", xMax-xMin,yMax-yMin);
    printf("HX                      : %8.5f \n",(xMax-xMin)/(double)xPanels);
    printf("HY                      : %8.5f \n\n",(yMax-yMin)/(double)yPanels);
    printf("Extension Factor        : %8.5f  \n\n",extensionFactor);

    printf("Error (Inf)  : %10.5e \n",error2d.normInf());
    printf("Error (L2)   : %10.5e \n",error2d.norm2());
    printf("\n");
    printf("Soln norm (Inf)  : %10.5e \n",potential2d.normInf());
    printf("Soln norm (L2)   : %10.5e \n",potential2d.norm2());
    printf("\n\n\n");

#ifdef _OUTPUT_PLOTS
    UCLAQ::GridFunction2dUtility gUtility2d;
    gUtility2d.outputToGNUplot(potential2d,"computed2d.dat");
    gUtility2d.outputToGNUplot(exactPotential2d,"exact2d.dat");
    gUtility2d.outputToGNUplot(error2d,"error2d.dat");
#endif

//
//  Increase extension factor to improve accuracy
//
    extensionFactor *= 2.0;
    inversePoissonOp.initialize(laplaceCoeff,xPanels,xMin,xMax,yPanels,yMin,yMax,extensionFactor);

    potential2d = source2d;

    inversePoissonOp.applyInverseOp(potential2d);

    error2d = potential2d - exactPotential2d;

    printf("############################################\n");
    printf("###         InversePoisson2d             ###\n");
    printf("### Infinite domain potential error      ###\n");
    printf("### using 2 X ExtensionFactor            ###\n\n");

    printf("Source differentiabiliy : %d \n\n",sourceDifferentiability);

    printf("Comp. Domain Size       : %8.5f X %-8.5f  \n", xMax-xMin,yMax-yMin);
    printf("HX                      : %8.5f \n",(xMax-xMin)/(double)xPanels);
    printf("HY                      : %8.5f \n\n",(yMax-yMin)/(double)yPanels);
    printf("New Extension Factor    : %8.5f  \n\n",extensionFactor);

    printf("Error (Inf)  : %10.5e \n",error2d.normInf());
    printf("Error (L2)   : %10.5e \n",error2d.norm2());
    printf("\n");
    printf("Soln norm (Inf)  : %10.5e \n",potential2d.normInf());
    printf("Soln norm (L2)   : %10.5e \n",potential2d.norm2());
    printf("\n\n\n");


    printf("XXXX Execution Complete XXXXX\n");
}

