#include <iostream>
#include <string>
#include <cstdio>
#include <functional>

#include "InterpolationNd/LegendreInterp1d.h"
#include "InterpolationNd/FixedLocationLegendreInterp1d.h"

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"

using namespace std;

#undef  _OUTPUT_PLOTS    // #define to output 1d data for gnuplot

int main(int argc, char* argv[])
{
    int interpolationOrder = 4;
	double xMin = -5.0;
	double xMax  = 5.0;
    long xPanels = 50;

    //#####################################################
    // Initialize a UCLAQ::GridFunction1d with values of
    // a Gaussian centered at xBar and variance sigma^2.
    //#####################################################


	double xBar    = -0.2;
	double sigma   =  1.0;

    std::function < double(double) > gaussianFunction1d = [xBar,sigma](double x)
    {
    return exp(-((x-xBar)/sigma)*((x-xBar)/sigma));
    };

	UCLAQ::GridFunction1d gaussian1d(xPanels,xMin,xMax);
	gaussian1d.setToValue(0.0);

	gaussian1d.specify(gaussianFunction1d);

    //#####################################################
    // Initialize LegendreInterp1d instance and test
	// the evaluation of the interpolant at xPos
    //#####################################################

    double xPos           =  1.13;


    LegendreInterp1d interp1d(interpolationOrder);

    double interp = interp1d.evaluateInterpolant(xPos,gaussian1d);
    double exact  = gaussianFunction1d(xPos);

    printf("############################################\n");
    printf("###         LegendreInterp1d             ###\n");
    printf("###        Interpolation error           ###\n");
    printf("### evaluateInterpolant(xPos,gaussian1d) ###\n\n");
    printf("Interpolation order : %d \n",interpolationOrder);
    printf("Evaluation point    : %10.5e \n\n",xPos);

    printf("Exact        : %10.5e \n",exact);
    printf("Interpolated : %10.5e \n\n",interp);
    printf("Error (Abs)  : %10.5e \n",abs(exact-interp));
    printf("Error (Rel)  : %10.5e \n",abs(exact-interp)/abs(exact));
    printf("\n\n\n");


    //#####################################################
    // Initialize FixedLocationLegendreInterp1d instance and test
	// the evaluation of the interpolant at xPos
    //#####################################################

    FixedLocationLegendreInterp1d fixedInterp(interpolationOrder,xPos,gaussian1d);

    // Evaluate the interpolant using data from the grid function

    interp = fixedInterp.getInterpolatedValue(gaussian1d);

    printf("#################################################\n");
    printf("###    FixedLocationLegendreInterp1d          ###\n");
    printf("###        Interpolation error                ###\n");
    printf("###    getInterpolatedValue(gaussian1d)       ###\n\n");
    printf("Interpolation order : %d \n",interpolationOrder);
    printf("Evaluation point    : %10.5e \n\n",xPos);

    printf("Exact        : %10.5e \n",exact);
    printf("Interpolated : %10.5e \n\n",interp);
    printf("Error (Abs)  : %10.5e \n",abs(exact-interp));
    printf("Error (Rel)  : %10.5e \n",abs(exact-interp)/abs(exact));
    printf("\n\n");


    //#####################################################
    // Encapsulate the LegendreInterp1d instance inside a
    // lambda function and use it to initialize the values
    // of a UCLAQ::GridFunction1d that has mesh size 1/2
    // the mesh size of the grid providing the interpolant
    // data.
    //#####################################################

    std::function < double(double) > interpFunction1d = [&interp1d, &gaussian1d](double x)
    {
    return interp1d.evaluateInterpolant(x,gaussian1d);
    };

    // Create a grid function with 1/2 the mesh size of the grid
    // used to provide the values being interpolated.


	UCLAQ::GridFunction1d refinedInterpolated1d(2*xPanels,xMin,xMax);
	refinedInterpolated1d.setToValue(0.0);

	UCLAQ::GridFunction1d refinedExact1d(2*xPanels,xMin,xMax);
	refinedExact1d.setToValue(0.0);

	UCLAQ::GridFunction1d refinedError(2*xPanels,xMin,xMax);
	refinedError.setToValue(0.0);

	refinedInterpolated1d.specify(interpFunction1d);
    refinedExact1d.specify(gaussianFunction1d);

    refinedError = refinedExact1d - refinedInterpolated1d;

    printf("############################################\n");
    printf("###         LegendreInterp1d             ###\n");
    printf("### Interpolation error of evaluation    ###\n");
    printf("###   at nodes of a 2X refined grid      ###\n\n");
    printf("Interpolation order : %d \n\n",interpolationOrder);

    printf("Error (Inf)  : %10.5e \n",refinedError.normInf());
    printf("Error (L2)   : %10.5e \n",refinedError.norm2());
    printf("\n\n\n");

#ifdef _OUTPUT_PLOTS
    UCLAQ::GridFunction1dUtility gUtility1d;
    gUtility1d.outputToGNUplot(refinedInterpolated1d,"interp1d.dat");
    gUtility1d.outputToGNUplot(refinedExact1d,"exact1d.dat");
    gUtility1d.outputToGNUplot(refinedError,"error1d.dat");
#endif

    printf("XXXX Execution Complete XXXXX\n");



}

