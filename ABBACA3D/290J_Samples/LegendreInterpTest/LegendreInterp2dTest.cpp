#include <iostream>
#include <string>
#include <cstdio>
#include <functional>

#include "InterpolationNd/LegendreInterp2d.h"
#include "InterpolationNd/FixedLocationLegendreInterp2d.h"

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"

using namespace std;

#undef  _OUTPUT_PLOTS // #define to output 2d data for gnuplot

int main(int argc, char* argv[])
{
    int interpolationOrder = 4;

	double xMin = -5.0;
	double xMax  = 5.0;

	double yMin = -3.0;
	double yMax =  2.0;

    long xPanels = 50;
    long yPanels = 45;

    //##################################################################
    // Initialize a UCLAQ::GridFunction2d with values of a
    // Gaussian centered at (xBar,yBar) and variance sigma^2.
    //##################################################################


	double xBar    = -0.2; double yBar = 0.1;
	double sigma   =  1.0;

    std::function < double(double, double) > gaussianFunction2d = [xBar,yBar,sigma](double x,double y)
    {
    return exp(-( (x-xBar)*(x-xBar) + (y-yBar)*(y-yBar) )/(sigma*sigma));
    };

	UCLAQ::GridFunction2d gaussian2d(xPanels,xMin,xMax,yPanels,yMin,yMax);
	gaussian2d.setToValue(0.0);

	gaussian2d.specify(gaussianFunction2d);

    //#####################################################
    // Initialize LegendreInterp2d instance and test
	// the evaluation of the interpolant at (xPos,yPos)
    //#####################################################

    double xPos           =  1.13;
    double yPos           =  0.13;


    LegendreInterp2d interp2d(interpolationOrder);

    double interp = interp2d.evaluateInterpolant(xPos,yPos,gaussian2d);
    double exact  = gaussianFunction2d(xPos,yPos);

    printf("#################################################\n");
    printf("###             LegendreInterp2d              ###\n");
    printf("###            Interpolation error            ###\n");
    printf("### evaluateInterpolant(xPos,yPos,gaussian2d) ###\n\n");
    printf("Interpolation order : %d \n",interpolationOrder);
    printf("Evaluation point    : (%10.5e, %10.5E) \n\n",xPos,yPos);

    printf("Exact        : %10.5e \n",exact);
    printf("Interpolated : %10.5e \n\n",interp);
    printf("Error (Abs)  : %10.5e \n",abs(exact-interp));
    printf("Error (Rel)  : %10.5e \n",abs(exact-interp)/abs(exact));
    printf("\n\n\n");



    //#####################################################
    // Initialize FixedLocationLegendreInterp2d instance and test
	// the evaluation of the interpolant at (xPos,yPos)
    //#####################################################

    FixedLocationLegendreInterp2d fixedInterp(interpolationOrder,xPos,yPos, gaussian2d);

    // Evaluate the interpolant using data from the grid function

    interp = fixedInterp.getInterpolatedValue(gaussian2d);

    printf("#################################################\n");
    printf("###    FixedLocationLegendreInterp2d          ###\n");
    printf("###        Interpolation error                ###\n");
    printf("###    getInterpolatedValue(gaussian2d)       ###\n\n");
    printf("Interpolation order : %d \n",interpolationOrder);
    printf("Evaluation point    : (%10.5e, %10.5E) \n\n",xPos,yPos);

    printf("Exact        : %10.5e \n",exact);
    printf("Interpolated : %10.5e \n\n",interp);
    printf("Error (Abs)  : %10.5e \n",abs(exact-interp));
    printf("Error (Rel)  : %10.5e \n",abs(exact-interp)/abs(exact));
    printf("\n\n");


    //#####################################################
    // Encapsulate the LegendreInterp2d instance inside a
    // lambda function and use it to initialize the values
    // of a UCLAQ::GridFunction2d that has mesh size 1/2
    // the mesh size of the grid providing the interpolant
    // data.
    //#####################################################

    std::function < double(double,double) > interpFunction2d = [&interp2d, &gaussian2d](double x,double y)
    {
    return interp2d.evaluateInterpolant(x,y,gaussian2d);
    };

    // Create a grid function with 1/2 the mesh size of the grid
    // used to provide the values being interpolated.


	UCLAQ::GridFunction2d refinedInterpolated2d(2*xPanels,xMin,xMax,2*yPanels,yMin,yMax);
	refinedInterpolated2d.setToValue(0.0);

	UCLAQ::GridFunction2d refinedExact2d(2*xPanels,xMin,xMax,2*yPanels,yMin,yMax);
	refinedExact2d.setToValue(0.0);

	UCLAQ::GridFunction2d refinedError(2*xPanels,xMin,xMax,2*yPanels,yMin,yMax);
	refinedError.setToValue(0.0);

	refinedInterpolated2d.specify(interpFunction2d);
    refinedExact2d.specify(gaussianFunction2d);

    refinedError = refinedExact2d - refinedInterpolated2d;

    printf("############################################\n");
    printf("###         LegendreInterp2d             ###\n");
    printf("### Interpolation error of evaluation    ###\n");
    printf("###   at nodes of a 2X refined grid      ###\n\n");
    printf("Interpolation order : %d \n\n",interpolationOrder);

    printf("Error (Inf)  : %10.5e \n",refinedError.normInf());
    printf("Error (L2)   : %10.5e \n",refinedError.norm2());
    printf("\n\n\n");

#ifdef _OUTPUT_PLOTS
    UCLAQ::GridFunction2dUtility gUtility2d;
    gUtility2d.outputToGNUplot(refinedInterpolated2d,"interp2d.dat");
    gUtility2d.outputToGNUplot(refinedExact2d,"exact2d.dat");
    gUtility2d.outputToGNUplot(refinedError,"error2d.dat");
#endif

    printf("XXXX Execution Complete XXXXX\n");



}

