#include <iostream>
#include <string>
#include <cstdio>
#include <functional>

#include "InterpolationNd/LegendreInterp3d.h"
#include "InterpolationNd/FixedLocationLegendreInterp3d.h"

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

using namespace std;

#undef  _OUTPUT_PLOTS // #define to output 3d data for VTK

int main(int argc, char* argv[])
{
    int interpolationOrder = 4;

	double xMin = -5.0;
	double xMax  = 5.0;

	double yMin = -3.0;
	double yMax =  2.0;

	double zMin = -2.0;
	double zMax =  5.0;

    long xPanels = 50;
    long yPanels = 45;
    long zPanels = 55;

    //##################################################################
    // Initialize a UCLAQ::GridFunction3d with values of a
    // Gaussian centered at (xBar,yBar,zBar) and variance sigma^2.
    //##################################################################


	double xBar    = -0.2; double yBar = 0.1; double zBar = 0.3;
	double sigma   =  1.0;

    std::function < double(double, double, double) > gaussianFunction3d = [xBar,yBar,zBar,sigma](double x,double y, double z)
    {
    return exp(-( (x-xBar)*(x-xBar) + (y-yBar)*(y-yBar) + (z-zBar)*(z-zBar) )/(sigma*sigma));
    };

	UCLAQ::GridFunction3d gaussian3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	gaussian3d.setToValue(0.0);

	gaussian3d.specify(gaussianFunction3d);

    //#####################################################
    // Initialize LegendreInterp3d instance and test
	// the evaluation of the interpolant at (xPos,yPos,zPos)
    //#####################################################

    double xPos           =  1.13;
    double yPos           =  0.13;
    double zPos           =  0.13;

    LegendreInterp3d interp3d(interpolationOrder);

    double interp = interp3d.evaluateInterpolant(xPos,yPos,zPos, gaussian3d);
    double exact  = gaussianFunction3d(xPos,yPos,zPos);

    printf("########################################################\n");
    printf("###                LegendreInterp3d                  ###\n");
    printf("###              Interpolation error                 ###\n");
    printf("###  evaluateInterpolant(xPos,yPos,zPos,gaussian3d)  ###\n\n");
    printf("Interpolation order : %d \n",interpolationOrder);
    printf("Evaluation point    : (%7.3e, %7.3e, %7.3e) \n\n",xPos,yPos,zPos);

    printf("Exact        : %10.5e \n",exact);
    printf("Interpolated : %10.5e \n\n",interp);
    printf("Error (Abs)  : %10.5e \n",abs(exact-interp));
    printf("Error (Rel)  : %10.5e \n",abs(exact-interp)/abs(exact));
    printf("\n\n\n");

    //#####################################################
    // Initialize FixedLocationLegendreInterp3d instance and test
	// the evaluation of the interpolant at (xPos,yPos,zPos)
    //#####################################################

    FixedLocationLegendreInterp3d fixedInterp(interpolationOrder,xPos,yPos,zPos,gaussian3d);

    // Evaluate the interpolant using data from the grid function

    interp = fixedInterp.getInterpolatedValue(gaussian3d);

    printf("########################################################\n");
    printf("###         FixedLocationLegendreInterp3d           ###\n");
    printf("###             Interpolation error                 ###\n");
    printf("###      getInterpolatedValue(gaussian3d)           ###\n\n");
    printf("Interpolation order : %d \n",interpolationOrder);
    printf("Evaluation point    : (%7.3e, %7.3e, %7.3e) \n\n",xPos,yPos,zPos);

    printf("Exact        : %10.5e \n",exact);
    printf("Interpolated : %10.5e \n\n",interp);
    printf("Error (Abs)  : %10.5e \n",abs(exact-interp));
    printf("Error (Rel)  : %10.5e \n",abs(exact-interp)/abs(exact));
    printf("\n\n");

    //#####################################################
    // Encapsulate the LegendreInterp3d instance inside a
    // lambda function and use it to initialize the values
    // of a UCLAQ::GridFunction3d that has mesh size 1/2
    // the mesh size of the grid providing the interpolant
    // data.
    //#####################################################

    std::function < double(double,double,double) > interpFunction3d = [&interp3d, &gaussian3d](double x, double y, double z)
    {
    return interp3d.evaluateInterpolant(x,y,z,gaussian3d);
    };

    // Create a grid function with 1/2 the mesh size of the grid
    // used to provide the values being interpolated.

	UCLAQ::GridFunction3d refinedInterpolated3d(2*xPanels,xMin,xMax,2*yPanels,yMin,yMax,2*zPanels,zMin,zMax);
	refinedInterpolated3d.setToValue(0.0);

	UCLAQ::GridFunction3d refinedExact3d(2*xPanels,xMin,xMax,2*yPanels,yMin,yMax,2*zPanels,zMin,zMax);
	refinedExact3d.setToValue(0.0);

	UCLAQ::GridFunction3d refinedError(2*xPanels,xMin,xMax,2*yPanels,yMin,yMax,2*zPanels,zMin,zMax);
	refinedError.setToValue(0.0);


    printf("********************************************\n");
    printf("Evaluating interpolant on 2X refined grid. \n");
    printf("In 3d this can a while ... ");
    fflush(stdout);

	refinedInterpolated3d.specify(interpFunction3d);
    refinedExact3d.specify(gaussianFunction3d);

    printf("done. \n");
    printf("********************************************\n\n");

    refinedError = refinedExact3d - refinedInterpolated3d;

    printf("############################################\n");
    printf("###         LegendreInterp3d             ###\n");
    printf("### Interpolation error of evaluation    ###\n");
    printf("###   at nodes of a 2X refined grid      ###\n\n");
    printf("Interpolation order : %d \n\n",interpolationOrder);

    printf("Error (Inf)  : %10.5e \n",refinedError.normInf());
    printf("Error (L2)   : %10.5e \n",refinedError.norm2());
    printf("\n\n\n");

#ifdef _OUTPUT_PLOTS
    UCLAQ::GridFunction3dUtility gUtility3d;
    gUtility3d.outputDataToVTKfile(refinedInterpolated3d,"interp3d.vtk","Interpolated");
    gUtility3d.outputDataToVTKfile(refinedExact3d,"exact3d.vtk","Exact");
    gUtility3d.outputDataToVTKfile(refinedError,"error3d.vtk","Error");
#endif

    printf("XXXX Execution Complete XXXXX\n");



}

