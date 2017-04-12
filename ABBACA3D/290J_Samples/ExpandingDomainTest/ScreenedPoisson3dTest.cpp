#include <iostream>
#include <string>
#include <cstdio>
#include <functional>


#include "ExpandingDomainNd/InversePoisson3d.h"
#include "ExpandingDomainNd/ScreenedSoln3d.h"
#include "MollifierNd/SmoothPolyMollifier3d.h"

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

using namespace std;

#undef  _OUTPUT_PLOTS // #define to output 3d data for VTK

int main(int argc, char* argv[])
{

	int threadCount = -1;
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

//
// The extension factor is the factor by which the
// x-coordinate domain (1st coordinate) is implicitly
// expanded in the InversePoisson3d instance.
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
    double extensionFactor = 2.0;

	double xMin = -1.0;
	double xMax  = 1.0;

	double yMin = -1.0;
	double yMax =  1.0;

    double zMin = -1.0;
	double zMax =  1.0;

    long xPanels = 50;
    long yPanels = 50;
    long zPanels = 50;

	double xBar     = -0.1;
	double yBar     =  0.2;
	double zBar     =  0.3;

	double sigma    =  0.4;
	double strength =  1.0;

    double laplaceCoeff            =  2.0;
    double screenCoeff             = -1.0;
	int    sourceDifferentiability =    7; // Source differentiability range = 0 to 9

	// Instantiate an exact potential with a specified source differentiability

	ScreenedSoln3d screenedPotFunction3d(xBar,yBar, zBar, sigma, strength, laplaceCoeff, screenCoeff);
	screenedPotFunction3d.setSourceDifferentiability(sourceDifferentiability);

	UCLAQ::GridFunction3d source3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	source3d.setToValue(0.0);

	UCLAQ::GridFunction3d potential3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	potential3d.setToValue(0.0);

    UCLAQ::GridFunction3d exactPotential3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	exactPotential3d.setToValue(0.0);

	UCLAQ::GridFunction3d error3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	error3d.setToValue(0.0);

	// Capture the exact potential and source values on GridFunction3d
	// instances

	exactPotential3d.specify(screenedPotFunction3d.getEvaluationPtr());
	source3d.specify(screenedPotFunction3d.getSourceEvaluationPtr());

    //#####################################################
    // Initialize InversePoisson3d instance
    //#####################################################

    InversePoisson3d inversePoissonOp;
    inversePoissonOp.initialize(laplaceCoeff,screenCoeff, xPanels,  xMin, xMax, yPanels, yMin, yMax, zPanels, zMin, zMax,extensionFactor);

    potential3d = source3d;

    inversePoissonOp.applyInverseOp(potential3d);

    error3d = potential3d - exactPotential3d;

    printf("######################################################\n");
    printf("###         InversePoisson3d                       ###\n");
    printf("### Infinite domain screened potential error       ###\n\n");

    printf("Source differentiabiliy : %d \n\n",sourceDifferentiability);

    printf("Comp. Domain Size       : %8.5f X %-8.5f X %-8.5f  \n", xMax-xMin,yMax-yMin,zMax-zMin);
    printf("HX                      : %8.5f \n",(xMax-xMin)/(double)xPanels);
    printf("HY                      : %8.5f \n\n",(yMax-yMin)/(double)yPanels);
    printf("HZ                      : %8.5f \n\n",(zMax-zMin)/(double)zPanels);
    printf("Extension Factor        : %8.5f  \n\n",extensionFactor);

    printf("Error (Inf)  : %10.5e \n",error3d.normInf());
    printf("Error (L2)   : %10.5e \n",error3d.norm2());
    printf("\n");
    printf("Soln norm (Inf)  : %10.5e \n",potential3d.normInf());
    printf("Soln norm (L2)   : %10.5e \n",potential3d.norm2());
    printf("\n\n\n");

#ifdef _OUTPUT_PLOTS
    UCLAQ::GridFunction3dUtility gUtility3d;
    gUtility3d.outputDataToVTKfile(potential3d,"computed3d.vtk","Potential");
    gUtility3d.outputDataToVTKfile(exactPotential3d,"exact3d.vtk","ExactPotential");
    gUtility3d.outputDataToVTKfile(error3d,"error3d.vtk","Error");
    gUtility3d.outputDataToVTKfile(source3d,"source3d.vtk","Error");
#endif

//
//  Increase extension factor to improve accuracy
//
    extensionFactor *= 2.0;
    inversePoissonOp.initialize(laplaceCoeff, screenCoeff, xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,extensionFactor);

    potential3d = source3d;

    inversePoissonOp.applyInverseOp(potential3d);

    error3d = potential3d - exactPotential3d;

    printf("######################################################\n");
    printf("###         InversePoisson3d                       ###\n");
    printf("### Infinite domain screened potential error       ###\n\n");
    printf("### using 2 X ExtensionFactor                      ###\n\n");

    printf("Source differentiabiliy : %d \n\n",sourceDifferentiability);

    printf("Comp. Domain Size       : %8.5f X %-8.5f X %-8.5f  \n", xMax-xMin,yMax-yMin,zMax-zMin);
    printf("HX                      : %8.5f \n",(xMax-xMin)/(double)xPanels);
    printf("HY                      : %8.5f \n",(yMax-yMin)/(double)yPanels);
    printf("HZ                      : %8.5f \n\n",(zMax-zMin)/(double)zPanels);
    printf("New Extension Factor    : %8.5f  \n\n",extensionFactor);

    printf("Error (Inf)  : %10.5e \n",error3d.normInf());
    printf("Error (L2)   : %10.5e \n",error3d.norm2());
    printf("\n");
    printf("Soln norm (Inf)  : %10.5e \n",potential3d.normInf());
    printf("Soln norm (L2)   : %10.5e \n",potential3d.norm2());
    printf("\n\n\n");


    printf("XXXX Execution Complete XXXXX\n");

}

