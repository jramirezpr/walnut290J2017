#include <iostream>
#include <cstdio>
using namespace std;

//
// Specifying directory prefix one need only specify
// a single include path to the 290J_2015 directory on the command line
// or in the makescript.
//
//
#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"
#include "GridFunctionNd/UCLAQ_GridFunction2d.h"          // Required for extracting slices of 3D grid function
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"   // Output of 2d slice


#include "MollifierNd/HighOrderPolyPotential3d.h"           // Mollified nuclear potential

//
// A sample program demonstrating the construction  
// and use of the HighOrderPolyPotential3d class.
// 
// The test program assumes the directories 290J_2015 and 290J_Samples are organized as follows
//
// Course Directory
// ---------------------------------------------------------------
//      |                |                |            |
// Project_1 Project_2   ....        290J_Samples   290J_2015
//
//
// where the directory 290J_2015 contains the
//
// 	DoubleVectorNd
// 	GridFunctionNd
// 	Mollifier
//
// source directories.
//
// The command line compilation command is
//
// g++ MollifiedNuclearPotTest.cpp -std=c++11 -I../../290J_2015 -o MollifiedNuclearPotTest.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/MollifiedNuclearPotTest is
//
// make -f MollifiedNuclearPotTest.mk release
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Oct. 27, 2015
//
int main()
{
	double xMin = -5.0; double xMax = 5.0;
    double yMin = -5.0; double yMax = 5.0;
    double zMin = -5.0; double zMax = 5.0;

    long xPanels = 50;
    long yPanels = 50;
    long zPanels = 50;

    // Set up grid function to hold the smoothed nuclear potential values

	UCLAQ::GridFunction3d nuclearPotential(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	nuclearPotential.setToValue(0.0);

	// The far field potential supplied by a HighOrderPolyPotential3d instance
	// is the potential
	//
	// V(r) = -(strength/(4*PI*|r - rPos|)
	//
	// where |r - rPos| = radial distance to rPos = (xPos,yPos,zPos)
	//
	// (For |r - rPos| < radius a smoothed value of this function is returned).
	//
	// In atomic units, the single particle Schroedinger equation for the hydrogen
	// atom has the potential
	//
	// -1/(|r - rPos|)
	//
	// so we set the strength to be 4*PI
	//

	double pi4      = 4.0*3.14159265358979323846;
	double strength = pi4;

	// (xPos,yPos,zPos) : location of potential

	double xPos   =  0.1;
	double yPos   =  0.01;
	double zPos   = -0.1;

	// radius : smoothing radius (epsilon)

	double radius =  1.5;

	// Instantiate and initialize the smoothed potential

	HighOrderPolyPotential3d smoothPotentialFun;

	smoothPotentialFun.initialize(xPos,yPos,zPos,radius,strength);

    // Evaluate the smoothed nuclear potential on a grid

    nuclearPotential.specify(smoothPotentialFun.getEvaluationPtr());

    // Output 2D slice 

    long zIndex = zPanels/2;

	UCLAQ::GridFunction2d zSliceFun =nuclearPotential.getConstantZslice(zIndex); // (x-y function) slice
    
    UCLAQ::GridFunction2dUtility gUtility2d;
    gUtility2d.outputToGNUplot(zSliceFun,"Potential_XY.dat","%10.5f");

    cout << "Constant Z function slice output to " << "Potential_XY.dat" << endl << endl;

    UCLAQ::GridFunction3dUtility gUtility3d;

    gUtility3d.outputDataToVTKfile(nuclearPotential,"Potential.vtk", "Nuclear_Potential");
    cout << "VTK output to " << "Potential.vtk" << endl << endl;
    //
    // Output file that can be viewed using mayavi2 (or another

	return 0;
}
