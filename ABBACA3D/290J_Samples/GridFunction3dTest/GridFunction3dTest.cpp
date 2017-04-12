#include <iostream>
#include <cstdio>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2d.h"          // Required for extracting slices of 3D grid function
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"   // Output of 2d slice

//
// A sample program demonstrating the construction  
// and use of the UCLAQ::GridFunction3d class. 
// 
// The test program assumes the directories 290J_2015 and 290J_Samples are organized as follows
//
// Course Directory
// ---------------------------------------------------------------
//      |                |                |            |
// Project_1 Project_2   ....        290J_Samples   290J_2015
//
//
// where the directory 290J_2015 contains the DoubleVectorNd and GridFunctionNd
// source directories.
//
// The command line compilation command is
//
// g++ GridFunction3dTest.cpp -std=c++11 -I../../290J_2015 -o GridFunction3dTest.exe
//
// Alternately, if one is using MakeScripts, the build command executed from within 290J_Samples/GridFunction3dTest is
//
// make -f GridFunction3dTest.mk release
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Oct. 27, 2015
//
int main()
{
	double xMin = -1.0; double xMax = 1.0;
    double yMin = -1.0; double yMax = 1.0;
    double zMin = -1.0; double zMax = 1.0;

    long xPanels = 10;
    long yPanels = 15;
    long zPanels = 20;

	UCLAQ::GridFunction3d Fgrid(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	UCLAQ::GridFunction3d Ggrid(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d FGgridSum(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

    // Initialize grid values to 0

    Fgrid.setToValue(0.0);
    Ggrid.setToValue(0.0);


    // Initialize grid functions using lambda functions (c++11 construct)

	std::function<double(double,double,double)> F      
      = [](double x,double y,double z) {return 2.0*x*x*(z+1.0)*(z+1.0);};

    // To pass an "external" value into the lambda function, 
    // insert it into the  [ ]'s

    double alpha = 1.0;

    std::function<double(double,double,double)> G      
      = [alpha](double x,double y,double z) {return alpha*y;};

    std::function<double(double,double,double)> FplusG 
     = [alpha](double x,double y,double z) {return 2.0*x*x*(z+1.0)*(z+1.0) + alpha*y;};


    Fgrid.specify(F);
    Ggrid.specify(G);
    FGgridSum.specify(FplusG);

    //
    // Output 2D slice 
    //
    long zIndex = zPanels/2;

	UCLAQ::GridFunction2d zSliceFun = Fgrid.getConstantZslice(zIndex); // (x-y function) slice
    
    UCLAQ::GridFunction2dUtility gUtility2d;
    gUtility2d.outputToGNUplot(zSliceFun,"Fgrid_XY.dat","%10.5f");

    cout << "Constant Z function slice output to " << "Fgrid_XY.dat" << endl << endl;

	//
    // Check addition
    //

	Fgrid     += Ggrid;
    FGgridSum -= Fgrid;

    cout << "L2 Norm of Fgrid + Ggrid (mesh weighted) " << Fgrid.norm2() << endl;
    cout << "L2 Error   Fgrid + Ggrid (mesh weighted) " << FGgridSum.norm2() << endl;

	return 0;
}
