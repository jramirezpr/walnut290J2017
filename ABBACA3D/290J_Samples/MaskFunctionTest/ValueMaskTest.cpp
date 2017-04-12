#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <functional>

using namespace std;
#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

#include "MollifierNd/SmoothPolyPotential1d.h"
#include "MollifierNd/SmoothPolyPotential2d.h"
#include "MollifierNd/SmoothPolyPotential3d.h"

#include "MaskNd/ValueMaskFun1d.h"
#include "MaskNd/ValueMaskFun2d.h"
#include "MaskNd/ValueMaskFun3d.h"

int main(int argc, char* argv[])
{

    UCLAQ::GridFunction1dUtility gUtility;

    SmoothPolyPotential1d potential1d(0.0,0.25,1.0);
    SmoothPolyPotential2d potential2d(0.0,0.0,0.25,1.0);
    SmoothPolyPotential3d potential3d(0.0,0.0,0.0,0.25,1.0);

    std::function<double(double)>                 sourcePotential = potential1d.getEvaluationPtr();
    std::function<double(double,double)>        sourcePotential2d = potential2d.getEvaluationPtr();
    std::function<double(double,double,double)> sourcePotential3d = potential3d.getEvaluationPtr();

    long   xPanel = 50;
    long   yPanel = 50;
    long   zPanel = 50;

    double xMin    = -1.0;
    double xMax    =  1.0;

    double yMin    = -1.0;
    double yMax    =  1.0;

    double zMin    = -1.0;
    double zMax    =  1.0;

	double funMin;
    double funMax;


    UCLAQ::GridFunction1d potFunMask(xPanel,xMin,xMax);
    UCLAQ::GridFunction1d potFun(xPanel,xMin,xMax);
    UCLAQ::GridFunction1d maskFun(xPanel,xMin,xMax);

    potFun.specify(sourcePotential);
    gUtility.outputToGNUplot(potFun,"potData.dat");


    ValueMaskFun1d valueMaskFun1d;
    funMin = 0.2;
    funMax = 0.4;

    valueMaskFun1d.createMaxMask(funMin,funMax,potFun,maskFun);
    gUtility.outputToGNUplot(maskFun,"maxMask.dat");


    valueMaskFun1d.createMinMask(funMin,funMax,potFun,maskFun);
    gUtility.outputToGNUplot(maskFun,"minMask.dat");


    UCLAQ::GridFunction2dUtility gUtility2d;
    ValueMaskFun2d valueMaskFun2d;

    UCLAQ::GridFunction2d potFunMask2d(xPanel,xMin,xMax,yPanel,yMin,yMax);
    UCLAQ::GridFunction2d potFun2d(xPanel,xMin,xMax,yPanel,yMin,yMax);
    UCLAQ::GridFunction2d maskFun2d(xPanel,xMin,xMax,yPanel,yMin,yMax);

    funMin = -0.1;
    funMax = 0.05;

    potFun2d.specify(sourcePotential2d);
    gUtility2d.outputToGNUplot(potFun2d,"potData2d.dat");


    valueMaskFun2d.createMaxMask(funMin,funMax,potFun2d,maskFun2d);
    gUtility2d.outputToGNUplot(maskFun2d,"maxMask2d.dat");


    valueMaskFun2d.createMinMask(funMin,funMax,potFun2d,maskFun2d);
    gUtility2d.outputToGNUplot(maskFun2d,"minMask2d.dat");


    UCLAQ::GridFunction3dUtility gUtility3d;
    ValueMaskFun3d valueMaskFun3d;

    UCLAQ::GridFunction3d potFunMask3d(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);
    UCLAQ::GridFunction3d potFun3d(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);
    UCLAQ::GridFunction3d maskFun3d(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);


    funMin = -0.3;
    funMax = -0.1;


    potFun3d.specify(sourcePotential3d);
    gUtility3d.outputDataToVTKfile(potFun3d,"potData3d.vtk","potData");

    valueMaskFun3d.createMaxMask(funMin,funMax,potFun3d,maskFun3d);
    gUtility3d.outputDataToVTKfile(maskFun3d,"maxMask3d.vtk","maxMask");


    valueMaskFun3d.createMinMask(funMin,funMax,potFun3d,maskFun3d);
    gUtility3d.outputDataToVTKfile(maskFun3d,"minMask3d.vtk","minMask");


    printf("XXXX Execution Complete XXXXX\n");

}


