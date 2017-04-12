#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <functional>

using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"
#include "MaskNd/RectMaskFun3d.h"

int main(int argc, char* argv[])
{
    double xPanel =  50;
    double xMin    = -2.5;
    double xMax    = 3.0;

    double yPanel =   50;
    double yMin    =  0.0;
    double yMax    =  2.0;

    double zPanel =   50;
    double zMin    =  -1.0;
    double zMax    =   2.0;

    double maskMinX            = xMin +  0.1;
    double maskMaxX            = xMax -  0.2;
    double transitionDistanceX = 0.3;

    double maskMinY            = yMin +  0.1;
    double maskMaxY            = yMax -  0.2;
    double transitionDistanceY = 0.3;

    double maskMinZ            = zMin +  0.1;
    double maskMaxZ            = zMax -  0.2;
    double transitionDistanceZ = 0.3;


    UCLAQ::GridFunction3d gridMask(xPanel,xMin,xMax,yPanel,yMin,yMax,zPanel,zMin,zMax);

    RectMaskFun3d rectMaskFun3d(maskMinX,maskMaxX,transitionDistanceX,
                                  maskMinY,maskMaxY,transitionDistanceY,
                                  maskMinZ,maskMaxZ,transitionDistanceZ);

    gridMask.specify(rectMaskFun3d.getEvaluationPtr());


    UCLAQ::GridFunction3dUtility gUtility;
    gUtility.outputDataToVTKfile(gridMask,"mask3d.vtk","mask");



    printf("XXXX Execution Complete XXXXX\n");

}

