#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <functional>

using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"
#include "MaskNd/RectMaskFun2d.h"

int main(int argc, char* argv[])
{
    double xPanel =  100;
    double xMin    = -2.5;
    double xMax    = 3.0;

    double yPanel =   50;
    double yMin    =  0.0;
    double yMax    =  2.0;

    double maskMinX            = xMin +  0.1;
    double maskMaxX            = xMax -  0.2;
    double transitionDistanceX = 0.3;


    double maskMinY            = yMin +  0.1;
    double maskMaxY            = yMax -  0.2;
    double transitionDistanceY = 0.3;

    UCLAQ::GridFunction2d gridMask(xPanel,xMin,xMax,yPanel,yMin,yMax);

    RectMaskFun2d rectMaskFun2d(maskMinX,maskMaxX,transitionDistanceX,
                                maskMinY,maskMaxY,transitionDistanceY);

    gridMask.specify(rectMaskFun2d.getEvaluationPtr());


    UCLAQ::GridFunction2dUtility gUtility;
    gUtility.outputToGNUplot(gridMask,"mask2d.dat");



    printf("XXXX Execution Complete XXXXX\n");

}

