#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <functional>

using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"
#include "MaskNd/RectMaskFun1d.h"

int main(int argc, char* argv[])
{
    double xPanel =  100;
    double xMin    = -2.5;
    double xMax    = 3.0;

    double maskMinX           = xMin +  0.1;
    double maskMaxX           = xMax -  0.2;
    double transitionDistance = 0.3;

    UCLAQ::GridFunction1d gridMask(xPanel,xMin,xMax);

    RectMaskFun1d rectMaskFun1d(maskMinX,maskMaxX,transitionDistance);
    gridMask.specify(rectMaskFun1d.getEvaluationPtr());

    UCLAQ::GridFunction1dUtility gUtility;
    gUtility.outputToGNUplot(gridMask,"mask1d.dat");



    printf("XXXX Execution Complete XXXXX\n");

}

