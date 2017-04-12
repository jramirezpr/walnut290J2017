#include <cmath>
#include <functional>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "MollifierNd/SmoothPolyStep.h"
#include "RectMaskFun1d.h"

#ifndef _RectMaskFun3d_
#define _RectMaskFun3d_

#define _DEFAULT_STEP_DIFFERENTIABLITY_ 6
//
// A product of 3 RectMaskFun1d instances, each of which implements a masking function
// with the structure:
//
//     0             0 <= mask <= 1                  1                 0 <= mask <= 1          0
// | ----  | --- transition distance --- | --- domain interior | --- transition distance --- | --- |
//    maskMin                                                                            maskMax
//
// for each dimension
//
/*
#############################################################################
#
# Copyright  2016 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

class RectMaskFun3d
{
public:

    RectMaskFun3d() : maskFunX(), maskFunY(),maskFunZ()
	{
    initialize();
	}

    RectMaskFun3d(double maskMinX, double maskMaxX, double transitionDistanceX,
                   double maskMinY, double maskMaxY, double transitionDistanceY,
                   double maskMinZ, double maskMaxZ, double transitionDistanceZ,
                   int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    initialize(maskMinX,maskMaxX,transitionDistanceX,
               maskMinY,maskMaxY,transitionDistanceY,
               maskMinZ,maskMaxZ,transitionDistanceZ,polyDiffOrder);
    }

    virtual  ~RectMaskFun3d(){};

    void initialize()
    {
    maskFunX.initialize();
    maskFunY.initialize();
    maskFunZ.initialize();
    }

    void initialize(double maskMinX, double maskMaxX, double transitionDistanceX,
                      double maskMinY, double maskMaxY, double transitionDistanceY,
                      double maskMinZ, double maskMaxZ, double transitionDistanceZ,
                      int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    maskFunX.initialize(maskMinX,maskMaxX,transitionDistanceX,polyDiffOrder);
    maskFunY.initialize(maskMinY,maskMaxY,transitionDistanceY,polyDiffOrder);
    maskFunZ.initialize(maskMinZ,maskMaxZ,transitionDistanceZ,polyDiffOrder);
    }

    void setDiffFactor(int diffOrder)
	{
    maskFunX.setDiffFactor(diffOrder);
    maskFunY.setDiffFactor(diffOrder);
    maskFunZ.setDiffFactor(diffOrder);
	}

	double operator()(double xPos,double yPos,double zPos) const
	{
    return maskFunX(xPos)*maskFunY(yPos)*maskFunZ(zPos);
	}

	//  Returns a std::function that is bound to the evaluation operator of *this

#if __cplusplus > 199711L
	std::function<double(double,double,double)> getEvaluationPtr() const
	{
	std::function<double(double,double,double)> F = [this](double x,double y,double z) {return this->operator()(x,y,z);};
	return std::move(F);
	}
#endif

	RectMaskFun1d maskFunX;
	RectMaskFun1d maskFunY;
    RectMaskFun1d maskFunZ;
};

#undef _DEFAULT_STEP_DIFFERENTIABLITY_
#endif
