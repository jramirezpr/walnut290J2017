#include <cmath>
#include <functional>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "MollifierNd/SmoothPolyStep.h"
#include "RectMaskFun1d.h"

#ifndef _RectMaskFun2d_
#define _RectMaskFun2d_

#define _DEFAULT_STEP_DIFFERENTIABLITY_ 6

//
// A product of 2 RectMaskFun1d instances, each of which implements a masking function
// with the structure:
//
//     0             0 <= mask <= 1                  1                 0 <= mask <= 1           0
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

class RectMaskFun2d
{
public:

    RectMaskFun2d() : maskFunX(), maskFunY()
	{
    initialize();
	}

    RectMaskFun2d(double maskMinX, double maskMaxX, double transitionDistanceX,
    double maskMinY, double maskMaxY, double transitionDistanceY,
    int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    initialize(maskMinX,maskMaxX,transitionDistanceX,
               maskMinY,maskMaxY,transitionDistanceY, polyDiffOrder);
    }

    virtual  ~RectMaskFun2d(){};

    void initialize()
    {
    maskFunX.initialize();
    maskFunY.initialize();
    }

    void initialize(double maskMinX, double maskMaxX, double transitionDistanceX,
                      double maskMinY, double maskMaxY, double transitionDistanceY,
                      int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    maskFunX.initialize(maskMinX,maskMaxX,transitionDistanceX,polyDiffOrder);
    maskFunY.initialize(maskMinY,maskMaxY,transitionDistanceY,polyDiffOrder);
    }

    void setDiffFactor(int diffOrder)
	{
    maskFunX.setDiffFactor(diffOrder);
    maskFunY.setDiffFactor(diffOrder);
	}


	double operator()(double xPos,double yPos) const
	{
    return maskFunX(xPos)*maskFunY(yPos);
	}

	//  Returns a std::function that is bound to the evaluation operator of *this

#if __cplusplus > 199711L
	std::function<double(double,double)> getEvaluationPtr() const
	{
	std::function<double(double,double)> F = [this](double x,double y) {return this->operator()(x,y);};
	return std::move(F);
	}
#endif

	RectMaskFun1d maskFunX;
	RectMaskFun1d maskFunY;

};

#undef _DEFAULT_STEP_DIFFERENTIABLITY_
#endif
