#include <cmath>
#include <functional>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "MollifierNd/SmoothPolyStep.h"

#ifndef _RectMaskFun1d_
#define _RectMaskFun1d_

#define _DEFAULT_STEP_DIFFERENTIABLITY_ 6

//
//     0             0 <= mask <= 1                  1                 0 <= mask <= 1          0
// | ----  | --- transition distance --- | --- domain interior | --- transition distance --- | --- |
//    maskMinX                                                                            maskMaxX

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
class RectMaskFun1d
{
public:

    RectMaskFun1d() : smoothStepFunction()
	{
    initialize();
	}

    RectMaskFun1d(double maskMinX, double maskMaxX, double transitionDistance,
    int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    initialize(maskMinX,maskMaxX,transitionDistance,polyDiffOrder);
    }

    virtual  ~RectMaskFun1d(){};

    void initialize()
    {
    polyDiffOrder = _DEFAULT_STEP_DIFFERENTIABLITY_;
    maskMinX = 0.0;
    maskMaxX = 1.0;
    transitionDistance = 0.0;
    smoothStepFunction.initialize();
    smoothStepFunction.setTransitionDistance(transitionDistance);
    }

    void initialize(double maskMinX, double maskMaxX, double transitionDistance,
    int polyDiffOrder = _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    this-> polyDiffOrder     = polyDiffOrder;
    this->maskMinX           = maskMinX;
    this->maskMaxX           = maskMaxX;
    this->transitionDistance = transitionDistance;

    smoothStepFunction.initialize();
    smoothStepFunction.setDifferentiability(polyDiffOrder);
    smoothStepFunction.setTransitionDistance(transitionDistance);
    }

    void setDiffFactor(int diffOrder)
	{
    polyDiffOrder = diffOrder;
    smoothStepFunction.setDifferentiability(polyDiffOrder);
	}


	double operator()(double xPos) const
	{
	double  dist;

	if((xPos-maskMinX) < (maskMaxX- xPos))
	{
		dist = xPos-maskMinX;
	}
	else
	{
		dist = maskMaxX-xPos;
	}
	return smoothStepFunction(dist - (transitionDistance/2.0));
	}

	//  Returns a std::function that is bound to the evaluation operator of *this

#if __cplusplus > 199711L
	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return std::move(F);
	}
#endif


    double                  maskMinX;
    double                  maskMaxX;
    double        transitionDistance;
    SmoothPolyStep smoothStepFunction;
    int polyDiffOrder;
};

#undef _DEFAULT_STEP_DIFFERENTIABLITY_
#endif
