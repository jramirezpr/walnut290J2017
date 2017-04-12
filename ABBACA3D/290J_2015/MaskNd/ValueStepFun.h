#include <cmath>
#include <functional>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "MollifierNd/SmoothPolyStep.h"

#ifndef _ValueStepFun_
#define _ValueStepFun_

#define _DEFAULT_STEP_DIFFERENTIABLITY_ 6


//
//
// The operator()(double v) of this class returns
//
//                 v <= funMin  : maskFunMin
//
//        funMin < v < funMax   : smooth transition from
//                                maskFunMin to maskFunMax
//
//              funMax <= v     : maskFunMax
//
//
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

class ValueStepFun
{
public:

    ValueStepFun() : smoothStepFunction()
	{
    initialize();
	}

    ValueStepFun(double funMin, double funMax, double maskFunMin, double maskFunMax, int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    initialize(funMin, funMax, maskFunMin,maskFunMax,polyDiffOrder);
    }

    virtual  ~ValueStepFun(){};


    void initialize()
    {
    polyDiffOrder = _DEFAULT_STEP_DIFFERENTIABLITY_;
    funMin             = 0.0;
    funMax             = 1.0;
    maskFunMin         = 0.0;
    maskFunMax         = 1.0;
    centValue          = 0.5;
    transitionDistance = 1.0;
    smoothStepFunction.initialize();
    smoothStepFunction.initialize(centValue,transitionDistance,1.0,false);
    }

    void initialize(double funMin, double funMax, double maskFunMin, double maskFunMax,
    int polyDiffOrder = _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    this->polyDiffOrder      = polyDiffOrder;
    this->maskFunMin         = maskFunMin;
    this->maskFunMax         = maskFunMax;
    this->funMin             = funMin;
    this->funMax             = funMax;

    this->centValue          = (funMax + funMin)/2.0;
    this->transitionDistance = abs(funMax-funMin);

    bool  stepDownFlag  = false;
    double stepStrength = 1.0;

    smoothStepFunction.initialize(centValue,transitionDistance,stepStrength,stepDownFlag);
    smoothStepFunction.setDifferentiability(polyDiffOrder);
    }

    void setDiffFactor(int diffOrder)
	{
    polyDiffOrder = diffOrder;
    smoothStepFunction.setDifferentiability(polyDiffOrder);
	}

	double operator()(double val) const
	{
    return maskFunMin + (maskFunMax-maskFunMin)*smoothStepFunction(val);
	}

	//  Returns a std::function that is bound to the evaluation operator of *this

#if __cplusplus > 199711L
	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return std::move(F);
	}
#endif

    double                      funMin;
    double                      funMax;
    double                  maskFunMin;
    double                  maskFunMax;

    double                   centValue;
    double          transitionDistance;
    SmoothPolyStep  smoothStepFunction;
    int                  polyDiffOrder;
};

#undef _DEFAULT_STEP_DIFFERENTIABLITY_
#endif   
