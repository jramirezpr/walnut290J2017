#include <cmath>
#include <functional>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "ValueStepFun.h"

#ifndef _ValueMaskFun2d_
#define _ValueMaskFun2d_

#define _DEFAULT_STEP_DIFFERENTIABLITY_ 6

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

class ValueMaskFun2d
{
public:

    ValueMaskFun2d(int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_) :  valueStepFun()
    {
    this-> polyDiffOrder = polyDiffOrder;
    }

    virtual  ~ValueMaskFun2d(){};


    void initialize(int polyDiffOrder =  _DEFAULT_STEP_DIFFERENTIABLITY_)
    {
    this-> polyDiffOrder     = polyDiffOrder;
    valueStepFun.initialize();
    }

    void setDiffFactor(int diffOrder)
	{
    polyDiffOrder = diffOrder;
	}

     /*!
     * This member function sets the values of maskFun to be
     * 0 when F < funMin and equal to 1 when  F > funMax.
     * For values funMin < F < funMax, maskFun is assigned values of a function that
     * transitions between 0 and 1 based on the values of F. This
     * transition is smooth if F is smooth.
     *
    */
    void  createMinMask(double funMin, double funMax,const UCLAQ::GridFunction2d& F,UCLAQ::GridFunction2d& maskFun)
    {
    	valueStepFun.initialize(funMin,funMax,0.0,1.0,polyDiffOrder);
    	maskFun.initialize(F);
    	maskFun.transformValues(valueStepFun.getEvaluationPtr());
    }

     /*!
     * This member function sets the values of maskFun to be
     * 1 when F < funMin and equal to 0 when  F > funMax.
     * For values funMin < F < funMax, maskFun is assigned values of a function that
     * transitions between 1 and 0 based on the values of F. This
     * transition is smooth if F is smooth.
     *
    */
    void  createMaxMask(double funMin, double funMax,const UCLAQ::GridFunction2d& F,UCLAQ::GridFunction2d& maskFun)
    {
    	valueStepFun.initialize(funMin,funMax,1.0,0.0,polyDiffOrder);
    	maskFun.initialize(F);
    	maskFun.transformValues(valueStepFun.getEvaluationPtr());
    }


    /*!
     * This member function sets the values of the input grid function F to
     * 0 when  F < funMin and is the identity when F > funMax.
     * For values funMin < F < funMax, F is assigned values of a function that
     * transitions between 0 and funMax based on the values of F. This
     * transition is smooth if F is smooth.
     *
     */
    void applyMinMask(double funMin, double funMax, UCLAQ::GridFunction2d& F)
    {
    	valueStepFun.initialize(funMin,funMax,0.0,1.0,polyDiffOrder);
    	mask.initialize(F);

    	// Create mask

    	mask.transformValues(valueStepFun.getEvaluationPtr());

    	// Apply mask to F

    	F *= mask;
    }

     /*!
     * This member function sets the values of the input grid function F to
     * 0 when  F > funMax and is the identity when F < funMin.
     * For values funMin < F < funMax, F is assigned values of a function that
     * transitions between funMin and 0 based on the values of F. This
     * transition is smooth if F is smooth.
     *
     */
    void applyMaxMask(double funMin, double funMax, UCLAQ::GridFunction2d& F)
    {
    	valueStepFun.initialize(funMin,funMax,1.0,0.0,polyDiffOrder);
    	mask.initialize(F);

    	// Create mask

    	mask.transformValues(valueStepFun.getEvaluationPtr());

    	// Apply mask to F

    	F *= mask;
    }

     /*!
     * This member function sets the values Fout to
     * 0 when  F < funMin and equal to F when  F > funMax.
     * For values funMin < F < funMax, Fout is assigned values of a function that
     * transitions between 0 and funMax based on the values of F. This
     * transition is smooth if F is smooth.
     *
     */
    void applyMinMask(double funMin, double funMax,const UCLAQ::GridFunction2d& F,UCLAQ::GridFunction2d& Fout)
    {
    	valueStepFun.initialize(funMin,funMax,0.0,1.0,polyDiffOrder);
    	Fout.initialize(F);

    	// Create mask

    	Fout.transformValues(valueStepFun.getEvaluationPtr());

    	// Apply mask to F

    	Fout *= F;
    }

     /*!
     * This member function sets the values of Fout to
     * 0 when  F > funMax and equal to F when F  < funMin.
     * For values funMin < F < funMax, Fout is assigned values of a function that
     * transitions between funMin and 0 based on the values of F. This
     * transition is smooth if F is smooth.
     *
     */
    void applyMaxMask(double funMin, double funMax, const UCLAQ::GridFunction2d& F,UCLAQ::GridFunction2d& Fout)
    {
    	valueStepFun.initialize(funMin,funMax,1.0,0.0,polyDiffOrder);
    	Fout.initialize(F);

    	// Create mask

    	Fout.transformValues(valueStepFun.getEvaluationPtr());

    	// Apply mask to F

    	Fout *= F;
    }



    UCLAQ::GridFunction2d       mask;
    ValueStepFun        valueStepFun;
    int polyDiffOrder;
};

#undef _DEFAULT_STEP_DIFFERENTIABLITY_
#endif
