#ifndef _SmoothedH2op_
#define _SmoothedH2op_

//##################################################################
//   SmoothedH2op:
//
//   A class whose apply member function implements the single particle
//   operator
//
//   (-1/2)*Delta_h -  V(r)
//
//   where V(r) is a high-order smooth approximation to the 
//   nuclear potential induced by two nucleii of charge 1
//.
//##################################################################

/*
#############################################################################
#
# Copyright 2015-16 Chris Anderson
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

#include <vector> 
using namespace std;


#include "GridFunctionNd/UCLAQ_GridFunction3d.h"        // Grid function class
#include "MollifierNd/HighOrderPolyPotential3d.h"       // To create smoothed nuclear potentials
#include "Qdomain/QdomainKineticOp.h"                   // Variable order Laplace operator approximation

// Default values that will be used to override the defaults
// specified in the HighOrderPolyPotential3d and QdomainKineticOp header 
// files. 


#define _POTENTIAL_MOMENT_ORDER  6
#define _POTENTIAL_SOURCE_DIFFERENTIABILITY  6
#define _LAPLACE_OP_ORDER 6

class SmoothedH2op
{

public:

	SmoothedH2op()
	{
		initialize();
	}

	SmoothedH2op(vector<double>& A, vector<double> &B, double epsilon)
	{
		initialize(A,B,epsilon);
	}

	void initialize()
	{
        this->A.clear();
        this->B.clear();
		this->epsilon = 0.0;

		potentialMomentOrder             = 0;
		potentialSourceDifferentiability = 0;

		nuclearPotential.clear();
		kineticEnergyOp.initialize();

		potential.initialize();
		Vtmp.initialize();
	}

	void initialize(vector<double>& A, vector<double> &B, double epsilon)
	{
		this->A       = A;
		this->B       = B;
		this->epsilon = epsilon;

		double strength = 4.0*3.141592653589793238;
		double opCoeff  = 1.0;

		potentialMomentOrder             = _POTENTIAL_MOMENT_ORDER;
		potentialSourceDifferentiability = _POTENTIAL_SOURCE_DIFFERENTIABILITY;

		nuclearPotential.resize(2);
	    nuclearPotential[0].initialize(A[0],A[1],A[2],epsilon,strength,opCoeff);
		nuclearPotential[0].setSourceDifferentiability(potentialSourceDifferentiability);
		nuclearPotential[0].setOrder(potentialMomentOrder);

	    nuclearPotential[1].initialize(B[0],B[1],B[2],epsilon,strength,opCoeff);
		nuclearPotential[1].setSourceDifferentiability(potentialSourceDifferentiability);
		nuclearPotential[1].setOrder(potentialMomentOrder);

		double laplaceCoeff    =  -0.5;
		laplaceOpOrder         =  _LAPLACE_OP_ORDER;

		kineticEnergyOp.initialize(laplaceOpOrder, laplaceCoeff);

		potential.initialize();
		Vtmp.initialize();
	}


	void setPotentialSmoothingRadius(double epsilon)
	{
		this->epsilon = epsilon;
		nuclearPotential[0].setRadius(epsilon);
		nuclearPotential[1].setRadius(epsilon);
	}

	void setPotentialMomentOrder(int potentialMomentOrder)
	{
		this->potentialMomentOrder = potentialMomentOrder;
		nuclearPotential[0].setOrder(potentialMomentOrder);
		nuclearPotential[1].setOrder(potentialMomentOrder);
	}

	void setPotentialSourceDifferentiability(int potentialSourceDifferentiability)
	{
		this->potentialSourceDifferentiability = potentialSourceDifferentiability;
		nuclearPotential[0].setSourceDifferentiability(potentialSourceDifferentiability);
		nuclearPotential[1].setSourceDifferentiability(potentialSourceDifferentiability);
	}

	void setLaplaceOrder(int laplaceOpOrder)
	{
		double laplaceCoeff    =  -0.5;
		this->laplaceOpOrder   =  laplaceOpOrder;
		kineticEnergyOp.initialize(laplaceOpOrder, laplaceCoeff);
	}

	void apply(UCLAQ::GridFunction3d& V)
	{
		// Update potential
		// (A no-op unless parameters defining the potential have changed)

		updatePotential(V);

		// Create potential component

		Vtmp.initialize(V);
		Vtmp *= potential;

		// Create kinetic energy component

		kineticEnergyOp.apply(V);
		V += Vtmp;
	}


	void updatePotential(UCLAQ::GridFunction3d& V)
	{
		// Change in grid structure (or initial call)

		if(not potential.isCoincident(V))
		{
		potential.initialize(V);
		potential.specify(nuclearPotential[0]);
		potential += nuclearPotential[1];
		}

		// Change in smoothing distance

		if(nuclearPotential[0].getOrder() != potentialMomentOrder)
		{
		potential.specify(nuclearPotential[0]);
		potential += nuclearPotential[1];
		}

		// Change in source differentiability

	    if(nuclearPotential[0].getSourceDifferentiablity() != potentialSourceDifferentiability)
		{
		potential.specify(nuclearPotential[0]);
		potential += nuclearPotential[1];
		}
	}


    vector<double> A;                                  // Location of nucleii
    vector<double> B;

	UCLAQ::QdomainKineticOp   kineticEnergyOp;         // Kinetic energy operator (Laplace operator)
	int                        laplaceOpOrder;         // Order of difference approximation to the Laplace operator: 2,4, 6, or 8

	vector<HighOrderPolyPotential3d> nuclearPotential;  // Potentials induced by a "high order" approximation to a delta function
	double                                    epsilon; // Smoothing radius
	int                          potentialMomentOrder; // Order of vanishing moments: 2,4,6
	int              potentialSourceDifferentiability; // Differentiability of source approximaiton to delta function: range from  1 to 10

    UCLAQ::GridFunction3d                Vtmp;
	UCLAQ::GridFunction3d           potential;
};

#undef  _POTENTIAL_MOMENT_ORDER
#undef  _POTENTIAL_SOURCE_DIFFERENTIABILITY
#undef  _LAPLACE_OP_ORDER

#endif
