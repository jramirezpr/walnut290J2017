#include <iostream>
#include <cmath>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"

/**
 *                       Class LaplaceOp1d
 *  A class whose apply member function implements the 3-point finite difference
 *  discrete Laplace operator with Dirichlet boundary conditions.
 *
 *  The apply operator of this class applies the difference operator
 *
 *          alpha*[ (D+D-)_x ]
 *
 *  all interior values and sets the boundary values to 0.0
 */

#ifndef _LaplaceOp1d_
#define _LaplaceOp1d_

class LaplaceOp1d
{
public:

	LaplaceOp1d()
	{
	this->alpha = 0.0;
	}

	virtual ~LaplaceOp1d()
	{}

	LaplaceOp1d(double alpha)
	{
	this->alpha   = alpha;
	}

	void initialize(double alpha)
	{
    this->alpha   = alpha;
	}
void initialize(){
alpha=0.0;
}
	/**
	This routine applies alpha times the 3-point discrete Laplace operator to
	the interior grid points associated with a uniform discretization
	of a 1D  domain.

	Input  : UCLAQ_GriddFunction1d class instance whose values are those of a uniform 1D grid.
	         The function values associated with this class are the values at both the interior
	         and boundary points of the discretization.

    Output : The interior values of the input GridFunction are overwritten
             with the 3-point difference approximation boundary values are
             set to 0.

    If _DEBUG is defined at compile time, bounds checking is performed.
	*/

	void apply(UCLAQ::GridFunction1d& V)
	{

	// Capture values since we are over-writing the input vector
    //
	// Note: The use of initialize instead of = here is because the initialize
	// member function is "smart" in the sense that if the existing instance
	// of Vtmp has identical size to the input V, then it just copies
	// over the values, and doesn't destroy and recreate a new instance
	// of the GridFunction2d.
	//

	Vtmp.initialize(V);

    // Extract grid information from V

    double hx = V.getHx();
    long  xPanels = V.getXpanelCount();

	// Interior grid points not adjacent to the edge

	for(long i = 1; i < xPanels; i++)
	{
	V(i) =  alpha*((Vtmp(i+1) - 2.0*Vtmp(i) + Vtmp(i-1))/(hx*hx));
	}

    V.setBoundaryValues(0.0);

	}

    double alpha;               // Coefficient of the discrete Laplace operator
	UCLAQ::GridFunction1d Vtmp; // Temporary
};

#endif



