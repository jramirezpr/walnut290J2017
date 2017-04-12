
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
using namespace std;

#ifndef  _UCLAQ_HermiteProductFunction2d
#define  _UCLAQ_HermiteProductFunction2d
//
//######################################################################
//                UCLAQ_HermiteProductFunction2d.h
//######################################################################
//
/*

The UCLAQ::HermiteProductFunction2d creates functions from R^3 -> R that
are products of three Hermite functions.

If gammaX, gammaY, and gammaZ are double values greater than  0 and nx,ny, and ny
are integer values all greater than or equal to zero, then these functions are products
of functions of the form

H_n(s) = C*exp(-((gamma*s)^2)/2)*P_n(s*gamma)

where C = sqrt(gamma/ [sqrt(pi)*(2^n)*n!] ), P_n is the nth
Hermite polynomial, and n >= 0.

When s = x, gamma = gammaX and n = nx
     s = y, gamma = gammaY and n = ny
     s = z, gamma = gammaZ and n = nz


These functions are orthonormal over [-oo, oo]X[-oo, oo]X[-oo, oo].

If alphaX, alphaY, and alphaZ and betaX, betaY and betaZ are real values such that
alphaX*betaX  < 0, alphaY*betaY  < 0, and alphaX*betaZ  < 0,  then these functions
are eigenfunctions of the separable operator

L = alphaX*d^/dx^2 + alphaY*d^/dy^2 + alphaX*d^/dz^2 + betaX*x^2 + betaY*y^2 + betaZ*z^2

when gammaX = (betaX/alphaX)^(1/4), gammaY = (betaY/alphaY)^(1/4) and gammaZ = (betaZ/alphaZ)^(1/4)

The eigenvalue is given by

lambda  = [ -(2*nX + 1)*sign(alphaX)*(|alphaX*betaX|)^(1/2) ] +
          [ -(2*nY + 1)*sign(alphaY)*(|alphaY*betaY|)^(1/2) ] +
          [ -(2*nZ + 1)*sign(alphaZ)*(|alphaZ*betaZ|)^(1/2) ]

e.g. the sum of the eigenvalues of each of the separable components.

By specifying shiftX, shiftY and shiftZ as additional parameters in the
constructor (or initializer) the center of the functions created can be
shifted to arbitrary locations in R^3.

*/

/*
#############################################################################
#
# Copyright  2015 Chris Anderson
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

#include "UCLAQ_PolyFun.h"
#include "UCLAQ_OrthoPoly.h"
#include "UCLAQ_HermiteFunction.h"

namespace UCLAQ
{

class HermiteProductFunction2d
{
    public :

    HermiteProductFunction2d()
	{
    initialize();
	}

    HermiteProductFunction2d(double gammaX, double gammaY)
    {
    initialize(gammaX,gammaY,shiftX,shiftY);
    }

    HermiteProductFunction2d(double gammaX, double gammaY, double shiftX, double shiftY)
    {
    initialize(gammaX,gammaY,shiftX,shiftY);
    }

    virtual ~HermiteProductFunction2d(){};

    void initialize()
	{
    gammaX = gammaY = 1.0;
    shiftX = shiftY = 0.0;
    orthoPoly.initialize(OrthoPoly::Hermite);
	}

    void initialize(double gammaX, double gammaY)
    {
    shiftX = shiftY = 0.0;
    initialize(gammaX,gammaY,shiftX,shiftY);
    }


    void initialize(double gammaX, double gammaY,  double shiftX, double shiftY)
    {
    this->gammaX = gammaX;
    this->shiftX = shiftX;

    this->gammaY = gammaY;
    this->shiftY = shiftY;

    orthoPoly.initialize(OrthoPoly::Hermite);
    }

    std::function< double(double,double) > getHermiteProductFunction2d(long nX,long nY) const
    {
    double sqrtPi  = sqrt(3.141592653589793238);
    PolyFun PX     = orthoPoly.getNthOrthoPolynomial(nX);
    PolyFun PY     = orthoPoly.getNthOrthoPolynomial(nY);

    double normConstant;

    normConstant  = gammaX/(sqrtPi*pow(2.0,nX)*exp(std::lgamma(nX+1)));
    PX *= sqrt(normConstant);

    normConstant  = gammaY/(sqrtPi*pow(2.0,nY)*exp(std::lgamma(nY+1)));
    PY *= sqrt(normConstant);

    std::function< double(double,double) > H = [=](double x,double y)
	{
    	double val = -0.5*(gammaX*(x-shiftX)*gammaX*(x-shiftX) + gammaY*(y-shiftY)*gammaY*(y-shiftY));
    	return exp(val)*PX(gammaX*(x-shiftX))*PY(gammaY*(y-shiftY));
	};

    return H;
    }

    private :

    double     gammaX,gammaY;
    double     shiftX,shiftY;
    OrthoPoly      orthoPoly;
};

}; // namespace UCLAQ

#endif


