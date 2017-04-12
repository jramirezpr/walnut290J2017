
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
using namespace std;

#ifndef  _UCLAQ_HermiteFunction
#define  _UCLAQ_HermiteFunction
//
//######################################################################
//                    UCLAQ_HermiteFunction.h
//######################################################################
//
/*

Instances of UCLAQ::HermiteFunction are the functions

H_n(x) = C*exp(-((gamma*x)^2)/2)*P_n(x*gamma)

where C = sqrt(gamma/ [sqrt(pi)*(2^n)*n!] ), P_n is the nth
Hermite polynomial, and n >= 0.

These functions are orthonormal over (-oo,+oo).

If alpha and beta are real values such that alpha*beta < 0, then
H_n(x) are eigenfunctions of the operator

L = alpha*d^2/dx^2 + beta*x^2

when gamma = (beta/alpha)^(1/4). The eigenvalue is given by

lambda  = -(2*n + 1)*sign(alpha)*(|alpha*beta|)^(1/2)

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
#include "UCLAQ_OrthoPolyUtility.h"

namespace UCLAQ
{

class HermiteFunction
{
    public :

//
//  Constructors
//
    HermiteFunction()
	{
    gamma         = 0.0;
	shift         = 0.0;
    orthoPoly.initialize(OrthoPoly::Hermite);
	}

    HermiteFunction(double gamma, double shift = 0.0)
    {
    this->gamma          = gamma;
    this->shift          = shift;
    orthoPoly.initialize(OrthoPoly::Hermite,1.0,shift);
    }

    virtual ~HermiteFunction(){};

    void initialize()
    {
    gamma         = 0.0;
	shift         = 0.0;
    orthoPoly.initialize(OrthoPoly::Hermite);
    }

    void initialize(double gamma, double shift = 0.0)
    {
    this->gamma          = gamma;
    this->shift          = shift;
    orthoPoly.initialize(OrthoPoly::Hermite);
    }

    void   setGamma(double gamma) {this->gamma = gamma;}
    double getGamma()             {return this->gamma;}

    void   setShift(double shift) {this->shift = shift;}
    double getshift()             {return shift;};


    std::function< double(double) > getNthHermiteFunction(long n) const
    {
    double sqrtPi = sqrt(3.141592653589793238);
    PolyFun P     = orthoPoly.getNthOrthoPolynomial(n);

    double normConstant = gamma/(sqrtPi*pow(2.0,n)*exp(std::lgamma(n+1)));
    P  *= sqrt(normConstant);

    double gamma  = this->gamma;
    double shift  = this->shift;

    std::function< double(double) > Hn = [P, gamma, shift](double x)
	{
    	return exp(-0.5*gamma*(x-shift)*gamma*(x-shift))*P(gamma*(x-shift));
	};

    return Hn;
    }


	vector< std::function< double(double) > >  getHermiteFunctionArray(long maxIndex)
    {
	vector< std::function< double(double) > > H;
    for(int i = 0; i < maxIndex; i++)
    {
    	H.push_back(getNthHermiteFunction(i));
    }

    return H;
    }

    private :

    double                    gamma;
    double                    shift;
    OrthoPoly             orthoPoly;
    OrthoPolyUtility   orthoUtility;
};

}; // namespace UCLAQ

#endif


