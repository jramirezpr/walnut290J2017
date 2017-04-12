/*
 * UCLAQ_OrthoPolyUtil.h
 *
 *  Created on: Dec 13, 2015
 *      Author: anderson
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

#ifndef _UCLAQ_OrthoPolyUtility
#define _UCLAQ_OrthoPolyUtility

#include <iostream>
#include <vector>
using namespace std;
//
//######################################################################
//                        UCLAQ_OrthoPolyUtility.h
//######################################################################
//
namespace UCLAQ
{

#include "UCLAQ_OrthoPoly.h"

class OrthoPolyUtility
{
    public :

void getNodesAndWeights(int polyType,long n, double alpha, vector<double>& x, vector<double>& wt) const
{
	switch(polyType)
    {
    case OrthoPoly::Legendre : getLegendreNodesAndWeights(n, x, wt);        break;
    case OrthoPoly::Laguerre : getLaguerreNodesAndWeights(n, x, wt, alpha); break;
    case OrthoPoly::Chebyshev: getChebyshevNodesAndWeights(n, x, wt);       break;
    case OrthoPoly::Hermite  : getHermiteNodesAndWeights(n, x, wt,alpha);   break;
    default: break;
    }

}
void getLegendreNodesAndWeights(long n, vector<double>& x, vector<double>& wt) const
{
//
// This routine generates n Gaussian quadrature weights for integrals of the
// form
//
//     --  +1
//     |
//     |
//     |    f(x) dx
//     |
//    --  -1
//
//        x  = nodes     vector<double>
//        wt = weights   vector<double>
//

    x.clear();
	wt.clear();

	x.resize(n,0.0);
	wt.resize(n,0.0);

    double eps   = 1.0e-14;
    double pi    = 3.141592653589793238;
	long   i,j,m;
	double z1,z,pp,p3,p2,p1;

    m=(n+1)/2;
	for (i=1; i<= m;i++)
    {
		z=cos(pi*(i-0.25)/(n+0.5));
		do
        {
        p1=1.0;
        p2=0.0;
        for (j=1; j<=n; j++)
        {
		p3=p2;
		p2=p1;
		p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
		}
		pp=n*(z*p1-p2)/(z*z-1.0);
		z1=z;
		z=z1-p1/pp;
		} while (abs(z-z1) > eps);

		x[i-1]=-z;
		x[n-i]= z;
		wt[i-1]=2.0/((1.0-z*z)*pp*pp);
		wt[n-i]=wt[i-1];
	}
}




void getLaguerreNodesAndWeights(long n, vector<double>& x, vector<double>& wt,double alpha) const
{
//
// This routine generates n Gaussian quadrature weights for integrals of the
// form
//
//     --  +oo
//     |
//     |        alpha  -x
//     |     x         e                      f(x) dx
//     |
//    --  0
//
//    x  = nodes     vector<double>
//    wt = weights   vector<double>
//
//    alpha = x exponent scale factor (double)
//
	x.clear();
	wt.clear();

	x.resize(n,0.0);
	wt.resize(n,0.0);

    double eps    = 1.0e-14;
    long   maxit  = 12;

	long i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++)
     {
		if (i == 1)
          {
			z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*n+1.8*alpha);
		}
          else if (i == 2)
          {
			z += (15.0+6.25*alpha)/(1.0+0.9*alpha+2.5*n);
		}
          else
          {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alpha/
				(1.0+3.5*ai))*(z-x[i-3])/(1.0+0.3*alpha);
		}
		for (its=1;its<=maxit;its++)
          {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++)
                {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alpha-z)*p2-(j-1+alpha)*p3)/j;
			}
			pp=(n*p1-(n+alpha)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= eps) break;
		}
		if (its > maxit)
          cerr << "Error in Computation of Laguerre Weights" << endl;
		x[i-1] =  z;
		wt[i-1] = -exp(lgamma(alpha+n)- lgamma(double(n)))/(pp*n*p2);
	}
}


void getChebyshevNodesAndWeights(long n, vector<double>& x, vector<double>& wt) const
{
//
// This routine generates n Gaussian quadrature weights for integrals of the
// form
//
//      --  +1
//      |
//      |
//      |     f(x)/(1.0 - x^2) dx
//      |
//     --  -1
//
//        x  = nodes     vector<double>
//        wt = weights   vector<double>
//

	x.clear();
	wt.clear();

	x.resize(n,0.0);
	wt.resize(n,0.0);

	double pi  = 3.141592653589793238;
	long     j;
	double theta;
	for(j = 0; j < n; j++)
	{
		theta = ((double(j+1) - 0.5)/double(n))*pi;
		x[j]   = cos(theta);
		wt[j]  = pi/double(n);
	}
}

void getHermiteNodesAndWeights(long n, vector<double>& x, vector<double>& wt, double alpha) const
{

//
// This routine generates n Gaussian quadrature weights for integrals of the
// form
//
//      --  +oo
//      |
//      |            -alpha^2*x^2
//      |      e                             f(x) dx
//      |
//     --  -oo
//
//        x  = nodes     vector<double>
//        wt = weights   vector<double>
//
// alpha = exponential scale factor (double)
//
    x.clear();
	wt.clear();

	x.resize(n,0.0);
	wt.resize(n,0.0);

    double eps         = 1.0e-14;
    double pi          = 3.141592653589793238;
    double pifactr     = sqrt(sqrt(1.0/pi));
    long    maxit        = 12;
    int     i,its,j,m;
	double p1,p2,p3,pp,z,z1;

	m=(n+1)/2;
	for (i=1;i<=m;i++)
    {
    //
    // initial guesses for roots
    //
    switch(i)
    {
	case  1 : z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667); break;
	case  2 : z -= 1.14*pow((double)n,0.426)/z; break;
	case  3 : z =  1.86*z-0.86*x[0]; break;
	case  4 : z =  1.91*z-0.91*x[1]; break;
	default : z =  2.0*z-x[i-3]; break;
    }
    //
    // refine using Newton interation
    //
	for (its=1;its<=maxit;its++)
    {
	p1=pifactr;
	p2=0.0;
	for (j=1;j<=n;j++)
	{
		p3=p2;
		p2=p1;
		p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
	}
	pp=sqrt((double)2*n)*p2;
	z1=z;
	z=z1-p1/pp;
	if (fabs(z-z1) <= eps) break;
	}
	if (its > maxit)
    cerr << "Error in Computation of Hermite Weights" << endl;
	x[i-1]=z;
	x[n-i] = -z;
	wt[i-1]=2.0/(pp*pp);
	wt[n-i]=wt[i-1];
	}
//
// scale for exponent
//
	for(i = 0; i < n; i++)
    {
    wt[i] = wt[i]/alpha;
    x[i]  = x[i]/alpha;
    }
}


};
} // Namespace UCLAQ
#endif 
