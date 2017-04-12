#ifndef  _UCLAQ_PolyFun
#define  _UCLAQ_PolyFun
//
//######################################################################
//                        UCLAQ_PolyFun.h
//######################################################################
//

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

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

namespace UCLAQ
{

class PolyFun
{

    public :

	PolyFun() : coeff()
	{
	    degree = -1;
	    zero   = 0.0;
	}

	PolyFun(const PolyFun& P)
	{
	    degree = P.degree;
	    coeff      = P.coeff;
	    zero   = P.zero;
	}

	PolyFun(int degree)
	{
	    this->degree = degree;
	    coeff.resize(degree+1,0.0);
	    zero = 0.0;
	}

	PolyFun(const vector<double>& coefficients)
	{
	    degree = coefficients.size()-1;
	    coeff      = coefficients;
	    zero   = 0.0;
	}

	~PolyFun()
	{}

	void initialize()
	{
	    coeff.clear();
	    degree = -1;
	}

	void initialize(const PolyFun& P)
	{
	    coeff      = P.coeff;
	    degree = P.degree;
	}

	void initialize(int degree)
	{
	    this->degree = degree;
	    coeff.resize(degree+1,0.0);
	}

	void initialize(const vector<double>& coefficients)
	{
	    degree = coefficients.size()-1;
	    coeff      = coefficients;
	}

	PolyFun integrate()
	{
	    int i;
	    PolyFun R(degree+1);

	    R.coeff[0] = 0.0;                     // pick arbitrary constant
	                                      // to be zero
	    for(i = 1; i <= degree + 1; i++)
	    {R.coeff[i] =  coeff[i-1]/double(i);}

	    return R;
	}

	double integrate(double Xmin, double Xmax)
	{
	    PolyFun R = this->integrate();
	    return (R(Xmax) - R(Xmin));
	}


	PolyFun differentiate()
	{
	    int i;
	    double rDegree;

	    if(degree > 0)
	    {rDegree = degree - 1;}
	    else
	    {rDegree = 0;}

	    PolyFun R(rDegree);

	    if(degree > 0)
	    {
	    for(i = 0; i <= degree - 1; i++)
	    {R.coeff[i] =  coeff[i+1]*double(i+1);}
	    }
	    else
	    {
	    R.coeff[0] = 0.0;
	    }

	    return R;
	}

	void operator=(const PolyFun& P)
	{
	    degree = P.degree;
	    coeff  = P.coeff;
	}

	PolyFun operator+(const PolyFun& P)
	{
	    PolyFun  R;
	    int i;

	    if(degree >= P.degree)
	    {
	      R = *this;
	      for(i = 0; i <= P.degree; i++)
	      {R.coeff[i] += P.coeff[i];}
	    }
	    else
	    {
	       R = P;
	       for(i = 0; i <= degree; i++)
	       {R.coeff[i] += coeff[i];}
	    }
	    return R;
	}

	PolyFun operator-(const PolyFun& P)
	{
	    PolyFun  R;
	    int i;

	    if(degree >= P.degree)
	    {
	      R = *this;
	      for(i = 0; i <= P.degree; i++)
	      {R.coeff[i] -= P.coeff[i];}
	    }
	    else
	    {
	       R  =    P;
	       R *= -1.0;
	       for(i = 0; i <= degree; i++)
	       {R.coeff[i] += coeff[i];}
	    }
	    return R;
	}
	PolyFun operator*(const PolyFun& P)
	{
	    int i; int k;
	    PolyFun R(degree+P.degree);

	    for(i = 0; i <= degree; i++)
	    {
	    for(k = i; k <= i+P.degree; k++)
	    {R.coeff[k] += coeff[i]*P.coeff[k-i];}
	    }

	    return R;
	}

	PolyFun operator-()
	{
	    PolyFun  R(*this);
	    int i;
	    for(i = 0; i <= degree; i++)
	    {R.coeff[i] = -coeff[i];}
	    return R;
	}

	double operator()(double x) const
	{
	//  Evaluate using Horner's method
	//
	    double result = 0.0;
	    int i;
	    for(i = degree; i>=0; i--)
	    {result = x*result + coeff[i];}
	    return result;
	}

#if __cplusplus > 199711L
	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return std::move(F);
	}
#endif

	const double&  operator[](long i) const
	{
	   if(i > degree) {return zero;}
	   return coeff[i];
	}

	double&  operator[](long i)
	{
	   if(i > degree)
	   {degree = i; coeff.resize(degree+1,0.0);}
	   return coeff[i];
	}

	PolyFun operator*(double alpha)
	{
	   PolyFun R(*this);
	   R *= alpha;
	   return R;
	}

	void operator*=(double alpha)
	{
	   for(int i =0; i <= degree; i++) coeff[i] *= alpha;
	}

	void operator/=(double alpha)
	{
	   for(int i =0; i <= degree; i++) coeff[i] /= alpha;
	}


	PolyFun operator/(double alpha)
	{
	   PolyFun R(*this);
	   R /= alpha;
	   return R;
	}

	friend PolyFun operator*(double alpha, const PolyFun& P)
	{
	   PolyFun R(P);
	   R *= alpha;
	   return R;
	}


	// Shift the independent variable x -> x - p

	PolyFun shift(double p)
	{
	     PolyFun R(degree);

	     PolyFun S(1);
	     S[0] = -p; S[1] = 1.0;      // S = (x - p)
	     PolyFun Q(S);

	     R[0] = coeff[0];
	     for(long i = 1; i <= degree; i++)
	     {
	     R = R + coeff[i]*Q;
	     Q = Q*S;
	     }
	     return R;
	}


	// Scale the independent variable x -> alpha*x

	PolyFun scale(double alpha)
	{
	     PolyFun R(*this);
	     double s = 1.0;
	     for(long i = 1; i <= R.degree; i++)
	     {
	     s = s*alpha;
	     R.coeff[i] = R.coeff[i]*s;
	     }
	     return R;
	}

	//

	friend ostream& operator <<(ostream& outStream, const PolyFun& P)
	{
	    outStream.setf(ios::fixed);
	    outStream.precision(2);
	    if(P.degree >= 0)
	    {outStream << setw(2) << P.coeff[0];}
	    int i;
	    if(P.degree >= 1)
	    {
	    for(i = 1; i <= P.degree; i++)
	    {
	    outStream << " + " << setw(2) << P.coeff[i] << "x^";
	    outStream.setf(ios::left);
	    outStream << setw(2) << i;
	    outStream.setf(ios::right);
	    }
	    }
	    return outStream;
	}


	int getDegree() const
	{
	     long i;
	     int deg = 0;
	     for(i = degree; i >= 0; i--)
	     {
	     if((deg == 0)&&(coeff[i] != 0.0)) deg = i;
	     }
	     return deg;
	}

	vector<double> getCoefficients() const
	{
	     return coeff;
	}


    bool operator==(const PolyFun& P)
	{
    	double epsTol     = 1.0e-14;
    	bool returnValue  = true;
    	int minDegree     = (coeff.size() < P.coeff.size()) ? coeff.size()-1 : P.coeff.size() -1;

    	for(int i = 0; i <= minDegree; i++)
    	{
    	if(abs(coeff[i] - P.coeff[i]) > epsTol) returnValue = false;
    	}

    	for(int i = minDegree+1; i < (int)coeff.size(); i++)
    	{
    	if(abs(coeff[i]) > epsTol) returnValue = false;
    	}

    	for(int i = minDegree+1; i < (int)P.coeff.size(); i++)
    	{
    	if(abs(P.coeff[i]) > epsTol) returnValue = false;
    	}

    	return returnValue;
	}

    bool operator!=(const PolyFun& P)
	{
    	return not this->operator==(P);
	}

    private :

    vector<double>     coeff;
    int               degree;
    double              zero;  // reference to 0 coefficient value
};

} // Namespace UCLAQ
#endif
  
