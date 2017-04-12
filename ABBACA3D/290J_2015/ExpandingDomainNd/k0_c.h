#include <cmath>
#include <cfloat>
using namespace std;

/* k0.f -- Translated by f2c (version 20100827).
   from http://www.netlib.org/specfun

Modifications:

removed static
doublereal <- double
integer    <- long 
inserted machine constants from cfloat 
(DBL_EPSILON,DBL_MAX)


Chris Anderson (C) UCLA July 8, 2015

This class has no dependencies other than C++ header files 
*/

/*
#############################################################################
#
# Modification of the f2c translation of k0.f by 
# W. J. Cody and Laura Stoltz 
#
# Copyright 2015 Chris Anderson
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
#ifndef _BesselK0_
#define _BesselK0_

class BesselK0
{

public :

/*   This function program computes approximate values for the */
/*   modified Bessel function of the second kind of order zero */
/*   for arguments 0.0 .LT. ARG .LE. XMAX (see comments heading */
/*   Authors: W. J. Cody and Laura Stoltz */
/*   Latest Modification: January 19, 1988 */

double besk0(double x) const
{
    double ret_val;
    long jint;
    double result;

    jint = 1;
    calck0_(&x, &result, &jint);
    ret_val = result;
    return ret_val;
}

/* This function program computes approximate values for the */
/* modified Bessel function of the second kind of order zero */
/* multiplied by the Exponential function, for arguments */
/* 0.0 .LT. ARG. */

/*  Authors: W. J. Cody and Laura Stoltz */
/*  Latest Modification: January 19, 1988 */

double besek0(double x) const
{
    double ret_val;
    long jint;
    double result;
    jint = 2;
    calck0_(&x, &result, &jint);
    ret_val = result;
    return ret_val;
}  

private: 

int calck0_(double *arg, double *result, long *jint) const
{
    /* Initialized data */

    double one = 1.;
    double pp[10] = { 113.94980557384778174,3683.258995734026794,
	    31075.408980684392399,105770.68948034021957,173988.67902565686251,
	    150976.46353289914539,71557.062783764037541,18321.525870183537725,
	    2344.4738764199315021,116.00249425076035558 };
    double qq[10] = { 200.13443064949242491,4432.9628889746408858,
	    31474.655750295278825,97418.829762268075784,151446.44673520157801,
	    126898.39587977598727,58824.616785857027752,14847.228371802360957,
	    1882.1890840982713696,92.556599177304839811 };
    double zero = 0.;
    double xsmall = DBL_EPSILON; // 1.11e-16;
    double xinf   = DBL_MAX;     // 1.79e308;
    double xmax = 705.342;
    double p[6] = { 5.85992214128261e-4,.1316605256498957185,
	    11.999463724910714109,468.50901201934832188,5916.9059852270512312,
	    2470.8152720399552679 };
    double q[2] = { -249.94418972832303646,21312.71430384912038 };
    double f[4] = { -1.64144528372990641,-296.01657892958843866,
	    -17733.784684952985886,-403203.40761145482298 };
    double g[3] = { -250.6497244587799273,29865.713163054025489,
	    -1612813.6304458193998 };

    /* Builtin functions */
    //double log(double), exp(double), sqrt(double);

    /* Local variables */
    long i__;
    double x, xx, temp, sumf, sumg, sump, sumq;

/* -------------------------------------------------------------------- */

/* This packet computes modified Bessel functions of the second kind */
/*   and order zero, K0(X) and EXP(X)*K0(X), for real */
/*   arguments X.  It contains two function type subprograms, BESK0 */
/*   and BESEK0, and one subroutine type subprogram, CALCK0. */
/*   the calling statements for the primary entries are */

/*                   Y=BESK0(X) */
/*   and */
/*                   Y=BESEK0(X) */

/*   where the entry points correspond to the functions K0(X) and */
/*   EXP(X)*K0(X), respectively.  The routine CALCK0 is */
/*   intended for internal packet use only, all computations within */
/*   the packet being concentrated in this routine.  The function */
/*   subprograms invoke CALCK0 with the statement */
/*          CALL CALCK0(ARG,RESULT,JINT) */
/*   where the parameter usage is as follows */

/*      Function                     Parameters for CALCK0 */
/*       Call              ARG                  RESULT          JINT */

/*     BESK0(ARG)   0 .LT. ARG .LE. XMAX       K0(ARG)           1 */
/*     BESEK0(ARG)     0 .LT. ARG           EXP(ARG)*K0(ARG)     2 */

/*   The main computation evaluates slightly modified forms of near */
/*   minimax rational approximations generated by Russon and Blair, */
/*   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461, */
/*   1969.  This transportable program is patterned after the */
/*   machine-dependent FUNPACK packet NATSK0, but cannot match that */
/*   version for efficiency or accuracy.  This version uses rational */
/*   functions that theoretically approximate K-SUB-0(X) to at */
/*   least 18 significant decimal digits.  The accuracy achieved */
/*   depends on the arithmetic system, the compiler, the intrinsic */
/*   functions, and proper selection of the machine-dependent */
/*   constants. */

/* ******************************************************************* */
/* ******************************************************************* */

/* Explanation of machine-dependent constants */

/*   laplaceCoeff   = Radix for the floating-point system */
/*   minexp = Smallest representable power of laplaceCoeff */
/*   maxexp = Smallest power of laplaceCoeff that overflows */
/*   XSMALL = Argument below which BESK0 and BESEK0 may */
/*            each be represented by a constant and a log. */
/*            largest X such that  1.0 + X = 1.0  to machine */
/*            precision. */
/*   XINF   = Largest positive machine number; approximately */
/*            laplaceCoeff**maxexp */
/*   XMAX   = Largest argument acceptable to BESK0;  Solution to */
/*            equation: */
/*               W(X) * (1-1/8X+9/128X**2) = laplaceCoeff**minexp */
/*            where  W(X) = EXP(-X)*SQRT(PI/2X) */


/*     Approximate values for some important machines are: */


/*                           laplaceCoeff       minexp       maxexp */

/*  CRAY-1        (S.P.)       2        -8193         8191 */
/*  Cyber 180/185 */
/*    under NOS   (S.P.)       2         -975         1070 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (S.P.)       2         -126          128 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (D.P.)       2        -1022         1024 */
/*  IBM 3033      (D.P.)      16          -65           63 */
/*  VAX D-Format  (D.P.)       2         -128          127 */
/*  VAX G-Format  (D.P.)       2        -1024         1023 */


/*                          XSMALL       XINF         XMAX */

/* CRAY-1        (S.P.)    3.55E-15   5.45E+2465    5674.858 */
/* Cyber 180/855 */
/*   under NOS   (S.P.)    1.77E-15   1.26E+322      672.788 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (S.P.)    5.95E-8    3.40E+38        85.337 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (D.P.)    1.11D-16   1.79D+308      705.342 */
/* IBM 3033      (D.P.)    1.11D-16   7.23D+75       177.852 */
/* VAX D-Format  (D.P.)    6.95D-18   1.70D+38        86.715 */
/* VAX G-Format  (D.P.)    5.55D-17   8.98D+307      706.728 */

/* ******************************************************************* */
/* ******************************************************************* */

/* Error returns */

/*  The program returns the value XINF for ARG .LE. 0.0, and the */
/*  BESK0 entry returns the value 0.0 for ARG .GT. XMAX. */

/*  Intrinsic functions required are: */

/*     EXP, LOG, SQRT */

/*  Latest modification: March 19, 1990 */

/*  Authors: W. J. Cody and Laura Stoltz */
/*           Mathematics and Computer Science Division */
/*           Argonne National Laboratory */
/*           Argonne, IL 60439 */

    x = *arg;
    if (x > zero) {
	if (x <= one) {
/* -------------------------------------------------------------------- */
/*     0.0 .LT.  ARG  .LE. 1.0 */
/* -------------------------------------------------------------------- */
	    temp = log(x);
	    if (x < xsmall) {
/* -------------------------------------------------------------------- */
/*     Return for small ARG */
/* -------------------------------------------------------------------- */
		*result = p[5] / q[1] - temp;
	    } else {
		xx = x * x;
		sump = ((((p[0] * xx + p[1]) * xx + p[2]) * xx + p[3]) * xx + 
			p[4]) * xx + p[5];
		sumq = (xx + q[0]) * xx + q[1];
		sumf = ((f[0] * xx + f[1]) * xx + f[2]) * xx + f[3];
		sumg = ((xx + g[0]) * xx + g[1]) * xx + g[2];
		*result = sump / sumq - xx * sumf * temp / sumg - temp;
		if (*jint == 2) {
		    *result *= exp(x);
		}
	    }
	} else if (*jint == 1 && x > xmax) {
/* -------------------------------------------------------------------- */
/*     Error return for ARG .GT. XMAX */
/* -------------------------------------------------------------------- */
	    *result = zero;
	} else {
/* -------------------------------------------------------------------- */
/*     1.0 .LT. ARG */
/* -------------------------------------------------------------------- */
	    xx = one / x;
	    sump = pp[0];
	    for (i__ = 2; i__ <= 10; ++i__) {
		sump = sump * xx + pp[i__ - 1];
/* L120: */
	    }
	    sumq = xx;
	    for (i__ = 1; i__ <= 9; ++i__) {
		sumq = (sumq + qq[i__ - 1]) * xx;
/* L140: */
	    }
	    sumq += qq[9];
	    *result = sump / sumq / sqrt(x);
	    if (*jint == 1) {
		*result *= exp(-x);
	    }
	}
    } else {
/* -------------------------------------------------------------------- */
/*     Error return for ARG .LE. 0.0 */
/* -------------------------------------------------------------------- */
	*result = xinf;
    }
/* -------------------------------------------------------------------- */
/*     Update error counts, etc. */
/* -------------------------------------------------------------------- */
    return 0;
/* ---------- Last line of CALCK0 ---------- */
} /* calck0_ */

};

#endif

