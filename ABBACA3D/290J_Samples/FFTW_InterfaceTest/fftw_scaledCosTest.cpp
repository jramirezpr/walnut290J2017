#include <iostream>
#include <cmath>
#include <functional>
#include <random>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_cos1d.h"

int main()
{
	long xPanels  = 9;

	double xMin = -1.0;
	double xMax =  2.0;
	double LX   = (xMax-xMin);

	double errCoeff;
    double dotCoeff;
    double err;

	double pi = 3.1415926535897932384;

	std::mt19937_64 gen(9265358979);  // to seed mersenne twister.
	std::uniform_real_distribution<> distr(0.0,1.0);


	std::function<double(double)> uRandom1d = [&distr,&gen](double x)
	{return  distr(gen);};

	UCLAQ::GridFunction1d            f(xPanels,xMin,xMax);
	UCLAQ::GridFunction1d        forig(xPanels,xMin,xMax);
	UCLAQ::GridFunction1d       fbasis(xPanels,xMin,xMax);
	f.specify(uRandom1d);
	f.enforcePeriodicity();

	forig = f;

	//
    // Tests coefficients obtained via an inner product with
    // an orthonormal basis w.r.t. the scaled inner product
    // (w.scaledDot(v)) and coefficients obtain by scaling the
	// transform values obtained with FFTW_InterfaceNd cos
	// transforms.
    //

    UCLAQ::fftw3_cos1d DFT;       // Discrete cos transform interface
    DFT.initialize(xPanels,LX);   // Specify the number of grid panels and domain size

    UCLAQ::DoubleVector1d fTransform(xPanels+1);
    fTransform.setToValue(0.0);

    DFT.fftw1d_cos_forward(f,fTransform);

    // Create basis functions and verify coefficient construction

	long  kx;

    std::function<double(double)> cos1dbasis= [&kx,pi,xMin,LX,xPanels] (double x)
	{
    double val = sqrt(2.0/LX)*cos((double(kx)*pi*(x-xMin))/LX);
    if(kx == 0)       return val/sqrt(2.0);
    if(kx == xPanels) return val/sqrt(2.0);
    return val;
	};

    errCoeff = 0.0;
    for(kx = 0; kx <= xPanels; kx++)
    {
    fbasis.specify(cos1dbasis);

	dotCoeff =  f.scaledDot(fbasis);
	err      =  abs(dotCoeff - fTransform(kx));
	errCoeff = (errCoeff < err) ? err : errCoeff;
    }

    printf("1D: Error in consistency of cos transform coefficient creation: %10.5e \n",errCoeff);


    // Check orthonormality

    UCLAQ::GridFunction1d       fbasisA(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d       fbasisB(xPanels,xMin,xMax);

    double orthoErr = 0.0;
    for(long ka = 0; ka <= xPanels; ka++)
    {
    for(long kb = 0; kb <= xPanels; kb++)
    {
    kx = ka;
    fbasisA.specify(cos1dbasis);

    kx = kb;
    fbasisB.specify(cos1dbasis);

	dotCoeff =  fbasisA.scaledDot(fbasisB);
	if(ka == kb) { err =  abs(dotCoeff - 1.0);}
	else         { err =  abs(dotCoeff);}
	orthoErr = (orthoErr < err) ? err : orthoErr;
    }}

    printf("1D: Error in orthonormality : %10.5e \n",orthoErr);



    // Compare inverse transform computation using inner products vs. fftw inverse transform

    DFT.fftw1d_cos_inverse(fTransform,f);

    errCoeff = (f - forig).normInf();
    printf("1D: Error in consistency of fftw inverse transform : %10.5e \n",errCoeff);




}

