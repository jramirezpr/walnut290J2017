#include <iostream>
#include <cmath>
#include <functional>
#include <random>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin1d.h"

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "GridFunctionNd/UCLAQ_GridFunction2dUtility.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin2d.h"

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin3d.h"


int main()
{
	long xPanels  = 10;
	long yPanels  = 8;
	long zPanels  = 9;

	double xMin = -1.0;
	double xMax =  2.0;
	double LX   = (xMax-xMin);

	double yMin = -1.0;
	double yMax =  3.0;
	double LY   = (yMax-yMin);

    double zMin = -3.0;
	double zMax = -2.4;
	double LZ   = (zMax-zMin);

	long kx;   long ky;  long kz;

	double errCoeff;
    double dotCoeff;
    double err;

	double pi = 3.1415926535897932384;

	std::mt19937_64 gen(9265358979);  // to seed mersenne twister.
	std::uniform_real_distribution<> distr(0.0,1.0);


	std::function<double(double)> uRandom1d = [&distr,&gen](double x)
	{
    return  distr(gen);
	};

    std::function<double(double,double)> uRandom2d = [&distr,&gen](double x,double y)
	{
    return  distr(gen);
	};

    std::function<double(double,double,double)> uRandom3d = [&distr,&gen](double x,double y,double z)
	{
    return  distr(gen);
	};

	UCLAQ::GridFunction1d            f(xPanels,xMin,xMax);
	UCLAQ::GridFunction1d        forig(xPanels,xMin,xMax);
	UCLAQ::GridFunction1d       fbasis(xPanels,xMin,xMax);
	f.specify(uRandom1d);
	f.setBoundaryValues(0.0);

	forig = f;


    // Initialize the FFTW3 interface routine

	//
    // Tests coefficients obtained via an inner product with
    // an orthonormal basis w.r.t. the scaled inner product
    // (w.scaledDot(v)) and coefficients obtain by scaling the
	// transform values obtained with FFTW_InterfaceNd sin
	// transforms.
    //

    UCLAQ::fftw3_sin1d DFT;    // Discrete sin transform interface
    DFT.initialize(xPanels,LX);   // Specify the number of grid panels!

    UCLAQ::DoubleVector1d fTransform(xPanels-1);
    fTransform.setToValue(0.0);

    DFT.fftw1d_sin_forward(f,fTransform);

    // Create basis functions and verify coefficient construction

    std::function<double(double)> sin1dbasis= [&kx,pi,xMin,LX] (double x)
	{
    return sqrt(2.0/LX)*sin((double(kx)*pi*(x-xMin))/LX);
	};

    errCoeff = 0.0;
    for(kx = 1; kx <= xPanels-1; kx++)
    {
    fbasis.specify(sin1dbasis);

	dotCoeff =  f.scaledDot(fbasis);
	err      =  abs(dotCoeff - fTransform(kx-1));
	errCoeff = (errCoeff < err) ? err : errCoeff;
    }

    printf("1D: Error in consistency of sin transform coefficient creation: %10.5e \n",errCoeff);

    // Compare inverse transform computation using inner products vs. fftw inverse transform

    DFT.fftw1d_sin_inverse(fTransform,f);

    errCoeff = (f - forig).normInf();
    printf("1D: Error in consistency of fftw inverse transform : %10.5e \n",errCoeff);

    // Evaluation using inner products

    double hx = (xMax-xMin)/(double)xPanels;
    double x;
    for(long i = 1; i <= xPanels-1; i++)
    {
    	x = double(i)*hx;
    	for(long k=1; k < xPanels; k++)
    	{
		fbasis(k-1) = sqrt(2.0/LX)*sin((double(k)*pi*x)/LX);
    	}
    	f(i) = fTransform.dot(fbasis);
    }

    errCoeff = (f - forig).normInf();
    printf("1D: Error in consistency of dot product inverse transform : %10.5e \n",errCoeff);


//   Two Dimensions

	UCLAQ::GridFunction2d            f2(xPanels,xMin,xMax,yPanels,yMin,yMax);
	UCLAQ::GridFunction2d       f2basis(xPanels,xMin,xMax,yPanels,yMin,yMax);
	f2.specify(uRandom2d);
	f2.setBoundaryValues(0.0);

    UCLAQ::fftw3_sin2d DFT2;                  // Discrete sin transform interface
    DFT2.initialize(xPanels,yPanels,LX,LY);   // Specify the number of grid panels and domain size

    UCLAQ::DoubleVector2d f2Transform(xPanels-1,yPanels-1);
    f2Transform.setToValue(0.0);

    DFT2.fftw2d_sin_forward(f2,f2Transform);

    // Create basis functions and verify coefficient construction

    std::function<double(double,double)> sin2dbasis= [&kx,&ky,pi,xMin,yMin,LX,LY] (double x,double y)
	{
    return sqrt(2.0/LX)*sin((double(kx)*pi*(x-xMin))/LX)*sqrt(2.0/LY)*sin((double(ky)*pi*(y-yMin))/LY);
	};

    errCoeff = 0.0;
    for(kx = 1; kx <= xPanels-1; kx++)
    {
    for(ky = 1; ky <= yPanels-1; ky++)
    {
    	f2basis.specify(sin2dbasis);
    	dotCoeff =  f2.scaledDot(f2basis);
    	err      =  abs(dotCoeff - f2Transform(kx-1,ky-1));
    	errCoeff = (errCoeff < err) ? err : errCoeff;
    }}

    printf("2D: Error in consistency of sin transform coefficient creation: %10.5e \n",errCoeff);

    // Compare inverse transform computation using inner products vs. fftw inverse transform

    DFT2.fftw2d_sin_inverse(f2Transform,f2);

    errCoeff = (f - forig).normInf();
    printf("2D: Error in consistency of fftw inverse transform : %10.5e \n",errCoeff);

//
//   ThreeDimensions
//
	UCLAQ::GridFunction3d            f3(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	UCLAQ::GridFunction3d        f3orig(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	UCLAQ::GridFunction3d       f3basis(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	f3.specify(uRandom3d);
	f3.setBoundaryValues(0.0);

    UCLAQ::fftw3_sin3d DFT3;                           // Discrete sin transform interface
    DFT3.initialize(xPanels,yPanels,zPanels,LX,LY,LZ); // Specify the number of grid panels and domain sizes


    UCLAQ::DoubleVector3d f3Transform(xPanels-1,yPanels-1,zPanels-1);
    f3Transform.setToValue(0.0);

    DFT3.fftw3d_sin_forward(f3,f3Transform);

    // Create basis functions and verify coefficient construction

    std::function<double(double,double,double)> sin3dbasis= [&kx,&ky,&kz,pi,xMin,yMin,zMin,LX,LY,LZ] (double x,double y,double z)
	{
    return sqrt(2.0/LX)*sin((double(kx)*pi*(x-xMin))/LX)*sqrt(2.0/LY)*sin((double(ky)*pi*(y-yMin))/LY)*sqrt(2.0/LZ)*sin((double(kz)*pi*(z-zMin))/LZ);
	};

    errCoeff = 0.0;
    for(kx = 1; kx <= xPanels-1; kx++)
    {
    for(ky = 1; ky <= yPanels-1; ky++)
    {
    for(kz = 1; kz <= zPanels-1; kz++)
    {
    	f3basis.specify(sin3dbasis);
    	dotCoeff =  f3.scaledDot(f3basis);
    	err      =  abs(dotCoeff - f3Transform(kx-1,ky-1,kz-1));
    	errCoeff = (errCoeff < err) ? err : errCoeff;
    }}}

    printf("3D: Error in consistency of sin transform coefficient creation: %10.5e \n",errCoeff);

    // Compare inverse transform computation using inner products vs. fftw inverse transform

    DFT3.fftw3d_sin_inverse(f3Transform,f3);

    errCoeff = (f - forig).normInf();
    printf("3D: Error in consistency of fftw inverse transform : %10.5e \n",errCoeff);


}

