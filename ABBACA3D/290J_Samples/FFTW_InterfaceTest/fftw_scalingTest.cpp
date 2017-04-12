#include <iostream>
#include <cmath>
#include <functional>
#include <random>
using namespace std;


#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_1d.h"

#include "GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector2d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_2d.h"

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "FFTW3_InterfaceNd/UCLAQ_fftw3_3d.h"


int main()
{
	long xPanels  = 10;
	long yPanels  =  8;
	long zPanels  =  8;

	double xMin =  -2.0;
	double xMax =   1.0;
	double LX   = (xMax-xMin);

	double yMin = -1.0;
	double yMax =  3.0;
	double LY   = (yMax-yMin);

    double zMin = -3.0;
	double zMax = -2.4;
	double LZ   = (zMax-zMin);

	long      kx;      long ky;      long kz;
	long kxIndex; long kyIndex; long kzIndex;

	double errCoeff;
    double dotRealCoeff; double dotImagCoeff;
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


    // 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

	UCLAQ::GridFunction1d    f1realData(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d    f1imagData(xPanels,xMin,xMax);

    UCLAQ::GridFunction1d    f1realDataOrig(xPanels,xMin,xMax);
    UCLAQ::GridFunction1d    f1imagDataOrig(xPanels,xMin,xMax);

	UCLAQ::GridFunction1d    realbasis1d(xPanels,xMin,xMax);
	UCLAQ::GridFunction1d    imagbasis1d(xPanels,xMin,xMax);

	f1realData.specify(uRandom1d);
	f1imagData.specify(uRandom1d);

	f1realData.enforcePeriodicity();
	f1imagData.enforcePeriodicity();

    f1realDataOrig = f1realData;
	f1imagDataOrig = f1imagData;

    UCLAQ::fftw3_1d DFT;          // Discrete sin transform interface
    DFT.initialize(xPanels,LX);   // Specify the number of grid panels and domain width

    UCLAQ::DoubleVector1d realTransform1Data(xPanels);
    UCLAQ::DoubleVector1d imagTransform1Data(xPanels);

    DFT.fftw1d_forward(f1realData,f1imagData,realTransform1Data,imagTransform1Data);

    // Create basis functions and verify coefficient construction

    std::function<double(double)> realFun1d= [&kx,pi,xMin,LX] (double x)
	{
    return sqrt(1.0/LX)*cos(double(kx)*2.0*pi*(x-xMin)/LX);
	};

    std::function<double(double)> imagFun1d= [&kx,pi,xMin,LX] (double x)
	{
    return sqrt(1.0/LX)*sin(double(kx)*2.0*pi*(x-xMin)/LX);
	};

    errCoeff = 0.0;
    for(kx = -(xPanels/2); kx <= (xPanels-1)/2; kx++)
    {
    realbasis1d.specify(realFun1d);
    imagbasis1d.specify(imagFun1d);

    kxIndex = kx + (xPanels/2);

	dotRealCoeff =  f1realData.scaledDot(realbasis1d) + f1imagData.scaledDot(imagbasis1d);
	err          =  abs(dotRealCoeff - realTransform1Data(kxIndex));
    dotImagCoeff =  -f1realData.scaledDot(imagbasis1d) + f1imagData.scaledDot(realbasis1d);
	err          =  abs(dotImagCoeff - imagTransform1Data(kxIndex));
	errCoeff     = (errCoeff < err) ? err : errCoeff;
    }

    printf("1D: Error in consistency of transform coefficient creation: %10.5e \n",errCoeff);

    // Check Parseval's realation

    double f1norm2 = f1realData.norm2squared() + f1imagData.norm2squared();
    f1norm2        = sqrt(f1norm2);

    double f1transNorm2 = sqrt(imagTransform1Data.dot(imagTransform1Data) + realTransform1Data.dot(realTransform1Data));

    printf("1D: Error in consistency of Parseval's relation : %10.5e \n",abs(f1norm2 - f1transNorm2));

    // Compare inverse transform computation using inner products vs. fftw inverse transform

    DFT.fftw1d_inverse(realTransform1Data,imagTransform1Data,f1realData,f1imagData);

    errCoeff = (f1realDataOrig - f1realData).normInf();
    printf("1D: Error in consistency of real component of fftw inverse transform : %10.5e \n",errCoeff);

    errCoeff = (f1imagDataOrig - f1imagData).normInf();
    printf("1D: Error in consistency of imag component of fftw inverse transform : %10.5e \n",errCoeff);

// 22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222


    UCLAQ::GridFunction2d    f2realData(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d    f2imagData(xPanels,xMin,xMax,yPanels,yMin,yMax);

    UCLAQ::GridFunction2d    f2realDataOrig(xPanels,xMin,xMax,yPanels,yMin,yMax);
    UCLAQ::GridFunction2d    f2imagDataOrig(xPanels,xMin,xMax,yPanels,yMin,yMax);

	UCLAQ::GridFunction2d    realbasis2d(xPanels,xMin,xMax,yPanels,yMin,yMax);
	UCLAQ::GridFunction2d    imagbasis2d(xPanels,xMin,xMax,yPanels,yMin,yMax);

	f2realData.specify(uRandom2d);
	f2imagData.specify(uRandom2d);

	f2realData.enforcePeriodicity();
	f2imagData.enforcePeriodicity();

    f2realDataOrig = f2realData;
	f2imagDataOrig = f2imagData;

    UCLAQ::fftw3_2d DFT2;                     // Discrete  transform interface
    DFT2.initialize(xPanels,yPanels,LX,LY);   // Specify the number of grid panels and domain sizes

    UCLAQ::DoubleVector2d realTransform2Data(xPanels,yPanels);
    UCLAQ::DoubleVector2d imagTransform2Data(xPanels,yPanels);

    DFT2.fftw2d_forward(f2realData,f2imagData,realTransform2Data,imagTransform2Data);

    // Create basis functions and verify coefficient construction

    std::function<double(double,double)> realFun2d= [&kx,&ky,pi,xMin,yMin,LX,LY] (double x, double y)
	{
    return sqrt(1.0/(LX*LY))*(
    		 cos(double(kx)*2.0*pi*(x-xMin)/LX)*cos(double(ky)*2.0*pi*(y-yMin)/LY)
           - sin(double(kx)*2.0*pi*(x-xMin)/LX)*sin(double(ky)*2.0*pi*(y-yMin)/LY)
		     );
	};

    std::function<double(double,double)> imagFun2d= [&kx,&ky,pi,xMin,yMin,LX,LY] (double x, double y)
	{
    return sqrt(1.0/(LX*LY))*(
    		 sin(double(kx)*2.0*pi*(x-xMin)/LX)*cos(double(ky)*2.0*pi*(y-yMin)/LY)
           + cos(double(kx)*2.0*pi*(x-xMin)/LX)*sin(double(ky)*2.0*pi*(y-yMin)/LY)
		    );
	};

    errCoeff = 0.0;
    for(kx = -(xPanels/2); kx <= (xPanels-1)/2; kx++)
    {
    for(ky = -(yPanels/2); ky <= (yPanels-1)/2; ky++)
    {
    realbasis2d.specify(realFun2d);
    imagbasis2d.specify(imagFun2d);

    kxIndex = kx + (xPanels/2);
    kyIndex = ky + (yPanels/2);

	dotRealCoeff =  f2realData.scaledDot(realbasis2d) + f2imagData.scaledDot(imagbasis2d);
	err          =  abs(dotRealCoeff - realTransform2Data(kxIndex,kyIndex));
    dotImagCoeff =  -f2realData.scaledDot(imagbasis2d) + f2imagData.scaledDot(realbasis2d);
	err          =  abs(dotImagCoeff - imagTransform2Data(kxIndex,kyIndex));
	errCoeff     = (errCoeff < err) ? err : errCoeff;
    }}

    printf("2D: Error in consistency of transform coefficient creation: %10.5e \n",errCoeff);

    // Check Parseval's realation

    double f2norm2 = f2realData.norm2squared() + f2imagData.norm2squared();
    f2norm2        = sqrt(f2norm2);

    double f2transNorm2 = sqrt(imagTransform2Data.dot(imagTransform2Data) + realTransform2Data.dot(realTransform2Data));

    printf("2D: Error in consistency of Parseval's relation : %10.5e \n",abs(f2norm2 - f2transNorm2));

    // Compare inverse transform computation using inner products vs. fftw inverse transform

    DFT2.fftw2d_inverse(realTransform2Data,imagTransform2Data,f2realData,f2imagData);


    errCoeff = (f2realDataOrig - f2realData).normInf();
    printf("2D: Error in consistency of real component of fftw inverse transform : %10.5e \n",errCoeff);

    errCoeff = (f2imagDataOrig - f2imagData).normInf();
    printf("2D: Error in consistency of imag component of fftw inverse transform : %10.5e \n",errCoeff);


    // 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333


    UCLAQ::GridFunction3d    f3realData(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d    f3imagData(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

    UCLAQ::GridFunction3d    f3realDataOrig(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
    UCLAQ::GridFunction3d    f3imagDataOrig(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

	UCLAQ::GridFunction3d    realbasis3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	UCLAQ::GridFunction3d    imagbasis3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

	f3realData.specify(uRandom3d);
	f3imagData.specify(uRandom3d);

	f3realData.enforcePeriodicity();
	f3imagData.enforcePeriodicity();

    f3realDataOrig = f3realData;
	f3imagDataOrig = f3imagData;

    UCLAQ::fftw3_3d DFT3;                                 // Discrete transform interface
    DFT3.initialize(xPanels,yPanels,zPanels,LX,LY,LZ);    // Specify the number of grid panels and domain sizes

    UCLAQ::DoubleVector3d realTransform3Data(xPanels,yPanels,zPanels);
    UCLAQ::DoubleVector3d imagTransform3Data(xPanels,yPanels,zPanels);

    DFT3.fftw3d_forward(f3realData,f3imagData,realTransform3Data,imagTransform3Data);


    // Create basis functions and verify coefficient construction

    std::function<double(double,double,double)> realFun3d= [&kx,&ky,&kz,pi,xMin,yMin,zMin,LX,LY,LZ] (double x, double y,double z)
	{
    return sqrt(1.0/(LX*LY*LZ))*(
    		   cos(double(kx)*2.0*pi*(x-xMin)/LX)*cos(double(ky)*2.0*pi*(y-yMin)/LY)*cos(double(kz)*2.0*pi*(z-zMin)/LZ)
			  -cos(double(kx)*2.0*pi*(x-xMin)/LX)*sin(double(ky)*2.0*pi*(y-yMin)/LY)*sin(double(kz)*2.0*pi*(z-zMin)/LZ)
	          -sin(double(kx)*2.0*pi*(x-xMin)/LX)*cos(double(ky)*2.0*pi*(y-yMin)/LY)*sin(double(kz)*2.0*pi*(z-zMin)/LZ)
			  -sin(double(kx)*2.0*pi*(x-xMin)/LX)*sin(double(ky)*2.0*pi*(y-yMin)/LY)*cos(double(kz)*2.0*pi*(z-zMin)/LZ)
		     );
	};

    std::function<double(double,double,double)> imagFun3d= [&kx,&ky,&kz,pi,xMin,yMin,zMin,LX,LY,LZ] (double x, double y,double z)
	{
    return sqrt(1.0/(LX*LY*LZ))*(
    		  sin(double(kx)*2.0*pi*(x-xMin)/LX)*cos(double(ky)*2.0*pi*(y-yMin)/LY)*cos(double(kz)*2.0*pi*(z-zMin)/LZ)
			 +cos(double(kx)*2.0*pi*(x-xMin)/LX)*cos(double(ky)*2.0*pi*(y-yMin)/LY)*sin(double(kz)*2.0*pi*(z-zMin)/LZ)
             +cos(double(kx)*2.0*pi*(x-xMin)/LX)*sin(double(ky)*2.0*pi*(y-yMin)/LY)*cos(double(kz)*2.0*pi*(z-zMin)/LZ)
			 -sin(double(kx)*2.0*pi*(x-xMin)/LX)*sin(double(ky)*2.0*pi*(y-yMin)/LY)*sin(double(kz)*2.0*pi*(z-zMin)/LZ)
		    );
	};

    errCoeff = 0.0;
    for(kx = -(xPanels/2); kx <= (xPanels-1)/2; kx++)
    {
    for(ky = -(yPanels/2); ky <= (yPanels-1)/2; ky++)
    {
    for(kz = -(zPanels/2); kz <= (zPanels-1)/2; kz++)
    {
    realbasis3d.specify(realFun3d);
    imagbasis3d.specify(imagFun3d);

    kxIndex = kx + (xPanels/2);
    kyIndex = ky + (yPanels/2);
    kzIndex = kz + (zPanels/2);

	dotRealCoeff =  f3realData.scaledDot(realbasis3d) + f3imagData.scaledDot(imagbasis3d);
	err          =  abs(dotRealCoeff - realTransform3Data(kxIndex,kyIndex,kzIndex));
    dotImagCoeff =  -f3realData.scaledDot(imagbasis3d) + f3imagData.scaledDot(realbasis3d);
	err          =  abs(dotImagCoeff - imagTransform3Data(kxIndex,kyIndex,kzIndex));
	errCoeff     = (errCoeff < err) ? err : errCoeff;
    }}}

    printf("3D: Error in consistency of transform coefficient creation: %10.5e \n",errCoeff);

    // Check Parseval's realation

    double f3norm2 = f3realData.norm2squared() + f3imagData.norm2squared();
    f3norm2        = sqrt(f3norm2);

    double f3transNorm2 = sqrt(imagTransform3Data.dot(imagTransform3Data) + realTransform3Data.dot(realTransform3Data));

    printf("3D: Error in consistency of Parseval's relation : %10.5e \n",abs(f3norm2 - f3transNorm2));

    // Compare inverse transform computation using inner products vs. fftw inverse transform

    DFT3.fftw3d_inverse(realTransform3Data,imagTransform3Data,f3realData,f3imagData);


    errCoeff = (f3realDataOrig - f3realData).normInf();
    printf("3D: Error in consistency of real component of fftw inverse transform : %10.5e \n",errCoeff);

    errCoeff = (f3imagDataOrig - f3imagData).normInf();
    printf("3D: Error in consistency of imag component of fftw inverse transform : %10.5e \n",errCoeff);
}

