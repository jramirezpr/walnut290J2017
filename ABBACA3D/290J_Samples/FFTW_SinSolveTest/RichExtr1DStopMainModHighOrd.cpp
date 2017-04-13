#include <iostream>
#include <cmath>
#include <functional>
#include<fstream>
#include<string>
using namespace std;

#ifdef _FFTW_OPENMP
#include <omp.h>
#endif
//Remote changes to richardson
//
//RichardsonExtr.cpp
//
// Solve Schroedingers equation with  RICHARDSON APPLIED TOabbacabba procedure.
//

// This program depends on an installation of the FFTW3 libraries. You may need to modify
// the include paths, library paths, and library to link to an specific installation of
// FFTW3. This program was created an tested on an Ubuntu linux system in which the FFTW3
// libraries were installed using the synaptic package manager.
//
// The test program assumes the directories 290J_2015 and 290J_Samples are organized as follows
//
// Course Directory
// ---------------------------------------------------------------
//      |                |                |            |
// Project_1 Project_2   ....        290J_Samples   290J_2015
//
//
// where the directory 290J_2015 contains the DoubleVectorNd, GridFunctionNd and 
// FFTW3_InterfaceNd source directories.
//
// The command line compilation command for single threaded execution is
//
//
// g++ RichExtr1DStopMainModHighOrd.cpp -std=c++11 -I../../290J_2015 -lfftw3 -o Rich1DTestTemplateRecursion.exe

//

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"

#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "RandOpNd/RandOpDirichlet1d.h"
#include "ShroedingOp3D1Par.cpp"

#include "DoubleVectorNd/UCLAQ_DoubleVector1d.h"
#include "RandOpNd/RandOpDirichlet1d.h"
#include "ShroedingOp1D.cpp"

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin3d.h"
#include "FFTW3_InterfaceNd/UCLAQ_FFT_Nvalues.h"
#include "MollifierNd/HighOrderPolyPotential3d.h"           // Mollified nuclear potential dont need it for 1D
#include "RayleighChebyshev/RayleighChebyshev.h"

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin1d.h"
//#include "MollifierNd/HighOrderPolyPotential1D.h"//     Does not exist for 1D
#include "MollifierNd/SmoothPolyMollifier1d.h"           // Try and see if this one works
 const double pi = 3.141592653589793238;

void heatEqStep1D(double tau, UCLAQ::fftw3_sin1d& DFT, UCLAQ::GridFunction1d&  uIn, UCLAQ::GridFunction1d&  uOut, double alpha = -1.0){
	//1D heatequation step
	double minuseigkx = 0;
	int xPanels = uIn.xPanels;
	double LX = uIn.xMax - uIn.xMin;
	UCLAQ::DoubleVector1d fTransform(xPanels - 1);
	fTransform.setToValue(0.0);
	DFT.fftw1d_sin_forward(uIn, fTransform); // Note: We can use a GridFunction1d as an argument since it's extended
	for (long kx = 1; kx <= xPanels - 1; kx++){
		minuseigkx = tau*alpha*(((kx*pi) / LX)*((kx*pi) / LX));
		fTransform(kx - 1) *= exp(minuseigkx);// its exp(-tau*alpha*eig_kx,ky)
	}
	DFT.fftw1d_sin_inverse(fTransform, uOut);      // Note: We can use a GridFunction3d as an argument since it's extended
}

void heatEqStep(double tau, int xPanels, int yPanels,int zPanels, UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,
	UCLAQ::GridFunction3d&  uOut,double LX,double LY,double LZ,double alpha  = -1.0){
	double threshold=exp(20);
	double minuseigkxyz=0;
	UCLAQ::DoubleVector3d fTransform(xPanels-1,yPanels-1,zPanels-1);
    fTransform.setToValue(0.0);
	long sizegrid=(xPanels-1)*(yPanels-1)*(zPanels-1);
    DFT.fftw3d_sin_forward(uIn,fTransform); // Note: We can use a GridFunction3d as an argument since it's extended
                                            // from a DoubleVector3d
    for(long kx = 1; kx <= xPanels-1; kx++){
		for(long ky = 1; ky <= yPanels-1; ky++){
			for(long kz = 1; kz <= zPanels-1; kz++){
				minuseigkxyz=tau*alpha*( (((kx*pi)/LX)*((kx*pi)/LX)) +  (((ky*pi)/LY)*((ky*pi)/LY)) +  (((kz*pi)/LZ)*((kz*pi)/LZ)) );
				fTransform(kx-1,ky-1,kz-1) *=exp(tau*alpha*( (((kx*pi)/LX)*((kx*pi)/LX)) +  (((ky*pi)/LY)*((ky*pi)/LY)) +  (((kz*pi)/LZ)*((kz*pi)/LZ)) ));
				// its exp(-tau*alpha*eig_kx,ky)
			}
		}
	}
    DFT.fftw3d_sin_inverse(fTransform,uOut);	// Note: We can use a GridFunction3d as an argument since it's extended
												// from a DoubleVector3d
}
void PotentialStep1D(double tau, UCLAQ::GridFunction1d&  uIn, UCLAQ::GridFunction1d&  uOut, UCLAQ::GridFunction1d& potential){
	for (long kx = 1; kx <= potential.xPanels - 1; kx++){
		uOut.Values(kx) = (uIn.Values(kx))* exp(-tau*potential.Values(kx));
	}
}

void PotentialStep(double tau, UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential){
	double max=0;
	int xm=-1,ym=-1,zm=-1;
	for(long kx = 1; kx <= potential.xPanels-1; kx++){
		for(long ky = 1; ky <= potential.yPanels-1; ky++){
			for(long kz = 1; kz <= potential.zPanels-1; kz++){
				uOut.Values(kx,ky,kz) =(uIn.Values(kx,ky,kz))* exp(-tau*potential.Values(kx,ky,kz));
				/*
				if (-tau*potential.Values(kx,ky,kz)>max){
				max=-tau*potential.Values(kx,ky,kz);
				xm=kx; ym=ky;zm=kz;
				*/
			}
		}
	}
	//cout<<"POTENTIAL MAX:"<<max<<" exp  :"<< exp(max)<<" Type 0 to continue \n";
	//cin>>max;
}
/*

//1D version
void ABBAstepwithPROJ1D(double tau,UCLAQ::fftw3_sin1d& DFT,UCLAQ::GridFunction1d&  uIn,    UCLAQ::GridFunction1d&  uOut,UCLAQ::GridFunction1d& potential){
 	int xPanels = uIn.xPanels;
	//AB step
	UCLAQ::GridFunction1d  AUX(xPanels,uIn.xMin,uIn.xMax); 
	AUX.setToValue(0.0);
	heatEqStep1D(tau/2.0,DFT,uIn,AUX);
	UCLAQ::GridFunction1d  AUX2(xPanels,uIn.xMin,uIn.xMax); 
	AUX2.setToValue(0.0);
	PotentialStep1D(tau,AUX, AUX2,potential);
	//BA step
	heatEqStep1D(tau/2.0,DFT,AUX2,uOut);
	//Cstep
	uOut/=uOut.norm2();
}
/////////////////////////

void ABBAstepwithPROJ(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax,
int zPanels,double zMin,double zMax, UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,
UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
	//AB step
	UCLAQ::GridFunction3d  AUX(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX.setToValue(0.0);
	heatEqStep(tau/2.0,xPanels,yPanels,zPanels,DFT,uIn,AUX,LX,LY,LZ);
	UCLAQ::GridFunction3d  AUX2(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX2.setToValue(0.0);
	PotentialStep(tau,AUX, AUX2,potential);
	//BA step
	heatEqStep(tau/2.0,xPanels,yPanels,zPanels,DFT,AUX2,uOut,LX,LY,LZ);
	//Cstep
	uOut/=uOut.norm2();
}
//--------------------------------------
//--------------------------------------

void ABBAstep1D(double tau, UCLAQ::fftw3_sin1d& DFT,UCLAQ::GridFunction1d&  uIn,    UCLAQ::GridFunction1d&  uOut,UCLAQ::GridFunction1d& potential){
 	int xPanels = uIn.xPanels;
	double LX = uIn.xMax - uIn.xMin;
	double xMin=uIn.xMin;
	double xMax=uIn.xMax;
	//AB step
	UCLAQ::GridFunction1d  AUX(xPanels,xMin,xMax); 
	AUX.setToValue(0.0);
	heatEqStep1D(tau/2.0,DFT,uIn,AUX);
	UCLAQ::GridFunction1d  AUX2(xPanels,xMin,xMax); 
	AUX2.setToValue(0.0);
	PotentialStep1D(tau,AUX, AUX2,potential);
	//BA step
	heatEqStep1D(tau/2.0,DFT,AUX2,uOut);
}
void Richardson1D(double tau,UCLAQ::fftw3_sin1d& DFT,UCLAQ::GridFunction1d&  uIn,    UCLAQ::GridFunction1d&  uOut,UCLAQ::GridFunction1d& potential){
 	int xPanels = uIn.xPanels;
	double LX = uIn.xMax - uIn.xMin;
	double xMin=uIn.xMin;
	double xMax=uIn.xMax;
	UCLAQ::GridFunction1d  AUX(xPanels,xMin,xMax); 
	AUX.setToValue(0.0);
	UCLAQ::GridFunction1d  AUX2(xPanels,xMin,xMax); 
	AUX.setToValue(0.0);
	uOut.setToValue(0.0);
	double alph=-1.0/7.0;
	double bet=8.0/7.0;
	ABBAstep1D(tau,DFT,uIn,AUX,potential);
	//two half steps
	ABBAstep1D(tau/2.0,DFT,uIn,AUX2,potential);
	ABBAstep1D(tau/2.0,DFT,AUX2,uOut,potential);
	uOut*=bet;
	AUX*=alph;
	uOut+=AUX;
	uOut/=uOut.norm2();
};
void Richardsonnoproj1D(double tau,UCLAQ::fftw3_sin1d& DFT,UCLAQ::GridFunction1d&  uIn,UCLAQ::GridFunction1d&  uOut,UCLAQ::GridFunction1d& potential){
 	int xPanels = uIn.xPanels;
	double LX = uIn.xMax - uIn.xMin;
	double xMin=uIn.xMin;
	double xMax=uIn.xMax;
	UCLAQ::GridFunction1d  AUX(xPanels,xMin,xMax); 
	AUX.setToValue(0.0);
	UCLAQ::GridFunction1d  AUX2(xPanels,xMin,xMax); 
	AUX.setToValue(0.0);
	uOut.setToValue(0.0);
	double alph=-1.0/7.0;
	double bet=8.0/7.0;
	ABBAstep1D(tau,DFT,uIn,AUX,potential);
	//two half steps
	ABBAstep1D(tau/2.0,DFT,uIn,AUX2,potential);
	ABBAstep1D(tau/2.0,DFT,AUX2,uOut,potential);
	uOut*=bet;
	AUX*=alph;
	uOut+=AUX;
};
*/
//--------------------------------------
//--------------------------------------
/*
void ABBAstep(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax,
UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
	//AB step
	UCLAQ::GridFunction3d  AUX(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX.setToValue(0.0);
	heatEqStep(tau/2.0,xPanels,yPanels,zPanels,DFT,uIn,AUX,LX,LY,LZ);
	UCLAQ::GridFunction3d  AUX2(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX2.setToValue(0.0);
	PotentialStep(tau,xPanels,yPanels,zPanels,AUX, AUX2,potential);
	//BA step
	heatEqStep(tau/2.0,xPanels,yPanels,zPanels,DFT,AUX2,uOut,LX,LY,LZ);
	//Cstep
	//uOut/=uOut.norm2();
}
//--------------------------------------
//--------------------------------------
//--------------------------------------
//--------------------------------------

void Richardson(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax,
UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
	UCLAQ::GridFunction3d  AUX(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX.setToValue(0.0);
	UCLAQ::GridFunction3d  AUX2(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX.setToValue(0.0);
	uOut.setToValue(0.0);
	double alph=-1.0/7.0;
	double bet=8.0/7.0;
	ABBAstep(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,uIn,AUX,potential,LX,LY,LZ);
	//two half steps
	ABBAstep(tau/2.0,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,uIn,AUX2,potential,LX,LY,LZ);
	ABBAstep(tau/2.0,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,AUX2,uOut,potential,LX,LY,LZ);
	uOut*=bet;
	AUX*=alph;
	uOut+=AUX;
	uOut/=uOut.norm2();
};
//---------------------------
void RichardsonNOPROJ(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax,
UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
	UCLAQ::GridFunction3d  AUX(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX.setToValue(0.0);
	UCLAQ::GridFunction3d  AUX2(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX.setToValue(0.0);  
	uOut.setToValue(0.0);
	double alph=-1.0/7.0;
	double bet=8.0/7.0;
	ABBAstep(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,uIn,AUX,potential,LX,LY,LZ);
	//two half steps
	ABBAstep(tau/2.0,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,uIn,AUX2,potential,LX,LY,LZ);
	ABBAstep(tau/2.0,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,AUX2,uOut,potential,LX,LY,LZ);
	uOut*=bet;
	AUX*=alph;
	uOut+=AUX;
};
void Neri4thorder(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax, 
UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
	double x_0=-pow(2.0,1.0/3.0)/(2.0-pow(2.0,1.0/3.0));
	double x_1= 1.0             /(2.0-pow(2.0,1.0/3.0));
	UCLAQ::GridFunction3d  AUX(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX.setToValue(0.0);
	UCLAQ::GridFunction3d  AUX2(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
	AUX2.setToValue(0.0);
	uOut.setToValue(0.0);
	ABBAstep(tau*x_1,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,uIn,AUX,potential,LX,LY,LZ);
	ABBAstep(tau*x_0,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,AUX,AUX2,potential,LX,LY,LZ);
	ABBAstep(tau*x_1,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,AUX2,uOut,potential,LX,LY,LZ);
	uOut/=uOut.norm2();
}
*/
#include "TemplateRichardson.cpp"

int main(){
	ofstream fout;
	//timestep
	//double tau=0.05;
	double tau=0.2;
	//double tau=0.3;
	cout<<"Stepsize?";
	cin>>tau;
	string nameOfFile=to_string(tau);
	string threadCountInput = "-1";
	#ifdef _FFTW_OPENMP
    int threadIndex;
    int threadCount;
    if(not threadCountInput.empty()){
		threadCount = atoi(threadCountInput.c_str());
    }
    if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
    if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);
	printf("\n");
    printf("#############\n");
	printf("############# Using OpenMP With %d Threads\n",omp_get_max_threads());
	printf("#############\n");
	printf("\n");
    #endif
	long xPanels  = 100;
	//long yPanels  = 100;
    //long zPanels  = 100;
	if(xPanels==150){	
		nameOfFile="RichardsonResultsTemp1Dim150x150x150pt"+nameOfFile+".txt";//1DIM
		fout.open(nameOfFile.c_str());
		cout<<"FILE 150OPEN";
	}
	if(xPanels==100){
		nameOfFile="RichardsonResultsTemp1Dim100x100x100pt"+nameOfFile+".txt";//1DIM
		fout.open(nameOfFile.c_str());
		cout<<"FILE OPEN";
	}
	if(xPanels==50){
		nameOfFile="RichardsonResultsTemp1Dim50x50x50"+nameOfFile+".txt";//1DIM
		fout.open(nameOfFile.c_str());
		cout<<"FILE 50 OPEN";
	}
    // Reset panel count so efficient FFT's are used; the new
    // panel count is the next larger value that's a product
    // of primes < 13

    UCLAQ::FFT_Nvalues fft_Nvalues;

    printf("Original xPanels = %ld ",xPanels);
    xPanels = fft_Nvalues.getFFT_N(xPanels);
    printf("::: New xPanels = %ld \n",xPanels);
	/*    //3D stuff
    printf("Original yPanels = %ld ",yPanels);
    yPanels = fft_Nvalues.getFFT_N(yPanels);
    printf("::: New yPanels = %ld \n",yPanels);

    printf("Original zPanels = %ld ",zPanels);
    zPanels = fft_Nvalues.getFFT_N(zPanels);
    printf("::: New zPanels = %ld \n",zPanels);
	*/ 


	//   Problem set up.
	//
	//   Use non-unit domain so that this program demonstrates
	//   how to scale the transforms of derivatives appropriately
 
	double xMin   = -6.0;
	double xMax   =  6.0;
	double LX     = (xMax-xMin);
	/*
	double yMin   = -6.0;
	double yMax   =  6.0;

	double zMin   = -6.0;
	double zMax   =  6.0;


	double LY     = (yMax-yMin);
	double LZ     = (zMax-zMin);
	*/
	double pi     =  3.141592653589793238;

	double alpha  = 1.0;  // Coefficient of laplace operator
	double beta  = 1;
	UCLAQ::GridFunction1d nuclearPotential1D(xPanels,xMin,xMax);
	nuclearPotential1D.setToValue(0.0);
	//here I need the 1D version
	// The far field potential supplied by a HighOrderPolyPotential3d instance
	// is the potential
	//
	// V(r) = -(strength/(4*PI*|r - rPos|)
	//
	// where |r - rPos| = radial distance to rPos = (xPos,yPos,zPos)
	//
	// (For |r - rPos| < radius a smoothed value of this function is returned).
	//
	// In atomic units, the single particle Schroedinger equation for the hydrogen
	// atom has the potential
	//
	// -1/(|r - rPos|)
	//
	// so we set the strength to be 4*PI
	//
	double pi4      = 4.0*3.14159265358979323846;
	double strength = pi4;
	// (xPos,yPos,zPos) : location of potential
	double xPos   =  0.1;
	//double yPos   =  0.01;
	//double zPos   = -0.1;

	// radius : smoothing radius (epsilon)

	double radius =  10.0/(double)xPanels;

	// Instantiate and initialize the smoothed potential

	SmoothPolyMollifier1d smoothPotentialFun;//try mollifier 1d

	smoothPotentialFun.initialize(xPos,radius,strength);

    // Evaluate the smoothed nuclear potential on a grid

    nuclearPotential1D.specify(smoothPotentialFun.getEvaluationPtr());
	cout<<"minof nuclearPotential"<<nuclearPotential1D.min()<<endl;
	cout<<"maxof nuclearPotential"<<nuclearPotential1D.max()<<endl;
	//----------------------------------------------------------------------
   
	//compute exact eigenvalues NOT DONE RIGHT NOW FOR 1D,maybe later?
	//double exactEigValue;

	//exactEigValue=-2.29586;
    // Initialize the FFTW3 interface routine
    //UCLAQ::fftw3_sin3d DFT;                    // Discrete sin transform interface
	UCLAQ::fftw3_sin1d DFT;//1D discrete FT
    #ifdef _FFTW_OPENMP
	DFT.initialize(xPanels,threadCount);
	#else
	DFT.initialize(xPanels);
	#endif
	//----------------------------------------------------------------------------
    //Solving Schroedinger Equation
    //choose random initial vector:
	RandOpDirichlet1d RandOp;
	UCLAQ::GridFunction1d  Uinit(xPanels,xMin,xMax); 
	UCLAQ::GridFunction1d  UBeforeInnLoop(xPanels,xMin,xMax); 
	RandOp.randomize(Uinit);
	cout<<"minof U"<<Uinit.min()<<endl;
	cout<<"maxof nuclearPotential"<<Uinit.max()<<endl;
	//Create output vector
	UCLAQ::GridFunction1d  Uout(xPanels,xMin,xMax); 
	Uout.setToValue(0.0);
	// create Laplace operator WE WILL NEED IT TO EVALUATE THE  CANDIDATE EIGENVALUE VIA RAYLEIGH QUOTIENT
	LaplaceOp1d Lop1d(alpha);
	//create ShroedingerOperator using Laplace Operator and Potential function
	ShroedingOp1D ShrOP1d( Lop1d,nuclearPotential1D);
	// Temporary grids to store eigenvalue:
	UCLAQ::GridFunction1d  uTemp2(xPanels,xMin,xMax);
	UCLAQ::GridFunction1d  uTemp3(xPanels,xMin,xMax);
	uTemp2=Uinit;
	uTemp2/=uTemp2.norm2();
	uTemp3=uTemp2;
	//ShrOP3d.apply(uTemp2);
	ShrOP1d.apply(uTemp2);
	double lambdaAUX=uTemp2.scaledDot(uTemp3)/(uTemp3.scaledDot(uTemp3));
	cout<<"ORIGINAL VALUE"<<lambdaAUX<<endl;
	//-----------------------
	//-----------------------
	//-----------------------
	//COMPUTATION with richardson Extrapolation
	int maxIt=24000;
	int innMaxIt=24000;
	int i=1;
	double subspaceTol         = 2.0e-6;
	double innTol       = 0.001;
	double inTolLowBd=subspaceTol/10000.0;
	double outRelerr=subspaceTol+1;
	double inRelerr=outRelerr;
	double outLambdaprev=lambdaAUX;
	double lambdaSteadyState;
	bool abovebedrock_flag=1;
	double inLambdaprev;
	//while ((i<maxIt)&&(outRelerr>subspaceTol)&&(abovebedrock_flag))
	while ((outRelerr>subspaceTol)&&(abovebedrock_flag)){
		inRelerr=innTol+1;
		inLambdaprev=lambdaAUX;
		int jInner=0;
		UBeforeInnLoop=Uinit;
		//while(((jInner<innMaxIt)&&(inRelerr>innTol))||(jInner<=1))// or jinner>3
		while((inRelerr>tau*tau*innTol)||(jInner<=1)){
			//Richardson(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uinit,Uout,nuclearPotential,LX,LY,LZ);3D operator
			RichardsonRecANDprojT<UCLAQ::GridFunction1d,UCLAQ::fftw3_sin1d>(tau,DFT,Uinit,Uout,nuclearPotential1D,heatEqStep1D,PotentialStep1D,2);//done tests with this one
			//ABBAstepwithPROJ(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uinit,Uout,nuclearPotential,LX,LY,LZ);
 			//ABBAstepwithPROJ1D(tau,DFT,Uinit,Uout,nuclearPotential1D);//1D ABBAstepwith PROJ
			Uinit=Uout;
			//compute and print eigenvalue
			uTemp2=Uout;
			uTemp2/=uTemp2.norm2();
			uTemp3=Uout;
			//ShrOP3d.apply(uTemp2);
			ShrOP1d.apply(uTemp2);
			//RichardsonNOPROJ(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uout,uTemp2,nuclearPotential,LX,LY,LZ);//uTemp2 has Ax
			lambdaAUX=uTemp2.scaledDot(uTemp3)/(uTemp3.scaledDot(uTemp3));// eigenvalue estimate x'AX/||x||^2
			if((i%2)==1)
			{
				cout<<"Iteration number:   "<<i<<"*\n";
				cout<<"lambda etAetB "<<lambdaAUX<<endl;
				fout<<lambdaAUX<<endl;
				cout<<" interior error "<<inRelerr<<endl;
			}
      			inRelerr= abs( lambdaAUX - inLambdaprev)/(1+abs(lambdaAUX))  ;
			inLambdaprev=lambdaAUX;
			i++;
			jInner++;
		}
	cout << "\n";
	cout << "\n";
	cout<<"inner iterations:"<< jInner<<endl;
	uTemp2=Uout;
	uTemp2/=uTemp2.norm2();
	//ShrOP3d.apply(uTemp2);
	ShrOP1d.apply(uTemp2);
	lambdaSteadyState=uTemp2.scaledDot(uTemp3)/(uTemp3.scaledDot(uTemp3));// eigenvalue estimate x'AX/||x||^2
	
	cout<<"lambda steady state:"<<lambdaSteadyState<<endl;
	outRelerr=abs(lambdaSteadyState-outLambdaprev)/(1+abs(lambdaSteadyState));
	cout << "outerRelerr:" << outRelerr << endl;
	cout << "\n";
	cout << "\n";
	//if((outRelerr<inRelerr)||(jInner<=5)){
		
		/*if(jInner<=5)
		{//reset vector.
			lambdaSteadyState=outLambdaprev;
			Uinit=UBeforeInnLoop; 
			outRelerr=subspaceTol+1;
			cout<<"tightening inner tolerance bound"<<endl;
			innTol/=10.0;
			innTol=max(innTol,inTolLowBd);
			cout<<"new inner tolerance:"<<innTol;
		}
	//}
	else{
		cout<<"timestep update"<<endl;
		tau/=2.0;
	}
*/
	cout<<"timestep update"<<endl;
	tau/=2.0;
 	outLambdaprev= lambdaSteadyState;
	cout<<"enter a number to continue"<<endl;
	int trashpause;
	cin>>trashpause;
	}
	cout<<"bedrock reached"<<endl;
	//reached bedrock of accuracy that is possible by halving t/2, at this point just iterate
	int jInner=0;
	double bedTol=2.0e-10;
	double bedrockRelErr=bedTol+1;
	/*
	while((i<maxIt)&&(bedrockRelErr>bedTol))
	{
		if(i>0){
	RichardsonT<UCLAQ::GridFunction1d,UCLAQ::fftw3_sin1d>(tau,DFT,Uinit,Uout,nuclearPotential1D,heatEqStep1D,PotentialStep1D);
			}
		else
			ABBAstepwithPROJ(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uinit,Uout,nuclearPotential,LX,LY,LZ);
		Uinit=Uout;
		//compute and print eigenvalue
		uTemp2=Uout;
		uTemp2/=uTemp2.norm2();
		uTemp3=Uout;
		ShrOP3d.apply(uTemp2);
		//RichardsonNOPROJ(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uout,uTemp2,nuclearPotential,LX,LY,LZ);//uTemp2 has Ax
		lambdaAUX=uTemp2.scaledDot(uTemp3)/(uTemp3.scaledDot(uTemp3));// eigenvalue estimate x'AX/||x||^2
		if((i%2)==1)
		{
			cout<<"Iteration number"<<i<<endl;
			cout<<"lambda etAetB "<<lambdaAUX<<endl;
			fout<<lambdaAUX<<endl;
			cout<<" interior error "<<inRelerr<<endl;
		}
      	inRelerr= abs( lambdaAUX - inLambdaprev)/(1+abs(lambdaAUX))  ;
	inLambdaprev=lambdaAUX;

	i++;
	jInner++;
	}
	*/
	//UCLAQ::RandOpDirichlet3d  randomOp;
	UCLAQ::RandOpDirichlet1d randomOp;
	// Allocate arrays for eigenvectors and eigenvalues
    //vector <UCLAQ::GridFunction3d>  eigVectors;
	vector<UCLAQ::GridFunction1d> eigVectors;
    vector <double>                  eigValues;
	// Declare an instance of the Raylegh-Chebyshev eigensystem procedure
    RayleighChebyshev < UCLAQ::GridFunction1d,ShroedingOp1D, UCLAQ::RandOpDirichlet1d > RCeigProcedure;
    //                    |                |                       |
    //                    |                |                       |
    //             vector class    linear operator class     randomize operator class

    RCeigProcedure.setEigDiagnosticsFlag();
    RCeigProcedure.setVerboseFlag();
	//UCLAQ::GridFunction3d vTmp(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);     // A temporary vector is required as input. This vector must
                                                                         // be a non-null instance of the vector class

	UCLAQ::GridFunction1d vTmp(xPanels,xMin,xMax);
    double dimension           = vTmp.getDimension();
    long subspaceIncrementSize = 1;
    long bufferSize            = 1;
    long eigCount              = dimension < 1 ? dimension : 1;
	//get the eigenvalues for Shroedinger Operator ShrOP3d
    RCeigProcedure.getMinEigenSystem(eigCount, subspaceTol, subspaceIncrementSize, bufferSize, vTmp,
	ShrOP1d, randomOp, eigValues, eigVectors);//CHange this to SHROP3d in 3D case
	//PRINT OUT RESULTS 
    printf("\n\nXXXX   RC_OperatorEig_Test Results XXXX\n\n");
    printf("Tolerance Specified : %10.5e\n\n",subspaceTol);
    printf("       Eigenvalue\n");
    for(long  k = 0; k < eigCount; k++ )
    {
    	printf("%-5ld %-10.5e \n", k+1, eigValues[k]);
    }
	printf("\nXXX Execution Completed XXXX\n");
	fout.close();
}
