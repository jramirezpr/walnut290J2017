#include <iostream>
#include <cmath>
#include <functional>
#include<fstream>
#include<string>
using namespace std;

#ifdef _FFTW_OPENMP
#include <omp.h>
#endif
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
// g++ RichardsonExtr.cpp -std=c++11 -I../../290J_2015 -lfftw3 -o Richfixedprojpt2.exe
//

//

#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "RandOpNd/RandOpDirichlet3d.h"
#include "ShroedingOp3D1Par.cpp"

#include "FFTW3_InterfaceNd/UCLAQ_fftw3_sin3d.h"
#include "FFTW3_InterfaceNd/UCLAQ_FFT_Nvalues.h"
#include "MollifierNd/HighOrderPolyPotential3d.h"           // Mollified nuclear potential
#include "RayleighChebyshev/RayleighChebyshev.h"

void heatEqStep(double tau, int xPanels, int yPanels,int zPanels, UCLAQ::fftw3_sin3d& DFT,       UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,double LX,double LY,double LZ,double alpha  = -1.0){
double threshold=exp(20);
int count=0;
int xm=-1,ym=-1,zm=-1;	

double pi     =  3.141592653589793238;
double minuseigkxyz=0;
UCLAQ::DoubleVector3d fTransform(xPanels-1,yPanels-1,zPanels-1);
    fTransform.setToValue(0.0);
long sizegrid=(xPanels-1)*(yPanels-1)*(zPanels-1);
    DFT.fftw3d_sin_forward(uIn,fTransform); // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d

    for(long kx = 1; kx <= xPanels-1; kx++)
    {
    for(long ky = 1; ky <= yPanels-1; ky++)
    {
    for(long kz = 1; kz <= zPanels-1; kz++)
    {
    	
 minuseigkxyz=tau*alpha*( (((kx*pi)/LX)*((kx*pi)/LX)) +  (((ky*pi)/LY)*((ky*pi)/LY)) +  (((kz*pi)/LZ)*((kz*pi)/LZ)) );
 /*if ( exp(minuseigkxyz)>threshold){ //UNCOMMENT THESE LINES WHEN THRESHOLDING
count++;
        //max= minuseigkxyz;
	xm=kx; ym=ky;zm=kz;
}
else{
fTransform(kx-1,ky-1,kz-1) *=exp(tau*alpha*( (((kx*pi)/LX)*((kx*pi)/LX)) +  (((ky*pi)/LY)*((ky*pi)/LY)) +  (((kz*pi)/LZ)*((kz*pi)/LZ)) ));// its exp(-tau*alpha*eig_kx,ky)
}
*/


//COMMENT THE FOLLOWING LINE WHEN THRESHOLDING
fTransform(kx-1,ky-1,kz-1) *=exp(tau*alpha*( (((kx*pi)/LX)*((kx*pi)/LX)) +  (((ky*pi)/LY)*((ky*pi)/LY)) +  (((kz*pi)/LZ)*((kz*pi)/LZ)) ));// its exp(-tau*alpha*eig_kx,ky)
    }}}
//cout<<"HEAT MAX:"<<max<<" exp  :"<< exp(max)<<" Type 0 to continue \n";
 //UNCOMMENT WHEN THRESHOLDING
if(count!=0){

//cout<<"NUMBER OF ENTRIES LARGER THAN THRESHOLD:"<<count<<" Type 0 to continue \n";
//cout<<"Percentage of entries larger than "<<(double)count/(double)(sizegrid)<<" Type 0 to continue \n";
//cin>>count;
}

    DFT.fftw3d_sin_inverse(fTransform,uOut);      // Note: We can use a GridFunction3d as an argument since it's extended
                                               // from a DoubleVector3d



}
void PotentialStep(double tau, int xPanels, int yPanels, int zPanels,        UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential){
double max=0;
int xm=-1,ym=-1,zm=-1;
for(long kx = 1; kx <= potential.xPanels-1; kx++)
    {
    for(long ky = 1; ky <= potential.yPanels-1; ky++)
    {
    for(long kz = 1; kz <= potential.zPanels-1; kz++)
    {
    	uOut.Values(kx,ky,kz) =(uIn.Values(kx,ky,kz))* exp(-tau*potential.Values(kx,ky,kz));
if (-tau*potential.Values(kx,ky,kz)>max){
max=-tau*potential.Values(kx,ky,kz);
	xm=kx; ym=ky;zm=kz;

}

    }}}


//cout<<"POTENTIAL MAX:"<<max<<" exp  :"<< exp(max)<<" Type 0 to continue \n";

//cin>>max;



}

void ABBAstepwithPROJ(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax, UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
 
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
uOut/=uOut.norm2();


}







//--------------------------------------
//--------------------------------------
//--------------------------------------
//--------------------------------------
void ABBAstep(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax, UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
 
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


void Richardson(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax, UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
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


void Neri4thorder(double tau, int xPanels,double xMin,double xMax, int yPanels,double yMin,double yMax, int zPanels,double zMin,double zMax, UCLAQ::fftw3_sin3d& DFT,UCLAQ::GridFunction3d&  uIn,    UCLAQ::GridFunction3d&  uOut,UCLAQ::GridFunction3d& potential,double LX,double LY,double LZ){
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



int main()
{
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
    if(not threadCountInput.empty())
    {
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
	long yPanels  = 100;
    long zPanels  = 100;
if(xPanels==150){
nameOfFile="RichardsonResults150x150x150pt"+nameOfFile+".txt";
fout.open(nameOfFile.c_str());
cout<<"FILE 150OPEN";

}
if(xPanels==100){
nameOfFile="RichardsonResults100x100x100pt"+nameOfFile+".txt";
fout.open(nameOfFile.c_str());
cout<<"FILE OPEN";

}
if(xPanels==50){
nameOfFile="RichardsonResults50x50x50"+nameOfFile+".txt";
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

    printf("Original yPanels = %ld ",yPanels);
    yPanels = fft_Nvalues.getFFT_N(yPanels);
    printf("::: New yPanels = %ld \n",yPanels);

    printf("Original zPanels = %ld ",zPanels);
    zPanels = fft_Nvalues.getFFT_N(zPanels);
    printf("::: New zPanels = %ld \n",zPanels);

//   Problem set up.
//
//   Use non-unit domain so that this program demonstrates
//   how to scale the transforms of derivatives appropriately
 
	double xMin   = -6.0;
	double xMax   =  6.0;

	double yMin   = -6.0;
	double yMax   =  6.0;

	double zMin   = -6.0;
	double zMax   =  6.0;

	double LX     = (xMax-xMin);
	double LY     = (yMax-yMin);
	double LZ     = (zMax-zMin);

	double pi     =  3.141592653589793238;

	double alpha  = 1.0;  // Coefficient of laplace operator
	double beta  = 1;
UCLAQ::GridFunction3d nuclearPotential(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
	nuclearPotential.setToValue(0.0);

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
	double yPos   =  0.01;
	double zPos   = -0.1;

	// radius : smoothing radius (epsilon)

	double radius =  10.0/(double)xPanels;

	// Instantiate and initialize the smoothed potential

	HighOrderPolyPotential3d smoothPotentialFun;

	smoothPotentialFun.initialize(xPos,yPos,zPos,radius,strength);

    // Evaluate the smoothed nuclear potential on a grid

    nuclearPotential.specify(smoothPotentialFun.getEvaluationPtr());
cout<<"minof nuclearPotential"<<nuclearPotential.min()<<endl;
cout<<"maxof nuclearPotential"<<nuclearPotential.max()<<endl;
//----------------------------------------------------------------------
    



//compute exact eigenvalues
double exactEigValue;


	exactEigValue=-2.29586;


    

    // Initialize the FFTW3 interface routine

    UCLAQ::fftw3_sin3d DFT;                    // Discrete sin transform interface


    #ifdef _FFTW_OPENMP
	DFT.initialize(xPanels,yPanels,zPanels,threadCount);
	#else
	DFT.initialize(xPanels,yPanels,zPanels);
	#endif

//----------------------------------------------------------------------------
    //Solving Schroedinger Equation

    ///choose random initial vector:
RandOpDirichlet3d RandOp;
UCLAQ::GridFunction3d  Uinit(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
RandOp.randomize(Uinit);
cout<<"minof U"<<Uinit.min()<<endl;
cout<<"maxof nuclearPotential"<<Uinit.max()<<endl;


//Create output vector
UCLAQ::GridFunction3d  Uout(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax); 
Uout.setToValue(0.0);

// create Laplace operator WE WILL NEED IT TO EVALUATE THE  CANDIDATE EIGENVALUE VIA RAYLEIGH QUOTIENT
   LaplaceOp3D Lop3d(alpha);
//create ShroedingerOperator using Laplace Operator and Potential function

ShroedingOp3D ShrOP3d( Lop3d,nuclearPotential);
// Temporary grids to store eigenvalue:
 UCLAQ::GridFunction3d  uTemp2(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
 UCLAQ::GridFunction3d  uTemp3(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);
uTemp2=Uinit;
uTemp2/=uTemp2.norm2();
uTemp3=uTemp2;
ShrOP3d.apply(uTemp2);
 double lambdaAUX=uTemp2.scaledDot(uTemp3)/(uTemp3.scaledDot(uTemp3));
cout<<"ORIGINAL VALUE"<<lambdaAUX<<endl;
//-----------------------
//-----------------------
//-----------------------
//COMPUTATION with richardson Extrapolation
for(int i=1;i<351;i++){
if(i>0)
{      //Neri4thorder(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uinit,Uout,nuclearPotential,LX,LY,LZ);
	Richardson(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uinit,Uout,nuclearPotential,LX,LY,LZ);
}
else
ABBAstepwithPROJ(tau,xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax,DFT,Uinit,Uout,nuclearPotential,LX,LY,LZ);
 
Uinit=Uout;
/*
if(((i%50)==0)      )
{tau=tau/2.0;
cout<<tau<<" HERE";
}
*/
if((i==60))
{tau=tau/2.0;
cout<<tau<<" HERE";
}

if((i==100))
{tau=tau/2.0;
cout<<tau<<" HERE";
}

if((i%2)==1)
	{//compute and print eigenvalue
uTemp2=Uout;
uTemp2/=uTemp2.norm2();
uTemp3=uTemp2;
ShrOP3d.apply(uTemp2);
double lambdaAUX=uTemp2.scaledDot(uTemp3)/(uTemp3.scaledDot(uTemp3));
cout<<"lambda"<<lambdaAUX<<endl;
fout<<lambdaAUX<<endl;


	}





		}

  UCLAQ::RandOpDirichlet3d  randomOp;
// Allocate arrays for eigenvectors and eigenvalues

    vector <UCLAQ::GridFunction3d>  eigVectors;
    vector <double>                  eigValues;
// Declare an instance of the Raylegh-Chebyshev eigensystem procedure

    RayleighChebyshev < UCLAQ::GridFunction3d,ShroedingOp3D, UCLAQ::RandOpDirichlet3d > RCeigProcedure;
    //                    |                |                       |
    //                    |                |                       |
    //             vector class    linear operator class     randomize operator class

    RCeigProcedure.setEigDiagnosticsFlag();
    RCeigProcedure.setVerboseFlag();

UCLAQ::GridFunction3d vTmp(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);     // A temporary vector is required as input. This vector must
                                                                         // be a non-null instance of the vector class

    double dimension           = vTmp.getDimension();
    double subspaceTol         = 2.0e-6;
    long subspaceIncrementSize = 1;
    long bufferSize            = 1;
    long eigCount              = dimension < 1 ? dimension : 1;
//get the eigenvalues for Shroedinger Operator ShrOP3d
    RCeigProcedure.getMinEigenSystem(eigCount, subspaceTol, subspaceIncrementSize, bufferSize, vTmp,
    		ShrOP3d, randomOp, eigValues, eigVectors);
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
