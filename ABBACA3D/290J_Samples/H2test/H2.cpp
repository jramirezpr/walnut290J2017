#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
using namespace std;
//
// H2.cpp
//
// A test code that demonstrates the use of a finite difference
// procedure to  create approximate eigenfunctions of the two particle
// Schroedinger operator with the potential induced by two nucleii
// located at A and B, points in R^3.
//
// The basis set consists of the first M eigenfunctions (M specified) 
// of the operator  
//
//  [ -(1/2)*DELTA  + V_epsilon(x) ]
//
// where the Laplace operator DELTA is a high order finite difference approximation
// with homogeneous boundary conditions and V(x) is a smoothed nuclear potential 
// of unit strength with centers at (xA,yA,zA) and (xB,yB,zB). The smoothing
// radius is epsilon. 
//
//
// In the FCI procedure used, the single particle operator is also taken to be 
// 
// [ -(1/2)*DELTA  + V_epsilon(x) ]
//
// but the electron-electron interaction potential is not mollified. 
//

#include "Timing/ClockIt.h"


// Problem specific classes

// Single particle operator combining a high order 
// approximatio of the Laplace operator with 
// a smoothed nuclear potential

#include "SmoothedH2op.h" 

// Classes for basis set construction as eigenfunctions of the
// single particle operator (Bare Atomic Molecular Orbitals)

#include "RandOpNd/RandOpDirichlet3d.h"
#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "RayleighChebyshev/RayleighChebyshev.h"

// Classes needed for orbital integral construction 

#include "ExpandingDomainNd/InversePoisson3d.h"                     // Infinite domain Poisson solver
#include "MultiParticleHamiltonian/RealSpatialOrbitalIntegrals.h"   // Orbital integral value data structure
#include "Qdomain/QdomainOrbitalIntegrals.h"

// Qdomain classes for using orbital integral data to construct a 
// Ritz representation of the multi-particle Hamiltonian in a Slater 
// determinant basis, find the lowest states, and make available the 
// energies and electron densities associated with those states.   
// 

#include "Qdomain/QdomainFCI.h"
#include "Qdomain/QdomainFCIdensity.h"


int main()
{
	 ClockIt clockIt;

	long panelCount               = 50;
	double domainWidth            = 5.0;
    double epsilon                = 0.5;

	int threadCount = 8;

	#ifdef _OPENMP
    if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
    if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);

	printf("\n");
    printf("#############\n");
	printf("############# Using OpenMP With %d Threads\n",omp_get_max_threads());
	printf("#############\n");
	printf("\n");
    #endif

	//   Problem set up.

	long xPanels = panelCount; long yPanels = panelCount; long zPanels = panelCount;

	double xMin = -domainWidth - 0.5; double xMax = domainWidth + 0.5;
	double yMin = -domainWidth; double yMax = domainWidth;
	double zMin = -domainWidth; double zMax = domainWidth;

	// Location of nuclear charges (using ground state separation)

	vector<double> A = {-0.70055,0.0,0.0};
	vector<double> B = { 0.70055,0.0,0.0};

//##################################################################
//   Initializing smoothed single particle H2 Op
//##################################################################

	SmoothedH2op sOp(A,B, epsilon);

//##################################################################
//   Creating basis elements
//##################################################################

	// Instantiate Random Operator

    UCLAQ::RandOpDirichlet3d  randomOp;

    // Allocate arrays for eigenvectors and eigenvalues

    vector <UCLAQ::GridFunction3d>  eigVectors;
    vector <double>                  eigValues;

    // Declare an instance of the Raylegh-Chebyshev eigensystem procedure

    RayleighChebyshev < UCLAQ::GridFunction3d, SmoothedH2op , UCLAQ::RandOpDirichlet3d > RCeigProcedure;
    //                    |                      |                       |
    //                    |                      |                       |
    //             vector class        linear operator class     randomize operator class

    RCeigProcedure.setEigDiagnosticsFlag();
    RCeigProcedure.setVerboseFlag();

    UCLAQ::GridFunction3d vTmp(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);     // A temporary vector is required as input. This vector must
                                                                                           // be a non-null instance of the vector class

    long basisSize             = 20;
    long dimension             = vTmp.getDimension();
    double subspaceTol         = 2.0e-5;
    long subspaceIncrementSize = basisSize;
    long bufferSize            = 5;
    long eigCount              = dimension < 10 ? dimension : 10;


    printf("Creating Eigensystem : ");
    fflush(stdout);

    clockIt.start();

    RCeigProcedure.getMinEigenSystem(eigCount, subspaceTol, subspaceIncrementSize, bufferSize, vTmp, sOp, randomOp, eigValues, eigVectors);

    clockIt.stop();
	printf("\n\n");
	printf("Done --- Time (sec) : %10.5f \n\n",clockIt.getSecElapsedTime());

    printf("\n\nXXXX   RC_OperatorEig_Test Results XXXX\n\n");
    printf("Tolerance Specified : %10.5e\n\n",subspaceTol);

    printf("       Eigenvalue    \n");
    for(long  k = 0; k < eigCount; k++ )
    {
    	printf("%-5ld %-10.5e   \n", k+1, eigValues[k]);
    }

    double eigNorm = 0.0;

    for(long  k = 0; k < eigCount; k++ )
    {
    	eigNorm = sqrt(abs(eigVectors[k].scaledDot(eigVectors[k])));
        eigVectors[k] /= eigNorm;
    }

//##################################################################
//   Creating orbital integral data
//##################################################################

    // Orbital integral data

    bool verboseFlag                 = true;
    bool timingConstantsFlag         = false;
    double invLapExtFactor           = 2.0;

    UCLAQ::RealSpatialOrbitalIntegrals realSpatialOrbitalIntegrals;

    UCLAQ::QdomainOrbitalIntegrals < InversePoisson3d , SmoothedH2op >    qDomainOrbitalIntegrals;

    qDomainOrbitalIntegrals.initialize();
    qDomainOrbitalIntegrals.setInvLaplaceExtFactor(invLapExtFactor);


    // laplaceFixUpCoefficient =  the coefficient of the Laplace operator so that when the 
    // INVERSE Laplace operator routine is called it returns -1/r convolved with arguments.

    double laplaceFixUpCoefficient      = 1.0/(4.0*3.141592653589793238);

    qDomainOrbitalIntegrals.setLaplaceOpCoefficient(laplaceFixUpCoefficient);

    qDomainOrbitalIntegrals.setVerbose(verboseFlag);
    qDomainOrbitalIntegrals.setTimingConstantsFlag(timingConstantsFlag);

    try{qDomainOrbitalIntegrals.attachBasisVectors(eigVectors);}
    catch (const runtime_error& e){throw e;}


    qDomainOrbitalIntegrals.attachSingleParticleOperator(&sOp);

    realSpatialOrbitalIntegrals.initialize();
    qDomainOrbitalIntegrals.createOrbitalIntegrals(realSpatialOrbitalIntegrals);

//
// Compute potential energy contriution (energyshift) to the distribution of nucleii; 
// this energy is added to the energies of electronic energy to determine the total energy. 
//
    double rdist       = sqrt((A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1])+ (A[2]-B[2])*(A[2]-B[2]));
    double energyShift = 1.0/rdist;


//##################################################################
//   Construct FCI matrix from orbital integral data
//   and then create eigensystem
//##################################################################

    // FCI output data

    vector < string >                       spinBlockIdentifiers;
    vector <vector <double> >                   spinBlockEnergies;
    vector < vector <UCLAQ::DoubleVector2d> >  spinBlockDensities;

     // FCI computation data members

    UCLAQ::QdomainFCI     qDomainFCI;

    bool verboseEigFlag    = verboseFlag;
    long eigCountMax       = 5;
    subspaceIncrementSize  = 5;
    bufferSize             = 3;

    qDomainFCI.initialize();
    qDomainFCI.clearVerbose();
    qDomainFCI.clearVerboseEigensystemComputation();
    qDomainFCI.setVerbose(verboseFlag);
    qDomainFCI.setVerboseEigensystemComputation(verboseEigFlag);

    qDomainFCI.setEigCountMax(eigCountMax);
    qDomainFCI.setSubspaceTol(subspaceTol);
    qDomainFCI.setSubspaceIncrementSize(subspaceIncrementSize);
    qDomainFCI.setBufferSize(bufferSize);

    spinBlockIdentifiers.clear();
    spinBlockEnergies.clear();
    spinBlockDensities.clear();


    int electronCount = 2;

    // Using all eigenfunctions that have been computed 

    long spatialOrbitalBasisCount = eigVectors.size();

    qDomainFCI.run(electronCount, spatialOrbitalBasisCount, energyShift, realSpatialOrbitalIntegrals,
    spinBlockIdentifiers,spinBlockEnergies,spinBlockDensities);


//##################################################################
//   Output energies
//##################################################################


    for(long k = 0; k < (long)spinBlockIdentifiers.size(); k++)
    {
    cout << spinBlockIdentifiers[k] << endl;
    for(long j = 0; j < (long)spinBlockEnergies[k].size(); j++)
    {
    cout << spinBlockIdentifiers[k] << " " << j << " : " << spinBlockEnergies[k][j] << endl;
    }
    }

//##################################################################
//   Create density
//##################################################################

   UCLAQ::QdomainFCIdensity qDomainFCIdensity;

   qDomainFCIdensity.initialize(electronCount, spinBlockIdentifiers,
   spinBlockEnergies, spinBlockDensities, eigVectors);

   UCLAQ::GridFunction3d H2_0;
   UCLAQ::GridFunction3d H2_1;
   UCLAQ::GridFunction3d H2_2;

   qDomainFCIdensity.getFCIdensity(0,"UD",H2_0);
   qDomainFCIdensity.getFCIdensity(1,"UD",H2_1);
   qDomainFCIdensity.getFCIdensity(2,"UD",H2_2);


//##################################################################
//   Output results
//##################################################################

   UCLAQ::GridFunction3dUtility gUtility;
   gUtility.outputDataToVTKfile(H2_0,"H2_0.vtk","H2_0");
   gUtility.outputDataToVTKfile(H2_1,"H2_1.vtk","H2_1");
   gUtility.outputDataToVTKfile(H2_2,"H2_2.vtk","H2_2");


	printf("XXXX Execution Complete XXXXX\n");

}

