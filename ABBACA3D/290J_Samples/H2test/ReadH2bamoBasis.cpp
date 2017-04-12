#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
using namespace std;
//
// ReadH2bamoBasis.cpp
//
// A test code that demonstrates reading a file containing a
// vector of orbital basis functions and reading a file
// containing the potential associated with a single particle
// Schroedinger operator.
//
// The orbitals basis are eigenfuncitons of a single particle
// Schroedinger operator with the potential induced by two
// nucleii located at A and B, points in R^3
// (Bare Atomic Molecular Orbitals).
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

#include "DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3d.h"
#include "GridFunctionNd/UCLAQ_GridFunction3dUtility.h"

#include "Qdomain/QdomainBasisFunctions.h"


int main()
{
	// Names of input files to be read

	string orbitalBasisFileName = "H2bamoBasis.dat";
    string potentialFilename    = "H2nuclearPotential.dat";

    //####################################
    //    Reading in basis functions
    //####################################

	long orbitalBasisSize       = -1;  // -1 indicates all basis vectors will be read
	bool verboseFlag            = true;

	FILE* dataFile = 0;
	if((dataFile = fopen(orbitalBasisFileName.c_str(), "rb" )) == 0 )
	{
	string msg ="\nXXX initializeOrbitalBasis(...) XXX \n";
	msg.append("External basisl file \"" + orbitalBasisFileName + "\" not found \n");
    throw std::runtime_error(msg);
	}

	UCLAQ::QdomainBasisFunctions basisFunctionUtility;

	if(verboseFlag)
	{
	basisFunctionUtility.setVerboseFlag();
	}

	vector<double>                orbitalBasisIndexValues;
	vector<UCLAQ::GridFunction3d> orbitalBasis;

	basisFunctionUtility.inputFromFile(orbitalBasisIndexValues,orbitalBasis,orbitalBasisFileName,orbitalBasisSize);

	//
	// Check norms of input basis functions
	//

   long basisSize = orbitalBasis.size();

   cout << " Basis Size " << basisSize << endl;
   for(long i = 0; i < basisSize; i++)
   {
   printf("Norm of %ld basis function : %10.5e \n",i, orbitalBasis[i].norm2());
   }

   //####################################
   //    Reading in potential
   //####################################

   UCLAQ::GridFunction3d H2nuclearPotential;


   int noFileFlag  = 0;
   UCLAQ::GridFunction3dUtility gUtility3d;
   gUtility3d.inputFromBinaryDataFile(H2nuclearPotential, potentialFilename, noFileFlag);

   if(noFileFlag)
   {
	string msg ="\nXXX inputFromBinaryDataFile(...) XXX \n";
	msg.append("Specified input file \"" +potentialFilename + "\" not found \n");
    throw std::runtime_error(msg);
   }

   // Print out data
   printf("\n");
   printf(" Potential Grid Structure \n\n");
   printf(" [xPanels, xMin, xMax] : %ld, %10.5f %10.5f \n",H2nuclearPotential.xPanels,H2nuclearPotential.xMin, H2nuclearPotential.xMax);
   printf(" [yPanels, yMin, yMax] : %ld, %10.5f %10.5f \n",H2nuclearPotential.yPanels,H2nuclearPotential.yMin, H2nuclearPotential.yMax);
   printf(" [zPanels, zMin, zMax] : %ld, %10.5f %10.5f \n",H2nuclearPotential.zPanels,H2nuclearPotential.zMin, H2nuclearPotential.zMax);
   printf("\n");
   printf(" Potential Norm_Inf: %20.15e   Norm_L2: %20.15e  \n\n",H2nuclearPotential.normInf(),H2nuclearPotential.norm2());

   printf("XXXX Execution Complete XXXXX\n");

}

