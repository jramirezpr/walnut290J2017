/*
#############################################################################
#
# Copyright 2015-16 Chris Anderson
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


//
//  RealSpatialOrbitalIntegrals.h 
//
#ifndef __RealSpatialOrbitalIntegrals__
#define __RealSpatialOrbitalIntegrals__

#include <sstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <stdexcept>
using namespace std;

namespace UCLAQ
{

class RealSpatialOrbitalIntegrals
{

public :

RealSpatialOrbitalIntegrals()
{
    orbitalEnergies           = 0;
    singleOrbitalIntegralData = 0;
    doubleOrbitalIntegralData = 0;

    orbitalBasisCount           = 0;
    singleOrbitalIntegralCount  = 0;
    doubleOrbitalIntegralCount  = 0;
};


RealSpatialOrbitalIntegrals(const RealSpatialOrbitalIntegrals& R)
{
    orbitalEnergies           = 0;
    singleOrbitalIntegralData = 0;
    doubleOrbitalIntegralData = 0;

    orbitalBasisCount           = R.orbitalBasisCount ;
    singleOrbitalIntegralCount  = R.singleOrbitalIntegralCount;
    doubleOrbitalIntegralCount  = R.doubleOrbitalIntegralCount;

    orbitalEnergies            = new double[orbitalBasisCount];
    singleOrbitalIntegralData  = new double[singleOrbitalIntegralCount];
    doubleOrbitalIntegralData  = new double[doubleOrbitalIntegralCount];

    long i;

    for(i = 0; i < orbitalBasisCount; i++)
    {orbitalEnergies[i] = R.orbitalEnergies[i];}

    for(i = 0; i < singleOrbitalIntegralCount; i++)
    {singleOrbitalIntegralData[i] = R.singleOrbitalIntegralData[i];}

    for(i = 0; i < doubleOrbitalIntegralCount; i++)
    {doubleOrbitalIntegralData[i] = R.doubleOrbitalIntegralData[i];}
};

// Creates instance with internal data structures to hold single and double orbital
// integral data associated with a set of orbitalBasisCount spatial orbitals.


RealSpatialOrbitalIntegrals(long orbitalBasisCount)
{
	orbitalEnergies           = 0;
    singleOrbitalIntegralData = 0;
    doubleOrbitalIntegralData = 0;
	initialize(orbitalBasisCount);
}

void initialize()
{
	destroyData();

	orbitalEnergies           = 0;
    singleOrbitalIntegralData = 0;
    doubleOrbitalIntegralData = 0;

    orbitalBasisCount           = 0;
    singleOrbitalIntegralCount  = 0;
    doubleOrbitalIntegralCount  = 0;
}


// Initializes an instance with internal data structures to hold single and double orbital
// integral data associated with a set of orbitalBasisCount spatial orbitals.

void initialize(long orbitalBasisCount)
{
	destroyData();

	this->orbitalBasisCount = orbitalBasisCount;

	long p                 = this->orbitalBasisCount;

	this->singleOrbitalIntegralCount = (p)*(p+1)/2;
	this->doubleOrbitalIntegralCount = (p*(p*p*p + 2*p*p + 3*p +2))/8;

	this->orbitalEnergies     = new double[this->orbitalBasisCount];
    singleOrbitalIntegralData = new double[this->singleOrbitalIntegralCount];
    doubleOrbitalIntegralData = new double[this->doubleOrbitalIntegralCount];

    // Assign artificial energies to each orbital

    for(long i = 0; i < this->orbitalBasisCount; i++)
    {
    orbitalEnergies[i] = (double)(i);
    }

    // Set single and double orbital integral values to 0.0

    for(long i = 0; i < this->singleOrbitalIntegralCount; i++)
    {
    	singleOrbitalIntegralData[i] = 0.0;
    }

    for(long i = 0; i < this->doubleOrbitalIntegralCount; i++)
    {
    	doubleOrbitalIntegralData[i] = 0.0;
    }
}

~RealSpatialOrbitalIntegrals() {destroyData();};


/**
// This function sets the single orbital integral data associated with
// the integral
//
// phi_i(r1)*H_singleParticle(r1)*phi_j(r1) dr1
//
// i and j can take on values from 0 to orbitalBasisCount -1
//
// This function assumes that the singleParticleOperator is a symmetric operator.
*/

void setSingleOrbitalIntegralData(long i, long j, double value)
{
#ifdef _DEBUG
   boundsCheck(i,orbitalBasisCount,1);
   boundsCheck(j,orbitalBasisCount,2);
#endif

   // Shift indices by +1 since external indexing starts at 0

   long p = i+1; long q = j+1;

   if(p >= q)
   {
	   singleOrbitalIntegralData[getSingleOrbitalDataIndex(p,q)] = value;
   }
   else
   {
	   singleOrbitalIntegralData[getSingleOrbitalDataIndex(q,p)] = value;
   }
}
/**
// This function returns the single orbital integral data associated with
// the integral
//
// phi_i(r1)*H_singleParticle(r1)*phi_j(r1) dr1
//
// i and j can take on values from 0 to orbitalBasisCount -1
//
// This function assumes that the singleParticleOperator is a symmetric operator.
*/

double getSingleOrbitalIntegralData(long i, long j) const
{
#ifdef _DEBUG
   boundsCheck(i,orbitalBasisCount,1);
   boundsCheck(j,orbitalBasisCount,2);
#endif

   // Shift indices by +1 since external indexing starts at 0

   long p = i+1; long q = j+1;

   if(p >= q)
   {
	   return singleOrbitalIntegralData[getSingleOrbitalDataIndex(p,q)];
   }
   else
   {
	   return singleOrbitalIntegralData[getSingleOrbitalDataIndex(q,p)];
   }
}
/**
// This function sets the value of the internal data structure for the integral of
//
// phi_m(r2)*phi_n(r2)*G(r1,r2)*phi_i(r1)*phi_j(r1) dr1 dr2
//
// Note!!! : In this indexing each pair of indices (i,j) and (m,n)
//           are associated with the same integration variable.
//
// All indices range from 0 to orbitalBasisCount-1
*/

void setDoubleOrbitalIntegralData(long i, long j, long m, long n, double value)
{

#ifdef _DEBUG
   boundsCheck(i,orbitalBasisCount,1);
   boundsCheck(j,orbitalBasisCount,2);
   boundsCheck(m,orbitalBasisCount,3);
   boundsCheck(n,orbitalBasisCount,4);
#endif

// Offset indices so they start with 1

    long p = i+1; long r = m+1;
    long q = j+1; long s = n+1;

    long iTmp;

//  In this indexing
//
//  (p,q) is associated with r1
//  (r,s) is associated with r2
//
//  The integral is obtained from the data file
//  with (p,q) and (r,s) ordered such that
//
//  p >= q  and r >= s, p >= r, and when p = r then q >= s
//

    // Re-order pairs so second index of each pair
    // is less than the first

    if(p < q){iTmp = p; p = q; q = iTmp;} //  => (p >= q)
    if(r < s){iTmp = r; r = s; s = iTmp;} //  => (r >= s)

    // Re-order pairs so p >= r

    if(p < r)
    {
     iTmp = p; p = r; r = iTmp;
     iTmp = q; q = s; s = iTmp;
    }

    // when p = r insure that q >= s

    if((p == r)&&(q < s)) {iTmp = q; q = s; s = iTmp;}

    doubleOrbitalIntegralData[getDoubleOrbitalDataIndex(p, q, r, s)] = value;
}

/**
// This function returns the value of the internal data structure for the integral of
//
// phi_m(r2)*phi_n(r2)*G(r1,r2)*phi_i(r1)*phi_j(r1) dr1 dr2
//
// Note!!! : In this indexing each pair of indices (i,j) and (m,n)
//           are associated with the same integration variable.
//
// All indices range from 0 to orbitalBasisCount-1
*/

double getDoubleOrbitalIntegralData(long i, long j, long m, long n) const
{

#ifdef _DEBUG
   boundsCheck(i,orbitalBasisCount,1);
   boundsCheck(j,orbitalBasisCount,2);
   boundsCheck(m,orbitalBasisCount,3);
   boundsCheck(n,orbitalBasisCount,4);
#endif

// Offset indices so they start with 1

    long p = i+1; long r = m+1;
    long q = j+1; long s = n+1;

    long iTmp;

//  In this indexing
//
//  (p,q) is associated with r1
//  (r,s) is associated with r2
//
//  The integral is obtained from the data file
//  with (p,q) and (r,s) ordered such that
//
//  p >= q  and r >= s, p >= r, and when p = r then q >= s
//

    // Re-order pairs so second index of each pair
    // is less than the first

    if(p < q){iTmp = p; p = q; q = iTmp;} //  => (p >= q)
    if(r < s){iTmp = r; r = s; s = iTmp;} //  => (r >= s)

    // Re-order pairs so p >= r

    if(p < r)
    {
     iTmp = p; p = r; r = iTmp;
     iTmp = q; q = s; s = iTmp;
    }

    // when p = r insure that q >= s

    if((p == r)&&(q < s)) {iTmp = q; q = s; s = iTmp;}

    return doubleOrbitalIntegralData[getDoubleOrbitalDataIndex(p, q, r, s)];
}


/// Returns the number of basis orbitals in the data set used in the
/// construction of the integrals. 


long   getOrbitalBasisCount() {return orbitalBasisCount;}


//
//  Returns the single particle energy of the basis orbitals.
//  The index takes on values from 1 to orbitalBasisCount
//
double getOrbitalEnergy(long i)
{
return orbitalEnergies[i];
}

//
// Single orbital access function
//
// Returns the value of the integral
//
// phi_i(r1)*H_singleParticle(r1)*phi_j(r1) dr1 
//
// where H_singleParticle(r1) is the single particle Hamiltonian.
//
// Each index takes on values from 0 to orbitalBasisCount -1
//
//
//
double operator()(long i, long j) const
{

#ifdef _DEBUG
   boundsCheck(i,orbitalBasisCount,1);
   boundsCheck(j,orbitalBasisCount,2);
#endif
   // shift indices by +1 because external indexing starts at 0

   if(i >= j)
   {
   return singleOrbitalIntegralData[getSingleOrbitalDataIndex(i+1,j+1)];
   }
   return singleOrbitalIntegralData[getSingleOrbitalDataIndex(j+1,i+1)];
}


//
// Double orbital access function.
//
// Returns the value of the integral
//
// phi_m(r1)*phi_n(r2)*G(r1,r2)*phi_i(r1)*phi_j(r2) dr1 dr2
//
// where G(r1,r2) is the Greens' function for the
// electrostatic potential operator and phi_* are single 
// particle spatial orbitals.
// 
// Each index takes on values from 0 to orbitalBasisCount-1
//
//
//
double operator()(long i, long j, long m, long n) const
{
    long iTmp; 

#ifdef _DEBUG
   boundsCheck(i,orbitalBasisCount,1);
   boundsCheck(j,orbitalBasisCount,2);
   boundsCheck(m,orbitalBasisCount,3);
   boundsCheck(n,orbitalBasisCount,4);
#endif

    
//  Create new indices that map into the data set 
//  (shifting by +1 because the external indexing starts at 0)

    long p = i+1; long q = m+1;  
    long r = j+1; long s = n+1; 
    
//
//  In the new indexing 
// 
//  (p,q) is associated with r1
//  (r,s) is associated with r2 
// 
//  The integral is obtained from the data file
//  with (p,q) and (r,s) ordered such that
//
//  p >= q  and r >= s, p >= r, and when p = r then q >= s
//

    // Re-order pairs so second index of each pair 
    // is less than the first

    if(p < q){iTmp = p; p = q; q = iTmp;} //  => (p >= q)
    if(r < s){iTmp = r; r = s; s = iTmp;} //  => (r >= s)

    // Re-order pairs so p >= r

    if(p < r) 
    {
     iTmp = p; p = r; r = iTmp;
     iTmp = q; q = s; s = iTmp;
    }

    // when p = r insure that q >= s

    if((p == r)&&(q < s)) {iTmp = q; q = s; s = iTmp;}

 
    return doubleOrbitalIntegralData[getDoubleOrbitalDataIndex(p, q, r, s)];
}
//
// Reads in data from file with name fileName
//
// Returns 0 : Data successfully read in 
//
// Returns 1 : File open failure or insufficient data in data file 
//
int inputOrbitalIntegralData(const char* fileName)
{
    FILE* Fin;

    if((Fin = fopen(fileName, "r")) == 0)
    {
    	string msg ="\nXXX QdomainExec XXX\n";
		msg       += "Orbital integrals data file specified \"" + string(fileName) + "\" not found.\n";
		throw std::runtime_error(msg);
    }

    // remove any previously read data 

    destroyData();          


    char     iLine[512];
    char  dateLine[512];
    char  fileLine[512];

    //
    // Read header 
    //

    char* fgetReturn;

    fgetReturn = fgets(&iLine[0],512,Fin); 
    fgetReturn = fgets(&iLine[0],512,Fin); 
    fgetReturn = fgets(&iLine[0],512,Fin); 
    fgetReturn = fgets(&iLine[0],512,Fin); 

    fgetReturn = fgets(&fileLine[0],512,Fin); 
    fgetReturn = fgets(&dateLine[0],512,Fin); 

    fgetReturn = fgets(&iLine[0],512,Fin); 
    fgetReturn = fgets(&iLine[0],512,Fin); 
    fgetReturn = fgets(&iLine[0],512,Fin); 


    // 
    // Get pointer to the file name and date strings 
    //
 
    char* fileLinePtr = strchr(fileLine,':') + 2;
    char* dateLinePtr = strchr(dateLine,':') + 2;

    strcpy(oIntegralProgramDataFile,fileLinePtr);
    strcpy(oIntegralProgramRunDate,dateLinePtr);

    //
    // Read in orbital integral data 
    //

    int  orbitalBasisCountInt;
    int  singleOrbitalIntegralCountInt;
    int  doubleOrbitalIntegralCountInt;
    
    int fscanfReturn;
 
    fscanfReturn = fscanf(Fin,"%d",&orbitalBasisCountInt);
    fscanfReturn = fscanf(Fin,"%d",&singleOrbitalIntegralCountInt);
    fscanfReturn = fscanf(Fin,"%d",&doubleOrbitalIntegralCountInt);
    

    orbitalBasisCount          = orbitalBasisCountInt;
    singleOrbitalIntegralCount = singleOrbitalIntegralCountInt;
    doubleOrbitalIntegralCount = doubleOrbitalIntegralCountInt;

    orbitalEnergies            = new double[orbitalBasisCount];
    singleOrbitalIntegralData  = new double[singleOrbitalIntegralCount];
    doubleOrbitalIntegralData  = new double[doubleOrbitalIntegralCount];

    long k;

    for(k = 0; k < orbitalBasisCount; k++)
    {
    fscanfReturn = fscanf(Fin,"%lf",&orbitalEnergies[k]);
    }

    for(k = 0; k < singleOrbitalIntegralCount; k++)
    {
    fscanfReturn = fscanf(Fin,"%lf",&singleOrbitalIntegralData[k]);
    }

    for(k = 0; k < doubleOrbitalIntegralCount; k++)
    {
    fscanfReturn = fscanf(Fin,"%lf",&doubleOrbitalIntegralData[k]);
    }

    // check for errors 

    if(feof(Fin) != 0) 
    {
    printf("##########   Error    #############\n"); 
	printf("Reading Integrals Data File \n\n"); 
	printf("Input File Specified: \n\n%s\n\n",fileName); 
	printf("###################################\n");
	return 1;
    }

    fclose(Fin);

    // To suppress unused variable warnings

    (void)fgetReturn;
    (void)fscanfReturn;

    return 0;
}


//
// Reads in binary data from file with the name fileName
//
// Returns 0 : Data successfully read in 
//
// Returns 1 : File open failure or insufficient data in data file 
//
int inputOrbitalIntegralDataBinary(const char* fileName)
{
    FILE* binaryDataFile = 0;
    if( (binaryDataFile = fopen(fileName, "rb" )) == 0 )
    {
    	string msg ="\nXXX QdomainExec XXX\n";
		msg       += "Orbital integrals data file specified \"" + string(fileName) + "\" not found.\n";
		throw std::runtime_error(msg);
	}

    // remove any previously read data 

    destroyData(); 

    long dataSize = 1;

    int  orbitalBasisCountInt;
    int  singleOrbitalIntegralCountInt;
    int  doubleOrbitalIntegralCountInt;
     
    int freadReturn;
 
    freadReturn = fread(&orbitalBasisCountInt, sizeof(int),          dataSize, binaryDataFile);
    freadReturn = fread(&singleOrbitalIntegralCountInt, sizeof(int), dataSize, binaryDataFile);
    freadReturn = fread(&doubleOrbitalIntegralCountInt, sizeof(int), dataSize, binaryDataFile);
    
    orbitalBasisCount           = orbitalBasisCountInt;
    singleOrbitalIntegralCount  = singleOrbitalIntegralCountInt;
    doubleOrbitalIntegralCount  = doubleOrbitalIntegralCountInt;
    
    //fread(&orbitalBasisCount, sizeof(long),          dataSize, binaryDataFile);
    //fread(&singleOrbitalIntegralCount, sizeof(long), dataSize, binaryDataFile);
    //fread(&doubleOrbitalIntegralCount, sizeof(long), dataSize, binaryDataFile);

    orbitalEnergies            = new double[orbitalBasisCount];
    singleOrbitalIntegralData  = new double[singleOrbitalIntegralCount];
    doubleOrbitalIntegralData  = new double[doubleOrbitalIntegralCount];

    freadReturn = fread(orbitalEnergies,           sizeof(double),orbitalBasisCount,binaryDataFile);
    freadReturn = fread(singleOrbitalIntegralData, sizeof(double),singleOrbitalIntegralCount,binaryDataFile);
    freadReturn = fread(doubleOrbitalIntegralData, sizeof(double),doubleOrbitalIntegralCount,binaryDataFile);

    // check for errors 

    if(feof(binaryDataFile) != 0) 
    {
    printf("##########   Error    #############\n"); 
	printf("Reading Integrals Data File \n\n"); 
	printf("Input File Specified: \n\n%s\n\n",fileName); 
	printf("###################################\n");
	return 1;
    }

    fclose(binaryDataFile);

    (void)freadReturn; // Suppress unused variable warning

    return 0;
}



void initialize(long orbitalBasisCount,double* orbitalEnergies,long singleOrbitalIntegralCount, 
double*  singleOrbitalIntegralData, long doubleOrbitalIntegralCount, double* doubleOrbitalIntegralData)
{
    destroyData();

    this->orbitalBasisCount           = orbitalBasisCount;
    this->singleOrbitalIntegralCount  = singleOrbitalIntegralCount;
    this->doubleOrbitalIntegralCount  = doubleOrbitalIntegralCount;

    this->orbitalEnergies            = new double[orbitalBasisCount];
    this->singleOrbitalIntegralData  = new double[singleOrbitalIntegralCount];

    if(this->doubleOrbitalIntegralCount > 0)
    {
    this->doubleOrbitalIntegralData  = new double[doubleOrbitalIntegralCount];
    }
    else
    {
    this->doubleOrbitalIntegralData  = 0;
    }

    long i;

    for(i = 0; i < orbitalBasisCount; i++) 
    {this->orbitalEnergies[i] = orbitalEnergies[i];}
    
    for(i = 0; i < singleOrbitalIntegralCount; i++) 
    {this->singleOrbitalIntegralData[i] = singleOrbitalIntegralData[i];}

    for(i = 0; i < doubleOrbitalIntegralCount; i++) 
    {this->doubleOrbitalIntegralData[i] = doubleOrbitalIntegralData[i];}
}
void outputOrbitalIntegralData(string& outputFileName)
{
    ostringstream s;
    FILE *Fout;

    struct tm when;
    time_t now; time(&now);
    when = *localtime( &now );

    Fout = fopen(outputFileName.c_str(), "w");
    fprintf(Fout,"%% \n");
    fprintf(Fout,"%% ############################################################## \n");
    fprintf(Fout,"%% \n");
    fprintf(Fout,"%% One and Two Electron Integrals Data File \n");
    fprintf(Fout,"%% Input File : %s \n","InternalDataConstruction");
    fprintf(Fout,"%% Run Date   : %s", asctime( &when ) );
    fprintf(Fout,"%% \n");
    fprintf(Fout,"%% ############################################################## \n");
    fprintf(Fout,"%% \n");

    fprintf(Fout,"%ld\n",orbitalBasisCount);
    fprintf(Fout,"%ld\n",singleOrbitalIntegralCount);
    fprintf(Fout,"%ld\n",doubleOrbitalIntegralCount);


    long k;

    for(k = 0; k < orbitalBasisCount; k++)
    {
    fprintf(Fout,"%20.15e\n",orbitalEnergies[k]);
    }

    for(k = 0; k < singleOrbitalIntegralCount; k++)
    {
    fprintf(Fout,"%20.15e\n",singleOrbitalIntegralData[k]);
    }

    for(k = 0; k < doubleOrbitalIntegralCount; k++)
    {
    fprintf(Fout,"%20.15e\n",doubleOrbitalIntegralData[k]);
    }

    fclose(Fout);

    printf("\n#######################################################\n\n");
    printf("Orbital Integral Data File (ascii format): %s \n",outputFileName.c_str());

    /*
    if ((selfInteractionFlag == 1)&&(nonLocalSelfInteractionFlag == 0))
    {
    printf("\nSingle Particle Self-Interaction Computation.\n\n");
    }
    else if ((selfInteractionFlag == 0)&&(nonLocalSelfInteractionFlag == 1))
    {
    printf("\nNon-local Single Particle Self-Interaction Computation.\n\n");
    }
    printf("\n#######################################################\n");
    */
}


void outputOrbitalIntegralDataBinary(string& outputFileName)
{
    ostringstream s;
    FILE *binaryDataFile;

    if( (binaryDataFile = fopen(outputFileName.c_str(), "wb" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",outputFileName.c_str());
      exit(1);
    }

    int  orbitalBasisCountInt          = orbitalBasisCount;
    int  singleOrbitalIntegralCountInt = singleOrbitalIntegralCount;
    int  doubleOrbitalIntegralCountInt = doubleOrbitalIntegralCount;

    fwrite(&orbitalBasisCountInt,            sizeof(int), 1, binaryDataFile);
    fwrite(&singleOrbitalIntegralCountInt,   sizeof(int), 1, binaryDataFile);
    fwrite(&doubleOrbitalIntegralCountInt,   sizeof(int), 1, binaryDataFile);

    //fwrite(&orbFunBasisCount,             sizeof(long), 1, binaryDataFile);
    //fwrite(&singleOrbitalIntegralCount,   sizeof(long), 1, binaryDataFile);
    //fwrite(&doubleOrbitalIntegralCount,   sizeof(long), 1, binaryDataFile);


    fwrite(&orbitalEnergies[0],sizeof(double),orbitalBasisCount,binaryDataFile);
    fwrite(&singleOrbitalIntegralData[0],     sizeof(double),singleOrbitalIntegralCount,binaryDataFile);
    fwrite(&doubleOrbitalIntegralData[0],     sizeof(double),doubleOrbitalIntegralCount,binaryDataFile);


    fclose(binaryDataFile);


    printf("\n#######################################################\n\n");
    printf("Orbital Integral Data File (binary format): %s   \n",outputFileName.c_str());
    printf("\n#######################################################\n");
    return;
    }

private:

//
// This private function returns the value of the data index for 
//
// phi_i(r1)*H_singleParticle(r1)*phi_j(r1) dr1 
//
// Note : This indexing is different from that used externally, e.g.
//        j must be less than or equal to i. This restriction is 
//        accounted for in the operator()(...) routine. 
//
#ifndef _DEBUG
inline long getSingleOrbitalDataIndex(long i, long j) const
{
    return i*(i-1)/2 + (j-1);
}
#else
long getSingleOrbitalDataIndex(long i, long j) const
{
    long index = i*(i-1)/2 + (j-1);
    boundsCheck(index,singleOrbitalIntegralCount,1);
    return index;
}
#endif

//
// This private function returns the value of the data index for 
//
// phi_m(r2)*phi_n(r2)*G(r1,r2)*phi_i(r1)*phi_j(r1) dr1 dr2
//
// Note : This indexing is different from that used externally
//        e.g. (i,j) is associated with r1 and and (m,n) is
//        associated with r2. This change in indexing is 
//        accounted for in the operator()(...) routine. 
//

#ifndef _DEBUG
inline long getDoubleOrbitalDataIndex(long i, long j, long m, long n) const
{
    return i*(i*i*i -2*i*i + 3*i - 2)/8 
           + ((j-1)*i*(i-1) + (j-1)*j)/2 
           + ((m-1)*m)/2 
           + (n-1);
}
#else
long getDoubleOrbitalDataIndex(long i, long j, long m, long n) const
{
    long index = i*(i*i*i -2*i*i + 3*i - 2)/8 
           + ((j-1)*i*(i-1) + (j-1)*j)/2 
           + ((m-1)*m)/2 
           + (n-1);
    boundsCheck(index,doubleOrbitalIntegralCount,1);
    return index;
}
#endif

void destroyData()
{
    if(orbitalEnergies != 0)           delete [] orbitalEnergies;
    if(doubleOrbitalIntegralData != 0) delete [] doubleOrbitalIntegralData;
    if(singleOrbitalIntegralData != 0) delete [] singleOrbitalIntegralData;

    orbitalEnergies          = 0;
    doubleOrbitalIntegralData = 0;
    singleOrbitalIntegralData = 0;
}


#ifdef _DEBUG 
    void boundsCheck(long i, long Size, int coordinate) const
    {
    if((i < 0)||(i  >= Size))
    {
    printf("RealSpatialOribital index %d out of bounds \n",coordinate);
    printf("Offending index value %ld : Acceptable Range [0, %ld] \n",i, Size-1);
    exit(1);
    }
    }
#else
void boundsCheck(long, long, int) const {}
#endif

public:

    long             orbitalBasisCount;
    double*            orbitalEnergies;

    long    singleOrbitalIntegralCount;
    double*  singleOrbitalIntegralData;

    long    doubleOrbitalIntegralCount;
    double*  doubleOrbitalIntegralData;

    char oIntegralProgramDataFile[512];
    char  oIntegralProgramRunDate[512];
};

} // UCLAQ namespace
#endif




 
