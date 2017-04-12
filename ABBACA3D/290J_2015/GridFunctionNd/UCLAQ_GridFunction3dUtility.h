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
#ifndef _UCLAQ_GridFunction3dUtility_
#define _UCLAQ_GridFunction3dUtility_
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <functional>
#include <cassert>
using namespace std;

// MS compilers generate warnings if fopen is used instead of fopen_s (a Microsoft specific language
// extension, so this macro implements the appropriate equivalent to fopen that MS wants when
// MS compilers are being used. In both versions, the
// macro returns a non-zero value if the open fails (e.g. a non-zero error code).
//
#ifndef _MSC_VER
#define OPENFILE(dataFile,filename,mode) ((dataFile = fopen(filename,  mode)) == NULL)
#else
#define OPENFILE(dataFile,fileName,mode) ((fopen_s(&dataFile,fileName, mode)) != 0)
#pragma warning(push)
#pragma warning(disable: 4996)
#endif


#include "../DoubleVectorNd/UCLAQ_DoubleVector3d.h"
#include "../GridFunctionNd/UCLAQ_GridFunction3d.h"


namespace UCLAQ
{
class  GridFunction3dUtility
{
public :

//
// If xScalingFactor > 0 then the x-direction (vertical) is scaled so that the
// the extent of the x-coordinates span a distance of xScalingFactor times the
// maximal transverse coordinate span. Typical factors for this are .25, .3, etc
// enough to give an idea that the vertical layer is thin with respect to the
// transverse size, but large enough so that it is easy to visualize in 3D.
//
//
void outputDataToVTKfile(const GridFunction3d& gridFun, const string& fileName, const string& dataLabel,
double xScalingFactor = -1.0)
{
	FILE* dataFile;

    double a  = gridFun.getXmin();  double b  = gridFun.getXmax();
    double c  = gridFun.getYmin();  double d  = gridFun.getYmax();
	double e  = gridFun.getZmin();  double f  = gridFun.getZmax();

    double hx = gridFun.getHx();
    double hy = gridFun.getHy();
	double hz = gridFun.getHz();

    long mPt = gridFun.getXpanelCount() + 1;
    long nPt = gridFun.getYpanelCount() + 1;
	long pPt = gridFun.getZpanelCount() + 1;

    long dataCount = mPt*nPt*pPt;

    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      exit(1);
    }

    double hxScaling = 1.0;
    double transverseSizeMax;

    if(xScalingFactor > 0.0)
    {
    transverseSizeMax = f-e;
    transverseSizeMax = (transverseSizeMax > (d-c)) ? transverseSizeMax : (d-c);

    hxScaling  = (xScalingFactor*transverseSizeMax)/(b-a);
    }

    long i; long j; long k;
    double xPos; double yPos; double zPos;
    //
    // output the regular positions
    //
	fprintf(dataFile, "# vtk DataFile Version 2.0\n");
    fprintf(dataFile, "Phi Sample \n");
    fprintf(dataFile, "ASCII\n");

    fprintf(dataFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(dataFile, "DIMENSIONS %ld %ld %ld \n",mPt,nPt,pPt);
    fprintf(dataFile, "X_COORDINATES %ld float \n",mPt);
    for(i = 0; i < mPt; i++)
    {
    xPos = i*hx*hxScaling + a*hxScaling;
    fprintf(dataFile, "%10.5e ",xPos);
    }
    fprintf(dataFile, "\n");
    fprintf(dataFile, "Y_COORDINATES %ld float \n",nPt);
    for(j = 0; j < nPt; j++)
    {
    yPos = j*hy + c;
    fprintf(dataFile, "%10.5e ",yPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "Z_COORDINATES %ld float \n",pPt);
    for(k = 0; k < pPt; k++)
    {
    zPos = k*hz + e;
    fprintf(dataFile, "%10.5e ",zPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "POINT_DATA %ld\n",dataCount);
    fprintf(dataFile, "SCALARS %s float\n",dataLabel.c_str());
    fprintf(dataFile, "LOOKUP_TABLE default\n");
    for(k = 0; k <  pPt; k++)
    {
    for(j = 0; j < nPt; j++)
    {
    for(i = 0; i < mPt; i++)
    {
    fprintf(dataFile, "%15.8e ",gridFun.Values(i,j,k));
    }
    fprintf(dataFile, "\n");
    }}

    fclose(dataFile);

}

void outputToDataFile(const GridFunction3d& gF, const string& fileName, const string& formatString = "%15.10e")
{
//
//  Create format string
//
    ostringstream s;
    s.str("");
    s << formatString << " ";
//
//  Open and then write to a file
//
    FILE* dataFile;

    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }

    long xPt = gF.getXpanelCount() + 1;
    long yPt = gF.getYpanelCount() + 1;
	long zPt = gF.getZpanelCount() + 1;

    double xMin  = gF.getXmin();
    double yMin  = gF.getYmin();
	double zMin  = gF.getZmin();

    double xMax    = gF.getXmax();
    double yMax    = gF.getYmax();
    double zMax    = gF.getZmax();


    fprintf(dataFile,"%ld \n", xPt);
	fprintf(dataFile,"%ld \n", yPt);
	fprintf(dataFile,"%ld \n", zPt);

    fprintf(dataFile,"%15.10e \n",xMin);
	fprintf(dataFile,"%15.10e \n",xMax);
	fprintf(dataFile,"%15.10e \n",yMin);
	fprintf(dataFile,"%15.10e \n",yMax);
    fprintf(dataFile,"%15.10e \n",zMin);
	fprintf(dataFile,"%15.10e \n",zMax);


    for(long k = 0; k < zPt; k++)
    {
    for(long j = 0; j < yPt; j++)
    {
    for(long i = 0; i < xPt; i++)
    {
    fprintf(dataFile, s.str().c_str(),gF.Values(i,j,k));
    }
    fprintf(dataFile, "\n");
    }}


    fclose(dataFile);
}

void inputFromDataFile(GridFunction3d& gF, FILE* dataFile)
{
    long xPt;
    long yPt;
	long zPt;

    double xMin;
    double yMin;
	double zMin;

    double xMax;
    double yMax;
    double zMax;



    fscanf(dataFile,"%ld", &xPt);
	fscanf(dataFile,"%ld", &yPt);
	fscanf(dataFile,"%ld", &zPt);

    fscanf(dataFile,"%lf",&xMin);
	fscanf(dataFile,"%lf",&xMax);
	fscanf(dataFile,"%lf",&yMin);
	fscanf(dataFile,"%lf",&yMax);
    fscanf(dataFile,"%lf",&zMin);
	fscanf(dataFile,"%lf",&zMax);

    gF.initialize(xPt-1,xMin,xMax,yPt-1,yMin,yMax,zPt-1,zMin,zMax);


	for(long k = 0; k < zPt; k++)
    {
    for(long j = 0; j < yPt; j++)
    {
    for(long i = 0; i < xPt; i++)
    {
    fscanf(dataFile,"%lf",&gF.Values(i,j,k));
    }
    }}

}

void inputFromDataFile(GridFunction3d& gF, const string& fileName)
{
//
//  Open input file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "r" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    inputFromDataFile(gF, dataFile);

	fclose(dataFile);
}


void outputToBinaryDataFile(const GridFunction3d& gF, FILE* dataFile)
{
    long dataSize;

    long xPt = gF.getXpanelCount() + 1;
    long yPt = gF.getYpanelCount() + 1;
	long zPt = gF.getZpanelCount() + 1;

    double xMin  = gF.getXmin();
    double yMin  = gF.getYmin();
	double zMin  = gF.getZmin();

    double hx = gF.getHx();
    double hy = gF.getHy();
	double hz = gF.getHz();

	//
	//  Write out the grid structure information
	//

	//
	//
	// To avoid problems with nonstandard sizes of integers in binary format
	// output the dimensions as real values with a small epsilon to insure
	// the conversion back to integers results in the correct integer value
	//

	double integerEps = 1.0e-07;

	double doubleXpt  = (double)xPt + integerEps;
	double doubleYpt  = (double)yPt + integerEps;
	double doubleZpt  = (double)zPt + integerEps;

    fwrite(&doubleXpt,  sizeof(double), 1, dataFile);
	fwrite(&doubleYpt,  sizeof(double), 1, dataFile);
	fwrite(&doubleZpt,  sizeof(double), 1, dataFile);

    fwrite(&hx,  sizeof(double), 1, dataFile);
	fwrite(&hy,  sizeof(double), 1, dataFile);
	fwrite(&hz,  sizeof(double), 1, dataFile);

    fwrite(&xMin,  sizeof(double), 1, dataFile);
	fwrite(&yMin,  sizeof(double), 1, dataFile);
	fwrite(&zMin,  sizeof(double), 1, dataFile);
//
//  Write ot the function values
//
    dataSize = xPt*yPt*zPt;
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}

void inputFromBinaryDataFile(GridFunction3d& gF, const string& fileName,
		   int& noFileFlag)
{
	//
	//  Open input file (remember to use the b mode to specify binary!!!!)
	//
	FILE* dataFile = 0;
	noFileFlag     = 0;
	if(OPENFILE(dataFile,fileName.c_str(), "rb" ))
	{
		  noFileFlag = 1;
	      return;
	}

	inputFromBinaryDataFile(gF,dataFile);

    if(ferror(dataFile))
	{
	      printf( "GridFunction3d could not be initialized from data file  %s \n",fileName.c_str());
	      exit(1);
	}

    fclose(dataFile);
}
void inputFromBinaryDataFile(GridFunction3d& gF, FILE* dataFile)
{
    long dataSize;

    long    xPt;    long yPt;    long zPt;
    double   hx;   double hy;   double hz;
    double xMin; double yMin; double zMin;
    double xMax; double yMax; double zMax;

	double doubleXpt;
	double doubleYpt;
	double doubleZpt;

	fread(&doubleXpt,  sizeof(double), 1, dataFile);
	fread(&doubleYpt,  sizeof(double), 1, dataFile);
	fread(&doubleZpt,  sizeof(double), 1, dataFile);

	fread(&hx,  sizeof(double), 1, dataFile);
	fread(&hy,  sizeof(double), 1, dataFile);
	fread(&hz,  sizeof(double), 1, dataFile);

	fread(&xMin,  sizeof(double), 1, dataFile);
	fread(&yMin,  sizeof(double), 1, dataFile);
	fread(&zMin,  sizeof(double), 1, dataFile);

	xPt = (long)doubleXpt;
	yPt = (long)doubleYpt;
	zPt = (long)doubleZpt;

	xMax = xMin + (xPt-1)*hx;
	yMax = yMin + (yPt-1)*hy;
	zMax = zMin + (zPt-1)*hz;

	gF.initialize(xPt-1,xMin,xMax,yPt-1,yMin,yMax,zPt-1,zMin,zMax);
	dataSize = xPt*yPt*zPt;
	fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}


void inputValuesFromBinaryDataFile(GridFunction3d& gF, FILE* dataFile)
{
    long dataSize;

    long xPt = gF.getXpanelCount() + 1;
    long yPt = gF.getYpanelCount() + 1;
	long zPt = gF.getZpanelCount() + 1;
	dataSize = xPt*yPt*zPt;
	fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}


void appendValuesToBinaryDataFile(const GridFunction3d& gF, FILE* dataFile)
{
	long dataSize;
    long xPt = gF.getXpanelCount() + 1;
    long yPt = gF.getYpanelCount() + 1;
	long zPt = gF.getZpanelCount() + 1;
//
//  Write ot the function values
//
    dataSize = xPt*yPt*zPt;
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}

//
// This routine opens up a new file and write GridFunction3d structure and data
// and then closes the file
//
void outputToBinaryDataFile(const GridFunction3d& gF, const string& fileName)
{
//
//  Open and then write to a file (remember to use the b mode to specify binary!!!!)
//
    FILE* dataFile;

    if(OPENFILE(dataFile,fileName.c_str(), "wb" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }
    outputToBinaryDataFile(gF, dataFile);
    fclose(dataFile);
}

/*

void GridFunction3dUtility::getDomainData(XML_ParameterListArray& paramList,
                      long& xPanels, double& xMin, double& xMax,
                      long& yPanels, double& yMin, double& yMax,
                      long& zPanels, double& zMin, double& zMax)
{
	if(paramList.isParameterList("ComputationalDomain") == 0)
    {
    string mesg = "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    mesg       += "Error initializing GridFunction3d class \n";
    mesg       += "ComputationalDomain parameter list was not found \n";
    mesg       += "in input XML_ParameterListArray \n";
    mesg       += "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    throw runtime_error(mesg);
    }

    if(paramList.isParameterList("GridParameters") == 0)
    {
    string mesg = "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    mesg       += "Error initializing GridFunction3d class \n";
    mesg       += "GridParameters parameter list was not found \n";
    mesg       += "in input XML_ParameterListArray \n";
    mesg       += "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    throw runtime_error(mesg);
    }

    xPanels  = paramList.getParameterValue("xPanels","GridParameters");
    yPanels  = paramList.getParameterValue("yPanels","GridParameters");
    zPanels  = paramList.getParameterValue("zPanels","GridParameters");

    xMin  = paramList.getParameterValue("xMin","ComputationalDomain");
	yMin  = paramList.getParameterValue("yMin","ComputationalDomain");
	zMin  = paramList.getParameterValue("zMin","ComputationalDomain");
	xMax  = paramList.getParameterValue("xMax","ComputationalDomain");
	yMax  = paramList.getParameterValue("yMax","ComputationalDomain");
	zMax  = paramList.getParameterValue("zMax","ComputationalDomain");
}

*/

};
}

#undef OPENFILE
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif
 
