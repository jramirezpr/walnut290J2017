
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

#ifndef _UCLAQ_GridFunction2dUtility_
#define _UCLAQ_GridFunction2dUtility_
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <functional>
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

#include "../GridFunctionNd/UCLAQ_GridFunction2d.h"
#include "../DoubleVectorNd/UCLAQ_DoubleVector2d.h"

namespace UCLAQ
{
class  GridFunction2dUtility
{
public :

// Adds the values of F to the input gFun

void addToValues(GridFunction2d& gFun, std::function<double(double,double)>& F)
{
    long i; long j;

    double xMin     = gFun.getXmin();
    double yMin     = gFun.getYmin();

    long xPanels = gFun.getXpanelCount();
    long yPanels = gFun.getYpanelCount();

    double  hx   = gFun.getHx();
    double  hy   = gFun.getHy();

    double x; double y;

    for(i = 0; i <= xPanels; i++)
    {
    x = xMin + i*hx;
    for(j = 0; j <= yPanels; j++)
    {
    y = yMin + j*hy;
    gFun.Values(i,j) += F(x,y);
    }}
}


void outputToGNUplot(GridFunction2d& gF, const string& fileName, const string& formatString = "%15.10e")
{
    ostringstream s;
    s.str("");
    s << formatString << "  " << formatString << "  " << formatString << " \n";
//
//  Open and then write to a file
//
    FILE* dataFile;

    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }

    UCLAQ::DoubleVector2d* V = gF.getValuesPointer();

    double hx     = gF.getHx();
    long mPanel    = gF.getXpanelCount();
    double xMin   = gF.getXmin();

    double hy     = gF.getHy();
    long nPanel   = gF.getYpanelCount();
    double yMin   = gF.getYmin();

    long i; long j;

    fprintf(dataFile,"# %ld %ld\n",mPanel+1, nPanel+1);

    double x;
    double y;

    for(i = 0;  i <= nPanel; i++)
    {
    //for(j = mPanel+1; j >= 1;   j--)
   // {
   // x = xMin + double(mPanel+ 1 - j)*hx;
	for(j = 0; j <= mPanel; j++)
	{
	x = xMin + double(j)*hx;
    y = yMin + double(i)*hy;
    fprintf(dataFile,(s.str()).c_str(),x,y,V->operator()(j,i));
    }
    fprintf(dataFile,"\n");
    }

    fclose(dataFile);
}

//
// In the case that the input GridFun2d is a null instance as indicated by
// mPanel=nPanel=0, this routine initializes gF based upon data contained in
// fileName.

 void inputFromGNUplot(GridFunction2d& gF, const string& fileName, int& noFileFlag)
{

	//
	//  Open input file
	//
	// ifstream dataFile(fileName);

    FILE* dataFile = 0;
	noFileFlag     = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "r" ))
    {
	  noFileFlag = 1;
      return;
    }

    char* strTmp = new char[1024];
    long mPtmp;
    long nPtmp;

    long i;      long j;
    long mPanel; long nPanel;
    double x;    double y;

    double xMin; double xMax;
    double yMin; double yMax;
    double valTmp;

    if((gF.getXpanelCount() == 0)&&(gF.getYpanelCount() == 0))
    {

    // Extract panel count

    fscanf(dataFile,"%s",strTmp);
    fscanf(dataFile,"%ld",&mPanel);
    fscanf(dataFile,"%ld",&nPanel);
    mPanel -= 1;
    nPanel -= 1;
    //
    //  Capture bounds
    //
    //	x = xMin + double(j)*hx;
    //  y = yMin + double(i)*hy;

    for(i = 0;  i <= nPanel; i++)
    {
	for(j = 0;  j <= mPanel; j++)
	{
	fscanf(dataFile,"%lf %lf %lf",&x,&y,&valTmp);
	if(j == 0)      xMin = x;
    if(j == nPanel) xMax = x;
    }
	if(i == 0)      yMin = y;
	if(i == mPanel) yMax = y;
    }
    gF.initialize(mPanel,xMin,xMax,nPanel,yMin,yMax);
    rewind(dataFile);
    }

    UCLAQ::DoubleVector2d* V;
    V = gF.getValuesPointer();

    mPanel = gF.getXpanelCount();
    nPanel = gF.getYpanelCount();

    fscanf(dataFile,"%s",strTmp);
    fscanf(dataFile,"%ld",&mPtmp);
    fscanf(dataFile,"%ld",&nPtmp);

    for(i = 0;  i <= nPanel; i++)
    {
	for(j = 0; j <= mPanel; j++)
	{
	fscanf(dataFile,"%lf %lf %lf",&x,&y,&(V->operator()(j,i)));
    }
    }
    delete [] strTmp;
    fclose(dataFile);
}

};

}

#undef OPENFILE
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#endif




 
