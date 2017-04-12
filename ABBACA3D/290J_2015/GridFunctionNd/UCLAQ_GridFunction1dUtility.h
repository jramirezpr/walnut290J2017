/*
 * UCLAQ_GridFunction1dUtility.h
 *
 *  Created on: Jun 28, 2015
 *      Author: anderson
 */

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

#ifndef _UCLAQ_GridFunction1dUtility_
#define _UCLAQ_GridFunction1dUtility_


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
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


#include "UCLAQ_GridFunction1d.h"


namespace UCLAQ
{

class GridFunction1dUtility
{
public:

void outputToGNUplot(const GridFunction1d& gF, const string& fileName, const string& formatString = "%15.10e")
{
//
//  Open and then write to a file
//
    ostringstream s;
    FILE* dataFile;

    s.str("");
    s << formatString << "  " << formatString << " \n";

    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }

    long i;

    double* xp;

    double hx     = gF.getHx();
    long mPanel   = gF.getXpanelCount();
    double xMin   = gF.getXmin();

    xp     = new double[mPanel+1];

    for(i = 0; i <= mPanel; i++)
    {
    xp[i] = xMin + i*hx;
    }

    for(i = 0;  i <= mPanel; i++)
    {
    fprintf(dataFile,(s.str()).c_str(),xp[i], gF.Values(i));
    }

    delete [] xp;

    fclose(dataFile);
}

void outputToMatlab(const GridFunction1d& gF,  const string& fileName, const string& formatString = "%15.10e")
{
//
//  Open and then write to a file
//
    FILE* dataFile;
    ostringstream s;
    s.str("");
    s << formatString << "  " << formatString << " \n";

    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }

    long i;

    double xp;

    double hx     = gF.getHx();
    long mPanel   = gF.getXpanelCount();
    double xMin   = gF.getXmin();

	fprintf(dataFile,"%ld \n",mPanel+1);

	// print out x values

    for(i = 0;  i <= mPanel; i++)
    {
    xp = xMin + i*hx;
    fprintf(dataFile,formatString.c_str(),xp);
	fprintf(dataFile,"\n");
    }

	// print out y values

	for(i = 0;  i <= mPanel; i++)
    {
    fprintf(dataFile,formatString.c_str(),gF.Values(i));
	fprintf(dataFile,"\n");
    }

    fclose(dataFile);
}

void outputFunction(const GridFunction1d& fun, const string& fileName, const string& outputFormat,
const string& formatString = "%15.10e")
{
    ostringstream s;
    s.str("");

    if(outputFormat.compare("GNUPLOT")==0)
    {
    s << fileName << ".dat";
    outputToGNUplot(fun, (s.str()).c_str(),formatString);
    }
    else if(outputFormat.compare("MATLAB")==0)
    {
    s << fileName << ".dat";
    outputToMatlab(fun, (s.str()).c_str(),formatString);
    }
	else if(outputFormat.compare("MATLAB_PLAIN")==0)
    {
    s << fileName << ".dat";
    outputToMatlab(fun, (s.str()).c_str(),formatString);
    }
    else
    {
    s << fileName << ".dat";
    outputToGNUplot(fun, (s.str()).c_str(),formatString);
    }
}


void appendToGNUplot(const GridFunction1d& gF,const string& fileName, const string& formatString = "%15.10e")
{
//
//  Open and then write to a file
//
    FILE* dataFile;
	ostringstream s;

    s.str("");
    s << formatString << "  " << formatString << " \n";


    if(OPENFILE(dataFile,fileName.c_str(), "a+" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }

    long i;

    double* xp;

    double hx     = gF.getHx();
    long mPanel   = gF.getXpanelCount();
    double xMin   = gF.getXmin();

    xp     = new double[mPanel+1];

    for(i = 0; i <= mPanel; i++)
    {
    xp[i] = xMin + i*hx;
    }

    fprintf(dataFile,"\n");
    for(i = 0;  i <= mPanel; i++)
    {
    fprintf(dataFile,(s.str()).c_str(),xp[i], gF.Values(i));
    }

    delete [] xp;

    fclose(dataFile);
}

void inputFromGNUplot(GridFunction1d& gF, const string& fileName, int& noFileFlag)
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

    UCLAQ::DoubleVector1d* V;

    long i;
    long mPanel;
    double x;

    V      = gF.getValuesPointer();
    mPanel = gF.getXpanelCount();
    for(i = 0;  i <= mPanel; i++)
    {
	fscanf(dataFile,"%lf",&x);
    fscanf(dataFile,"%lf",&(V->operator()(i)));
    }
	fclose(dataFile);
}

//
// Outputs the data to a file in gnuplot format, and then creates a very basic Veusz plot file that
// is linked to that data.
//
void outputToVeusz(const GridFunction1d& gF, const string& fileName, const string& formatString = "%15.10e")
{
//	 Output data file in gnuplot format

	outputToGNUplot(gF, fileName,formatString);

	int lastindex    = fileName.find_last_of(".");
    string baseName  = fileName.substr(0, lastindex);
    string veuszName = baseName;
    veuszName.append(".vsz");

    FILE* dataFile;
    if(OPENFILE(dataFile,fileName.c_str(), "w" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }

	ostringstream s;
    s.str("");

	fprintf(dataFile,"%s","#Veusz saved document (version 1.23.1)\n");
	fprintf(dataFile,"%s","AddImportPath(u'./')\n");

	s << "ImportFile(u'./" << fileName << "', u'', ignoretext=True, linked=True, prefix=u'" << baseName << "')\n";
	fprintf(dataFile,"%s",(s.str()).c_str());

	fprintf(dataFile,"%s","Add('page', name='page1', autoadd=False)\n");
	fprintf(dataFile,"%s","To('page1')\n");
	fprintf(dataFile,"%s","Add('graph', name='graph1', autoadd=False)\n");
	fprintf(dataFile,"%s","To('graph1')\n");
	fprintf(dataFile,"%s","Add('axis', name='x', autoadd=False)\n");
	fprintf(dataFile,"%s","Add('axis', name='y', autoadd=False)\n");
	fprintf(dataFile,"%s","To('y')\n");
	fprintf(dataFile,"%s","Set('direction', 'vertical')\n");
	fprintf(dataFile,"%s","To('..')\n");
	fprintf(dataFile,"%s","Add('xy', name='xy1', autoadd=False)\n");
	fprintf(dataFile,"%s","To('xy1')\n");

	s.str("");
	s << "Set('xData', u'" << baseName << "1" << "')\n";
	fprintf(dataFile,"%s",(s.str()).c_str());
	s.str("");
	s << "Set('yData', u'" << baseName << "2" << "')\n";
	fprintf(dataFile,"%s",(s.str()).c_str());

	fprintf(dataFile,"%s","Set('markerSize', u'1.5pt')\n");
	fprintf(dataFile,"%s","To('..')\n");
    fprintf(dataFile,"%s","To('..')\n");
    fprintf(dataFile,"%s","To('..')\n");

    fclose(dataFile);
}

//
// Veusz plot used to determine veusz plot data structure for the above member function
//
/*
# Veusz saved document (version 1.23.1)
# Saved at 2015-06-29T22:42:43.183012

AddImportPath(u'./')
ImportFile(u'./Mollifier.dat', u'', ignoretext=True, linked=True, prefix=u'M')
Add('page', name='page1', autoadd=False)
To('page1')
Add('graph', name='graph1', autoadd=False)
To('graph1')
Add('axis', name='x', autoadd=False)
Add('axis', name='y', autoadd=False)
To('y')
Set('direction', 'vertical')
To('..')
Add('xy', name='xy1', autoadd=False)
To('xy1')
Set('xData', u'M1')
Set('yData', u'M2')
Set('markerSize', u'1.5pt')
To('..')
To('..')
To('..')
 */
};
} // namespace


#undef  OPENFILE
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif /* UCLAQ_GRIDFUNCTION_UCLAQ_GRIDFUNCTION1DUTILITY_H_ */
