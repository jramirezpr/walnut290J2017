#include <iostream>
#include <cstdio>
using namespace std;

#include "XML_ParameterList/XML_ParameterListArray.h"

//
// A sample program demonstrating the construction  
// and use of an XML_ParameterListArray.
//
// The test program assumes the directories 290J_2015 and 290J_Samples are organized as follows
//
// Course Directory
// ---------------------------------------------------------------
//      |                |                |            |
// Project_1 Project_2   ....        290J_Samples   290J_2015
//
//
// where the directory 290J_2015 contains the source directories :
//
//    XML_ParameterList
//
// The command line compilation command is
//
// g++ ParameterListTest.cpp  ../../290J_2015/XML_ParameterList/*.cpp  -I../../290J_2015  -o ParameterListTest.exe
//
// or, if one has created an XML_ParameterList library (see below)
//
// g++ ParameterListTest.cpp   -I../../290J_2015 -L../../290J_2015/XML_ParameterList/lib -lXML_ParameterList -o ParameterListTest.exe
//
//
// If using the supplied makefiles (*.mk),  one first creates the
// the XML_ParameterList library file by cd'ing to the XML_ParameterList
// directory and executing
//
// make -f XML_ParameterListLib.mk release
//
// This is then followed by cd'ing to 290J_Samples/XML_ParameterListTest and executing
//
// make -f XML_ParameterList.mk release
//
// Or if one wants to use the makefile that both builds the library and then the executable one
// uses the command
//
// make release
//
//
// Created for Math 290J UCLA Dept. of Mathematics
//
// Chris Anderson : Oct. 28, 2015
//
int main()
{
	XML_ParameterListArray paramList("ParameterListTest.xml");

    // Set to abort if an error is made using the XML_ParameterListArray

    paramList.setAbortOnErrorFlag();

    // Echo X domain size 

    printf("X-Domain Size \n");
    
    double xMin = paramList.getParameterValue("xMin","ComputationalDomain");
    double xMax = paramList.getParameterValue("xMax","ComputationalDomain");

    printf("xMin = %10.5f\n", xMin);
    printf("xMax = %10.5f\n", xMax);
    printf("\n");

    //
    // Reset panel counts if internal instance and then print out 
    //

    paramList.setParameter(10,"xPanels","GridParameters");
    paramList.setParameter( 7,"yPanels","GridParameters");
    paramList.setParameter(11,"zPanels","GridParameters");

    cout << paramList << endl;


    // Clear abort on error and then detect the fetch of a
    // non-existent data item by checking for a null instance
    // as the type of the XML_dataType object returned.

    paramList.clearAbortOnErrorFlag();

    XML_dataType paramValue = paramList.getParameterValue("minimumX","ComputationalDomain");

    if(paramValue.isNull())
    {
    printf("Attempt to access non-existent element \n");
    printf("%s",paramList.errorMessage.c_str());
    }


	return 0;
}
