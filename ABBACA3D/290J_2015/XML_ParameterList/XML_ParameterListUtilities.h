/*
 * XML_ParameterListUtilities.h
 *
 *  Created on: Mar 22, 2013
 *      Author: anderson
 */
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

#ifndef _XML_ParameterListUtilities_
#define _XML_ParameterListUtilities_

#include "XML_ParameterListArray.h"

class XML_ParameterListUtilities
{
public:

XML_ParameterListUtilities(){};


//
// This member function takes the input paramList and inserts all parameterLists
// whose name is not an element of the exclusions vector of names into
// the outputList.
//
// The outputList must be initialized before calling this routine, or
// created with an call to createParameterListArray(...)
//

void createDuplicateEntries(XML_ParameterListArray& paramList,
vector<string>& exclusions,    XML_ParameterListArray& outputList)
{
    vector<string> paramListNames;
    vector<string> paramNames;
    vector<string> paramChildNames;
    bool excludeFlag;
    bool abortOnErrorFlag  = paramList.abortOnErrorFlag;
    XML_dataType XMLvalue;

    vector<string>::iterator it;

    long      paramInstanceIndex;
    paramList.getParameterListNames(paramListNames);
    for(long k = 0; k < (long)paramListNames.size(); k++)
    {

    excludeFlag = true;
    it = std::find(exclusions.begin(),exclusions.end(),paramListNames[k]);
    if(it == exclusions.end()){excludeFlag = false;}

    if(not excludeFlag)
    {

    outputList.addParameterList(paramListNames[k].c_str());

    paramNames.clear();
    paramList.getParameterNames(paramListNames[k].c_str(),paramNames);

    // Remove duplicate parameter names

    std::sort(paramNames.begin(),paramNames.end());
    it = std::unique(paramNames.begin(), paramNames.end());
    paramNames.resize(it - paramNames.begin());

    for(long j = 0; j < (long)paramNames.size(); j++)
    {
    for(paramInstanceIndex = 0;
    paramInstanceIndex < (long)paramList.getParameterInstanceCount(paramNames[j].c_str(),paramListNames[k].c_str());
    paramInstanceIndex++)
    {

    paramChildNames.clear();
    paramList.getParameterChildNames(paramInstanceIndex, paramNames[j].c_str(),  paramListNames[k].c_str(), paramChildNames);
    if(paramChildNames.size() == 0)
    {
    outputList.addParameter(paramList.getParameterValue(paramNames[j].c_str(),paramListNames[k].c_str()),
                            paramNames[j].c_str(), paramListNames[k].c_str());
    }
    else // add parameter and child parameter values
    {
    outputList.addParameter(paramNames[j].c_str(), paramListNames[k].c_str());


    // Add parameter value if it is specified

    paramList.clearAbortOnErrorFlag();
    XMLvalue = paramList.getParameterValue(paramNames[j].c_str(), paramListNames[k].c_str());
    if(not XMLvalue.isNull())
    {outputList.setParameter(XMLvalue,paramNames[j].c_str(), paramListNames[k].c_str());}
    else
    {paramList.clearError();}
    if(abortOnErrorFlag) paramList.setAbortOnErrorFlag();


    // Add child parameter values

    for(long i = 0; i < (int)paramChildNames.size(); i++)
    {
    outputList.addParameterInstanceChild(paramList.getParameterInstanceChildValue(paramInstanceIndex,
                                         paramChildNames[i].c_str(),
                                         paramNames[j].c_str(),
                                         paramListNames[k].c_str()),
                                         paramInstanceIndex,
                                         paramChildNames[i].c_str(), paramNames[j].c_str(), paramListNames[k].c_str());
    }

    }}

    } // parameter instances

    }}
}

//
// This routine copies the parameter list from the input XML_ParameterListArray. An
// error message is generated if the parameter does not exist.
//
void copyParameterList(XML_ParameterListArray& paramListArray,const char* paramList, XML_ParameterListArray& outputList)
{
    bool errorFlag = true;

    vector<string> paramListArrayNames;
    vector<string> paramNames;
    vector<string> paramChildNames;
    bool abortOnErrorFlag  = paramListArray.abortOnErrorFlag;
    XML_dataType XMLvalue;

    vector<string>::iterator it;

    long      paramInstanceIndex;
    paramListArray.getParameterListNames(paramListArrayNames);
    for(long k = 0; k < (long)paramListArrayNames.size(); k++)
    {

    if(paramListArrayNames[k].compare(paramList) == 0)
    {
    errorFlag = false;
    outputList.addParameterList(paramListArrayNames[k].c_str());

    paramNames.clear();
    paramListArray.getParameterNames(paramListArrayNames[k].c_str(),paramNames);

    // Remove duplicate parameter names

    std::sort(paramNames.begin(),paramNames.end());
    it = std::unique(paramNames.begin(), paramNames.end());
    paramNames.resize(it - paramNames.begin());

    for(long j = 0; j < (long)paramNames.size(); j++)
    {
    for(paramInstanceIndex = 0;
    paramInstanceIndex < (long)paramListArray.getParameterInstanceCount(paramNames[j].c_str(),paramListArrayNames[k].c_str());
    paramInstanceIndex++)
    {

    paramChildNames.clear();
    paramListArray.getParameterChildNames(paramInstanceIndex, paramNames[j].c_str(),  paramListArrayNames[k].c_str(), paramChildNames);
    if(paramChildNames.size() == 0)
    {
    outputList.addParameter(paramListArray.getParameterValue(paramNames[j].c_str(),paramListArrayNames[k].c_str()),
                            paramNames[j].c_str(), paramListArrayNames[k].c_str());
    }
    else // add parameter and child parameter values
    {
    outputList.addParameter(paramNames[j].c_str(), paramListArrayNames[k].c_str());


    // Add parameter value if it is specified

    paramListArray.clearAbortOnErrorFlag();
    XMLvalue = paramListArray.getParameterValue(paramNames[j].c_str(), paramListArrayNames[k].c_str());
    if(not XMLvalue.isNull())
    {outputList.setParameter(XMLvalue,paramNames[j].c_str(), paramListArrayNames[k].c_str());}
    else
    {paramListArray.clearError();}
    if(abortOnErrorFlag) paramListArray.setAbortOnErrorFlag();


    // Add child parameter values

    for(long i = 0; i < (int)paramChildNames.size(); i++)
    {
    outputList.addParameterInstanceChild(paramListArray.getParameterInstanceChildValue(paramInstanceIndex,
                                         paramChildNames[i].c_str(),
                                         paramNames[j].c_str(),
                                         paramListArrayNames[k].c_str()),
                                         paramInstanceIndex,
                                         paramChildNames[i].c_str(), paramNames[j].c_str(), paramListArrayNames[k].c_str());
    }

    }}

    } // parameter instances

    }}

    //
    // Error checking : Abort if input paramListArray has abort on error flag set.
    //
	if((abortOnErrorFlag)&&(errorFlag))
	{
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("XML_ParameterListArrayUtilities Class Error \n");
	printf("copyParameterList(...) \n\n");
	printf("Copying a parameter list that does not exist in input XML_ParameterListArray. \n");
	printf("Input XML_ParameterListArray : ");
	printf("%s \n",paramListArray.getParameterListArrayName());
	printf("Offending parameterListName  : ");
	printf("%s \n",paramList);
	printf("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("\n");
	exit(0);
	}
}


//
// This routine assigns the values of all parameters of the parameters lists specified in the input parameter list
// to the corresponding parameter values of the identically named parameter lists in the specified output paramListArray.
//
//  An error message is generated if the parameter list, or parameter list parameter,
//  does not exist in the outputList.
//
void assignParameters(XML_ParameterListArray& paramListArray,  XML_ParameterListArray& outputList)
{
    vector<string> paramListArrayNames;
    paramListArray.getParameterListNames(paramListArrayNames);
    for(long k = 0; k < (long)paramListArrayNames.size(); k++)
    {
    assignParameters(paramListArray,paramListArrayNames[k].c_str(), outputList);
    }
}

//
// This routine assigns the values of all parameters specified in the input parameter list
// contained in the input paramListArray to the corresponding values of the identically named
// parameter list in the specified output paramListArray.
//
//  An error message is generated if the parameter list, or parameter list parameter,
//  does not exist in the outputList.
//
void assignParameters(XML_ParameterListArray& paramListArray,const char* paramList, XML_ParameterListArray& outputList)
{
    bool errorFlag = true;

    vector<string> paramListArrayNames;
    vector<string> paramNames;
    vector<string> paramChildNames;
    bool abortOnErrorFlag  = paramListArray.abortOnErrorFlag;
    XML_dataType XMLvalue;

    vector<string>::iterator it;

    long      paramInstanceIndex;
    paramListArray.getParameterListNames(paramListArrayNames);
    for(long k = 0; k < (long)paramListArrayNames.size(); k++)
    {

    if(paramListArrayNames[k].compare(paramList) == 0)
    {
    errorFlag = false;
    paramNames.clear();
    paramListArray.getParameterNames(paramListArrayNames[k].c_str(),paramNames);

    // Remove duplicate parameter names

    std::sort(paramNames.begin(),paramNames.end());
    it = std::unique(paramNames.begin(), paramNames.end());
    paramNames.resize(it - paramNames.begin());

    for(long j = 0; j < (long)paramNames.size(); j++)
    {
    for(paramInstanceIndex = 0;
    paramInstanceIndex < (long)paramListArray.getParameterInstanceCount(paramNames[j].c_str(),paramListArrayNames[k].c_str());
    paramInstanceIndex++)
    {

    paramChildNames.clear();
    paramListArray.getParameterChildNames(paramInstanceIndex, paramNames[j].c_str(),  paramListArrayNames[k].c_str(), paramChildNames);
    if(paramChildNames.size() == 0)
    {
    outputList.setParameter(paramListArray.getParameterValue(paramNames[j].c_str(),paramListArrayNames[k].c_str()),
    paramNames[j].c_str(), paramListArrayNames[k].c_str());
    }
    else // set parameter and child parameter values
    {
    paramListArray.clearAbortOnErrorFlag();
    XMLvalue = paramListArray.getParameterValue(paramNames[j].c_str(), paramListArrayNames[k].c_str());
    if(not XMLvalue.isNull())
    {outputList.setParameter(XMLvalue,paramNames[j].c_str(), paramListArrayNames[k].c_str());}
    else
    {paramListArray.clearError();}
    if(abortOnErrorFlag) paramListArray.setAbortOnErrorFlag();


    // Add child parameter values

    for(long i = 0; i < (int)paramChildNames.size(); i++)
    {
    outputList.setParameterInstanceChildValue(paramListArray.getParameterInstanceChildValue(paramInstanceIndex,
                                         paramChildNames[i].c_str(),
                                         paramNames[j].c_str(),
                                         paramListArrayNames[k].c_str()),
                                         paramInstanceIndex,
                                         paramChildNames[i].c_str(), paramNames[j].c_str(), paramListArrayNames[k].c_str());
    }

    }}

    } // parameter instances

    }}

    //
    // Error checking : Abort if input paramListArray has abort on error flag set.
    //
	if((abortOnErrorFlag)&&(errorFlag))
	{
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("XML_ParameterListArrayUtilities Class Error \n");
	printf("assignParameters(...) \n\n");
	printf("Attempting to assign parameters using a parameter list \n");
    printf("that does not exist in the input XML_ParameterListArray. \n\n");
	printf("Input XML_ParameterListArray : ");
	printf("%s \n",paramListArray.getParameterListArrayName());
	printf("Offending parameterListName  : ");
	printf("%s \n",paramList);
	printf("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("\n");
	exit(0);
	}
}

//
// Given a specified parameter child (named parameterChileName) whose
// value is a file, this routine replaces that file name with a full
// path file name based on the extraction of the base path using the
// realpath system call on the specified file. If the file does
// not exist then this routine is a no-op on that file name.
//
void insertFullPathNames(XML_ParameterListArray& paramList,
const char* parameterChildName, const char* parameterName,
const char*  parameterListName)
{
    string fileName;

    string basePath;
    string baseName;
    string fullFileName;

	for(long i = 0; i < paramList.getParameterInstanceCount(parameterName, parameterListName); i++)
	{
	if(paramList.isParameterInstanceChildValue(i,parameterChildName,parameterName, parameterListName))
	{
	fileName = (string)paramList.getParameterInstanceChildValue(i,parameterChildName,
	parameterName, parameterListName);
	basePath = getBasePath(fileName);
	baseName = getBaseName(fileName);

	if(basePath.empty())
	{
	return;
	}
    else
    {
	fullFileName = basePath;
	fullFileName.append("/");
	fullFileName.append(baseName);
	}
	paramList.setParameterInstanceChildValue(fullFileName,i,parameterChildName,parameterName, parameterListName);
	}
	}
}

//
// Given a specified parameter child (named parameterChileName) whose
// value is a file, this routine replaces that file name with a full
// path file name whose base path is specified.
//
void insertFullPathNames(const char* basePath, XML_ParameterListArray& paramList,
const char* parameterChildName, const char* parameterName,
const char*  parameterListName)
{
    string fileName;
    string baseName;
    string fullFileName;

	for(long i = 0; i < paramList.getParameterInstanceCount(parameterName, parameterListName); i++)
	{
	if(paramList.isParameterInstanceChildValue(i,parameterChildName,parameterName, parameterListName))
	{
	fileName = (string)paramList.getParameterInstanceChildValue(i,parameterChildName,parameterName, parameterListName);
	baseName = getBaseName(fileName);
	fullFileName = basePath;
	fullFileName.append("/");
	fullFileName.append(baseName);
	paramList.setParameterInstanceChildValue(fullFileName,i,parameterChildName,parameterName, parameterListName);
	}
	}
}

//
// If a file with name fileName does not exist then an empty string is returned.
//
string getBasePath(string fileName)
{
	const char *symlinkpath = fileName.c_str();
	char *actualpath;

	string actualPath;
	string   basePath;
    actualpath = realpath(symlinkpath, NULL);
    if (actualpath != NULL)
    {
    	actualPath.assign(actualpath);
    	free(actualpath);
    	basePath = actualPath.substr(0,actualPath.find_last_of("/\\"));
    	return basePath;
   }
   return basePath;
}

string getCWD()
{
	char *actualpath;
	string actualPath;
    actualpath = realpath("./", NULL);
    if (actualpath != NULL)
    {
    	actualPath.assign(actualpath);
    	free(actualpath);
   }
   return actualPath;
}
string getBaseName(string fileName)
{
	const char *symlinkpath = fileName.c_str();
	char *actualpath;

	string actualPath;
	string   baseName;
    actualpath = realpath(symlinkpath, NULL);
    if (actualpath != NULL)
    {
    	actualPath.assign(actualpath);
    	free(actualpath);
    	baseName = actualPath.substr(actualPath.find_last_of("/\\")+1);
    	return baseName;
   }
   else // Just strip away all directory modifiers
   {
	baseName = fileName.substr(fileName.find_last_of("/\\")+1);
   }
   return baseName;
}

bool fileExists(string fileName)
{
	if(getBasePath(fileName).empty()) return 0;
	return 1;
}
};



#endif /* _XML_ParameterListUtilities_ */
