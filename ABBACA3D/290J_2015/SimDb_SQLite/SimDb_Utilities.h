#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
using namespace std;

#include "sqlite3.h"

// Author: Chris Anderson
// (C) UCLA 2012

// Class SimDb_Utilities: Miscellaneous member functions for working
// with the simulation database.
//
// Motivation:
//
// By deciding that all non-blob simulation data values will be stored as strings has the
// benefit of avoiding any problems with unwanted implicit conversions taking place when
// one stores or retrieves information from a simulation database. The downside to using
// all string values is that type information has to be encoded with each value
// and data conversion has to be carried out externally.
//
// The majority of the utility functions in this class have been
// created to remove some of the burden of data transformation.
//

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
#ifndef _SimDb_Utilities_
#define _SimDb_Utilities_

class SimDb_Utilities
{
	public:

	SimDb_Utilities(){};

//
// This routine removes the all but one or two trailing zeros from a string containing
// a floating point number in scientific format. All zeros aren't removed to
// improve the readability of the resulting string.
//
	void removeScientificFormatTrailingZeros(string inputString, string& outputString)
	{
	vector < string > stringParts;
	string delimiters = "Ee.";
	createTokens(inputString,stringParts,delimiters);


    outputString.clear();
    outputString.append(stringParts[0]);
    outputString.append(".");

    if(stringParts.size() == 1)
    {
    outputString.append("0");
    return;
    }

    long itrim = stringParts[1].length();

    for(long i = stringParts[1].length()-1; i >= 1 ; i--)
    {
    if(stringParts[1].at(i) != '0') break;
    itrim --;
    }
    itrim += 1;
    if(itrim > (long)stringParts[1].length()) itrim = stringParts[1].length();

    outputString.append(stringParts[1].substr(0,itrim));
    if(stringParts.size() > 2)
    {
    outputString.append("e");
    outputString.append(stringParts[2]);
    }
    return;
	}

	string getFileNameExtension(const string& fileName)
	{
	string extension;
	vector < string > fileParts;
	createTokens(fileName,fileParts,".");
	if(fileParts.size() > 1)
	{
	extension = fileParts.back();
	}
	return extension;
	}

    //
    // This construction allows for single .'s embedded in the
    // input file name (multiple consecutive dots are not allowed).
    //
    string getFileNameBase(const string& fileName)
	{
	string fileBase;
	vector < string > fileParts;
	createTokens(fileName,fileParts,".");
	if(fileParts.size() <= 1) return fileName;

	fileBase.assign(fileParts[0]);
	for(long i = 1; i < (long)fileParts.size()-1; i++)
	{
	fileBase.append(".");
	fileBase.append(fileParts[i]);
	}
    return fileBase;
	}
    //
    // SQLite functions require UTF8 (= char) strings, so this utility is
    // used for constructing a char* string from the characters in the inputString
    //
    void createUTF8characterEncoding(const string& inputString, vector < char >& charArray)
    {
    charArray.clear();
    charArray.resize(inputString.size()+1);
    for(long i = 0; i < (long)inputString.size(); i++)
    {
    charArray[i] = (char)inputString.at(i);
    }
    charArray[inputString.size()] = '\0';
    }


    bool fileExists(const char* fileName)
    {
	struct stat fileInfo;
	int statReturn;
	statReturn = stat(fileName,&fileInfo);
	if(statReturn  == 0) {return true;}
	return false;
    }

    bool dirExists(const char* dirName)
    {
	struct stat fileInfo;
	int statReturn;
	statReturn = stat(dirName,&fileInfo);
	if(statReturn  == 0) {return true;}
	return false;
    }

    string getCurrentWorkingDirectory()
    {
    char cwdBuffer[4096];
    char* workingDirPtr = getcwd(cwdBuffer,4096);
    string workingDir;
    if(workingDirPtr != 0) workingDir.assign(workingDirPtr);
    return workingDir;
    }

    string trim(const string & inputString)
	{
    string whitespace = " \t";
    string returnString;
    size_t beginStr = inputString.find_first_not_of(whitespace);
    if (beginStr == std::string::npos) {return returnString;}

    size_t endStr = inputString.find_last_not_of(whitespace);
    size_t range  = endStr - beginStr + 1;
    returnString = inputString.substr(beginStr,range);
    return returnString;
	}
	
	string trimAndRemoveNewLines(const string & inputString)
	{
    string whitespace = " \t\n";
    string returnString;
    size_t beginStr = inputString.find_first_not_of(whitespace);
    if (beginStr == std::string::npos) {return returnString;}

    size_t endStr = inputString.find_last_not_of(whitespace);
    size_t range  = endStr - beginStr + 1;
    returnString = inputString.substr(beginStr,range);
    return returnString;
	}

    string getPID()
    {
    int PIDval       =  getpid();
    string PIDstring =  getSQdataString(PIDval);
    return PIDstring;
    }

    const string& getSQdataString(const string& V)
    {
    return V;
    }

    string getSQdataString(string& V)
    {
    string Vreturn(V);
    return Vreturn;
    }

    string getSQdataString(const int& val)
    {
    string SQdata;
    ostringstream convert;
    convert << val; SQdata.assign(convert.str());
    return SQdata;
    }

    string getSQdataString(const long& val)
    {
    string SQdata;
    ostringstream convert;
    convert << val; SQdata.assign(convert.str());
    return SQdata;
    }

    string getSQdataString(const float& val)
    {
    string SQdata;
    ostringstream convert;
    string inputString;
    string cleanString;
    convert << scientific << setprecision(15) << val;
    inputString.assign(convert.str());
    removeScientificFormatTrailingZeros(inputString, cleanString);
    SQdata.assign(cleanString);
    return SQdata;
    }


    string getSQdataString(const double& val)
    {
    string SQdata;
    ostringstream convert;
    string inputString;
    string cleanString;
    convert << scientific << setprecision(15) << val;
    inputString.assign(convert.str());
    removeScientificFormatTrailingZeros(inputString, cleanString);
    SQdata.assign(cleanString);
    return SQdata;
    }

    int isInArray(const vector <string>& sArray, const string& s)
    {
	if (std::find(sArray.begin(), sArray.end(), s) != sArray.end())
    {return 1;}
    return 0;
    }

    bool allUnique(const vector <string>& inputArray)
    {
    vector <string > sArray = inputArray;
    vector<string>::iterator it;
    std::sort( sArray.begin(),sArray.end());
    it = std::unique(sArray.begin(), sArray.end());
    if(it ==  sArray.end()){return true;}
    return false;
    }

    void removeDuplicates(vector <string>& sArray)
    {
    vector<string>::iterator it;
    std::sort( sArray.begin(),sArray.end());
    it = std::unique (sArray.begin(), sArray.end());
    sArray.resize(it - sArray.begin());
    }


    void createTokens(const string& str,vector<string>& tokens, const string& delimiters)
    {
    tokens.clear();
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
   }

   // This routine scans the directory specified for file names with extension ext and
   // add them to fileNames if they don't already exist.
   //
   // This is a utility routine for capturing files that are created during a run,
   //
   void addAdditionalDirFiles(const string& directory, vector <string >& fileNames,const string& ext)
   {
	   DIR *dir;
	   struct dirent *dirEntry;
	   string fileName;
	   dir = opendir(directory.c_str());
	   if (dir != NULL)
	   {
		   while ((dirEntry  = readdir (dir)) != NULL)
		   {
			   fileName.assign(dirEntry->d_name);
			   if(ext.compare(getFileNameExtension(fileName)) == 0)
			   {
				   if(isInArray(fileNames, fileName) == 0)
				   {
					   fileNames.push_back(fileName);
				   }
				 }
           	   }
       closedir (dir);
       }
       return;
   	 }
};

#endif
