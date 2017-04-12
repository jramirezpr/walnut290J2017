#include <iostream>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <vector>
using namespace std;

//
//******************************************************************************
//                 XML_ParameterListArray.h
//******************************************************************************
//
//********************************************************************************
//               Chris Anderson (C) UCLA  May 30, 2011
//********************************************************************************
//
// When initialized from an existing xml file, a non-explcitly typed value
// in the set of any case version of  ["true","yes","false","no"] are
// aliased as boolean values. To specify the actual string values of values
// in this set, one must specify the type explicitly as string.
//
// When initializing within a program, or setting an existing parameter,
// the automatic conversion from the set of any case version of
// ["true","yes",0,"false","no",1] also occur.
// To specify boolean one specifies the value with a bool
// variable.
//
// When evaluating parameter instance counts, if the specified parameterList
// doesn't exist, then the instance count returned is 0 and an error
// indicator is not set. This decision was made to avoid the need to specify
// empty parameter lists to avoid errors being thrown. Accessing an instance
// value of a parameter list that doesn't exist will induce the error flag to
// be set.
//
//


#ifndef _XML_ParameterListArray_
#define _XML_ParameterListArray_

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include "XML_dataType.h"
#include "tinyxml.h"

class XML_ParameterListArray
{

public:

	XML_ParameterListArray()
	{
		parameterArrayDocPtr = 0;
		abortOnErrorFlag     = true;
		initialize();
	}

	XML_ParameterListArray(const XML_ParameterListArray& P)
	{
		parameterArrayDocPtr = 0;
		abortOnErrorFlag     = true;
		initialize(P);
	}

	XML_ParameterListArray(const char* fileName)
	{
		parameterArrayDocPtr = 0;
		abortOnErrorFlag     = true;
		initialize(fileName);
	}

	~XML_ParameterListArray()
	{
	destroyData();
	}


////////////////////////////////////////////////////////////////////////
//
// For internal construction and manipulation of ParameterListArray
//
//
// To Do : Add parameter range checks. Write up documentation
//
////////////////////////////////////////////////////////////////////////

    void createParameterListArray(const char* listArrayName)
    {
    //
    // Set up header
    //

    string listFileName;
    listFileName.assign(listArrayName);
    listFileName.append(".xml");
    initialize();
    parameterArrayDocPtr           = new TiXmlDocument(listFileName.c_str());
    TiXmlDeclaration * declaration = new TiXmlDeclaration( "1.0", "", "" );
	TiXmlComment * comment         = new TiXmlComment();
	comment->SetValue(" XML_ParameterListArray ");


	parameterArrayDocPtr->LinkEndChild( declaration );
	parameterArrayDocPtr->LinkEndChild( comment );

    TiXmlElement * root = new TiXmlElement(listArrayName);
	parameterArrayDocPtr->LinkEndChild( root );

	abortOnErrorFlag     =  true;
    }

    const char* getParameterListArrayName()
    {
    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* root = docHandle.ToElement();
	return root->Value();
    }


    // Creates a duplicate

    void operator=(const XML_ParameterListArray& P)
    {
    initialize(P);
    }

    //
    // Check to make sure parameterListName doesn't exist more than once
    //

    void addParameterList(const char* parameterListName)
    {
    if(parameterListInstanceCount(parameterListName) != 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameterList(const char* parameterListName)\n\n");
	errorMessage.append("Adding a parameterList that already exists.\n");
	errorMessage.append("Multiple parameterLists in a single\nparameterListArray is not supported.\n\n");
	errorMessage.append("Offending parameterListName : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* root = docHandle.ToElement();
    TiXmlElement* parameterListElement = new TiXmlElement(parameterListName);
    root->LinkEndChild( parameterListElement );
    }


    void addParameter(const char* parameterName, const char* parameterListName)
    {

    if(isParameterList(parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameter(const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("ParameterList specified parameterListName does not exist.\n\n");
	errorMessage.append("Offending parameterListName : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterListElement    = docHandle.FirstChild(parameterListName).ToElement();
	TiXmlElement* parameterElement        = new TiXmlElement(parameterName);
    parameterListElement->LinkEndChild( parameterElement );
    }


    void addParameterNoTypeSpec(const char* value, const char* parameterName, const char* parameterListName)
    {

    if(isParameterList(parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameter(XML_dataType value, const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("ParameterList specified parameterListName does not exist.\n\n");
	errorMessage.append("Offending parameterListName : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterListElement    = docHandle.FirstChild(parameterListName).ToElement();
	TiXmlElement* parameterElement        = new TiXmlElement(parameterName);
    parameterListElement->LinkEndChild( parameterElement );

    parameterElement->SetAttribute("value",value);
    }



    void addParameter(XML_dataType value, const char* parameterName, const char* parameterListName)
    {

    if(isParameterList(parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameter(XML_dataType value, const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("ParameterList specified parameterListName does not exist.\n\n");
	errorMessage.append("Offending parameterListName : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterListElement    = docHandle.FirstChild(parameterListName).ToElement();
	TiXmlElement* parameterElement        = new TiXmlElement(parameterName);
    parameterListElement->LinkEndChild( parameterElement );

    parameterElement->SetAttribute("value", value.toString().c_str());
    if(value.isString()){parameterElement->SetAttribute("type","string");}
    else
    {parameterElement->SetAttribute("type",getDataType(value.toString().c_str()));}
    }

    void addParameterChild(const char* parameterChildName, const char* parameterName, const char* parameterListName)
    {

    if(isParameterList(parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameterChild(const char* parameterChildName, const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("ParameterList specified parameterListName does not exist.\n\n");
	errorMessage.append("Offending parameterListName : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameterChild(const char* parameterChildName, const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter specified by parameterName does not exist.\n\n");
	errorMessage.append("ParameterName : ");
	errorMessage.append(parameterName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter             = docHandle.FirstChild(parameterListName).FirstChild(parameterName).ToElement();
	TiXmlElement* parameterChild     = new TiXmlElement(parameterChildName);
    parameter->LinkEndChild( parameterChild );
    }

    //
    // check to make sure parameterList and parameter exist
    //

    void addParameterChild(XML_dataType value, const char* parameterChildName,
    const char* parameterName, const char* parameterListName)
    {
    if(isParameterList(parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameterChild(XML_dataType value, const char* parameterChildName,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("ParameterList specified parameterListName does not exist.\n\n");
	errorMessage.append("Offending parameterListName : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameterChild(XML_dataType value, const char* parameterChildName,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter specified by parameterName does not exist.\n\n");
	errorMessage.append("Parameter : ");
	errorMessage.append(parameterName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}


    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter             = docHandle.FirstChild(parameterListName).FirstChild(parameterName).ToElement();
	TiXmlElement* parameterChild     = new TiXmlElement(parameterChildName);
    parameter->LinkEndChild( parameterChild );
    parameterChild->SetAttribute("value", value.toString().c_str());
    if(value.isString())
    {parameterChild->SetAttribute("type","string");}
    else
    {parameterChild->SetAttribute("type",getDataType(value.toString().c_str()));}
    }


    void addParameterChildNoTypeSpec(XML_dataType value, const char* parameterChildName,
    const char* parameterName, const char* parameterListName)
    {
    if(isParameterList(parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameterChild(XML_dataType value, const char* parameterChildName,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("ParameterList specified parameterListName does not exist.\n\n");
	errorMessage.append("Offending parameterListName : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void addParameterChild(XML_dataType value, const char* parameterChildName,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter specified by parameterName does not exist.\n\n");
	errorMessage.append("Parameter : ");
	errorMessage.append(parameterName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}


    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter             = docHandle.FirstChild(parameterListName).FirstChild(parameterName).ToElement();
	TiXmlElement* parameterChild     = new TiXmlElement(parameterChildName);
    parameter->LinkEndChild( parameterChild );
    parameterChild->SetAttribute("value", value.toString().c_str());
    }



    void addParameterInstanceChild(XML_dataType value, int instanceIndex, const char* parameterChildName,
    const char* parameterName, const char* parameterListName)
    {

    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void setParameterChildValue(XML_dataType value,const char* childParameter,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter or ParameterList specified do not exist.\n\n");
	errorMessage.append("Parameter     : ");
	errorMessage.append(parameterName);
    errorMessage.append("\nParameterList : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

	TiXmlHandle   docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter          = docHandle.FirstChild(parameterListName).ToElement();
    TiXmlElement* parameterChild     = new TiXmlElement(parameterChildName);
	TiXmlNode*    node;

	long count = 0;
	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName))
	{
		if(count == instanceIndex)
		{
        node->ToElement()->LinkEndChild( parameterChild );
        parameterChild->SetAttribute("value", value.toString().c_str());
		if(value.isString())
		{parameterChild->SetAttribute("type","string");}
		else
		{parameterChild->SetAttribute("type",getDataType(value.toString().c_str()));}
		}
		count++;
	}

    }


    void setParameterOrIgnore(XML_dataType value, const char* parameterName, const char* parameterListName)
    {
    	if(this->isParameter(parameterName, parameterListName))
    	{
    		this->setParameter(value,parameterName, parameterListName);
    		return;
    	}
    	return;
    }

    void setParameter(XML_dataType value, const char* parameterName, const char* parameterListName)
    {

    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("setParameter(XML_dataType value, const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter or ParameterList specified do not exist.\n\n");
	errorMessage.append("Parameter     : ");
	errorMessage.append(parameterName);
    errorMessage.append("\nParameterList : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterElement    = docHandle.FirstChild(parameterListName).FirstChild(parameterName).ToElement();


    // Force a data conversion if necessary

    const char* inputDataType = getDataType(value.toString().c_str());
    const char* dataType      = parameterElement->Attribute("type");

	if(dataType == 0)
	{
	dataType = getDataType(parameterElement->Attribute("value"));
	parameterElement->SetAttribute("type",dataType);
    }

    std::string stringTemp;
    if((dataType != 0)&&(strcmp(dataType,inputDataType) != 0))
    {
    if(strcmp(dataType,"string") == 0) parameterElement->SetAttribute("value", XML_dataType(string(value)).toString().c_str());
    if(strcmp(dataType,"bool") == 0)   parameterElement->SetAttribute("value", XML_dataType(bool(value)).toString().c_str());
    if(strcmp(dataType,"long") == 0)   parameterElement->SetAttribute("value", XML_dataType(long(value)).toString().c_str());
    if(strcmp(dataType,"int") == 0)    parameterElement->SetAttribute("value", XML_dataType(int(value)).toString().c_str());
    if(strcmp(dataType,"float") == 0)  parameterElement->SetAttribute("value", XML_dataType(float(value)).toString().c_str());
    if(strcmp(dataType,"double") == 0) parameterElement->SetAttribute("value", XML_dataType(double(value)).toString().c_str());
    if((strcmp(dataType,"bool") == 0)&&(strcmp(inputDataType,"string")== 0))
    {
    stringTemp = XML_dataType(string(value)).toString();
    std::transform(stringTemp.begin(), stringTemp.end(),stringTemp.begin(), ::toupper);
	if(stringTemp.compare("FALSE") == 0)  parameterElement->SetAttribute("value","false");
    if(stringTemp.compare("TRUE")  == 0)  parameterElement->SetAttribute("value","true");
    if(stringTemp.compare("NO")    == 0)  parameterElement->SetAttribute("value","false");
	if(stringTemp.compare("YES")   == 0)  parameterElement->SetAttribute("value","true");
    }
    }
    else
    {
    parameterElement->SetAttribute("value", value.toString().c_str());
    if(dataType == 0)
    {
    if(value.isString())
	{parameterElement->SetAttribute("type","string");}
	else
	{parameterElement->SetAttribute("type",getDataType(value.toString().c_str()));}
	}
    }
    }


    void setParameterType(const char* typeValue, const char* parameterName, const char* parameterListName)
    {

    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("setParameter(XML_dataType value, const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter or ParameterList specified do not exist.\n\n");
	errorMessage.append("Parameter     : ");
	errorMessage.append(parameterName);
    errorMessage.append("\nParameterList : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterElement    = docHandle.FirstChild(parameterListName).FirstChild(parameterName).ToElement();

    {parameterElement->SetAttribute("type",typeValue);}
    }

    void setParameterChildValue(XML_dataType value,const char* childParameter,
    const char* parameterName, const char* parameterListName)
    {
    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void setParameterChildValue(XML_dataType value,const char* childParameter,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter or ParameterList specified do not exist.\n\n");
	errorMessage.append("Parameter     : ");
	errorMessage.append(parameterName);
    errorMessage.append("\nParameterList : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

    TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterChild  = docHandle.FirstChild(parameterListName).FirstChild(parameterName).FirstChild(childParameter).ToElement();

    const char* inputDataType;
    const char* dataType;

	if(parameterChild)
	{
	// Force a data conversion if necessary

    inputDataType = getDataType(value.toString().c_str());
    dataType      = parameterChild->Attribute("type");

    if(dataType == 0)
	{
	dataType = getDataType(parameterChild->Attribute("value"));
	parameterChild->SetAttribute("type",dataType);
    }

    std::string stringTemp;

    if((dataType != 0)&&(strcmp(dataType,inputDataType) != 0))
    {
    if(strcmp(dataType,"string") == 0) parameterChild->SetAttribute("value", XML_dataType(string(value)).toString().c_str());
    if(strcmp(dataType,"bool") == 0)   parameterChild->SetAttribute("value", XML_dataType(bool(value)).toString().c_str());
    if(strcmp(dataType,"long") == 0)   parameterChild->SetAttribute("value", XML_dataType(long(value)).toString().c_str());
    if(strcmp(dataType,"int") == 0)    parameterChild->SetAttribute("value", XML_dataType(int(value)).toString().c_str());
    if(strcmp(dataType,"float") == 0)  parameterChild->SetAttribute("value", XML_dataType(float(value)).toString().c_str());
    if(strcmp(dataType,"double") == 0) parameterChild->SetAttribute("value", XML_dataType(double(value)).toString().c_str());
    if((strcmp(dataType,"bool") == 0)&&((strcmp(inputDataType,"string"))==0))
    {
    stringTemp = XML_dataType(string(value)).toString();;
    std::transform(stringTemp.begin(),    stringTemp.end(),stringTemp.begin(), ::toupper);
	if(stringTemp.compare("FALSE") == 0)  parameterChild->SetAttribute("value","false");
    if(stringTemp.compare("TRUE")  == 0)  parameterChild->SetAttribute("value","true");
    if(stringTemp.compare("NO")    == 0)  parameterChild->SetAttribute("value","false");
	if(stringTemp.compare("YES")   == 0)  parameterChild->SetAttribute("value","true");
    }
    }
    else
    {
    parameterChild->SetAttribute("value", value.toString().c_str());

    if(dataType == 0)
    {
    if(value.isString())
	{parameterChild->SetAttribute("type","string");}
	else
	{parameterChild->SetAttribute("type",getDataType(value.toString().c_str()));}
	}

    }

    }
    else
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void setParameterChildValue(XML_dataType value,const char* childParameter,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter child value specifed  by childParameter does not exist.\n\n");
	errorMessage.append("Parameter Child     : ");
	errorMessage.append(childParameter);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}
    }


    void setParameterInstanceChildValue(XML_dataType value,int instanceIndex,const char* childParameter,
    const char* parameterName, const char* parameterListName)
    {
    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void setParameterChildValue(XML_dataType value,const char* childParameter,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter or ParameterList specified do not exist.\n\n");
	errorMessage.append("Parameter     : ");
	errorMessage.append(parameterName);
    errorMessage.append("\nParameterList : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).ToElement();

	TiXmlNode* node;

    errorFlag = true;

    const char* inputDataType;
    const char* dataType;

	long count = 0;
	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName) )
	{
		if(count == instanceIndex)
		{
			if(node->FirstChild(childParameter))
			{

			inputDataType = getDataType(value.toString().c_str());
            dataType      = (node->FirstChild(childParameter))->ToElement()->Attribute("type");
		    if((dataType != 0)&&(strcmp(dataType,inputDataType) != 0))
            {
                if(strcmp(dataType,"string") == 0)  (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(string(value)).toString().c_str());
                if(strcmp(dataType,"bool") == 0)    (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(bool(value)).toString().c_str());
		    	if(strcmp(dataType,"long") == 0)    (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(long(value)).toString().c_str());
		    	if(strcmp(dataType,"int") == 0)     (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(int(value)).toString().c_str());
		    	if(strcmp(dataType,"float") == 0)   (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(float(value)).toString().c_str());
		    	if(strcmp(dataType,"double") == 0)  (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(double(value)).toString().c_str());
            }
		    else
		    {
		    	(node->FirstChild(childParameter))->ToElement()->SetAttribute("value", value.toString().c_str()); // original

		    	if(dataType == 0)
		    	{
		    		if(value.isString())
		    		{(node->FirstChild(childParameter))->ToElement()->SetAttribute("type","string");}
		    		else
		    		{(node->FirstChild(childParameter))->ToElement()->SetAttribute("type",getDataType(value.toString().c_str()));}
		    	}

		    }

			errorFlag = false;
			}
			else
			{
			errorFlag = true;
			errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			errorMessage.append("XML_ParameterListArray Class Error \n");
			errorMessage.append("void setParameterInstanceChildValue(XML_dataType value,int instanceIndex,\n");
			errorMessage.append("const char* childParameter,const char* parameterName, const char* parameterListName)\n\n");
			errorMessage.append("Parameter child value specifed  by childParameter does not exist.\n\n");
			errorMessage.append("Parameter Child     : ");
			errorMessage.append(childParameter);
			errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			errorMessage.append("\n");
			}
			break;
		}
		count++;
	}
	if(abortOnErrorFlag){checkErrorAndAbort();}
	}


	void setParameterInstanceChildValue(XML_dataType value,const char* instanceName,const char* childParameter,
    const char* parameterName, const char* parameterListName)
    {
    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void setParameterChildValue(XML_dataType value,const char* childParameter,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter or ParameterList specified do not exist.\n\n");
	errorMessage.append("Parameter     : ");
	errorMessage.append(parameterName);
    errorMessage.append("\nParameterList : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}

	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).ToElement();
    TiXmlElement* parameterChild;
	TiXmlNode* node;
    string instanceString;
    string indexString(instanceName);
    errorFlag          = true;
    bool instanceError = true;

    const char* inputDataType;
    const char* dataType;

	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName) )
	{
		if(node->FirstChild("name"))
		{
		    parameterChild = node->FirstChild("name")->ToElement();
		    parameterChild->QueryStringAttribute("value",&instanceString);
		    if(instanceString.compare(indexString) == 0)
		    {
		    instanceError = false;
			if(node->FirstChild(childParameter))
			{

		    inputDataType = getDataType(value.toString().c_str());
            dataType      = (node->FirstChild(childParameter))->ToElement()->Attribute("type");
		    if((dataType != 0)&&(strcmp(dataType,inputDataType) != 0))
            {
                if(strcmp(dataType,"string") == 0)  (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(string(value)).toString().c_str());
                if(strcmp(dataType,"bool") == 0)    (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(bool(value)).toString().c_str());
		    	if(strcmp(dataType,"long") == 0)    (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(long(value)).toString().c_str());
		    	if(strcmp(dataType,"int") == 0)     (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(int(value)).toString().c_str());
		    	if(strcmp(dataType,"float") == 0)   (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(float(value)).toString().c_str());
		    	if(strcmp(dataType,"double") == 0)  (node->FirstChild(childParameter))->ToElement()->SetAttribute("value", XML_dataType(double(value)).toString().c_str());
            }
		    else
		    {
			(node->FirstChild(childParameter))->ToElement()->SetAttribute("value", value.toString().c_str()); // original

			if(dataType == 0)
		    {
		    	if(value.isString())
		    	{(node->FirstChild(childParameter))->ToElement()->SetAttribute("type","string");}
		    	else
		    	{(node->FirstChild(childParameter))->ToElement()->SetAttribute("type",getDataType(value.toString().c_str()));}
		   }


		   }

			errorFlag = false;
			}
			else
			{
			errorFlag = true;
			errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			errorMessage.append("XML_ParameterListArray Class Error \n");
			errorMessage.append("void setParameterInstanceChildValue(XML_dataType value,const char* instanceName,\n");
			errorMessage.append("const char* childParameter,const char* parameterName, const char* parameterListName)\n\n");
			errorMessage.append("Parameter child value specified  by childParameter does not exist.\n\n");
		    errorMessage.append("Parameter Instance Name     : ");
			errorMessage.append(indexString);
			errorMessage.append("\nParameter Child             : ");
			errorMessage.append(childParameter);

			errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			errorMessage.append("\n");
			}
			break;
		    }
	}
	}
    if(instanceError)
    {
    errorFlag = true;
    errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    errorMessage.append("XML_ParameterListArray Class Error \n");
    errorMessage.append("Parameter instance name specified does not exist.\n\n");
    errorMessage.append("Parameter List : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\n");
	errorMessage.append("Parameter      : ");
	errorMessage.append(parameterName);
	errorMessage.append("\n");
    errorMessage.append("ParameterChild : ");
	errorMessage.append(childParameter);
	errorMessage.append("\n");
    errorMessage.append("instance name  : ");
    errorMessage.append(indexString);
    errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    errorMessage.append("\n");
			}


	if(abortOnErrorFlag){checkErrorAndAbort();}
	}


	XML_dataType getParameterChildValueOrDefault(const char* childParameter,
    const char* parameterName, const char* parameterListName,XML_dataType defaultValue)
    {
	if(this->isParameterInstanceChildValue(0,childParameter,parameterName,parameterListName))
	{
		return this->getParameterChildValue(childParameter,parameterName,parameterListName);
	}
    return defaultValue;
    }

	XML_dataType getParameterInstanceChildValueOrDefault(int instanceIndex,const char* childParameter,
    const char* parameterName, const char* parameterListName,XML_dataType defaultValue)
    {
	if(this->isParameterInstanceChildValue(instanceIndex,childParameter,parameterName,parameterListName))
	{
		return this->getParameterInstanceChildValue(instanceIndex,childParameter,parameterName,parameterListName);
	}
    return defaultValue;
    }


	XML_dataType getParameterInstanceChildValue(const char* instanceName,const char* childParameter,
    const char* parameterName, const char* parameterListName)
    {
    XML_dataType returnValue;

    if(not isParameterList(parameterListName))
	{
		errorFlag = true;
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		errorMessage.append("XML_ParameterListArray Class Error \n");
		errorMessage.append("Parameter List not specified \n");
		errorMessage.append("Parameter List requested : ");
		errorMessage.append(parameterListName);
		errorMessage.append("\n");
		errorMessage.append("\n");
		if(not abortOnErrorFlag) return returnValue;
	}

	if(abortOnErrorFlag){checkErrorAndAbort();}

    if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("void setParameterChildValue(XML_dataType value,const char* childParameter,\n");
    errorMessage.append("const char* parameterName, const char* parameterListName)\n\n");
	errorMessage.append("Parameter or ParameterList specified do not exist.\n\n");
	errorMessage.append("Parameter     : ");
	errorMessage.append(parameterName);
    errorMessage.append("\nParameterList : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
	if(not abortOnErrorFlag) return returnValue;
    }

    if(abortOnErrorFlag){checkErrorAndAbort();}

	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).ToElement();
    TiXmlElement* parameterChild;
	TiXmlNode* node;
    string instanceString;
    string indexString(instanceName);
    errorFlag          = true;
    bool instanceError = true;

	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName) )
	{
		if(node->FirstChild("name"))
		{
		    parameterChild = node->FirstChild("name")->ToElement();
		    parameterChild->QueryStringAttribute("value",&instanceString);
		    if(instanceString.compare(indexString) == 0)
		    {
		    instanceError = false;
		    if(node->FirstChild(childParameter))
		    {
		    return getParameterValue(node->FirstChild(childParameter)->ToElement(),parameterName,  parameterListName);
		    }
		    else
		    {
		    errorFlag = true;
		    }
		    }
	     }
	}


    errorFlag = true;
    if(instanceError)
    {
    errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    errorMessage.append("XML_ParameterListArray Class Error \n");
    errorMessage.append("Parameter instance name specified does not exist.\n\n");
    errorMessage.append("Parameter List : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\n");
	errorMessage.append("Parameter      : ");
	errorMessage.append(parameterName);
	errorMessage.append("\n");
    errorMessage.append("ParameterChild : ");
	errorMessage.append(childParameter);
	errorMessage.append("\n");
    errorMessage.append("instance name  : ");
    errorMessage.append(indexString);
    errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    errorMessage.append("\n");
	if(abortOnErrorFlag){checkErrorAndAbort();}
	else return returnValue;
	}

    errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("Parameter Instance Child Parameter Not found.\n");
	errorMessage.append("Parameter List : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\n");
	errorMessage.append("Parameter      : ");
	errorMessage.append(parameterName);
	errorMessage.append("\n");
    errorMessage.append("ParameterChild : ");
	errorMessage.append(childParameter);
	errorMessage.append("\n");
	if(abortOnErrorFlag){checkErrorAndAbort();}

	return returnValue;
	}


	void getParameterListNames(vector < string >& paramListNames) const
	{
	paramListNames.clear();
	if(parameterArrayDocPtr == 0) return;
	TiXmlHandle   docHandle(parameterArrayDocPtr->RootElement());
	TiXmlNode* parameterList = docHandle.Element();
	TiXmlNode* child = 0;
    while((child = parameterList->IterateChildren( child )))
    {
	if(child->Type() != TiXmlNode::TINYXML_COMMENT)
	{
    paramListNames.push_back(child->Value());
	}
    }
	}


	void getParameterChildNames(long instanceIndex, const char* parameterName, const char* parameterListName, vector < string >& paramChildNames) const
	{
	paramChildNames.clear();
	if(not isParameter(parameterName, parameterListName)) return;

	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterList = docHandle.FirstChild(parameterListName).ToElement();

	TiXmlNode* node = 0;
    TiXmlNode* child = 0;

	long count = 0;
	for( node = parameterList->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName) )
	{
	    if(count == instanceIndex)
	    {
	    	child = 0;
	    	while((child = node->IterateChildren( child )))
	    	{
	    		if(child->Type() != TiXmlNode::TINYXML_COMMENT)
	    		{
                paramChildNames.push_back(child->Value());
	    		}
	    	}
        }
        count++;
    }
	}


    void getParameterNames(const char* parameterListName, vector < string >& paramNames) const
	{
	paramNames.clear();
	if(not isParameterList(parameterListName)) return;

	paramNames.clear();
	if(parameterArrayDocPtr == 0) return;
	TiXmlHandle   docHandle(parameterArrayDocPtr->RootElement());
	TiXmlNode* parameterList = docHandle.FirstChild(parameterListName).Element();
	TiXmlNode* child = 0;
    while((child = parameterList->IterateChildren( child )))
    {
    if(child->Type() != TiXmlNode::TINYXML_COMMENT)
	{
    paramNames.push_back(child->Value());
	}
    }
	}

    long parameterListCount() const
	{
	if(parameterArrayDocPtr == 0) return 0;
	TiXmlHandle   docHandle(parameterArrayDocPtr->RootElement());
	TiXmlNode* parameterList = docHandle.Element();
	TiXmlNode* child = 0;
	long count = 0;
    while((child = parameterList->IterateChildren( child )))
    {
    if(child->Type() != TiXmlNode::TINYXML_COMMENT)
    {
    count++;
    }
    }
    return count;
	}


    long parameterListInstanceCount(const char* parameterListName) const
	{
	if(parameterArrayDocPtr == 0) return 0;
	TiXmlHandle   docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterList = docHandle.FirstChild(parameterListName).ToElement();

	if(parameterList == 0){return 0;}
	TiXmlElement* node;
	long count = 0;
	for( node = parameterList;
		 node;
		 node = parameterList->NextSiblingElement(parameterListName))
	{
		if(node->Type() != TiXmlNode::TINYXML_COMMENT)
	    {
		count++;
	    }
	}

	return count;
	}

    /*  These need to be fixed */
    /*
    long parameterCount(const char* parameterListName)
    {
    if(parameterArrayDocPtr == 0) return 0;
    TiXmlHandle   docHandle(parameterArrayDocPtr->RootElement());

	TiXmlElement* parameterList = docHandle.FirstChild(parameterListName).ToElement();
	if(parameterList == NULL) return 0;

    TiXmlNode* parameter = parameterList->FirstChild();
    if(parameter == NULL) return 0;

    long count = 0;
	for(parameter; parameter; parameter = parameter->NextSiblingElement())
	{
	count++;
	}
	return count;
    }

    vector<string> getParameterNames(const char* parameterListName)
    {
    vector<string> paramNames;
    if(parameterArrayDocPtr == 0) return paramNames;

    TiXmlHandle   docHandle(parameterArrayDocPtr->RootElement());

	TiXmlElement* parameterList = docHandle.FirstChild(parameterListName).ToElement();
	if(parameterList == NULL) return paramNames;

    TiXmlNode* parameter = parameterList->FirstChild();
    if(parameter == NULL) return paramNames;

	for(parameter; parameter; parameter = parameter->NextSiblingElement())
	{
	paramNames.push_back(parameter->ToElement()->Value());
	}

	return paramNames;
    }
    */

////////////////////////////////////////////////////////////////////////
//
// Internal construction and manipulation of ParameterListArray
//
////////////////////////////////////////////////////////////////////////


	void initialize()
	{
	destroyData();
	}

	void initialize(const XML_ParameterListArray& P)
	{
	destroyData();
	parameterArrayDocPtr = new TiXmlDocument(*(P.parameterArrayDocPtr));
	abortOnErrorFlag     = P.abortOnErrorFlag;
	}

	void initialize(const char* fileName)
	{
		destroyData();
		FILE* dataFile;
		if((dataFile = fopen(fileName, "r" )) == NULL)
		{
		errorFlag = true;
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		errorMessage.append("XML_ParameterListArray Class Error \n");
		errorMessage.append("Error loading XML data from file : \n");
		errorMessage.append(fileName);
		errorMessage.append("\n");
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		}


		if(not errorFlag)
		{
		fclose(dataFile);
		parameterArrayDocPtr = new TiXmlDocument(fileName);
		if(parameterArrayDocPtr->LoadFile() == false)
		{
		errorFlag = true;
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		errorMessage.append("XML_ParameterListArray Class Error \n");
		errorMessage.append("Error loading XML data from file : \n");
		errorMessage.append(fileName);
		errorMessage.append("\n");
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		}}

		if(abortOnErrorFlag){checkErrorAndAbort();}

	}

	friend void operator >>(istream& in_stream, XML_ParameterListArray& P)
	{
    //
	//  Output using pretty printing
	//
		P.initialize();
		P.parameterArrayDocPtr = new TiXmlDocument();
		in_stream >> *P.parameterArrayDocPtr;
	}

	void destroyData()
	{
		if(parameterArrayDocPtr != 0){delete parameterArrayDocPtr;}
		parameterArrayDocPtr = 0;
		errorFlag            = false;
		errorMessage.clear();
	}

	friend ostream&  operator <<(ostream& out_stream, const XML_ParameterListArray& P)
	{
    //
	//  Output using pretty printing
	//
		TiXmlPrinter printer;
		P.parameterArrayDocPtr->Accept( &printer );
		out_stream << printer.CStr();
		return out_stream;
	}

	void checkErrorAndAbort()
	{
    if(errorFlag)
    {
    	printf("\n%s\n",getErrorMessage());
    	exit(1);
    }
	}

	void        setAbortOnErrorFlag(){abortOnErrorFlag   = true;}
	void        clearAbortOnErrorFlag(){abortOnErrorFlag = false;}
	bool        getErrorFlag()   {return errorFlag;}
	const char* getErrorMessage(){return errorMessage.c_str();}


	void clearError()
	{
	errorFlag = false;
	errorMessage.clear();
	}

	int isParameterList(const char* parameterListName) const
	{
	if(parameterArrayDocPtr == 0) return 0;
	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameterList = docHandle.FirstChild(parameterListName).ToElement();
	if(parameterList) {return 1;}
	return 0;
	}

	int isParameter(const char* parameterName, const char* parameterListName) const
	{
	if(isParameterList(parameterListName) == 0)
    {return 0;}

	if(parameterArrayDocPtr == 0) return 0;
	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).FirstChild(parameterName).ToElement();
	if(parameter) {return 1;}
	return 0;
	}

	long getParameterInstanceCount(const char* parameterName, const char* parameterListName)
	{
    if(isParameterList(parameterListName) == 0) {return 0;}

	if(parameterArrayDocPtr == 0) return 0;
	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).ToElement();

	TiXmlNode* node;

	long count = 0;
	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName))
	{
		count++;
	}
	return count;
	}

	//
	// Indexing of the instance starts at 0 to facilitate extracting data into C type arrays
	//
	XML_dataType getParameterInstanceValue(int instanceIndex,const char* parameterName, const char* parameterListName)
	{
    XML_dataType returnValue;
	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).ToElement();

	TiXmlNode* node;

	long count = 0;
	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName) )
	{
		if(count == instanceIndex)
		{
			if(node->ToElement())
			{
			returnValue = getParameterValue(node->ToElement(),parameterName,  parameterListName);
			}
			else
			{
			errorFlag   = true;
		    errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	     	errorMessage.append("XML_ParameterListArray Class Error \n");
	     	errorMessage.append("Parameter Instance Not found.\n");
	     	errorMessage.append("Parameter List : ");
	     	errorMessage.append(parameterListName);
	     	errorMessage.append("\n");
	     	errorMessage.append("Parameter      : ");
	     	errorMessage.append(parameterName);
	     	errorMessage.append("\n");
			}
		}
		count++;
	}

	if(returnValue.isNull())
	{
	errorFlag = true;
	}

	if(errorFlag)
	{
	errorMessage.append("Parameter Instance : ");
	stringStream.str("");
	stringStream << instanceIndex;
	errorMessage.append(stringStream.str());
	errorMessage.append("\n");
	}

	if(abortOnErrorFlag){checkErrorAndAbort();}
	return returnValue;
	}


	//
	// Indexing of the instance starts at 0 to facilitate extracting data into C type arrays
    //

	XML_dataType getParameterInstanceChildValue(int instanceIndex,const char* childParameter,
			const char* parameterName, const char* parameterListName)
	{
	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).ToElement();

	TiXmlNode* node;
	XML_dataType returnValue;

	long count = 0;
	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName) )
	{
	    if(count == instanceIndex)
	    {
	    	if(node->FirstChild(childParameter))
	    	{
	    	returnValue =  getParameterValue(node->FirstChild(childParameter)->ToElement(),parameterName,  parameterListName);
	    	}
	    	else
	    	{
	    	errorFlag   = true;
	     	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	     	errorMessage.append("XML_ParameterListArray Class Error \n");
	     	errorMessage.append("Parameter Instance Child Value Not found.\n");
	     	errorMessage.append("Parameter List : ");
	     	errorMessage.append(parameterListName);
	     	errorMessage.append("\n");
	     	errorMessage.append("Parameter      : ");
	     	errorMessage.append(parameterName);
	     	errorMessage.append("\n");
	    	}
	    }
	    count++;
	}


	if((returnValue.isNull())&&(not errorFlag))
	{
	errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("Parameter Instance Not found.\n");
	errorMessage.append("Parameter List : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\n");
	errorMessage.append("Parameter      : ");
	errorMessage.append(parameterName);
	errorMessage.append("\n");
	}

	if(errorFlag)
	{
	    	errorMessage.append("Parameter Instance : ");
	    	stringStream.str("");
	    	stringStream << instanceIndex;
	    	errorMessage.append(stringStream.str());
	    	errorMessage.append("\n");
	        errorMessage.append("Parameter Child  : ");
	    	errorMessage.append(childParameter);
	    	errorMessage.append("\n");
	}
	if(abortOnErrorFlag){checkErrorAndAbort();}
	return returnValue;
	}


    bool isParameterInstanceChildValue(int instanceIndex,const char* childParameter,
    const char* parameterName, const char* parameterListName) const
	{
	if(parameterArrayDocPtr == 0) return false;

	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).ToElement();

	TiXmlNode* node;

	long count = 0;
	for( node = parameter->FirstChild(parameterName);
		 node;
		 node = node->NextSibling(parameterName) )
	{
	    if(count == instanceIndex)
	    {
	    	if(node->FirstChild(childParameter))
	    	{
	    	return true;
	    	}
	    	else
	    	{
	    	return false;
	    	}
	    }
	    count++;
	}
    return false;
	}

    XML_dataType getParameterChildValue(const char* childParameter, const char* parameterName, const char* parameterListName)
	{
    XML_dataType returnValue;
	if(not isParameterList(parameterListName))
	{
		errorFlag = true;
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		errorMessage.append("XML_ParameterListArray Class Error \n");
		errorMessage.append("Parameter List not specified \n");
		errorMessage.append("Parameter List requested : ");
		errorMessage.append(parameterListName);
		errorMessage.append("\n");
		errorMessage.append("\n");
		if(not abortOnErrorFlag) return returnValue;
	}
	if(abortOnErrorFlag){checkErrorAndAbort();}

	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).FirstChild(parameterName).FirstChild(childParameter).ToElement();

	returnValue = getParameterValue(parameter, parameterName,  parameterListName);

	if(errorFlag)
	{
	errorMessage.append("Parameter Child  : ");
	errorMessage.append(childParameter);
	errorMessage.append("\n");
	if(not abortOnErrorFlag) {returnValue.initialize(); return returnValue;}
	}

	if(abortOnErrorFlag){checkErrorAndAbort();}
	return returnValue;
	}


   //
   // If the specified parameter exists, then this routine returns the value, otherwise it returns the
   // defaultValue specified.
   //
   XML_dataType getParameterValueOrDefault(const char* parameterName, const char* parameterListName, XML_dataType defaultValue )
   {

   if(this->isParameter(parameterName,parameterListName))
   {
    return this->getParameterValue(parameterName,parameterListName);
   }
   return defaultValue;
   }
	//
	// If one wants to extract a boolean value as a string one must append the toString() member function, e.g.
	//
	// string boolAsString = paramList.getParameterValue(...).toString();
	//


	XML_dataType getParameterValue(const char* parameterName, const char* parameterListName)
	{
	XML_dataType returnValue;

	if(not isParameterList(parameterListName))
	{
		errorFlag = true;
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		errorMessage.append("XML_ParameterListArray Class Error \n");
		errorMessage.append("Parameter List not specified \n");
		errorMessage.append("Parameter List requested : ");
		errorMessage.append(parameterListName);
		errorMessage.append("\n");
		errorMessage.append("\n");
		if(not abortOnErrorFlag) return returnValue;
	}
	if(abortOnErrorFlag){checkErrorAndAbort();}

	if(isParameter(parameterName,parameterListName) == 0)
    {
    errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("Parameter specified by parameterName does not exist.\n\n");
	errorMessage.append("ParameterName : ");
	errorMessage.append(parameterName);
	errorMessage.append("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("\n");
	if(not abortOnErrorFlag) return returnValue;
    }
    if(abortOnErrorFlag){checkErrorAndAbort();}


	TiXmlHandle  docHandle(parameterArrayDocPtr->RootElement());
	TiXmlElement* parameter = docHandle.FirstChild(parameterListName).FirstChild(parameterName).ToElement();

	returnValue = getParameterValue(parameter, parameterName,  parameterListName);

	if(abortOnErrorFlag){checkErrorAndAbort();}
	return returnValue;
	}

	XML_dataType getParameterValue(const TiXmlElement* parameter,const char* parameterName, const char* parameterListName)
	{
	int           intParam;
	double     doubleParam;
	float       floatParam;
	bool         boolParam;
	string     stringParam;
	string      stringTemp;
	bool     explicitType;



	const char* dataType = 0;
	if(parameter)
	{
		//
		// Check for values attribute
		//
		if(parameter->Attribute("value") == 0)
		{
		errorFlag = true;
		errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		errorMessage.append("XML_ParameterListArray Class Error \n");
		errorMessage.append("Parameter Value not specified using value= \"...\"\n");
		errorMessage.append("Parameter List : ");
		errorMessage.append(parameterListName);
		errorMessage.append("\n");
		errorMessage.append("Parameter      : ");
		errorMessage.append(parameterName);
		errorMessage.append("\n");
		return XML_dataType();
		}
		//
		// Check for explicitly typed parameters
		//
		explicitType = true;
		dataType = parameter->Attribute("type");

		//
		// If not explicitly typed, then discover the type by interrogating the string representation
		//

		if(dataType == 0)
		{
		 stringParam = parameter->Attribute("value");
		 dataType = getDataType(stringParam.c_str());
		 explicitType = false;
		}

		if(strcmp(dataType,"bool") == 0)
		{
			if(parameter->QueryBoolAttribute("value",&boolParam)  == TIXML_SUCCESS) {return XML_dataType(boolParam);}
		}
		if(strcmp(dataType,"int") == 0)
		{
			if(parameter->QueryIntAttribute("value",&intParam)  == TIXML_SUCCESS) {return XML_dataType(intParam);}
		}
		if(strcmp(dataType,"long") == 0)
		{
			if(parameter->QueryIntAttribute("value",&intParam)  == TIXML_SUCCESS) {return XML_dataType((long)intParam);}
		}
		if(strcmp(dataType,"float") == 0)
		{
			if(parameter->QueryFloatAttribute("value",&floatParam)  == TIXML_SUCCESS) {return XML_dataType(floatParam);}
		}
		if(strcmp(dataType,"double") == 0)
		{
			if(parameter->QueryDoubleAttribute("value",&doubleParam)  == TIXML_SUCCESS) {return XML_dataType(doubleParam);}
		}
		//
		// For strings, if not explicitly typed, check to see if it's a boolean variable first
		// as designated by any upper or lower case collection of letters spelling true, false
		// or yes or no.
		//
		// I don't use tinyXML QueryBoolAttribute since it returns incorrectly if one specifies the
		// string None.

		if(strcmp(dataType,"string") == 0)
		{
			//if(parameter->QueryBoolAttribute("value",&boolParam)  == TIXML_SUCCESS)     {return XML_dataType(boolParam);}

            if(parameter->QueryStringAttribute("value",&stringParam)  == TIXML_SUCCESS)
            {
            if(explicitType){return XML_dataType(stringParam);}
			else // interrogate first for conversion to boolean based on false, true, yes or no,
			{
			stringTemp = stringParam;
			std::transform(stringTemp.begin(), stringTemp.end(),stringTemp.begin(), ::toupper);
			if(stringTemp.compare("FALSE") == 0) {boolParam = false; return XML_dataType(boolParam);}
			if(stringTemp.compare("TRUE")  == 0) {boolParam = true;  return XML_dataType(boolParam);}
			if(stringTemp.compare("NO") == 0)    {boolParam = false; return XML_dataType(boolParam);}
			if(stringTemp.compare("YES")  == 0)  {boolParam = true;  return XML_dataType(boolParam);}
			return XML_dataType(stringParam);
			}}
		}
	}

	errorFlag = true;
	errorMessage.append("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	errorMessage.append("XML_ParameterListArray Class Error \n");
	errorMessage.append("Parameter Not found.\n");
	errorMessage.append("Parameter List : ");
	errorMessage.append(parameterListName);
	errorMessage.append("\n");
	errorMessage.append("Parameter      : ");
	errorMessage.append(parameterName);
	errorMessage.append("\n");

	return XML_dataType();
	}



//
// Interrogates string and returns one of
// bool, double, long, string
//
//
// Case independent versions of {true,yes} and {false,no}
// are returned as type bool
//
const char* getDataType(const char* sIn)
{
    //
    // remove white space before and after variables
    //
    string sInString(sIn);
    string sTrimmed(trim(sInString));
    const char* s   = sTrimmed.c_str();

    const char* boolType = "bool";

    std::string stringTemp(s);
	std::transform(stringTemp.begin(), stringTemp.end(),stringTemp.begin(), ::toupper);
	if(stringTemp.compare("FALSE") == 0)  return boolType;
	if(stringTemp.compare("TRUE")  == 0)  return boolType;
	if(stringTemp.compare("NO")    == 0)  return boolType;
	if(stringTemp.compare("YES")   == 0)  return boolType;

	int sLength   = strlen(s);
	int firstChar = (int)s[0];
	int scndChar;
	int thrdChar;

	const char* doubleType = "double";
	const char* longType   = "long";
	const char* stringType = "string";

	//
	// check starting digit first 0-9 (48-57)
	//
	if((firstChar >= 48)&&(firstChar <= 57))
	{
	if(strchr(s,'.'))   {return doubleType;}
	else                {return longType;}
	}

	// check starting with (.)

	if(firstChar == 46)
	{
    if(sLength == 1)                        {return stringType;}
    scndChar = (int)s[1];
    if((scndChar >= 48)&&(scndChar <= 57))  {return doubleType;}
    else                                    {return stringType;}
	}

    //
	// check starting with a (+) (-)
	//

    if((firstChar == 43)||(firstChar == 45))
	{
	if(sLength == 1){return stringType;}

    scndChar = (int)s[1];
	if((scndChar >= 48)&&(scndChar <= 57))
	{
	if(strchr(s,'.'))   {return doubleType;}
	else                {return longType;}
	}

	if(scndChar == 46)
	{
    if(sLength == 2)                        {return stringType;}
    thrdChar = (int)s[2];
    if((thrdChar >= 48)&&(thrdChar <= 57))  {return doubleType;}
    else                                    {return stringType;}
	}
	}

	//
	// Check for string type
	//

	if((firstChar == 39)||(firstChar == 34)) {return stringType;}
	return stringType;
}


// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}


	public:

	TiXmlDocument* parameterArrayDocPtr;
	string         errorMessage;
	bool              errorFlag;
	ostringstream  stringStream;
	bool       abortOnErrorFlag;

};

#endif

