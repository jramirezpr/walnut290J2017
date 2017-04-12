//
// XML_dataType.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
#ifndef _XML_dataType_
#define _XML_dataType_
#include <string>
#include <iostream>
#include <cstdlib>
using namespace std;

class  XML_dataType
{

public :

    bool    b;
    int     i;
    long    l;
    float   f;
    double  d;
    string  s;
    void*   v;
    int  type_name;

public :

//
//  Type Enumerations
//

    enum {TYPE_NULL, TYPE_BOOL, TYPE_INT ,TYPE_LONG ,TYPE_FLOAT ,TYPE_DOUBLE ,TYPE_STRING ,TYPE_BOOL_PTR, TYPE_INT_PTR ,TYPE_LONG_PTR ,TYPE_FLOAT_PTR ,TYPE_DOUBLE_PTR ,TYPE_STRING_PTR  };
//
//  Constructors
//
    XML_dataType();
    XML_dataType(const XML_dataType& A);

    XML_dataType(bool A );
    XML_dataType(int A );
    XML_dataType(long A );
    XML_dataType(float A );
    XML_dataType(double A );
    XML_dataType(const string& A );

    XML_dataType(const char* A);

    //XML_dataType(char* A);
    XML_dataType(int* A );
    XML_dataType(bool* A );
    XML_dataType(long* A );
    XML_dataType(float* A );
    XML_dataType(double* A );
    XML_dataType(string* A );
//
//  Destructor
//
    ~XML_dataType();

//
//  Initializers
//
    void initialize();
//
//  Assignment
//
    XML_dataType& operator = ( const XML_dataType& A);
//
//  Output
//
    friend ostream& operator <<(ostream& out_stream, const XML_dataType& A);
//
//  conversions
//
	//operator const char*();
	operator int();
	operator bool();
	operator long();
	operator float();
	operator double();
	operator string();
	friend istream&  operator >>(istream& in_stream, XML_dataType& A);
//
//  assignments
//
	void  operator =(bool A);
    void  operator =(int A);
    void  operator =(long A);
    void  operator =(float A);
    void  operator =(double A);
    void  operator =(const string& A);
    void  operator =(const char* A);

    //void  operator =(const char* A);
//
//  Utility for output
//
    string toString();
//
//  Typing mutatators
//

    bool  isString();
    void  setType(int XML_dataTypeType);
    int   getType();

    int isNull(){if(getType() == TYPE_NULL){return 1;} return 0;}
//
//  Error Handling
//
    static void  illegalConversion();
    static void  illegalConversion(int typeA, int typeB);
    static void  nullOperand();
    static void  illegalAssignment();
    static void  illegalAssignment(int typeA, int typeB);

};
#endif

