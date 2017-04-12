#include "XML_dataType.h"
#include <sstream>
//
//******************************************************************************
//                      XML_DataType.cpp
//******************************************************************************
//
//********************************************************************************
//               Chris Anderson (C) UCLA 2011
//********************************************************************************
//
//********************************************************************************
//                     Class XML_dataType
//********************************************************************************
//
/*
XML_dataType is a class which provides a minimal wrapping for a subset
of the basic data types available in C++. The data types supported are
bool, int, long, float, double and strings. The class is a
container class which holds the value of the data and it's type.
The constructors and conversions associated with this class have been
implemented so that functions which return a XML_dataType object can
"communicate" with variables of the standard data types. This class is
purely a mechanism for implementing functions with identical signature
but different return types. Algebraic operations (other than assignment)
are not implemented for the objects of this class.
*/
XML_dataType::XML_dataType(): s()
{
	b = false;
	i = 0;
	l = 0;
	f = 0.0;
	d = 0.0;
	v = 0;
	type_name = TYPE_NULL;
}

XML_dataType::XML_dataType( const XML_dataType& A)
{
	b = A.b;
	i = A.i;
	l = A.l;
	f = A.f;
	d = A.d;
	s = A.s;
	v = A.v;
	type_name = A.type_name;

}

XML_dataType::XML_dataType(bool A ): s()
{
	b = A;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = 0;
    type_name = TYPE_BOOL;
}



XML_dataType::XML_dataType(int A ): s()
{
	b = false;
    i = A;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = 0;
    type_name = TYPE_INT;
}

XML_dataType::XML_dataType(long A ): s()
{
	b = false;
    i = 0;
    l = A;
    f = 0.0;
    d = 0.0;
    v = 0;
    type_name = TYPE_LONG;
}

XML_dataType::XML_dataType(float A ): s()
{
	b = false;
    i = 0;
    l = 0;
    f = A;
    d = 0.0;
    v = 0;
    type_name = TYPE_FLOAT;
}

XML_dataType::XML_dataType(double A ): s()
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = A;
    v = 0;
    type_name = TYPE_DOUBLE;
}

XML_dataType::XML_dataType(const string& A ): s(A)
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = 0;
    type_name = TYPE_STRING;
}


XML_dataType::XML_dataType(const char* A): s(A)
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = 0;
    type_name = TYPE_STRING;
}


XML_dataType::XML_dataType(bool* A ): s()
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = A;
    type_name = TYPE_BOOL_PTR;
}

XML_dataType::XML_dataType(int* A ): s()
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = A;
    type_name = TYPE_INT_PTR;
}

XML_dataType::XML_dataType(long* A ): s()
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = A;
    type_name = TYPE_LONG_PTR;
}

XML_dataType::XML_dataType(float* A ): s()
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = A;
    type_name = TYPE_FLOAT_PTR;
}

XML_dataType::XML_dataType(double* A ): s()
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = A;
    type_name = TYPE_DOUBLE_PTR;
}

XML_dataType::XML_dataType(string* A): s()
{
	b = false;
    i = 0;
    l = 0;
    f = 0.0;
    d = 0.0;
    v = (void*)A;
    type_name = TYPE_STRING_PTR;
}
//
//********************************************************************************
//                    DESTRUCTOR
//********************************************************************************
//
XML_dataType::~XML_dataType()
{
}

//
//********************************************************************************
//                    Initialize
//********************************************************************************
//
void XML_dataType::initialize()
{
	b = false;
	i = 0;
	l = 0;
	f = 0.0;
	d = 0.0;
	v = 0;
	s.clear();
	type_name = TYPE_NULL;
}

//
//********************************************************************************
//                    ASSIGNMENT
//********************************************************************************
//
XML_dataType&  XML_dataType::operator =( const XML_dataType& A)
{
//
//  Memberwise Assignment
//
	b = A.b;
    i = A.i;
    l = A.l;
    f = A.f;
    d = A.d;
    s = A.s;
    v = A.v;
	type_name = A.type_name;

    return *this;
}

//
//********************************************************************************
//                    OUTPUT
//********************************************************************************
//
ostream&  operator <<(ostream& out_stream, const XML_dataType& A)
{
//
//  Default Output  : Memberwise Output
//
     switch(A.type_name)
     {
     case XML_dataType::TYPE_NULL         : out_stream << "Null Type"; break;
     case XML_dataType::TYPE_BOOL         : out_stream << A.b; break;
     case XML_dataType::TYPE_INT          : out_stream << A.i; break;
     case XML_dataType::TYPE_LONG         : out_stream << A.l; break;
     case XML_dataType::TYPE_FLOAT        : out_stream << A.f; break;
     case XML_dataType::TYPE_DOUBLE       : out_stream << A.d; break;
     case XML_dataType::TYPE_STRING       : out_stream << A.s; break;
//
     case XML_dataType::TYPE_BOOL_PTR     : out_stream << *((bool*)(A.v)); break;
     case XML_dataType::TYPE_INT_PTR      : out_stream << *((int*)(A.v)); break;
     case XML_dataType::TYPE_LONG_PTR     : out_stream << *((long*)(A.v)); break;
     case XML_dataType::TYPE_FLOAT_PTR    : out_stream << *((float*)(A.v)); break;
     case XML_dataType::TYPE_DOUBLE_PTR   : out_stream << *((double*)(A.v)); break;
     case XML_dataType::TYPE_STRING_PTR   : out_stream << *((string*)(A.v)); break;

     }
     return out_stream;
}
istream&  operator >>(istream& in_stream, XML_dataType& A)
{
//
//   Input
//
     switch(A.type_name)
     {
	 case XML_dataType::TYPE_NULL     : XML_dataType::nullOperand(); break;
	 case XML_dataType::TYPE_BOOL     : in_stream >> A.b; break;
     case XML_dataType::TYPE_INT      : in_stream >> A.i; break;
     case XML_dataType::TYPE_LONG     : in_stream >> A.l; break;
     case XML_dataType::TYPE_FLOAT    : in_stream >> A.f; break;
     case XML_dataType::TYPE_DOUBLE   : in_stream >> A.d; break;
     case XML_dataType::TYPE_STRING   : in_stream >> A.s; break;
//
     case XML_dataType::TYPE_BOOL_PTR     : in_stream >> *((bool*)(A.v)); break;
     case XML_dataType::TYPE_INT_PTR      : in_stream >> *((int*)(A.v)); break;
     case XML_dataType::TYPE_LONG_PTR     : in_stream >> *((long*)(A.v)); break;
     case XML_dataType::TYPE_FLOAT_PTR    : in_stream >> *((float*)(A.v)); break;
     case XML_dataType::TYPE_DOUBLE_PTR   : in_stream >> *((double*)(A.v)); break;
     case XML_dataType::TYPE_STRING_PTR   : in_stream >> *((string*)(A.v)); break;

     }
     return in_stream;
}
//
//********************************************************************************
//                    MEMBER_FUNCTIONS
//********************************************************************************
//

XML_dataType::operator bool()
{
     bool b_return = 0;

     switch(type_name)
     {
     case TYPE_BOOL     : b_return =     b; break;
     case TYPE_INT      : if(i == 0){b_return = false;} else {b_return = true;} break;
     case TYPE_LONG     : if(l == 0){b_return = false;} else {b_return = true;} break;
     case TYPE_FLOAT    : XML_dataType::illegalConversion(TYPE_FLOAT,   TYPE_BOOL); break;
     case TYPE_DOUBLE   : XML_dataType::illegalConversion(TYPE_DOUBLE,  TYPE_BOOL); break;
     case TYPE_STRING   : XML_dataType::illegalConversion(TYPE_STRING,  TYPE_BOOL); break;
     case TYPE_NULL     : XML_dataType::nullOperand(); break;

     case TYPE_BOOL_PTR     : b_return = *((bool*)(v)); break;
     case TYPE_INT_PTR      : if(*((int*)(v))  == 0){b_return = false;} else {b_return = true;} break;
     case TYPE_LONG_PTR     : if(*((long*)(v)) == 0){b_return = false;} else {b_return = true;} break;
     case TYPE_FLOAT_PTR    : XML_dataType::illegalConversion(); break;
     case TYPE_DOUBLE_PTR   : XML_dataType::illegalConversion(); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalConversion(); break;
     }
     return b_return;
}

XML_dataType::operator int()
{
     int i_return = 0;

     switch(type_name)
     {
     case TYPE_BOOL     : i_return = int(i); break;
     case TYPE_INT      : i_return =     i; break;
     case TYPE_LONG     : i_return = int(l); break;
     case TYPE_FLOAT    : i_return = int(f); break;
     case TYPE_DOUBLE   : i_return = int(d); break;
     case TYPE_STRING   : XML_dataType::illegalConversion(TYPE_STRING, TYPE_INT); break;
     case TYPE_NULL     : XML_dataType::nullOperand(); break;

     case TYPE_BOOL_PTR     : i_return = *((bool*)(v)); break;
     case TYPE_INT_PTR      : i_return = *((int*)(v)); break;
     case TYPE_LONG_PTR     : i_return = int(*((long*)(v))); break;
     case TYPE_FLOAT_PTR    : i_return = int(*((float*)(v))); break;
     case TYPE_DOUBLE_PTR   : i_return = int(*((double*)(v))); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalConversion(); break;
     }
     return i_return;
}

XML_dataType::operator long()
{
     long l_return = 0;

     switch(type_name)
     {
     case TYPE_BOOL     : l_return = long(i); break;
     case TYPE_INT      : l_return = long(i); break;
     case TYPE_LONG     : l_return =       l; break;
     case TYPE_FLOAT    : l_return = long(f); break;
     case TYPE_DOUBLE   : l_return = long(d); break;
     case TYPE_STRING   : XML_dataType::illegalConversion(TYPE_STRING,TYPE_LONG); break;
     case TYPE_NULL     : XML_dataType::nullOperand(); break;

     case TYPE_BOOL_PTR     : l_return = long(*((bool*)(v))); break;
     case TYPE_INT_PTR      : l_return = long(*((int*)(v))); break;
     case TYPE_LONG_PTR     : l_return = (*((long*)(v))); break;
     case TYPE_FLOAT_PTR    : l_return = long(*((float*)(v))); break;
     case TYPE_DOUBLE_PTR   : l_return = long(*((double*)(v))); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalConversion(); break;

     }
     return l_return;
}

XML_dataType::operator float()
{
     float f_return =0.0;

     switch(type_name)
     {
     case TYPE_BOOL     : XML_dataType::illegalConversion(TYPE_BOOL,TYPE_FLOAT); break;
     case TYPE_INT      : f_return = float(i); break;
     case TYPE_LONG     : f_return = float(l); break;
     case TYPE_FLOAT    : f_return =        f; break;
     case TYPE_DOUBLE   : f_return = float(d); break;
     case TYPE_STRING   : XML_dataType::illegalConversion(TYPE_STRING,TYPE_FLOAT); break;
     case TYPE_NULL     : XML_dataType::nullOperand(); break;

     case TYPE_BOOL_PTR     : XML_dataType::illegalConversion(); break;
     case TYPE_INT_PTR      : f_return = float(*((int*)(v))); break;
     case TYPE_LONG_PTR     : f_return = float(*((long*)(v))); break;
     case TYPE_FLOAT_PTR    : f_return = (*((float*)(v))); break;
     case TYPE_DOUBLE_PTR   : f_return = float(*((double*)(v))); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalConversion(); break;

     }
     return f_return;
}

XML_dataType::operator double()
{
     double d_return = 0.0;

     switch(type_name)
     {
     case TYPE_BOOL     : XML_dataType::illegalConversion(TYPE_BOOL,TYPE_DOUBLE); break;
     case TYPE_INT      : d_return = double(i); break;
     case TYPE_LONG     : d_return = double(l); break;
     case TYPE_FLOAT    : d_return = double(f); break;
     case TYPE_DOUBLE   : d_return =        d ; break;
     case TYPE_STRING   : XML_dataType::illegalConversion(TYPE_STRING,TYPE_DOUBLE); break;
     case TYPE_NULL     : XML_dataType::nullOperand(); break;

     case TYPE_BOOL_PTR     : XML_dataType::illegalConversion(); break;
     case TYPE_INT_PTR      : d_return = double(*((int*)(v))); break;
     case TYPE_LONG_PTR     : d_return = double(*((long*)(v))); break;
     case TYPE_FLOAT_PTR    : d_return = double(*((float*)(v))); break;
     case TYPE_DOUBLE_PTR   : d_return = (*((double*)(v))); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalConversion(); break;

     }
     return d_return;
}

XML_dataType::operator string()
{

     string string_return;

     switch(type_name)
     {
     case TYPE_BOOL     : XML_dataType::illegalConversion(TYPE_BOOL,   TYPE_STRING); break;
     case TYPE_INT      : XML_dataType::illegalConversion(TYPE_INT,    TYPE_STRING); break;
     case TYPE_LONG     : XML_dataType::illegalConversion(TYPE_LONG,   TYPE_STRING); break;
     case TYPE_FLOAT    : XML_dataType::illegalConversion(TYPE_FLOAT,  TYPE_STRING); break;
	 case TYPE_DOUBLE   : XML_dataType::illegalConversion(TYPE_DOUBLE, TYPE_STRING); break;
	 case TYPE_STRING   : string_return = s; break;
	 case TYPE_NULL     : XML_dataType::nullOperand(); break;

     case TYPE_BOOL_PTR  : XML_dataType::illegalConversion(); break;
     case TYPE_INT_PTR   : XML_dataType::illegalConversion(); break;
     case TYPE_LONG_PTR  : XML_dataType::illegalConversion(); break;
     case TYPE_FLOAT_PTR : XML_dataType::illegalConversion(); break;
     case TYPE_DOUBLE_PTR: XML_dataType::illegalConversion(); break;
     case TYPE_STRING_PTR: string_return = (*((string*)(v))); break;
     }
     return string_return;
}

/*
XML_dataType::operator const char*()
{
	 const char* charStar_return = 0;

     switch(type_name)
	 {
     case TYPE_BOOL     : XML_dataType::illegalConversion(); break;
	 case TYPE_INT      : XML_dataType::illegalConversion(); break;
	 case TYPE_LONG     : XML_dataType::illegalConversion(); break;
     case TYPE_FLOAT    : XML_dataType::illegalConversion(); break;
     case TYPE_DOUBLE   : XML_dataType::illegalConversion(); break;
	 case TYPE_STRING   : charStar_return = s.c_str(); break;
     case TYPE_NULL     : XML_dataType::nullOperand(); break;

     case TYPE_BOOL_PTR  : XML_dataType::illegalConversion(); break;
     case TYPE_INT_PTR   : XML_dataType::illegalConversion(); break;
     case TYPE_LONG_PTR  : XML_dataType::illegalConversion(); break;
     case TYPE_FLOAT_PTR : XML_dataType::illegalConversion(); break;
	 case TYPE_DOUBLE_PTR: XML_dataType::illegalConversion(); break;
	 case TYPE_STRING_PTR: charStar_return = ((string*)(v))->c_str(); break;
     }
	 return charStar_return;
}
*/
//
//******************************************************************************
//                    MEMBER_FUNCTIONS
//******************************************************************************
//


void  XML_dataType::operator =(bool A)
{
     if(type_name == TYPE_NULL) type_name = TYPE_BOOL;

     switch(type_name)
     {
     case TYPE_BOOL     : b = A; break;
     case TYPE_INT      : i = int(A); break;
     case TYPE_LONG     : l = long(A); break;
     case TYPE_FLOAT    : XML_dataType::illegalConversion(TYPE_BOOL,TYPE_FLOAT); break;
     case TYPE_DOUBLE   : XML_dataType::illegalConversion(TYPE_BOOL,TYPE_DOUBLE); break;
     case TYPE_STRING   : XML_dataType::illegalAssignment(TYPE_BOOL,TYPE_STRING); break;

     case TYPE_BOOL_PTR     : *((bool*)(v))   = A; break;
     case TYPE_INT_PTR      : *((int*)(v))    = int(A); break;
     case TYPE_LONG_PTR     : *((long*)(v))   = long(A); break;
     case TYPE_FLOAT_PTR    : XML_dataType::illegalConversion(); break;
     case TYPE_DOUBLE_PTR   : XML_dataType::illegalConversion(); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalAssignment(); break;

     }
}


void  XML_dataType::operator =(int A)
{
     if(type_name == TYPE_NULL) type_name = TYPE_INT;

     switch(type_name)
     {
     case TYPE_BOOL     : if(A == 0){b = false;} else {b = true;} break;
     case TYPE_INT      : i = A; break;
     case TYPE_LONG     : l = A; break;
     case TYPE_FLOAT    : f = (float)A; break;
     case TYPE_DOUBLE   : d = A; break;
     case TYPE_STRING   : XML_dataType::illegalAssignment(TYPE_STRING,TYPE_INT); break;

     case TYPE_BOOL_PTR     : if(A == 0){*((bool*)(v)) = false;} else {*((bool*)(v)) = true;} break;
     case TYPE_INT_PTR      : *((int*)(v))    = A; break;
     case TYPE_LONG_PTR     : *((long*)(v))   = long(A); break;
     case TYPE_FLOAT_PTR    : *((float*)(v))  = float(A); break;
     case TYPE_DOUBLE_PTR   : *((double*)(v)) = double(A); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalAssignment(); break;

     }
}

void  XML_dataType::operator =(long A)
{
     if(type_name == TYPE_NULL) type_name = TYPE_LONG;

     switch(type_name)
     {
     case TYPE_BOOL     : if(A == 0){b = false;} else {b = true;} break;
     case TYPE_INT      : i = int(A); break;
     case TYPE_LONG     : l = A; break;
     case TYPE_FLOAT    : f = float(A); break;
     case TYPE_DOUBLE   : d = A; break;
     case TYPE_STRING   : XML_dataType::illegalAssignment(TYPE_STRING,TYPE_LONG); break;

     case TYPE_BOOL_PTR     : if(A == 0){*((bool*)(v)) = false;} else {*((bool*)(v)) = true;} break;
     case TYPE_INT_PTR      : *((int*)(v))    = int(A); break;
     case TYPE_LONG_PTR     : *((long*)(v))   = A; break;
     case TYPE_FLOAT_PTR    : *((float*)(v))  = float(A); break;
     case TYPE_DOUBLE_PTR   : *((double*)(v)) = double(A); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalAssignment(); break;

     }
}

void  XML_dataType::operator =(float A)
{
     if(type_name == TYPE_NULL) type_name = TYPE_FLOAT;

     switch(type_name)
     {
     case TYPE_BOOL     : XML_dataType::illegalAssignment(TYPE_FLOAT,TYPE_BOOL); break;
     case TYPE_INT      : i = int(A); break;
     case TYPE_LONG     : l = long(A); break;
     case TYPE_FLOAT    : f = A; break;
     case TYPE_DOUBLE   : d = A; break;
     case TYPE_STRING   : XML_dataType::illegalAssignment(TYPE_FLOAT,TYPE_STRING); break;

     case TYPE_BOOL_PTR     :  XML_dataType::illegalAssignment(); break;
     case TYPE_INT_PTR      : *((int*)(v))    = int(A); break;
     case TYPE_LONG_PTR     : *((long*)(v))   = long(A); break;
     case TYPE_FLOAT_PTR    : *((float*)(v))  = A; break;
     case TYPE_DOUBLE_PTR   : *((double*)(v)) = double(A); break;
     case TYPE_STRING_PTR   : XML_dataType::illegalAssignment(); break;

     }
}

void  XML_dataType::operator =(double A)
{
     if(type_name == TYPE_NULL) type_name = TYPE_DOUBLE;

     switch(type_name)
     {
     case TYPE_BOOL     : XML_dataType::illegalAssignment(TYPE_DOUBLE,TYPE_BOOL); break;
     case TYPE_INT      : i = int(A); break;
     case TYPE_LONG     : l = long(A); break;
     case TYPE_FLOAT    : f = float(A); break;
     case TYPE_DOUBLE   : d = A; break;
     case TYPE_STRING   : XML_dataType::illegalAssignment(TYPE_DOUBLE,TYPE_STRING); break;

     case TYPE_BOOL_PTR     :  XML_dataType::illegalAssignment(); break;
     case TYPE_INT_PTR      : *((int*)(v))    = int(A); break;
     case TYPE_LONG_PTR     : *((long*)(v))   = long(A); break;
     case TYPE_FLOAT_PTR    : *((float*)(v))  = float(A); break;
     case TYPE_DOUBLE_PTR   : *((double*)(v)) = A; break;
     case TYPE_STRING_PTR   : XML_dataType::illegalAssignment(); break;

     }
}

void  XML_dataType::operator =(const string& A)
{
     if(type_name == TYPE_NULL) type_name = TYPE_STRING;

     switch(type_name)
     {
     case TYPE_BOOL     : XML_dataType::illegalAssignment(TYPE_BOOL,TYPE_STRING); break;
     case TYPE_INT      : XML_dataType::illegalAssignment(TYPE_INT,TYPE_STRING); break;
     case TYPE_LONG     : XML_dataType::illegalAssignment(TYPE_LONG,TYPE_STRING); break;
     case TYPE_FLOAT    : XML_dataType::illegalAssignment(TYPE_FLOAT,TYPE_STRING); break;
     case TYPE_DOUBLE   : XML_dataType::illegalAssignment(TYPE_DOUBLE,TYPE_STRING ); break;
     case TYPE_STRING   : s = A; break;

     case TYPE_BOOL_PTR     :  XML_dataType::illegalAssignment(); break;
     case TYPE_INT_PTR      : XML_dataType::illegalAssignment(); break;
	 case TYPE_LONG_PTR     : XML_dataType::illegalAssignment(); break;
	 case TYPE_FLOAT_PTR    : XML_dataType::illegalAssignment(); break;
	 case TYPE_DOUBLE_PTR   : XML_dataType::illegalAssignment(); break;
     case TYPE_STRING_PTR   : *((string*)(v)) = A; break;

     }
}


void  XML_dataType::operator =(const char* A)
{
     if(type_name == TYPE_NULL) type_name = TYPE_STRING;

     switch(type_name)
	 {
     case TYPE_BOOL     : XML_dataType::illegalAssignment(TYPE_BOOL,TYPE_STRING); break;
     case TYPE_INT      : XML_dataType::illegalAssignment(TYPE_INT,TYPE_STRING); break;
     case TYPE_LONG     : XML_dataType::illegalAssignment(TYPE_LONG,TYPE_STRING); break;
     case TYPE_FLOAT    : XML_dataType::illegalAssignment(TYPE_FLOAT,TYPE_STRING); break;
     case TYPE_DOUBLE   : XML_dataType::illegalAssignment(TYPE_DOUBLE,TYPE_STRING ); break;
	 case TYPE_STRING   : s = A; break;

	 case TYPE_BOOL_PTR     :  XML_dataType::illegalAssignment(); break;
     case TYPE_INT_PTR      : XML_dataType::illegalAssignment(); break;
     case TYPE_LONG_PTR     : XML_dataType::illegalAssignment(); break;
	 case TYPE_FLOAT_PTR    : XML_dataType::illegalAssignment(); break;
	 case TYPE_DOUBLE_PTR   : XML_dataType::illegalAssignment(); break;
     case TYPE_STRING_PTR   : *((string*)(v)) = A; break;
     }
}
 

//
//******************************************************************************
//                    MEMBER_FUNCTIONS
//******************************************************************************
//
/*
      case TYPE_BOOL     : XML_dataType::illegalAssignment(); break;
     case TYPE_INT      : i = int(A); break;
     case TYPE_LONG     : l = long(A); break;
     case TYPE_FLOAT    : f = A; break;
     case TYPE_DOUBLE   : d = A; break;
     case TYPE_STRING   : XML_dataType::illegalAssignment(); break;
 */

string XML_dataType::toString()
{
	 std::ostringstream oss(ostringstream::out);
	 string string_return = "NULL";
     switch(type_name)
     {
     case TYPE_BOOL     : if(b == true) {string_return = "true"; } else {string_return = "false";}  break;
     case TYPE_INT      : oss << i; string_return = oss.str(); break;
     case TYPE_LONG     : oss << l; string_return = oss.str(); break;
     case TYPE_FLOAT    : oss.setf(ios::scientific); oss.precision(8);   oss << f; string_return = oss.str(); break;
     case TYPE_DOUBLE   : oss.setf(ios::scientific); oss.precision(15);  oss << d; string_return = oss.str(); break;
     case TYPE_STRING   : string_return = s; break;
     }
     return string_return;
}
void  XML_dataType::setType(int XML_dataTypeType)
{
    type_name = XML_dataTypeType;
}

int  XML_dataType::getType()
{
    return type_name;
}

bool XML_dataType::isString()
{
	if(type_name == TYPE_STRING) return true;
	return false;
}
//
//******************************************************************************
//              	ERROR HANDLING
//******************************************************************************
//
void XML_dataType::illegalConversion()
{
	const char* ErrMsg =
	"XML_parameterList Class Error :\n\nIllegal Data Type Conversion";
	cerr << ErrMsg << endl << endl << endl;
	cerr << " Fatal Error " << endl;
	exit(1);
}
void XML_dataType::illegalConversion(int typeA, int typeB)
{
	string typeAstring;
	string typeBstring;
	switch(typeA)
	{
     case TYPE_BOOL     : typeAstring.assign("bool")  ; break;
     case TYPE_INT      : typeAstring.assign("int")   ; break;
     case TYPE_LONG     : typeAstring.assign("long")  ; break;
     case TYPE_FLOAT    : typeAstring.assign("float") ; break;
     case TYPE_DOUBLE   : typeAstring.assign("double"); break;
     case TYPE_STRING   : typeAstring.assign("string"); break;
     case TYPE_NULL     : typeAstring.assign("null" ) ; break;
    }
    switch(typeB)
	{
     case TYPE_BOOL     : typeBstring.assign("bool")   ; break;
     case TYPE_INT      : typeBstring.assign("int")    ; break;
     case TYPE_LONG     : typeBstring.assign("long")   ; break;
     case TYPE_FLOAT    : typeBstring.assign( "float") ; break;
     case TYPE_DOUBLE   : typeBstring.assign("double") ; break;
     case TYPE_STRING   : typeBstring.assign("string") ; break;
     case TYPE_NULL     : typeBstring.assign( "null")  ; break;
    }
    string ErrMsg =
	"XML_parameterList Class Error :\n\nIllegal Data Type Conversion\n";
	ErrMsg.append("Attempting to convert a " + typeAstring + " to a " + typeBstring +"\n");
	cerr << ErrMsg.c_str() << endl << endl << endl;
	cerr << " Fatal Error " << endl;
	exit(1);
}

void XML_dataType::nullOperand()
{
	const char* ErrMsg =
	"XML_parameterList Class Error :\n\nNull Operand ";
	cerr << ErrMsg << endl << endl << endl;
	cerr << " Fatal Error " << endl;
	exit(1);
}

void XML_dataType::illegalAssignment()
{
	const char* ErrMsg =
	"XML_parameterList Class Error :\n\nIllegal Assignment ";
	cerr << ErrMsg << endl << endl << endl;
	cerr << " Fatal Error " << endl;
	exit(1);
}


void XML_dataType::illegalAssignment(int typeA, int typeB)
{
	string typeAstring;
	string typeBstring;
	switch(typeA)
	{
     case TYPE_BOOL     : typeAstring.assign("bool")  ; break;
     case TYPE_INT      : typeAstring.assign("int")   ; break;
     case TYPE_LONG     : typeAstring.assign("long")  ; break;
     case TYPE_FLOAT    : typeAstring.assign("float") ; break;
     case TYPE_DOUBLE   : typeAstring.assign("double"); break;
     case TYPE_STRING   : typeAstring.assign("string"); break;
     case TYPE_NULL     : typeAstring.assign("null" ) ; break;
    }
    switch(typeB)
	{
     case TYPE_BOOL     : typeBstring.assign("bool")   ; break;
     case TYPE_INT      : typeBstring.assign("int")    ; break;
     case TYPE_LONG     : typeBstring.assign("long")   ; break;
     case TYPE_FLOAT    : typeBstring.assign( "float") ; break;
     case TYPE_DOUBLE   : typeBstring.assign("double") ; break;
     case TYPE_STRING   : typeBstring.assign("string") ; break;
     case TYPE_NULL     : typeBstring.assign( "null")  ; break;
    }
    string ErrMsg =
	"XML_parameterList Class Error :\n\nIllegal Assignment \n";
	ErrMsg.append("Attempting to assign a  " + typeAstring + " to a " + typeBstring +"\n");
	cerr << ErrMsg.c_str() << endl << endl << endl;
	cerr << " Fatal Error " << endl;
	exit(1);
}
