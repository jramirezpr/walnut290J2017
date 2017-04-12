/*
 * UCLAQ_LongArray2d.h
 *
 *  Created on: Jun 7, 2016
 *      Author: anderson
 *
 *  A minimal 2d long array class with move semantic implementation. 
*/
/*
#############################################################################
#
# Copyright 2015 Chris Anderson
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


#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
using namespace std;

#ifdef  _DEBUG
#include <iostream>
#include <cstdio>
#else
#define _NDEBUG
#endif
#include <cassert>

#undef _VERBOSE_OPS_

#ifndef _UCLAQ_LongArray2d_
#define _UCLAQ_LongArray2d_

namespace UCLAQ
{
class LongArray2d
{
	public:

    LongArray2d()
	{
	dataPtr    = nullptr;
	index1Size = 0;
	index2Size = 0;
	}

	explicit LongArray2d(long m, long n)
	{
	dataPtr    = new long[m*n];
	index1Size = m;
	index2Size = n;
	}

   LongArray2d(const LongArray2d& V)
   {
      #ifdef _VERBOSE_OPS_
      cout << "Standard Copy " << endl;
      #endif

      if(V.dataPtr == nullptr) {dataPtr = nullptr; index1Size = 0; index2Size = 0; return;}

      dataPtr     = new long[V.index1Size*V.index2Size];
      index1Size = V.index1Size;
      index2Size = V.index2Size;
#ifdef _MSC_VER
      std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size*index2Size);
#else
      std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif
   }

    LongArray2d(LongArray2d&& V)
    {
      #ifdef _VERBOSE_OPS_
      cout << "Move Copy " << endl;
      #endif

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }

   	virtual ~LongArray2d()
	{
	if(dataPtr != nullptr) delete [] dataPtr;
	}

	void initialize()
	{
	if(dataPtr != nullptr) delete [] dataPtr;
	dataPtr    = nullptr;
	index1Size = 0;
	index2Size = 0;
	}

	void initialize(long m, long n)
	{
		if((index1Size != m)||(index2Size != n))
    	{
    		if(dataPtr != nullptr) delete [] dataPtr;
			dataPtr    = new long[n*m];
			index1Size = m;
			index2Size = n;
    	}
	}

    void initialize(const LongArray2d& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      index2Size = 0;
      return;
      }

      if((index1Size != V.index1Size)||(index2Size != V.index2Size))
      {
	  if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr     = new long[V.index1Size*V.index2Size];
      index1Size = V.index1Size;
      index2Size = V.index2Size;
      }

#ifdef _MSC_VER
      std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size*index2Size);
#else
      std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif
    }

    void initialize(LongArray2d&& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      index2Size = 0;
      return;
      }

      if(dataPtr != nullptr) delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }

    // Assignment operators : Being careful with nullptr instances

    LongArray2d& operator=(const LongArray2d& V)
    {
      #ifdef _VERBOSE_OPS_
      cout << "Standard Assignment" << endl;
      #endif

      if (this != &V)
      {
         if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
         {
         index1Size  = V.index1Size;
         index2Size  = V.index2Size;
         dataPtr     = new long[index1Size*index2Size];
#ifdef _MSC_VER
         std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size*index2Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif
         }
         else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; index2Size = 0; return *this;}
         else
         {
         assert(sizeCheck(this->index1Size,V.index1Size,1));
         assert(sizeCheck(this->index2Size,V.index2Size,2));
#ifdef _MSC_VER
         std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size*index2Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif
         }
      }
      return *this;
    }

	LongArray2d& operator=(LongArray2d&& V)
	{
    #ifdef _VERBOSE_OPS_
    cout << "Move Assignment" << endl;
    #endif

	if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
    {
      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }
    else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; index2Size = 0; return *this;}
    else
    {
      assert(sizeCheck(this->index1Size,V.index1Size,1));
      assert(sizeCheck(this->index2Size,V.index2Size,2));

      // Remove existing data

      delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }
    return *this;
    }

/*!  Sets all values of the vector to d. */
	void setToValue(long d)
	{
	for(long i = 0; i < index1Size*index2Size; i++)
	{dataPtr[i] =  d;}
	}

    void addValue(long d)
	{
	for(long i = 0; i < index1Size*index2Size; i++)
	{dataPtr[i] += d;}
	}


/*!  Returns the dimension of the vector */

virtual long getDimension()
{
    return index1Size*index2Size;
}


//
//###################################################################
//      Bounds checking enforced by defining  _DEBUG
//###################################################################
//

#ifdef _DEBUG
    long&  operator()(long i1, long i2)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +  i2 + i1*index2Size);
    };

    const long&  operator()(long i1, long i2) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +   i2  + i1*index2Size);
    };
#else
    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
    */
    inline long&  operator()(long i1, long i2)
    {
    return *(dataPtr +  i2 + i1*index2Size);
    };

    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
     */
    inline const long&  operator()(long i1, long i2) const
    {
    return *(dataPtr +  i2  + i1*index2Size);
    };
#endif


    long getIndex1Size()  const {return index1Size;}
    long getIndex2Size()  const {return index2Size;}


    long* getDataPointer(){return dataPtr;}
    const  long* getDataPointer() const {return dataPtr;};

	long*       dataPtr;
	long index1Size;
    long index2Size;


//###################################################################
//                      Bounds Checking
//###################################################################
//
#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end, int coordinate) const
        {
        if((i < begin)||(i  > end))
        {
        cerr << "UCLAQ::LongArray2d index " << coordinate << " out of bounds " << endl;
        cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << endl;
        return false;
        }
        return true;
        }
#else
        bool boundsCheck(long, long, long, int) const {return true;}
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2, int coordinate)
    {
    if(size1 != size2)
    {
    cerr << "UCLAQ::LongArray2d sizes are incompatible : " << size1 << " != " << size2;
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2, int coordinate) const
    {
    if(size1 != size2)
    {
    cerr << "UCLAQ::LongArray2d sizes are incompatible : " << size1 << " != " << size2;
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long, int) {return true;}
    bool sizeCheck(long, long, int) const{return true;}
#endif

};
}

#endif /* UCLAQ_LongArray2d_ */
