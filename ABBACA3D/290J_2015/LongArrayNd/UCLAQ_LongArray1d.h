/*
 * UCLAQ_LongArray1d.h
 *
 *  Created on: Jun 7, 2016
 *      Author: anderson
 *
 *
 *  A minimal 1d long array class with move semantic implementation.
 *
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

#ifndef _UCLAQ_LongArray1d_
#define _UCLAQ_LongArray1d_

namespace UCLAQ
{
class LongArray1d
{
	public:

    LongArray1d()
	{
	dataPtr    = nullptr;
	index1Size = 0;
	}

	LongArray1d(long n)
	{
	dataPtr     = new long[n];
	index1Size = n;
	}

    LongArray1d(const LongArray1d& V)
    {
      #ifdef _VERBOSE_OPS_
      cout << "Standard Copy " << endl;
      #endif

      if(V.dataPtr == nullptr)
      {dataPtr = nullptr; index1Size = 0; return;}

      dataPtr     = new long[V.index1Size];
      index1Size = V.index1Size;
#ifdef _MSC_VER
      std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size);
#else
      std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif
   }

    LongArray1d(LongArray1d&& V)
    {
      #ifdef _VERBOSE_OPS_
      cout << "Move Copy " << endl;
      #endif

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;;
    }

   	virtual ~LongArray1d()
	{
	if(dataPtr != nullptr) delete [] dataPtr;
	}

	void initialize()
	{
	if(dataPtr != nullptr) delete [] dataPtr;
	dataPtr    = nullptr;
	index1Size = 0;
	}

	void initialize(long n)
	{
      if(index1Size != n)
      {
	  if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr     = new long[n];
      index1Size = n;
      }
	}

    void initialize(const LongArray1d& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      return;
      }

      if(index1Size != V.index1Size)
      {
	  if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr     = new long[V.index1Size];
      index1Size = V.index1Size;
      }

#ifdef _MSC_VER
      std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size);
#else
      std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif
     }

    void initialize(LongArray1d&& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      return;
      }

      if(dataPtr != nullptr) delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = 0;
      V.index1Size = 0;
    }

    // Assignment operators : Being careful with nullptr instances

    LongArray1d& operator=(const LongArray1d& V)
    {
      #ifdef _VERBOSE_OPS_
      cout << "Standard Assignment" << endl;
      #endif

      if (this != &V)
      {
         if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
         {
         index1Size  = V.index1Size;
         dataPtr     = new long[index1Size];
#ifdef _MSC_VER
         std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif
         }
         else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; return *this;}
         else
         {
         assert(sizeCheck(this->index1Size,V.index1Size));
#ifdef _MSC_VER
         std::memcpy(dataPtr, V.dataPtr, (sizeof(long))*index1Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif
         }
      }
      return *this;
    }

	LongArray1d& operator=(LongArray1d&& V)
	{
    #ifdef _VERBOSE_OPS_
    cout << "Move Assignment" << endl;
    #endif

	if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
    {
      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
    }
    else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; return *this;}
    else
    {
      assert(sizeCheck(this->index1Size,V.index1Size));

      // Remove existing data

      delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
    }
    return *this;
    }

/*!  Returns the dimension of the vector */

	virtual long getDimension()
	{
    return index1Size;
	}

#ifndef _NDEBUG
    long&  operator()(long i1)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    return *(dataPtr +  i1);
    };

    const long&  operator()(long i1) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    return *(dataPtr +  i1);
    };

#else
    inline long&  operator()(long i1)
    {
    return *(dataPtr + i1);
    };

    inline const long&  operator()(long i1) const
    {
    return *(dataPtr + i1);
    };
#endif


    long getSize()         const {return index1Size;}
    long getIndex1Size()  const {return index1Size;}

    long* getDataPointer(){return dataPtr;}
    const  long* getDataPointer()  const  {return dataPtr;}

	long*   dataPtr;
	long index1Size;


//###################################################################
//                      Bounds Checking
//###################################################################
//
#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end, int coordinate) const
        {
        if((i < begin)||(i  > end))
        {
        cerr << "UCLAQ::LongArray1d index " << coordinate << " out of bounds " << endl;
        cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << endl;
        return false;
        }
        return true;
        }
#else
        bool boundsCheck(long, long, long, int) const {return true;}
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2)
    {
    if(size1 != size2)
    {
    cerr << "UCLAQ::LongArray1d sizes are incompatible : " << size1 << " != " << size2;
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2) const
    {
    if(size1 != size2)
    {
    cerr << "UCLAQ::LongArray1d sizes are incompatible : " << size1 << " != " << size2;
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long) {return true;}
    bool sizeCheck(long, long) const{return true;}
#endif

};
}

#endif /* UCLAQ_LongArray1d_ */
