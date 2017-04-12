/*
 * GridFunction3dSliceOp: A utility class to facilitate thread-safe insertion and extraction
 * of coordinate slices of data associated with GridFunction3d instances.
 *
 * The only data member of this class is a pointer to a UCLAQ::GridFunction3d instance
 * to which one wants to extract coordinate slices or insert coordinate slices.
 *
 * To create thread-safe code one typically creates a vector of instances of this class that
 * are all associated with a single UCLAQ::GridFunction3d instance. One then multi-threads
 * an insertion or extraction loop using the instance with index equal to the thread index
 * inside the loop. As long as the overall insertion or extraction loop is
 * indexed disjoint data of the UCLAQ::GridFunction2d, the extraction or insertion
 * will be thread-safe.
 *
 *
 * Warning: When inserting data you have to be a bit careful that one doesn't insert
 * data that "crosses", e.g. insert data slices corresponding to different coordinate
 * directions.
 *
 *
 *
 * UCLAQ_GridFunction3dSliceOp.h
 *
 *  Created on: Feb 25, 2016
 *      Author: anderson
 */

 /*
#############################################################################
#
# Copyright  2016 Chris Anderson
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
###

*/

#ifndef _UCLAQ_GridFunction3dSliceOp_
#define _UCLAQ_GridFunction3dSliceOp_

#include "UCLAQ_GridFunction1d.h"
#include "UCLAQ_GridFunction2d.h"
#include "UCLAQ_GridFunction3d.h"

using namespace std;

namespace UCLAQ
{
class GridFunction3dSliceOp
{

public:

GridFunction3dSliceOp()
{
	gridFunPtr= nullptr;
};

GridFunction3dSliceOp(GridFunction3d& R)
{
	gridFunPtr= &R;
};

GridFunction3dSliceOp(const GridFunction3dSliceOp& R)
{
	gridFunPtr= R.gridFunPtr;
};

void initialize()
{
	gridFunPtr= nullptr;
};

void initialize(GridFunction3d& R)
{
	gridFunPtr= &R;
};

void initialize(const GridFunction3dSliceOp& R)
{
	gridFunPtr= R.gridFunPtr;
};


GridFunction2d getConstantZslice(long zIndex) const //(x-y function)
{
    GridFunction2d R(gridFunPtr->xPanels, gridFunPtr->xMin, gridFunPtr->xMax, gridFunPtr->yPanels, gridFunPtr->yMin, gridFunPtr->yMax);
    for(long i = 0; i <=  gridFunPtr->xPanels; i++)
    {
    for(long j = 0; j <=  gridFunPtr->yPanels; j++)
    {
    R.Values(i,j) = gridFunPtr->Values(i,j,zIndex);
    }}
    return R;
}

void getConstantZslice(long zIndex, GridFunction2d& R) const //(x-y function)
{
    for(long i = 0; i <=  gridFunPtr->xPanels; i++)
    {
    for(long j = 0; j <=  gridFunPtr->yPanels; j++)
    {
    R.Values(i,j) = gridFunPtr->Values(i,j,zIndex);
    }}
}

GridFunction2d getConstantYslice(long yIndex) const //(x-z function)
{
	GridFunction2d R(gridFunPtr->xPanels, gridFunPtr->xMin, gridFunPtr->xMax, gridFunPtr->zPanels, gridFunPtr->zMin, gridFunPtr->zMax);

    for(long i = 0; i <= gridFunPtr->xPanels; i++)
    {
    for(long k = 0; k <= gridFunPtr->zPanels; k++)
    {
    R.Values(i,k) = gridFunPtr->Values(i,yIndex,k);
    }}
    return R;
}

void getConstantYslice(long yIndex, GridFunction2d& R) const //(x-z function)
{
    for(long i = 0; i <= gridFunPtr->xPanels; i++)
    {
    for(long k = 0; k <= gridFunPtr->zPanels; k++)
    {
    R.Values(i,k) = gridFunPtr->Values(i,yIndex,k);
    }}
}


GridFunction2d getConstantXslice(long xIndex) const //(y-z function)
{
	GridFunction2d R(gridFunPtr->yPanels, gridFunPtr->yMin, gridFunPtr->yMax, gridFunPtr->zPanels, gridFunPtr->zMin, gridFunPtr->zMax);

    for(long j = 0; j <= gridFunPtr->yPanels; j++)
    {
    for(long k = 0; k <= gridFunPtr->zPanels; k++)
    {
    R.Values(j,k) = gridFunPtr->Values(xIndex,j,k);
    }}
    return R;
}

void getConstantXslice(long xIndex,GridFunction2d& R) const //(y-z function)
{
    for(long j = 0; j <= gridFunPtr->yPanels; j++)
    {
    for(long k = 0; k <= gridFunPtr->zPanels; k++)
    {
    R.Values(j,k) = gridFunPtr->Values(xIndex,j,k);
    }}
}


GridFunction1d getConstantYZslice(long yIndex, long zIndex) const  // (x function)
{
	GridFunction1d R(gridFunPtr->xPanels, gridFunPtr->xMin, gridFunPtr->xMax);

	for(long i = 0; i <= gridFunPtr->xPanels; i++)
	{
	R.Values(i) = gridFunPtr->Values(i,yIndex,zIndex);
	}
	return R;
}

void getConstantYZslice(long yIndex, long zIndex,GridFunction1d& R) const  // (x function)
{
	for(long i = 0; i <= gridFunPtr->xPanels; i++)
	{
	R.Values(i) = gridFunPtr->Values(i,yIndex,zIndex);
	}
}

GridFunction1d getConstantXZslice(long xIndex, long zIndex) const  //( y function)
{
	GridFunction1d R(gridFunPtr->yPanels,gridFunPtr->yMin,gridFunPtr->yMax);

	for(long j = 0; j <= gridFunPtr->yPanels; j++)
	{
	R.Values(j) = gridFunPtr->Values(xIndex,j,zIndex);
	}
	return R;
}

void getConstantXZslice(long xIndex, long zIndex,GridFunction1d& R) const  //( y function)
{
	for(long j = 0; j <= gridFunPtr->yPanels; j++)
	{
	R.Values(j) = gridFunPtr->Values(xIndex,j,zIndex);
	}
}


GridFunction1d getConstantXYslice(long xIndex, long yIndex) const  //( z function)
{
	GridFunction1d R(gridFunPtr->zPanels,gridFunPtr->zMin,gridFunPtr->zMax);

	for(long k = 0; k <= gridFunPtr->zPanels; k++)
	{
	R.Values(k) = gridFunPtr->Values(xIndex,yIndex,k);
	}
	return R;
}

void  getConstantXYslice(long xIndex, long yIndex, GridFunction1d& R) const  //( z function)
{
	for(long k = 0; k <= gridFunPtr->zPanels; k++)
	{
	R.Values(k) = gridFunPtr->Values(xIndex,yIndex,k);
	}
}

///
/// Insertion
///

void setConstantZslice(long zIndex, const GridFunction2d& R) const //(x-y function)
{
    for(long i = 0; i <=  gridFunPtr->xPanels; i++)
    {
    for(long j = 0; j <=  gridFunPtr->yPanels; j++)
    {
    gridFunPtr->Values(i,j,zIndex) = R.Values(i,j);
    }}
}


void setConstantYslice(long yIndex, const GridFunction2d& R) const //(x-z function)
{
    for(long i = 0; i <= gridFunPtr->xPanels; i++)
    {
    for(long k = 0; k <= gridFunPtr->zPanels; k++)
    {
    gridFunPtr->Values(i,yIndex,k) = R.Values(i,k);
    }}
}

void setConstantXslice(long xIndex,const GridFunction2d& R) const //(y-z function)
{
    for(long j = 0; j <= gridFunPtr->yPanels; j++)
    {
    for(long k = 0; k <= gridFunPtr->zPanels; k++)
    {
    gridFunPtr->Values(xIndex,j,k) = R.Values(j,k);
    }}
}


void setConstantYZslice(long yIndex, long zIndex,const GridFunction1d& R) const  // (x function)
{
	for(long i = 0; i <= gridFunPtr->xPanels; i++)
	{
	gridFunPtr->Values(i,yIndex,zIndex) = R.Values(i);
	}
}

void setConstantXZslice(long xIndex, long zIndex,const GridFunction1d& R) const  //( y function)
{
	for(long j = 0; j <= gridFunPtr->yPanels; j++)
	{
	gridFunPtr->Values(xIndex,j,zIndex) = R.Values(j);
	}
}


void  setConstantXYslice(long xIndex, long yIndex, const GridFunction1d& R) const  //( z function)
{
	for(long k = 0; k <= gridFunPtr->zPanels; k++)
	{
	gridFunPtr->Values(xIndex,yIndex,k) = R.Values(k);
	}
}

UCLAQ::GridFunction3d* gridFunPtr;
};

} // Namespace UCLAQ
#endif /* CLAQ_GridFunction3dSliceOp */
