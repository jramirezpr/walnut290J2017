/*
 * GridFunction2dSliceOp: A utility class to facilitate thread-safe insertion and extraction
 * of coordinate slices of data associated with GridFunction2d instances.
 *
 * The only data member of this class is a pointer to a UCLAQ::GridFunction2d instance
 * to which one wants to extract coordinate slices or insert coordinate slices.
 *
 * To create thread-safe code one typically creates a vector of instances of this class that
 * are all associated with a single UCLAQ::GridFunction2d instance. One then multi-threads
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
 * UCLAQ_GridFunction2dSliceOp3d.h
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

#ifndef _UCLAQ_GridFunction2dSliceOp_
#define _UCLAQ_GridFunction2dSliceOp_

#include "UCLAQ_GridFunction1d.h"
#include "UCLAQ_GridFunction2d.h"
#include "UCLAQ_GridFunction2d.h"

using namespace std;

namespace UCLAQ
{
class GridFunction2dSliceOp
{

public:

GridFunction2dSliceOp()
{
	gridFunPtr= nullptr;
};

GridFunction2dSliceOp(GridFunction2d& R)
{
	gridFunPtr= &R;
};

GridFunction2dSliceOp(const GridFunction2dSliceOp& R)
{
	gridFunPtr= R.gridFunPtr;
};

void initialize()
{
	gridFunPtr= nullptr;
};

void initialize(GridFunction2d& R)
{
	gridFunPtr= &R;
};

void initialize(const GridFunction2dSliceOp& R)
{
	gridFunPtr= R.gridFunPtr;
};


GridFunction1d getConstantYslice(long yIndex) const  // (x function)
{
	GridFunction1d R(gridFunPtr->xPanels, gridFunPtr->xMin, gridFunPtr->xMax);

	for(long i = 0; i <= gridFunPtr->xPanels; i++)
	{
	R.Values(i) = gridFunPtr->Values(i,yIndex);
	}
	return R;
}

void getConstantYslice(long yIndex, GridFunction1d& R) const  // (x function)
{
	for(long i = 0; i <= gridFunPtr->xPanels; i++)
	{
	R.Values(i) = gridFunPtr->Values(i,yIndex);
	}
}

GridFunction1d getConstantXslice(long xIndex) const  //( y function)
{
	GridFunction1d R(gridFunPtr->yPanels,gridFunPtr->yMin,gridFunPtr->yMax);

	for(long j = 0; j <= gridFunPtr->yPanels; j++)
	{
	R.Values(j) = gridFunPtr->Values(xIndex,j);
	}
	return R;
}

void getConstantXslice(long xIndex, GridFunction1d& R) const  //( y function)
{
	for(long j = 0; j <= gridFunPtr->yPanels; j++)
	{
	R.Values(j) = gridFunPtr->Values(xIndex,j);
	}
}


///
/// Insertion
///


void setConstantYslice(long yIndex, const GridFunction1d& R) const  // (x function)
{
	for(long i = 0; i <= gridFunPtr->xPanels; i++)
	{
	gridFunPtr->Values(i,yIndex) = R.Values(i);
	}
}

void setConstantXslice(long xIndex, const GridFunction1d& R) const  //( y function)
{
	for(long j = 0; j <= gridFunPtr->yPanels; j++)
	{
	gridFunPtr->Values(xIndex,j) = R.Values(j);
	}
}


UCLAQ::GridFunction2d* gridFunPtr;
};

} // Namespace UCLAQ
#endif /* CLAQ_GridFunction2dSliceOp */
