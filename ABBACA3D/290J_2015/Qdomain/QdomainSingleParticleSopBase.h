/*
 * QdomainSingleParticleSopBase.h
 *
 * Virtual base class for single particle Schroedinger operator
 *
 *  Created on: Jul 3, 2014
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright 2015-2016 Chris Anderson
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

#ifndef  _QdomainSingleParticleSopBase_
#define  _QdomainSingleParticleSopBase_

#include "UCLAQ::GridFunction3d.h"

namespace UCLAQ
{

class QdomainSingleParticleSopBase
{
    public:
	QdomainSingleParticleSopBase(){};
	QdomainSingleParticleSopBase(const QdomainSingleParticleSopBase& Qop){};
	virtual ~QdomainSingleParticleSopBase(){};
	virtual void applyForwardOp(UCLAQ::GridFunction3d& F) = 0;
};


} // UCLAQ namespace

#endif /* _QdomainSingleParticleSopBase_ */
