/*
 * QdomainSchroedingerOp.h
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

#ifndef  _QdomainSingleParticleSop_
#define  _QdomainSingleParticleSop_

#include "QdomainSingleParticleSopBase.h"
#include "QdomainKineticOp.h"
#include "GridFunction3d.h"

namespace UCLAQ
{

class QdomainSingleParticleSop : public QdomainSingleParticleSopBase
{
    public:

	QdomainSingleParticleSop()
	{
	initialize();
	}

	QdomainSingleParticleSop(const QdomainKineticOp& kineticOp, const GridFunction3d& potential)
	{
	initialize(kineticOp,potential);
	}

	QdomainSingleParticleSop(const QdomainSingleParticleSop& Qop)
	{
	initialize(Qop);
	}

	void initialize()
	{
	kineticOp.initialize();
	potential.initialize();
	opTemp.initialize();
	kineticOnlyFlag = false;
	}

	void initialize(const QdomainKineticOp& kineticOp, const GridFunction3d& potential)
	{
	this->kineticOp.initialize(kineticOp);
	this->potential.initialize(potential);
	this->opTemp.initialize(potential);
	kineticOnlyFlag = false;
	}

	void initialize(const QdomainSingleParticleSop& Qop)
	{
	initialize(Qop.kineticOp,Qop.potential);
	}

    void setKineticOnly()
    {
    kineticOnlyFlag = true;
    }

    void clearKineticOnly()
    {
    kineticOnlyFlag = false;
    }

	virtual void applyForwardOp(GridFunction3d& F)
	{
	if(kineticOnlyFlag)
	{
	kineticOp.applyForwardOp(F);
	return;
	}

	opTemp  = F;
	opTemp *= potential;
	kineticOp.applyForwardOp(F);
	F     += opTemp;
	}

    QdomainKineticOp     kineticOp;
	GridFunction3d       potential;
	GridFunction3d         opTemp;
	bool          kineticOnlyFlag;
};

} // UCLAQ namespace

#endif /* _QdomainSingleParticleSop_ */
