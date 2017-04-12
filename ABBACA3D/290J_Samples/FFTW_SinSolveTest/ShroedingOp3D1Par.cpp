/*
 * ShroedingOp3D1Par.cpp
 *
 *  Created on: October 7, 2015
 *      Author: Juan Carlos Ramirez, based on andersons 2D version
 */

//#include "UCLAQ_GridFunction3d.h"


using namespace UCLAQ;


#ifndef _ShroedingOp3D_
#define _ShroedingOp3D_
#include "LaplaceOp3D.cpp"

/**

 *                       Class ShroedingOp3D
 *  A class whose apply member function implements the 7-point finite difference
 *  discete  Schoedinger operator with homogeneous boundary conditions AND potential V given by a gridfunction.
 */

class ShroedingOp3D
{
public:

    ShroedingOp3D()
    {
    initialize();
    }

    virtual ~ ShroedingOp3D()
    {}

    ShroedingOp3D(LaplaceOp3D LapOp3D,GridFunction3d Potential)//, long xPanels, double xMax,double xMin,long yPanels,double yMax,double yMin,long zPanels,double zMax, double zMin)
    {
    initialize(LapOp3D,Potential);//, xPanels xMax,xMin,yPanels, yMax, yMin, zPanels,zMax,zMin);
    }

    void initialize(LaplaceOp3D LapOp3D,GridFunction3d Potential)//, long xPanels, double xMax,double xMin,long yPanels,double yMax,double yMin,long zPanels,double zMax, double zMin)
    {
this->LapOp3D  = LapOp3D;
this->Potential=Potential;
/*    
this->xPanels = xPanels;
    this->nx  = xPanels-1;
	
	
	// quantities in the y direction

    this->yPanels = yPanels;
    this->ny  = yPanels-1;
		//quantities in the z direction
 
    this->zPanels = zPanels;
    this->nz  = zPanels-1;


    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);
//boundaries
this->xMax=xMax;
this->xMin=xMin;

this->yMax=yMax;
this->yMin=yMin;

this->zMax=zMax;
this->zMin=zMin;

    Vtmp.initialize(xPanels,xMin,xMax,yPanels,yMin, yMax,zPanels,zMin,zMax);
*/
    }

    void initialize()
    {
LapOp3D.initialize();
Potential.initialize();
/*
    long xPanels =  V.getXpanelCount();;
    this->nx  = 0;
	this->hy=0.0;
	this->yPanels=0;
	this->ny=0;
	this->hz=0.0;
	this->zPanels=0;
	this->nz=0;

this->xMax=0;
this->xMin=0;


this->yMax=0;
this->yMin=0;

this->zMax=0;
this->zMin=0;
    Vtmp.initialize();*/
    }

    /**
    This routine applies the 7-point discrete SCHROEDINGER operator to
    the interior grid points associated with a uniform discretization
    with m INTERIOR points that are spaced hx apart in the x-direction and 
	 hy in the y-direction 
	.

    Input  : M270::Vector2D a vector of size m whose values are those of
             associated with the interior values of the discretization

    Output : The input vector is overwritten with the result

    If _DEBUG is defined at compile time, bounds checking is performed
    on operations associated on the M270::Vector1D class, which is why there is no explicit
    bounds checking in this routine.
    */

    void apply(GridFunction3d& V)
    {




    GridFunction3d VtmpPot = V;
	long i;
	long j;
	long k;

 double hx = V.getHx();
  double hy = V.getHy();
  double hz = V.getHy();

    long  xPanels = V.getXpanelCount();
    long  yPanels = V.getYpanelCount();
  long zPanels=V.getZpanelCount();

    long nx  = xPanels-1;
     long ny  = yPanels-1;
long nz= zPanels-1;


    // Interior grid points not adjacent to the edge, Compute potential
	for(j = 1; j <= ny; j++){
            for(i = 1; i <= nx; i++){
		for(k=1;k<=nz;k++){

    			VtmpPot(i,j,k)=Potential(i,j,k)*V(i,j,k);
    				    
					 }
				}

    }
//V is now alpha Lap(V)
LapOp3D.apply(V);
//Schoedinger is Potential-alpha Lap(V)
VtmpPot-=V;
V=VtmpPot;
	
}
     GridFunction3d Potential;
    LaplaceOp3D       LapOp3D; // Coefficient of the Laplace Operator

/*
    double            hx; // Grid spacing in the x-direction
    long         xPanels; // Number of panels in the x-direction
    long          nx; // Number of unknowns in the x-direction
	

	// global variables in y direction 
	double	hy      ; // Grid spacing in the y-direction
    long	yPanels ; // Number of panels in the y-dir
    long	ny  ;  // Number of unknowns in the y-direction
    UCLAQ::GridFunction3d Vtmp; // Temporary 2D vector

	// global variables in y direction 
	double	hz      ; // Grid spacing in the y-direction
    long	zPanels ; // Number of panels in the y-dir
    long	nz  ;  // Number of unknowns in the y-direction

//boundary limits
double xMin,xMax,yMin,yMax,zMin,zMax;*/

};




#endif /* LAPLACEOP1D_H_ */
