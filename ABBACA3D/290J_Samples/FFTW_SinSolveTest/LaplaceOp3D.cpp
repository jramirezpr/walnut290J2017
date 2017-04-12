/*
 * LaplaceOp2D.h
 *
 *  Created on: October 7, 2015
 *      Author: Juan Carlos Ramirez, based on andersons 2D version
 */

//#include "UCLAQ_GridFunction3d.h"
//#include "UCLAQ_GridFunction2d.h"
using namespace UCLAQ;


#ifndef _LaplaceOp3D_
#define _LaplaceOp3D_

/**

 *                       Class LaplaceOp3D
 *  A class whose apply member function implements the 7-point finite difference
 *  discete  Laplace operator with homogeneous boundary conditions.
 */

class LaplaceOp3D
{
public:

    LaplaceOp3D()
    {
    initialize();
    }

    virtual ~LaplaceOp3D()
    {}

    LaplaceOp3D(double alpha)//, long xPanels, double xMax,double xMin,long yPanels,double yMax,double yMin,long zPanels,double zMax, double zMin)
    {
    initialize(alpha);//, xPanels xMax,xMin,yPanels, yMax, yMin, zPanels,zMax,zMin);
    }

    void initialize(double alpha)//, long xPanels, double xMax,double xMin,long yPanels,double yMax,double yMin,long zPanels,double zMax, double zMin)
    {
    this->alpha   = alpha;
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
alpha   = 0.0;
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
    This routine applies the 7-point discrete Laplace operator to
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



    GridFunction3d Vtmp = V;
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
/* Boundary cases temporarily ignored
    // Grid points adjacent to edge: START WITH CORNERS

    long i;
	long j;

	//1st corner
    i = 0;
    j =0;
	V(i,j)   = alpha*(((-2.0*Vtmp(i,j) + Vtmp(i+1,j))/(hx*hx))+((-2.0*Vtmp(i,j) + Vtmp(i,j+1))/(hy*hy)));

	//2nd corner 
    i = nx-1;

    V(i,j) =   alpha*(((-2.0*Vtmp(i,j) + Vtmp(i-1,j))/(hx*hx))+((-2.0*Vtmp(i,j) + Vtmp(i,j+1))/(hy*hy)));


	// upper corners
	//3rd corner
	i=0;
	j=ny-1;
	V(i,j)   = alpha*(((-2.0*Vtmp(i,j) + Vtmp(i+1,j))/(hx*hx))+((-2.0*Vtmp(i,j) + Vtmp(i,j-1))/(hy*hy)));

	//4th corner
	i=nx-1;

	 V(i,j) =   alpha*(((-2.0*Vtmp(i,j) + Vtmp(i-1,j))/(hx*hx))+((-2.0*Vtmp(i,j) + Vtmp(i,j-1))/(hy*hy)));


	 //non-corner edges

	 //edges along the x axis
	 
	 
	 //1st (bottom) edge
	 j=0;
	   for(i = 1; i < nx-1; i++)
    {
    V(i,j) =  alpha*(((Vtmp(i+1,j) - 2.0*Vtmp(i,j) + Vtmp(i-1,j))/(hx*hx))+((-2.0*Vtmp(i,j) + Vtmp(i,j+1))/(hy*hy)));
    }

	   //2nd (top) edge
	    j=ny-1;
			   for(i = 1; i < nx-1; i++)
    {
    V(i,j) =  alpha*(((Vtmp(i+1,j) - 2.0*Vtmp(i,j) + Vtmp(i-1,j))/(hx*hx))+((-2.0*Vtmp(i,j) + Vtmp(i,j-1))/(hy*hy)));
    }

			   //now edges along y axis
	 //3rd left  (left)
	 i=0;
  for(j = 1; j < ny-1; j++)
    {
    V(i,j) =  alpha*(((-2.0*Vtmp(i,j) + Vtmp(i+1,j))/(hx*hx))+ ((Vtmp(i,j+1) - 2.0*Vtmp(i,j) + Vtmp(i,j-1))/(hy*hy)));
    }

  //4rth (right) edge
  i=nx-1;
    for(j = 1; j < ny-1; j++)
    {
    V(i,j) =  alpha*(((-2.0*Vtmp(i,j) + Vtmp(i-1,j))/(hx*hx))+ ((Vtmp(i,j+1) - 2.0*Vtmp(i,j) + Vtmp(i,j-1))/(hy*hy)));
    }

	*/

    // Interior grid points not adjacent to the edge
	for(j = 1; j <= ny; j++){
            for(i = 1; i <= nx; i++){
		for(k=1;k<=nz;k++){

    			V(i,j,k)=  alpha*(((Vtmp(i+1,j,k) - 2.0*Vtmp(i,j,k) + Vtmp(i-1,j,k))/(hx*hx))+((Vtmp(i,j+1,k) - 2.0*Vtmp(i,j,k) + Vtmp(i,j-1,k))/(hy*hy))
					+((Vtmp(i,j,k+1) - 2.0*Vtmp(i,j,k) + Vtmp(i,j,k-1))/(hz*hz)));
    				    
					 }
				}

    }
	
}

    double         alpha; // Coefficient of the Laplace Operator

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