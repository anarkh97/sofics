#ifndef _BARF_COROTATOR_H_
#define _BARF_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class BarFCorotator : public Corotator {
     int n1; 			// node number 1
     int n2; 			// node number 2
     double em;                 // elastic modulus
     double a0;                 // initial cross-sectional area
     double l0;			// initial length
     double preload;		// preload
     double damage;		// damage in the element
     double lambda;		// damage growth parameter
     double Ucrit;		// Stretch where Damage Initiates
     double Uf;			// Failure Stretch for the Yarn
     int op;			// Material Option

   public:

     // Constructor
     BarFCorotator(int node1, int node2, double em,
			double lambda, double area, 
			int op, double h, double d,
			double Uc, double Uf,
			int np, int Nf, double dlambda, int Seed,
			double preload, CoordSet &cs);

     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt, double t);
     void   getExternalForce(GeomState &geomState, CoordSet &cs, double *f);

     void   formInternalForce(double t[3], double p, double *f);

     void   getDExternalForceDu(GeomState &geomState, CoordSet &cs, 
    				     FullSquareMatrix &elK, double *locF);

     void   AssignMicroScaleProp(double ef, double h, double d, 
			int np, int Nf, double lambda_g);

     void   MicroCalcLambda(int np, double *pDamage, double *pU_vec, 
			double lambda_g);

     double RndNorm(double mean, double stdev);

     void   formGeometricStiffness(GeomState &gs, CoordSet &cs, 
                                   FullSquareMatrix &elk, double *f);

     void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                              int &nlflag);

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                 double *vlr);

     double getElementEnergy(GeomState &gs, CoordSet &cs);

     double getDamage() { return damage; }

};

#endif
