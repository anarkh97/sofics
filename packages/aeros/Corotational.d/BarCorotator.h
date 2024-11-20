#ifndef _BAR_COROTATOR_H_
#define _BAR_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class BarCorotator : public Corotator {
     int n1; 			// node number 1
     int n2; 			// node number 2
     double preload;		// preload
     StructProp *prop;
   public:

     // Constructor
     BarCorotator(int node1, int node2, StructProp *_prop, double preload, CoordSet &cs);
     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt, double t);
     void   getInternalForce(GeomState &geomState, CoordSet &cs, 
                                  FullSquareMatrix &elK, double *f, double dt, double t);
     void   getExternalForce(GeomState &geomState, CoordSet &cs, double *f);

     void   formInternalForce(double t[3], double p, double *f);

     void   getDExternalForceDu(GeomState &geomState, CoordSet &cs, 
    				     FullSquareMatrix &elK, double *locF);

     void   formTangentStiffness(double t[3], double p, double ld, double l0, 
                                 double kt[6][6]);

     void   formGeometricStiffness(GeomState &gs, CoordSet &cs, 
                                   FullSquareMatrix &elk, double *f);

     void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                              int &nlflag);

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                 double *vlr);

     double getElementEnergy(GeomState &gs, CoordSet &cs);
     void   getLocalDisp(GeomState& , CoordSet &, Vector&);	
     void   getGlobalDisp(GeomState& , CoordSet &, Vector&);

};

#endif
