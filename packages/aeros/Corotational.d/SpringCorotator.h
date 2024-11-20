#ifndef _SPRING_COROTATOR_H_
#define _SPRING_COROTATOR_H_

#include <Corotational.d/Corotator.h>
#include <Math.d/FullSquareMatrix.h>

class SpringCorotator : public Corotator {
     int n1; 			// node number 1
     int n2; 			// node number 2
     double em;                 // elastic modulus
     double a0;                 // initial cross-sectional area
     double l0;			// initial length
     int fitAlg;                // fit algorithm to use
     FullSquareMatrix & origK;  // Original Stiffness matrix, K
     double u[3], v[3], w[3];   // Original spring element coordinate system vectors
     double Kx, Ky, Kz;         // Spring Stiffness in Each coordinate direction
     enum SpringType {translational, rotational1, rotational2} type;
   public:

     // Constructors -- NOTE: Order of arguments 
     // signifies the type of spring(i.e. translational, corotational)
     // Constructor for translational corotator 
     SpringCorotator(CoordSet &, int *nn, double, double,
                        double, FullSquareMatrix&);

     // Constructor for rotational corotator
     SpringCorotator(int *nn, CoordSet &, int, double, double,
			double,  FullSquareMatrix&);

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt, double t);

     double * getOriginalStiffness() { return origK.data(); }

     void extractDeformations(GeomState &,
                              CoordSet &, double *, int &nlflag);

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                 double *vlr);

     void formTransformationTensor(double x0[2][3], double t[3][3]);

     void formInternalForce(double t[3][3], double p, double *f);

     double getElementEnergy(GeomState &gs, CoordSet &cs); 

};

#endif
