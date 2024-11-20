#ifndef _QUADTHERMAL_COROTATOR_H_
#define _QUADTHERMAL_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class QuadThermalCorotator : public Corotator {
     int n1; 			// node number 1
     int n2; 			// node number 2
     int n3;                    // node number 3
     int n4;                    // node number 4
     double eps;                // emissivity of the body
     double sigma;              // Stefan's constant 
     double Tr;                 // Temperature of the enclosure receiving the radiation

   public:

     // Constructor
     QuadThermalCorotator(int node1, int node2, int node3, int node4, double epsilon, double sigma, double Tr, CoordSet &cs);

     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &ts, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void   getInternalForce(GeomState &ts, CoordSet &cs,
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void   formInternalForce(double xl[4], double yl[4], double xn[4], double eps, double sigma, double Tr, double f[4]);  

     void   formTangentStiffness(double xl[4], double yl[4], double xn[4], double eps, double sigma, double kt[4][4]);

};

#endif
