#ifndef _TET_COROTATOR_H_
#define _TET_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class MFTTData;

class TetCorotator : public Corotator {
     int nodeNum[4];
     double em;                 // elastic modulus
     double nu;                 // Poisson's ratio
     double Tref;               // ambient temperature
     double alpha;              // thermal expansion coefficient
     MFTTData *ymtt;
     MFTTData *ctett;
   public:

     // Constructor
     TetCorotator(int nn[4], double, double, CoordSet &, double, double, MFTTData *, MFTTData *);
     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void   getInternalForce(GeomState &gs, CoordSet &cs,
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void   extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                                int &nlflag);

     void   extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                   double *vlr);

     void   getNLVonMises(Vector& stress, Vector& weight, GeomState &curState,
                          GeomState *refState, CoordSet& c0, int strIndex, int surface = 0,
                          double ylayer = 0, double zlayer = 0, int avgnum = 0, int measure = -1);

     void   getNLAllStress(FullM &stress, Vector &weight, GeomState &curState,
                           GeomState *refState, CoordSet &c0, int strInd, int surface = 0,
                           int measure = -1);

     double computeShapeGrad(CoordSet &nodes, double nGrad[4][3]);
 
     void   computeStrainGrad(GeomState &, double nGrad[4][3], double dedU[12][6]);

     void   computePiolaStress(GeomState &, CoordSet &cs,
                               double stress[4][7], double strain[4][7]);

     double getElementEnergy(GeomState &gs, CoordSet &cs);

};

#endif
