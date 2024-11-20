#ifndef _BEAM_COROTATOR_H_
#define _BEAM_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class Node;
class NodeState;

class BeamCorotator : public Corotator {
   int n1;			// node number 1
   int n2;			// node number 2
   double* zVec;
   double origK[12][12];	// element stiffness matrix
   int fitAlg;			// fit algorithm
 public:
   // Constructor
   BeamCorotator(int node1, int node2, double [3], FullSquareMatrix &kel, 
                 int fitAlgBeam);

   double * getOriginalStiffness() override { return (double*)origK; }

   void localCoord(double zVec[3], double zVecL[2][3], 
                double (* rot[2])[3][3],
                double x0[2][3], double xn[2][3], double t0[3][3], 
                double t0n[3][3],
                double xl0[2][3], double xln[2][3]);

   void getStiffAndForce(GeomState &, CoordSet &, FullSquareMatrix &, double *, double, double) override;
   
   void formCorrectGeometricStiffness(
                            double pmat[12][12],double gmat[3][12], 
			    double f[12], double stiffGeo1[12][12], 
                            double stiffGeo2[12][12],double fe[12],
			    Node &node1, Node &node2, NodeState &ns1,
			    NodeState &ns2);
   
   void getInternalForce(GeomState &geomState, CoordSet &cs, 
                         FullSquareMatrix &elK, double *f, double dt, double t) override;

   void formGeometricStiffness(GeomState &, CoordSet &, 
                               FullSquareMatrix &, double *) override;

   void extractDefDisp(Node &nd1, Node &nd2, NodeState &ns1, NodeState &ns2,
                       double zVecL[2][3], double xl0[2][3], double xln[2][3],
                       double t0[3][3], double t0n[3][3], double vld[12]);

   void corotStiffGeo(double zVecL[2][3],
                      double xln[2][3], double pmat[12][12], double f[12],
                      double stiffGeo1[12][12], double stiffGeo2[12][12]);

   void spinAxialAndMoment (double f[], double fnm[][3]);

   void spinAxial(double f[], double fn[][3]);
    
   void gradLocRot(double len, double zVecL[2][3], 
                   double gmat[3][12]);

   void gradDefDisp(double zVecL[][3],double xln[][3],
                    double pmat[12][12], double gmat[3][12]);

   void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                            int &nlflag) override;

   void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                               double *vlr) override;

   double getElementEnergy(GeomState &, CoordSet &) override;
 
   void reBuildorigK(FullSquareMatrix &) override;
   
   void localCoord(double zVec[3], double dzVec[3], double zVecL[2][3],
                  double dzVecL[2][3], double (* rot[2])[3][3], 
                  double (* drot[2])[3][3], double x0[2][3],
                  double dx0[2][3], double xn[2][3], double dxn[2][3],
                  double t0[3][3], double dt0[3][3],double t0n[3][3], 
                  double dt0n[3][3], double xl0[2][3], double dxl0[2][3],
                  double xln[2][3], double dxln[2][3]);

 void  dgradLocRot(double len, double dlen, double zVecL[2][3],
                   double dzVecL[2][3], double gmat[3][12],
                   double dgmat[3][12]);
};


#endif
