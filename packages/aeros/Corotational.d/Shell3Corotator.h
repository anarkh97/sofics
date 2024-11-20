#ifndef _SHELL3_COROTATOR_H_
#define _SHELL3_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class Node;
class NodeState;

class Shell3Corotator : public Corotator {
protected:
	int n1, n2, n3;
	double origK[18][18];
	int fitAlg;
public:
	Shell3Corotator();
	Shell3Corotator(int, int, int, FullSquareMatrix &, int fitAlgShell);

	void getStiffAndForce(GeomState &, CoordSet &, FullSquareMatrix &, double *, double dt, double t) override;

	void getDExternalForceDu(GeomState &geomState, CoordSet &cs,
							 FullSquareMatrix &elK, double *locF) override;

	void getInternalForce(GeomState &, CoordSet &, FullSquareMatrix &, double *, double, double) override;

	void getExternalForce(GeomState &,CoordSet &, double*) override;

	void formGeometricStiffness(GeomState &, CoordSet &,
								FullSquareMatrix &, double *) override;

	double * getOriginalStiffness() override { return (double *)origK; }

	void extractDefDisp(Node &nd1, Node &nd2, Node &nd3, NodeState &ns1,
						NodeState &ns2, NodeState &ns3,
						double xl0[3][3], double xln[3][3],
						double t0[3][3], double t0n[3][3], double vld[18]);

	void getGlobalDisp(GeomState& , CoordSet&, Vector& ) override;

	void formCorrectGeometricStiffness(double rotvar[3][3][3],
									   double xln[3][3], double pmat[18][18],
									   double gmat[3][18], double f[18],
									   double stiffGeo1[18][18],
									   double stiffGeo2[18][18], double fe[18],
									   double t0n[3][3]);

	void spinAxialAndMoment(double f[], double fnm[][3]);

	void spinAxial(double f[], double fn[][3]);

	void formRotationGradientMatrix(double xdij[3][3],
									double ydij[3][3], double xln[3][3], double gmat[3][18]);

	void gradDefDisp(double xl0[][3], double xln[][3],
					 double pmat[18][18], double gmat[3][18]);

	void localCoord(double x0[3][3], double xn[3][3],
					double t0[3][3], double t0n[3][3], double xl0[3][3],
					double xln[3][3]);

	void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
							 int &nlflag) override;
	void extractDeformationsDisplacementSensitivity(GeomState &geomState, CoordSet &cs, double *dvld) override;

	void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
								double *vlr) override;

	double getElementEnergy(GeomState &, CoordSet &) override;

	void reBuildorigK(FullSquareMatrix &) override;

};

#endif
