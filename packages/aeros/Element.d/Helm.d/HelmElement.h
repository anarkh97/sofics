#ifndef _HELMELEM_H_
#define _HELMELEM_H_

#include <Element.d/Element.h>
#include <Math.d/ComplexD.h>
#include <Math.d/Vector.h>
#include <Driver.d/PolygonSet.h>


class HelmElement {
public:
	virtual void getHelmForce(CoordSet&, ComplexVector &vc, ComplexVector &force);
	virtual void addFaces(PolygonSet *pset);
	virtual void edgeShapeFunctions(int n1, int n2, int *ng,
									double **gw, double **N);
	virtual FullSquareMatrix acousticm(CoordSet&, double *kel);
	virtual void wErrors(CoordSet&,
						 double *l2e, double *h1e, double *l2, double *h1,
						 ComplexD *u, double kappa, double *waveDir);
	virtual int isFluid() { return 1; }
	virtual void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*);
	virtual void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
								double kappa, double *waveDir);
};

#endif
