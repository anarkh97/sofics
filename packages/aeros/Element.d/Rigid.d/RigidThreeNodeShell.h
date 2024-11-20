#ifndef _RIGIDTHREENODESHELL_H_
#define _RIGIDTHREENODESHELL_H_

#include <Element.d/SuperElement.h>

class RigidThreeNodeShell : public SuperElement
{
	PressureBCond *pbc;

public:
	explicit RigidThreeNodeShell(int*);
	int getElementType() const override { return 73; }
	int getTopNumber() const override { return 108; }
	bool isRigidElement() const override { return true; }
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront *mf) override;

	int getMassType() const override { return 0; }
	FullSquareMatrix massMatrix(const CoordSet&, double* mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;
	void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
									 GeomState *gs) override;

	void             computeDisp(CoordSet&, State &, const InterpPoint &,
								 double*, GeomState *gs) override;
	void             getFlLoad(CoordSet &, const InterpPoint &,
							   double *flF, double *resF, GeomState *gs) override;

	void setPressure(PressureBCond *_pbc) override { pbc = _pbc; }
	PressureBCond* getPressure() override { return pbc; }
	void computePressureForce(CoordSet&, Vector& elPressureForce,
							  GeomState *gs, int cflg, double t) override;
};

#endif
