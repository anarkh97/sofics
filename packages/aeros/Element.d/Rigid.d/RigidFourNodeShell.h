#ifndef _RIGIDFOURNODESHELL_H_
#define _RIGIDFOURNODESHELL_H_

#include <Element.d/SuperElement.h>

class GeomState;
class MultiFront;
class NLMaterial;
class ExpMat;

class RigidFourNodeShell : public SuperElement
{
	ExpMat *expmat;
	PressureBCond *pbc;

public:
	explicit RigidFourNodeShell(int*);

	int getElementType() const override { return 76; }
	int getTopNumber() const override { return 188; }
	bool isRigidElement() const override { return true; }
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront *mf) override;

	int getMassType() const override { return 0; }
	FullSquareMatrix massMatrix(const CoordSet&, double* mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;
	void getGravityForce(CoordSet&, double* gravity, Vector&, int gravflg,
						 GeomState *gs) override;

	void setPressure(PressureBCond *_pbc) override { pbc = _pbc; }
	PressureBCond* getPressure() override { return pbc; }
	void computePressureForce(CoordSet&, Vector& elPressureForce,
							  GeomState* gs, int cflg, double t) override;

	void computeDisp(CoordSet&, State&, const InterpPoint&, double*,
					 GeomState*) override;
	void getFlLoad(CoordSet&, const InterpPoint&, double*, double *,
				   GeomState*) override;
};

#endif

