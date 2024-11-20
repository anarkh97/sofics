#ifndef _RIGIDEIGHTNODEBRICK_H_
#define _RIGIDEIGHTNODEBRICK_H_

#include <Element.d/SuperElement.h>

class RigidEightNodeBrick : public SuperElement
{
public:
	explicit RigidEightNodeBrick(int*);
	int getElementType() const override { return 70; }
	int getTopNumber() const override { return 117; }
	bool isRigidElement() const override { return true; }
	bool isSafe() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;

	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const override;
	double getMass(const CoordSet& cs) const override;
	void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
						 GeomState *gs) override;
};

#endif

