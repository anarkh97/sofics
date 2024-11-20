#ifndef _RIGIDTWONODETRUSS_H_
#define _RIGIDTWONODETRUSS_H_

#include <Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.h>

class RigidTwoNodeTruss : public ConstantDistanceConstraint
{
public:
	explicit RigidTwoNodeTruss(int*);
	int getTopNumber() const override { return 101; }
	bool isRigidElement() const override { return true; }
	bool isSafe() const override { return false; }
	PrioInfo examine(int sub, MultiFront*) override;
};

class RigidTwoNodeTrussWithMass : public ConstantDistanceConstraint
{
public:
	explicit RigidTwoNodeTrussWithMass(int*);

	int getElementType() const override { return 65; }
	int getTopNumber() const override { return 101; }
	bool isRigidElement() const override { return true; }
	bool isSafe() const override { return false; }
	PrioInfo examine(int sub, MultiFront*) override;

	int getMassType() const override { return 2; } // both consistent and lumped
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double  getMass(const CoordSet& cs) const override;
	void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg,
						 GeomState *gs) override;
};

#endif
