#ifndef _NONLINEARTORSIONALSPRING_H_
#define _NONLINEARTORSIONALSPRING_H_

#include <Element.d/MpcElement.d/AngleType1ConstraintElement.h>

// this element is a torsional spring for large displacements and rotations, implemented using
// the penalized constraint method, in which the penalty parameter is the value of the spring stiffness.
// The constraint which is penalized in the relative rotation of node A w.r.t node B about a specified axis.
// i.e f(x) = arccos(c1.dot(c2))
// To handle arbitriarily large rotations and also avoid singularities dynamic reparameterization and 
// history variables are employed

class NonlinearTorsionalSpring : public AngleType1ConstraintElement
{
	int m_axis1, m_axis2;
	int propIndex; // 0: use StructProp::k1, 1: use StructProp::k2, 2: use StructProp::k3
	double offset2;
	int quadrant;

public:
	NonlinearTorsionalSpring(int*, int, int, int=0, int=0, int=1);

	int getElementType() const override { return this->type == 0 ? 202 : 205; }
	Category getCategory() const override { return Category::Structural; }
	void setProp(StructProp *p, bool _myProp) override;
	void update(GeomState *refState, GeomState& gState, CoordSet& cs, double t) override;

	int numStates() override;
	void initStates(double *) override;
	void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0, double dt) override;

	bool isSpring() const override { return true; }
	bool hasRot() const override { return true; }
};

#endif
