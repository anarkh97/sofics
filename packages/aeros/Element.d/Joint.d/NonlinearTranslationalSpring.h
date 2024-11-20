#ifndef _NONLINEARTRANSLATIONALSPRING_H_
#define _NONLINEARTRANSLATIONALSPRING_H_

#include <Element.d/MpcElement.d/DotType2ConstraintElement.h>

// this element is a translational spring for large displacements and rotations, implemented using
// the penalized constraint method, in which the penalty parameter is the value of the spring stiffness.
// The constraint which is penalized in the scalar projection of vector d (from node A to node B)
// onto the unit vector c1 (either the x,y or z axis of the local frame attached to node A) minus 
// the constant value of same scalar projection in the undeformed configuration.
// i.e f(x) = d.dot(c1) - sp0

class NonlinearTranslationalSpring : public DotType2ConstraintElement
{
	double sp0; // scalar projection of d onto c0 in the undeformed configuration
	int propIndex; // 0: use StructProp::k1, 1: use StructProp::k2, 2: use StructProp::k3

public:
	NonlinearTranslationalSpring(int*, int, int=0, int=0, int=1);

	int getElementType() const override { return 201; }
	Category getCategory() const override { return Category::Structural; }
	void setProp(StructProp *p, bool _myProp) override;
	void buildFrame(CoordSet&) override;
	void update(GeomState *refState, GeomState& gState, CoordSet& cs, double) override;

	bool isSpring() const override { return true; }
	bool hasRot() const override { return true; }
};

class NonlinearTranslationalSpring204 : public NonlinearTranslationalSpring {
public:
	using NonlinearTranslationalSpring::NonlinearTranslationalSpring;

	int getElementType() const override { return 204; }
};
#endif
