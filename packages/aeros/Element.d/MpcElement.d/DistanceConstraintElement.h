#ifndef _DISTANCECONSTRAINTELEMENT_H_
#define _DISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/DistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DistanceConstraintElement : public ConstraintFunctionElement<Simo::DistanceConstraintFunction>
{
public:
	DistanceConstraintElement(int* _nn, double f0, int type = 0);
	Category getCategory() const override { return Category::Structural; }
	double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
	double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;

protected:
	double f0;
	void getConstants(const CoordSet & cs, Eigen::Array<double,4,1>& sconst, Eigen::Array<int,0,1>&,
					  const GeomState* = nullptr) const override;
};

#endif
