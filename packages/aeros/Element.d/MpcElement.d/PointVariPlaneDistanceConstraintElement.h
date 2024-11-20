#ifndef _POINTVARIPLANEDISTANCECONSTRAINTELEMENT_H_
#define _POINTVARIPLANEDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointVariPlaneDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPlaneDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointVariPlaneDistanceConstraintFunction>
{
public:
	PointVariPlaneDistanceConstraintElement(int* _nn);

	int getElementType() const override { return 179; }
protected:
	void getConstants(const CoordSet & cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst,
	                  const GeomState* = nullptr) const;
};

#endif
