#ifndef _POINTVARIPLANESEGMENTDISTANCECONSTRAINTELEMENT_H_
#define _POINTVARIPLANESEGMENTDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointVariPlaneSegmentDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPlaneSegmentDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointVariPlaneSegmentDistanceConstraintFunction>
{
public:
	PointVariPlaneSegmentDistanceConstraintElement(int* _nn);

	int getElementType() const override { return 279; }
protected:
	void getConstants(const CoordSet & cs,
	                  Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst,
	                  const GeomState* = nullptr) const override;
};

class PointVariPlaneSegmentDistanceConstraintElement379 : public PointVariPlaneSegmentDistanceConstraintElement {
public:
	using PointVariPlaneSegmentDistanceConstraintElement::PointVariPlaneSegmentDistanceConstraintElement;

	int getElementType() const override { return 379; }
};
#endif
