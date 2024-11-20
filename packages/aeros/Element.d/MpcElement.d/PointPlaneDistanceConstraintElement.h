#ifndef _POINTPLANEDISTANCECONSTRAINTELEMENT_H_
#define _POINTPLANEDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointPlaneDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointPlaneDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointPlaneDistanceConstraintFunction>
{
	double x1[3], x2[3], x3[3]; // coordinates of the 3 points defining the plane

public:
	PointPlaneDistanceConstraintElement(int* _nn);
	void setFrame(EFrame *) override;
	FunctionType functionType() override { return LINEAR; }
	int getElementType() const override { return 79; }
protected:
	void getConstants(const CoordSet & cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst,
	                  const GeomState* = nullptr) const override;
};

#endif
