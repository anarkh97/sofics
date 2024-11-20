#ifndef _POINTLINEDISTANCECONSTRAINTELEMENT_H_
#define _POINTLINEDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointLineDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointLineDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointLineDistanceConstraintFunction>
{
	double x1[3], x2[3]; // coordinates of the 2 points defining the line

public:
	PointLineDistanceConstraintElement(int* _nn);
	void setFrame(EFrame *) override;
	int getElementType() const override { return 78; }
protected:
	void getConstants(const CoordSet & cs, Eigen::Array<double,14,1>& sconst, Eigen::Array<int,1,1>& iconst,
					  const GeomState* = nullptr) const override;
};

#endif
