#ifndef _SEGMENTSEGMENTDISTANCECONSTRAINTELEMENT_H_
#define _SEGMENTSEGMENTDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/SegmentSegmentDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class SegmentSegmentDistanceConstraintElement : public ConstraintFunctionElement<Simo::SegmentSegmentDistanceConstraintFunction>
{
	double x1[3], x2[3]; // coordinates of the 2 points defining the line

public:
	SegmentSegmentDistanceConstraintElement(int* _nn);
	void setFrame(EFrame *) override;
	int getElementType() const override { return 174; }
protected:
	void getConstants(const CoordSet & cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst,
					  const GeomState* = nullptr) const override;
};

#endif
