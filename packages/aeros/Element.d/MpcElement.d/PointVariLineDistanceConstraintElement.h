#ifndef _POINTVARILINEDISTANCECONSTRAINTELEMENT_H_
#define _POINTVARILINEDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointVariLineDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariLineDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointVariLineDistanceConstraintFunction>
{
  public:
    PointVariLineDistanceConstraintElement(int* _nn);

	int getElementType() const override { return 178; }
  protected:
    void getConstants(const CoordSet & cs, Eigen::Array<double,14,1>& sconst, Eigen::Array<int,1,1>& iconst,
                      const GeomState* = nullptr) const override;
};

#endif
