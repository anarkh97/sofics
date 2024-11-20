#ifndef _SEGVARISEGDISTANCECONSTRAINTELEMENT_H_
#define _SEGVARISEGDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/SegVariSegDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class SegVariSegDistanceConstraintElement : public ConstraintFunctionElement<Simo::SegVariSegDistanceConstraintFunction>
{

  public:
    SegVariSegDistanceConstraintElement(int* _nn);

	int getElementType() const override { return 173; }
  protected:
    void getConstants(const CoordSet & cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst,
                      const GeomState* = nullptr) const override;
};

#endif
