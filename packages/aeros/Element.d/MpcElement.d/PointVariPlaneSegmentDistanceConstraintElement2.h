#ifndef _POINTVARIPLANESEGMENTDISTANCECONSTRAINTELEMENT2_H_
#define _POINTVARIPLANESEGMENTDISTANCECONSTRAINTELEMENT2_H_

#include <Element.d/Function.d/Constraint.d/PointVariPlaneSegmentDistanceConstraintFunction2.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPlaneSegmentDistanceConstraintElement2 : public ConstraintFunctionElement<Simo::PointVariPlaneSegmentDistanceConstraintFunction2>
{
  public:
    PointVariPlaneSegmentDistanceConstraintElement2(int* _nn); 
  protected:
    void getConstants(const CoordSet & cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst,
                      const GeomState* = nullptr) const override;
};

#endif
