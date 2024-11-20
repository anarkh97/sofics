#ifndef _DOTTYPE1CONSTRAINTELEMENT_H_
#define _DOTTYPE1CONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/DotType1ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType1ConstraintElement : public ConstraintFunctionElement<Simo::DotType1ConstraintFunction>
{
  protected:
    double (*C0)[3]; // initial frame (axes stored row-wise)
    int axis1, axis2;
    double d0;

  public:
    DotType1ConstraintElement(int*, int, int, double=0); 
    ~DotType1ConstraintElement();
    void buildFrame(CoordSet&) override;
    void setFrame(EFrame *) override;
    double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
    double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;

  protected:
    void getConstants(const CoordSet &, Eigen::Array<double,7,1>& sconst, Eigen::Array<int,0,1>&,
                      const GeomState *gs = NULL) const override;
};

#endif
