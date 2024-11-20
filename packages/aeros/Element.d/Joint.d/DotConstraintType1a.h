#ifndef _DOTCONSTRAINTTYPE1A_H_
#define _DOTCONSTRAINTTYPE1A_H_

#include <Element.d/MpcElement.d/DotType1ConstraintElement.h>

class DotConstraintType1a : public DotType1ConstraintElement
{
    int axis1_copy;
    double t_reparam, offset;

  public:
    DotConstraintType1a(int*, int, int);
    void buildFrame(CoordSet& cs) override;
    void update(GeomState*, GeomState&, CoordSet&, double) override;
    double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
    double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0);
};

#endif
