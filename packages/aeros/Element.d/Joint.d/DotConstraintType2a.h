#ifndef _DOTCONSTRAINTTYPE2A_H_
#define _DOTCONSTRAINTTYPE2A_H_

#include <Element.d/MpcElement.d/DotType2ConstraintElement.h>

class DotConstraintType2a : public DotType2ConstraintElement
{
  public:
    DotConstraintType2a(int*, int);
    void buildFrame(CoordSet& cs) override;
    void update(GeomState*, GeomState&, CoordSet&, double) override;
    double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
    double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0);
};

#endif
