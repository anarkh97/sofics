#ifndef _NON_LIN_QSTATIC_H_
#define _NON_LIN_QSTATIC_H_

#include <Problems.d/NonLinStatic.h>

class NonLinQStatic : public NonLinStatic
{
    Vector &rhs;
    GeomState *geomState;

  public:
    NonLinQStatic(Domain *d, Vector &_rhs, GeomState *_geomState)
     : rhs(_rhs), geomState(_geomState), NonLinStatic(d) {}

    void getRHS(Vector &_rhs) { _rhs = rhs; }
    GeomState* createGeomState() { return new GeomState(*geomState); }
    void staticOutput(GeomState *, double, Vector &, Vector &, GeomState *) {}
};

#endif
