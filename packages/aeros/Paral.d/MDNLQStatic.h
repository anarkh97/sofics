#ifndef _MD_NL_QSTATIC_H_
#define _MD_NL_QSTATIC_H_

#include <Paral.d/MDNLStatic.h>

class MDNLQStatic : public MDNLStatic
{
    DistrVector &rhs;
    DistrGeomState *geomState;

  public:
    MDNLQStatic(Domain *d, DecDomain *dd, DistrVector &_rhs, DistrGeomState *_geomState)
     : rhs(_rhs), geomState(_geomState), MDNLStatic(d,dd) {}

    void getRHS(DistrVector &_rhs) { _rhs = rhs; }
    DistrGeomState* createGeomState() { return new DistrGeomState(*geomState); }
    void staticOutput(DistrGeomState *, double, DistrVector &, DistrVector &, DistrGeomState *) {}
};

#endif

