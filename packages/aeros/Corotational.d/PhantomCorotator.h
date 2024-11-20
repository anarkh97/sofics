#ifndef _PHANTOMCOROTATOR_H_
#define _PHANTOMCOROTATOR_H_

#include <Corotational.d/Corotator.h>

class PhantomCorotator : public Corotator {
  public:
    void getStiffAndForce(GeomState &, CoordSet &, 
                          FullSquareMatrix &, double *, double, double) override;

    void formGeometricStiffness(GeomState &, CoordSet &, 
                                FullSquareMatrix &, double *) override;

    double* getOriginalStiffness() override;

    void extractDeformations(GeomState &geomState, CoordSet &cs, 
                             double *vld, int &nlflag) override;
    void extractDeformations(GeomState &geomState, CoordSet &cs,
                             DComplex *vld, int &nlflag) override;

    double getElementEnergy(GeomState &, CoordSet &) override;
 
    void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                double *vlr) override;
};

#endif
