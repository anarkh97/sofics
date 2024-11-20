#ifndef _FORCEFUNCTIONELEMENT_H_
#define _FORCEFUNCTIONELEMENT_H_

#include <Element.d/Force.d/BoundaryElement.h>

class DofSet;
class GeomState;

template<template <typename S> class VectorValuedFunctionTemplate>
class ForceFunctionElement : public BoundaryElement
{
public:
    ForceFunctionElement(int, DofSet, int*);
    ForceFunctionElement(int, DofSet*, int*);
    ForceFunctionElement(int, DofSet*, DofSet*, int*);

    FullSquareMatrix stiffness(const CoordSet&, double*, int = 1) const override;
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
    void getStiffAndForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
    void getInternalForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0.0) override;
private:
    void getJacobian(const GeomState *refState, const GeomState &c1, const CoordSet& c0, FullM& B, double t) const;

protected:
    virtual void getConstants(const CoordSet&,
                              Eigen::Array<typename VectorValuedFunctionTemplate<double>::ScalarConstantType,
                                      VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
                              Eigen::Array<int,
                                      VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
                              const GeomState* = nullptr, double = 0) const {}

    virtual void getInputs(Eigen::Matrix<double,VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
                           const CoordSet& c0, const GeomState *c1 = nullptr, const GeomState *refState = nullptr) const;
};

#ifdef _TEMPLATE_FIX_
#include <Element.d/Force.d/ForceFunctionElementImpl.h>
#endif

#endif
