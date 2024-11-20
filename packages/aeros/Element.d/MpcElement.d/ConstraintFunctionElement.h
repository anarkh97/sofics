#ifndef _CONSTRAINTFUNCTIONELEMENT_H_
#define _CONSTRAINTFUNCTIONELEMENT_H_

#include <Element.d/MpcElement.d/MpcElement.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

class DofSet;
class GeomState;

template<template <typename S> class ConstraintFunctionTemplate>
class ConstraintFunctionElement : public MpcElement
{
  private:
    CoordSet *c0;
  protected:
#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
    template <typename S>
      using ConstraintJacobian = Simo::Jacobian<S,ConstraintFunctionTemplate>;
    template <typename S>
      using ConstraintJacobianVectorProduct = Simo::JacobianVectorProduct<S,ConstraintFunctionTemplate>;
#endif
  public:
    ConstraintFunctionElement(int, DofSet, int*, int);
    ConstraintFunctionElement(int, DofSet*, int*, int);

    void buildFrame(CoordSet&) override;
    void setProp(StructProp *p, bool _myProp) override;
    void update(GeomState*, GeomState&, CoordSet&, double) override;
    void getHessian(const GeomState*, const GeomState&, const CoordSet&, FullSquareMatrix&, double) const override;
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0.0) override;
    double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
    double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&s, double) override;

    FunctionType functionType() override { return NONLINEAR; }

  protected:
    virtual void getConstants(const CoordSet&,
                              Eigen::Array<typename ConstraintFunctionTemplate<double>::ScalarConstantType,
                                           ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
                              Eigen::Array<int,
                                           ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
                              const GeomState* = nullptr) const {}

   virtual void getInputs(Eigen::Matrix<double,ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
                          const CoordSet& c0, const GeomState *c1 = NULL, const GeomState *refState = NULL) const;
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/MpcElement.d/ConstraintFunctionElement.C>
#endif

#endif
