#ifndef _MIXEDFINITEELEMENT_H_
#define _MIXEDFINITEELEMENT_H_

#include <Element.d/Force.d/BoundaryElement.h>

class DofSet;
class GeomState;
class NLMaterial;

template<template <typename S> class ScalarValuedFunctionTemplate>
class MixedFiniteElement : public BoundaryElement
{
    enum {
      N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,
      NumberOfNodes = ScalarValuedFunctionTemplate<double>::NumberOfNodes,
      NumberOfDimensions = ScalarValuedFunctionTemplate<double>::NumberOfDimensions
    };
    bool first_time;
    Eigen::Matrix<double,N,1> q_copy;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> CinvBt;
    Eigen::Matrix<double,Eigen::Dynamic,1> Cinvg;
    double epsilon; // penalty parameter for incompressible materials

  protected:
    int materialType;
    NLMaterial *mat;
    int nIV; // number of internal variables (excluding lagrange multipliers for incompressible materials)
    int nLM; // number of lagrange multiplier variables for incompressible materials

  public:
    MixedFiniteElement(int, DofSet, int*);
    MixedFiniteElement(int, DofSet*, int*);

    void setMaterial(NLMaterial *_mat);
    bool isSafe() const override { return true; }

    FullSquareMatrix stiffness(const CoordSet&, double*, int = 1) const override;
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
    void getStiffAndForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
    void getHessian(GeomState *refState, GeomState *c1, CoordSet& c0, Eigen::Matrix<double,
                    ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,
                    ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates>& H, double t);

    int numStates() { return nIV+nLM; }
    void initStates(double *states) { for(int i=0; i<nIV+nLM; ++i) states[i] = 0; }

    void initMultipliers(GeomState& c1) override;
    double getError(GeomState& c1) override;
    void updateMultipliers(GeomState& c1) override;

    virtual int getQuadratureOrder() = 0;

  private:
    void getConstants(const CoordSet &cs, Eigen::Array<typename ScalarValuedFunctionTemplate<double>::ScalarConstantType,
                      ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
                      Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
                      const GeomState *gs = nullptr) const;
    void getInputs(Eigen::Matrix<double,ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
                   const CoordSet& c0, const GeomState *c1 = nullptr) const;
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/Force.d/MixedFiniteElement.C>
#endif

#endif
