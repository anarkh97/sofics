#ifndef _POTENTIALFUNCTIONELEMENT_H_
#define _POTENTIALFUNCTIONELEMENT_H_

#include <Element.d/Force.d/BoundaryElement.h>

class DofSet;
class GeomState;

template<template <typename S> class ScalarValuedFunctionTemplate>
class PotentialFunctionElement : public BoundaryElement
{
public:
	PotentialFunctionElement(int, DofSet, int*);
	PotentialFunctionElement(int, DofSet*, int*);

	FullSquareMatrix stiffness(const CoordSet&, double*, int) const override;
	void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
	void getStiffAndForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
	void getInternalForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
	void computePressureForce(CoordSet&, Vector& elPressureForce,
							  GeomState *gs, int cflg, double t) override;
private:
	void getHessian(const GeomState *refState, const GeomState &c1, const CoordSet &c0,
					FullSquareMatrix &B, double t) const;

protected:
	virtual void getConstants(const CoordSet&,
							  Eigen::Array<typename ScalarValuedFunctionTemplate<double>::ScalarConstantType,
									  ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
							  Eigen::Array<int,
									  ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
							  const GeomState* = nullptr, double = 0) const {}

	virtual void getInputs(Eigen::Matrix<double,ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
						   const CoordSet& c0, const GeomState *c1 = nullptr, const GeomState *refState = nullptr) const;
};

#endif
