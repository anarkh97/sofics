#if !defined(_NLTETRAHEDRAL_H_) && defined(USE_EIGEN3)
#define _NLTETRAHEDRAL_H_

#include <Element.d/NonLinearity.d/SolidElementTemplate.h>
#include <Element.d/Function.d/Shape.d/Tet4LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Tet10LagrangePolynomial.h>

typedef SolidElementTemplate<Tet4LagrangePolynomialShapeFunction,4,1> NLTetrahedral4;
typedef SolidElementTemplate<Tet10LagrangePolynomialShapeFunction,10,15> NLTetrahedral10;

#endif
