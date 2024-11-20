#if !defined(_NLPENTAHEDRAL_H_) && defined(USE_EIGEN3)
#define _NLPENTAHEDRAL_H_

#include <Element.d/NonLinearity.d/SolidElementTemplate.h>
#include <Element.d/Function.d/Shape.d/Wedge6LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Wedge15LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Wedge26LagrangePolynomial.h>

typedef SolidElementTemplate<Wedge6LagrangePolynomialShapeFunction,6,6> NLPentahedral6;
typedef SolidElementTemplate<Wedge15LagrangePolynomialShapeFunction,15,9> NLPentahedral15;
typedef SolidElementTemplate<Wedge26LagrangePolynomialShapeFunction,26,18> NLPentahedral26;

#endif
