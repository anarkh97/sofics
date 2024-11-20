#if !defined(_NLHEXAHEDRAL_H_) && defined(USE_EIGEN3)
#define _NLHEXAHEDRAL_H_

#include <Element.d/NonLinearity.d/SolidElementTemplate.h>
#include <Element.d/Function.d/Shape.d/Hex8LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Hex20LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Hex32LagrangePolynomial.h>

class NLHexahedral8 : public SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>
{
public:
	using SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>::SolidElementTemplate;
	int getElementType() const override { return 8808; }
};
class NLHexahedral20 : public SolidElementTemplate<Hex20LagrangePolynomialShapeFunction,20,27>
{
public:
	using SolidElementTemplate<Hex20LagrangePolynomialShapeFunction,20,27>::SolidElementTemplate;
	int getElementType() const override { return 8820; }
};
class NLHexahedral32 : public SolidElementTemplate<Hex32LagrangePolynomialShapeFunction,32,64>
{
public:
	using SolidElementTemplate<Hex32LagrangePolynomialShapeFunction,32,64>::SolidElementTemplate;
	int getElementType() const override { return 8832; }
};

// Nonlinear 8-node brick element with 2x2x2 gauss points
class NLHexahedral : public SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>
{
    bool linearKinematics;
  public:
    NLHexahedral(int *nd, bool isLinKin) : linearKinematics(isLinKin), SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>(nd) {}
    StrainEvaluator* getStrainEvaluator() const override;
    int getTopNumber() const override { return 117; }
    int getElementType() const override { return 1117; } // TODO Update???
};

#endif
