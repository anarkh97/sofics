#ifndef _SOLIDELEMENTTEMPLATE_H_
#define _SOLIDELEMENTTEMPLATE_H_

#include <Element.d/NonLinearity.d/ShapeFunction.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>

class NLMaterial;
class ShapeFunction;
class StrainEvaluator;

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes>
class AutoShapeFunction : public ShapeFunction
{
public:
	AutoShapeFunction() : ShapeFunction(3*NumberOfNodes) {}
	void getLocalDerivatives(Tensor *localDerivatives, double xi[3]) override;
	void getValues(Tensor *val, double xi[3]) override {}
	Tensor *getValInstance() override { return 0; }
	double interpolateScalar(double *_q, double _xi[3]) override;
};

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
class SolidElementTemplate : public GaussIntgElement
{
private:
	static const double nodeRefCoords[NumberOfNodes][3];
	static const AutoShapeFunction<ShapeFunctionTemplate,NumberOfNodes> shapeFunction;

protected:
	int n[NumberOfNodes];
	NLMaterial *material;

	int getNumGaussPoints() const override;
	void getGaussPointAndWeight(int i, double *point, double &weight) const override;
	void getLocalNodalCoords(int i, double *coords) override;
	ShapeFunction *getShapeFunction() const override;
	StrainEvaluator *getStrainEvaluator() const override;
	NLMaterial *getMaterial() const override;
	void getNodeRefCoords(double (*nodeRefCoords)[3]);

public:
	enum { numGaussPoints = NumIntgPts };
	SolidElementTemplate(int *nd);

	int getElementType() const { return -1; }
	int numNodes() const override;
	int numDofs() const override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int* nodes(int * = 0) const override;
	void setMaterial(NLMaterial *) override;
};

#ifdef _TEMPLATE_FIX_
#include <Element.d/NonLinearity.d/SolidElementTemplate.C>
#endif

#endif
