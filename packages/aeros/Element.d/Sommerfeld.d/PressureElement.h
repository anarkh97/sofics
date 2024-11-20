#ifndef _PRESSUREELEMENT_H_
#define _PRESSUREELEMENT_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class DofSet;

class GeomState;

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
class PressureElement : public SommerElement {
protected:
	int nNodes;               // number of nodes
	int *nn;                  // node numbers
	std::vector<BCond> terms;
	int nterms;
	PressureBCond *pbc;

	void addTerms(DofSet);

public:
	PressureElement(int *, PressureBCond *);

	~PressureElement();

	int getElementType() const override { return eType; }
	void renum(const int *) override;

	void renum(EleRenumMap &) override;

	int numNodes() const override;

	int *nodes(int * = 0) const override;

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override;

	int *dofs(DofSetArray &, int *) const override;

	void markDofs(DofSetArray &) const override;

	int findAndSetEle(const CoordSet &cs, Elemset &eset, const Connectivity *nodeToEle, int *eleTouch, int *eleCount,
	                  int myNum,
	                  int it) override;

	PressureBCond *getPressure() override { return pbc; }

	void neumVector(CoordSet &, Vector &, int pflag = 0, GeomState * = 0, double t = 0) override;

	void neumVectorJacobian(CoordSet &, FullSquareMatrix &, int pflag = 0, GeomState * = 0, double t = 0) override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;
};

#ifdef _TEMPLATE_FIX_

#include <Element.d/Sommerfeld.d/PressureElement.C>

#endif

#endif
