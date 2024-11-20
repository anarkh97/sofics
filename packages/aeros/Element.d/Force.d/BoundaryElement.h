#ifndef _BOUNDARYELEMENT_H_
#define _BOUNDARYELEMENT_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <vector>

class DofSet;

class BoundaryElement : public Element, public Corotator
{
protected:
	int nNodes;               // number of nodes
	int *nn;                  // node numbers
	std::vector<BCond> terms;
	int nterms;

	void addTerms(DofSet);
	void addTerms(DofSet*);
	void addTerms(DofSet*, DofSet*);

	std::vector<int> inputs;
	std::vector<int> outputs;

public:
	BoundaryElement(int, DofSet, int*);
	BoundaryElement(int, DofSet*, int*);
	BoundaryElement(int, DofSet*, DofSet*, int*);
	~BoundaryElement() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	int numNodes() const override;
	int *nodes(int *) const override;

	int numDofs() const override;
	int* dofs(DofSetArray&, int*) const override;
	void markDofs(DofSetArray&) const override;

	Corotator* getCorotator(CoordSet&, double*, int, int) override;
	double getElementEnergy(GeomState&, CoordSet&) override { return 0; }

	bool isSafe() const override { return false; }

	int getTopNumber() const override { return 101; }

	void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
							 int &nlflag) override { nlflag = 2; }

	void getNLVonMises(Vector&, Vector& weight,
					   GeomState &, CoordSet &, int);

	void getGravityForce(CoordSet&, double*, Vector&, int, GeomState*) override;
};
#endif
