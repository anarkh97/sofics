#ifndef _MPCELEMENT_H_
#define _MPCELEMENT_H_

#include <Driver.d/Mpc.h>
#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <map>

class DofSet;

class MpcElement : public Element, public Corotator, public LMPCons
{
protected:
	int nNodes;              // number of nodes (not including "internal node")
	int *nn;                 // node numbers
	// for direct elimination / coordinate split / state space  set prop->lagrangeMult to false and prop->penalty to 0
	// for lagrange multipliers method set prop->lagrangeMult to true and prop->penalty to 0
	// for penalty method set prop->lagrangeMult to false and prop->penalty to some large number
	// for augmented lagrangian method set prop->lagrangeMult to true and prop->penalty to some large number
	std::map<int,std::vector<int> > rotation_indices;
	std::map<int,std::vector<double> > rotation_coefs;

	void addTerms(DofSet);
	void addTerms(DofSet*);

public:
	MpcElement(int, DofSet, int*);
	MpcElement(int, DofSet*, int*);
	MpcElement(LMPCons *mpc, bool nlflag);
	~MpcElement() override;

	int getNumMPCs() override;
	LMPCons** getMPCs() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	int numNodes() const override;
	int* nodes(int*) const override;

	int numInternalNodes() const override;
	void setInternalNodes(int*) override;

	int numDofs() const override;
	int* dofs(DofSetArray&, int*) const override;
	void markDofs(DofSetArray &dsa) const override;

	bool hasRot() const override;

	FullSquareMatrix stiffness(const CoordSet&, double*, int = 1) const override;

	void getGravityForce(CoordSet&, double*, Vector& f, int, GeomState*) override;

	bool isMpcElement() override { return true; }

	Corotator* getCorotator(CoordSet&, double*, int, int) override;
	void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
	void getStiffAndForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
	void getInternalForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
	void getResidualCorrection(GeomState& c1, double* r) override;
	double getElementEnergy(GeomState&, CoordSet&) override;

	virtual void update(GeomState*, GeomState&, CoordSet&, double);
	virtual void getHessian(const GeomState*, const GeomState&, const CoordSet&, FullSquareMatrix&, double) const;
	virtual double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double);
	virtual double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double);

	PrioInfo examine(int sub, MultiFront *mf) override;
	int getElementType() const override;
	bool isSafe() const override { return false; }

	void computePressureForce(CoordSet&, Vector& elPressureForce,
	                          GeomState *gs, int cflg, double t) override;
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override { return 101; }

	void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
	                         int &nlflag) override { nlflag = 2; }

	void getNLVonMises(Vector&, Vector&, GeomState&, CoordSet&, int);
	void getNLAllStress(FullM&, Vector&, GeomState&, CoordSet&, int);

	void initMultipliers(GeomState& c1) override;
	void updateMultipliers(GeomState& c1) override;

	double getError(GeomState& c1) override;

	enum FunctionType { LINEAR=0, QUADRATIC, NONLINEAR };
	virtual FunctionType functionType() { return NONLINEAR; }

	bool isFreeplayElement() const override { return type == 1 && prop->penalty != 0; }
};
#endif
