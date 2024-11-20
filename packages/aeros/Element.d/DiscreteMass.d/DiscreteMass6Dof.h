#ifndef _DISCRETEMASS6DOF_H_
#define _DISCRETEMASS6DOF_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>

#include <Eigen/Core>

// Discrete mass+inertia with offset, no gravity force

class DiscreteMass6Dof : public Element, public Corotator
{
	int nn[1];
	Eigen::Matrix3d *C0;
	Eigen::Vector3d f0; // m*g;

public:
	explicit DiscreteMass6Dof(int*);
	~DiscreteMass6Dof() override;

	int getElementType() const override { return 131; }
	Category getCategory() const override { return Category::Structural; }
	void setFrame(EFrame *elemframe) override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	int numNodes() const override { return 1; }
	int *nodes(int *) const override;
	int numDofs() const override { return 6; }
	int* dofs(DofSetArray&, int*) const override;
	void markDofs(DofSetArray&) const override;
	bool hasRot() const override { return true; }
	int getTopNumber() const override { return 506; }

	FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg) const override;
	int getMassType() const override { return 0; } // lumped only
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double  getMass(const CoordSet& cs) const override;
	void getGravityForce(CoordSet&, double *g, Vector &f, int gravflg, GeomState *gs) override;
	void computePressureForce(CoordSet&, Vector& elPressureForce, GeomState *, int, double) override;

	// nonlinear functions
	Corotator* getCorotator(CoordSet&, double*, int, int) override { return this; }
	void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
	void getStiffAndForce(GeomState *refState, GeomState &curState, CoordSet &c0,
						  FullSquareMatrix &elk, double *f, double dt, double t) override;
	void getInternalForce(GeomState *, GeomState &curState, CoordSet &,
						  FullSquareMatrix&, double *f, double, double) override;
	double getElementEnergy(GeomState&, CoordSet&) override;
	bool useDefaultInertialStiffAndForce() override { return false; }
	void getInertialStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0,
								  FullSquareMatrix &elK, double *f, double dt, double t,
								  double beta, double gamma, double alphaf, double alpham) override;
	void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
							 int &nlflag) override { nlflag = 2; }
	void getNLVonMises(Vector&, Vector& weight, GeomState &, CoordSet &, int);
};

#endif
