#ifndef _RIGIDBEAM_H_
#define _RIGIDBEAM_H_

#include <Element.d/SuperElement.h>

class RigidBeam : public SuperElement
{
	EFrame *elemframe;
	double c0[3][3];
	int variant;
	double length;
public:
	RigidBeam(int*, int=0);

	int getElementType() const override { return 132; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override { return 106; }
	bool isRigidElement() const override { return true; }
	bool hasRot() const override { return true; }
	bool isSafe() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;

	void buildFrame(CoordSet&) override;

private:
	void getLength(CoordSet&, double &length);
};

class RigidBeam133 : public RigidBeam {
public:
	using RigidBeam::RigidBeam;

	int getElementType() const override { return 133; }
};

class RigidBeamWithMass : public SuperElement
{
	EFrame *elemframe;
	double c0[3][3];
	int variant;
	double length;
public:
	RigidBeamWithMass(int*, int=0);

	int getElementType() const override { return 66; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override { return 106; }
	bool isRigidElement() const override { return true; }
	bool hasRot() const override { return true; }
	bool isSafe() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;

	void buildFrame(CoordSet&) override;
	void setProp(StructProp *p, bool _myProp) override;
	int getMassType() const override { return 0; } // lumped
	double computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs, GeomState *gs,
	                                double stable_tol, int stable_maxit) override;

private:
	void getLength(CoordSet&, double &length);
};

class RigidBeamWithMass106 : public RigidBeamWithMass {
public:
	using RigidBeamWithMass::RigidBeamWithMass;

	int getElementType() const override { return 106; }
};

#endif
