#ifndef _NLMEMBRANE_H_
#define _NLMEMBRANE_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/NonLinearity.d/ShapeFunction.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Math.d/TTensor.h>

template <int n>
struct TwoDTensorTypes
{
	static const int ndofs=n;
	typedef Stress2D StressTensor;
	typedef Stress2D StrainTensor;
	typedef SimpleTensor<StrainTensor, n> BTensor;
	typedef SimpleTensor<SimpleTensor<double,3>,2> GradUTensor;
	typedef SimpleTensor<GradUTensor, n> GradUDerivTensor;
	typedef SimpleTensor<BTensor,n> DBTensor;
	// Material tangent matrix tensor type
	typedef SymTensor<StressTensor,2> DTensor;
};

// Double contraction operator
template <int n>
SimpleTensor<Stress2D, n>
operator||(SimpleTensor<SimpleTensor<Stress2D,n>,n> &, SimpleTensor<Stress2D, n> &);

template <int n>
class LinearStrain2D : public GenStrainEvaluator<TwoDTensorTypes<n> >
{
	void getE(typename TwoDTensorTypes<n>::StrainTensor &e,
			  typename TwoDTensorTypes<n>::GradUTensor &gradU);
	void getEandB(typename TwoDTensorTypes<n>::StrainTensor &e,
				  typename TwoDTensorTypes<n>::BTensor &B,
				  typename TwoDTensorTypes<n>::GradUTensor &gradU,
				  typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
	void getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e,
					typename TwoDTensorTypes<n>::BTensor &B,
					typename TwoDTensorTypes<n>::DBTensor &DB,
					typename TwoDTensorTypes<n>::GradUTensor &gradU,
					typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
};

template <int n>
class GLStrain2D : public GenStrainEvaluator<TwoDTensorTypes<n> >
{
public:
	void getE(typename TwoDTensorTypes<n>::StrainTensor &e,
			  typename TwoDTensorTypes<n>::GradUTensor &gradU);
	void getEandB(typename TwoDTensorTypes<n>::StrainTensor &e,
				  typename TwoDTensorTypes<n>::BTensor &B,
				  typename TwoDTensorTypes<n>::GradUTensor &gradU,
				  typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
	void getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e,
					typename TwoDTensorTypes<n>::BTensor &B,
					typename TwoDTensorTypes<n>::DBTensor &DB,
					typename TwoDTensorTypes<n>::GradUTensor &gradU,
					typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
	bool isNonLinear() { return true; }
};

class TriMembraneShapeFunct : public GenShapeFunction< TwoDTensorTypes<9> >
{
public:
	TriMembraneShapeFunct() : GenShapeFunction< TwoDTensorTypes<9> >() {}
	void getGlobalGrads(Grad2D &gradU, Grad2DDeriv9 &dGradUdqk, double *_jac,
						Node *nodes, double xi[3], Vector &disp) override;
	void getGradU(Grad2D &gradU,
				  Node *nodes, double xi[3], Vector &disp) override;
	double interpolateScalar(double *_q, double _xi[3]) override;
};

class NLMembrane : public GenGaussIntgElement<TwoDTensorTypes<9> >
{
	int n[3];
    NLMaterial *material, *linearMaterial;
	double *cCoefs;
	double *cFrame;
	bool useDefaultMaterial;
	PressureBCond *pbc;

protected:
	int getNumGaussPoints() const override;
	void getGaussPointAndWeight(int i, double *point, double &weight) const override;
    void getLocalNodalCoords(int, double *) override;
	GenShapeFunction< TwoDTensorTypes<9> > *getShapeFunction() const override;
	GenStrainEvaluator<TwoDTensorTypes<9> > *getGenStrainEvaluator() const override;
	const NLMaterial *getMaterial() const override;
    NLMaterial *getLinearMaterial() const override;
    void rotateCFrame(const CoordSet &cs, double *T, double *Tinv) const;

public:
	explicit NLMembrane(int *nd);
	~NLMembrane() override;

	int getElementType() const override { return 129; }
	void setPreLoad(std::vector<double> &_preload) override { preload = _preload; }
	std::vector<double> getPreLoad() override { return preload; }
	int numNodes() const override { return 3; }
	int numDofs() const override { return 9; }
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int* nodes(int *) const override;
	void updateStates(Node *nodes, double *states, double *un, double *unp) {}
	void setProp(StructProp *p, bool _myProp) override;
    void buildFrame(CoordSet &cs) override;
	void setCompositeData(int, int, double *, double *coefs, double *frame) override;
	double* setCompositeData2(int, int, double *, double *coefs, CoordSet &cs, double theta) override;
	void setMaterial(NLMaterial *) override;
	void setPressure(PressureBCond *_pbc) override { pbc = _pbc; }
	PressureBCond* getPressure() override { return pbc; }
	void computePressureForce(CoordSet &cs, Vector& elPressureForce,
							  GeomState *gs, int cflg, double t) override;

	Corotator* getCorotator(CoordSet &, double *, int , int) override;
	int getTopNumber() const override { return 104; }
	FullSquareMatrix  stiffness(const CoordSet& cs, double *k, int flg) const override;
#ifdef USE_EIGEN3
	int getMassType() const override { return 2; } // both consistent and lumped
#else
	int getMassType() const override { return 0; } // lumped only
#endif
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg) const override;
	double getMass(const CoordSet&) const override;

	void computeDisp(CoordSet &cs, State &state, const InterpPoint &, double *res, GeomState *gs) override;
	void getFlLoad(CoordSet &, const InterpPoint &,  double *flF, double *resF, GeomState *gs) override;

	PrioInfo examine(int sub, MultiFront *mf) override;
};

#include <Element.d/SuperElement.h>

class NLMembrane4 : public SuperElement
{
public:
	explicit NLMembrane4(int *nodenums);

	int getElementType() const override { return 128; }
	int getTopNumber() const override;
	void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs) override;
	void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs) override;
	PrioInfo examine(int sub, MultiFront *mf) override;
};

#endif
