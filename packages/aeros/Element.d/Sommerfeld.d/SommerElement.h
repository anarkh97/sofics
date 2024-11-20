#ifndef _SOMMERELEMENT_H_
#define _SOMMERELEMENT_H_

#include <Element.d/Element.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/ComplexD.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/FullMatrix.h>
#include <Math.d/FullRectMatrix.h>
#include <Driver.d/Mpc.h>
#include <gsl/span>

class Domain;

class GeomState;

class SommerElement : public Element {

private:
	SommerElement(const SommerElement &);

protected:
	static bool first; //HB for debugging coupled pb

public:
	Element *el;//Adjacent 3D element
	int iEle;   // and it's index in the CoordSet
	Element *el2;
	Domain *dom;
	complex<double> soundSpeed; // for wet scattering rhs
        bool sFlag;
	SommerElement(Element *_el = 0, Domain *_dom = 0) {
		el = _el;
		dom = _dom;
		el2 = 0;
		iEle = -1;
	}
	Category getCategory() const override { return Category::Acoustic; }
	virtual int getNode(int) const = 0;

	virtual const int *getNodes() const = 0;
	virtual int *getNodes() = 0;

	gsl::span<const int> getNodeSpan() const {
		auto nd = getNodes();
		return { nd, nd+numNodes() };
	}

	int *nodes(int * = 0) const override;

	virtual int numWetDofs() const;

	virtual int numSolidDofs() const;

	int dim() const override;

	virtual int *wetDofs(DofSetArray &, int *p = 0) const;

	virtual int *solidDofs(DofSetArray &, int *p = 0) const;

	void markDofs(DofSetArray &) const override;

	void renum(const int *) override;

	void renum(EleRenumMap &) override;

	SommerElement *clone() override;

	virtual int nFaceCorners() const { return 0; }

	virtual int *faceCorners() const { return nullptr; }

	virtual void flipNormal();

	virtual void checkAndSetNormal(CoordSet &cs);

	virtual void center(CoordSet &cs, double *c);

	virtual void findBothEle(const Connectivity *nodeToElem, int *eleTouch,
	                         int *eleCount, int myNum, const Elemset *eset, int *ie);

	virtual int findAndSetBothEle(const CoordSet &cs, Elemset &eset,
	                              const Connectivity *nodeToElem, int *eleTouch, int *eleCount, int myNum);

	virtual int findEle(const Connectivity *nodeToEle, int *eleTouch,
	                    int *eleCount, int myNum, Elemset *eset = 0,
	                    int it = 0);

	virtual int findAndSetEle(const CoordSet &cs, Elemset &eset,
	                          const Connectivity *nodeToEle, int *eleTouch, int *eleCount, int myNum,
	                          int it = 0);

	virtual FullSquareMatrix sommerMatrix(CoordSet &) const;

	virtual FullSquareMatrix turkelMatrix(CoordSet &) const;

	virtual FullSquareMatrix refinedSommerMatrix(CoordSet &);

	FullSquareMatrix HSommerMatrix(const CoordSet &) const;

	//virtual FullSquareMatrix HKSommerMatrix(CoordSet&);
	virtual FullSquareMatrix interfMatrixConsistent(CoordSet &);

	virtual FullSquareMatrix interfMatrixLumped(CoordSet &);

	virtual FullSquareMatrix sommerMatrix(CoordSet &, double *) const =0;

	virtual GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &,
	                                                        double *);

	virtual void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd);

	virtual void ffp(CoordSet &cs, int numFFP, double *dirFFP,
	                 complex<double> *sol, complex<double> *ffpv, bool direction = true);

	virtual FullRectMatrix transferMatrix(CoordSet &, double *);

	virtual FullSquareMatrix turkelMatrix(CoordSet &, double *) const;

	virtual FullSquareMatrix refinedSommerMatrix(CoordSet &, double *);

	//virtual FullSquareMatrix surfStiffMatrix(CoordSet&, double *);
	virtual FullSquareMatrix HSommerMatrix(const CoordSet &cs, double *d) const;

	//virtual FullSquareMatrix HKSommerMatrix(CoordSet&, double *);
	virtual FullSquareMatrix interfMatrixConsistent(CoordSet &, double *);

	virtual FullSquareMatrix interfMatrixLumped(CoordSet &, double *);

	virtual void BT2(CoordSet &cs, double *e, double *f, double *g,
	                 double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d);

	virtual void BT2n(CoordSet &cs, double *e, double *f, double *g,
	                  double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d, int n);

	virtual void sphereBT2(CoordSet &cs, double r, double k, ComplexD *d);

	virtual void ellipsoidBT2(CoordSet &cs, double a, double b, double k, ComplexD *d);

	virtual void sommerMatrixEllipsoid(CoordSet &cs, double kappa, double H[3], double K[3], ComplexD *d);

	virtual void neumVector(CoordSet &, Vector &, int pflag = 0, GeomState * = 0, double t = 0);

	virtual void neumVectorJacobian(CoordSet &, FullSquareMatrix &, int pflag = 0, GeomState * = 0, double t = 0);

	virtual void neumVector(CoordSet &, ComplexVector &,
	                        double, double, double, double, int pflag = 0);

// RT: obsolete? and I need it for something else
//        virtual void neumVector(CoordSet&,ComplexVector&,
//                                double,double,double,double,int);
	virtual void neumVectorDeriv(CoordSet &, ComplexVector &,
	                             double, double, double, double, int,
	                             int pflag = 0);

	virtual void wetInterfaceVector(CoordSet &, ComplexVector &,
	                                double, double, double, double, int, int);

	virtual void wetInterfaceVector(CoordSet &, ComplexVector &,
	                                complex<double> (*)[3], complex<double> *);

	virtual void wetInterfaceVectorDeriv(CoordSet &, ComplexVector &,
	                                     complex<double> (*)[3], complex<double> *,
	                                     complex<double> *, int);

	virtual void sommerVector(CoordSet &, ComplexVector &, ComplexVector &);

	virtual void btVector(CoordSet &, ComplexVector &, ComplexVector &);

	virtual void ffpNeum(int, ComplexD *, CoordSet &, ComplexD *,
	                     double, double(*)[3], double *);

	virtual void ffpDir(int, ComplexD *, CoordSet &, ComplexD *, ComplexD *,
	                    double, double(*)[3], double *);

	virtual void ffpDir(int, ComplexD *, CoordSet &, ComplexD *,
	                    double, double(*)[3], double *);

	virtual ComplexD ffpCoef(double k) const;

	virtual void getNormal(const CoordSet &, double[3]) const;

	virtual double getSize(CoordSet &);

	virtual FullSquareMatrixC turkelMatrix(CoordSet &, double, int);

	virtual FullSquareMatrixC turkelMatrix(CoordSet &, double, int, DComplex *);

	bool isSommerElement() const override { return true; }

	bool isPhantomElement() override { return false; }

	FullSquareMatrix stiffness(const CoordSet &, double *d, int flg = 1) const override;

	FullSquareMatrix massMatrix(const CoordSet &, double *mel, int cmflg = 1) const override;

	FullSquareMatrix dampingMatrix(CoordSet &, double *cel, int cmflg = 1) override;
	//virtual FullSquareMatrix imStiffness(CoordSet&, double *d, int flg = 1);

	// TODO Get rid of this here!!! Pure math functions on arguments have nothing to do here.
	int ludcmp(FullSquareMatrix A, int n, int *indx) const;// LU decomposition
	void lubksb(FullSquareMatrix A, int n, int *indx, double *b) const;//LU factorisation
	void invert(FullSquareMatrix A, FullSquareMatrix B) const;

	int getAdjElementIndex() { return iEle; }

	void setAdjElementIndex(int _iEle) { iEle = _iEle; }
};

#endif

