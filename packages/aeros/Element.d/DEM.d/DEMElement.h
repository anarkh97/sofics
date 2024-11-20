#ifndef _DEMELEMENT_H_
#define _DEMELEMENT_H_

#include <cmath>

#ifndef SANDIA

#include <Element.d/Element.h>
#include <Math.d/ComplexD.h>

#else
#include <complex>
using std::complex;
#endif

class DEMElement;

class DEMLM {
public:
	DEMElement *e1,*e2;
	int nOffset;
	virtual int type()const =0;
	virtual int nDofs()const=0;
	virtual void init()=0;
	int setNodeOffset(int j) { nOffset = j; return(0); }
	int getNodeOffset() { return nOffset; }
};

// 
//  | A  B^T | | a | = | d |
//  | B   C  | | c |   | f |
//  C is the condensed part
// default setup condenses out all enrichment dofs

class DEMMatrices {
public:
	DEMMatrices() { B=C=f=0; }
	complex<double>* B;
	complex<double>* C;
	complex<double>* f;
};


class DEMCoreElement {
public:
	int *nn; // geometrical nodes
	bool storeMatricesF;
	bool condensedF;
	DEMLM **lm;
	int *bc;
	virtual ~DEMCoreElement() {
		if (bc) delete[] bc;
//   if (lm) delete[] lm;
		if (nn) delete[] nn;
		if (demm.B) delete[] demm.B;
		if (demm.C) delete[] demm.C;
		if (demm.f) delete[] demm.f;
	}
	virtual int defaultLMType() const =0;
	virtual bool condensedFlag() const =0;
	virtual bool dgmFlag() const = 0;
	virtual bool storeMatricesFlag()=0;
	virtual int nPolynomialDofs() const = 0;
	virtual int nEnrichmentDofs() const = 0;
	virtual int nGeomNodes() const = 0;
	virtual int nFaces() const =0;
	virtual int nFaceCorners(int fi) const =0;
	virtual int * faceCorners(int fi)=0;
	int nLagrangeDofs() const;
	virtual void createM(complex<double>*)=0;
	virtual void createRHS(complex<double>*)=0;
	virtual void createSol(double *xyz, complex<double>*, complex<double>*)=0;
	virtual void systemMatrix(complex<double>*);
	virtual void interfMatrix(int fi, DEMElement*,complex<double>*);
	virtual void systemRHS(complex<double>*);
	virtual void enrichmentF(double *x, complex<double> *f);
	virtual void polynomialF(double *x, double *f);

	DEMMatrices demm;
	virtual void condensedDofs(int &nc, int *&condensed, int *&kept);
	virtual int createCondensedMatrices(complex<double>* kk, complex<double>* A);

	int *ipiv; // used with one version of static condensation to store pivots
	void staticCondensationLHS(int ne, int nl, complex<double> *kee,
	                           complex<double> *kel, complex<double> *kll);
	void staticCondensationRHS(int ne, int nl, int nr,
	                           complex<double> *kee, complex<double> *kel,
	                           complex<double> *re, complex<double>* rl);
	void staticCondensationSolution(int ne, int nl,
	                                complex<double> *kee, complex<double> *kel,
	                                complex<double> *pl,
	                                complex<double> *re, complex<double> *pe);
};

class DEMElement: public Element, public DEMCoreElement {
public:
	int ne; // enrichment node offset
	virtual int nodesE(int no);
	virtual void nodalSolution(CoordSet &cs, complex<double>* sol,
	                           complex<double>(*nodalSol)[8]);
	complex<double> *forceVector;
// Element functions
	Category getCategory() const override { return Category::Acoustic; }

	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	int numDofs() const override {
		return ( (dgmFlag())?0:nPolynomialDofs() ) +
		       ( (condensedFlag())?0:nEnrichmentDofs() ) +
		       nLagrangeDofs() ;
	}
	virtual int polyDofType() const = 0;
	virtual int polyDofsPerNode() const = 0;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numNodes() const override {
		return ( (dgmFlag())?0:nGeomNodes() ) +
		       ( (condensedFlag())?0:nEnrichmentDofs() )+
		       ( nLagrangeDofs() );
	}

	int* nodes(int * = 0) const override;
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;

// Interface
	virtual double getOmega() const;
	virtual double *getWaveDirection() const;
	virtual complex<double> *getField() { return forceVector; }
	virtual void getNodalCoord(int n, int *nn, double *xyz) const;
	virtual double getRho() const { return prop->rho; }
	virtual double getSpeedOfSound() const;
	virtual double getE() const { return prop->E; }
	virtual double getNu() const { return prop->nu; }
	virtual PMLProps* getPMLProps() const { return &prop->fp; }
};

class DEMCoreInterfaceElement {
public:
	int fi;
	DEMElement *deme, *deme2;
	virtual void systemMatrix(complex<double>*);
};

#ifdef SANDIA
class DEMInterfaceElement: public DEMCoreInterfaceElement {
 DEMInterfaceElement(DEMElement *_deme, DEMElement *_deme2, int _fi);
};

#else

class DEMInterfaceElement: public Element, public DEMCoreInterfaceElement {
public:
	DEMInterfaceElement(DEMElement *_deme, DEMElement *_deme2, int _fi);

	int getElementType() const override { return -1; }
	Category getCategory() const override { return Category::Undefined; }
	void renum(const int *) override;

	void renum(EleRenumMap&) override;

	int numDofs() const override { return deme->numDofs()-deme->nLagrangeDofs()+
	                                      deme2->numDofs()-deme2->nLagrangeDofs(); }

	int* dofs(DofSetArray &, int *p) const override;

	void markDofs(DofSetArray &) const override;
	int numNodes() const override { return deme->numNodes()-deme->nLagrangeDofs()+
				deme2->numNodes()-deme2->nLagrangeDofs(); }
	virtual int* nodes(int *) const override;
};
#endif

#ifdef _TEMPLATE_FIX_
//#include <Element.d/Helm.d/DEMElementT.C>
#endif

#endif

