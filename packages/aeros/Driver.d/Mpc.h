#ifndef _LMPC_H_
#define _LMPC_H_

#include <iostream>
#include <iomanip>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif
#include <Utils.d/MyComplex.h>
#include <Utils.d/resize_array.h>
#include <Corotational.d/GeomState.h>
#include <Element.d/Element.h>
#include <Utils.d/Connectivity.h>
#include <vector>

struct RealOrComplex
{
	double r_value;
	DComplex c_value;
};

class LMPCTerm
{
public:
	bool isComplex;
	int nnum;           // node number
	int dofnum;         // dof number (0-6)
	RealOrComplex coef; // value of coefficient

	// real constructor
	LMPCTerm(int _nnum, int _dofnum, double _coef);

	// complex constructor
	LMPCTerm(int _nnum, int _dofnum, double _rcoef, double _icoef);

	// copy constructors
	LMPCTerm(const LMPCTerm &t);
	LMPCTerm(const LMPCTerm &t, bool _isComplex);

// JLchange 
	LMPCTerm(LMPCTerm &t, double weight);

	// default constructor
	LMPCTerm(bool _isComplex = false);

	// convert a real term to a complex term if it is not already complex
	void makeComplex();

	template <class Scalar>
	Scalar val() {
		if(isComplex) {
			Scalar s;
			ScalarTypes::copy(s, coef.c_value);
			return s;
		}
		else
			return coef.r_value;
	}

	bool isNull() {
		return isComplex ? coef.c_value == 0.0 : coef.r_value == 0.0;
	}
};


template<class Scalar>
class GenLMPCTerm
{
public:
	int nnum;           // node number
	int dofnum;         // dof number (0-6)
	Scalar coef;        // value of coefficient
	int dof;            // local subdomain dof (dsa numbering), -1 if doesn't exist
	int cdof;           // local subdomain dof (c_dsa numbering), -1 if constrained or doesn't exist
	int ccdof;          // local subdomain dof (cc_dsa numbering), -1 if corner, constrained or doesn't exist
	GenLMPCTerm() { dof = cdof = ccdof = -1; }
};

namespace mpc {
enum ConstraintType { Equality, Inequality };
enum ConstraintSource { Lmpc, NodalContact, TiedSurfaces, ContactSurfaces, Rigid, Joint, FetiPrimal, FetiDual, Undefined };
}

/** Linear Multi-Point Constraint class */
class LMPCons
{
private:
	mpc::ConstraintType m_type;
	mpc::ConstraintSource m_source;

public:
	bool isComplex;
	RealOrComplex rhs;              // right hand side of mpc
	RealOrComplex original_rhs;
	std::vector<LMPCTerm> terms;    // terms of the mpc (node, dof & coef)
	union {
		int lmpcnum;                  // id number of the mpc from input
		int fluid_node;
		int pairNb;
	};
	int nterms;                     // number of terms in mpc

	// TODO avoid having state variables.
	mutable int type;               // 0: equality constraint
	// 1: inequality constraint (contact)

	int psub, nsub;                 // subdomains involved in Feti constraint.
	// psub will have coef +1.0, nsub will have coef -1.0

	int lagrangeMult;               // -1 (undefined) or 0/1 (multiplier on/off)
	double penalty;

	std::pair<int,int> id;          // unique identifier
	// currently defined for source == ContactSurfaces as follows:
	// id.first = mortar handler id, id.second = slave node global id
#ifdef USE_EIGEN3
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> H;
#endif
	// real constructor
	LMPCons(int _lmpcnum, double _rhs, LMPCTerm *term0 = 0);

	// complex constructor
	LMPCons(int _lmpcnum, double rrhs, double irhs, LMPCTerm *term0 = 0);

	// templated data access functions
	template<class Scalar> Scalar getRhs();
	template<class Scalar> GenLMPCTerm<Scalar> getTerm(int i);

	// make rhs and all terms complex if they aren't already
	void makeComplex();

	// add term
	void addterm(LMPCTerm *term);

	// memory
	long mem() { return nterms*(4+1+1)+1; } //HB

	bool isPrimalMPC() { return ((type == 3) || (type == 4) || (type == 5)); }
	bool isBoundaryMPC() { return ((type == 2) || (type == 5)); }

	void print();

	/** remove the zero terms in this constraint */
	void removeNullTerms();

	void setType(mpc::ConstraintType _type) { m_type = _type; }
	mpc::ConstraintType getType() { return m_type; }
	void setSource(mpc::ConstraintSource _source) { m_source = _source; }
	mpc::ConstraintSource getSource() const { return m_source; }

	void normalize();
};

template<> double LMPCons::getRhs<double>();
template<> DComplex LMPCons::getRhs<DComplex>();
template<> GenLMPCTerm<double> LMPCons::getTerm<double>(int i);
template<> GenLMPCTerm<DComplex> LMPCons::getTerm<DComplex>(int i);


// PJSA & HB: create specific class for subdomain LMPC to handle the lmpc stiffness scaling
// -> store for each term of the SUBDOMAIN lmpc, its index (position) in the ORIGINAL DOMAIN lmpc
//    store the ORIGINAL number of terms of the DOMAIN lmpc the SUBDOMAIN lmpc comes from
//    store for each term of the SUBDOMAIN lmpc, the diag stiffness value of the associated dof
//    store for each term of the SUBDOMAIN lmpc, the sum of the diag stiffness values (from the shared dofs)

template<class Scalar>
class SubLMPCons
{
private:
	mpc::ConstraintType m_type;
	mpc::ConstraintSource m_source;
	SubLMPCons(const SubLMPCons<Scalar> &) { std::cerr << "SubLMPCons copy constructor is not implemented \n"; }
public:
	Scalar rhs;
	Scalar original_rhs;
	ResizeArray<GenLMPCTerm<Scalar> > terms;
	union {
		int lmpcnum;  // for distributed input using the CU_Feti interface, this is the global mpc number
		int fluid_node;
		int pairNb;
	};
	int nterms;

	int type;                       // 0: equality constraint
	// 1: inequality constraint (contact)

	std::pair<int,int> id;          // unique identifier
	// currently defined for source == ContactSurfaces as follows:
	//   id.first = mortar handler id, id.second = slave node global id

	bool active;                    // defines active set, used for contact
	// in FETI-DP active means than the lagrange multiplier associated with the mpc is constrained to be zero (ie the dual constraint is active)

	ResizeArray<int> gi;      // index of mpc term before it is distributed to subd
	int gsize;                // number of terms in mpc before distributing to subds
	ResizeArray<Scalar> k;    // diag stiffness for term dof, used for stiffness scaling
	Scalar *ksum;             // sum of diag stiffness for term dof, used for stiffness scaling

	SubLMPCons(int _lmpcnum, Scalar _rhs, GenLMPCTerm<Scalar> term0, int _gsize, int _gi) : terms(term0, 1), gi(0,1), k(1.0,1)
	{
		lmpcnum = _lmpcnum;
		original_rhs = _rhs;
		rhs = _rhs;
		nterms = 1;
		ksum = 0;
		gsize = _gsize;
		gi[0] = _gi;
		ScalarTypes::initScalar(k[0], 1.0);
		type = 0;
		active = false;
		m_type = mpc::Equality;
		m_source = mpc::Undefined;
	}

	virtual ~SubLMPCons() { if(ksum) delete [] ksum; }

	void addterm(GenLMPCTerm<Scalar> term, int _gi)
	{
		int i=0;
		// --- Verify if term already exists
		while((i<nterms) &&
		      ((terms[i].nnum!=term.nnum) |
		       (terms[i].dofnum!=term.dofnum)))
			i++;
		// if term not already implied, add term
		if(i==nterms) {
			terms[nterms++]=term;
			gi[i] = _gi;
			ScalarTypes::initScalar(k[i], 1.0);  // default value (used for topo scaling, i.e. A=I)
			// overwritten for stiffness scaling (i.e. A=1/diag(K))
		}
			// if term already implied, add coefficients
		else terms[i].coef += term.coef;
	}

	void initKsum()
	{
		ksum = new Scalar[nterms];
		for(int i=0; i<nterms; i++) ScalarTypes::initScalar(ksum[i], 1.0);
	}

	void print()
	{
		std::cerr << "lmpcnum = " << lmpcnum << ", rhs = " << rhs << ", nterms = " << nterms << std::endl;
		for(int i=0; i<nterms; ++i)
			std::cerr << "  term " << i+1 << ": node " << terms[i].nnum+1 << "  dof "
			          << terms[i].dofnum << "  coef " << terms[i].coef << std::endl;
	}

	void setType(mpc::ConstraintType _type) { m_type = _type; }
	mpc::ConstraintType getType() { return m_type; }
	void setSource(mpc::ConstraintSource _source) { m_source = _source; }
	mpc::ConstraintSource getSource() { return m_source; }
};

inline void LMPCons::removeNullTerms() {
	std::vector<LMPCTerm>::iterator i = terms.begin();
	while(i != terms.end()) {
		if(i->isNull())
			i = terms.erase(i);
		else
			++i;
	}
	nterms = terms.size();
}

template<>
class SetAccess<LMPCons>
{
	int numLMPC;
	//ResizeArray<LMPCons *>&lmpc;
	LMPCons **lmpc;
	std::map<std::pair<int,int>, int> &dofID;
public:
	//SetAccess(int n, ResizeArray<LMPCons *>&l, std::map<std::pair<int,int>, int> &d)
	SetAccess(int n, LMPCons **l, std::map<std::pair<int,int>, int> &d)
			: numLMPC(n), lmpc(l), dofID(d) {}
	int size() const {
		return numLMPC;
	}
	int numNodes(int i) const {
		return lmpc[i]->nterms;
	}
	void nodes(int i, int *nd) const {
		for(int j = 0; j < lmpc[i]->nterms; ++j)
			nd[j] = dofID[std::pair<int,int>(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum)];
	}
};

#endif
