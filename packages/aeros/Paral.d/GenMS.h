#ifndef _GEN_M_S_H_
#define _GEN_M_S_H_

#include<Math.d/CuCSparse.h>

/** \brief Set of non owning pointers to matrices that are part of the coarse problem */
template<class Scalar>
class GenMultiSparse : public GenSparseMatrix<Scalar> 
{
   GenSparseMatrix<Scalar> *K, *Kii;
   GenAssembledFullM<Scalar> *Kcc;
   GenCuCSparse<Scalar>      *Krc;
   GenCuCSparse<Scalar>      *Kib;
   GenDBSparseMatrix<Scalar> *Kbb;
  public:
   GenMultiSparse(GenSparseMatrix<Scalar> *_K, GenSparseMatrix<Scalar> *_Kii, 
                  GenDBSparseMatrix<Scalar> *_Kbb, GenCuCSparse<Scalar> *_Kib) 
   { K = _K; Kii = _Kii; Kbb = _Kbb; Kib = _Kib; Kcc = 0; Krc = 0; }

   GenMultiSparse(GenSparseMatrix<Scalar> *_K, GenSparseMatrix<Scalar> *_Kii, 
                  GenDBSparseMatrix<Scalar> *_Kbb, GenCuCSparse<Scalar> *_Kib,
                  GenCuCSparse<Scalar> *_Krc, GenAssembledFullM<Scalar> *_Kcc)
   { K = _K; Kii = _Kii; Kbb = _Kbb; Kib = _Kib; Krc = _Krc; Kcc = _Kcc; }

   GenMultiSparse(GenSparseMatrix<Scalar> *_K, GenCuCSparse<Scalar> *_Krc,
                  GenAssembledFullM<Scalar> *_Kcc)
   { K = _K; Krc = _Krc; Kcc = _Kcc; Kii = 0; Kbb = 0; Kib = 0; }
   
   ~GenMultiSparse() {};
    
   void add(const FullSquareMatrix & kel, const int *dofs) override;
   void add(const FullSquareMatrixC & kel, const int *dofs) override;
   void addImaginary(const FullSquareMatrix & kel, const int *dofs) override;
   void addDiscreteMass(int dof, Scalar mass) override;
   Scalar diag(int i) const override { return K->diag(i); }
   Scalar &diag(int i) override { return K->diag(i); }
   int dim() const override { return K->dim() ; }
   int neqs() const override { return K->neqs() ; }
   void zeroAll() override;
};

template<class Scalar>
void
GenMultiSparse<Scalar>::zeroAll() 
{
 if(K) K->zeroAll();
 if(Kbb) Kbb->zeroAll();
 if(Kib) Kib->zeroAll();
 if(Kii) Kii->zeroAll();
 if(Kcc) Kcc->zero();
 if(Krc) Krc->zeroAll();
}

template<class Scalar>
void
GenMultiSparse<Scalar>::add(const FullSquareMatrix & kel, const int *dofs)
{
 if(K)     K->add(kel, dofs);
 if(Krc) Krc->add(kel, dofs);
 if(Kcc) Kcc->add(kel, { dofs, kel.numRow() });
 if(Kbb) Kbb->add(kel, dofs);
 if(Kib) Kib->add(kel, dofs);
 if(Kii) Kii->add(kel, dofs);
}


template<class Scalar>
void
GenMultiSparse<Scalar>::add(const FullSquareMatrixC & kel, const int *dofs)
{
 if(K)     K->add(kel, dofs);
 if(Krc) Krc->add(kel, dofs);
 if(Kcc) Kcc->add(kel, dofs);
 if(Kbb) Kbb->add(kel, dofs);
 if(Kib) Kib->add(kel, dofs);
 if(Kii) Kii->add(kel, dofs);
}

template<class Scalar>
void
GenMultiSparse<Scalar>::addImaginary(const FullSquareMatrix & kel, const int *dofs)
{
 if(K)     K->addImaginary(kel, dofs);
 if(Krc) Krc->addImaginary(kel, dofs);
 if(Kcc) Kcc->addImaginary(kel, dofs);
 if(Kbb) Kbb->addImaginary(kel, dofs);
 if(Kib) Kib->addImaginary(kel, dofs);
 if(Kii) Kii->addImaginary(kel, dofs);
}

template<class Scalar>
void
GenMultiSparse<Scalar>::addDiscreteMass(int dof, Scalar mass)
{
// Check on if DOF is unconstrained must be done
// before this is called.
 if(K)     K->addDiscreteMass(dof, mass);
 if(Kbb) Kbb->addDiscreteMass(dof, mass);
 if(Kib) Kib->addDiscreteMass(dof, mass);
 if(Kii) Kii->addDiscreteMass(dof, mass);
 if(Kcc) Kcc->addDiscreteMass(dof, mass);
 if(Krc) Krc->addDiscreteMass(dof, mass);
}

//-----------------------------------------------------------------

#include <Math.d/DBSparseMatrix.h>
template <typename Scalar>
class SolverWrapper : public GenSparseMatrix<Scalar> {
public:

	SolverWrapper(const Connectivity *nToN, const DofSetArray *dsa, const ConstrainedDSA *c_dsa, GenSparseMatrix<Scalar> *wrapped) :
		multMatrix(nToN, dsa, c_dsa), wrapped(wrapped) {}

	void add(const FullSquareMatrix & kel, const int *dofs) override {
		multMatrix.add(kel, dofs);
		wrapped->add(kel, dofs);
	}
	void add(const FullSquareMatrixC & kel, const int *dofs) override {
		multMatrix.add(kel, dofs);
		wrapped->add(kel, dofs);
	}
	void addImaginary(const FullSquareMatrix & kel, const int *dofs) override {
		multMatrix.addImaginary(kel, dofs);
		wrapped->addImaginary(kel, dofs);
	}
	void addDiscreteMass(int dof, Scalar mass) override {
		multMatrix.addDiscreteMass(dof, mass);
		wrapped->addDiscreteMass(dof, mass);
	}
	void mult(const Scalar *rhs, Scalar *result) const override {
		multMatrix.mult(rhs, result);
	}
	Scalar diag(int i) const override { return multMatrix.diag(i); }
	Scalar &diag(int i) override { return multMatrix.diag(i); }
	int dim() const override { return multMatrix.dim() ; }
	int neqs() const override { return multMatrix.neqs() ; }
	void zeroAll() override {
		multMatrix.zeroAll();
		wrapped->zeroAll();
	}
private:
	GenDBSparseMatrix<Scalar> multMatrix;
    GenSparseMatrix<Scalar> *wrapped;
};
#endif
