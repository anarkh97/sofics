
/*****************************************************************************
 *                   Copyright (C) 1999 CMSoft                               *
 *                                                                           *
 *  These lines of code and declarations contain unpublished proprietary     *
 *  information of CMSoft. They may not be copied or duplicated in whole     *
 *  or part without prior authorization from CMSoft.                         *
 *                                                                           *
 *****************************************************************************/

#ifndef FS_FULL_MATRIX_H_
#define FS_FULL_MATRIX_H_

#include <Utils.d/MyComplex.h>

// ======================================================== //
//                                                          //
// Author: K. H. Pierson                                    //
//                                                          //
// ======================================================== //
// Description:                                             //
//                                                          //
//             FSFullMatrix = Feti Solver Full Matrix class //
//             stores an mxn full matrix, certain member    //
//             functions only work for the square full      //
//             matrix case (nxn)                            //
// ======================================================== //

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
template <class Scalar> class GenStackVector;
typedef GenStackVector<double> StackVector;

template<class Scalar>
class GenFSFullMatrix {

 protected:
   int nrow;	// number of rows
   int ncolumn; // number of columns
   int ld;      // leading dimension
   Scalar *v;   // pointer to matrix data
 public:

   // Constructors
   GenFSFullMatrix();                               // Creates an empty full matrix
   GenFSFullMatrix(int nr);                         // Creates M(nr,nr)
   GenFSFullMatrix(int nr, int nc);                 // Creates M(nr,nc)
   GenFSFullMatrix(int nr, int nc, Scalar initVal); // Creates M(nr,nc)=initVal
   GenFSFullMatrix(int nr, int nc, Scalar *array);  // Creates M(nr,nc)=array(i,j)

   // Copy Constructors
   GenFSFullMatrix(const GenFSFullMatrix<Scalar> &A); // Creates a new full matrix equal to A
   GenFSFullMatrix(const GenFSFullMatrix<Scalar> &A, int nr, int sr, int nc, int sc);

   // Destructor
   ~GenFSFullMatrix();

   void copy(Scalar *array);
   void reSize(int nr, int nc, Scalar d=0.0); // resize to M(nr,nc)=0.0

   // Operators
   Scalar* operator[](int i)const;
   void  operator = (const GenFSFullMatrix<Scalar> &B);
   void  operator = (const Scalar c);
   void  operator *= (const Scalar c);

   // matrix-matrix multiplication
   GenFSFullMatrix<Scalar> operator *(GenFSFullMatrix<Scalar> &B); // product C=A*B
   GenFSFullMatrix<Scalar> operator ^(GenFSFullMatrix<Scalar> &B); // product C=A^T*B
   GenFSFullMatrix<Scalar> operator %(GenFSFullMatrix<Scalar> &B); // product C=A*B^T

   GenFSFullMatrix<Scalar> transpose();

   // matrix-vector multiplication routines (calls BLAS-3 routines)
   // y=alpha*A*x + beta*y
   void mult(Scalar *x, Scalar *y, Scalar alpha=1.0,Scalar beta=0.0);
   // y=alpha*A^T*x + beta*y
   void trMult(const Scalar *x, Scalar *y, Scalar alpha=1.0,Scalar beta=0.0 );
   // y=alpha*A*x + beta*y
   void mult(  const GenVector<Scalar> &x, GenVector<Scalar> &y, Scalar alpha=1.0,Scalar beta=0.0);
   // y=alpha*A^T*x + beta*y
   void trMult(const GenVector<Scalar> &x, GenVector<Scalar> &y, Scalar alpha=1.0,Scalar beta=0.0);

   void mult( const ComplexVector &x, StackVector &y) {}

   // general matrix-matrix multiplication routines (calls BLAS-3 routines)
   // C=alpha*A*B + beta*C
   void mult(    const GenFSFullMatrix<Scalar> &B, GenFSFullMatrix<Scalar> &C,
                 Scalar alpha=1.0, Scalar beta=0.0);
   // C=alpha*A^T*B + beta*C
   void trMult(  const GenFSFullMatrix<Scalar> &B, GenFSFullMatrix<Scalar> &C,
                 Scalar alpha=1.0, Scalar beta=0.0);
   // C=alpha*A*B^T + beta*C
   void multTr(  const GenFSFullMatrix<Scalar> &B, GenFSFullMatrix<Scalar> &C,
                 Scalar alpha=1.0, Scalar beta=0.0);
   // C=alpha*A^T*B^T + beta*C
   void trMultTr(const GenFSFullMatrix<Scalar> &B, GenFSFullMatrix<Scalar> &C,
                 Scalar alpha=1.0, Scalar beta=0.0);


   GenFSFullMatrix<Scalar> invert();      // form A^-1
   double max();               // maximum matrix entry
   void luFactor();            // perform LU factorization
   void symLuFactor();            // perform LU factorization
   int  symLuFactor(int *perm, double tol, Scalar *origDiag, int *sing=0);
   void Lm1Mult(Scalar *a, int nc, int lda);
   void Um1Mult(Scalar *a, int nc, int lda);
   void Um1TMult(Scalar *a, int nc, int lda);
   void Lm1Mult(Scalar *a, int nc, int lda, int *perm);
   void Um1Mult(Scalar *a, int nc, int lda, int *perm);
   void Um1TMult(Scalar *a, int nc, int lda, int *perm);
   //void solve(Scalar *rhs);    // perform forward/backward [Deactivated since buggy]
   void zero();                // zero all matrix entries
   void identity();            // set square full matrix to identity matrix
   void print(const char *msg = "",const char *msg2=""); // print the full matrix

   void add(GenFSFullMatrix<Scalar>& subMatrix, int sr, int sc);

   // inline functions
   int dim()      const { return nrow;    }
   int ldim()     const { return ld;      }
   int numRow()   const { return nrow;    }
   int numCol()   const { return ncolumn; }
   Scalar* data() const { return v;       }
};

template<class Scalar>
inline
Scalar *
GenFSFullMatrix<Scalar>::operator[](int i) const
 { return v+i*ld; }

template<class Scalar>
class GenStackFSFullMatrix : public GenFSFullMatrix<Scalar> {
 public:
   GenStackFSFullMatrix(const GenFSFullMatrix<Scalar> &matrix);
   GenStackFSFullMatrix(const GenFSFullMatrix<Scalar> &matrix, int nr, int nc,
          int rOff, int cOff);
   GenStackFSFullMatrix(int nr, Scalar *data);
   GenStackFSFullMatrix(int nr, int nc, Scalar *data);
   GenStackFSFullMatrix(int nr, int nc, int ldim, Scalar *data);
   ~GenStackFSFullMatrix() { this->v = 0; }
   GenFSFullMatrix<Scalar>&  operator = (const Scalar c)
      { GenFSFullMatrix<Scalar>::operator=(c); return *this; }
};

template<class Scalar>
inline
GenStackFSFullMatrix<Scalar>::GenStackFSFullMatrix(const GenFSFullMatrix<Scalar> &mat)
{
 this->nrow    = mat.numRow();
 this->ncolumn = mat.numCol();
 this->ld      = mat.ldim();
 this->v       = mat.data();
}

template<class Scalar>
inline
GenStackFSFullMatrix<Scalar>::GenStackFSFullMatrix(const GenFSFullMatrix<Scalar> &mat,
          int nr, int nc, int rOff, int cOff)
{
 this->nrow    = nr;
 this->ncolumn = nc;
 this->ld      = mat.ldim();
 this->v       = mat.data() + rOff*this->ld + cOff;
}


template<class Scalar>
inline
GenStackFSFullMatrix<Scalar>::GenStackFSFullMatrix(int nr, Scalar *data)
{
 this->nrow    = nr;
 this->ncolumn = nr;
 this->ld      = nr;
 this->v       = data;
}

template<class Scalar>
inline
GenStackFSFullMatrix<Scalar>::GenStackFSFullMatrix(int nr, int nc, Scalar *data)
{
 this->nrow    = nr;
 this->ncolumn = nc;
 this->ld      = nc;
 this->v       = data;
}

template<class Scalar>
inline
GenStackFSFullMatrix<Scalar>::GenStackFSFullMatrix(int nr, int nc, int ldim, Scalar *data)
{
 this->nrow    = nr;
 this->ncolumn = nc;
 this->ld      = ldim;
 this->v       = data;
}

typedef GenFSFullMatrix<double> FSFullMatrix;
typedef GenStackFSFullMatrix<double> StackFSFullMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/FullMatrix.C>
#endif

#endif
