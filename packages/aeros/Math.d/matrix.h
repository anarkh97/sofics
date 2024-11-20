#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cstdio>
#include <vector>
#include <gsl/span>
#include <Utils.d/MyComplex.h>

// FullM = Full Matrix class
//         stores an mxn matrix
//         certain member functions only work for
//         the square matrix case (nxn)

template <class Scalar> class GenStackFullM;
typedef GenStackFullM<double> StackFullM;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;

// TODO Separate the storage aspect and the ?GECP factorization aspects.
// TODO Merge storage with storage of GenFullSquareMatrix
template<class Scalar>
class GenFullM {
protected:
    int nrow;	// number of rows
    int ncolumn; // number of columns
    Scalar *v;   // pointer to matrix data
    int *iprow; //!< Row pivoting indices from ?gecp routine
    int *ipcol; //!< Column pivoting indices from ?gecp routine
    int ndef;

    void subAddProd(int i, int nThreads, GenFullM<Scalar> *A, GenFullM<Scalar> *B,
                    Scalar coefC, Scalar coefAB);
    void subAddUTProd(int i, int nThreads, GenFullM<Scalar> *A, GenFullM<Scalar> *B,
                      Scalar coefC, Scalar coefAB);
    void subAddTTProd(int i, int nThreads, GenFullM<Scalar> *A, GenFullM<Scalar> *B,
                      Scalar coefC, Scalar coefAB);
    void subZero(int i, int nThreads);

public:

    // constructors
    GenFullM();               // Creates an empty matrix
    explicit GenFullM(int nr);         // Creates an NxN matrix
    GenFullM(int nr, int nc); // Creates an NxM matrix
    GenFullM(int nr, int nc, Scalar init_val); // Creates an NxM matrix initialized a[i][j] = init_val
    GenFullM(const GenFullM<Scalar> &m, int nr, int sr, int nc, int sc);
    GenFullM(const GenFullM<Scalar> &m);
    GenFullM(const GenFullM<Scalar> &m, int nr, int *rows, int nc, int *cols);
    GenFullM(const GenFullM<Scalar> &m, int nr, int *rows, int nc, int sc);
    GenFullM(const GenFullM<Scalar> &m, int nr, int sr, int nc, int *cols);
    GenFullM(const GenFullM<Scalar> &m, double scale_factor);
    GenFullM(Scalar *data, int nr, int nc, int flag = 0);

    // destructor
    ~GenFullM();

    void copy(Scalar *array);
    void setNewSize(int nr, int nc, Scalar d=0.0);
    void setNewSize(int nr, Scalar d=0.0);

    // OPERATORS
    void  operator =  (const GenFullM<Scalar> &);
    void  operator =  (const Scalar c);
    void  operator *=  (const Scalar c);
    GenFullM<Scalar>  operator *(const GenFullM<Scalar>&) const;   // product A*B
    GenFullM<Scalar>  operator *=(const GenFullM<Scalar>&);   // product A=A*B
    GenFullM<Scalar>  operator ^=(const GenFullM<Scalar>&);   // product A=B*A
    GenVector<Scalar> operator *(const GenVector<Scalar>&) const; // product A*x
    GenFullM<Scalar>  operator ^(const GenFullM<Scalar>&) const;   // product A^T*B
    GenVector<Scalar> operator ^(const GenVector<Scalar>&) const; // product A^T*x
    GenFullM<Scalar>  operator %(const GenFullM<Scalar>&) const;   // product A*B^T
    GenFullM<Scalar>  operator +(const GenFullM<Scalar>&) const;   // A+B = [A B]^T
    GenFullM<Scalar>  operator +=(const GenFullM<Scalar>&);
    GenFullM<Scalar>  operator -=(const GenFullM<Scalar>&);

    void  mult(  Scalar *x, Scalar *y, Scalar alpha=1.0, Scalar beta=0.0);
    void  trMult(Scalar *x, Scalar *y, Scalar alpha=1.0, Scalar beta=0.0);

    GenFullM<Scalar> invert();
    GenFullM<Scalar> Invert(double tol=1.0e-6); // invert using gaussian elimination with full pivoting
    GenFullM<Scalar> transpose();

    int dim() const { return nrow;    }
    int numRow() const { return nrow;    }
    int numCol() const { return ncolumn; }
    int getNdef() { return ndef; }

    Scalar *operator[](int i) const;
    Scalar* data() const { return v; }
    Scalar* getData() { return v; }
    Scalar* Column(int i) const;

    Scalar & operator() (int row, int col) { return (*this)[row][col]; }
    const Scalar operator() (int row, int col) const { return (*this)[row][col]; }

    double max();
    double min();
    double maxAbs();
    void symmetrize_from_uptriag(); // Symmetrize by   lower diagonal terms = upper diag terms
    void print(const char *msg = "", const char *msg2="", FILE* f=stderr);
    void Print();
    void factor();
    void Factor(double tol=1.0e-6, bool print_ndef = true);
    void reSolve(Scalar *d);
    void ReSolve(Scalar *d); // resolve using factorization by gaussian elimination with full pivoting
    void zero();
    void clean_up();
    void add(const GenFullM<Scalar>&, int, int);
    void add(GenVector<Scalar>&, int, int);
    void add(const GenFullM<Scalar>&, int *);
    void addrows(GenFullM<Scalar>&, int *);
    void subtract(GenFullM<Scalar>&, int, int);
    void negate();

    void paralAddProd(GenFullM<Scalar> &A, GenFullM<Scalar> &B, Scalar coefC, Scalar coefAB);
    void paralAddUTProd(GenFullM<Scalar> &A, GenFullM<Scalar> &B, Scalar coefC, Scalar coefAB);
    void paralAddTTProd(GenFullM<Scalar> &A, GenFullM<Scalar> &B, Scalar coefC, Scalar coefAB);
    void paralZero();
    void transposeAssign(GenFullM<Scalar>&);
    void transposeMult(GenFullM<Scalar>&, GenFullM<Scalar> &);
    void transposeMultD(GenFullM<Scalar> &, GenVector<Scalar> &, GenFullM<Scalar> &);

    Scalar diag(int i) { return (*this)[i][i]; }
};

template<class Scalar>
inline
Scalar *
GenFullM<Scalar>::operator[](int i) const
{ return v+i*ncolumn; }

template<class Scalar>
class GenStackFullM : public GenFullM<Scalar> {
public:
    GenStackFullM(int nr, Scalar *data);
    GenStackFullM(Scalar *data, int nr);
    GenStackFullM(int nr, int nc, Scalar *data);
    GenStackFullM(Scalar *data, int nr, int nc);
    ~GenStackFullM() { this->v=0; }
};

template<class Scalar>
inline
GenStackFullM<Scalar>::GenStackFullM(int nr,Scalar *data)
{
 this->nrow    = nr;
 this->ncolumn = nr;
 this->v       = data;
}

template<class Scalar>
inline
GenStackFullM<Scalar>::GenStackFullM(Scalar *data,int nr)
{
 this->nrow    = nr;
 this->ncolumn = nr;
 this->v       = data;
}

template<class Scalar>
inline
GenStackFullM<Scalar>::GenStackFullM(int nr, int nc, Scalar *data)
{
 this->nrow    = nr;
 this->ncolumn = nc;
 this->v       = data;
}

template<class Scalar>
inline
GenStackFullM<Scalar>::GenStackFullM(Scalar *data, int nr, int nc)
{
 this->nrow    = nr;
 this->ncolumn = nc;
 this->v       = data;
}


template<class Scalar>
class GenAssembledFullM : public GenFullM<Scalar>
{
	std::vector<int> rowMap;
	std::vector<int> colMap;
public:
	GenAssembledFullM(int nr, gsl::span<const int> rMap);
	GenAssembledFullM(int nr, int *rMap, int nc, int *cMap);
	~GenAssembledFullM() { clean_up(); };
	/** \brief Add an elemental contribution.
	 *
	 * @param kel Elemental matrix.
	 * @param dofs Dofs of the elemental matrix.
	 */
	void add(const FullSquareMatrix &kel, gsl::span<const int> dofs);
	void add(const FullSquareMatrixC &kel, gsl::span<const int> dofs) {
		add(kel, dofs.data());
	}
	void add(const FullSquareMatrixC &kel, const int *dofs);
	void addImaginary(const FullSquareMatrix &kel, const int *dofs);
	void add(Scalar **matrix);
	void addDiscreteMass(int dof, Scalar dmass);
	void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
	void clean_up();
};

template<class Scalar>
inline
GenAssembledFullM<Scalar>::GenAssembledFullM(int nr, gsl::span<const int> rMap) :
	rowMap(rMap.begin(), rMap.end()), colMap(rMap.begin(), rMap.end())
{
 this->nrow    = nr;
 this->ncolumn = nr;
 this->v       = new Scalar[this->nrow*this->ncolumn];
 this->zero();
}

template<class Scalar>
inline
GenAssembledFullM<Scalar>::GenAssembledFullM(int nr, int *rMap, int nc, int *cMap) :
	rowMap(rMap, rMap+nr), colMap(cMap, cMap+nc)
{
 this->nrow    = nr;
 this->ncolumn = nc;
 this->v       = new Scalar[this->nrow*this->ncolumn];
 this->zero();
}

template<class Scalar>
inline
void
GenAssembledFullM<Scalar>::clean_up()
{
 GenFullM<Scalar>::clean_up();
}

typedef GenFullM<double> FullM;
typedef GenStackFullM<double> StackFullM;
typedef GenAssembledFullM<double> AssembledFullM;
typedef GenFullM<DComplex> FullMC;
typedef GenAssembledFullM<DComplex> AssembledFullMC;

#endif
