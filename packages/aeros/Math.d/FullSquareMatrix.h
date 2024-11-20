#ifndef	_FULLSQUAREMATRIX_H_
#define	_FULLSQUAREMATRIX_H_
#include <Utils.d/NodeSpaceArray.h>
#include <cstdio>
#include <Utils.d/MyComplex.h>
#include <Math.d/Vector.h>

template<class Scalar>
class GenFullSquareMatrix {
private:
  int	size;
  int     myval;
  Scalar*	value;
  int     length;
public:
  GenFullSquareMatrix(int, Scalar *l=0);
  GenFullSquareMatrix(); 
  GenFullSquareMatrix(GenFullSquareMatrix<Scalar> &m, Scalar s);
  ~GenFullSquareMatrix();

  // WARNING: Default copy constructor and assignement operator are used
  
  // Behavior is likely to be faulty if the original matrix is responsible
  // for managing its own memory (That is, when myval == true)
  
  // GenFullSquareMatrix(const GenFullSquareMatrix &) = default;
  // const GenFullSquareMatrix & operator=(const GenFullSquareMatrix &) = default;

  void setSize(int); 
  void changeSize(int i, int numMax);
  void reSize(int newSize); // Resize and preserve data

  void setMyval(int _myval) { myval = _myval; }

  GenFullSquareMatrix<Scalar> operator * (Scalar v);
  GenFullSquareMatrix<Scalar> operator / (Scalar v);
  GenFullSquareMatrix<Scalar> operator + (const GenFullSquareMatrix<Scalar> &M2);
  GenFullSquareMatrix<Scalar> operator - (const GenFullSquareMatrix<Scalar> &M2);
  GenFullSquareMatrix<Scalar> &operator *= (Scalar v);
  GenFullSquareMatrix<Scalar> &operator /= (Scalar v);
  GenFullSquareMatrix<Scalar> &operator += (const GenFullSquareMatrix<Scalar> &M2); 
  GenFullSquareMatrix<Scalar> &operator -= (const GenFullSquareMatrix<Scalar> &M2);
  GenFullSquareMatrix<Scalar> &operator += (const Tensor_d2s0 &t);
  Scalar *operator[] (int row);
  const Scalar *operator[] (int row) const;
  int dim() const   { return size; }
  int numRow() const { return size; }
  int numCol() const { return size; }
  void multiply(GenFullSquareMatrix<Scalar> &res, double d);
  void zero();
  void unitary();
  void unitaryDiag();
  void copy(const GenFullSquareMatrix<Scalar> &m);

  void symmetrize();
  void print(const char *msg = "", const char *msg2="");
  void printDiagonals();

  Scalar* data() { return value; }
  const Scalar* data() const { return value; }
  Scalar getValue(int i) const { return value[i]; }
  void setValue(int i, Scalar s) { value[i] = s; }
  void copy(Scalar *d);

  void add(GenFullSquareMatrix<Scalar> &m, int *rc);
 
  enum TransposeFlag {
    NORMAL = 0,
    TRANSPOSED
  };

  template<class Scalar1, class Scalar2, class Scalar3> void multiply(GenVector<Scalar1>& a, GenVector<Scalar2>& b, Scalar3 c, TransposeFlag transpose = NORMAL);
  
  void multiply(GenFullSquareMatrix<Scalar> &M2, GenFullSquareMatrix<Scalar> &result);
  void eigenVals(Scalar*);
  void eigenV(Scalar*);

  // Computes b <- M * a
  void mult(const GenVector<Scalar> &a, GenVector<Scalar> &b) const;

  // Computes b <- b + M * a
  void multAdd(const GenVector<Scalar> &a, GenVector<Scalar> &b) const;

  void linAdd(const Scalar &, const GenFullSquareMatrix<Scalar> &);
};

template<class Scalar>
inline
Scalar *
GenFullSquareMatrix<Scalar>::operator[](int row)
 { return value+size*row; }

template<class Scalar>
inline
const Scalar *
GenFullSquareMatrix<Scalar>::operator[](int row) const
 { return value+size*row; }


typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;

template<class Scalar> template<class Scalar1, class Scalar2, class Scalar3>
void GenFullSquareMatrix<Scalar>::multiply(GenVector<Scalar1>& a, GenVector<Scalar2>& b, Scalar3 c, typename GenFullSquareMatrix<Scalar>::TransposeFlag transposed)
{
	// multiply  b += c * matrix * a
	// Warning: no check for correct dimensions
	Scalar1* ap = a.data();
	Scalar2* bp = b.data();
	if (transposed == NORMAL) {
		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++)
			{
				bp[i] += c*(value+size*i)[j]*ap[j];
			}
		}
	} else {
		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++)
			{
				bp[j] += c*(value+size*i)[j]*ap[i];
			}
		}
	}
	return;
}

#endif

