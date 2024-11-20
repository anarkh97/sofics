#ifndef SYMFULLMATRIX_H_
#define SYMFULLMATRIX_H_

#include <Utils.d/MyComplex.h> 

// Symmetric Full Matrix class
// Only the lower triangular part of the matrix is stored in a packed format 

template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;

template<class Scalar> 
class GenSymFullMatrix {
public :
  GenSymFullMatrix(int dim = 0);
  ~GenSymFullMatrix();
  
  int dim() const { return dim_; };
  void reSize(int newDim);
  void zero();

  // Caution: Must have first index greter or equal to second index
  Scalar *operator[](int row) { return v + rowStart(row); }             // Allows A[i][j] indexing, i >= j
  const Scalar *operator[](int row) const { return v + rowStart(row); } // Allows A[i][j] indexing, i >= j

  Scalar *data() { return v; }
  const Scalar *data() const { return v; }

  void add(const GenFullM<Scalar>&, int, int);
  void print() const;

protected:
  static int rowStart(int row) { return (row * (row + 1)) / 2; }
  int arrayLength() const { return rowStart(dim_); }

private:
  int dim_;
  Scalar *v;
};

typedef GenSymFullMatrix<double> SymFullMatrix;
typedef GenSymFullMatrix<DComplex> SymFullMatrixC;

#ifdef _TEMPLATE_FIX_
#include <Math.d/SymFullMatrix.C>
#endif

#endif
