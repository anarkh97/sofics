#include <cstdio>
#include <Math.d/matrix.h>
#include <Math.d/SymFullMatrix.h>

template<class Scalar>
GenSymFullMatrix<Scalar>::GenSymFullMatrix(int dim)
{
  dim_ = dim;
  v = new Scalar[this->arrayLength()];
}

template<class Scalar>
GenSymFullMatrix<Scalar>::~GenSymFullMatrix()
{
  delete[] v;
}

template<class Scalar>
void
GenSymFullMatrix<Scalar>::reSize(int newDim)
{
  if (dim() == newDim)
    return;

  int newLength = rowStart(newDim);
  int copyCount = std::min(newLength, this->arrayLength());

  const Scalar * temp = v;
  dim_ = newDim;
  v = new Scalar[newLength];
  std::copy(temp, temp + copyCount, v);

  delete[] temp;
}

template<class Scalar>
void
GenSymFullMatrix<Scalar>::zero()
{
  int size = this->arrayLength();
  for (int i = 0; i < size; ++i)
    v[i] = 0.0;
}

template<class Scalar>
void
GenSymFullMatrix<Scalar>::add(const GenFullM<Scalar> &mat, int fRow, int fCol)
{
  int nrow = mat.numRow();
  int ncol = mat.numCol();

  int iCol, iRow;

  for(iCol = 0; iCol < ncol; ++iCol) {
    int rowStop = std::min(nrow, fCol+iCol-fRow+1);
    for(iRow = 0; iRow < rowStop; ++iRow)
      (*this)[fCol+iCol][fRow+iRow] += mat[iRow][iCol];
  }
}

template<class Scalar>
void
GenSymFullMatrix<Scalar>::print() const
{
  for(int i = 0; i < dim_; ++i)
     for (int j = 0; j <= i; ++j)
        fprintf(stderr, " sym[%2d][%2d] = %14.7e\n",i,j,(*this)[i][j]);
}

