#include "DynamState.h"
#include <Math.d/SparseMatrix.h>

namespace Pita { namespace Old {

/* Constructors & destructor */

template <typename Scalar>
DynamState<Scalar>::DynamState(const DynamState<Scalar> & ds, const GenSparseMatrix<Scalar> * K, const GenSparseMatrix<Scalar> * M)
{
  init(ds.vectorSize_);
  const_cast<GenSparseMatrix<double>*>(K)->mult(ds.dataArray_, dataArray_); 
  const_cast<GenSparseMatrix<double>*>(M)->mult(ds.dataArray_ + vectorSize_, dataArray_ + vectorSize_);
}

/* Memory management */

template <typename Scalar>
void DynamState<Scalar>::init(int vectorSize)
{
  if (vectorSize <= 0)
  {
    ownData_ = false;
    vectorSize_ = 0;
    dataArray_ = NULL;
    disp_ = NULL;
    vel_ = NULL;
  }
  else
  {
    ownData_ = true;
    vectorSize_ = vectorSize;
    dataArray_ = new Scalar[2 * vectorSize_];
    disp_ = new VectorType(vectorSize_, dataArray_, false);
    vel_ = new VectorType(vectorSize_, dataArray_ + vectorSize_, false);
  }
}

template <typename Scalar>
void DynamState<Scalar>::init(int vectorSize, Scalar * dataArray)
{
  ownData_ = false;
  if (vectorSize <= 0)
  {
    vectorSize_ = 0;
    dataArray_ = NULL;
    disp_ = NULL;
    vel_ = NULL;
  }
  else
  {
    vectorSize_ = vectorSize;
    dataArray_ = dataArray;
    disp_ = new VectorType(vectorSize_, dataArray_, false);
    vel_ = new VectorType(vectorSize_, dataArray_ + vectorSize_, false);
  }
}

template <typename Scalar>
void DynamState<Scalar>::setState(const VectorType & disp, const VectorType & vel)
{
  if (vectorSize_ != disp.size() || vectorSize_ !=  vel.size())
  {
    std::cerr << "Warning -- In void DynamState<Scalar>::setState(const VectorType &, const VectorType &) : Non matching lenghts\n";
    return;
  }
  disp_->copy(disp);
  vel_->copy(vel);
}

template <typename Scalar>
void DynamState<Scalar>::erase()
{
  if (vectorSize_ > 0)
  {
    delete vel_;
    delete disp_;
    if (ownData_)
    {
      delete[] dataArray_;
    }
  }
}

/* Operators */

template <typename Scalar>
DynamState<Scalar> & DynamState<Scalar>::operator=(Scalar value)
{
  std::fill(dataArray_, dataArray_ + 2 * vectorSize_, value);
  return *this;
}


template <typename Scalar>
DynamState<Scalar> & DynamState<Scalar>::operator+=(const DynamState<Scalar> & ds)
{
  int imax = 2 * vectorSize_;
  for (int i = 0; i < imax; ++i)
  {
    dataArray_[i] += ds.dataArray_[i];
  }
  return *this;
}

template <typename Scalar>
DynamState<Scalar> & DynamState<Scalar>::operator-=(const DynamState<Scalar> & ds)
{
  int imax = 2 * vectorSize_;
  for (int i = 0; i < imax; ++i)
  {
    dataArray_[i] -= ds.dataArray_[i];
  }
  return *this;
}

template <typename Scalar>
DynamState<Scalar> & DynamState<Scalar>::operator*=(Scalar coef)
{
  int imax = 2 * vectorSize_;
  for (int i = 0; i < imax; ++i)
  {
    dataArray_[i] *= coef;
  }
  return *this;
}

template <typename Scalar>
Scalar DynamState<Scalar>::operator*(const DynamState<Scalar> & ds) const 
{
  return ((*disp_) * (*(ds.disp_))) + ((*vel_) * (*(ds.vel_)));
}

template <typename Scalar>
void DynamState<Scalar>::linAdd(Scalar coef, const DynamState<Scalar> & ds)
{
  int imax = 2 * vectorSize_;
  for (int i = 0; i < imax; ++i)
  {
    dataArray_[i] += coef * ds.dataArray_[i];
  }
}

} /* end namespace Old */ } /* end namespace Pita */
