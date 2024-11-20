#ifndef PITA_OLD_DYNAMSTATE_H
#define PITA_OLD_DYNAMSTATE_H

#include <algorithm>

template <typename Scalar> class GenVector;
template <typename Scalar> class GenSparseMatrix;

namespace Pita { namespace Old {

template <typename Scalar>
class DynamState
{
public:
  typedef Scalar DataType;
  typedef GenVector<Scalar> VectorType;

  // Constructors & destructor
  DynamState(int vectorSize = 0) { init(vectorSize); }
  DynamState(int vectorSize, Scalar initValue) { init(vectorSize, initValue); }
  DynamState(int vectorSize, Scalar * dataArray) { init(vectorSize, dataArray); }
  DynamState(const DynamState<Scalar> & ds) { init(ds.vectorSize_); copyData(ds.dataArray_);  }
  DynamState(const DynamState<Scalar> &, const GenSparseMatrix<Scalar> * K, const GenSparseMatrix<Scalar> * M);
  ~DynamState() { erase(); }

  // Global operations
  void reset(int vectorSize = 0) { erase(); init(vectorSize); }
  void reset(int vectorSize, Scalar initValue) { erase(); init(vectorSize, initValue); }
  void reset(int vectorSize, Scalar * dataArray) { erase(); init(vectorSize, dataArray); }
  void clear() { reset(0); }

  // Accessors
  const VectorType & disp() const { return *disp_; }
  VectorType & disp() { return *disp_; }
  const VectorType & vel() const { return *vel_; }
  VectorType & vel() { return *vel_; }
  int vectorSize() const { return vectorSize_; }
  const Scalar * data() const { return dataArray_; }
    
  // Mutator
  void setState(const VectorType & disp, const VectorType & vel);

  // Operators
  DynamState<Scalar> & operator=(const DynamState<Scalar> & ds) { reset(ds.vectorSize_); copyData(ds.dataArray_); return *this; }
  DynamState<Scalar> & operator=(Scalar);
  DynamState<Scalar> & operator+=(const DynamState<Scalar> &);
  DynamState<Scalar> & operator-=(const DynamState<Scalar> &);
  DynamState<Scalar> & operator*=(Scalar);
  Scalar operator*(const DynamState<Scalar> &) const;
  void linAdd(Scalar, const DynamState<Scalar> &);
 
  // Raw data manipulation
  Scalar * data() { return dataArray_; }
  void setRaw(const Scalar * buffer) { copyData(buffer); }
  void getRaw(Scalar * buffer) const { std::copy(dataArray_, dataArray_ + 2 * vectorSize_, buffer); }

protected:
  int vectorSize_;
  Scalar * dataArray_;
  VectorType * disp_;
  VectorType * vel_;
  bool ownData_;

  // Memory management
  void init(int vectorSize);
  void init(int vectorSize, Scalar initValue) { init(vectorSize); operator=(initValue); }
  void init(int vectorSize, Scalar * dataArray);
  void copyData(const Scalar * dataArray) { std::copy(dataArray, dataArray + 2 * vectorSize_, dataArray_); }
  void erase();
};

} /* end namespace Old */ } /* end namespace Pita */

#ifdef _TEMPLATE_FIX_
#include "DynamState.C"
#endif

#endif
