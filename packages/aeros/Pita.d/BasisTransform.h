#ifndef PITA_BASISTRANSFORM_H
#define PITA_BASISTRANSFORM_H

#include "Fwk.h"
#include "DynamStatePlainBasis.h"

#include "SimpleBuffer.h"

namespace Pita {
  
// Common interface

class BasisTransform : public Fwk::PtrInterface<BasisTransform> {
public:
  EXPORT_PTRINTERFACE_TYPES(BasisTransform);

  size_t vectorSize() const { return output_->vectorSize(); }
  
  virtual void inputValueIs(const DynamStateBasis & iv) = 0; // Pass-by-value

  const DynamStatePlainBasis * output() const { return output_.ptr(); }
  DynamStatePlainBasis * output() { return output_.ptr(); }
  void outputIs(DynamStatePlainBasis * o) { output_ = o; }

protected:
  explicit BasisTransform(size_t vectorSize) :
    output_(DynamStatePlainBasis::New(vectorSize))
  {}

private:
  DynamStatePlainBasis::Ptr output_;

  DISALLOW_COPY_AND_ASSIGN(BasisTransform);
};

// Implementations

class BasisCopy : public BasisTransform {
public:
  EXPORT_PTRINTERFACE_TYPES(BasisCopy);

  explicit BasisCopy(size_t vectorSize) :
    BasisTransform(vectorSize)
  {}

  virtual void inputValueIs(const DynamStateBasis & iv); // overriden
};

class GramSchmidtOrtho : public BasisTransform {
public:
  EXPORT_PTRINTERFACE_TYPES(GramSchmidtOrtho);

  explicit GramSchmidtOrtho(size_t vectorSize, double energyTol);

  virtual void inputValueIs(const DynamStateBasis & iv); // overriden

private:
  double energyTol_;
};

class CholeskyOrtho : public BasisTransform {
public:
  EXPORT_PTRINTERFACE_TYPES(CholeskyOrtho);

  explicit CholeskyOrtho(size_t vectorSize, double energyTol);

  virtual void inputValueIs(const DynamStateBasis & iv); // overriden

private:
  double energyTol_;
};


class QROrtho : public BasisTransform {
public:
  EXPORT_PTRINTERFACE_TYPES(QROrtho);

  explicit QROrtho(size_t vectorSize, double energyTol);

  virtual void inputValueIs(const DynamStateBasis & iv); // overriden

private:
  double energyTol_;
  SimpleBuffer<double> data_;
  SimpleBuffer<int> permutation_;
  SimpleBuffer<double> tau_;
  
  mutable SimpleBuffer<double> workspace_;
};

} /* end namespace Pita */

#endif /* PITA_BASISTRANSFORM_H */
