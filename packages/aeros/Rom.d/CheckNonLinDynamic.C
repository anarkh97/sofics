#include "CheckNonLinDynamic.h"

#include "VecBasis.h"
#include "VecBasisFile.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "BasisOps.h"

#include <Driver.d/Domain.h>

#include <cstddef>

// DEBUG
#include <iostream>

namespace Rom {

class CheckNonLinDynamic::Impl {
public:
  explicit Impl(Domain *);

  void lastStateSnapIs(const GeomState &state);

private:
  Domain * domain_;
  VecBasis projectionBasis_;

  Vector displacement_;
  Vector projectionError_;
  Vector components_;

  int skipCounter_;

  // Disallow copy & assignment
  Impl(const Impl &);
  const Impl &operator=(const Impl &);
};

CheckNonLinDynamic::Impl::Impl(Domain *d) :
  domain_(d),
  skipCounter_(0)
{
  FileNameInfo fileInfo;
  VecNodeDof6Conversion converter(*domain_->getCDSA());
  BasisInputStream<6> podFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), converter);
  podFile >> projectionBasis_;

  displacement_.reset(projectionBasis_.vectorSize());
  projectionError_.reset(projectionBasis_.vectorSize());
  components_.reset(projectionBasis_.vectorCount());
}

void
CheckNonLinDynamic::Impl::lastStateSnapIs(const GeomState &state) {
  ++skipCounter_;
  if (skipCounter_ >= domain_->solInfo().skipPodRom) {
    const_cast<GeomState &>(state).get_tot_displacement(displacement_);
    reduce(projectionBasis_, displacement_, components_);
    expand(projectionBasis_, components_, projectionError_);
    projectionError_ -= displacement_;

    // DEBUG
    std::cout << "(Projection Error) / Displacement = " << projectionError_.norm() / displacement_.norm() << "\n";

    skipCounter_ = 0;
  }
}

CheckNonLinDynamic::CheckNonLinDynamic(Domain *domain) :
  NonLinDynamic(domain),
  impl_(nullptr)
{}

CheckNonLinDynamic::~CheckNonLinDynamic() {
  // Nothing to do
}

void
CheckNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  impl_.reset(new CheckNonLinDynamic::Impl(this->domain));
}

void
CheckNonLinDynamic::handleStateSnap(const GeomState &state) {
  impl_->lastStateSnapIs(state);
}

} /* end namespace Rom */
