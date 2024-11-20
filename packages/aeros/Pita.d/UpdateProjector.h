#ifndef PITA_UPDATEPROJECTOR_H
#define PITA_UPDATEPROJECTOR_H

#include "Fwk.h"
#include "Seed.h"

namespace Pita {

class UpdateProjector : public Fwk::PtrInterface<UpdateProjector> {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdateProjector);
  typedef Fwk::GenManagerInterface<UpdateProjector *, String> Manager;

  virtual size_t reducedBasisSize() const = 0;

  /* Results */
  const Vector & reducedJump() const { return reducedJump_; }
  const Vector & updateComponents() const { return updateComponents_; }

  /* Sources */
  const Seed * jump() const { return jump_.ptr(); }
  const ReducedSeed * reducedCorrection() const { return reducedCorrection_.ptr(); }
  
  virtual void jumpIs(const Seed * j) = 0;
  virtual void reducedCorrectionIs(const ReducedSeed * rc) = 0;

protected:
  UpdateProjector() :
    reducedJump_(),
    updateComponents_(),
    jump_(NULL),
    reducedCorrection_(NULL)
  {}

  void setJump(const Seed * j) { jump_ = j; }
  void setReducedCorrection(const ReducedSeed * rc) { reducedCorrection_ = rc; }

private:
  Vector reducedJump_;
  Vector updateComponents_;
  Seed::PtrConst jump_;
  ReducedSeed::PtrConst reducedCorrection_;

  DISALLOW_COPY_AND_ASSIGN(UpdateProjector);
};

} /* end namespace Pita */

#endif /* PITA_UPDATEPROJECTOR_H */
