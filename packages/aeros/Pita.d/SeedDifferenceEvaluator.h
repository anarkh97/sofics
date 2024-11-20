#ifndef PITA_SEEDDIFFERENCEEVALUATOR_H
#define PITA_SEEDDIFFERENCEEVALUATOR_H

#include "Fwk.h"
#include "Seed.h"
#include "DynamOps.h"

#include "PitaNonLinDynam.h"

namespace Pita {

class SeedDifferenceEvaluator : public SharedState<DynamState>::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(SeedDifferenceEvaluator);

  virtual void onIteration(); // overriden

  const Seed * seedDifference() const { return notifier(); } 
  const Seed * referenceSeed() const { return referenceSeed_.ptr(); }

  virtual void referenceSeedIs(const Seed * rs) { setReferenceSeed(rs); }

protected:
  SeedDifferenceEvaluator(const Seed * seedDifference, const Seed * referenceSeed);
  
  void setReferenceSeed(const Seed * rs) { referenceSeed_ = rs; }

  virtual double getRefMagnitude() const = 0;
  virtual double getDiffMagnitude() const = 0;

private:
  Seed::PtrConst referenceSeed_;
  double firstDiffMagnitude_;

  DISALLOW_COPY_AND_ASSIGN(SeedDifferenceEvaluator);
};

class LinSeedDifferenceEvaluator : public SeedDifferenceEvaluator {
public:
  EXPORT_PTRINTERFACE_TYPES(LinSeedDifferenceEvaluator);

  class Manager;

  const DynamOps * metric() const { return metric_.ptr(); }

protected:
  LinSeedDifferenceEvaluator(const DynamOps * metric, const Seed * seedError, const Seed * referenceSeed);

  virtual double getRefMagnitude() const;
  virtual double getDiffMagnitude() const;
  
  friend class Manager;

private:
  DynamOps::PtrConst metric_;
};


class LinSeedDifferenceEvaluator::Manager : public Fwk::GenManager<LinSeedDifferenceEvaluator, const Seed *> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  const DynamOps * defaultMetric() const { return defaultMetric_.ptr(); }

  static Ptr New(const DynamOps * defaultMetric) {
    return new Manager(defaultMetric);
  }

protected:
  explicit Manager(const DynamOps * defaultMetric);

  virtual LinSeedDifferenceEvaluator::Ptr createNewInstance(const KeyType & key); // overriden

private:
  DynamOps::PtrConst defaultMetric_;
};

class NonLinSeedDifferenceEvaluator : public SeedDifferenceEvaluator {
public:
  EXPORT_PTRINTERFACE_TYPES(NonLinSeedDifferenceEvaluator);

  class Manager;

  NonLinSeedDifferenceEvaluator(PitaNonLinDynamic * probDesc,
                                const Seed * seedDifference,
                                const Seed * referenceSeed,
                                Seconds refTime = Seconds(0.0));
protected:
  virtual double getRefMagnitude() const;
  virtual double getDiffMagnitude() const;
  
private:
  PitaNonLinDynamic * probDesc_;
  Seconds refTime_;
};


class NonLinSeedDifferenceEvaluator::Manager : public Fwk::GenManager<NonLinSeedDifferenceEvaluator, const Seed *> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  static Ptr New(PitaNonLinDynamic * probDesc) {
    return new Manager(probDesc);
  }

protected:
  explicit Manager(PitaNonLinDynamic * probDesc) :
    probDesc_(probDesc)
  {}

  virtual NonLinSeedDifferenceEvaluator::Ptr createNewInstance(const KeyType & key); // overriden

private:
  PitaNonLinDynamic * probDesc_;
};

} /* end namespace Pita */

#endif /* PITA_SEEDDIFFERENCEEVALUATOR_H */
