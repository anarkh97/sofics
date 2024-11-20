#ifndef PITA_POSTPROCESSINGMANAGER_H
#define PITA_POSTPROCESSINGMANAGER_H

#include "Fwk.h"
#include "PostProcessor.h"
#include "DynamPropagator.h"
#include "DynamTimeIntegrator.h"

#include <map>

namespace Pita {

namespace PostProcessing {

// class IntegratorReactor

class IntegratorReactor : public DynamTimeIntegrator::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(IntegratorReactor);

  class Builder : public Fwk::PtrInterface<Builder> {
  public:
    EXPORT_PTRINTERFACE_TYPES(Builder);

    virtual PostProcessorRoot * postProcessor() const = 0;
    
    virtual IntegratorReactor * reactorNew(const DynamTimeIntegrator * notifier, PostProcessor::FileSetId outputFileSet) const = 0;

  protected:
    Builder() {}

  private:
    DISALLOW_COPY_AND_ASSIGN(Builder);
  };

  PostProcessor::FileSetId outputFileSet() const { return outputFileSet_; }

protected:
  IntegratorReactor(const DynamTimeIntegrator * notifier, PostProcessor::FileSetId outputFileSet) :
    DynamTimeIntegrator::NotifieeConst(notifier),
    outputFileSet_(outputFileSet)
  {}

private:
  PostProcessor::FileSetId outputFileSet_;

  DISALLOW_COPY_AND_ASSIGN(IntegratorReactor);
};

// class PropagatorReactor

class PropagatorReactor : public DynamPropagator::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagatorReactor);

  class Manager : public Fwk::PtrInterface<Manager> {
  public:
    EXPORT_PTRINTERFACE_TYPES(Manager);
    
    const IntegratorReactor::Builder * reactorBuilder() const { return reactorBuilder_.ptr(); }

    PostProcessor::FileSetId outputFileSet(const DynamPropagator * observedPropagator) const;
    void outputFileSetIs(const DynamPropagator * observedPropagator, PostProcessor::FileSetId fileSet);

    static Ptr New(IntegratorReactor::Builder * reactorBuilder) { 
      return new Manager(reactorBuilder);
    }

  protected:
    explicit Manager(IntegratorReactor::Builder * reactorBuilder);
    
    friend class PropagatorReactor;

  private:
    IntegratorReactor::Builder::Ptr reactorBuilder_;

    typedef std::map<const DynamPropagator *, PropagatorReactor::Ptr> PropagatorReactorContainer;
    PropagatorReactorContainer propagatorReactor_;

    DISALLOW_COPY_AND_ASSIGN(Manager);
  };

  virtual void onInitialState();
  virtual void onFinalState();

  PostProcessor::FileSetId fileSet() const { return fileSet_; }
  void fileSetIs(PostProcessor::FileSetId fs) { fileSet_ = fs; }

protected:
  PropagatorReactor(const DynamPropagator * notifier,
      Manager * parent,
      PostProcessor::FileSetId fileSet);
  
  friend class Manager;

private:
  const DynamPropagator * notifier_; // Need a specialized DynamPropagator
  Manager * parent_;
  PostProcessor::FileSetId fileSet_;

  DynamTimeIntegrator::NotifieeConst::Ptr integratorReactor_; // Collects the output data
};

typedef PropagatorReactor::Manager Manager;

// class  IntegratorReactorImpl

template <typename PostProcessorType>
class IntegratorReactorImpl : public IntegratorReactor {
public:
  EXPORT_PTRINTERFACE_TYPES(IntegratorReactorImpl);

  typedef typename PostProcessorType::Integrator Integrator;

  class Builder : public IntegratorReactor::Builder {
  public:
    EXPORT_PTRINTERFACE_TYPES(Builder);
  
    PostProcessorType * postProcessor() const { return postProcessor_.ptr(); }

    virtual IntegratorReactorImpl<PostProcessorType> * reactorNew(
        const DynamTimeIntegrator * integrator,
        PostProcessor::FileSetId fileSet) const {
      return new IntegratorReactorImpl<PostProcessorType>(
          dynamic_cast<const typename PostProcessorType::Integrator*>(integrator),
          postProcessor_.ptr(),
          fileSet);
    }

    static Ptr New(PostProcessorType * postProcessor) {
      return new Builder(postProcessor);
    }

  protected:
    Builder(PostProcessorType * postProcessor) :
      postProcessor_(postProcessor)
    {}

  private:
    typename PostProcessorType::Ptr postProcessor_;
  };

  const PostProcessorType * target() { return target_; }

  virtual void onInitialCondition() {
    // Reset file(s) before writing initial state
    this->target_->fileStatusIs(this->outputFileSet(), PostProcessorRoot::CLOSED);
    this->target_->fileStatusIs(this->outputFileSet(), PostProcessorRoot::OPEN);
    this->performOutput();
  }

  virtual void onCurrentCondition() {
    this->performOutput();
  }

  virtual void notifierIs(const DynamTimeIntegrator * notifier) {
    // If notifier is not actually a PostProcessorType::Integrator, do not accept notifications
    this->notifier_ = dynamic_cast<const Integrator *>(notifier);
    IntegratorReactor::notifierIs(notifier_);
  }

protected:
  IntegratorReactorImpl(const Integrator * n, 
                        PostProcessorType * t, 
                        PostProcessor::FileSetId ofs);

  void performOutput() const {
    this->target_->outputNew(this->outputFileSet(), this->notifier_);
  }

  friend class Builder;

private:
  const Integrator * notifier_; // Need a specialized DynamTimeIntegrator
  PostProcessorType * target_;
  PostProcessor::FileSetId outputFileSet_;
};

template <typename PostProcessorType>
IntegratorReactorImpl<PostProcessorType>::IntegratorReactorImpl(
    const IntegratorReactorImpl<PostProcessorType>::Integrator * notifier,
    PostProcessorType * target,
    PostProcessor::FileSetId outputFileSet) :
  IntegratorReactor(NULL, outputFileSet),
  target_(target)
{
  this->notifierIs(notifier);
  if (!this->target() || this->target()->fileStatus(this->outputFileSet()) == PostProcessorRoot::NO_FILE) {
    throw Fwk::RangeException("in PostProcessing::IntegratorReactorImpl::IntegratorReactorImpl");
  }
}

} // end namespace PostProcessing

} // end namespace Pita

#endif /* PITA_POSTPROCESSINGMANAGER_H */
