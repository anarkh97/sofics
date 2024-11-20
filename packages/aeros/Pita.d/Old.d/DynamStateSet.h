#ifndef PITA_OLD_DYNAMSTATESET_H
#define PITA_OLD_DYNAMSTATESET_H

#include "DynamState.h"

#include <vector>

template <typename Scalar> class GenSparseMatrix;

namespace Pita { namespace Old {

template <typename Scalar>
class DynamStateSet
{

  public:
  
    typedef Scalar DataType; 
    typedef DynamState<DataType> State;
    typedef GenVector<DataType> VectorType;
 
  public:

    // Constructors & destructor
    explicit DynamStateSet(int vectorSize = 0, int maxNumStates = 0) { init_(vectorSize, maxNumStates); }
    DynamStateSet(const DynamStateSet<Scalar> & dss) : vectorSize_(dss.vectorSize_), stateSet_(dss.stateSet_) {}
    ~DynamStateSet() {}

    // Accessors
    int numStates() const { return stateSet_.size(); }
    int maxNumStates() const { return stateSet_.capacity(); }
    int vectorSize() const { return vectorSize_; }
    State & operator[](int i) { return stateSet_[i]; }
    const State & operator[](int i) const { return stateSet_[i]; }

    // Global operations
    DynamStateSet<Scalar> & operator=(const DynamStateSet<Scalar> & dss) { copy(dss); return *this;}
    void reset(int vectorSize, int maxNumStates) { stateSet_.clear(); init_(vectorSize, maxNumStates); }
    void copy(const DynamStateSet<Scalar> &dss) { vectorSize_ = dss.vectorSize_; stateSet_ = dss.stateSet_; }
    void merge(const DynamStateSet<Scalar> &dss);
    void resize(int newMaxNumStates) { stateSet_.reserve(newMaxNumStates); }
    void clear() { stateSet_.clear(); }

    // Add states
    void addState(const DynamState<Scalar> & newState) { stateSet_.push_back(newState); }
    void addState() { stateSet_.push_back(DynamState<Scalar>(vectorSize_)); }
    
    // Extract to raw data array
    void getRaw(Scalar * buffer) const;

  private:
    // Data Fields
    int vectorSize_;
    std::vector<State> stateSet_; 

    // Memory management
    void init_(int vectorSize, int maxNumStates);
};

} /* end namespace Old */ } /* end namespace Pita */

#ifdef _TEMPLATE_FIX_
#include "DynamStateSet.C"
#endif

#endif
