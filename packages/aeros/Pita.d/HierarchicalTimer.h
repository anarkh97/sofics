#ifndef PITA_HIERARCHICALTIMER_H
#define PITA_HIERARCHICALTIMER_H

#include "Fwk.h"

#include <deque>
#include <map>
#include <string>
#include <ostream>

namespace Pita {

// Class HierarchicalTimer: Hierarchical timer with a tree-like dynamic structure.
// The structure expands when a timer that does not exist yet is started.
class HierarchicalTimer : public Fwk::PtrInterface<HierarchicalTimer> {
public:
  EXPORT_PTRINTERFACE_TYPES(HierarchicalTimer);

  bool running() const;
  HierarchicalTimer & start(); // Start global timer
  HierarchicalTimer & startAll(); // Start global timer
  HierarchicalTimer & start(const String & name); // Start subtimer (create if necessary)
  HierarchicalTimer & swap(const String & name); // Switch to subtimer (create if necessary)
  HierarchicalTimer & stop(); // Stop current subtimer
  HierarchicalTimer & stop(const String & name);
  HierarchicalTimer & stopAll(); // Stop timer & all active subtimers
  HierarchicalTimer & reset(); // Reset current subtimer, destroy its subtimers
  HierarchicalTimer & resetIteration(); // Reset current iteration
  HierarchicalTimer & resetAll(); // Reset (destroy) the whole structure
  HierarchicalTimer & newIteration(); // Start new iteration (Make current sutimer iterative if necessary)
  HierarchicalTimer & collapseIterations(bool flag = true); // Only count iterations & add times
  int iterationNumber() const;
  double time() const; // Get current time from current subtimer
  String status() const; // Get chain of currently running subtimers
  void print(std::ostream &) const; // Get complete current subtimer report
  void printAll(std::ostream &) const; // Get complete timer report
  
  HierarchicalTimer();
  explicit HierarchicalTimer(const String & name);

private:
  class TData;
    
  // Labelled subtimer
  struct TNode {
    String name;                     // Subtimer identifyer
    std::deque<TData> timeData;   // Time values (if iterative, one for each iteration)
    bool iterative;                       // If true, subtimer is iterative
    bool collapseIterations;              // If true, sums all iteration timers

    explicit TNode(const String & n = "Partial time");
  };
    
  // Subtimer data
  struct TData {
    double timeValue;                              // Timer value (in natural time units
    std::map<String, int> mapping;            // Name <-> Subtimer#
    std::deque<TNode> subs;                // Subtimers

    explicit TData();
  };
   
  static const double timeConversionRatio_; // Conversion milliseconds (natural time unit for FEM) -> seconds 
  
  // Subfunctions for print
  static void printSkeleton(std::ostream &, const std::deque<bool> &);
  static void printTimerNode(std::ostream &, const TNode &, std::deque<TNode *> &, std::deque<bool> &);
  static void printTimerData(std::ostream &, const TData &, std::deque<TNode *> &, std::deque<bool> &);
      
  TNode root_;                    // Main subtimer
  std::deque<TNode *> running_;   // Chain of timers currently running

  DISALLOW_COPY_AND_ASSIGN(HierarchicalTimer);
};

// Output (C++ std I/O)
std::ostream &
operator<<(std::ostream & out, const HierarchicalTimer & th);

// Inline functions
// ------------------------------------------------------------------

inline
bool
HierarchicalTimer::running() const {
  return !running_.empty();
}

inline
HierarchicalTimer &
HierarchicalTimer::startAll() {
  return start();
}

} /* end namespace Pita */

#endif /* PITA_HIERARCHICALTIMER_H */
