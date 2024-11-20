#ifndef PITA_OLD_NICETIMER_H
#define PITA_OLD_NICETIMER_H

#include <deque>
#include <map>
#include <string>
#include <ostream>

namespace Pita { namespace Old {

/* Class NiceTimerHandler: Hierarchical timer with a tree-like dynamic structure.
   The structure expands when a timer that does not exist yet is started.
*/
class NiceTimerHandler
{
  public:
  explicit NiceTimerHandler(const std::string & name = "Total time");  
  explicit NiceTimerHandler(const char * name = "Total time");
  bool running() const;
  NiceTimerHandler & start(); // Start global timer
  NiceTimerHandler & startAll(); // Start global timer
  NiceTimerHandler & start(const std::string & name); // Start subtimer (create if necessary)
  NiceTimerHandler & start(const char * name);
  NiceTimerHandler & swap(const std::string & name); // Switch to subtimer (create if necessary)
  NiceTimerHandler & swap(const char * name);
  NiceTimerHandler & stop(); // Stop current subtimer
  NiceTimerHandler & stop(const std::string & name);
  NiceTimerHandler & stop(const char * name);
  NiceTimerHandler & stopAll(); // Stop timer & all active subtimers
  NiceTimerHandler & reset(); // Reset current subtimer, destroy its subtimers
  NiceTimerHandler & resetIteration(); // Reset current iteration
  NiceTimerHandler & resetAll(); // Reset (destroy) the whole structure
  NiceTimerHandler & newIteration(); // Start new iteration (Make current sutimer iterative if necessary)
  NiceTimerHandler & collapseIterations(bool flag = true); // Only count iterations & add times
  int iterationNumber() const;
  double time() const; // Get current time from current subtimer
  std::string status() const; // Get chain of currently running subtimers
  void print(std::ostream &) const; // Get complete current subtimer report
  void printAll(std::ostream &) const; // Get complete timer report
  
  private:
  class NiceTimerData;
    
  // Labelled subtimer
  class NiceTimerNode
  {
    public:
    std::string name;                     // Subtimer identifyer
    std::deque<NiceTimerData> timeData;   // Time values (if iterative, one for each iteration)
    bool iterative;                       // If true, subtimer is iterative
    bool collapseIterations;              // If true, sums all iteration timers

    explicit NiceTimerNode(const std::string & n = "Partial time");
  };
    
  // Subtimer data
  class NiceTimerData
  {
    public:
    double timeValue;                              // Timer value (in natural time units
    std::map<std::string, int> mapping;            // Name <-> Subtimer#
    std::deque<NiceTimerNode> subs;                // Subtimers

    explicit NiceTimerData();
  };
   
  static const double timeConversionRatio_; // Conversion milliseconds (natural time unit for FEM) -> seconds 
  
  // Subfunctions for print
  static void printSkeleton(std::ostream &, const std::deque<bool> &);
  static void printTimerNode(std::ostream &, const NiceTimerNode &, std::deque<NiceTimerNode *> &, std::deque<bool> &);
  static void printTimerData(std::ostream &, const NiceTimerData &, std::deque<NiceTimerNode *> &, std::deque<bool> &);
      
  NiceTimerNode root_;                    // Main subtimer
  std::deque<NiceTimerNode *> running_;   // Chain of timers currently running
};

// Output (C++ std I/O)
std::ostream & operator<<(std::ostream & out, const NiceTimerHandler & th);

// Inline functions
// ------------------------------------------------------------------

inline bool NiceTimerHandler::running() const
{
  return !running_.empty();
}

inline NiceTimerHandler & NiceTimerHandler::startAll()
{
  return start();
}

// Functions allowing C-style strings
inline NiceTimerHandler & NiceTimerHandler::start(const char * name)
{
  return start(std::string(name));
}

inline NiceTimerHandler & NiceTimerHandler::swap(const char * name)
{
  return swap(std::string(name));
}

inline NiceTimerHandler & NiceTimerHandler::stop(const char * name)
{
  return stop(std::string(name));
}

} /* end namespace Old */ } /* end namespace Pita */

#endif
