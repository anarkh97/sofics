#include <sstream>
#include <Timers.d/GetTime.h>
#include "HierarchicalTimer.h"

namespace Pita {

const double HierarchicalTimer::timeConversionRatio_ = 1.0e-3;

HierarchicalTimer::TNode::TNode(const String & n) :
   name(n),
   timeData(1),
   iterative(false),
   collapseIterations(false)
{}

HierarchicalTimer::TData::TData() :
  timeValue(0.0)
{}

HierarchicalTimer::HierarchicalTimer() :
  root_("Total time")
{}

HierarchicalTimer::HierarchicalTimer(const String & name) :
  root_(name) 
{}

HierarchicalTimer &
HierarchicalTimer::start() {
  if(!running()) {
    running_.push_back(&root_);
    running_.back()->timeData.back().timeValue -= getTime();
  }
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::start(const String & id) {
  if (running()) { 
    std::map<String, int> & mapping = running_.back()->timeData.back().mapping;
    std::map<String, int>::iterator it = mapping.find(id);
    int num;
    if (it != mapping.end()) {
      num = it->second;
    } else {
      num = running_.back()->timeData.back().subs.size();
      running_.back()->timeData.back().subs.push_back(TNode(id));
      mapping[id] = num;
    }
    running_.push_back(&(running_.back()->timeData.back().subs.at(num)));
    running_.back()->timeData.back().timeValue -= getTime();
  } else {
    if (id == root_.name) {
      return start();
    }
  }
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::swap(const String & id) {
  stop();
  return start(id);
}

HierarchicalTimer &
HierarchicalTimer::stop() {
  if (running()) {
    running_.back()->timeData.back().timeValue += getTime();
    running_.pop_back();
  }
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::stop(const String & id) {
  if (running()) {
    std::deque<TNode *>::reverse_iterator it_target;
    int numTimersToStop = 0;
    for (it_target = running_.rbegin(); it_target != running_.rend(); ++it_target) {
      ++numTimersToStop;
      if ((*it_target)->name == id) break;
    }
    if (it_target != running_.rend()) {
      for (int i = 0; i < numTimersToStop; ++i) {
        stop();
      }
    }
  }
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::stopAll() {
  while (running()) stop();
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::reset() {
  if (running()) {
    running_.back()->timeData.clear();
    running_.back()->timeData.push_back(TData());
    running_.back()->timeData.back().timeValue -= getTime();
  } else {
    return resetAll();
  }
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::resetIteration() {
  if (running()) {
    running_.back()->timeData.back().subs.clear();
    running_.back()->timeData.back().timeValue = - getTime();
  } else {
    return resetAll();
  }
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::resetAll() {
  root_.timeData.clear();
  root_.timeData.resize(1);
  running_.clear();
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::newIteration() {
  if (running()) {
    if (running_.back()->iterative) {
      running_.back()->timeData.back().timeValue += getTime();
      running_.back()->timeData.push_back(TData());
      running_.back()->timeData.back().timeValue -= getTime();
    } else {
      running_.back()->iterative = true;
    }
  } else {
    if (root_.iterative) {
      root_.timeData.push_back(TData());
    } else {
      root_.iterative = true;
    }
  }
  return *this;
}

HierarchicalTimer &
HierarchicalTimer::collapseIterations(bool flag) {
  if (running()) {
    running_.back()->collapseIterations = flag;
  } else {
    root_.collapseIterations = flag;
  }
  return *this;
}

int
HierarchicalTimer::iterationNumber() const {
  if (running()) {
    return running_.back()->iterative ? running_.back()->timeData.size() : 0;
  } else {
    return root_.iterative ? root_.timeData.size() : 0;
  }
}

double
HierarchicalTimer::time() const {
  return running() ? running_.back()->timeData.back().timeValue + getTime() : root_.timeData.back().timeValue;
}

String
HierarchicalTimer::status() const {
  std::stringstream s;
  if (running()) {
    std::deque<TNode *>::const_iterator it = running_.begin();
    s << (*it)->name;
    for (++it; it != running_.end(); ++it) {
      s << " -> " << (*it)->name;
    }
    s << " (running)";
  } else {
    s << root_.name << " (stopped)";
  }
  return s.str();
}

void
HierarchicalTimer::printSkeleton(OStream & out, const std::deque<bool> & skel) {
  for (std::deque<bool>::const_iterator it = skel.begin(); it != skel.end(); ++it) {
    out << (*it ? '|' : ' ') << ' ';
  }
}

void
HierarchicalTimer::printTimerData(OStream & out, const TData & td, std::deque<TNode*> & run, std::deque<bool> & skel) {
  double actualTime;
  if (run.empty()) {
    actualTime = td.timeValue * timeConversionRatio_;
  } else {
    actualTime = (td.timeValue + getTime()) * timeConversionRatio_;
    run.pop_front();
  }
  out << ": " << actualTime << " s" << std::endl;
  for(std::deque<TNode>::const_iterator it_node = td.subs.begin(); it_node != td.subs.end(); ++it_node) {
    printSkeleton(out, skel);
    bool lastSub = (it_node + 1 == td.subs.end());
    out << "|-";
    skel.push_back(!lastSub);
    printTimerNode(out, *it_node, run, skel);
    skel.pop_back();
  }
}

void
HierarchicalTimer::printTimerNode(OStream & out, const TNode & tn, std::deque<TNode*> & run, std::deque<bool> & skel)
{
  out << tn.name;
  bool runFlag = (!run.empty() && run.front() == &tn);
  std::deque<TNode*> runEmpty;
  if (tn.iterative) {
    if (tn.collapseIterations) {
      out << ": ";
      double totalTime = 0.0;
      for (std::deque<TData>::const_iterator it_iter = tn.timeData.begin(); it_iter != tn.timeData.end(); ++it_iter) {
        totalTime += it_iter->timeValue;
      }
      double actualTime = (runFlag ? totalTime + getTime() : totalTime) * timeConversionRatio_;
      out << actualTime << " s (For " << tn.timeData.size() << " iterations)" << std::endl;
    } else {
      out << std::endl;
      int numIter = 1;
      for(std::deque<TData>::const_iterator it_iter = tn.timeData.begin(); it_iter != tn.timeData.end(); ++it_iter) {
        printSkeleton(out, skel);
        out << "|-Iteration # " << numIter++;
        skel.push_back(it_iter + 1 != tn.timeData.end());
        printTimerData(out, *it_iter, (runFlag && !skel.back()) ? run : runEmpty, skel);
        skel.pop_back();
      }
    }
  } else {
    printTimerData(out, tn.timeData.back(), (runFlag ? run : runEmpty), skel);
  }
}

void
HierarchicalTimer::print(OStream & out) const {
  if (running()) {
    std::deque<TNode*> run(1, running_.back());
    std::deque<bool> skel;
    printTimerNode(out, *(running_.back()), run, skel);
  } else {
    printAll(out);
  }
}

void
HierarchicalTimer::printAll(OStream & out) const {
  std::deque<TNode*> run(running_);
  std::deque<bool> skel;
  printTimerNode(out, root_, run, skel);
}

OStream &
operator<<(OStream & out, const HierarchicalTimer & th) {
  th.printAll(out);
  return out;
}

} /* end namespace Pita */
