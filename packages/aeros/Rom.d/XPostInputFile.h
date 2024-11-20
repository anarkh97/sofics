#ifndef ROM_XPOSTINPUTFILE_H
#define ROM_XPOSTINPUTFILE_H

#include <string>
#include <cstdio>
#include <vector>

namespace Rom {

template <int> class NodeDofBuffer;
typedef NodeDofBuffer<6> NodeDof6Buffer;

class XPostInputFile {
public:
  // Static parameters
  const std::string &fileName() const { return fileName_; }
  
  int stateCount() const { return stateCount_; }
  
  int nodeCount() const { return nodeCount_; }
  typedef std::vector<int>::const_iterator NodeIdIterator;
  NodeIdIterator nodeIdBegin() const { return nodeIndices_.begin(); }
  NodeIdIterator nodeIdEnd() const { return nodeIndices_.end(); }

  // Iteration and retrieval
  int currentStateIndex() const { return currentStateIndex_; }
  double currentStateHeaderValue() const { return currentStateHeaderValue_; }
  // Node ordering determined by NodeDof6Buffer, file node ordering ignored
  const NodeDof6Buffer &currentStateBuffer(NodeDof6Buffer &target) const;
  
  void currentStateIndexInc(); // Must have called currentStateBuffer() at least once before calling currentStateIndexInc()

  bool validCurrentState() const { return currentStateIndex() < stateCount(); }

  // Ctor & dtor
  explicit XPostInputFile(const std::string &fileName);

  ~XPostInputFile();

private:
  const std::string fileName_;
  
  int nodeCount_;
  std::vector<int> nodeIndices_;

  int stateCount_;
  
  FILE *stream_;
  
  int currentStateIndex_;
  double currentStateHeaderValue_;
  mutable bool currentStateRead_;
  std::fpos_t currentStatePosition_;

  void readNodeIndices();
  void readCurrentStateHeader();
  void positionAtStateStart() const;

  // Disallow copy & assignment
  XPostInputFile(const XPostInputFile&);
  XPostInputFile& operator=(const XPostInputFile&);
};

} /* end namespace Rom */

#endif /* ROM_BASISINPUTFILE_H */
