#ifndef ROM_BASISBINARYFILE_H
#define ROM_BASISBINARYFILE_H

#include "NodeDof6Buffer.h"

#include <Utils.d/BinaryResultFile.h>

#include <string>
#include <algorithm>
#include <iterator>
#include <cassert>

namespace Rom {

class BasisBinaryFile {
protected:
  typedef BinFileHandler::OffType OffsetType;
  
  static const double VERSION;
  static const std::string DESC;
  static const int NODAL_DATA_FLAG;
  static const int DEFAULT_DOFS_PER_NODE;
};

class BasisBinaryOutputFile : public BasisBinaryFile {
public:
  const std::string &fileName() const { return binFile_.pathName(); }

  int nodeCount() const { return binFile_.itemCount(); }

  typedef BinaryResultOutputFile::ItemIdIterator NodeIdIterator;
  NodeIdIterator nodeIdBegin() const { return binFile_.itemIdBegin(); }
  NodeIdIterator nodeIdEnd() const { return binFile_.itemIdEnd(); }

  int stateCount() const { return binFile_.stateCount(); }

  // Node ordering determined by NodeDofBuffer, file node ordering ignored
  template<int DOFS_PER_NODE>
  void stateAdd(const NodeDofBuffer<DOFS_PER_NODE> &data);
  template<int DOFS_PER_NODE>
  void stateAdd(const NodeDofBuffer<DOFS_PER_NODE> &data, double headValue);
  
  // Node indices: [0, nodeCount) 
  BasisBinaryOutputFile(const std::string &fileName, int nodeCount, bool restart,
                        int dofs_per_node = DEFAULT_DOFS_PER_NODE);
  
  // Node indices: [first, last)
  template <typename NodeIdIt>
  BasisBinaryOutputFile(const std::string &fileName, NodeIdIt first, NodeIdIt last, bool restart,
                        int dofs_per_node = DEFAULT_DOFS_PER_NODE);

private:
  BinaryResultOutputFile binFile_;

  // Disallow copy & assignment
  BasisBinaryOutputFile(const BasisBinaryOutputFile&);
  BasisBinaryOutputFile& operator=(const BasisBinaryOutputFile&);
};

template <typename NodeIdIt>
BasisBinaryOutputFile::BasisBinaryOutputFile(const std::string &fileName, NodeIdIt first, NodeIdIt last, bool restart,
                                             int dofs_per_node) :
  binFile_(fileName, NODAL_DATA_FLAG, DESC, std::distance(first, last), dofs_per_node, 0, first, last, VERSION, restart)
{}

template<int DOFS_PER_NODE>
void
BasisBinaryOutputFile::stateAdd(const NodeDofBuffer<DOFS_PER_NODE> &data) {
  stateAdd(data, static_cast<double>(stateCount() + 1));
}


class BasisBinaryInputFile : public BasisBinaryFile {
public:
  // Static parameters
  const std::string &fileName() const { return binFile_.pathName(); }
  
  int nodeCount() const { return binFile_.itemCount(); }
  
  typedef BinaryResultInputFile::ItemIdIterator NodeIdIterator;
  NodeIdIterator nodeIdBegin() const { return binFile_.itemIdBegin(); }
  NodeIdIterator nodeIdEnd() const { return binFile_.itemIdEnd(); }
  
  int stateCount() const { return binFile_.stateCount(); }

  // Iteration and retrieval
  int currentStateIndex() const { return binFile_.stateRank(); }
  double currentStateHeaderValue() const { return currentStateHeaderValue_; }
  
  // Node ordering determined by NodeDofBuffer, file node ordering ignored
  template<int DOFS_PER_NODE>
  const NodeDofBuffer<DOFS_PER_NODE> &currentStateBuffer(NodeDofBuffer<DOFS_PER_NODE> &target);
  // Must have called currentStateBuffer() at least once before calling currentStateIndexInc()
  // NOTE: This restriction could be lifted if necessary
 
  void currentStateIndexInc();

  bool validCurrentState() const { return currentStateIndex() < stateCount(); }

  // Constructor
  explicit BasisBinaryInputFile(const std::string &fileName);

private:
  void cacheStateHeaderValue();

  double currentStateHeaderValue_;

protected:
  BinaryResultInputFile binFile_;

private:
  // Disallow copy & assignment
  BasisBinaryInputFile(const BasisBinaryInputFile&);
  BasisBinaryInputFile& operator=(const BasisBinaryInputFile&);
};

} /* end namespace Rom */

#endif /* ROM_BASISBINARYFILE_H */
