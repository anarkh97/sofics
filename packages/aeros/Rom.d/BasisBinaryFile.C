#include "BasisBinaryFile.h"

#include <sys/types.h>

#include <stdexcept>
#include <cassert>
#include <iostream>

namespace Rom {

const double BasisBinaryFile::VERSION = 2.0;
const std::string BasisBinaryFile::DESC = "rob";
const int BasisBinaryFile::NODAL_DATA_FLAG = 1;
const int BasisBinaryFile::DEFAULT_DOFS_PER_NODE = 6;


BasisBinaryOutputFile::BasisBinaryOutputFile(const std::string &fileName, int nodeCount, bool restart,
                                             int dofs_per_node) :
  binFile_(fileName, NODAL_DATA_FLAG, DESC, nodeCount, dofs_per_node, VERSION, restart)
{}

template<int DOFS_PER_NODE>
void
BasisBinaryOutputFile::stateAdd(const NodeDofBuffer<DOFS_PER_NODE> &data, double headValue) {
  assert(nodeCount() == data.size());

  // Dump all information in one pass
  binFile_.stateAdd(headValue, data.array());
}


BasisBinaryInputFile::BasisBinaryInputFile(const std::string &fileName) :
  binFile_(fileName)
{
  if (binFile_.version() != VERSION) {
    throw std::runtime_error("Incompatible binary file version");
  }

  if (binFile_.dataType() != NODAL_DATA_FLAG) {
    throw std::runtime_error("Non-nodal data");
  }
 
  if (binFile_.description() != DESC) {
    throw std::runtime_error("Incorrect description");
  }
 
  assert(currentStateIndex() == 0);

  cacheStateHeaderValue();
}

template<int DOFS_PER_NODE>
const NodeDofBuffer<DOFS_PER_NODE> &
BasisBinaryInputFile::currentStateBuffer(NodeDofBuffer<DOFS_PER_NODE> &target) {
  assert(validCurrentState());
  assert(nodeCount() == target.size());
  
  // Retrieve all information in one pass
  binFile_.state(target.array());

  return target;
}

void
BasisBinaryInputFile::currentStateIndexInc() {
  if (validCurrentState()) {
    binFile_.stateRankInc();
    cacheStateHeaderValue();
  }
}

void BasisBinaryInputFile::cacheStateHeaderValue() {
  if (validCurrentState()) {
    currentStateHeaderValue_ = binFile_.stateStamp();
  }
}

template
void
BasisBinaryOutputFile::stateAdd(const NodeDofBuffer<6> &data, double headValue);

template
const NodeDofBuffer<6> &
BasisBinaryInputFile::currentStateBuffer(NodeDofBuffer<6> &target);

template
void
BasisBinaryOutputFile::stateAdd(const NodeDofBuffer<1> &data, double headValue);

template
const NodeDofBuffer<1> &
BasisBinaryInputFile::currentStateBuffer(NodeDofBuffer<1> &target);

} /* end namespace Rom */
