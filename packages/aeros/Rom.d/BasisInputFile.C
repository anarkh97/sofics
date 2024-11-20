#include "BasisInputFile.h"

#include "NodeDof6Buffer.h"

#include <cstddef>
#include <stdexcept>

#include <cassert>

namespace Rom {

BasisInputFile::BasisInputFile(const std::string &fileName) :
  fileName_(fileName),
  nodeCount_(0),
  stateCount_(0),
  stream_(NULL),
  currentStateIndex_(0),
  currentStateHeaderValue_(),
  currentStateRead_(true)
{
  stream_ = std::fopen(fileName_.c_str(), "rt");
  
  if (!stream_) {
   throw std::runtime_error("Cannot open input file");
  }

  int info = std::fscanf(stream_, "%d %d", &stateCount_, &nodeCount_);
  
  if (info != 2) {
   throw std::runtime_error("Input file has unrecognized format");
  }

  readNodeIndices();

  if (stateCount_ > 0) {
    readCurrentStateHeader();
  }
}

BasisInputFile::~BasisInputFile() {
  std::fclose(stream_);
}

const NodeDof6Buffer &
BasisInputFile::currentStateBuffer(NodeDof6Buffer &target) const {
  assert(nodeCount() == target.size());

  positionAtStateStart();

  for (int iNode = 0; iNode != nodeCount(); ++iNode) {
    double *nodeBuffer = target[iNode];
    const int info = std::fscanf(stream_, "%le %le %le %le %le %le",
                                 &nodeBuffer[0], &nodeBuffer[1], &nodeBuffer[2],
                                 &nodeBuffer[3], &nodeBuffer[4], &nodeBuffer[5]);
    assert(info == 6);
  }

  currentStateRead_ = true;

  return target;
}

void
BasisInputFile::readNodeIndices() {
  nodeIndices_.resize(nodeCount());
  for (int iNode = 0; iNode != nodeCount(); ++iNode) { 
    const int info = std::fscanf(stream_, "%d", &nodeIndices_[iNode]);
    assert(info == 1);
    nodeIndices_[iNode] -= 1;
  }
}

void
BasisInputFile::currentStateIndexInc() {
  assert(validCurrentState());
  assert(currentStateRead_);

  if (currentStateIndex() + 1 < stateCount()) {
    readCurrentStateHeader();
  } else {
    currentStateHeaderValue_ = double();
    currentStateRead_ = false;
  }
  ++currentStateIndex_;
}

void
BasisInputFile::readCurrentStateHeader() {
  int info;
  info = std::fscanf(stream_, "%le", &currentStateHeaderValue_);
  assert(info == 1);

  info = std::fgetpos(stream_, &currentStatePosition_);
  assert(info == 0);

  currentStateRead_ = false;
}

void
BasisInputFile::positionAtStateStart() const {
  assert(validCurrentState());
  
  if (!currentStateRead_) {
    const int info = std::fsetpos(stream_, &currentStatePosition_);
    assert(info == 0);

    currentStateRead_ = false;
  }
}

} /* end namespace Rom */
