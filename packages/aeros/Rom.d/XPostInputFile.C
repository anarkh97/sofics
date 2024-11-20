#include "XPostInputFile.h"

#include "NodeDof6Buffer.h"

#include <cstddef>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cassert>

namespace Rom {

XPostInputFile::XPostInputFile(const std::string &fileName) :
  fileName_(fileName),
  nodeCount_(0),
  stateCount_(0),
  stream_(NULL),
  currentStateIndex_(0),
  currentStateHeaderValue_(),
  currentStateRead_(true)
{
  stream_ = std::fopen(fileName_.c_str(), "r");
  
  if (!stream_) {
   throw std::runtime_error("Cannot open input file");
  }

  // Read in number of modes and number of nodes
  char buf[80];
  char *str = fgets(buf, sizeof buf, stream_);
  bool b;
  if(strncmp("Vector", buf, 6) == 0) {
    int count = fscanf(stream_, "%d", &nodeCount_);
    b = (count != 1);
    stateCount_ = 12; // TODO in this case, keep reading until end of file is reached
    std::ifstream fin(fileName_.c_str());
    int linesCount=std::count(std::istreambuf_iterator<char>(fin), 
                              std::istreambuf_iterator<char>(), '\n');
    stateCount_ = (linesCount-2)/(nodeCount_+1);
  }
  else {
    int count1 = sscanf(buf, "%d", &stateCount_);
    int count2 = fscanf(stream_, "%d", &nodeCount_);
    b = ((count1 != 1) || (count2 != 1));
  }

  if(b) {
    throw std::runtime_error("Input file has unrecognized format");
  }

  if (stateCount_ > 0) {
    readCurrentStateHeader();
  }

  readNodeIndices();
}

XPostInputFile::~XPostInputFile() {
  std::fclose(stream_);
}

const NodeDof6Buffer &
XPostInputFile::currentStateBuffer(NodeDof6Buffer &target) const {
  assert(nodeCount() == target.size());

  positionAtStateStart();

  int nodeIndex;
  for (int iNode = 0; iNode != nodeCount(); ++iNode) {
    double *nodeBuffer = target[iNode];
    const int info = std::fscanf(stream_, "%d %le %le %le %le %le %le", &nodeIndex,
                                 &nodeBuffer[0], &nodeBuffer[1], &nodeBuffer[2],
                                 &nodeBuffer[3], &nodeBuffer[4], &nodeBuffer[5]);
    assert(info == 7);
  }

  currentStateRead_ = true;

  return target;
}

void
XPostInputFile::readNodeIndices() {
  nodeIndices_.resize(nodeCount());

  double nodeBuffer[6];
  for (int iNode = 0; iNode != nodeCount(); ++iNode) {
    const int info = std::fscanf(stream_, "%d %le %le %le %le %le %le", &nodeIndices_[iNode],
                                 &nodeBuffer[0], &nodeBuffer[1], &nodeBuffer[2],
                                 &nodeBuffer[3], &nodeBuffer[4], &nodeBuffer[5]);
    nodeIndices_[iNode]--;
    assert(info == 7);
  }
}

void
XPostInputFile::currentStateIndexInc() {
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
XPostInputFile::readCurrentStateHeader() {
  int info;
  info = std::fscanf(stream_, "%le", &currentStateHeaderValue_);
  assert(info == 1);

  info = std::fgetpos(stream_, &currentStatePosition_);
  assert(info == 0);

  currentStateRead_ = false;
}

void
XPostInputFile::positionAtStateStart() const {
  assert(validCurrentState());
  
  if (!currentStateRead_) {
    const int info = std::fsetpos(stream_, &currentStatePosition_);
    assert(info == 0);

    currentStateRead_ = false;
  }
}

} /* end namespace Rom */
