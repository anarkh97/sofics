#include "BasisOutputFile.h"

#include "NodeDof6Buffer.h"
#include "SimpleBuffer.h"

#include <string>
#include <cstdio>

#include <stdexcept>

#include <cassert>

namespace Rom {

const int
BasisOutputFile::STATE_COUNT_LENGTH = 10;

BasisOutputFile::BasisOutputFile(const std::string &fileName, int nodeCount, bool) :
  fileName_(fileName),
  nodeCount_(nodeCount),
  width_(23),
  precision_(15),
  stateCount_(0),
  stream_(NULL),
  stateCountOnFile_(0)
{
  stream_ = std::fopen(fileName_.c_str(), "wt");

  if (!stream_) {
   throw std::runtime_error("Cannot open output file"); 
  }

  writeStateCount();
  writeNodeCount();

  // TODO: Counting iterator [0, nodeCount)
  SimpleBuffer<int> nodeIndices(nodeCount_);
  for (int i = 0; i < nodeCount_; ++i) {
    nodeIndices[i] = i;
  }
  writeIndexMapping(nodeIndices.array(), nodeIndices.array() + nodeCount_);
}

BasisOutputFile::~BasisOutputFile() {
  if (stateCountStatus() != UP_TO_DATE) {
    try {
      rewindAndWriteStateCount();
    } catch (std::runtime_error &) {
      // Log error and swallow exception
      std::fprintf(stderr, "WARNING: Corrupted output file %s\n", fileName_.c_str());
    }
  }

  std::fclose(stream_);
}

void
BasisOutputFile::stateAdd(const NodeDof6Buffer &data, double headValue) {
  assert(nodeCount() == data.size());

  writeStateHeader(headValue);
  
  for (int iNode = 0; iNode < nodeCount(); iNode++)  {
    std::fprintf(stream_, " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                 width_, precision_, data[iNode][0], width_, precision_, data[iNode][1],
                 width_, precision_, data[iNode][2], width_, precision_, data[iNode][3],
                 width_, precision_, data[iNode][4], width_, precision_, data[iNode][5]);
  }

  stateCount_++;
}

void
BasisOutputFile::updateStateCountStatus() {
  if (stateCountStatus() == UP_TO_DATE) {
    return;
  }

  rewindAndWriteStateCount();

  std::fflush(stream_);
  std::fseek(stream_, 0, SEEK_END);
}

void
BasisOutputFile::rewindAndWriteStateCount() {
  std::fflush(stream_);
  std::rewind(stream_);

  writeStateCount();
}

void
BasisOutputFile::writeStateCount() {
  const int actualLength = std::fprintf(stream_, "%-*d", STATE_COUNT_LENGTH, stateCount_);
  if (actualLength != STATE_COUNT_LENGTH) {
    throw std::runtime_error("Corrupted state count");
  }
  std::fprintf(stream_, "\n");
  stateCountOnFile_ = stateCount_;
}

void
BasisOutputFile::writeNodeCount() {
  std::fprintf(stream_, "%d\n", nodeCount_);
}

void
BasisOutputFile::writeStateHeader(double value) {
  std::fprintf(stream_, "  % *.*E\n", width_, precision_, value);
}

} /* end namespace Rom */
