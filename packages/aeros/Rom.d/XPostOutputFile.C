#include "XPostOutputFile.h"

#include "NodeDof6Buffer.h"
#include "SimpleBuffer.h"

#include <string>
#include <cstdio>

#include <stdexcept>

#include <cassert>

namespace Rom {

template<int NumColumns>
XPostOutputFile<NumColumns>::XPostOutputFile(const std::string &fileName, int nodeCount, bool) :
  fileName_(fileName),
  nodeCount_(nodeCount),
  width_(22),
  precision_(15),
  stateCount_(0),
  stream_(NULL)
{
  stream_ = std::fopen(fileName_.c_str(), "wt");

  if (!stream_) {
   throw std::runtime_error("Cannot open output file"); 
  }

  writeHeader();
  writeNodeCount();

  // TODO: Counting iterator [0, nodeCount)
  SimpleBuffer<int> nodeIndices(nodeCount_);
  for (int i = 0; i < nodeCount_; ++i) {
    nodeIndices[i] = i;
  }
}

template<int NumColumns>
XPostOutputFile<NumColumns>::~XPostOutputFile() {
  std::fclose(stream_);
}

template<int NumColumns>
void
XPostOutputFile<NumColumns>::stateAdd(const NodeDof6Buffer &data, double headValue) {
  assert(nodeCount() == data.size());

  writeStateHeader(headValue);

  switch(NumColumns) {
    default:
    case 3:
      for (int iNode = 0; iNode < nodeCount(); iNode++) {
        std::fprintf(stream_, " % *.*E % *.*E % *.*E\n",
                     width_, precision_, data[iNode][0], width_, precision_, data[iNode][1],
                     width_, precision_, data[iNode][2]);
      }
      break;
    case 6:
      for (int iNode = 0; iNode < nodeCount(); iNode++) {
        std::fprintf(stream_, " %d % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n", iNode+1,
                     width_, precision_, data[iNode][0], width_, precision_, data[iNode][1],
                     width_, precision_, data[iNode][2], width_, precision_, data[iNode][3],
                     width_, precision_, data[iNode][4], width_, precision_, data[iNode][5]);
      }
      break;
  }

  stateCount_++;
}

template<int NumColumns>
void
XPostOutputFile<NumColumns>::writeHeader() {
  switch(NumColumns) {
    default:
    case 3:
      std::fprintf(stream_, "Vector MODE under Modal for nodeset\n");
      break;
    case 6:
      std::fprintf(stream_, "Vector MODE6 under Modal for nodeset\n");
      break;
  }
}

template<int NumColumns>
void
XPostOutputFile<NumColumns>::writeNodeCount() {
  std::fprintf(stream_, "%d\n", nodeCount_);
}

template<int NumColumns>
void
XPostOutputFile<NumColumns>::writeStateHeader(double value) {
  std::fprintf(stream_, "  % *.*E\n", width_, precision_, value);
}

} /* end namespace Rom */
