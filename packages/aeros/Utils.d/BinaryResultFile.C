#include "BinaryResultFile.h"

#include <cstddef>
#include <cassert>

BinaryResultOutputFile::BinaryResultOutputFile(const std::string &pathName, int dataType, const std::string &description, int itemCount, int itemDimension, double version, bool restart) :
  pathName_(pathName),
  dataType_(dataType),
  description_(description),
  itemCount_(itemCount),
  itemDimension_(itemDimension),
  stateCount_(0),
  itemIds_(itemCount),
  localOffset_(0),
  binHandler_(pathName.c_str(), (!restart) ? "w" : "w+", version)
{
  int iNode = 0;
  for (std::vector<int>::iterator it = itemIds_.begin(), it_end = itemIds_.end(); it != it_end; ++it) {
    *it = iNode++;
  }

  if(!restart) writePrelude();
  else {
    int rawItemCount;
    BinFileHandler binHandler(pathName.c_str(), "r", version);
    readHeaderFromBinaryOutputFile(binHandler, dataType_, description_, rawItemCount, itemDimension_, stateCount_);
  }
}

void
BinaryResultOutputFile::writePrelude() {
  if (isMaster()) {
    writeHeaderToBinaryOutputFile(binHandler_, dataType(), description(), itemCount(), itemDimension());
  }
  if (!itemIds_.empty()) {
    writePartialRangeToBinaryOutputFile(binHandler_, &itemIds_[0], localItemCount(), localOffset_, descriptionSize());
  }
}

void
BinaryResultOutputFile::stateAdd(double stamp, const double *buffer) {
  writePartialResultToBinaryOutputFile(binHandler_,
                                       buffer, localDataSize(), localOffset_ * itemDimension(),
                                       isMaster(),
                                       stateCount() + 1, stamp, stateSize(), itemCount(),
                                       descriptionSize());
  ++stateCount_;
  updateStateCount();
}

void
BinaryResultOutputFile::updateStateCount() {
  // Already up-to-date, nothing to do
}


BinaryResultInputFile::BinaryResultInputFile(const std::string &pathName) :
  pathName_(pathName),
  binHandler_(pathName.c_str(), "r"),
  stateRank_(0)
{
  int rawItemCount;
  readHeaderFromBinaryOutputFile(binHandler_, dataType_, description_, rawItemCount, itemDimension_, stateCount_);

  itemIds_.resize(rawItemCount);
  if (rawItemCount != 0) {
    readRangeFromBinaryOutputFile(binHandler_, &itemIds_[0], itemIds_.size(), descriptionSize());
  }
}

void
BinaryResultInputFile::stateRankInc() {
  assert(stateRank() < stateCount()); 
  
  ++stateRank_;
}

double
BinaryResultInputFile::stateStamp() {
  assert(stateRank() < stateCount()); 
  
  double result;
  readResultTimeStampFromBinaryOutputFile(binHandler_, stateRank_ + 1, result, itemCount(), itemDimension(), descriptionSize());
  return result;
}

const double *
BinaryResultInputFile::state(double *buffer) {
  assert(stateRank() < stateCount()); 
  
  readPartialResultFromBinaryOutputFile(binHandler_, stateRank_ + 1,
                                        buffer, 0, itemCount(),
                                        itemCount(), itemDimension(),
                                        descriptionSize());
  return buffer;
}

