#ifndef UTILS_BINARYRESULTFILE_H
#define UTILS_BINARYRESULTFILE_H

#include "BinFileHandler.h"
#include "BinaryOutputFile.h"

#include <string>
#include <vector>
#include <iostream>

class BinaryResultOutputFile {
public:
  const std::string &pathName() const { return pathName_; }

  double version() const { return binHandler_.getVersion(); }
  
  // Prelude
  const std::string &description() const { return description_; }
  int dataType() const { return dataType_; }
  int itemCount() const { return itemCount_; }
  int itemDimension() const { return itemDimension_; }
  int stateCount() const { return stateCount_; }
  
  typedef std::vector<int>::const_iterator ItemIdIterator;
  ItemIdIterator itemIdBegin() const { return itemIds_.begin(); }
  ItemIdIterator itemIdEnd() const { return itemIds_.end(); }

  // Data
  int localItemCount() const { return itemIds_.size(); }
  int localDataSize() const { return localItemCount() * itemDimension(); }
  // Writes (localDataSize) x scalars
  void stateAdd(double stamp, const double *data);

  // Constructors
  BinaryResultOutputFile(const std::string &pathName, int dataType, const std::string &description,
                         int itemCount, int itemDimension,
                         double version, bool restart);
  
  template <typename IdxInpIt>
  BinaryResultOutputFile(const std::string &pathName, int dataType, const std::string &description,
                         int itemCount, int itemDimension,
                         int localOffset, IdxInpIt localIdBegin, IdxInpIt localIdEnd,
                         double version, bool restart);

  void writePrelude();

private:
  int stateSize() const { return itemCount() * itemDimension(); }
  bool isMaster() const { return localOffset_ == 0; }

  void updateStateCount();

  std::vector<int>::size_type descriptionSize() const { return description_.size(); }

  std::string pathName_;

  int dataType_;
  std::string description_;
  int itemCount_;
  int itemDimension_;
  int stateCount_;
  
  std::vector<int> itemIds_;
  int localOffset_;

  BinFileHandler binHandler_;

  // Disallow copy and assignment
  BinaryResultOutputFile(const BinaryResultOutputFile &);
  BinaryResultOutputFile &operator=(const BinaryResultOutputFile &);
};

template <typename IdxInpIt>
BinaryResultOutputFile::BinaryResultOutputFile(const std::string &pathName, int dataType, const std::string &description,
                                               int itemCount, int itemDimension,
                                               int localOffset, IdxInpIt localIdBegin, IdxInpIt localIdEnd,
                                               double version, bool restart) :
  pathName_(pathName),
  dataType_(dataType),
  description_(description),
  itemCount_(itemCount),
  itemDimension_(itemDimension),
  stateCount_(0),
  itemIds_(localIdBegin, localIdEnd),
  localOffset_(localOffset),
#if defined(_POSIX_THREAD_SAFE_FUNCTIONS) && defined(_AEROS_ASYNCHRONOUS_IO)
  // need to make sure all of the master thread opens the file first
  // see Rom.d/DistrBasisFile.h
  binHandler_(pathName.c_str(), (isMaster() && !restart) ? "wb" : "rb+", version)
#else
  binHandler_(pathName.c_str(), (isMaster() && !restart) ? "ws" : "ws+", version)
#endif
{
  if(!restart) writePrelude();
  else {
    int rawItemCount;
    BinFileHandler binHandler(pathName.c_str(), "r", version);
    readHeaderFromBinaryOutputFile(binHandler, dataType_, description_, rawItemCount, itemDimension_, stateCount_);
  }
}

class BinaryResultInputFile {
public:
  const std::string &pathName() const { return pathName_; }
 
  double version() const { return binHandler_.getVersion(); }

  // Prelude
  const std::string &description() const { return description_; }
  int dataType() const { return dataType_; }
  int itemCount() const { return itemIds_.size(); }
  int itemDimension() const { return itemDimension_; }
  int stateCount() const { return stateCount_; }
  
  typedef std::vector<int>::const_iterator ItemIdIterator;
  ItemIdIterator itemIdBegin() const { return itemIds_.begin(); }
  ItemIdIterator itemIdEnd() const { return itemIds_.end(); } 

  // Data
  int stateRank() const { return stateRank_; }
  double stateStamp();
  const double* state(double *buffer);

  void stateRankInc();

  // Constructor
  explicit BinaryResultInputFile(const std::string &pathName);

private:
  int stateSize() const { return itemCount() * itemDimension(); }
  
  std::vector<int>::size_type descriptionSize() const { return description_.size(); }
  
  std::string pathName_;

  int dataType_;
  std::string description_;
  int itemDimension_;
  int stateCount_;
  
  std::vector<int> itemIds_;

  BinFileHandler binHandler_;

  int stateRank_;

  // Disallow copy and assignment
  BinaryResultInputFile(const BinaryResultInputFile &);
  BinaryResultInputFile &operator=(const BinaryResultInputFile &);
};

#endif /*UTILS_BINARYRESULTFILE_H */
