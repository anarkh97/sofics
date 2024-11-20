#include "BinaryOutputFile.h"

#include <cassert>
#include <vector>
#include <iostream>

// Auxiliary functions

inline
BinFileHandler::OffType headerBytes(BinFileHandler::OffType headerNameBytes) {
  return headerNameBytes + 5 * sizeof(int);
}

// Output functions

void writeHeaderToBinaryOutputFile(BinFileHandler &binFile, int dataType, const std::string &description, int itemCount, int itemDimension)
{
  // Data type: int x 1
  binFile.write(&dataType, 1);
  // Header length: int x 1
  const int descriptionLength = description.size();
  binFile.write(&descriptionLength, 1);
  // Header description: char x (header length)
  binFile.write(description.data(), descriptionLength);

  // Number of items in data of cluster: int x 1
  binFile.write(&itemCount, 1);
  // Dimension of data: int x 1
  binFile.write(&itemDimension, 1);

  // Number of results: int x 1
  const int dummy = 0;
  binFile.write(&dummy, 1);
}

void writePartialRangeToBinaryOutputFile(BinFileHandler &binFile, const int *indices, int indexCount, int indexOffset, int headerNameBytes)
{
  const BinFileHandler::OffType beginOffset = headerBytes(headerNameBytes) + indexOffset * sizeof(int);
  binFile.seek(beginOffset);

  // Range: int x (nData) 
  binFile.write(indices, indexCount);
}

void writePartialResultToBinaryOutputFile(BinFileHandler &binFile, const double *data, int dataSize, int dataOffset, bool doSerialPart, 
                                          int resultRank, double timeStamp, int inStateDataCount, int clusterItemCount,
                                          int headerNameBytes) {
  const BinFileHandler::OffType outputRangeBytes = clusterItemCount * sizeof(int);
  const BinFileHandler::OffType firstResultOffset = headerBytes(headerNameBytes) + outputRangeBytes;
  const BinFileHandler::OffType stateBytes = (inStateDataCount + 1) * sizeof(double);

  const BinFileHandler::OffType resultOffset = firstResultOffset + stateBytes * (resultRank - 1);

  //assert(doSerialPart == (dataOffset == 0));
  if (doSerialPart) { // Must be done only once
    // Update number of results
    const BinFileHandler::OffType resultCountOffset = headerNameBytes + 4 * sizeof(int);
    binFile.seek(resultCountOffset);
    binFile.write(&resultRank, 1);
    // Timestamp: double x 1
    binFile.seek(resultOffset);
    binFile.write(&timeStamp, 1);
  }

  // Data: double x (totLen)
  // Write only part attributed to subdomain
  const BinFileHandler::OffType resultDataOffset = resultOffset + (dataOffset + 1) * sizeof(double);
  binFile.seek(resultDataOffset);
  binFile.write(data, dataSize);
}


// Input functions

void readHeaderFromBinaryOutputFile(BinFileHandler &binFile, int &dataType, std::string &description, int &itemCount, int &itemDimension, int &resultCount) {
  // Rewind to beginning of file
  binFile.seek(BinFileHandler::OffType(0));
  
  // Data type: int x 1
  binFile.read(&dataType, 1);
  // Header length: int x 1
  int descriptionLength;
  binFile.read(&descriptionLength, 1);
  // Header description: char x (header length)
  std::vector<char> rawDescription(descriptionLength);
  if (descriptionLength) {
    binFile.read(&rawDescription[0], descriptionLength);
  }
  description.assign(rawDescription.begin(), rawDescription.end());

  // Number of items in data of cluster: int x 1
  binFile.read(&itemCount, 1);
  // Dimension of data: int x 1
  binFile.read(&itemDimension, 1);

  // Number of results: int x 1
  binFile.read(&resultCount, 1);
}

void readRangeFromBinaryOutputFile(BinFileHandler &binFile, int *indices, int itemCount, int headerNameBytes) {
  const BinFileHandler::OffType beginOffset = headerBytes(headerNameBytes);
  binFile.seek(beginOffset);

  // Range: int x (nData)
  binFile.read(indices, itemCount);
}

inline
BinFileHandler::OffType resultStampOffset(int resultRank, int itemCount, int itemDimension, int headerNameBytes) {
  const BinFileHandler::OffType outputRangeBytes = itemCount * sizeof(int);
  const BinFileHandler::OffType firstResultOffset = headerBytes(headerNameBytes) + outputRangeBytes;
  const BinFileHandler::OffType stateBytes = (itemCount * itemDimension + 1) * sizeof(double);

  return firstResultOffset + stateBytes * (resultRank - 1);
}

void readResultTimeStampFromBinaryOutputFile(BinFileHandler &binFile, int resultRank,
                                             double &timeStamp,
                                             int itemCount, int itemDimension,
                                             int headerNameBytes) {
  const BinFileHandler::OffType resultOffset = resultStampOffset(resultRank, itemCount, itemDimension, headerNameBytes);
  
  // Timestamp: double x 1
  binFile.seek(resultOffset);
  binFile.read(&timeStamp, 1);
}

void readPartialResultFromBinaryOutputFile(BinFileHandler &binFile, int resultRank,
                                           double &timeStamp,
                                           double *data, int dataItemFirst, int dataItemCount,
                                           int itemCount, int itemDimension,
                                           int headerNameBytes) {
  const BinFileHandler::OffType resultOffset = resultStampOffset(resultRank, itemCount, itemDimension, headerNameBytes);
  
  // Timestamp: double x 1
  binFile.seek(resultOffset);
  binFile.read(&timeStamp, 1);

  readPartialResultFromBinaryOutputFile(binFile, resultRank,
                                        data, dataItemFirst, dataItemCount,
                                        itemCount, itemDimension,
                                        headerNameBytes);
}

void readPartialResultFromBinaryOutputFile(BinFileHandler &binFile, int resultRank,
                                           double *data, int dataItemFirst, int dataItemCount,
                                           int itemCount, int itemDimension,
                                           int headerNameBytes) {
  const BinFileHandler::OffType resultOffset = resultStampOffset(resultRank, itemCount, itemDimension, headerNameBytes);
  const BinFileHandler::OffType resultDataOffset = resultOffset + (dataItemFirst * itemDimension + 1) * sizeof(double);
 
  // Data: double x (dataItemCount * itemDimension) 
  binFile.seek(resultDataOffset);
  binFile.read(data, dataItemCount * itemDimension);
}

