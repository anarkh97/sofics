#ifndef UTILS_BINARYOUTPUTFILE_H
#define UTILS_BINARYOUTPUTFILE_H

#include <Utils.d/BinFileHandler.h>

#include <string>

// Ouput

void writeHeaderToBinaryOutputFile(BinFileHandler &binFile, int dataType, const std::string &description, int itemCount, int itemDimension);

void writePartialRangeToBinaryOutputFile(BinFileHandler &binFile, const int *indices, int indexCount, int indexOffset, int headerNameBytes);

void writePartialResultToBinaryOutputFile(BinFileHandler &binFile, const double *data, int dataSize, int dataOffset, bool doSerialPart, 
                                          int resultRank, double timeStamp, int inStateDataCount, int clusterItemCount,
                                          int headerNameBytes);

// Input

void readHeaderFromBinaryOutputFile(BinFileHandler &binFile, int &dataType, std::string &description, int &itemCount, int &itemDimension, int &resultCount);

// Indices must point to an array of at least (itemCount) elements
void readRangeFromBinaryOutputFile(BinFileHandler &binFile, int *indices, int itemCount, int headerNameBytes);

void readResultTimeStampFromBinaryOutputFile(BinFileHandler &binFile, int resultRank,
                                             double &timeStamp,
                                             int itemCount, int itemDimension,
                                             int headerNameBytes);

// Data must point to an array of at least (itemCount * dataItemCount) elements
void readPartialResultFromBinaryOutputFile(BinFileHandler &binFile, int resultRank,
                                           double &timeStamp,
                                           double *data, int dataItemFirst, int dataItemCount,
                                           int itemCount, int itemDimension,
                                           int headerNameBytes);

// Data must point to an array of at least (itemCount * dataItemCount) elements
void readPartialResultFromBinaryOutputFile(BinFileHandler &binFile, int resultRank,
                                           double *data, int dataItemFirst, int dataItemCount,
                                           int itemCount, int itemDimension,
                                           int headerNameBytes);

#endif /* UTILS_BINARYOUTPUTFILE_H */
