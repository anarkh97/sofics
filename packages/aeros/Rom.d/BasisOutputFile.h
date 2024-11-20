#ifndef ROM_BASISOUTPUTFILE_H
#define ROM_BASISOUTPUTFILE_H

#include <string>
#include <cstdio>
#include <stdexcept>
#include <iterator>

namespace Rom {

template <int> class NodeDofBuffer;
typedef NodeDofBuffer<6> NodeDof6Buffer;

class BasisOutputFile {
public:
  const std::string &fileName() const { return fileName_; }

  int nodeCount() const  { return nodeCount_;  }
  int stateCount() const { return stateCount_; }

  enum StateCountStatus { UP_TO_DATE, OUTDATED };
  StateCountStatus stateCountStatus() const;
  void updateStateCountStatus();

  // Node ordering determined by NodeDof6Buffer, file node ordering ignored
  void stateAdd(const NodeDof6Buffer &data);
  void stateAdd(const NodeDof6Buffer &data, double headValue);
  
  BasisOutputFile(const std::string &fileName, int nodeCount, bool);
  
  template <typename NodeIdIt>
  BasisOutputFile(const std::string &fileName, NodeIdIt first, NodeIdIt last, bool);
 
  ~BasisOutputFile();

private:
  const std::string fileName_;
  const int nodeCount_;
  const int width_;
  const int precision_;

  int stateCount_;
  
  FILE *stream_;

  static const int STATE_COUNT_LENGTH;
  int stateCountOnFile_;
  void writeStateCount();
  void rewindAndWriteStateCount(); 

  void writeNodeCount();
  template <typename NodeIdIt>
  void writeIndexMapping(NodeIdIt first, NodeIdIt last);

  void writeStateHeader(double value);

  // Disallow copy & assignment
  BasisOutputFile(const BasisOutputFile&);
  BasisOutputFile& operator=(const BasisOutputFile&);
};

template <typename NodeIdIt>
BasisOutputFile::BasisOutputFile(const std::string &fileName, NodeIdIt first, NodeIdIt last, bool) :
  fileName_(fileName),
  nodeCount_(std::distance(first, last)),
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
  writeIndexMapping(first, last);
}

template <typename NodeIdIt>
void
BasisOutputFile::writeIndexMapping(NodeIdIt first, NodeIdIt last) {
  for (NodeIdIt it = first; it != last; ++it) {
    std::fprintf(stream_, " %d\n", *it + 1);
  }
}

inline
BasisOutputFile::StateCountStatus
BasisOutputFile::stateCountStatus() const {
  return (stateCount() == stateCountOnFile_) ? UP_TO_DATE : OUTDATED;
}

inline
void
BasisOutputFile::stateAdd(const NodeDof6Buffer &data) {
  stateAdd(data, static_cast<double>(stateCount() + 1));
}

} /* end namespace Rom */

#endif /* ROM_BASISOUTPUTFILE_H */
