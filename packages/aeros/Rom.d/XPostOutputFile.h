#ifndef ROM_XPOSTOUTPUTFILE_H
#define ROM_XPOSTOUTPUTFILE_H

#include <string>
#include <cstdio>
#include <stdexcept>
#include <iterator>

namespace Rom {

template <int> class NodeDofBuffer;
typedef NodeDofBuffer<6> NodeDof6Buffer;

template <int NumColumns=3>
class XPostOutputFile {
public:
  const std::string &fileName() const { return fileName_; }

  int nodeCount() const  { return nodeCount_;  }
  int stateCount() const { return stateCount_; }

  // Node ordering determined by NodeDof6Buffer, file node ordering ignored
  void stateAdd(const NodeDof6Buffer &data);
  void stateAdd(const NodeDof6Buffer &data, double headValue);
  
  XPostOutputFile(const std::string &fileName, int nodeCount, bool);

  template <typename NodeIdIt>
  XPostOutputFile(const std::string &fileName, NodeIdIt first, NodeIdIt last, bool);
  
  ~XPostOutputFile();

private:
  const std::string fileName_;
  const int nodeCount_;
  const int width_;
  const int precision_;

  int stateCount_;
  
  FILE *stream_;

  void writeHeader();
  void writeNodeCount();
  void writeStateHeader(double value);

  // Disallow copy & assignment
  XPostOutputFile(const XPostOutputFile&);
  XPostOutputFile& operator=(const XPostOutputFile&);
};

template <int NumColumns>
template <typename NodeIdIt>
XPostOutputFile<NumColumns>::XPostOutputFile(const std::string &fileName, NodeIdIt first, NodeIdIt last, bool) :
  fileName_(fileName),
  nodeCount_(std::distance(first, last)),
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
}

template<int NumColumns>
inline
void
XPostOutputFile<NumColumns>::stateAdd(const NodeDof6Buffer &data) {
  stateAdd(data, static_cast<double>(stateCount() + 1));
}

} /* end namespace Rom */

#ifdef _TEMPLATE_FIX_
#include <Rom.d/XPostOutputFile.C>
#endif

#endif /* ROM_XPOSTOUTPUTFILE_H */
