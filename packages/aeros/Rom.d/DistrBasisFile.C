#include "DistrBasisFile.h"

#include "RenumberingUtils.h"

#include <algorithm>
#include <iterator>
#include <stdexcept>

#include <cassert>

namespace Rom {

template<int DOFS_PER_NODE>
void
DistrBasisOutputFile::stateAdd(const DistrNodeDofBuffer<DOFS_PER_NODE> &data) {
  const double defaultHeader = static_cast<double>(stateCount() + 1);
  stateAdd(data,defaultHeader);
}

template<int DOFS_PER_NODE>
void
DistrBasisOutputFile::stateAdd(const DistrNodeDofBuffer<DOFS_PER_NODE> &data, double headValue) {
  assert(data.localNodeCount() == binFile_->localItemCount());
  SimpleBuffer<double> buffer(binFile_->localDataSize());

  double *target = buffer.array();
  for (BinaryResultOutputFile::ItemIdIterator it     = binFile_->itemIdBegin(),
                                              it_end = binFile_->itemIdEnd();
                                              it    != it_end;
                                              ++it) {
    const double *origin = data[*it];
    std::copy(origin, origin + DOFS_PER_NODE, target);
    target += DOFS_PER_NODE;
  }

  binFile_->stateAdd(headValue, buffer.array());
}


template<int DOFS_PER_NODE>
DistrBasisInputFileTemplate<DOFS_PER_NODE>::DistrBasisInputFileTemplate(const std::string &fileName) :
  BasisBinaryInputFile(fileName),
  fileNodeIds_(),
  fileBuffer_(nodeCount())
{
  inverse_numbering(nodeIdBegin(), nodeIdEnd(), std::inserter(fileNodeIds_, fileNodeIds_.end()));
  assert(fileNodeIds_.size() <= nodeCount());
}

template<int DOFS_PER_NODE>
const DistrNodeDofBuffer<DOFS_PER_NODE> &
DistrBasisInputFileTemplate<DOFS_PER_NODE>::currentStateBuffer(DistrNodeDofBuffer<DOFS_PER_NODE> &target) {
  // Retrieve all information in the internal buffer
  // TODO: More economical approach
  BasisBinaryInputFile::currentStateBuffer(fileBuffer_);
  
  typedef typename DistrNodeDofBuffer<DOFS_PER_NODE>::NodeItConst NodeIt;
  const NodeIt it_end = target.globalNodeIndexEnd();

  for (NodeIt it = target.globalNodeIndexBegin(); it != it_end; ++it) {
    const int iNode = *it; // Id of the node requested by the local buffer

    // Location of requested node in internal buffer
    std::map<int, int>::const_iterator it_loc = fileNodeIds_.find(iNode);
    if (it_loc == fileNodeIds_.end()) {
      throw std::runtime_error("Requested nodal data missing from file");
    }
    const double *nodeBuffer = fileBuffer_[it_loc->second];

    std::copy(nodeBuffer, nodeBuffer + DOFS_PER_NODE, &target[iNode][0]);
  }

  return target;
}

template<int DOFS_PER_NODE>
void
DistrBasisInputFileTemplate<DOFS_PER_NODE>::currentStateBuffer(std::map<int,double> &target) {
  BasisBinaryInputFile::currentStateBuffer(fileBuffer_);

  for(NodeIdIterator it = binFile_.itemIdBegin(); it != binFile_.itemIdEnd(); ++it){
    double *muVal = fileBuffer_[*it];
    if(*muVal > 1e-16 || *muVal < -1e-16){
      target.insert(std::make_pair(*it,*muVal));
    }
  }
}

template
void DistrBasisOutputFile::stateAdd(const DistrNodeDofBuffer<6> &data);

template
void DistrBasisOutputFile::stateAdd(const DistrNodeDofBuffer<6> &data, double headValue);

template
DistrBasisInputFileTemplate<6>::DistrBasisInputFileTemplate(const std::string &fileName);

template
const DistrNodeDofBuffer<6> &
DistrBasisInputFileTemplate<6>::currentStateBuffer(DistrNodeDofBuffer<6> &target);

template
void DistrBasisOutputFile::stateAdd(const DistrNodeDofBuffer<1> &data);

template
void DistrBasisOutputFile::stateAdd(const DistrNodeDofBuffer<1> &data, double headValue);

template
DistrBasisInputFileTemplate<1>::DistrBasisInputFileTemplate(const std::string &fileName);

template
const DistrNodeDofBuffer<1> &
DistrBasisInputFileTemplate<1>::currentStateBuffer(DistrNodeDofBuffer<1> &target);

template
void
DistrBasisInputFileTemplate<1>::currentStateBuffer(std::map<int,double> &target);

} /* end namespace Rom */
