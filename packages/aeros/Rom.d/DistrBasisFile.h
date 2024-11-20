#ifndef ROM_DISTRBASISFILE_H
#define ROM_DISTRBASISFILE_H

#include "BasisBinaryFile.h"

#include "NodeDof6Buffer.h"
#include "DistrNodeDof6Buffer.h"

#include <Comm.d/Communicator.h>

#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <numeric>
#include <memory>
#include <cstddef>

#include <cassert>

namespace Rom {

class DistrBasisOutputFile : public BasisBinaryFile {
public:
  const std::string &fileName() const { return binFile_->pathName(); }

  int nodeCount() const  { return binFile_->itemCount(); }
  int stateCount() const { return binFile_->stateCount(); }

  template<int DOFS_PER_NODE>
  void stateAdd(const DistrNodeDofBuffer<DOFS_PER_NODE> &data);
  template<int DOFS_PER_NODE>
  void stateAdd(const DistrNodeDofBuffer<DOFS_PER_NODE> &data, double headValue);

  // Constructor
  // Collective operation (must be called by all processes sharing the communicator)
  template <typename IdxIt>
  DistrBasisOutputFile(const std::string &fileName, int globalNodeCount, IdxIt localIdxBegin, IdxIt localIdxEnd, Communicator *comm,
                       bool restart, int dofs_per_node = DEFAULT_DOFS_PER_NODE);

private:
  template <typename IdxIt>
  void
  resetHandler(const std::string &fileName, int globalNodeCount, int localOffset, IdxIt localIdxBegin, IdxIt localIdxEnd, bool restart, int dofs_per_node) {
    binFile_.reset(new BinaryResultOutputFile(fileName, NODAL_DATA_FLAG, DESC, globalNodeCount, dofs_per_node, localOffset, localIdxBegin, localIdxEnd, VERSION, restart));
  }

  std::unique_ptr<BinaryResultOutputFile> binFile_;

  // Disallow copy and assigment
  DistrBasisOutputFile(const DistrBasisOutputFile &);
  DistrBasisOutputFile &operator=(const DistrBasisOutputFile &);
};

template <typename Scalar>
Scalar
distr_exclusive_partial_sum(Scalar v, Communicator *comm) {
  const int cpuCount = comm->numCPUs();
  std::vector<Scalar> a(cpuCount, Scalar());
  
  typename std::vector<Scalar>::iterator it(a.begin() + comm->myID());
  *it = v;
  comm->globalSum(a.size(), &a[0]);

  return std::accumulate(a.begin(), it, Scalar());
}

template <typename Scalar, typename Scalar2>
Scalar
asserting_distr_exclusive_partial_sum(Scalar v, Communicator *comm, Scalar2 sum) {
  const int cpuCount = comm->numCPUs();
  std::vector<Scalar> a(cpuCount, Scalar());
  
  typename std::vector<Scalar>::iterator it(a.begin() + comm->myID());
  *it = v;
  comm->globalSum(a.size(), &a[0]);

  assert(std::accumulate(a.begin(), a.end(), Scalar()) == sum);

  return std::accumulate(a.begin(), it, Scalar());
}

template <typename IdxIt>
DistrBasisOutputFile::DistrBasisOutputFile(const std::string &fileName, int globalNodeCount, IdxIt localIdxBegin, IdxIt localIdxEnd, Communicator *comm, bool restart,
                                           int dofs_per_node) :
  binFile_(nullptr)
{
  const int localOffset = asserting_distr_exclusive_partial_sum(std::distance(localIdxBegin, localIdxEnd), comm, globalNodeCount);
  if (localOffset == 0) {
    // Master process goes first, creates/truncates the file without interference
    resetHandler(fileName, globalNodeCount, localOffset, localIdxBegin, localIdxEnd, restart, dofs_per_node);
    comm->sync();
  } else {
    // Other processes wait until master process has created the file
    comm->sync();
    resetHandler(fileName, globalNodeCount, localOffset, localIdxBegin, localIdxEnd, restart, dofs_per_node);
  }
}

template<int DOFS_PER_NODE>
class DistrBasisInputFileTemplate : public BasisBinaryInputFile {
public:
  explicit DistrBasisInputFileTemplate(const std::string &fileName);
  
  const DistrNodeDofBuffer<DOFS_PER_NODE> &currentStateBuffer(DistrNodeDofBuffer<DOFS_PER_NODE> &target);
  void currentStateBuffer(std::map<int,double> &target);
  using BasisBinaryInputFile::currentStateBuffer; // Do not hide inherited member function

private:
  std::map<int, int> fileNodeIds_;
  NodeDofBuffer<DOFS_PER_NODE> fileBuffer_; 
};

typedef DistrBasisInputFileTemplate<6> DistrBasisInputFile;

} /* end namespace Rom */

#endif /* ROM_DISTRBASISFILE_H */
