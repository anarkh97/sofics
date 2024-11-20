#ifndef ROM_VECBASIS_H
#define ROM_VECBASIS_H

#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>
#include <Utils.d/dofset.h>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Eigen/Sparse>
#endif
#include <memory>
#include <algorithm>
#include <vector>

template <typename Scalar> class GenVector;

namespace Rom {

template <typename Scalar, template <typename> class GenVecType>
struct VecTraits {
  typedef GenVecType<Scalar> Type;
  typedef typename Type::InfoType InfoType;
  typedef InfoType InternalInfoType;
  
  static InfoType defaultInfo() { return InfoType(); }
  static int length(InfoType info) { return info; }
  static bool equals(InfoType i, InfoType j) { return i == j; }
  static bool not_equals(InfoType i, InfoType j) { return i != j; }
};

template <typename Scalar, template <typename> class GenVecType = GenVector>
class GenVecBasis : private std::allocator<GenVecType<Scalar> > {
private:
  typedef VecTraits<Scalar, GenVecType> Traits;

public:
  typedef typename Traits::Type VecType;
  typedef typename Traits::InfoType InfoType;

  // Common vector information
  InfoType vectorInfo() const { return vectorInfo_; }

  // Individual vector size (and compatibility aliases)
  int vectorSize() const { return Traits::length(vectorInfo_); }
  int globalVectorSize() const { return vectors_[0].size();}
  int size()       const { return vectorSize(); }
 
  Scalar *data()    const { return buffer_; }

  // Vector count (and compatibility aliases)
  int vectorCount() const { return blockCols_; }
  int numVec()      const { return blockCols_; }
  int numVectors()  const { return blockCols_; }

  // Iteration
  typedef const VecType *const_iterator;
  const_iterator begin() const { return vectors_ ; }
  const_iterator end() const   { return vectors_ + vectorCount_; }
  
  typedef VecType *iterator;
  iterator begin() { return vectors_ ; }
  iterator end()   { return vectors_ + vectorCount_; }

  // Unchecked direct individual vector read access
  const VecType &operator[](int i) const { return vectors_[startCol_+i]; }

  // Unchecked direct individual vector write access
  // (must take care to NOT reallocate underlying memory)
  VecType &operator[](int i) { return vectors_[startCol_+i]; }

  VecType & expand(VecType &, VecType &, bool=true) const;
  VecType & fullExpand(VecType &, VecType &) const;
  VecType & expand(std::vector<Scalar> &, VecType &) const;
  VecType & reduce(VecType &, VecType &, bool=true) const;
  VecType & reduceAll(VecType &, VecType &) const;
  VecType & compressedVecReduce(VecType &, VecType &) const;
  VecType & sparseVecReduce(VecType &, VecType &) const;
  Scalar  & sparseVecReduce(VecType &, Scalar *) const;
  VecType & expand2(VecType &, VecType &) const;
  VecType & addLocalPart(VecType &, VecType &) const;
  
  void makeSparseBasis(const std::vector<std::vector<int> > &, DofSetArray **);
  void makeSparseBasis(const std::vector<std::vector<std::pair<int, int> > > &, DofSetArray **);
  void makeSparseBasis(const std::vector<std::pair<int, int> > &, DofSetArray *);
  void makeSparseBasis(const std::vector<std::map<int,std::vector<int> > > & nodeVec);
  void makeSparseBasis(const std::map<int,std::vector<int> > &nodeVec);
  void makeSparseBasis(const std::vector<int> &, DofSetArray *); 
  void makeSparseBasis(const std::vector<int> & nodeVec, DofSetArray *dsa, DofSetArray *reduced_dsa);
  void makeSparseBasis2(const std::vector<std::vector<int> > &, DofSetArray **);
  void makeSparseBasis2(const std::vector<int> &, DofSetArray *);

  timespec tS1, tS2;
  double time1, time2, time3, time4, time5, time6, counter;

  // Constructors
  GenVecBasis();
  GenVecBasis(int vCount, InfoType vInfo);
 
  // Copy, assignment and swap 
  GenVecBasis(const GenVecBasis &);
  GenVecBasis &operator=(const GenVecBasis &);
  
  // Reshaping
  void dimensionIs(int vCount, InfoType vInfo);

  // Local bases
  void localBasisIs(int startCol, int blockCols);
  int startCol() const { return startCol_; }

  ~GenVecBasis();

private:
  typedef std::allocator<VecType> Allocator;

  void placeVectors();
  void cleanUp();
  void copyBufferContent(const GenVecBasis &other);

  typename Traits::InternalInfoType vectorInfo_;
  int vectorCount_;
  int startCol_, blockCols_; // local bases

  Scalar *buffer_;
  VecType *vectors_;
#ifdef USE_EIGEN3
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> basis_;
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> compressedBasis_, compressedBasis2_;
  std::vector<int> compressedKey_, compressedKey2_;

public:
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& basis() const { return basis_; }
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& compressedBasis() const { return compressedBasis_; }
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& compressedBasis2() const { return compressedBasis2_; }
  const std::vector<int>& compressedKey() const { return compressedKey_; }
#endif
};

template <typename Scalar, template <typename> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::placeVectors() {
#ifdef USE_EIGEN3
  basis_.resize(vectorSize(), vectorCount_);
  buffer_ = basis_.data();
#else
  buffer_ = new Scalar[vectorSize() * vectorCount_];
#endif
  vectors_ = Allocator::allocate(vectorCount_);
  for (int iVec = 0; iVec < vectorCount_; ++iVec) {
    new(vectors_ + iVec) VecType(vectorInfo_, buffer_ + (iVec * vectorSize()), false);
  }
}

template <typename Scalar, template <typename> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::cleanUp() {
  int iVec = vectorCount_;
  while (iVec--) {
    Allocator::destroy(vectors_ + iVec);
  }

  Allocator::deallocate(vectors_, vectorCount_);
#ifndef USE_EIGEN3
  delete[] buffer_;
#endif
}

template <typename Scalar, template <typename> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::copyBufferContent(const GenVecBasis &other) {
  std::copy(other.buffer_, other.buffer_ + vectorSize() * vectorCount_, buffer_);
}

template <typename Scalar, template <typename> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis() :
 vectorInfo_(Traits::defaultInfo()),
 vectorCount_(0)
{
  placeVectors();
}

template <typename Scalar, template <typename> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(int vCount, InfoType vInfo) :
 vectorInfo_(vInfo),
 vectorCount_(vCount),
 startCol_(0),
 blockCols_(vCount)
{
  placeVectors();
}

template <typename Scalar, template <typename> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(const GenVecBasis &other) :
 vectorInfo_(other.vectorInfo_),
 vectorCount_(other.vectorCount_),
 startCol_(other.startCol_),
 blockCols_(other.blockCols_)
{
  placeVectors();
  copyBufferContent(other);
}

template <typename Scalar, template <typename> class GenVecType>
GenVecBasis<Scalar, GenVecType> &
GenVecBasis<Scalar, GenVecType>::operator=(const GenVecBasis &other) {
  if (this != &other) {
    if (vectorCount_ != other.vectorCount_ || Traits::not_equals(vectorInfo_, other.vectorInfo_)) {
      cleanUp();
      vectorCount_ = other.vectorCount_;
      startCol_ = other.startCol_;
      blockCols_ = other.blockCols_;
      typename Traits::InternalInfoType temp(other.vectorInfo_);
      std::swap(vectorInfo_, temp);
      placeVectors();
    }
    copyBufferContent(other);
  }

  return *this;
}

template <typename Scalar, template <typename> class GenVecType>
void
GenVecBasis<Scalar, GenVecType>::dimensionIs(int vCount, InfoType vInfo) {
  if (vCount != vectorCount_ || Traits::not_equals(vInfo, vectorInfo_)) {
    cleanUp();
    vectorCount_ = vCount;
    startCol_ = 0;
    blockCols_ = vCount;
    typename Traits::InternalInfoType temp(vInfo);
    std::swap(vectorInfo_, temp);
    placeVectors();
  }
}

template <typename Scalar, template <typename> class GenVecType>
GenVecBasis<Scalar, GenVecType>::~GenVecBasis() {
  cleanUp();
}

template <typename Scalar, template <typename> class GenVecType>
void
GenVecBasis<Scalar, GenVecType>::localBasisIs(int startCol, int blockCols)
{
  startCol_ = startCol;
  blockCols_ = blockCols;
}

typedef GenVecBasis<double> VecBasis;

} /* end namespace Rom */

#endif /* ROM_VECBASIS_H */
