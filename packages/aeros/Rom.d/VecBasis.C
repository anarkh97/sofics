#include "VecBasis.h"
#include "Utils.d/DistHelper.h"
#include "Math.d/Vector.h"

#include <set>
namespace Rom {

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::expand(GenDistrVector<double> &x, GenDistrVector<double> &_result,
                                            bool useCompressedBasis) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data()+startCol_, blockCols_);
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
    result.setZero();

    if(useCompressedBasis && compressedKey_.size() > 0) {
      Eigen::VectorXd resultBuffer(compressedKey_.size());
      resultBuffer = compressedBasis_*GenCoordinates;
      for(int i = 0; i < compressedKey_.size(); i++)
        result(compressedKey_[i]) = resultBuffer(i);
    }
    else {
      result = basis_.block(0,startCol_,basis_.rows(),blockCols_)*GenCoordinates;
    }
  }
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::fullExpand(GenDistrVector<double> &x, GenDistrVector<double> &_result) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data(), x.size());
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
    result = basis_*GenCoordinates;
  }
#endif
  return _result;
}

template <>
GenVector<double> &
GenVecBasis<double, GenVector>::expand(GenVector<double> &x, GenVector<double> &_result,
                                       bool useCompressedBasis) const {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data()+startCol_, blockCols_);
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());

  if(useCompressedBasis && compressedKey_.size() > 0) {
    result.setZero();
    Eigen::VectorXd resultBuffer(compressedKey_.size());
    resultBuffer = compressedBasis_*GenCoordinates;
    for(int i = 0; i < compressedKey_.size(); i++)
      result(compressedKey_[i]) = resultBuffer(i);
  }
  else {
    result = basis_.block(0,startCol_,basis_.rows(),blockCols_)*GenCoordinates;
  }
#endif
  return _result;
}

template <>
GenVector<double> &
GenVecBasis<double, GenVector>::fullExpand(GenVector<double> &x, GenVector<double> &_result) const {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());

  result = basis_*GenCoordinates;
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::expand2(GenDistrVector<double> &x, GenDistrVector<double> &_result) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data()+startCol_, blockCols_);
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
    result.setZero();

    if(compressedKey2_.size() > 0) {
      Eigen::VectorXd resultBuffer(compressedKey2_.size());
      resultBuffer = compressedBasis2_.block(0,startCol_,compressedBasis2_.rows(),blockCols_)*GenCoordinates;
      for(int i = 0; i < compressedKey2_.size(); i++)
        result(compressedKey2_[i]) = resultBuffer(i);
    }
    else {
      result = basis_.block(0,startCol_,basis_.rows(),blockCols_)*GenCoordinates;
    }
  }
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::expand(std::vector<double> &x, GenDistrVector<double> &_result) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data()+startCol_, blockCols_);
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());

    result = basis_.block(0,startCol_,basis_.rows(),blockCols_)*GenCoordinates;
  }
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::compressedVecReduce(GenDistrVector<double> &x, GenDistrVector<double> &_result) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data()+startCol_, blockCols_);
    result.setZero();

    Eigen::VectorXd coordBuffer(compressedKey_.size());
    for(int i = 0; i < compressedKey_.size(); i++)
      coordBuffer(i) = FullCoordinates(compressedKey_[i]); 

    result = compressedBasis_.transpose()*coordBuffer;
    //each process gets a copy of reduced coordinates
    if(structCom)
      structCom->globalSum(result.size(), result.data());
  }
  else if(structCom) {
    Eigen::Matrix<double,Eigen::Dynamic, 1> result(compressedBasis_.cols());
    result.setZero();
    structCom->globalSum(result.size(), result.data());
  }
#endif 
  return _result;
}

template <>
GenVector<double> &
GenVecBasis<double, GenVector>::compressedVecReduce(GenVector<double> &x, GenVector<double> &_result) const {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data()+startCol_, blockCols_);
  result.setZero();

  Eigen::VectorXd coordBuffer(compressedKey_.size());
  for(int i = 0; i < compressedKey_.size(); i++)
    coordBuffer(i) = FullCoordinates(compressedKey_[i]);

  result = compressedBasis_.transpose()*coordBuffer;
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::sparseVecReduce(GenDistrVector<double> &x, GenDistrVector<double> &_result) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data()+startCol_, blockCols_);
    Eigen::SparseVector<double> sparsef(FullCoordinates.rows());
    result.setZero();

    for(int i = 0; i < FullCoordinates.rows(); i++)
      if(FullCoordinates(i) != 0)
        sparsef.insert(i) = FullCoordinates(i);

    result = basis_.block(0,startCol_,basis_.rows(),blockCols_).transpose()*sparsef;
    //each process gets a copy of reduced coordinates
    if(structCom)
      structCom->globalSum(result.size(), result.data());
  }
  else if(structCom) {
    Eigen::Matrix<double,Eigen::Dynamic, 1> result(blockCols_);
    result.setZero();
    structCom->globalSum(result.size(), result.data());
  }
#endif
  return _result;
}

template <>
double &
GenVecBasis<double, GenDistrVector>::sparseVecReduce(GenDistrVector<double> &x, double *_result) const {
#ifdef USE_EIGEN3
  if(_result) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result, blockCols_);
    Eigen::SparseVector<double> sparsef(FullCoordinates.rows());
    result.setZero();

    for(int i = 0; i < FullCoordinates.rows(); i++)
      if(FullCoordinates(i) != 0)
        sparsef.insert(i) = FullCoordinates(i);

    result = basis_.block(0,startCol_,basis_.rows(),blockCols_).transpose()*sparsef;
    //each process gets a copy of reduced coordinates
    if(structCom)
      structCom->globalSum(result.size(), result.data());
  }
  else if(structCom) {
    Eigen::Matrix<double,Eigen::Dynamic, 1> result(blockCols_);
    result.setZero();
    structCom->globalSum(result.size(), result.data());
  }
#endif
  return *_result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::reduce(GenDistrVector<double> &x, GenDistrVector<double> &_result, bool useCompressedBasis) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) { 
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data()+startCol_, blockCols_);
    result.setZero();

    if(useCompressedBasis && compressedKey_.size() > 0) {
      /*Eigen::VectorXd coordBuffer(compressedKey_.size());
      for(int i = 0; i < compressedKey_.size(); i++)
        coordBuffer(i) = FullCoordinates(compressedKey_[i]);

      result = compressedBasis_.transpose()*coordBuffer;*/
      sparseVecReduce(x, _result);
    }
    else {
      result = basis_.block(0,startCol_,basis_.rows(),blockCols_).transpose()*FullCoordinates;
      //each process gets a copy of reduced coordinates
      if(structCom)
        structCom->globalSum(result.size(), result.data());
    }
  }
  else if(structCom) {
    Eigen::Matrix<double,Eigen::Dynamic, 1> result(blockCols_);
    result.setZero();
    structCom->globalSum(result.size(), result.data());
  }
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::reduceAll(GenDistrVector<double> &x, GenDistrVector<double> &_result) const {
#ifdef USE_EIGEN3
  if(_result.size() > 0) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
    result = basis_.transpose()*FullCoordinates;

    //each process gets a copy of reduced coordinates
    if(structCom)
      structCom->globalSum(result.size(), result.data());
  }
  else if(structCom) {
    Eigen::Matrix<double,Eigen::Dynamic, 1> result(basis_.cols());
    result.setZero();
    structCom->globalSum(result.size(), result.data());
  }
#endif
  return _result;
}

template <>
GenVector<double> &
GenVecBasis<double, GenVector>::reduce(GenVector<double> &x, GenVector<double> &_result, bool useCompressedBasis) const {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data()+startCol_, blockCols_);

  result.setZero();
  if(useCompressedBasis && compressedKey_.size() > 0) {
    Eigen::VectorXd coordBuffer(compressedKey_.size());
    for(int i = 0; i < compressedKey_.size(); i++)
      coordBuffer(i) = FullCoordinates(compressedKey_[i]);

    result = compressedBasis_.transpose()*coordBuffer;
  }
  else {
    result = basis_.block(0,startCol_,basis_.rows(),blockCols_).transpose()*FullCoordinates;
  }
#endif
  return _result;
}

template <>
GenVector<double> &
GenVecBasis<double, GenVector>::reduceAll(GenVector<double> &x, GenVector<double> &_result) const {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());

  result = basis_.transpose()*FullCoordinates;
#endif
  return _result;
}

template<>
void
GenVecBasis<double, GenDistrVector>::makeSparseBasis(const std::vector<std::vector<int> > & nodeVec, DofSetArray **dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  compressedKey_.clear();
  for(int n = 0; n < nodeVec.size(); n++) {
    for(int i = 0; i < nodeVec[n].size(); i++) {
      dof1 = dsa[n]->firstdof(nodeVec[n][i]);
      numdofs = dsa[n]->weight(nodeVec[n][i]);
      for(int j = 0; j < numdofs; j++) {
        compressedKey_.push_back(vectors_[0].subOffset(n)+dof1+j);
      }
    }
  }

  compressedBasis_.resize(compressedKey_.size(), blockCols_);

  for(int i = 0; i < compressedKey_.size(); i++) {
    compressedBasis_.row(i) = basis_.row(compressedKey_[i]).segment(startCol_, blockCols_);
  }
#endif
}

template<>
void
GenVecBasis<double, GenDistrVector>::makeSparseBasis(const std::vector<std::vector<std::pair<int, int> > > & nodeVec, DofSetArray **dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  compressedKey_.clear();
  for(int n = 0; n < nodeVec.size(); n++) {
    const std::vector<std::pair<int,int> > &subNVMap = nodeVec[n];
    for(std::vector<std::pair<int,int> >::const_iterator it = subNVMap.begin(); it != subNVMap.end(); it++) {
      dof1 = dsa[n]->firstdof(it->first);
      compressedKey_.push_back(vectors_[0].subOffset(n)+dof1+(it->second));
    }
  }

  new (&compressedBasis_) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey_.size(), vectorCount_);

  for(int i = 0; i < compressedKey_.size(); i++) {
    compressedBasis_.row(i) = basis_.row(compressedKey_[i]);
  }
#endif
}

template<>
void
GenVecBasis<double, GenVector>::makeSparseBasis(const std::vector<std::pair<int, int> > & nodeVec, DofSetArray *dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  compressedKey_.clear();
   for(std::vector<std::pair<int,int> >::const_iterator it = nodeVec.begin(); it != nodeVec.end(); it++) {
     dof1 = dsa->firstdof(it->first);
     compressedKey_.push_back(dof1+(it->second));
   }

  new (&compressedBasis_) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey_.size(), vectorCount_);

  for(int i = 0; i < compressedKey_.size(); i++) {
    compressedBasis_.row(i) = basis_.row(compressedKey_[i]);
  }
#endif
}

template<>
void
GenVecBasis<double, GenDistrVector>::makeSparseBasis(const std::vector<std::map<int,std::vector<int> > > & nodeVec)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;
  std::vector<int> dummyKey;
  std::set<int> rowSet;

  compressedKey_.clear();
  for(int n = 0; n < nodeVec.size(); n++) {//loop over subdomains
    const std::map<int,std::vector<int> > &subNVMap = nodeVec[n];
    for(std::map<int,std::vector<int> >::const_iterator it = subNVMap.begin(); it != subNVMap.end(); it++) { //loop over elements
      const std::vector<int> &ColMap = it->second;
      for(std::vector<int>::const_iterator it2 = ColMap.begin(); it2 != ColMap.end(); it2++){//loop over DOFS
      dummyKey.push_back(vectors_[0].subOffset(n)+(*it2));
      rowSet.insert(*it2);
      }
    }
  }

  compressedKey_.resize(dummyKey.size());
  for(int row = 0; row != dummyKey.size(); row++)
     compressedKey_[dummyKey[row]] = row;

  new (&compressedBasis_) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey_.size(), vectorCount_);

  int row = 0;
  for(std::set<int>::const_iterator it = rowSet.begin(); it != rowSet.end(); it++) {
    compressedBasis_.row(row) = basis_.row(*it);
    row++;
  }
#endif
}

template<>
void
GenVecBasis<double, GenVector>::makeSparseBasis(const std::map<int,std::vector<int> > &nodeVec)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;
  std::vector<int> dummyKey;
  std::set<int> rowSet;

  compressedKey_.clear();
  for(std::map<int,std::vector<int> >::const_iterator it = nodeVec.begin(); it != nodeVec.end(); it++) { //loop over elements
    const std::vector<int> &ColMap = it->second;
    for(std::vector<int>::const_iterator it2 = ColMap.begin(); it2 != ColMap.end(); it2++){//loop over DOFS
    dummyKey.push_back(*it2);
    rowSet.insert(*it2);
    }
  }

  compressedKey_.resize(dummyKey.size());
  for(int row = 0; row != dummyKey.size(); row++)
     compressedKey_[dummyKey[row]] = row;

  new (&compressedBasis_) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey_.size(), vectorCount_);

  int row = 0;
  for(std::set<int>::const_iterator it = rowSet.begin(); it != rowSet.end(); it++) {
    compressedBasis_.row(row) = basis_.row(*it);
    row++;
  }
#endif
}

template<>
void
GenVecBasis<double, GenDistrVector>::makeSparseBasis2(const std::vector<std::vector<int> > & nodeVec, DofSetArray **dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  compressedKey2_.clear();
  for(int n = 0; n < nodeVec.size(); n++) {
    for(int i = 0; i < nodeVec[n].size(); i++) {
      dof1 = dsa[n]->firstdof(nodeVec[n][i]);
      numdofs = dsa[n]->weight(nodeVec[n][i]);
      for(int j = 0; j < numdofs; j++) {
        compressedKey2_.push_back(vectors_[0].subOffset(n)+dof1+j);
      }
    }
  }

  new (&compressedBasis2_) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey2_.size(), vectorCount_);

  for(int i = 0; i < compressedKey2_.size(); i++) {
    compressedBasis2_.row(i) = basis_.row(compressedKey2_[i]);
  }
#endif
}

template<>
void
GenVecBasis<double, GenVector>::makeSparseBasis(const std::vector<int> & nodeVec, DofSetArray *dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  compressedKey_.clear();
  for(int i = 0; i < nodeVec.size(); i++) {
    dof1 = dsa->firstdof(nodeVec[i]);
    numdofs = dsa->weight(nodeVec[i]);
    for(int j = 0; j < numdofs; j++) {
      compressedKey_.push_back(dof1+j);
    }
  }

  compressedBasis_.resize(compressedKey_.size(), blockCols_);

  for(int i = 0; i < compressedKey_.size(); i++) {
    compressedBasis_.row(i) = basis_.row(compressedKey_[i]).segment(startCol_, blockCols_);
  }
#endif
}

template<>
void
GenVecBasis<double, GenVector>::makeSparseBasis(const std::vector<int> & nodeVec, DofSetArray *dsa, DofSetArray *reduced_dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs, numdofs2;

  compressedKey_.clear();
  for(int i = 0; i < nodeVec.size(); i++) {
    dof1 = dsa->firstdof(nodeVec[i]);
    numdofs = dsa->weight(nodeVec[i]);
    numdofs2 = reduced_dsa->weight(i);
    int *rdofs = new int[numdofs];
    reduced_dsa->number(i, (*dsa)[nodeVec[i]], rdofs);
    for(int j = 0; j < numdofs; j++) {
      // note: in this version we only push back if dof is in reduced dsa
      if(rdofs[j] >= 0)
        compressedKey_.push_back(dof1+j);
    }
    delete [] rdofs;
  }

  new (&compressedBasis_) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey_.size(), vectorCount_);

  for(int i = 0; i < compressedKey_.size(); i++) {
    compressedBasis_.row(i) = basis_.row(compressedKey_[i]);
  }
#endif
}

template<>
void
GenVecBasis<double, GenVector>::makeSparseBasis2(const std::vector<int> & nodeVec, DofSetArray *dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  compressedKey2_.clear();
  for(int i = 0; i < nodeVec.size(); i++) {
    dof1 = dsa->firstdof(nodeVec[i]);
    numdofs = dsa->weight(nodeVec[i]);
    for(int j = 0; j < numdofs; j++) {
      compressedKey2_.push_back(dof1+j);
    }
  }

  new (&compressedBasis2_) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey2_.size(), vectorCount_);

  for(int i = 0; i < compressedKey2_.size(); i++) {
    compressedBasis2_.row(i) = basis_.row(compressedKey2_[i]);
  }
#endif
}

template <>
GenVector<double> &
GenVecBasis<double, GenVector>::expand2(GenVector<double> &x, GenVector<double> &_result) const {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());

  //full coordinates distributed over MPI processes
  if(compressedKey2_.size() > 0) {
    Eigen::VectorXd resultBuffer(compressedKey2_.size());
    resultBuffer = compressedBasis2_*GenCoordinates;
    for(int i = 0; i < compressedKey2_.size(); i++)
      result(compressedKey2_[i]) = resultBuffer(i);
  }
  else {
    result = basis_*GenCoordinates;
  }
#endif
  return _result;
}

template<>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::addLocalPart(GenDistrVector<double> & _x, GenDistrVector<double> & _y) const{
#ifdef USE_EIGEN3
  if(_x.size() > 0 ) {
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > x(_x.data(), _x.size());
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > y(_y.data(), _y.size());

    y.segment(startCol_,blockCols_) += x.segment(startCol_,blockCols_);
  }
#endif
  return _y;
}

template <>
GenVector<double> &
GenVecBasis<double, GenVector>::addLocalPart(GenVector<double> &_x, GenVector<double> &_y) const {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > x(_x.data(), _x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > y(_y.data(), _y.size());

  y.segment(startCol_,blockCols_) += x.segment(startCol_,blockCols_);
#endif
  return _y;
}

}
