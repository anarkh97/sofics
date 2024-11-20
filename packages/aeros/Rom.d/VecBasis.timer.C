#include "VecBasis.h"

namespace Rom{

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::project(GenDistrVector<double> &x, GenDistrVector<double> &_result) {
  std::cout << "\n counter = " << counter << std::endl;
  counter += 1;

  double dummy = 0;

  dummy -= getTime();
  Eigen::Matrix<double, Eigen::Dynamic, 1> GenCoordinates;
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > f(x.data(), x.size());
  Eigen::SparseVector<double> sparsef(f.size());
  dummy += getTime();
//  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > basis(vectors_[0].data(), vectors_[0].size(), numVectors());

  time1 += dummy;

  std::cout << " \n             matrix allocation time   = " << time1/counter << std::endl;

  long dummycount = 0;
  long dummycount2 = 0;
  dummy = 0;
  dummy -= getTime();
  for (int i = 0; i < f.size(); i++) {
   if(f(i) != 0.) {
    sparsef.insert(i) = f(i);
    dummycount += 1;
    }
   dummycount2 += 1;
  }
  dummy += getTime();

  time2 += dummy;

  std::cout << "             sparse vector allocation = " << time2/counter << std::endl;

//  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned, Eigen::OuterStride<> > basisT(vectors_[0].data(), vectors_[0].size(), vectorCount(), Eigen::OuterStride<>(size()));
  
  dummy = 0;
  dummy -= getTime();
  GenCoordinates = basis.transpose()*sparsef;
  dummy += getTime();

  time3 += dummy;

  std::cout << "       basis^T * sparseforce vector t = " << time3/counter << std::endl;

  dummy = 0;
  dummy -= getTime();
  GenCoordinates = basis.transpose()*f;
  dummy += getTime();
 
  time4 += dummy;

  std::cout << "             basis^T * force vector t = " << time4/counter << std::endl;

  dummy = 0;
  dummy -= getTime();
  if(structCom)
    structCom->globalSum(GenCoordinates.size(), GenCoordinates.data());
  dummy += getTime();

  time5 += dummy;

  std::cout << "             global sum reduction time= " << time5/counter << std::endl;

  dummy = 0;
  dummy -= getTime();
  result = basis*GenCoordinates;
  dummy += getTime();

  time6 += dummy;

  std::cout << "             Basis * gencoordinates   = " << time6/counter << std::endl;
  
  std::cout << "                   percent nodes left = " << 100*double(dummycount)/double(dummycount2) << std::endl;

   return _result;
}

}
