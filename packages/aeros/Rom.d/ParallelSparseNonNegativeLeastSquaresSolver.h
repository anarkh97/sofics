#ifndef ROM_PARALLELSPARSENONNEGATIVELEASTSQUARESSOLVER_H
#define ROM_PARALLELSPARSENONNEGATIVELEASTSQUARESSOLVER_H

#include "SparseNonNegativeLeastSquaresSolver.h"
#include "SimpleBuffer.h"
#include <vector>
#include <utility>
#include <list>

struct long_int {
  long index;
  int sub;
};

namespace Rom {

class ParallelSparseNonNegativeLeastSquaresSolver {
public:

  int subdomainCount() const { return nsub_; }
  SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> * subdomainSolver(int i) { return sd_[i]; }

  // Problem size
  long equationCount() const { return equationCount_; }
  long unknownCount() const { return unknownCount_; }

  void problemSizeIs(long eqnCount, long unkCount);

  // Options
  double relativeTolerance() const { return relativeTolerance_; }
  void relativeToleranceIs(double relTol) { relativeTolerance_ = relTol; }

  bool verboseFlag() const { return verboseFlag_; }
  void verboseFlagIs(bool verFlg) { verboseFlag_ = verFlg; }

  bool scalingFlag() const { return scalingFlag_; }
  void scalingFlagIs(bool scaFlg) { scalingFlag_ = scaFlg; }

  bool centerFlag() const { return centerFlag_; }
  void centerFlagIs(bool scaFlg) { centerFlag_ = scaFlg; }

  bool projectFlag() const { return projectFlag_; }
  void projectFlagIs(bool scaFlg) { projectFlag_ = scaFlg; }

  bool positivityFlag() const { return positivity_; }
  void positivityIs(bool posFlg) { positivity_ = posFlg; }

  int solverType() const { return solverType_; }
  void solverTypeIs(int solTyp) { solverType_ = solTyp; }

  double maxSizeRatio() const { return maxSizeRatio_; }
  void maxSizeRatioIs(double maxSze) { maxSizeRatio_ = maxSze; }

  double maxIterRatio() const { return maxIterRatio_; }
  void maxIterRatioIs(double maxIte) { maxIterRatio_ = maxIte; }

  void useHotStart(bool hs) { hotStart_ = hs; }

  // Rhs buffer: [equationCount]
  const double * rhsBuffer() const { return rhsBuffer_.array(); }
  double * rhsBuffer() { return rhsBuffer_.array(); }

  // Error magnitude
  double errorMagnitude() const { return errorMagnitude_; }

  // Scalapack LH solver arguments
  int npMax() const { return npMax_; }
  void npMaxIs(int npMax) { npMax_ = npMax; }

  int scpkMB() const { return scpkMB_; }
  void scpkMBIs(int scpkMB) { scpkMB_ = scpkMB; }

  int scpkNB() const { return scpkNB_; }
  void scpkNBIs(int scpkNB) { scpkNB_ = scpkNB; }

  int scpkMP() const { return scpkMP_; }
  void scpkMPIs(int scpkMP) { scpkMP_ = scpkMP; }

  int scpkNP() const { return scpkNP_; }
  void scpkNPIs(int scpkNP) { scpkNP_ = scpkNP; }


  // Perform solve
  void solve();

  // Constructor
  ParallelSparseNonNegativeLeastSquaresSolver(int nsub, SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> **sd);

private:

  long equationCount_;
  long unknownCount_;

  SimpleBuffer<double> rhsBuffer_;

  std::list<std::pair<int,long_int> > hotIndices;

  double relativeTolerance_;
  double errorMagnitude_;
  bool verboseFlag_;
  bool scalingFlag_;
  bool centerFlag_;
  bool projectFlag_;
  bool positivity_;
  bool hotStart_;
  int solverType_;
  double maxSizeRatio_;
  double maxIterRatio_;

  int nsub_;
  SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> **sd_;

  // Variables for Scalapack LH solver
  int npMax_;  // Maximum sparse mesh size.
  int scpkMB_; // Scalapack row block size.
  int scpkNB_; // Scalapack column block size.
  int scpkMP_; // Scalapack row processor grid size.
  int scpkNP_; // Scalapack column processor grid size.

  // Disallow copy & assignment
  ParallelSparseNonNegativeLeastSquaresSolver(const ParallelSparseNonNegativeLeastSquaresSolver &);
  ParallelSparseNonNegativeLeastSquaresSolver &operator=(const ParallelSparseNonNegativeLeastSquaresSolver &);
};

} // end namespace Rom

#endif /* ROM_PARALLELSPARSENONNEGATIVELEASTSQUARESSOLVER_H */

