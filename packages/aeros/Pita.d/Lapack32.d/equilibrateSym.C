#include <algorithm>
#include <cmath>

namespace Pita {  

void equilibrateSym(int matrixSize, const double * data, double relTol, double * scaling, double * scond, double * amax) {
  int diagShift = matrixSize + 1;
  const double * diagEnd = data + (matrixSize * diagShift);
 
  double diagMax = data[0];
  double diagMin = data[0];
  for (const double * diagAddr = data + diagShift; diagAddr < diagEnd; diagAddr += diagShift) {
    diagMax = std::max(diagMax, *diagAddr);
    diagMin = std::min(diagMin, *diagAddr);
  }
  const double absTol = diagMax * relTol;

  double * scalingAddr = scaling;
  for (const double * diagAddr = data; diagAddr < diagEnd; diagAddr += diagShift) {
    if (*diagAddr >= absTol) {
      *scalingAddr = 1.0 / std::sqrt(*diagAddr);
    } else {
      *scalingAddr = 0.0;
    }
    ++scalingAddr;
  }

  *amax = diagMax;
  *scond = std::sqrt(diagMin) / std::sqrt(diagMax);
}

} // end namespace Pita
