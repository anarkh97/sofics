#include "ElementSamplingDriver.C"

namespace Rom {

template
void
outputFullWeights<Vector,std::vector<int> >(const Vector &weights, const std::vector<int> &elemIds, int j);

template
void
outputFullWeights<std::vector<double>,std::vector<int> >(const std::vector<double> &weights, const std::vector<int> &elemIds, int j);

template
int
ElementSamplingDriver<std::vector<double>,size_t>
::elementCount() const;

template
int
ElementSamplingDriver<std::vector<double>,size_t>
::vectorSize() const;

template
ElementSamplingDriver<std::vector<double>,size_t>
::ElementSamplingDriver(Domain *d);

template
ElementSamplingDriver<std::vector<double>,size_t>
::~ElementSamplingDriver();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::assembleTrainingData(VecBasis &podBasis, int podVectorCount, VecBasis &displac,
                       VecBasis *veloc, VecBasis *accel, VecBasis *ndscfgCoords,
                       VecBasis *ndscfgMassOrthogonalBases, int j);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::assembleTrainingData(std::vector<StackVector> &podBasis, int podVectorCount, std::vector<StackVector> &displac,
                       std::vector<StackVector> *veloc, std::vector<StackVector> *accel, std::vector<StackVector> *ndscfgCoords,
                       std::vector<StackVector> *ndscfgMassOrthogonalBases, int j);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::solve();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::computeSolution(Vector &solution, double tol, bool verboseFlag);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::postProcess(Vector &solution, bool verboseFlag);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::preProcess();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::buildDomainCdsa();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::clean();

} // end namespace Rom
