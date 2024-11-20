#include <Driver.d/BinaryOutputInclude.C>

#define BINARYOUTPUT_INSTANTIATION_HELPER(Scalar,dim) \
template void GeoSource::writeNodeVectorToFile<Scalar,dim>(SVec<Scalar, dim>&, int, int, int, int, int, double, int, int, int*);

BINARYOUTPUT_INSTANTIATION_HELPER(double,11);
BINARYOUTPUT_INSTANTIATION_HELPER(double,6);
BINARYOUTPUT_INSTANTIATION_HELPER(complex<double>,11);
BINARYOUTPUT_INSTANTIATION_HELPER(complex<double>,6);
