#include <Driver.d/EFrameData.C>
#include <complex>

#define NFRAMEDATA_INSTANTIATION_HELPER(Scalar) \
template void NFrameData::transformVector3<Scalar>(Scalar*); \
template void NFrameData::invTransformVector3<Scalar>(Scalar*); \
template void NFrameData::transformVector6<Scalar>(Scalar*); \
template void NFrameData::invTransformVector6<Scalar>(Scalar*); \
template void NFrameData::transformMatrix3<Scalar>(Scalar*); \
template void NFrameData::invTransformMatrix3<Scalar>(Scalar*); \
template void NFrameData::transformMatrix6<Scalar>(Scalar*); \
template void NFrameData::invTransformMatrix6<Scalar>(Scalar*); \
template void NFrameData::transformSymMatrix3<Scalar>(Scalar*); \
template void NFrameData::invTransformSymMatrix3<Scalar>(Scalar*);

NFRAMEDATA_INSTANTIATION_HELPER(double);
NFRAMEDATA_INSTANTIATION_HELPER(std::complex<double>);
