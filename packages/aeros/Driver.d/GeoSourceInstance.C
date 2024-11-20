#include <Driver.d/GeoSource.C>

template
void GeoSource::distributeBCs<double>(GenSubDomain<double>*&, int*, int*);
template
void GeoSource::distributeBCs<complex<double> >(GenSubDomain<complex<double> >*&, int*, int*);

#ifdef DISTRIBUTED
template
void GeoSource::distributeOutputNodes<double>(GenSubDomain<double>*&, int*, int*);
template
void GeoSource::distributeOutputNodes<complex<double> >(GenSubDomain<complex<double> >*&, int*, int*);

template
void GeoSource::distributeOutputNodesX<double>(GenSubDomain<double>*, Connectivity*);
template
void GeoSource::distributeOutputNodesX<complex<double> >(GenSubDomain<complex<double> >*, Connectivity*);
#endif

template
std::vector<GenSubDomain<double> *> GeoSource::readDistributedInputFiles<double>(gsl::span<const int>);
template
std::vector<GenSubDomain<complex<double> > *> GeoSource::readDistributedInputFiles<complex<double> >(gsl::span<const int>);

#ifdef USE_EIGEN3
template
void GeoSource::outputSensitivityScalars<double>(int, Eigen::Matrix<double, Eigen::Dynamic, 1>*, double,
                                                 Eigen::Matrix<double, Eigen::Dynamic, 1> *);
template
void GeoSource::outputSensitivityScalars<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic, 1>*,
                                                           double, Eigen::Matrix<double, Eigen::Dynamic, 1> *);

template
void GeoSource::outputSensitivityVectors<double>(int, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*, double,
                                                 Eigen::Matrix<double, Eigen::Dynamic, 1> *);
template
void GeoSource::outputSensitivityVectors<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>*,
                                                           double, Eigen::Matrix<double, Eigen::Dynamic, 1> *);

template
void GeoSource::outputSensitivityDispVectors<double>(int, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>**,
                                                     double, int, int);

template
void GeoSource::outputSensitivityDispVectors<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic,
                                                               Eigen::Dynamic>**, double, int, int);

template
void GeoSource::outputSensitivityAdjointStressVectors<double>(int, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*,
                                                              double *, double, int, std::vector<int>,
                                                              Eigen::Matrix<double, Eigen::Dynamic, 1> *);
template
void GeoSource::outputSensitivityAdjointStressVectors<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic,
                                                                        Eigen::Dynamic>*, complex<double> *, double, int,
                                                                        std::vector<int>,Eigen::Matrix<double, Eigen::Dynamic, 1> *);

template
void GeoSource::outputSensitivityAdjointDispVectors<double>(int, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>**,
                                                            double *, double, int, std::vector<DispNode>,
                                                            Eigen::Matrix<double, Eigen::Dynamic, 1> *);
template
void GeoSource::outputSensitivityAdjointDispVectors<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic,
                                                                      Eigen::Dynamic>**, complex<double> *, double, int,
                                                                      std::vector<DispNode>,Eigen::Matrix<double,
                                                                      Eigen::Dynamic, 1> *);

template
void GeoSource::outputSensitivityDispVectors<double>(int, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*,
                                                     double, int);

template
void GeoSource::outputSensitivityDispVectors<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic,
                                                               Eigen::Dynamic>*,double, int);
#endif

#define GEOSOURCE_INSTANTIATION_HELPER(dim) \
template void GeoSource::outputNodeVectors<dim>(int, double (*)[dim], int, double);\
template void GeoSource::outputNodeVectors<dim>(int, DComplex (*)[dim], int, double);
#define GEOSOURCE_INSTANTIATION_HELPER6(dim) \
template void GeoSource::outputNodeVectors6<dim>(int, double (*)[dim], int, double);\
template void GeoSource::outputNodeVectors6<dim>(int, DComplex (*)[dim], int, double);
#define GEOSOURCE_INSTANTIATION_HELPER9(dim) \
template void GeoSource::outputNodeVectors9<dim>(int, double (*)[dim], int, double);
#define GEOSOURCE_INSTANTIATION_HELPER4(dim) \
template void GeoSource::outputNodeVectors4<dim>(int, double (*)[dim], int, double);

GEOSOURCE_INSTANTIATION_HELPER(3);
GEOSOURCE_INSTANTIATION_HELPER(4);
GEOSOURCE_INSTANTIATION_HELPER(6);
GEOSOURCE_INSTANTIATION_HELPER(9);
GEOSOURCE_INSTANTIATION_HELPER(11);

GEOSOURCE_INSTANTIATION_HELPER4(4);
GEOSOURCE_INSTANTIATION_HELPER4(6);
GEOSOURCE_INSTANTIATION_HELPER4(9);
GEOSOURCE_INSTANTIATION_HELPER4(11);

GEOSOURCE_INSTANTIATION_HELPER6(6);
GEOSOURCE_INSTANTIATION_HELPER6(9);
GEOSOURCE_INSTANTIATION_HELPER6(11);

GEOSOURCE_INSTANTIATION_HELPER9(9);
GEOSOURCE_INSTANTIATION_HELPER9(11);
