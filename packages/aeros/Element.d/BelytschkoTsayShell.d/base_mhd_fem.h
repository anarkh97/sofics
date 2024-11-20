#ifdef USE_EIGEN3
#include <Eigen/Core>

void rotprojbt1(const Eigen::Matrix<double,3,4> &ecurnloc, double *_evecloc0,
                double *_evecloc);
void rotprojbt2(const Eigen::Matrix<double,3,4> &ecurnloc, double *_evecloc0,
                double *_evecloc);
void getrotpmatbt(int inode, const Eigen::Matrix<double,3,4> &ecurnloc,
                  Eigen::Matrix<double,3,3> &rotpmat);
void getrotpmatnbt(int inode, const Eigen::Matrix<double,3,4> &ecurnloc,
                   Eigen::Matrix<double,3,3> &rotpmat, Eigen::Matrix<double,3,1> &e3vec);
#endif
