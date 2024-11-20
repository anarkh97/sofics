#ifdef USE_EIGEN3
#ifndef _TETRAELEMENTSTRESSWRTNODALCOORDINATESENSITIVITY_H_
#define _TETRAELEMENTSTRESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Tetra.d/TetraElementTemplate.cpp>

template<typename Scalar>
class TetraElementStressWRTNodalCoordinateSensitivity : public VectorValuedFunction<12,4,Scalar,14,1,double>
{
  public:
    TetraElementTemplate<Scalar> ele;
    Eigen::Array<Scalar,12,1> globalu; // element displacements
    Scalar E, nu; // material properties
    int surface;

  public:
    TetraElementStressWRTNodalCoordinateSensitivity(const Eigen::Array<double,14,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalu = sconst.segment<12>(0).cast<Scalar>();
      E = Scalar(sconst[12]);
      nu = Scalar(sconst[13]);
      surface = iconst[0];  
    }

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,12,1>& q, Scalar)
    {
      Eigen::Array<Scalar,4,1> globalx; 
      Eigen::Array<Scalar,4,1> globaly; 
      Eigen::Array<Scalar,4,1> globalz; 
      globalx << q[0], q[3], q[6], q[9];
      globaly << q[1], q[4], q[7], q[10];
      globalz << q[2], q[5], q[8], q[11];

      int maxgus = 4;
      int maxstr = 7;
      int outerr = 6;
      int elm = 1;
//      bool vmflg = true;
//      bool strainflg = false;

      Eigen::Array<Scalar,7,4> stress;
      Eigen::Array<Scalar,7,4> strain;
      ele.sands23(elm, globalx.data(), globaly.data(), globalz.data(), E, nu, globalu.data(),
                  stress.data(), strain.data(), 4, 7, 1, outerr, 1, 0);

      // return value:
      // von mises stresses at nodes  
      Eigen::Matrix<Scalar,4,1> v;
      v << stress(6,0), stress(6,1), stress(6,2), stress(6,3);

      return v; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
