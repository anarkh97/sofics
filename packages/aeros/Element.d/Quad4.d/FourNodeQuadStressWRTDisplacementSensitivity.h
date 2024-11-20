#ifdef USE_EIGEN3
#ifndef _FOURNODEQUADSTRESSWRTDISPLACEMENTSENSITIVITY_H_
#define _FOURNODEQUADSTRESSWRTDISPLACEMENTSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Quad4.d/QuadElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class FourNodeQuadStressWRTDisplacementSensitivity : public VectorValuedFunction<8,4,Scalar,17,1,double>
{
  public:
    QuadElementTemplate<Scalar> ele;
    Eigen::Array<Scalar,4,1> globalx, globaly; // nodal coordinates
    Scalar A, E, nu, alpha, Ta; // material properties
    Eigen::Array<Scalar,4,1> ndTemps;

  public:
    FourNodeQuadStressWRTDisplacementSensitivity(const Eigen::Array<double,17,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<4>(0).cast<Scalar>();
      globaly = sconst.segment<4>(4).cast<Scalar>();
      E = sconst[8];
      A = sconst[9];
      nu = sconst[10]; 
      alpha = sconst[11];
      Ta = sconst[12];
      ndTemps[0] = sconst[13];
      ndTemps[1] = sconst[14];
      ndTemps[2] = sconst[15];
      ndTemps[3] = sconst[16];
    }

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,8,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      Eigen::Matrix<Scalar,8,1> globalu = q;

      Eigen::Array<Scalar,7,4> elStress;
      Eigen::Array<Scalar,7,4> elStrain;
      int numel = 1;
      int elm = 1;
      int maxgus = 4;
      int maxstr = 7;
      bool vmflg = true;
      bool strainFlg = false;
      char escm[7] = "extrap";  // ... stress extrapolation from gauss points
      Scalar c[9];
      ele.getcmt(A, E, nu, c);

      Scalar tc = E*alpha/(1.0-nu);

      ele.sands2(escm, globalx.data(), globaly.data(), c,
                 globalu.data(),elStress.data(),
                 elStrain.data(), maxgus, maxstr, elm, numel,
                 vmflg, strainFlg, tc, Ta, ndTemps.data());

      // return value:

      Eigen::Matrix<Scalar,4,1> v;
      v[0] = elStress(6,0); 
      v[1] = elStress(6,1);
      v[2] = elStress(6,2);
      v[3] = elStress(6,3);

      return v; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
