#ifndef _FOURFIELDSTRAINENERGYFUNCTION_H_
#define _FOURFIELDSTRAINENERGYFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/StrainEnergyDensity.d/IsoLinearElasticStrainEnergyDensityFunction.h>
#include <Element.d/Function.d/StrainEnergyDensity.d/SaintVenantKirchhoffStrainEnergyDensityFunction.h>
#include <Element.d/Function.d/StrainEnergyDensity.d/NeoHookeanStrainEnergyDensityFunction.h>
#include <Element.d/Function.d/StrainEnergyDensity.d/MooneyRivlinStrainEnergyDensityFunction.h>
#include <Element.d/Function.d/StrainEnergyDensity.d/OgdenStrainEnergyDensityFunction.h>

#include <cassert>
#include <iostream>

// References
// [1] J.C. Simo, R.L. Taylor and K.S. Pister, "Variational and projection methods for the volume constraint in finite
//     deformation plasticity", Computer Methods in Applied Mechanics and Engineering, 51 (1985) 177-208.
// [2] J.C. Simo and R.L. Taylor, "Quasi-incompressible finite elasticity in principal stretches. Continuum basis and
//     numerical algorithms", Computer Methods in Applied Mechanics and Engineering, 8.5 (1991) 273-310.

namespace Simo {

template<typename Scalar, 
         template <typename S> class ShapeFunctionTemplate,
         template <typename S> class ShapeFunctionTemplate2,
         template <typename S> class ShapeFunctionTemplate3,
         template <typename S> class ShapeFunctionTemplate4,
         typename QuadratureRule>
class FourFieldStrainEnergyFunction
 : public ScalarValuedFunction<ShapeFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates*ShapeFunctionTemplate<Scalar>::NumberOfValues
                               +ShapeFunctionTemplate2<Scalar>::NumberOfValues
                               +ShapeFunctionTemplate3<Scalar>::NumberOfValues,
                               Scalar,
                               ShapeFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates*ShapeFunctionTemplate<Scalar>::NumberOfValues+4
                               +ShapeFunctionTemplate2<Scalar>::NumberOfValues
                               +ShapeFunctionTemplate3<Scalar>::NumberOfValues
                               +ShapeFunctionTemplate4<Scalar>::NumberOfValues,
                               2,
                               double>
{
  public:
    enum {
      NumberOfNodes  = ShapeFunctionTemplate<Scalar>::NumberOfValues,
      MaximumNumberOfMaterialConstants = 6,
      NumberOfNodes2 = ShapeFunctionTemplate2<Scalar>::NumberOfValues,
      NumberOfNodes3 = ShapeFunctionTemplate3<Scalar>::NumberOfValues,
      NumberOfNodes4 = ShapeFunctionTemplate4<Scalar>::NumberOfValues,
      NumberOfDimensions = ShapeFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      NumberOfGeneralizedCoordinates = NumberOfNodes*NumberOfDimensions+NumberOfNodes2+NumberOfNodes3,
      NumberOfScalarConstants = NumberOfGeneralizedCoordinates+MaximumNumberOfMaterialConstants+NumberOfNodes4+1,
      NumberOfIntegerConstants = 2
    };

  private:
    Eigen::Matrix<double,NumberOfNodes,NumberOfDimensions> X;
    Eigen::Array<double,MaximumNumberOfMaterialConstants,1> c;
    Eigen::Matrix<double,NumberOfNodes2,1> Theta0;  // dilatation reference configuration
    Eigen::Matrix<double,NumberOfNodes3,1> p0;      // pressure reference configuration
    Eigen::Matrix<double,NumberOfNodes4,1> lambda;  // Lagrange multiplier
    double epsilon; // penalty parameter
    int deg;  // quadrature rule degree
    int wopt; // strain energy density function

  public:
    FourFieldStrainEnergyFunction(const Eigen::Array<double,NumberOfScalarConstants,1>& sconst,
                                  const Eigen::Array<int,NumberOfIntegerConstants,1>& iconst)
    {
      for(int inode = 0; inode < NumberOfNodes; ++inode) {
        X.row(inode) = Eigen::Map<Eigen::Matrix<double,NumberOfDimensions,1>>(const_cast<double*>(sconst.data())+NumberOfDimensions*inode);
      }
      c       = sconst.template segment<MaximumNumberOfMaterialConstants>(NumberOfDimensions*NumberOfNodes);
      Theta0  = sconst.template segment<NumberOfNodes2>(NumberOfDimensions*NumberOfNodes+MaximumNumberOfMaterialConstants);
      p0      = sconst.template segment<NumberOfNodes3>(NumberOfDimensions*NumberOfNodes+MaximumNumberOfMaterialConstants+NumberOfNodes2);
      lambda  = sconst.template segment<NumberOfNodes4>(NumberOfDimensions*NumberOfNodes+MaximumNumberOfMaterialConstants+NumberOfNodes2+NumberOfNodes3);
      epsilon = sconst[NumberOfScalarConstants-1];

      deg  = iconst[0];
      wopt = iconst[1];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t)
    {
      // q[0] = x translation of node 1
      // q[1] = y translation of node 1
      // q[2] = z translation of node 1
      // q[3] = x translation of node 2
      // q[4] = y translation of node 2
      // q[5] = z translation of node 2
      // etc...

      // return value: strain energy
      Scalar V = 0.0;

      using std::abs;
      using std::sqrt;
      using std::pow;

      // Set current configuration
      Eigen::Matrix<Scalar,NumberOfNodes,NumberOfDimensions> x;
      Eigen::Matrix<Scalar,NumberOfNodes2,1> Theta;
      Eigen::Matrix<Scalar,NumberOfNodes3,1> p;
      for(int inode = 0; inode < NumberOfNodes; ++inode) {
        x.row(inode) = X.row(inode).template cast<Scalar>() + q.template segment<NumberOfDimensions>(inode*NumberOfDimensions).transpose();
      }
      Theta = Theta0.template cast<Scalar>() + q.template segment<NumberOfNodes2>(NumberOfDimensions*NumberOfNodes);
      p     = p0.template cast<Scalar>()     + q.template segment<NumberOfNodes3>(NumberOfDimensions*NumberOfNodes+NumberOfNodes2);

      // strain energy density
      StrainEnergyDensityFunction<Scalar> *W;
      switch(wopt) {
        case 0 :
          W = new IsoLinearElasticStrainEnergyDensityFunction<Scalar>(c[0],c[1]);
          break;
        case 1 :
          W = new SaintVenantKirchhoffStrainEnergyDensityFunction<Scalar>(c[0],c[1]);
          break;
        case 3 :
          W = new NeoHookeanStrainEnergyDensityFunction<Scalar>(c[1]);
          break;
        case 4 :
          W = new MooneyRivlinStrainEnergyDensityFunction<Scalar>(c[0],c[1]);
          break;
        case 5 : 
          W = new OgdenStrainEnergyDensityFunction<Scalar>(c[0],c[1],c[2],c[3],c[4],c[5]);
          break;
        default : {
          std::cerr << "Error: unsupported strain energy density in FourFieldStrainEnergyFunction\n";
          exit(-1);
        }
      }

      // shape function and derivatives
      ShapeFunctionTemplate2<double> N2;
      ShapeFunctionTemplate3<double> N3;
      ShapeFunctionTemplate4<double> N4;
      Jacobian<double, ShapeFunctionTemplate> dN(Eigen::Array<double,0,1>::Zero(), Eigen::Array<int,0,1>::Zero());

      // quadrature rule
      QuadratureRule qrule(deg);

      // local variables to be computed at each integration point
      Eigen::Matrix<double,NumberOfDimensions,1> xi;
      double weight;
      Eigen::Matrix<double,NumberOfNodes,NumberOfDimensions> dNdxi;
      Eigen::Matrix<double,NumberOfDimensions,NumberOfDimensions> J;
      Eigen::Matrix<double,NumberOfNodes,NumberOfDimensions> nGrad;
      Eigen::Matrix<Scalar,NumberOfDimensions,NumberOfDimensions> F, Fbar;
      Scalar Theta_G, p_G;
      double lambda_G;

      // Loop over the integration points
      for(int ip = 0; ip < qrule.getN(); ++ip) {

        // get the integration point abscissa and weight
        qrule.getAbscissaAndWeight(ip, xi, weight);

        // compute shape function derivatives
        dNdxi = dN(xi, 0.);

        // compute the jacobian of the isoparametric mapping and the deformation gradient
        J = dNdxi.transpose()*X;
        nGrad = dNdxi*J.transpose().inverse();
        F = x.transpose().lazyProduct(nGrad);

        // interpolate the scalar fields
        Theta_G  = N2(xi,0.).template cast<Scalar>().dot(Theta); assert(Theta_G > 0);
        p_G      = N3(xi,0.).template cast<Scalar>().dot(p);
        lambda_G = N4(xi,0.).dot(lambda);

        // compute the modified deformation gradient
        Scalar j = F.determinant(); assert(j > 0);
        Fbar = pow(Theta_G/j,1/3.)*F;

        // evaluate the augmented strain energy density function at the integration point
        V += (abs(J.determinant())*weight)*((*W)(Fbar) + p_G*(j-Theta_G) + epsilon*0.5*((Theta_G-1)*(Theta_G-1)) + lambda_G*(Theta_G-1));
      }

      delete W;

      return V;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo

#endif
