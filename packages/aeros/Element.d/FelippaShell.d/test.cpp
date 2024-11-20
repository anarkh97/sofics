#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/FelippaShell.d/ShellElementStressWRTThicknessSensitivity.h>
#include <Element.d/FelippaShell.d/ShellElementStressWRTDisplacementSensitivity.h>
#include <Element.d/FelippaShell.d/ShellElementGravityForceWRTThicknessSensitivity.h>
#include <Element.d/FelippaShell.d/ShellElementStiffnessWRTThicknessSensitivity.h>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unsupported/Eigen/NumericalDiff>
#include <iostream>
#include <cmath>

int main()
{
  const int num_test = 10000;
  const bool verbose = false;
  const double eps = 1e-10; // used for finite difference approximation
  double max_error_ad = 0, max_error_fd = 0;
  for(int i=0; i<num_test; ++i) {
    // generate random inputs
    Eigen::Array<double,3,1> g = Eigen::Array<double,3,1>::Random()*10,
                             x = Eigen::Array<double,3,1>::Random(),
                             y = Eigen::Array<double,3,1>::Random(),
                             z = Eigen::Array<double,3,1>::Random();
    Eigen::Array<double,18,1> u = Eigen::Array<double,18,1>::Random()*0.001;
    double E = std::abs(Eigen::Array<double,1,1>::Random()[0])*1e11;
    double nu = std::abs(Eigen::Array<double,1,1>::Random()[0])*0.5;
    double rho = std::abs(Eigen::Array<double,1,1>::Random()[0])*1e4;
    double eh = std::abs(Eigen::Array<double,1,1>::Random()[0])*0.01;
    double Ta = 0;
    double W = 0;
    Eigen::Array<double,9,1> eframe = Eigen::Array<double,9,1>::Random();
    eframe.segment<3>(0) = eframe.segment<3>(0).matrix().normalized();
    eframe.segment<3>(3) = (eframe.matrix().segment<3>(0).cross(eframe.matrix().segment<3>(6))).normalized();
    eframe.segment<3>(6) = (eframe.matrix().segment<3>(0).cross(eframe.matrix().segment<3>(3))).normalized();
    Eigen::Array<double,42,1> coefs = Eigen::Array<double,42,1>::Random()*1e6; // XXX
    Eigen::Array<double,3,1> ndtemps = Eigen::Array<double,3,1>::Random()*100;
    int surface = i%3+1;
    int ctyp = 0;
    int gravflg = i%2;
    int flg = i%2;
    int sflg = 1;
    int tflg = 1;

    double nsm = (ctyp == 2 || ctyp == 3) ? rho : 0;

    ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle> elem;
    ShellMaterial<double> *mat;
    if(ctyp == 0) 
      mat = new ShellMaterialType0<double>(E, eh, nu, rho, Ta, W);
    else if(ctyp == 1)
      mat = new ShellMaterialType1<double>(coefs.data(), eframe.data(), rho, eh, Ta);

    // #1 test stress wrt thickness sensitivity
    {
      // real parameters
      Eigen::Array<double,86,1> dconst;
      dconst << x, y, z, u, E, nu, rho, eframe, coefs, Ta, W, ndtemps;

      // integer parameters
      Eigen::Array<int,3,1> iconst;
      iconst << surface, ctyp, sflg;

      // inputs
      Eigen::Matrix<double,1,1> q;
      q[0] = mat->GetShellThickness();

      // function evaluation
      ShellElementStressWRTThicknessSensitivity<double> foo(dconst,iconst);
      Eigen::Matrix<double,3,1> S = foo(q, 0);
      if(verbose) std::cerr << "#1 S = " << S.transpose() << std::endl;

      // Jacobian evaluation using hand-coded routine
      Eigen::Matrix<double,3,1> J;
      elem.andesvmsWRTthic(0, nu, x.data(), y.data(), z.data(), u.data(), J.data(), ctyp, mat, surface, sflg, ndtemps.data());
      if(verbose) std::cerr << "#1 J = " << J.transpose() << std::endl;

      // Jacobian evaluation by automatic differentiation
      Simo::Jacobian<double,ShellElementStressWRTThicknessSensitivity> dSdh(dconst,iconst);
      Eigen::Matrix<double,3,1> J_ad = dSdh(q, 0);
      double error_ad = (J.norm() == 0) ? (J-J_ad).norm() : (J-J_ad).norm()/J.norm();
      if(verbose) std::cerr << "#1 J_ad = " << J_ad.transpose() << std::endl;
      if(verbose || error_ad > 1e-12) std::cerr << "#1 error_ad = " << error_ad << std::endl;
      max_error_ad = std::max(max_error_ad, error_ad);

      // Jacobian approximation by finite difference
      Simo::SpatialView<double,ShellElementStressWRTThicknessSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementStressWRTThicknessSensitivity>,Eigen::Central> cd(sv, eps);
      Eigen::Matrix<double,3,1> J_fd;
      cd.df(q, J_fd);
      double error_fd = (J.norm() == 0) ? (J-J_fd).norm() : (J-J_fd).norm()/J.norm();
      if(verbose) std::cerr << "#1 J_fd = " << J_fd.transpose() << std::endl;
      if(verbose) std::cerr << "#1 error_fd = " << error_fd << std::endl;
      max_error_fd = std::max(max_error_fd, error_fd);
    }

    // #2 test stress wrt displacement sensitivity
    {
      // real parameters
      Eigen::Array<double,69,1> dconst;
      dconst << x, y, z, E, nu, rho, eh, eframe, coefs, Ta, W, ndtemps;

      // integer parameters
      Eigen::Array<int,3,1> iconst;
      iconst << surface, ctyp, sflg;

      // inputs
      Eigen::Matrix<double,18,1> q = u;

      // function evaluation
      ShellElementStressWRTDisplacementSensitivity<double> foo(dconst,iconst);
      Eigen::Matrix<double,3,1> S = foo(q, 0);
      if(verbose) std::cerr << "#1 S = " << S.transpose() << std::endl;

      // Jacobian evaluation using hand-coded routine
      Eigen::Matrix<double,3,18> J;
      elem.andesvmsWRTdisp(0, nu, x.data(), y.data(), z.data(), u.data(), J.data(), ctyp, mat, surface, sflg, ndtemps.data());
      if(verbose) std::cerr << "#2 J = " << J.transpose() << std::endl;

      // Jacobian evaluation by automatic differentiation
      Simo::Jacobian<double,ShellElementStressWRTDisplacementSensitivity> dSdh(dconst,iconst);
      Eigen::Matrix<double,3,18> J_ad = dSdh(q, 0);
      double error_ad = (J.norm() == 0) ? (J-J_ad).norm() : (J-J_ad).norm()/J.norm();
      if(verbose) std::cerr << "#2 J_ad = " << J_ad.transpose() << std::endl;
      if(verbose || error_ad > 1e-12) std::cerr << "#2 error_ad = " << error_ad << std::endl;
      max_error_ad = std::max(max_error_ad, error_ad);

      // Jacobian approximation by finite difference
      Simo::SpatialView<double,ShellElementStressWRTDisplacementSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementStressWRTDisplacementSensitivity>,Eigen::Central> cd(sv, eps);
      Eigen::Matrix<double,3,18> J_fd;
      cd.df(q, J_fd);
      double error_fd = (J.norm() == 0) ? (J-J_fd).norm() : (J-J_fd).norm()/J.norm();
      if(verbose) std::cerr << "#2 J_fd = " << J_fd.transpose() << std::endl;
      if(verbose) std::cerr << "#2 error_fd = " << error_fd << std::endl;
      max_error_fd = std::max(max_error_fd, error_fd);
    }

    // #3 test gravity force wrt thickness sensitivity
    {
      // real parameters
      Eigen::Array<double,14,1> dconst;
      dconst << g, x, y, z, mat->GetAreaDensity()/mat->GetShellThickness(), nsm;

      // integer parameters
      Eigen::Array<int,1,1> iconst;
      iconst << gravflg;

      // inputs
      Eigen::Matrix<double,1,1> q;
      q[0] = mat->GetShellThickness();

      // function evaluation
      ShellElementGravityForceWRTThicknessSensitivity<double> foo(dconst,iconst);
      Eigen::Matrix<double,18,1> S = foo(q, 0);
      if(verbose) std::cerr << "#3 S = " << S.transpose() << std::endl;

      // Jacobian evaluation using hand-coded routine
      Eigen::Matrix<double,18,1> J;
      elem.andesgf(0, x.data(), y.data(), z.data(), J.data(), g.data(), gravflg, mat->GetAreaDensity()-nsm);
      J /= mat->GetShellThickness();
      if(verbose) std::cerr << "#3 J = " << J.transpose() << std::endl;

      // Jacobian evaluation by automatic differentiation
      Simo::Jacobian<double,ShellElementGravityForceWRTThicknessSensitivity> dSdh(dconst,iconst);
      Eigen::Matrix<double,18,1> J_ad = dSdh(q, 0);
      double error_ad = (J.norm() == 0) ? (J-J_ad).norm() : (J-J_ad).norm()/J.norm();
      if(verbose) std::cerr << "#3 J_ad = " << J_ad.transpose() << std::endl;
      if(verbose || error_ad > 1e-12) std::cerr << "#3 error_ad = " << error_ad << std::endl;
      max_error_ad = std::max(max_error_ad, error_ad);

      // Jacobian approximation by finite difference
      Simo::SpatialView<double,ShellElementGravityForceWRTThicknessSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementGravityForceWRTThicknessSensitivity>,Eigen::Central> cd(sv, eps);
      Eigen::Matrix<double,18,1> J_fd;
      cd.df(q, J_fd);
      double error_fd = (J.norm() == 0) ? (J-J_fd).norm() : (J-J_fd).norm()/J.norm();
      if(verbose) std::cerr << "#3 J_fd = " << J_fd.transpose() << std::endl;
      if(verbose) std::cerr << "#3 error_fd = " << error_fd << std::endl;
      max_error_fd = std::max(max_error_fd, error_fd);
    }

    // #4 test stiffness wrt thickness sensitivity
    {
      // real parameters
      Eigen::Array<double,68,1> dconst;
      dconst << x, y, z, E, nu, rho, eframe, coefs, Ta, W, ndtemps;

      // integer parameters
      Eigen::Array<int,3,1> iconst;
      iconst << ctyp, flg, tflg;

      // inputs
      Eigen::Matrix<double,1,1> q;
      q[0] = mat->GetShellThickness();

      // function evaluation
      ShellElementStiffnessWRTThicknessSensitivity<double> foo(dconst,iconst);
      Eigen::Matrix<double,18,18> S = foo(q, 0);
      if(verbose) std::cerr << "#4 S =\n" << S << std::endl;

      // Jacobian evaluation using hand-coded routine
      Eigen::Matrix<double,18,18> J;
      elem.andesstfWRTthick(0, J.data(), nu, x.data(), y.data(), z.data(), ctyp, mat, flg, tflg, ndtemps.data());
      if(verbose) std::cerr << "#4 J = \n" << J << std::endl;

      // Jacobian evaluation by automatic differentiation
      Simo::FirstPartialSpaceDerivatives<double,ShellElementStiffnessWRTThicknessSensitivity> dSdh(dconst,iconst);
      Eigen::Matrix<double,18,18> J_ad = dSdh(q, 0)[0];
      double error_ad = (J.norm() == 0) ? (J-J_ad).norm() : (J-J_ad).norm()/J.norm();
      if(verbose) std::cerr << "#4 J_ad = \n" << J_ad << std::endl;
      if(verbose || error_ad > 1e-12) std::cerr << "#4 error_ad = " << error_ad << std::endl;
      max_error_ad = std::max(max_error_ad, error_ad);

      // Jacobian approximation by finite difference
      Simo::SpatialView<double,ShellElementStiffnessWRTThicknessSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementStiffnessWRTThicknessSensitivity>,Eigen::Central> cd(sv, eps);
      Eigen::Matrix<double,324,1> j_fd;
      cd.df(q, j_fd);
      Eigen::Matrix<double,18,18> J_fd = Eigen::Map<Eigen::Matrix<double,18,18>>(j_fd.data());
      double error_fd = (J.norm() == 0) ? (J-J_fd).norm() : (J-J_fd).norm()/J.norm();
      if(verbose) std::cerr << "#4 J_fd = \n" << J_ad << std::endl;
      if(verbose) std::cerr << "#4 error_fd = " << error_fd << std::endl;
      max_error_fd = std::max(max_error_fd, error_fd);
    }

  }
  std::cerr << "max_error_ad = " << max_error_ad << std::endl;
  std::cerr << "max_error_fd = " << max_error_fd << std::endl;
}

