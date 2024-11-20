#ifndef _SHELLMATERIAL_CPP_
#define _SHELLMATERIAL_CPP_

#ifdef USE_EIGEN3
#include <cmath>
#include <stdexcept>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <Element.d/Function.d/AutoDiffScalarPlugin.h>

template<typename doublereal>
Eigen::Matrix<doublereal,3,3>
ShellMaterial<doublereal>::andesinvt(doublereal *_eframe, doublereal *_aframe, doublereal thetaf)
{
  // Initialized data 
  doublereal zero = 0.;
  doublereal one = 1.;
  Eigen::Matrix<doublereal,3,1> r; r << 1, 1, 2;

  // Builtin functions 
  using std::atan;
  using std::cos;
  using std::sin;

  // Local variables 
  int i;
  doublereal norm1, norm2, normref, proj1, proj2;
  doublereal cosine1, cosine2;
  doublereal theta, thetad, pi, twopi;
  doublereal costheta, sintheta;
  Eigen::Map<Eigen::Matrix<doublereal,3,3> > eframe(_eframe), aframe(_aframe);
  Eigen::Matrix<doublereal,3,1> refvec, orifiber;
  Eigen::Matrix<doublereal,3,3> T, invT;

  try {

    pi = acos(-1.);
    twopi = pi * 2.;

// .....SET THE REFERENCE VECTOR FOR ORIENTING THE FIBERS OF THE LAYER 
// .....(ALWAYS TAKE THE FIRST VECTOR OF THE FRAME) 

    refvec = aframe.col(0);

// .....PROJECT THE REFERENCE VECTOR INTO THE PLANE OF 
// .....THE ELEMENT TO GET THE FIBER ORIENTATION VECTOR 

    norm1 = eframe.col(0).norm();
    norm2 = eframe.col(1).norm();
    normref = refvec.norm();
    proj1 = eframe.col(0).dot(refvec);
    proj2 = eframe.col(1).dot(refvec);

    if (normref == zero) {
        cosine1 = one;
        cosine2 = zero;
    } else {
        cosine1 = proj1 / (norm1 * normref);
        cosine2 = proj2 / (norm2 * normref);
    }

    for (i = 0; i < 3; ++i) {
        orifiber[i] = cosine1 * eframe(i, 0) + cosine2 * eframe(i, 1);
    }

// .....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL 
// .....COORDINATE SYSTEM TO THE REFERENCE ORIENTATION VECTOR 

    proj1 = eframe.col(0).dot(orifiber);
    proj2 = eframe.col(1).dot(orifiber);
    thetad = zero;

    if (proj1 == zero) {
        if (proj2 == zero) {
            throw 700;
        }
        if (proj2 > zero) {
            thetad = pi * .5;
        }
        if (proj2 < zero) {
            thetad = pi * 1.5;
        }
    } else {
        if (proj2 == zero) {
            if (proj1 == zero) {
                throw 700;
            }
            if (proj1 > zero) {
                thetad = zero;
            }
            if (proj1 < zero) {
                thetad = pi;
            }
        } else {
            thetad = atan(proj2 / proj1);
        }
    }

    if (thetad < zero) {
        thetad += twopi;
    }

    if (thetad < zero || thetad > twopi) {
        throw 800;
    }

// .....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL 
// .....COORDINATE SYSTEM TO THE DIRECTION OF THE FIBER 

    theta = thetad + thetaf;

    if (theta > twopi) {
        theta -= twopi;
    }

// .....INITIALIZE THE ROTATION MATRIX FROM THE FIBER COORDINATE 
// .....SYSTEM {1;2} TO THE ELEMENT TRIANGULAR SYSTEM {x;y} 

    costheta = cos(theta);
    sintheta = sin(theta);

    T(0, 0) = costheta * costheta;
    T(0, 1) = sintheta * sintheta;
    T(0, 2) = costheta * 2. * sintheta;

    T(1, 0) = sintheta * sintheta;
    T(1, 1) = costheta * costheta;
    T(1, 2) = costheta * -2. * sintheta;

    T(2, 0) = -costheta * sintheta;
    T(2, 1) = costheta * sintheta;
    T(2, 2) = costheta * costheta - sintheta * sintheta;

// .....COMPUTE THE INVERSE OF [T]: 
// .....[invT] = inverse(diag[R]) * [T]^t * diag[R] 

    invT.noalias() = r.asDiagonal().inverse()*T.transpose()*r.asDiagonal();

  }

  catch (int err) {

    switch(err) {

// .....ERROR-MESSAGE IF THE REFERENCE ORIENTATION IS BUGGY 
      case 700:
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellMaterial::andesinvt ***\n"
          "*** The reference orientation vector        ***\n"
          "*** is parallel to the two in-plane         ***\n"
          "*** and orthogonal local frames!            ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE          ***\n");
        break;

// .....ERROR-MESSAGE IF THE ORIENTATION ANGLE IS OUT-OF-BOUNDS 
      case 800:
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellMaterial::andesinvt ***\n"
          "*** The angle from the local [x]            ***\n"
          "*** axis of the triangular coordinate       ***\n"
          "*** system to the reference direction       ***\n"
          "*** is out-of-bounds: it must be            ***\n"
          "*** within the range 0-2pi radians          ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE          ***\n");
        break;
    }
  }

  return invT;
}

template
Eigen::Matrix<double,3,3>
ShellMaterial<double>
::andesinvt(double *, double *, double);

#endif

#endif
