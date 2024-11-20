// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#ifndef ContactAlgorithm_h
#define ContactAlgorithm_h

#include <cmath>
#include <iostream>
#include "Contact_Defines.h"

namespace acme {

static const Real colinearity_tolerance = 0.99;

inline Real colinearityTolerance() 
{
  return colinearity_tolerance;
}

static const Real determinant_tolerance = 1.0e-30;

inline Real determinantTolerance() {
  return determinant_tolerance;
};

/*!
 * Compute the cross product of a and b.
 */
template<typename DataType>
inline void Cross(const DataType a[3], const DataType b[3], DataType c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/*!
 * Cross product where the first 3 real values represent the vector a
 * and the second set of three real values represents the vector b.
 */
template<typename DataType>
inline void Cross(const DataType dx1, const DataType dy1, const DataType dz1,
                  const DataType dx2, const DataType dy2, const DataType dz2,
		  DataType normal[3])
{
  normal[0] = dy1*dz2-dz1*dy2;
  normal[1] = dz1*dx2-dx1*dz2;
  normal[2] = dx1*dy2-dy1*dx2;
}

/*!
 * Compute the scalar triple product = (a X b).c
 */
template<typename DataType>
inline DataType ScalarTripleProduct(const DataType a[3], const DataType b[3], const DataType c[3])
{
  return  (a[1]*b[2] - a[2]*b[1])*c[0]
        + (a[2]*b[0] - a[0]*b[2])*c[1]
        + (a[0]*b[1] - a[1]*b[0])*c[2];
}

namespace optimized {

  /*!
   * Compute the cross product of a and b.
   */
  template<typename DataType>
  inline void Cross(const DataType a[3], const DataType b[3], DataType c[3])
  {
    const DataType u1 = a[0];
    const DataType u2 = a[1];
    const DataType u3 = a[2];

    const DataType v1 = b[0];
    const DataType v2 = b[1];
    const DataType v3 = b[2];

    const DataType t1 = u1 - u2;
    const DataType t2 = v2 + v3;
    const DataType t3 = u1*v3;
    const DataType t4 = t1*t2 - t3;

    c[0] = v2*(t1-u3) - t4;
    c[1] = u3*v1 - t3;
    c[2] = t4 - u2 *(v1 - t2);
  }

} // end namespace optimized

template<typename DataType>
inline bool isZero(const DataType a[3])
{
  return (   a[0] == 0.0
          && a[1] == 0.0
	  && a[2] == 0.0);
}

template<typename DataType>
inline void Zero(DataType a[3])
{
  a[0] = 0.0;
  a[1] = 0.0;
  a[2] = 0.0;
}

template<typename DataType>
inline void Scale(DataType a[3], const DataType scale)
{
  a[0] *= scale;
  a[1] *= scale;
  a[2] *= scale;
}

template<typename DataType>
inline void Copy(const DataType source[3], DataType destination[3])
{
  destination[0] = source[0];
  destination[1] = source[1];
  destination[2] = source[2];
}

/*!
 * Compute the sum of all elements of a 3D vector.
 */
template<typename DataType>
inline DataType Sum(const DataType a[3])
{
  return (a[0]+a[1]+a[2]);
}

/*!
 * Compute the sum of the absolute value of all elements of a 3D vector.
 */
template<typename DataType>
inline DataType AbsSum(const DataType a[3])
{
  return (std::fabs(a[0])+std::fabs(a[1])+std::fabs(a[2]));
}

/*!
 * Compute the dot product.
 */
template<typename DataType>
inline DataType Dot(const DataType a[3])
{
  return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

/*!
 * Compute the dot product between two 3D vectors.
 */
template<typename DataType>
inline DataType Dot(const DataType a[3], const DataType b[3])
{
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

/*!
 * Compute the magnitude of a 3D vector.
 */
template<typename DataType>
inline DataType Magnitude(const DataType a[3])
{
  using std::sqrt;
  return sqrt(Dot(a));
}

/*!
 * Compute the squared magnitude of a 3D vector.
 */
template<typename DataType>
inline DataType Magnitude2(const DataType a[3])
{
  return Dot(a);
}

/*!
 * Compute the distance between vectors a and b.
 */
template<typename DataType>
inline DataType Distance(const DataType a[3], const DataType b[3])
{
  return std::sqrt((b[0]-a[0])*(b[0]-a[0]) +
              (b[1]-a[1])*(b[1]-a[1]) +
              (b[2]-a[2])*(b[2]-a[2]));
}

/*!
 * Compute the squared distance between vectors a and b.
 */
template<typename DataType>
inline DataType Distance2(const DataType a[3], const DataType b[3])
{
  return ((b[0]-a[0])*(b[0]-a[0]) +
          (b[1]-a[1])*(b[1]-a[1]) +
          (b[2]-a[2])*(b[2]-a[2]));
}

/*!
 * Normalize a 3D vector. v = v /||v||
 *
 * Return the magnitude of the input vector.
 */
template<typename DataType>
inline DataType Normalize(DataType v[3])
{
  DataType mag = Magnitude(v);
  if ( mag != 0.0 ) {
    DataType invMag = 1.0 / mag;
    v[0] *= invMag;
    v[1] *= invMag;
    v[2] *= invMag;
  }
  return mag;
}

/*!
 * Compute the area of the triangle formed from the 2 vectors a and b.
 */
template<typename DataType>
inline DataType Area(const DataType a[3], const DataType b[3])
{
  DataType c[3];
  Cross(a, b, c);
  return 0.5*Magnitude(c);
}

template<typename DataType>
inline DataType det2x2(const DataType A[2][2])
{
  // A = | a b |
  //   = | c d |
  return (A[0][0]*A[1][1] - A[1][0]*A[0][1]);
}

template<typename DataType>
inline DataType det2x2(DataType a, DataType b,
                   DataType c, DataType d)
{
  // A = | a b |
  //   = | c d |
  return (a*d - b*c);
}

template<typename DataType>
inline void Solve_2x2_Matrix(const DataType A[2][2], const DataType x[2], DataType y[2])
{
  DataType a = A[0][0], b = A[0][1];
  DataType c = A[1][0], d = A[1][1];

  DataType det = det2x2(a,b,c,d);

  if ( det != 0.0 ) {
    DataType inv_det = 1.0 / det;
    DataType x0_det = x[0] * d - x[1] * b;
    DataType x1_det = a * x[0] - c * x[1];
    y[0] = x0_det*inv_det;
    y[1] = x1_det*inv_det;
  } else {
    // error! singular matrix.
  }
}

template<typename DataType>
inline void MatrixVectorProduct(const DataType M[3][3], const DataType x[3], DataType y[3])
{ 
  DataType M11 = M[0][0], M12 = M[0][1], M13 = M[0][2];
  DataType M21 = M[1][0], M22 = M[1][1], M23 = M[1][2];
  DataType M31 = M[2][0], M32 = M[2][1], M33 = M[2][2];

  DataType x1 = x[0];
  DataType x2 = x[1];
  DataType x3 = x[2];

  y[0] = M11*x1 + M12*x2 + M13*x3;
  y[1] = M21*x1 + M22*x2 + M23*x3; 
  y[2] = M31*x1 + M32*x2 + M33*x3;
}

template<typename DataType>
inline void SymMatrixVectorProduct(const DataType M[3][3], const DataType x[3], DataType y[3])
{ 
  DataType M11 = M[0][0], M12 = M[0][1], M13 = M[0][2];
  DataType                M22 = M[1][1], M23 = M[1][2];
  DataType                               M33 = M[2][2];

  DataType x1 = x[0], x2 = x[1], x3 = x[2];

  y[0] = M11*x1 + M12*x2 + M13*x3;
  y[1] = M12*x1 + M22*x2 + M23*x3; 
  y[2] = M13*x1 + M23*x2 + M33*x3;
}

template<typename DataType>
inline DataType Invert_3x3_Sym_Matrix(const DataType M [3][3], DataType M_inv [3][3], const char *msg="")
{
  DataType M11 = M[0][0], M12 = M[0][1], M13 = M[0][2];
  DataType M21 = M[1][0], M22 = M[1][1], M23 = M[1][2];
  DataType M31 = M[2][0], M32 = M[2][1], M33 = M[2][2];
  //
  // Cofactors
  //
  DataType A11 = M22*M33 - M23*M32;
  DataType A12 = M13*M32 - M12*M33;
  DataType A13 = M12*M23 - M13*M22;
  DataType A22 = M11*M33 - M13*M31;
  DataType A23 = M13*M21 - M11*M23;
  DataType A33 = M11*M22 - M12*M21;
  //
  // Determinant
  //
  DataType det = M11*A11 + M12*A12 + M13*A13 ;
  if( std::fabs(det) < determinant_tolerance ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }
  DataType det_inv = 1.0/det;
  M_inv[0][0] = A11*det_inv ;
  M_inv[1][0] = A12*det_inv ;
  M_inv[2][0] = A13*det_inv ;

  M_inv[0][1] = M_inv[1][0];
  M_inv[1][1] = A22*det_inv ;
  M_inv[2][1] = A23*det_inv ;

  M_inv[0][2] = M_inv[2][0];
  M_inv[1][2] = M_inv[2][1];
  M_inv[2][2] = A33*det_inv ;

  return det;
}

template<typename DataType>
inline bool isEqual(const DataType x, const DataType y)
{
  const double epsilon = 1e-5; /* some small number such as 1e-5 */;
  return std::abs(x - y) <= epsilon * std::abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

template<typename DataType>
inline bool isNormal(const DataType a[3])
{
  //
  // KHP: Will floating point precision haunt me on this comparison?
  //
  // return ( std::fabs(Magnitude(a) - 1.0 ) < 1.0e-8 );
  //
  return isEqual( Magnitude(a), 1.0 );
}

template<typename DataType>
inline void print(const DataType a[3], const char *name="a")
{
  std::cout << name << "(" << a[0] << "," << a[1] << "," << a[2] << ")\n";
}

template<typename DataType>
inline void printMag(const DataType a[3], const char *name="a")
{
  std::cout << "(" << a[0] << "," << a[1] << "," << a[2] << "), ||" << name << "||= " << Magnitude(a) << std::endl;
}

template<typename DataType>
inline DataType Invert_2x2_Matrix(const DataType M[2][2], DataType M_inv[2][2], const char *msg="")
{
  DataType a = M[0][0], b = M[0][1];
  DataType c = M[1][0], d = M[1][1];

  DataType det = a*d - b*c;
  if( std::fabs(det) < determinant_tolerance) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }

  DataType inv_det = 1.0 / det;

  M_inv[0][0] =  d * inv_det; M_inv[0][1] = -b * inv_det;
  M_inv[1][0] = -c * inv_det; M_inv[1][1] =  a * inv_det;

  return det;
}

template<typename DataType>
inline DataType Invert_2x2_Sym_Matrix(const DataType M[2][2], DataType M_inv[2][2], const char *msg="")
{
  //
  // if ( M[0][1] != M[1][0] ) error, should have called regular
  // Invert_2x2_Matrix method.
  //
  DataType a = M[0][0], b = M[0][1];
  DataType              d = M[1][1];

  DataType det = a*d - b*b;
  if( std::fabs(det) < determinant_tolerance ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }

  DataType inv_det = 1.0 / det;

  M_inv[0][0] =  d * inv_det; M_inv[0][1] = -b * inv_det;
  M_inv[1][0] =  M_inv[0][1]; M_inv[1][1] =  a * inv_det;

  return det;
}

template<typename DataType>
inline DataType Invert_3x3Matrix(const DataType M [3][3], DataType M_inv [3][3], const char *msg="")
{
  //
  // Cofactors
  //
  DataType cofM [3] [3];
  cofM[0][0] =   M[1][1]*M[2][2] - M[1][2]*M[2][1] ;
  cofM[1][0] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) ;
  cofM[2][0] =   M[0][1]*M[1][2] - M[0][2]*M[1][1] ;
  cofM[0][1] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) ;
  cofM[1][1] =   M[0][0]*M[2][2] - M[0][2]*M[2][0] ;
  cofM[2][1] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) ;
  cofM[0][2] =   M[1][0]*M[2][1] - M[1][1]*M[2][0] ;
  cofM[1][2] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) ;
  cofM[2][2] =   M[0][0]*M[1][1] - M[0][1]*M[1][0] ;
  //
  // Determinant
  //
  DataType det = M[0][0]*cofM[0][0] + M[0][1]*cofM[0][1] + M[0][2]*cofM[0][2] ;
  if( std::fabs(det) < determinant_tolerance ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }
  DataType det_inv = 1.0/det  ;
  M_inv[0][0] = cofM[0][0]*det_inv ;
  M_inv[1][0] = cofM[0][1]*det_inv ;
  M_inv[2][0] = cofM[0][2]*det_inv ;
  M_inv[0][1] = cofM[1][0]*det_inv ;
  M_inv[1][1] = cofM[1][1]*det_inv ;
  M_inv[2][1] = cofM[1][2]*det_inv ;
  M_inv[0][2] = cofM[2][0]*det_inv ;
  M_inv[1][2] = cofM[2][1]*det_inv ;
  M_inv[2][2] = cofM[2][2]*det_inv ;
  return det;
}

//
//  Compute ACME edge curvature from a pair of connected face normals
//  Curvature is an artifical construct where 0.0 is flat, and 2.0 is
//  fully opposed as:
//
//    curv = 1.0 - Theta = 1.0 - norm1.norm2
//
template<typename DataType>
inline DataType ComputeCurvatureFromNormals(const DataType *norm1, const DataType *norm2) {
  return (1.0 - Dot(norm1, norm2));
}




} // end namespace acme

#endif
