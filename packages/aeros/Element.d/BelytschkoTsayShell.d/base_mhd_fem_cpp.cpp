#ifdef USE_EIGEN3
#include <Element.d/BelytschkoTsayShell.d/base_mhd_fem.h>
#include <Eigen/Geometry>
#include <iostream>

void 
rotprojbt1(const Eigen::Matrix<double,3,4> &ecurn, double *_evec0,
           double *_evec) {
  //=======================================================================
  //  rotprojbt1 = drilling rotation projection
  //
  //  reference: Belytschko, Ted, and Itai Leviathan. "Projection schemes for one-point quadrature
  //  shell elements." Computer methods in applied mechanics and engineering 115.3 (1994): 277-286.
  //
  //  arguments description
  //  ---------------------
  //  input
  //  -----
  //  ecurn(3,4) : element nodal positions
  //
  //  evec0(24,1) : force or velocity vector
  //
  //  output:
  //  ------
  //  evec(24,1) : force or velocity vector
  //
  // ======================================================================

  // ====================================
  // subroutine argument
  // ===================
  Eigen::Map<Eigen::Matrix<double,24,1> > evec0(_evec0);
  Eigen::Map<Eigen::Matrix<double,24,1> > evec(_evec);

  // ====================================
  // local variable
  // ==============
  Eigen::Matrix<double,3,3> rotpmat;

  int i0,i3;
  // ====================================

  // loop over nodes
  for(int inode=0; inode<4; ++inode) {

     // address
     i0 = 6*inode;
     i3 = i0+3;

     // -----------------------------------------------------
     // force
     evec.segment<3>(i0) = evec0.segment<3>(i0); // f_x, f_y, f_z

     // -----------------------------------------------------
     // moment

     // get rotation projection matrix
     getrotpmatbt(inode, ecurn, rotpmat);
        // input : inode,ecurn
        // output : rotpmat

     // rotation projection:
     evec.segment<3>(i3) = rotpmat*evec0.segment<3>(i3); // m_x, m_y, m_z

  }

  return;
}

void 
rotprojbt2(const Eigen::Matrix<double,3,4> &ecurn, double *_evec0,
           double *_evec) {
  //=======================================================================
  //  rotprojbt2 = Coupled projection of rigid body rotation and drill rotation
  //
  //  reference: Belytschko, Ted, and Itai Leviathan. "Projection schemes for one-point quadrature
  //  shell elements." Computer methods in applied mechanics and engineering 115.3 (1994): 277-286.
  //
  //  arguments description
  //  ---------------------
  //  input
  //  -----
  //  ecurn(3,4) : element nodal positions
  //
  //  evec0(24,1) : force or velocity vector
  //
  //  output:
  //  ------
  //  evec(24,1) : force or velocity vector
  //
  // ======================================================================

  // ====================================
  // subroutine argument
  // ===================
  Eigen::Map<Eigen::Matrix<double,24,1> > evec0(_evec0);
  Eigen::Map<Eigen::Matrix<double,24,1> > evec(_evec);

  // ====================================
  // local variable
  // ==============
  Eigen::Matrix<double,3,3> rotpmat;
  Eigen::Matrix<double,3,1> nI;
  Eigen::Matrix<double,24,7> R;
  Eigen::Matrix<double,7,7> RtR;

  int i0,i3;
  // ====================================

  // loop over nodes
  for(int inode=0; inode<4; ++inode) {

     // get rotation projection matrix and node normal
     getrotpmatnbt(inode, ecurn, rotpmat, nI);
       // input : inode,ecurn
       // output : rotpmat,n

     // rigid body modes
     const double &xI = ecurn(0,inode), &yI = ecurn(1,inode), &zI = ecurn(2,inode);
     R.block<6,3>(6*inode,0) <<  0,   zI, -yI,
                                -zI,   0,  xI,
                                 yI, -xI,   0,
                                  1,   0,   0,
                                  0,   1,   0,
                                  0,   0,   1;
     // drilling modes
     R.block<6,4>(6*inode,3).setZero();
     R.block<3,1>(6*inode+3,3+inode) = nI;
  }

  RtR = R.transpose()*R;
  evec = evec0 - R*RtR.ldlt().solve(R.transpose()*evec0);

  return;
}

void
getrotpmatbt(int inode, const Eigen::Matrix<double,3,4> &ecurn,
             Eigen::Matrix<double,3,3> &rotpmat) {
  //=======================================================================
  //  getrotpmatbt = compute 3x3 rotation projection matrix
  //
  //                note:
  //                ----
  //                theta^bar = theta - ( theta . e3 ) e3
  //
  //  arguments description
  //  ---------------------
  //  input
  //  -----
  //  inode : current node number
  //
  //  ecurn(3,4) : current element nodal coordinate
  //
  //  output:
  //  ------
  //  rotpmat(3,3) : rotation projection matrix
  //
  // ======================================================================

  // ====================================
  // local variable
  // ==============
  Eigen::Matrix<int,4,3> ijkndx;

  Eigen::Matrix<double,3,1> nodi, nodj, nodk;
  Eigen::Matrix<double,3,1> e1vec, e2vec, e3vec;
  // ====================================

  // initialize
  ijkndx << 0,1,3, 1,2,0, 2,3,1, 3,0,2;

  //   k
  //    o
  //    |
  //    |
  //    |
  // e2 ^
  //    |
  //    |
  //    |
  //    o-------->--------o
  //   i        e1         j
  
  // define node
  nodi = ecurn.col(ijkndx(inode,0)); // note: inode should be 0,1,2,3
  nodj = ecurn.col(ijkndx(inode,1));
  nodk = ecurn.col(ijkndx(inode,2));

  // element edge vector: e1, e2
  e1vec = (nodj - nodi);
  e2vec = (nodk - nodi);

  // e3 vector: e3= e1 x e2 / ||e1 x e2||
  // ---------
  e3vec = e1vec.cross(e2vec).normalized();

  // set rotation projection matrix
  rotpmat(0,0) = 1.0 - e3vec[0]*e3vec[0];
  rotpmat(0,1) =     - e3vec[0]*e3vec[1];
  rotpmat(0,2) =     - e3vec[0]*e3vec[2];

  rotpmat(1,0) =     - e3vec[1]*e3vec[0];
  rotpmat(1,1) = 1.0 - e3vec[1]*e3vec[1];
  rotpmat(1,2) =     - e3vec[1]*e3vec[2];

  rotpmat(2,0) =     - e3vec[2]*e3vec[0];
  rotpmat(2,1) =     - e3vec[2]*e3vec[1];
  rotpmat(2,2) = 1.0 - e3vec[2]*e3vec[2];

  return;
}

void
getrotpmatnbt(int inode, const Eigen::Matrix<double,3,4> &ecurn,
              Eigen::Matrix<double,3,3> &rotpmat, Eigen::Matrix<double,3,1> &e3vec) {
  //=======================================================================
  //  getrotpmatnbt = compute 3x3 rotation projection matrix and node normal
  //
  //                note:
  //                ----
  //                theta^bar = theta - ( theta . e3 ) e3
  //
  //  arguments description
  //  ---------------------
  //  input
  //  -----
  //  inode : current node number
  //
  //  ecurn(3,4) : current element nodal coordinate
  //
  //  output:
  //  ------
  //  rotpmat(3,3) : rotation projection matrix
  //
  // ======================================================================

  // ====================================
  // local variable
  // ==============
  Eigen::Matrix<int,4,3> ijkndx;

  Eigen::Matrix<double,3,1> nodi, nodj, nodk;
  Eigen::Matrix<double,3,1> e1vec, e2vec;
  // ====================================

  // initialize
  ijkndx << 0,1,3, 1,2,0, 2,3,1, 3,0,2;

  //   k
  //    o
  //    |
  //    |
  //    |
  // e2 ^
  //    |
  //    |
  //    |
  //    o-------->--------o
  //   i        e1         j
  
  // define node
  nodi = ecurn.col(ijkndx(inode,0)); // note: inode should be 0,1,2,3
  nodj = ecurn.col(ijkndx(inode,1));
  nodk = ecurn.col(ijkndx(inode,2));

  // element edge vector: e1, e2
  e1vec = (nodj - nodi);
  e2vec = (nodk - nodi);

  // e3 vector: e3= e1 x e2 / ||e1 x e2||
  // ---------
  e3vec = e1vec.cross(e2vec).normalized();

  // set rotation projection matrix
  rotpmat(0,0) = 1.0 - e3vec[0]*e3vec[0];
  rotpmat(0,1) =     - e3vec[0]*e3vec[1];
  rotpmat(0,2) =     - e3vec[0]*e3vec[2];

  rotpmat(1,0) =     - e3vec[1]*e3vec[0];
  rotpmat(1,1) = 1.0 - e3vec[1]*e3vec[1];
  rotpmat(1,2) =     - e3vec[1]*e3vec[2];

  rotpmat(2,0) =     - e3vec[2]*e3vec[0];
  rotpmat(2,1) =     - e3vec[2]*e3vec[1];
  rotpmat(2,2) = 1.0 - e3vec[2]*e3vec[2];

  return;
}
#endif
