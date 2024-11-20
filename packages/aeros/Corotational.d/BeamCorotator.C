#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/BLAS.h>
#include <Element.d/Element.h>
#include <Corotational.d/BeamCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Utils.d/linkfc.h>

BeamCorotator::BeamCorotator(int _n1, int _n2, double z[3], 
                             FullSquareMatrix &_origK,
                             int fitAlgBeam)
{
 // Set Node Numbers
 n1 = _n1;
 n2 = _n2;

 // Copy the original stiffness matrix
 int i,j;
 for(i=0; i<12; ++i)
   for(j=0; j<12; ++j)
     origK[i][j] = _origK[i][j];

 // Copy the element normal vector
 zVec =  z;

 // Set the fit algorithm
 fitAlg  = fitAlgBeam;
}

void
BeamCorotator::getStiffAndForce(GeomState &geomState, CoordSet &cs, 
                                FullSquareMatrix &elK, double *f, double dt, double t)
{
 int i, j, inode;

 // Get Nodes original coordinates
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Get Nodes current coordinates
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 double xl0[2][3],xln[2][3],t0[3][3],t0n[3][3],vld[12],locF[12],zVecL[2][3];
 
 // Extract deformational displacement
 extractDefDisp(node1,node2,ns1,ns2,zVecL,xl0,xln,t0,t0n,vld);

 // copy origK (original stiffness matrix) to elK
  for(i=0; i<12; ++i)
   for(j=0; j<12; ++j)
     elK[i][j] = origK[i][j];

 // compute locF (unprojected internal local Force) as origK*vld
 _FORTRAN(dgemv)('N',12,12,1.0,(double *)elK.data(),12,vld,1,0.0,locF,1);


 // Compute gradients of nodal deformational pseudorotations
 // and correct stiffness and force

 double rotvar[2][3][3];
 
 for(inode = 0; inode < 2; ++inode)
   pseudorot_var(vld+inode*6+3, rotvar[inode]);

  leftmult_rotvar( 2, 1, rotvar, elK);
 rightmult_rotvar( 2, 0, rotvar, elK);

 double fe[12];

 for(inode = 0; inode < 2; ++inode)
   for(i = 0; i < 3; ++i) {
     fe[6*inode+i]   = locF[6*inode+i];
     fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                       rotvar[inode][1][i]*locF[6*inode+4] +
                       rotvar[inode][2][i]*locF[6*inode+5];
   }

 // Add second variation of pseudorotations contracted with the
 // nodal moment to the diagonal blocks of the stiffness

 for(inode = 0; inode < 2; ++inode) {
   pseudorot_2var(vld+inode*6+3, locF+inode*6+3, rotvar[inode]);
   for(i = 0; i < 3; ++i)
     for(j = 0; j < 3; ++j)
       elK[6*inode+3+i][6*inode+3+j] += rotvar[inode][i][j];
 }


 // Compute nonlinear projector matrix relative to deformed element
 // and correct stiffness and force
 // pmat = projector matrix
 double pmat[12][12];
 double gmat[3][12];
 gradDefDisp( zVecL, xln, pmat, gmat);

 double scrstiff[12][12];

 // Form: [K] = [P'][K]

 _FORTRAN(dgemm)('N','T',12,12,12,1.0,elK.data(),12,
                   (double*)pmat,12,0.0,(double*)scrstiff,12);
		   
 // Form:  [K] = [K*][P]

 _FORTRAN(dgemm)('N','N',12,12,12,1.0,(double*)pmat,12,
                   (double*)scrstiff,12,0.0,elK.data(),12);		   

 // Form: {f} = [P']{fe}

  _FORTRAN(dgemv)('N',12,12,1.0,(double *)pmat,12,fe,1,0.0,f,1);

 // Form geometric stiffness from internal force and material stiffness

 double stiffGeo1[12][12], stiffGeo2[12][12];

 // corotStiffGeo(zVecL,xln,pmat,f,stiffGeo1,stiffGeo2);
 formCorrectGeometricStiffness(pmat, gmat, f, stiffGeo1, 
                               stiffGeo2, fe, node1, node2, ns1, ns2);

 // Sum geometric stiffness contributions into elK first

 for(i=0; i<12; ++i)
   for(j=0; j<12; ++j)
     elK[i][j] += stiffGeo1[i][j] + stiffGeo2[i][j];

 // Transform internal force and stiffness matrix to global coordinate system 

 tran_fsl(f,elK,t0n,2);

 // The skew symmetric load stiffness matrix due to axial external moments is
 // added separately (see Domain::getFollowerForce in Driver.d/NLStatic.C)
 // For now, it can be removed from elK by symmetrizing to avoid doubling.
 elK.symmetrize();
}

//----------------------------------------------------------------------

void 
BeamCorotator::getInternalForce(GeomState &geomState, CoordSet &cs, 
                                FullSquareMatrix &elK, double *f, double dt, double t)
{
 // Get Nodes original coordinates
 Node &node1 = cs.getNode( n1 );
 Node &node2 = cs.getNode( n2 );

 // Get Nodes current coordinates
 NodeState &ns1 = geomState[ n1 ];
 NodeState &ns2 = geomState[ n2 ];

 double xl0[2][3], xln[2][3], t0[3][3], t0n[3][3], vld[12], locF[12], zVecL[2][3];

 // C0    = initial configuration
 // C0n   = nth configuration
 // xl0   = C0 local coordinates
 // xln   = C0n local coordinates
 // t0n   = transformation matrix between C0n and C0
 // vld   = local deformation vector
 // locF  = local unprojected internal force
 // origK = original stiffness matrix

 // f = T' P' H' K v

 // Extract deformational displacement from C0 to C0n configurations

 extractDefDisp(node1,node2, ns1,ns2, zVecL, xl0,xln, t0,t0n, vld);
 
 // copy origK (original stiffness matrix) to elK

 int i,j;
 for(i=0; i<12; ++i) {
   locF[i] = 0.0;
   for(j=0; j<12; ++j)
     elK[i][j] = origK[i][j];
  }
 
 // compute locF (unprojected internal local Force) as origK*vld

 _FORTRAN(dgemv)('N',12,12,1.0,(double *)elK.data(),12,vld,1,0.0,locF,1);

 // Compute gradients of nodal deformational pseudorotations
 // and correct stiffness and force

 double rotvar[2][3][3];
 
 int inode;
 for(inode=0; inode<2; ++inode)
   pseudorot_var( vld+inode*6+3, rotvar[inode] );
  
 double fe[12];

 for(inode=0; inode<2; ++inode)
   for(i=0; i<3; ++i) {
     fe[6*inode+i]   = locF[6*inode+i];
     fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                       rotvar[inode][1][i]*locF[6*inode+4] +
                       rotvar[inode][2][i]*locF[6*inode+5];
   }

 // Compute nonlinear projector matrix relative to deformed element
 // and correct stiffness and force

 double pmat[12][12], gmat[3][12];

 gradDefDisp(zVecL, xln, pmat, gmat);

 // Form: {f} = [P']{fe}

 _FORTRAN(dgemv)('N',12,12,1.0,(double *)pmat,12,fe,1,0.0,f,1);

 // Transform internal force and stiffness matrix to global coordinate system 

 tran_force(f, t0n, 2);

}				  

//----------------------------------------------------------------------------

void
BeamCorotator::formGeometricStiffness(GeomState &geomState, CoordSet &cs, 
                                FullSquareMatrix &elK, double *f)
{
 int i, j, inode;

 // Get Nodes original coordinates
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Get Nodes current coordinates
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 double xl0[2][3],xln[2][3],t0[3][3],t0n[3][3],vld[12],locF[12],zVecL[2][3];

 // Extract deformational displacement
 extractDefDisp(node1,node2,ns1,ns2,zVecL,xl0,xln,t0,t0n,vld);

 // copy origK (original stiffness matrix) to elK
  for(i=0; i<12; ++i)
   for(j=0; j<12; ++j)
     elK[i][j] = origK[i][j];

 // Form unprojected internal forces and initialize stiffness matrix
 _FORTRAN(dgemv)('N',12,12,1.0,(double *)origK,12,vld,1,0.0,locF,1);

 // Compute gradients of nodal deformational pseudorotations
 // and correct stiffness and force

 double rotvar[2][3][3];

 for(inode = 0; inode < 2; ++inode)
   pseudorot_var(vld+inode*6+3, rotvar[inode]);

 double fe[12];
 for(inode = 0; inode < 2; ++inode)
   for(i = 0; i < 3; ++i) {
     fe[6*inode+i]   = locF[6*inode+i];
     fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                       rotvar[inode][1][i]*locF[6*inode+4] +
                       rotvar[inode][2][i]*locF[6*inode+5];
   }

 // Compute nonlinear projector matrix relative to deformed element
 // and correct stiffness and force
 // pmat = projector matrix
 double pmat[12][12];
 double gmat[3][12];
 gradDefDisp( zVecL, xln, pmat, gmat);

 // Form: {f} = [P']{fe}

  _FORTRAN(dgemv)('N',12,12,1.0,(double *)pmat,12,fe,1,0.0,f,1);

// Form geometric stiffness from internal force and material stiffness

 double stiffGeo1[12][12], stiffGeo2[12][12];
 //corotStiffGeo(zVecL,xln,pmat,f,stiffGeo1,stiffGeo2);
 formCorrectGeometricStiffness(pmat, gmat, f, stiffGeo1,
                               stiffGeo2, fe, node1, node2, ns1, ns2);

 for(i = 0; i < 12; ++i)
   for(j =0; j < 12; ++j)
     elK[i][j] = stiffGeo1[i][j] + stiffGeo2[i][j];

 // transform geometric stiffness to global coordinate system
 tran_stiff(elK, t0n);

 // The skew symmetric load stiffness matrix due to axial external moments is
 // added separately (see Domain::getFollowerForce in Driver.d/NLStatic.C)
 // For now, it can be removed from elK by symmetrizing to avoid doubling.
 elK.symmetrize();
}

//-----------------------------------------------------------------------

void
BeamCorotator::extractDefDisp(Node &nd1, Node &nd2, NodeState &ns1, 
                            NodeState &ns2, double zVecL[2][3], 
                            double xl0[2][3], double xln[2][3], 
                            double t0[3][3], double t0n[3][3], double vld[12])
/****************************************************************
 *
 *  Purpose:
 *     Form the deformational displacement vector for an element
 *     by subtracting the rigid body displacements.
 *     The deformational nodal rotations are supposed to have
 *     vector properties.
 * 
 *     routine corot_defdisp in C++ programming
 *
 *  Output:
 *     zvecl:  local coordinate z axis for the nodes ( beam only )
 *     xl0  :  xl0[i][j] is coordinate component j of node i.
 *             Local coordinate system coordinates of initial element.
 *     xln  :  xln[i][j] is coordinate component j of node i.
 *             Local coordinate system coordinates of deformed element.
 *     t0   :  t0[i][j] is vector component j of basevector i
 *             for the initial element. i.e.
 *             transformation matrix to global.
 *     t0n  :  t0n[i][j] is vector component j of basevector i
 *             for the corotated and deformed element. i.e.
 *             transformation matrix to global.
 *     vld  :  deformational displacement vector for the element in
 *             local coordinate system.
 *
 *  Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
 int i, j, inode;
 double x0[2][3], xn[2][3];
 double (* rot[2])[3][3];

 x0[0][0]  = nd1.x; //x coordinate of node 1
 x0[0][1]  = nd1.y; //y coordinate of node 1
 x0[0][2]  = nd1.z; //z coordinate of node 1

 x0[1][0]  = nd2.x; //x coordinate of node 2
 x0[1][1]  = nd2.y; //y coordinate of node 2
 x0[1][2]  = nd2.z; //z coordinate of node 2

 xn[0][0]  = ns1.x; //x coordinate of node state 1
 xn[0][1]  = ns1.y; //y coordinate of node state 1
 xn[0][2]  = ns1.z; //z coordinate of node state 1

 xn[1][0]  = ns2.x; //x coordinate of node state 2
 xn[1][1]  = ns2.y; //y coordinate of node state 2
 xn[1][2]  = ns2.z; //z coordinate of node state 2

 rot[0]    = &(ns1.R); // rotation tensor of node state 1
 rot[1]    = &(ns2.R); // rotation tensor of node state 2

 // convert beam coordinates to local coordinate system
 localCoord(zVec, zVecL, rot, x0, xn, t0, t0n, xl0, xln);

 // translation part of the deformation vector for node 1
 vld[0] = xln[0][0] - xl0[0][0];
 vld[1] = xln[0][1] - xl0[0][1];
 vld[2] = xln[0][2] - xl0[0][2];

 // translation part of the deformation vector for node 2
 vld[6] = xln[1][0] - xl0[1][0];
 vld[7] = xln[1][1] - xl0[1][1];
 vld[8] = xln[1][2] - xl0[1][2];

 // Create rotation part of the deformation vector
 double rot0[3][3], dr[3][3];

 for(inode = 0; inode < 2; ++inode) {
    for(i = 0; i < 3; ++i)
       for(j = 0; j < 3; ++j)
         rot0[i][j] = (*rot[inode])[i][0]*t0[j][0]
                     +(*rot[inode])[i][1]*t0[j][1]
                     +(*rot[inode])[i][2]*t0[j][2];

    for(i = 0; i < 3; ++i)
       for(j = 0; j < 3; ++j)
         dr[i][j] = t0n[i][0]*rot0[0][j]
                   +t0n[i][1]*rot0[1][j]
                   +t0n[i][2]*rot0[2][j];
    mat_to_vec(dr, vld + inode*6 + 3);
 }

}
	
//-----------------------------------------------------------------------

void
BeamCorotator::corotStiffGeo(double zVecL[2][3],
              double xln[2][3], double pmat[12][12], double f[12], 
              double stiffGeo1[12][12], double stiffGeo2[12][12])
/****************************************************************
 *
 *  Purpose:
 *     Form geometric stiffness in local coordinate system for
 *     a corotational element based on current material stiffness
 *     and internal force.
 *
 *  Input:
 *     zvecl  : local coordinate z axis for nodes ( beam only )
 *     xl0    : local coordinates of undeformed shadow element
 *              xl0[i][j] is local coordinate component j of node i
 *              element is assumed best fit in xy-plane. i.e.
 *              z coordinate is zero for 3-nod element.
 *     xln    : local coordinates of deformed element
 *     pmat   : nonlinear projector matrix for current configuration
 *     fint   : internal force
 *
 *  Output:
 *     stiffGeo1  : geometric stiffness, rotation contribution
 *     stiffGeo2  : geometric stiffness, equilibrium contribution
 *
 *  Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int i, j, k;
   double fspin[12][3], fproj[12][3], gmat[3][12];

// Fspin with both axial and moment contributions
   spinAxialAndMoment( f, fspin );

// Gradients of local rotations with respect to local element dofs
   double len = xln[1][0] - xln[0][0];
   gradLocRot( len, zVecL, gmat);

// Geometric stiffness contribution Kgeo1 = -Fspin*Gmat
   for( i=0; i<12; i++ ) {
      for( j=0; j<12; j++ ) {
         stiffGeo1[i][j] = -( fspin[i][0]*gmat[0][j]
                             +fspin[i][1]*gmat[1][j]
                             +fspin[i][2]*gmat[2][j] );
      }
   }

// Geometric stiffness contribution Kgeo2 = -Gmat'*Fspin'*Pmat
// where Fspin does not contain any moment contributions

// Fspin with only axial contributions
   spinAxial( f, fspin );

// Compute Fproj' = Fspin'*Pmat
   for( i=0; i<3; i++ ) {
      for( j=0; j<12; j++ ) {
         fproj[j][i] = 0.0;
         for( k=0; k<12; k++ ) {
            fproj[j][i] += fspin[k][i]*pmat[k][j];
         }
      }
   }

   for( i=0; i<12; i++ ) {
      for( j=0; j<12; j++ ) {
         stiffGeo2[i][j] = -( gmat[0][i]*fproj[j][0]
                             +gmat[1][i]*fproj[j][1]
                             +gmat[2][i]*fproj[j][2] );
      }
   }

}

//-------------------------------------------------------------------

void 
BeamCorotator::spinAxialAndMoment ( double f[], double fnm[][3])
/****************************************************************
 *
 *  Purpose:
 *     Form Fnm matrix based on internal force vector f
 *     routine form_fnm in c  programming
 *
 *  Input:
 *     f      : internal force vector
 *
 *  Output:
 *     fnm    : matrix of | Spin(n) | for each  node ordered as columns
 *                        | Spin(m) |
 *
 *  Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int dof, nod;

// Fspin with both axial and moment contributions
   for( nod=0; nod<2; nod++ ) {
   // Zero diagonal elements
      for( dof=0; dof<3; dof++ ) {
         fnm[nod*6    +dof][dof] = 0.0;
         fnm[nod*6 +3 +dof][dof] = 0.0;
      }
   // Nonzero axial force
      fnm[nod*6 +1][2] = -f[nod*6 +0];
      fnm[nod*6   ][2] =  f[nod*6 +1];
      fnm[nod*6   ][1] = -f[nod*6 +2];
      fnm[nod*6 +2][1] =  f[nod*6 +0];
      fnm[nod*6 +2][0] = -f[nod*6 +1];
      fnm[nod*6 +1][0] =  f[nod*6 +2];
   // Nonzero moment force
      fnm[nod*6 +4][2] = -f[nod*6 +3];
      fnm[nod*6 +3][2] =  f[nod*6 +4];
      fnm[nod*6 +3][1] = -f[nod*6 +5];
      fnm[nod*6 +5][1] =  f[nod*6 +3];
      fnm[nod*6 +5][0] = -f[nod*6 +4];
      fnm[nod*6 +4][0] =  f[nod*6 +5];
   }

}


void
BeamCorotator::gradLocRot(double len, double zVecL[2][3], 
                                                  double gmat[3][12])
/***********************************************************************
 *
 *   Compute rotation gradient matrix for a 3 node element
 *   routine gmat_3node in c programs
 *
 *   Input:
 *     fitalg : == 1 x axis along side 1-2, both for t0 and t0n
 *              == 2 x axis rotated from side 1-2 in so that
 *                   sum of angles between xl0 and xln side edges
 *                   are zero.
 *              == 3 x axis rotated so that xl0 is rotated relative
 *                   to xln with the continuum mechanics definition
 *                   of rotations
 *   x       :  delta x-coordinates  ( element assumed in xy plane )
 *   y       :  delta y-coordinates  ( delta coordinates of Cn )
 *
 *   Output:
 *   gmat :  rotation gradients of rot_x rot_y rot_z in local system
 *           with respect to the local coord. nodal degree of freedoms
 *           The nodal dof ordering is : tx ty tz rx ry rz for each node
 *
 *   Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 **********************************************************************/
{
   int    i, j;
   double zsum[3];

// Compute sum of nodal z axis vectors
   zsum[0] = zVecL[0][0] + zVecL[1][0];
   zsum[1] = zVecL[0][1] + zVecL[1][1];
   zsum[2] = zVecL[0][2] + zVecL[1][2];

// Initialize all entries in gmat
   for( i=0; i<3; i++ )
     for( j=0; j<12; j++ ) 
       gmat[i][j] = 0.0;

   double coef = 1.0/len;

// Theta_y = ( w1 - w2 )/l
   gmat[1][2] =  coef;
   gmat[1][8] = -coef;

// Theta_z = ( -v1 + v2 )/l
   gmat[2][1] = -coef;
   gmat[2][7] =  coef;


// Fitalg == 1: Rotation of Z axis from node 1
   if (fitAlg == 1) {
      for( i=0; i<3; i++ ) {
         gmat[0][i]   =   gmat[2][i]*zVecL[0][0]/zVecL[0][2];
         gmat[0][6+i] = gmat[2][6+i]*zVecL[0][0]/zVecL[0][2];
      }
      gmat[0][3]  = 1.0;                            /* theta_x node 1 */
      gmat[0][4]  = 0.0;                            /* theta_y node 1 */
      gmat[0][5]  = -zVecL[0][0]/(zVecL[0][2]);     /* theta_z node 1 */
      gmat[0][9]  = 0.0;                            /* theta_x node 2 */
      gmat[0][10] = 0.0;                            /* theta_y node 2 */
      gmat[0][11] = 0.0;                            /* theta_z node 2 */
   }
   else {
      for( i=0; i<3; i++ ) {
         gmat[0][i]   =   gmat[2][i]*zsum[0]/zsum[2];
         gmat[0][6+i] = gmat[2][6+i]*zsum[0]/zsum[2];
      }
      gmat[0][3]  =  zVecL[0][2]/zsum[2];              /* theta_x node 1 */
      gmat[0][4]  =  0.0;                              /* theta_y node 1 */
      gmat[0][5]  = -zVecL[0][0]/zsum[2];              /* theta_z node 1 */
      gmat[0][9]  =  zVecL[1][2]/zsum[2];              /* theta_x node 2 */
      gmat[0][10] =  0.0;                              /* theta_y node 2 */
      gmat[0][11] = -zVecL[1][0]/zsum[2];              /* theta_z node 2 */
   }

}



void
BeamCorotator::spinAxial( double f[], double fn[][3] )
/****************************************************************
 *
 *  Purpose:
 *     Form the Fn matrix based on the internal force vector force
 *
 *  Input:
 *     fint   : internal force vector
 *
 *  Output:
 *     fnm    : matrix of | Spin(n) | for each  node ordered as columns
 *                        |    0    |
 *
 *  Coded by: Bjorn Haugen;  Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int dof, nod;

// Fspin with only axial contributions

   for( nod=0; nod<2; nod++ ) {
   // Zero diagonal elements
      for( dof=0; dof<3; dof++ ) {
         fn[nod*6    +dof][dof] = 0.0;
         fn[nod*6 +3 +dof][dof] = 0.0;
      }
   // Nonzero axial force
      fn[nod*6 +1][2] = -f[nod*6 +0];
      fn[nod*6   ][2] =  f[nod*6 +1];
      fn[nod*6   ][1] = -f[nod*6 +2];
      fn[nod*6 +2][1] =  f[nod*6 +0];
      fn[nod*6 +2][0] = -f[nod*6 +1];
      fn[nod*6 +1][0] =  f[nod*6 +2];
   // Zero moment force
      fn[nod*6 +4][2] = 0.0;
      fn[nod*6 +3][2] = 0.0;
      fn[nod*6 +3][1] = 0.0;
      fn[nod*6 +5][1] = 0.0;
      fn[nod*6 +5][0] = 0.0;
      fn[nod*6 +4][0] = 0.0;
   }

} 


void
BeamCorotator::gradDefDisp( double zVecL[][3], double xln[][3],
                            double pmat[12][12], double gmat[3][12])
/***********************************************************************
 *
 *   Compute the gradients of the deformational displacement
 *   vector with respect to the visible dofs for an element
 *   Nonlinear version
 *
 *   Input:
 *     fitalg :  as defined in gmat_3nod and gmat_quad
 *     zVecL  :  local coordinate z axis for the nodes ( beam only )
 *     xl0    :  local coordinates of shadow element C0n
 *     xln    :  local coordinates of deformed element Cn
 *               xl0 and xln are assumed
 *               relative to element centroid and
 *               best fit in x-y coordinate system )
 *   Output:
 *     pmat   :  gradient matrix of deformational displacement vector
 *               with resect to the local coord. nodal degree of freedoms
 *               The nodal dof ordering is : tx ty tz rx ry rz for each node
 *
 *   Coded by Bjorn Haugen;  Adjusted for C++ by Teymour Manzouri
 **********************************************************************/
{
   int    nod, inod, jnod, dof;
   double fac,  len;

// Gmat( gradient of local rotations) dependent part of projector matrix

   len = xln[1][0] - xln[0][0];
   gradLocRot( len, zVecL, gmat);

   for( nod=0; nod<2; nod++ ) {
      for( dof=0; dof<12; dof++ ) {

      // Translation part
         pmat[nod*6  ][dof] = -xln[nod][2]*gmat[1][dof]
                              +xln[nod][1]*gmat[2][dof];
         pmat[nod*6+1][dof] =  xln[nod][2]*gmat[0][dof]
                              -xln[nod][0]*gmat[2][dof];
         pmat[nod*6+2][dof] = -xln[nod][1]*gmat[0][dof]
                              +xln[nod][0]*gmat[1][dof];
      // Rotation part
         pmat[nod*6+3][dof] = -gmat[0][dof];
         pmat[nod*6+4][dof] = -gmat[1][dof];
         pmat[nod*6+5][dof] = -gmat[2][dof];
      }
   }

   // Constant part of projector matrix
   // NOTE: fac = -1.0/number_of_nodes
   fac = -0.5;

   for( inod=0; inod<2; inod++ ) {

   // Diagonal translation part
      pmat[inod*6   ][inod*6   ] += 1.0;
      pmat[inod*6 +1][inod*6 +1] += 1.0;
      pmat[inod*6 +2][inod*6 +2] += 1.0;

   // Rotation part
      pmat[inod*6 +3][inod*6 +3] += 1.0;
      pmat[inod*6 +4][inod*6 +4] += 1.0;
      pmat[inod*6 +5][inod*6 +5] += 1.0;

      for( jnod=0; jnod<2; jnod++ ) {

      // Nondiagonal translation part
         pmat[inod*6   ][jnod*6   ] += fac;
         pmat[inod*6 +1][jnod*6 +1] += fac;
         pmat[inod*6 +2][jnod*6 +2] += fac;
      }
   }

}

#include <cmath>
#include <Corotational.d/utilities.h>

void
BeamCorotator::localCoord(double zvec0[3], double zvecl[2][3], 
               double (* rot[2])[3][3], double x0[2][3], double xn[2][3], 
      double t0[3][3],double t0n[3][3],double xl0[2][3], double xln[2][3] )
/*****************************************************************
 *
 *  Purpose:
 *     Form the transformation matrix to local coordinate system
 *     for a 2 noded element in deformed and undeformed configuration.
 *     Form deformed and undeformed coordinates for the element.
 *
 *  Input;
 *     zvec0    : local coordinate system z-axis in undeformed configuration
 *     rot      : nodal rotation tensors
 *     x0[i][j] : global coordinate component j of node i in undeformed 
 *                configuration
 *     xn[i][j] : global coordinate component j of node i in deformed 
 *                configuration
 *
 *  Output:
 *     zvecl    : nodal z-axis in deformed configuration in local 
 *                coordinate system
 *     t0       : transformation matrix to local undeformed coordinates
 *     t0n      : transformation matrix to local deformed coordinates
 *     xl0[i][j]: local coordinate component j of node i in undef. conf.
 *                also best fit of undeformed element on deformed element.
 *     xln[i][j]: local coordinate component j of node i in deformed conf.
 *
 *  Coded by: Bjorn Haugen; adjusted for C++ by Teymour Manzouri
 *******************************************************************/
{
   double ln, rvec[3]; 

// Local X-axis for C0 configuration
   t0[0][0] = x0[1][0] - x0[0][0];
   t0[0][1] = x0[1][1] - x0[0][1];
   t0[0][2] = x0[1][2] - x0[0][2];

   double l0 = sqrt( t0[0][0]*t0[0][0] + t0[0][1]*t0[0][1] + t0[0][2]*t0[0][2]);
   t0[0][0] /= l0;
   t0[0][1] /= l0;
   t0[0][2] /= l0;

// Local Y-axis for the C0 configuration as cross product of z-axis and x-axis
   crossprod( zvec0, t0[0], t0[1] );
   normalize( t0[1] );

// Local Z-axis for the C0 configuration as cross product of x-axis and y-axis
   crossprod( t0[0], t0[1], t0[2] );
   normalize( t0[2] );

// Compute nodal rotated Z-axis in global coordinate system

   int i,nod;
   for(nod=0; nod<2; nod++ ) {
      for(i=0; i<3; i++ ) {
         zvecl[nod][i] = (*rot[nod])[i][0]*t0[2][0]
                        +(*rot[nod])[i][1]*t0[2][1]
                        +(*rot[nod])[i][2]*t0[2][2];
      }
   }

/* Fitalg 1: Z-axis from node 1 */

   if (fitAlg == 1) {
      t0n[2][0] = zvecl[0][0];
      t0n[2][1] = zvecl[0][1];
      t0n[2][2] = zvecl[0][2];
   }

/* Fitalg .ne. 1: Z-axis as sum of nodal z-axis */
   else {
      t0n[2][0] = zvecl[0][0] + zvecl[1][0];
      t0n[2][1] = zvecl[0][1] + zvecl[1][1];
      t0n[2][2] = zvecl[0][2] + zvecl[1][2];
   }

/* X-axis along element in Cn */
   t0n[0][0]  = xn[1][0] - xn[0][0];
   t0n[0][1]  = xn[1][1] - xn[0][1];
   t0n[0][2]  = xn[1][2] - xn[0][2];

   ln = sqrt(  t0n[0][0]*t0n[0][0]+t0n[0][1]*t0n[0][1]+t0n[0][2]*t0n[0][2] );
   t0n[0][0] /= ln;
   t0n[0][1] /= ln;
   t0n[0][2] /= ln;

/* Y-axis as cross product between z and x */
   crossprod( t0n[2], t0n[0], t0n[1]);
   normalize( t0n[1] );

/* Z-axis as cross product between x and y */
   crossprod( t0n[0], t0n[1], t0n[2]);
   normalize( t0n[2] );

/* Compute nodal rotated Z-axis in local coordinate system */

   for( nod=0; nod<2; nod++ ) {
      for( i=0; i<3; i++ ) {
         rvec[i] = (*rot[nod])[i][0]*t0[2][0]
                  +(*rot[nod])[i][1]*t0[2][1]
                  +(*rot[nod])[i][2]*t0[2][2];
      }
      for( i=0; i<3; i++ ) {
         zvecl[nod][i] = t0n[i][0]*rvec[0]
                        +t0n[i][1]*rvec[1]
                        +t0n[i][2]*rvec[2];
      }
   }

  double halfLength0 = 0.5*l0;
  double halfLengthn = 0.5*ln;

  xl0[0][0] = -halfLength0;
  xl0[0][1] =  0.0;
  xl0[0][2] =  0.0;

  xl0[1][0] =  halfLength0;
  xl0[1][1] =  0.0;
  xl0[1][2] =  0.0;

  xln[0][0] = -halfLengthn;
  xln[0][1] =  0.0;
  xln[0][2] =  0.0;

  xln[1][0] =  halfLengthn;
  xln[1][1] =  0.0;
  xln[1][2] =  0.0;
}

void
BeamCorotator::extractDeformations(GeomState &geomState, CoordSet &cs, 
                                   double *vld, int &nlflag)
{
 // Set Flag to Use Linear Routines for Stress
 nlflag = 1;

 // Get Nodes original coordinates
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Get Nodes current coordinates
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 double xl0[2][3], xln[2][3], t0[3][3], t0n[3][3], zVecL[2][3];

 // Extract deformational displacement
 extractDefDisp(node1,node2,ns1,ns2,zVecL,xl0,xln,t0,t0n,vld);

 // transform element displacement vector from local to global coordinates
 tran_force( vld, t0, 2 );
}

//
// KHPXXX: EXTRACTION OF RIGID BODY MOTION AT ELEMENT LEVEL
//         NECESSARY FOR THE RIGID BODY TRANFORMATION OF THE
//         FLUID MESH TO AVOID ELEMENT DISTORTION

void
BeamCorotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                      double *vlr)
{
 double vld[12];
 int nlflag=1;
 extractDeformations(geomState, cs, vld, nlflag);

 // Get Nodes current coordinates
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];
 
 vlr[0]  = ns1.x - vld[0];
 vlr[1]  = ns1.y - vld[1];
 vlr[2]  = ns1.z - vld[2];
 vlr[3]  = 0.0;
 vlr[4]  = 0.0;
 vlr[5]  = 0.0;

 vlr[6]  = ns2.x - vld[6];
 vlr[7]  = ns2.y - vld[7];
 vlr[8]  = ns2.z - vld[8];
 vlr[9]  = 0.0;
 vlr[10] = 0.0;
 vlr[11] = 0.0;
}

double
BeamCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
// Computes Internal Energy of Element in Given State

 // Get Nodes original coordinates
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Get Nodes current coordinates
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 int i,j;
 double xl0[2][3], xln[2][3], t0[3][3], t0n[3][3], zVecL[2][3];
 double vld[12], tmp[12];

 // Extract deformational displacement
 extractDefDisp(node1,node2,ns1,ns2,zVecL,xl0,xln,t0,t0n,vld);

 // transform element displacement vector from local to global coordinates
 //why do they do this?
 //tran_force( vld, t0, 2 );
 
 // Multiply by the original stiffness matrix
 for(i=0; i<12; ++i) {
   tmp[i] = 0.0;
   for(j=0; j<12; ++j) {
     tmp[i] += origK[i][j]*vld[j];
   }
 }

 // Compute Energy 
 double Energy = 0.0;
 for(i=0; i<12; ++i)
   Energy += tmp[i]*vld[i];

 Energy *= 0.5;
 return Energy;
}

//------------------------------------------------------------------------------

void
BeamCorotator::formCorrectGeometricStiffness( 
                            double pmat[12][12],double gmat[3][12], 
			    double f[12], double stiffGeo1[12][12], 
                            double stiffGeo2[12][12],double fe[12],
			    Node &node1, Node &node2, NodeState &ns1,
			    NodeState &ns2)
{
  int i, j, k, m, inod;
  double fspin[12][3], fproj[12][3];
  //zero stiffness matricies
  for(i=0; i<12; i++)
    for(j=0; j<12; j++) {
      if(j <= 2) {
        fspin[i][j] = 0.0;
        fproj[i][j] = 0.0;
      }
      stiffGeo1[i][j]=0.0;
      stiffGeo2[i][j]=0.0;
    }

  //First Compute Kgeo1
    
  // Fspin with both axial and moment contributions
  spinAxialAndMoment(f, fspin);

  // Geometric stiffness contribution Kgeo1 = -Fspin*Gmat ,Fspin(P'H'f)
  for(i=0; i<12; ++i) {
    for(j=0; j<12; ++j) {
      stiffGeo1[i][j] = -( fspin[i][0]*gmat[0][j]
                          +fspin[i][1]*gmat[1][j]
                          +fspin[i][2]*gmat[2][j] );
    }
  }

  // Geometric stiffness contribution Kgeo2  
  // 
  // Fspin with only axial contributions
  spinAxial(fe, fspin);

  // Compute Fproj' = Fspin'*Pmat
  for(i=0; i<3; ++i) {
    for(j=0; j<12; ++j) {
      fproj[j][i] = 0.0;
      for( k=0; k<12; ++k )
        fproj[j][i] += fspin[k][i]*pmat[k][j];
    } 
  }

  for(i=0; i<12; ++i)
    for(j=0; j<12; ++j)
      stiffGeo2[i][j] = -( gmat[0][i]*fproj[j][0]
                          +gmat[1][i]*fproj[j][1]
                          +gmat[2][i]*fproj[j][2] );

  //get contribution from the variation of G
  double PvarPartG[12][12];

  double x0[2][3], dx0[2][3], xn[2][3], dxn[2][3];
  double (*  rot[2])[3][3];
  double (* drot[2])[3][3];

  double drot1[3][3], drot2[3][3]; // PJSA
  drot[0] = &drot1;                // PJSA
  drot[1] = &drot2;                // PJSA

  x0[0][0]  = node1.x; //x coordinate of node 1
  x0[0][1]  = node1.y; //y coordinate of node 1
  x0[0][2]  = node1.z; //z coordinate of node 1
  dx0[0][0] = 0.0;
  dx0[0][1] = 0.0;
  dx0[0][2] = 0.0;
  x0[1][0]  = node2.x; //x coordinate of node 2
  x0[1][1]  = node2.y; //y coordinate of node 2
  x0[1][2]  = node2.z; //z coordinate of node 2
  dx0[1][0] = 0.0;
  dx0[1][1] = 0.0;
  dx0[1][2] = 0.0;

  xn[0][0]  = ns1.x; //x coordinate of node state 1
  xn[0][1]  = ns1.y; //y coordinate of node state 1
  xn[0][2]  = ns1.z; //z coordinate of node state 1

  xn[1][0]  = ns2.x; //x coordinate of node state 2
  xn[1][1]  = ns2.y; //y coordinate of node state 2
  xn[1][2]  = ns2.z; //z coordinate of node state 2

  rot[0]    = &(ns1.R); // rotation tensor of node state 1
  rot[1]    = &(ns2.R); // rotation tensor of node state 2

  //loop over each variation
  double temp[12];
  Vector temp_rot(3, 0.0);
  double rot_zero[3];
  for (i=0;i<3;i++)
    rot_zero[i] =0.0;
  double dummy[3][3];
  double dzVec[3] = { 0, 0, 0 };
  double zVecL[2][3], dzVecL[2][3];
  double t0[3][3], dt0[3][3];
  double t0n[3][3], dt0n[3][3];
  double xln[2][3], dxln[2][3];
  double xl0[2][3], dxl0[2][3];
  double Gvar[3][12];
  double tempstiff[12][12];
  double len, dlen;
	
  for(i=0; i<12; i++) {
	
    //zero variational terms	      
    for(j=0; j<2; j++)
      for(k=0; k<3; k++)
        dxn[j][k] = 0.0;
    for(j=0; j<2; j++)
      for(k=0; k<3; k++)
        for(m=0; m<3; m++)
          (*drot[j])[k][m] = 0.0; 
	    
    //write if statements to fill up the dxn, drot
    if(i<3)
      dxn[0][i] = 1.0;
    else if(i<6) {
      temp_rot.zero();
      temp_rot[i-3] = 1.0;
      dvec_to_mat(rot_zero, temp_rot.data(), dummy, *drot[0]);
      mat_mult_mat(*drot[0], *rot[0], dummy, 0);
      for(j=0; j<3; ++j) for(k=0; k<3; ++k) (*drot[0])[j][k] = dummy[j][k];
    }
    else if(i<9)
      dxn[1][i-6] = 1.0;
    else {
      temp_rot.zero();
      temp_rot[i-9] = 1.0;
      dvec_to_mat(rot_zero, temp_rot.data(), dummy, *drot[1]);
      mat_mult_mat(*drot[1], *rot[1], dummy, 0);
      for(j=0; j<3; ++j) for(k=0; k<3; ++k) (*drot[1])[j][k] = dummy[j][k];
    }

    localCoord(zVec, dzVec, zVecL, dzVecL, rot, drot, x0, dx0, xn, dxn, t0, dt0,
               t0n, dt0n, xl0, dxl0, xln, dxln);

    len =  xln[1][0] -  xln[0][0];
    dlen = dxln[1][0] - dxln[0][0];
    dgradLocRot(len, dlen, zVecL, dzVecL, gmat, Gvar);

    //create partial derivative of P wrt to G
    int nod,dof;
    for(nod=0; nod<2; nod++) {
      for(dof=0; dof<12; dof++) {

        // Translation part
        PvarPartG[nod*6  ][dof] = -xln[nod][2]*Gvar[1][dof]
                                  +xln[nod][1]*Gvar[2][dof];
        PvarPartG[nod*6+1][dof] =  xln[nod][2]*Gvar[0][dof]
                                  -xln[nod][0]*Gvar[2][dof];
        PvarPartG[nod*6+2][dof] = -xln[nod][1]*Gvar[0][dof]
                                  +xln[nod][0]*Gvar[1][dof];
        // Rotation part
        PvarPartG[nod*6+3][dof] = -Gvar[0][dof];
        PvarPartG[nod*6+4][dof] = -Gvar[1][dof];
        PvarPartG[nod*6+5][dof] = -Gvar[2][dof];
      }
    }

    _FORTRAN(dgemv)('N',12,12,1.0,(double*)PvarPartG,12,fe,1,0.0,temp,1);
    for(j=0; j<12; j++) tempstiff[j][i] = temp[j]; 
  }

  double t0nt[3][3];
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      t0nt[i][j] = t0n[j][i];

  for(i=0; i<12; i++) {
    for(j=0; j<12; j++) temp[j] = tempstiff[i][j];
    tran_force(temp, t0nt, 2);
    for(j=0; j<12; j++) stiffGeo2[i][j] += temp[j];
  }

}

//-----------------------------------------------------------------------------

void
BeamCorotator::reBuildorigK(FullSquareMatrix &_origK)
{
 // Copy Element's stiffness matrix to origK (original stiffness matrix)
 int i, j;
 for(i=0; i<12; ++i)
   for(j=0; j<12; ++j)
     origK[i][j] = _origK[i][j];
}


//----------------------------------------------------------------------------

void
BeamCorotator::dgradLocRot(double len, double dlen, double zVecL[2][3],
                           double dzVecL[2][3], double gmat[3][12],
			   double dgmat[3][12])
{
   int    i, j;
   double zsum[3], dzsum[3];

// Compute sum of nodal z axis vectors
   zsum[0] =  zVecL[0][0] +  zVecL[1][0];
   zsum[1] =  zVecL[0][1] +  zVecL[1][1];
   zsum[2] =  zVecL[0][2] +  zVecL[1][2];
  dzsum[0] = dzVecL[0][0] + dzVecL[1][0];
  dzsum[1] = dzVecL[0][1] + dzVecL[1][1];
  dzsum[2] = dzVecL[0][2] + dzVecL[1][2];
  
// Initialize all entries in gmat
   for( i=0; i<3; i++ )
     for( j=0; j<12; j++ ){ 
       gmat[i][j] = 0.0;
       dgmat[i][j] = 0.0;
       }
       
   double  coef = 1.0/len;
   double dcoef = -dlen/(len*len);

// Theta_y = ( w1 - w2 )/l
   gmat[1][2] =  coef;
   gmat[1][8] = -coef;
   dgmat[1][2] =  dcoef;
   dgmat[1][8] = -dcoef;

// Theta_z = ( -v1 + v2 )/l
   gmat[2][1] = -coef;
   gmat[2][7] =  coef;
   dgmat[2][1] = -dcoef;
   dgmat[2][7] =  dcoef;

// Fitalg == 1: Rotation of Z axis from node 1
   if (fitAlg == 1) {
      for( i=0; i<3; i++ ) {
         gmat[0][i]   =   gmat[2][i]*zVecL[0][0]/zVecL[0][2];
         gmat[0][6+i] = gmat[2][6+i]*zVecL[0][0]/zVecL[0][2];
	 dgmat[0][i]   =   dgmat[2][i]*zVecL[0][0]/zVecL[0][2]
	                  +gmat[2][i]*dzVecL[0][0]/zVecL[0][2]
			  -gmat[2][i]*zVecL[0][0]*dzVecL[0][2]/
			   (zVecL[0][2]*zVecL[0][2]);
         dgmat[0][6+i] = dgmat[2][6+i]*zVecL[0][0]/zVecL[0][2]
	                  +gmat[2][6+i]*dzVecL[0][0]/zVecL[0][2]
			  -gmat[2][6+i]*zVecL[0][0]*dzVecL[0][2]/
			   (zVecL[0][2]*zVecL[0][2]);
      }
      gmat[0][3]  = 1.0;                            /* theta_x node 1 */
      gmat[0][4]  = 0.0;                            /* theta_y node 1 */
      gmat[0][5]  = -zVecL[0][0]/(zVecL[0][2]);     /* theta_z node 1 */
     dgmat[0][5]  = -dzVecL[0][0]/(zVecL[0][2])+zVecL[0][0]*dzVecL[0][2]/(zVecL[0][2]*zVecL[0][2]);
      gmat[0][9]  = 0.0;                            /* theta_x node 2 */
      gmat[0][10] = 0.0;                            /* theta_y node 2 */
      gmat[0][11] = 0.0;                            /* theta_z node 2 */
   }
   else {
      for( i=0; i<3; i++ ) {
         gmat[0][i]   =   gmat[2][i]*zsum[0]/zsum[2];
	dgmat[0][i]   =   dgmat[2][i]*zsum[0]/zsum[2] + gmat[2][i]*dzsum[0]/zsum[2]
	                  -gmat[2][i]*zsum[0]*dzsum[2]/(zsum[2]*zsum[2]);
         gmat[0][6+i] = gmat[2][6+i]*zsum[0]/zsum[2];
	dgmat[0][6+i]   =   dgmat[2][6+i]*zsum[0]/zsum[2] + gmat[2][6+i]*dzsum[0]/zsum[2]
	                  -gmat[2][6+i]*zsum[0]*dzsum[2]/(zsum[2]*zsum[2]); 
      }
      gmat[0][3]  =  zVecL[0][2]/zsum[2];              /* theta_x node 1 */
     dgmat[0][3]  =  dzVecL[0][2]/zsum[2]-zVecL[0][2]*dzsum[2]/(zsum[2]*zsum[2]); 
      gmat[0][4]  =  0.0;                              /* theta_y node 1 */
      gmat[0][5]  = -zVecL[0][0]/zsum[2];              /* theta_z node 1 */
     dgmat[0][5]  = -dzVecL[0][0]/zsum[2]+zVecL[0][0]*dzsum[2]/(zsum[2]*zsum[2]);
      gmat[0][9]  =  zVecL[1][2]/zsum[2];              /* theta_x node 2 */
     dgmat[0][9]  =  dzVecL[1][2]/zsum[2]-zVecL[1][2]*dzsum[2]/(zsum[2]*zsum[2]);
      gmat[0][10] =  0.0;                              /* theta_y node 2 */
      gmat[0][11] = -zVecL[1][0]/zsum[2];              /* theta_z node 2 */
     dgmat[0][11] = -dzVecL[1][0]/zsum[2]+zVecL[1][0]*dzsum[2]/(zsum[2]*zsum[2]);
   }

}

//------------------------------------------------------------------------------------------------------

void
BeamCorotator::localCoord(double zVec[3], double dzVec[3], double zVecL[2][3], 
                 	  double dzVecL[2][3], double (* rot[2])[3][3], double (* drot[2])[3][3], 
			  double x0[2][3],
		 	  double dx0[2][3], double xn[2][3], double dxn[2][3], 
			  double t0[3][3], double dt0[3][3],double t0n[3][3], 
			  double dt0n[3][3], double xl0[2][3], double dxl0[2][3], 
		       	  double xln[2][3], double dxln[2][3])

{
   double ln, rvec[3]; 
   double dln, drvec[3];

// Local X-axis for C0 configuration
   t0[0][0] =  x0[1][0] -  x0[0][0];
   t0[0][1] =  x0[1][1] -  x0[0][1];
   t0[0][2] =  x0[1][2] -  x0[0][2];
  dt0[0][0] = dx0[1][0] - dx0[0][0];
  dt0[0][1] = dx0[1][1] - dx0[0][1];
  dt0[0][2] = dx0[1][2] - dx0[0][2];

   double l0  = sqrt( t0[0][0]*t0[0][0] + t0[0][1]*t0[0][1] + t0[0][2]*t0[0][2]);
   double dl0 = ( dt0[0][0]*t0[0][0] + dt0[0][1]*t0[0][1] + dt0[0][2]*t0[0][2])/l0;
  dt0[0][0] = dt0[0][0]/l0 - t0[0][0]*dl0/(l0*l0);
  dt0[0][1] = dt0[0][1]/l0 - t0[0][1]*dl0/(l0*l0);
  dt0[0][2] = dt0[0][2]/l0 - t0[0][2]*dl0/(l0*l0);
   t0[0][0] /= l0;
   t0[0][1] /= l0;
   t0[0][2] /= l0;

// Local Y-axis for the C0 configuration as cross product of z-axis and x-axis
   dcrossprod( zVec, t0[0], t0[1], dzVec, dt0[0], dt0[1] );
   dnormalize( t0[1], dt0[1] );

// Local Z-axis for the C0 configuration as cross product of x-axis and y-axis
   dcrossprod( t0[0], t0[1], t0[2], dt0[0], dt0[1], dt0[2]);
   dnormalize( t0[2], dt0[2] );

// Compute nodal rotated Z-axis in global coordinate system

   int i,nod;
   for(nod=0; nod<2; nod++ ) {
      for(i=0; i<3; i++ ) {
         zVecL[nod][i] = (*rot[nod])[i][0]*t0[2][0]
                        +(*rot[nod])[i][1]*t0[2][1]
                        +(*rot[nod])[i][2]*t0[2][2];
         dzVecL[nod][i] = (*rot[nod])[i][0]*dt0[2][0] + (*drot[nod])[i][0]*t0[2][0]
                         +(*rot[nod])[i][1]*dt0[2][1] + (*drot[nod])[i][1]*t0[2][1]
                         +(*rot[nod])[i][2]*dt0[2][2] + (*drot[nod])[i][2]*t0[2][2];
      }
   }

/* Fitalg 1: Z-axis from node 1 */

   if (fitAlg == 1) {
      t0n[2][0] = zVecL[0][0];
      t0n[2][1] = zVecL[0][1];
      t0n[2][2] = zVecL[0][2];
     dt0n[2][0] = dzVecL[0][0];
     dt0n[2][1] = dzVecL[0][1];
     dt0n[2][2] = dzVecL[0][2];
   }

/* Fitalg .ne. 1: Z-axis as sum of nodal z-axis */
   else {
      t0n[2][0] =  zVecL[0][0] +  zVecL[1][0];
      t0n[2][1] =  zVecL[0][1] +  zVecL[1][1];
      t0n[2][2] =  zVecL[0][2] +  zVecL[1][2];
     dt0n[2][0] = dzVecL[0][0] + dzVecL[1][0];
     dt0n[2][1] = dzVecL[0][1] + dzVecL[1][1];
     dt0n[2][2] = dzVecL[0][2] + dzVecL[1][2];
   }

/* X-axis along element in Cn */
   t0n[0][0]  =  xn[1][0] -  xn[0][0];
   t0n[0][1]  =  xn[1][1] -  xn[0][1];
   t0n[0][2]  =  xn[1][2] -  xn[0][2];
  dt0n[0][0]  = dxn[1][0] - dxn[0][0];
  dt0n[0][1]  = dxn[1][1] - dxn[0][1];
  dt0n[0][2]  = dxn[1][2] - dxn[0][2];

   ln = sqrt(  t0n[0][0]*t0n[0][0]+t0n[0][1]*t0n[0][1]+t0n[0][2]*t0n[0][2] );
   dln = (dt0n[0][0]*t0n[0][0]+dt0n[0][1]*t0n[0][1]+dt0n[0][2]*t0n[0][2])/ln;
   dt0n[0][0] = dt0n[0][0]/ln-t0n[0][0]*dln/(ln*ln);
   dt0n[0][1] = dt0n[0][1]/ln-t0n[0][1]*dln/(ln*ln);
   dt0n[0][2] = dt0n[0][2]/ln-t0n[0][2]*dln/(ln*ln);
   t0n[0][0] /= ln;
   t0n[0][1] /= ln;
   t0n[0][2] /= ln;

/* Y-axis as cross product between z and x */
   dcrossprod( t0n[2], t0n[0], t0n[1], dt0n[2], dt0n[0], dt0n[1]);
   dnormalize( t0n[1], dt0n[1] );

/* Z-axis as cross product between x and y */
   dcrossprod( t0n[0], t0n[1], t0n[2], dt0n[0], dt0n[1], dt0n[2]);
   dnormalize( t0n[2], dt0n[2] );

/* Compute nodal rotated Z-axis in local coordinate system */

   for( nod=0; nod<2; nod++ ) {
      for( i=0; i<3; i++ ) {
         rvec[i] = (*rot[nod])[i][0]*t0[2][0]
                  +(*rot[nod])[i][1]*t0[2][1]
                  +(*rot[nod])[i][2]*t0[2][2];
	drvec[i] = (*rot[nod])[i][0]*dt0[2][0] + (*drot[nod])[i][0]*t0[2][0]
                  +(*rot[nod])[i][1]*dt0[2][1] + (*drot[nod])[i][1]*t0[2][1]
                  +(*rot[nod])[i][2]*dt0[2][2] + (*drot[nod])[i][2]*t0[2][2];	  
      }
      for( i=0; i<3; i++ ) {
         zVecL[nod][i] = t0n[i][0]*rvec[0]
                        +t0n[i][1]*rvec[1]
                        +t0n[i][2]*rvec[2];
	 dzVecL[nod][i] = dt0n[i][0]*rvec[0] + t0n[i][0]*drvec[0]
	                 +dt0n[i][1]*rvec[1] + t0n[i][1]*drvec[1]
                         +dt0n[i][2]*rvec[2] + t0n[i][2]*drvec[2];		
      }
   }

  double  halfLength0 = 0.5*l0;
  double dhalfLength0 = 0.5*dl0;
  double  halfLengthn = 0.5*ln;
  double dhalfLengthn = 0.5*dln;

  xl0[0][0] = -halfLength0;
 dxl0[0][0] = -dhalfLength0;
  xl0[0][1] =  0.0;
 dxl0[0][1] =  0.0;
  xl0[0][2] =  0.0;
 dxl0[0][2] =  0.0;
 
  xl0[1][0] =  halfLength0;
 dxl0[1][0] =  dhalfLength0;
  xl0[1][1] =  0.0;
 dxl0[1][1] =  0.0;
  xl0[1][2] =  0.0;
 dxl0[1][2] =  0.0;
 
  xln[0][0] = -halfLengthn;
 dxln[0][0] = -dhalfLengthn;
  xln[0][1] =  0.0;
 dxln[0][1] =  0.0;
  xln[0][2] =  0.0;
 dxln[0][2] =  0.0;
 
  xln[1][0] =  halfLengthn;
 dxln[1][0] =  dhalfLengthn;
  xln[1][1] =  0.0;
 dxln[1][1] =  0.0;
  xln[1][2] =  0.0;
 dxln[1][2] =  0.0;
}
