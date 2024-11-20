#include <cstdio>
#include <cmath>

#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Utils.d/linkfc.h>
#include <Element.d/Element.h>
#include <Corotational.d/SpringCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>

extern "C" {
 void _FORTRAN(mstf21)(double*,double&, double&,double&,double&,
                               double&, double&,double&,double&,
                               double&, double&,double&,double&);

 void _FORTRAN(mstf22)(double*,double&, double&,double&,double&,
                               double&, double&,double&,double&,
                               double&, double&,double&,double&);
}

//-----------------------------------------------------------
// SpringCorotator Constructor for Translational Springs
SpringCorotator::SpringCorotator(CoordSet &cs, int *nn,
	double _Kx, double _Ky, double _Kz, 
	FullSquareMatrix & _origK) : origK(_origK)  {


  n1 = nn[0];
  n2 = nn[1];

  Kx = _Kx;
  Ky = _Ky;
  Kz = _Kz;

  type = translational;
}

//-----------------------------------------------------------
// SpringCorotator Constructor for Torsional Springs
SpringCorotator::SpringCorotator(int *nn, CoordSet &cs,
	int _type, double _Kx, double _Ky, double _Kz,
	 FullSquareMatrix & _origK) : origK(_origK)  {

  n1 = nn[0];

  if (_type == 1)  {
    n2 = -1;  //signifies 1-node torsion spring
    type = rotational1;
  }
  else if (_type == 2)  {
    n2 = nn[1];
    type = rotational2;  // signifies 2-node torsion spring
  }

  // set rotational stiffness coefficients
  Kx = _Kx;
  Ky = _Ky;
  Kz = _Kz;
}

//-----------------------------------------------------------

void
SpringCorotator::getStiffAndForce(GeomState &geomState, CoordSet & cs, 
                                  FullSquareMatrix &elK, double *f, double dt, double t)
/*******************************************************************
 *
 * Purpose :
 *  Compute Internal force vector and Tangent Stiffness matrix 
 *  for Co-rotated spring element in current configuration.
 *
 * Input Variables:
 *  geomState : current configuration 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * dofs      : number of total free degrees of freedom
 * nDim      : number of free degrees of freedom for 1-node
 * R	     : Global Rotation Matrix
 *
 * Output :
 *  f        : internal force vector in current configuration
 *  elK      : element tangent stiffness matrix
 *
 *****************************************************************/
{
  int i,j,k;
  int dofs;

  // get global rotation tensor
  double R[3][3];
  geomState.getGlobalRot(R);

  // form global stiffness Matrix: Kg = R*klocal*R(transpose)
  // mult. klocal * R(transpose)

  int nDim = 3;
  int *dim = new int[nDim];
  k = 0;
  for (i = 0; i < 3; i++)
    dim[k++] = i;

  // mult. klocal*R(trans) and put result in K, expanded stiff mat
  double K[3][3];
  for (i = 0; i < 3; i++)  {
    K[0][i] = Kx*R[i][0];
    K[1][i] = Ky*R[i][1];
    K[2][i] = Kz*R[i][2];
  }

  // mult. R*K to get expanded global stiff mat
  double resMat[3][3];
  mat_mult_mat(R, K, resMat, 0);

  double force[3];  // expanded force vector

/*
  fprintf(stderr,"Global Rot Mat: \n");
  for (i = 0; i < 3; i++)
    fprintf(stderr,"%f, %f, %f\n", R[i][0], R[i][1], R[i][2]);
*/

  // Case of 1-node Torsional Spring Element
  if(type == rotational1) {
    
    dofs = nDim;
    // compute local node rotation tensor: r(loc) = R(trans)*r(glob)
    double rTen[3][3];
    mat_mult_mat(R, geomState[n1].R, rTen, 1); 
    
    //compute rotation vector
    double rVec[3];
    mat_to_vec(rTen, rVec);

    // mult. Klocal * rVec
    double res[3];
    res[0] = Kx * rVec[0];
    res[1] = Ky * rVec[1];
    res[2] = Kz * rVec[2];

    // mult. R * res to get force vector
    for (i = 0; i < 3; i++)
      force[i] = R[i][0] * res[0] + R[i][1] * res[1] + R[i][2] * res[2];
    for (i = 0; i < dofs; i++)
      f[i] = force[dim[i]];
    
    for (i = 0; i < dofs; i++)
      for (j = 0; j < dofs; j++) 
        elK[i][j] = resMat[dim[i]][dim[j]];
  }

  // stiffmat computation for 2-node translation and rotation 
  // spring are the same
  /* Local Stiffness Matrix is:
        kx   0   0   -kx   0   0
        0    ky  0    0   -ky  0
        0    0   kz   0    0  -kz
       -kx   0   0    kx   0   0
        0   -ky  0    0    ky  0
        0    0  -kz   0    0   kz
  */
  else {
    dofs = 2 * nDim;

    // put results in element stiff mat.
    for (i = 0; i < nDim; i++)  {
      for (j = 0; j < nDim; j++)  {
        elK[i][j] = resMat[dim[i]][dim[j]];
        elK[i+3][j+3] = resMat[dim[i]][dim[j]];

        // compute off diag blocks
        elK[i+3][j] = -resMat[dim[i]][dim[j]];
        elK[i][j+3] = -resMat[dim[j]][dim[i]];
      }
    }

    // force computations are different for 2-node
    // torsion and translation springs
    if (type == rotational2)  {
      // case of 2-node torsional spring
      //Compute Forces: F = R*Klocal*R(trans)*rotationVector

      // Transpose R
      // compute local node rotation tensor: r(loc) = R(trans)*r(glob)
      double rTen[3][3];
      double rVec[2][3];

      // mult R(transpose)*local rot. tensor of node
      mat_mult_mat(R, geomState[n1].R, rTen, 1);
   
      //compute rotation vector
      mat_to_vec(rTen, rVec[0]);

      // Repeat for 2nd node
      mat_mult_mat(R, geomState[n2].R, rTen, 1);
      mat_to_vec(rTen, rVec[1]);
    
      // For Forces --> mult. Klocal * (rVec0 - rVec1)
      double res[3];
      res[0] = Kx * (rVec[0][0] - rVec[1][0]);
      res[1] = Ky * (rVec[0][1] - rVec[1][1]);
      res[2] = Kz * (rVec[0][2] - rVec[1][2]);

      // mult. R * res to get force vector
      for (i = 0; i < 3; i++)
        force[i] = R[i][0] * res[0] + R[i][1] * res[1] + R[i][2] * res[2];
    }

    // else case of 2-node Translational spring element
    else {

      double dx = cs[n1]->x - cs[n2]->x;
      double dy = cs[n1]->y - cs[n2]->y;
      double dz = cs[n1]->z - cs[n2]->z;

      // compute deformation
      double defNode[3];
      defNode[0] = geomState[n1].x - geomState[n2].x;
      defNode[1] = geomState[n1].y - geomState[n2].y;
      defNode[2] = geomState[n1].z - geomState[n2].z;

      // force = Kg(u1-u2)
      defNode[0] -= R[0][0]*dx + R[0][1]*dy + R[0][2]*dz;
      defNode[1] -= R[1][0]*dx + R[1][1]*dy + R[1][2]*dz;
      defNode[2] -= R[2][0]*dx + R[2][1]*dy + R[2][2]*dz;

      // compute forces
      for (i = 0; i < 3; i++)
        force[i] = resMat[i][0] * defNode[0]
                 + resMat[i][1] * defNode[1] 
                 + resMat[i][2] * defNode[2];
    }

    // remove forces for non-free dofs 
    for (i = 0; i < nDim; i++)  {
      f[i] = force[dim[i]];
      f[i+3] = -f[i];
    }
  }
    

/*
  // print global stiffness matrix
  fprintf(stderr,"\nStiffMat:\n");
    for (i = 0; i < dofs; i++)  {
      for (j = 0; j < dofs; j++)
        fprintf(stderr,"%e, ", elK[i][j]);
      fprintf(stderr,"\n");
    }
  // print forces
  for (i = 0; i < dofs; i++)
    fprintf(stderr,"Forces: %e, ",f[i]);
  fprintf(stderr,"\n");
*/
  delete [] dim;

  return;
}

//----------------------------------------------------------------------

void
SpringCorotator::extractDeformations(GeomState & geomState, 
                                     CoordSet & cs, double *vld, int &nlflag)
{
 // Set Flag for No Stress Recovery
 nlflag = 0;

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 double xn[2][3], x0[2][3];

 x0[0][0] = node1.x;
 x0[0][1] = node1.y;
 x0[0][2] = node1.z;

 x0[1][0] = node2.x;
 x0[1][1] = node2.y;
 x0[1][2] = node2.z;

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of Cn configuration
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Extract Deformations
 vld[0] = xn[0][0] - x0[0][0];
 vld[1] = xn[0][1] - x0[0][1];
 vld[2] = xn[0][2] - x0[0][2];

 vld[3] = xn[1][0] - x0[1][0];
 vld[4] = xn[1][1] - x0[1][1];
 vld[5] = xn[1][2] - x0[1][2];

}

void
SpringCorotator::formTransformationTensor(double x0[2][3], double t[3][3])
{

/*******************************************************************
 *
 * Purpose :
 *  Form transformation tensor from global to local for 2 node
 *  bar element.
 *
 * Static variables:
 *  unit : 3x3 identity matrix
 *
 * Input :
 *  x0[i][j] : global coordinate component j of node i in
 *             reference configuration
 * 
 * Output :
 *  t     : 3x3 rotation tensor from global to local
 *
 ******************************************************************/
     int    min;
     static double unit[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
 
  // Get local x-axis along bar
     t[0][0] = x0[1][0] - x0[0][0];
     t[0][1] = x0[1][1] - x0[0][1];
     t[0][2] = x0[1][2] - x0[0][2];

  // Get smallest component of local x-axis
     min = (fabs(t[0][0])   < fabs(t[0][1])) ?  0   : 1;
     min = (fabs(t[0][min]) < fabs(t[0][2])) ? min  : 2;

  // Get local y-axis as cross product of x and global min axis
     t[1][0] =  t[0][1]*unit[min][2] - unit[min][1]*t[0][2];
     t[1][1] = -t[0][0]*unit[min][2] + unit[min][0]*t[0][2];
     t[1][2] =  t[0][0]*unit[min][1] - unit[min][0]*t[0][1];

  // Normalize x-axis and y-axis vectors to unit vectors
     normalize( t[1] );
     normalize( t[0] );

  // Get local z-axis as cross product of x and y-axis
     crossprod( t[0], t[1], t[2]);
}


void
SpringCorotator::formInternalForce(double t[3][3], double p, double *f)
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector for Co-rotated spring/bar
 *  element in current configuration.
 *
 * Input :
 *  t[i][j]  : transformation matrix
 *  p        : axial force in local coordinates
 *
 * Output :
 *  f        : internal force vector in current configuration
 *
 *****************************************************************/
{
     // Compute internal force in local system and store in f

     f[0] = p;
     f[1] = 0.0;
     f[2] = 0.0;

     // Transform to global coordinate by Fg = T'*Fl and store in f

     f[3] = t[0][0]*f[0] + t[1][0]*f[1] + t[2][0]*f[2];
     f[4] = t[0][1]*f[0] + t[1][1]*f[1] + t[2][1]*f[2];
     f[5] = t[0][2]*f[0] + t[1][2]*f[1] + t[2][2]*f[2];

     f[0] = -f[3];
     f[1] = -f[4];
     f[2] = -f[5];
}

void
SpringCorotator::extractRigidBodyMotion(GeomState &geomState,CoordSet &cs,
                                        double *vlr)
{
}

double
SpringCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
// Computes Internal Energy of Element in Given State

  int i;
  double Energy = 0;

  // get global rotation tensor
  double R[3][3];
  geomState.getGlobalRot(R);

  // Case of 1-node Torsional Spring Element
  if(type == rotational1) {

    // compute local node rotation tensor: r(loc) = R(trans)*r(glob)
    double rTen[3][3];
    mat_mult_mat(R, geomState[n1].R, rTen, 1);

    //compute rotation vector
    double rVec[3];
    mat_to_vec(rTen, rVec);

    // mult. Klocal * rVec
    double res[3];
    res[0] = Kx * rVec[0];
    res[1] = Ky * rVec[1];
    res[2] = Kz * rVec[2];

    Energy +=  res[0]*rVec[0] + res[1]*rVec[1] + res[2]*rVec[2];

  } else {

    if (type == rotational2)  {
      // case of 2-node torsional spring

      // Transpose R
      // compute local node rotation tensor: r(loc) = R(trans)*r(glob)
      double rTen[3][3];
      double rVec[2][3];

      // mult R(transpose)*local rot. tensor of node
      mat_mult_mat(R, geomState[n1].R, rTen, 1);

      //compute rotation vector
      mat_to_vec(rTen, rVec[0]);

      // Repeat for 2nd node
      mat_mult_mat(R, geomState[n2].R, rTen, 1);
      mat_to_vec(rTen, rVec[1]);

      // For Forces --> mult. Klocal * (rVec0 - rVec1)
      double res[3], res2[3];
      res[0] = (rVec[0][0] - rVec[1][0]);
      res[1] = (rVec[0][1] - rVec[1][1]);
      res[2] = (rVec[0][2] - rVec[1][2]);
      res2[0] = Kx * res[0];
      res2[1] = Ky * res[1];
      res2[2] = Kz * res[2];

      Energy +=  res[0]*res2[0] + res[1]*res2[1] + res[2]*res2[2];

    }

    // else case of 2-node Translational spring element
    else {

      // mult. klocal*R(trans) and put result in K, expanded stiff mat
      double K[3][3];
      for (i = 0; i < 3; i++)  {
        K[0][i] = Kx*R[i][0];
        K[1][i] = Ky*R[i][1];
        K[2][i] = Kz*R[i][2];
      }

      // mult. R*K to get expanded global stiff mat
      double resMat[3][3];
      mat_mult_mat(R, K, resMat, 0);

      double dx = cs[n1]->x - cs[n2]->x;
      double dy = cs[n1]->y - cs[n2]->y;
      double dz = cs[n1]->z - cs[n2]->z;

      // compute deformation
      double defNode[3];
      defNode[0] = geomState[n1].x - geomState[n2].x;
      defNode[1] = geomState[n1].y - geomState[n2].y;
      defNode[2] = geomState[n1].z - geomState[n2].z;

      // force = Kg(u1-u2)
      defNode[0] -= R[0][0]*dx + R[0][1]*dy + R[0][2]*dz;
      defNode[1] -= R[1][0]*dx + R[1][1]*dy + R[1][2]*dz;
      defNode[2] -= R[2][0]*dx + R[2][1]*dy + R[2][2]*dz;

      // compute forces
      double force[3];
      for (i = 0; i < 3; i++)
        force[i] = resMat[i][0] * defNode[0]
                 + resMat[i][1] * defNode[1]
                 + resMat[i][2] * defNode[2];

      Energy +=  force[0]*defNode[0]+force[1]*defNode[1]+force[2]*defNode[2];
    }

  }
  Energy *= 0.5;
  return Energy;
}
