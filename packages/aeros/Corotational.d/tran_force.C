#include <cstdio>
#include <cmath>

void tran_force( double* force, double tmat[3][3], int num_nodes, int num_dofs_per_node )
/*****************************************************************
 *
 *  Purpose:
 *     Transform the force vector to global system
 *
 *  Method:
 *            force = T'*force
 *
 *  Input:
 *     force : local coordinate force vector
 *     tmat  : transformation matrix T to global system
 *
 *  Output:
 *     force : now with respect to global coordinate system
 *
 *  Coded by: Bjorn Haugen
 *******************************************************************/
{

   int dof, numdofs;
   double tf0, tf1, tf2;

// Tranform groups of 3 degrees of freedom

   numdofs = num_nodes*num_dofs_per_node;

   for( dof=0; dof<(numdofs-2); dof +=3 ) {
      tf0 = tmat[0][0]*force[dof  ]
           +tmat[1][0]*force[dof+1]
           +tmat[2][0]*force[dof+2];
      tf1 = tmat[0][1]*force[dof  ]
           +tmat[1][1]*force[dof+1]
           +tmat[2][1]*force[dof+2];
      tf2 = tmat[0][2]*force[dof  ]
           +tmat[1][2]*force[dof+1]
           +tmat[2][2]*force[dof+2];
      force[dof  ] = tf0;
      force[dof+1] = tf1;
      force[dof+2] = tf2;
   }

}

