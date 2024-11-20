#include <Driver.d/EFrameData.h>
#include <Driver.d/Mpc.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/utilities.hpp>
#include <Utils.d/dofset.h>
#include <Element.d/Element.h>
#include <Driver.d/ControlLawInfo.h>
#include <Utils.d/SolverInfo.h>
#include <Math.d/Vector.h>
#ifdef USE_EIGEN3
#include <Element.d/Function.d/Rotation.d/LeftCompositionRule.h>
#include <Element.d/Function.d/Rotation.d/RightCompositionRule.h>
#endif

//#define COMPUTE_GLOBAL_ROTATION
//#define NEW_PULL_BACK

extern SolverInfo &solInfo;

GeomState::GeomState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs, Elemset *elems, double *ndTemps)
 : X0(&cs)
/****************************************************************
 *
 *  Purpose: determine geometric state of nodal coordinates
 *
 *  Input:
 *     DofSetArray : Constrained Degree of freedom set array 
 *     CoordSet    : Coordinate set 
 *
 *  Output:
 *     ns   :  node state updated  
 *
 *  Coded by: Michel Lesoinne and Teymour Manzouri
 ***************************************************************/
{
  numnodes = dsa.numNodes();
  ns.resize(numnodes);
  loc.resize(numnodes);
  flag.resize(numnodes);
  haveRot = false;

  for(int i = 0; i < numnodes; ++i) {

    loc[i].resize(6);
    // Store location of each degree of freedom
    loc[i][0] = cdsa.locate( i, DofSet::Xdisp );
    loc[i][1] = cdsa.locate( i, DofSet::Ydisp );
    loc[i][2] = cdsa.locate( i, DofSet::Zdisp );
    loc[i][3] = cdsa.locate( i, DofSet::Xrot  );
    loc[i][4] = cdsa.locate( i, DofSet::Yrot  );
    loc[i][5] = cdsa.locate( i, DofSet::Zrot  );
    if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5]) haveRot = true;
 
    // Get Node i from the Coordinate (Node) set
    Node *node_i = cs[i];

    if (node_i)  {

      // Set the ith node's coordinates
      ns[i].x = node_i->x;
      ns[i].y = node_i->y;
      ns[i].z = node_i->z;
 
      // Set ith node's rotation tensor equal to identity
      ns[i].R[0][0] = 1.0;
      ns[i].R[0][1] = 0.0;
      ns[i].R[0][2] = 0.0;
      ns[i].R[1][0] = 0.0;
      ns[i].R[1][1] = 1.0;
      ns[i].R[1][2] = 0.0;
      ns[i].R[2][0] = 0.0;
      ns[i].R[2][1] = 0.0;
      ns[i].R[2][2] = 1.0;

      if(dsa[i].list() != 0) {
        flag[i] = 1;
      }
      else 
        flag[i] = 0;

    }
    else  {
      // Set ith node's coordinates equal to zero
      ns[i].x = 0.0;
      ns[i].y = 0.0;
      ns[i].z = 0.0;

      // Set ith node's rotation tensor equal to identity
      ns[i].R[0][0] = 1.0;
      ns[i].R[0][1] = 0.0;
      ns[i].R[0][2] = 0.0;
      ns[i].R[1][0] = 0.0;
      ns[i].R[1][1] = 1.0;
      ns[i].R[1][2] = 0.0;
      ns[i].R[2][0] = 0.0;
      ns[i].R[2][1] = 0.0;
      ns[i].R[2][2] = 1.0;

      int dof;
      if((dof = cdsa.locate( i, DofSet::LagrangeE )) > -1) {
        loc[i][0] = dof;
        flag[i] = 0;
      }
      else if((dof = cdsa.locate( i, DofSet::LagrangeI )) > -1) {
        loc[i][0] = dof;
        flag[i] = -1;
      }
      else 
        flag[i] = 0;
    }
    
  }

  // Initialize Global Rotation Matrix to Identity
  double zeroRot[3] = {0.0, 0.0, 0.0};
  computeRotMat(zeroRot, gRot);
  computeCG(refCG);

  numelems = 0;
  if(elems) {
    int last = elems->last();
    for(int i = 0; i < last; ++i) {
      if((*elems)[i]->numStates()) numelems++;
    }
    es.resize(numelems);
    numelems = 0;
    for(int i = 0; i < last; ++i) {
      int numStates = (*elems)[i]->numStates();
      if(numStates > 0) {
        es[numelems].numInternalStates = numStates;
        es[numelems].internalStates = new double[numStates];
        for(int j = 0; j < numStates; ++j) es[numelems].internalStates[j] = 0;
        (*elems)[i]->setStateOffset(0);
        (*elems)[i]->initStates(es[numelems].internalStates);
        emap[(*elems)[i]->getGlNum()] = numelems;
        numelems++;
      }
    }
  }

  numnodesFixed = numnodes;
  setNodalTemperatures(ndTemps);
}

CoordSet emptyCoord;

GeomState::GeomState() : ns(0), numnodes(0), numnodesFixed(0), loc(0), X0(&emptyCoord), flag(0), numelems(0), haveRot(false)
{
}

GeomState::GeomState(const CoordSet &cs) : X0(&cs), numelems(0), haveRot(false)
{
  numnodes = numnodesFixed = cs.size(); // Number of nodes
  ns.resize(numnodes);

  for(int i = 0; i < numnodes; ++i) {
    // Get Node i from the Coordinate (Node) set
    Node *node_i = cs[i];

    if (node_i)  {
      // Set the ith node's coordinates
      ns[i].x = node_i->x;
      ns[i].y = node_i->y;
      ns[i].z = node_i->z;
    }
    else  {
      // Set the ith node's coordinates equal to zero
      ns[i].x = 0.0;
      ns[i].y = 0.0;
      ns[i].z = 0.0;
    }

    // Set ith node's rotation tensor equal to identity
    ns[i].R[0][0] = 1.0;
    ns[i].R[0][1] = 0.0;
    ns[i].R[0][2] = 0.0;
    ns[i].R[1][0] = 0.0;
    ns[i].R[1][1] = 1.0;
    ns[i].R[1][2] = 0.0;
    ns[i].R[2][0] = 0.0;
    ns[i].R[2][1] = 0.0;
    ns[i].R[2][2] = 1.0;
  }
}

GeomState::~GeomState()
{
}

void
GeomState::resize(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs, Elemset *elems)
{
  multiplier_nodes.clear();
  X0 = &cs;
  if(true/*dsa.numNodes() != numnodesFixed*/) { // XXX
    numnodes = dsa.numNodes();
    ns.resize(numnodes);
    loc.resize(numnodes);
    flag.resize(numnodes);

    for(int i = numnodesFixed; i < numnodes; ++i) {

      loc[i].resize(6);
      // Store location of each degree of freedom
      loc[i][0] = cdsa.locate( i, DofSet::Xdisp );
      loc[i][1] = cdsa.locate( i, DofSet::Ydisp );
      loc[i][2] = cdsa.locate( i, DofSet::Zdisp );
      loc[i][3] = cdsa.locate( i, DofSet::Xrot  );
      loc[i][4] = cdsa.locate( i, DofSet::Yrot  );
      loc[i][5] = cdsa.locate( i, DofSet::Zrot  );
      if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5]) haveRot = true;
 
      // Get Node i from the Coordinate (Node) set
      Node *node_i = cs[i];

      if (node_i)  {

        // Set the ith node's coordinates
        ns[i].x = node_i->x;
        ns[i].y = node_i->y;
        ns[i].z = node_i->z;
 
        // Set ith node's rotation tensor equal to identity
        ns[i].R[0][0] = 1.0;
        ns[i].R[0][1] = 0.0;
        ns[i].R[0][2] = 0.0;
        ns[i].R[1][0] = 0.0;
        ns[i].R[1][1] = 1.0;
        ns[i].R[1][2] = 0.0;
        ns[i].R[2][0] = 0.0;
        ns[i].R[2][1] = 0.0;
        ns[i].R[2][2] = 1.0;

        if(dsa[i].list() != 0)  {
          flag[i] = 1;
        }
        else 
          flag[i] = 0;

      }
      else  {
        // Set ith node's coordinates equal to zero
        ns[i].x = 0.0;
        ns[i].y = 0.0;
        ns[i].z = 0.0;

        // Set ith node's rotation tensor equal to identity
        ns[i].R[0][0] = 1.0;
        ns[i].R[0][1] = 0.0;
        ns[i].R[0][2] = 0.0;
        ns[i].R[1][0] = 0.0;
        ns[i].R[1][1] = 1.0;
        ns[i].R[1][2] = 0.0;
        ns[i].R[2][0] = 0.0;
        ns[i].R[2][1] = 0.0;
        ns[i].R[2][2] = 1.0;

        int dof;
        if((dof = cdsa.locate( i, DofSet::LagrangeE )) > -1) {
          loc[i][0] = dof;
          flag[i] = 0;
        }
        else if((dof = cdsa.locate( i, DofSet::LagrangeI )) > -1) {
          loc[i][0] = dof;
          flag[i] = -1;
        }
        else 
          flag[i] = 0;
      }
    
    }

    int last = elems->last();
    for(int i = 0; i < last; ++i) {
      if((*elems)[i] && (*elems)[i]->isMpcElement() && (*elems)[i]->numInternalNodes() == 1) {
        LMPCons *lmpc = dynamic_cast<LMPCons*>((*elems)[i]);
        if(lmpc && lmpc->getSource() == mpc::ContactSurfaces) {
          int numNodes = (*elems)[i]->numNodes();
          int *nodes = (*elems)[i]->nodes();
          multiplier_nodes[lmpc->id] = nodes[numNodes-1];
          delete [] nodes;
        }
      }
    }
  }

  // XXX don't need to resize element states because added elements are constraint
  //     elements which don't have any internal variables.
}

void
GeomState::clearMultiplierNodes()
{
  numnodes -= multiplier_nodes.size();
  ns.resize(numnodes);
  multiplier_nodes.clear();
}

void
GeomState::resizeLocAndFlag(DofSetArray &cdsa)
{
  // note ns has already been resized
  // now we need to resize and update loc and flag
  loc.resize(ns.size());
  flag.resize(ns.size());

  for(int i = numnodes-multiplier_nodes.size(); i < ns.size(); ++i) {

    loc[i].resize(6); for(int j=0; j<6; ++j) loc[i][j] = -1;
    int dof;
    if((dof = cdsa.locate( i, DofSet::LagrangeE )) > -1) {
      loc[i][0] = dof;
      flag[i] = 0;
    }
    else if((dof = cdsa.locate( i, DofSet::LagrangeI )) > -1) {
      loc[i][0] = dof;
      flag[i] = -1;
    }
    else {
      flag[i] = 0;
    }
  }
}

void
GeomState::print()
{
  // Prints nodal coordinates and associated rotation tensor
  for(int i = 0; i < numnodes; ++i) 
    printNode(i);
}

void
GeomState::printNode(int i)
{
  fprintf(stderr,"inode\tx\ty\tz\n");
  fprintf(stderr,"#%d\t%e\t%e\t%e\n",i,ns[i].x,ns[i].y,ns[i].z);
  fprintf(stderr,"Rotation Tensor\n");
  fprintf(stderr,"% e % e % e\n",ns[i].R[0][0],ns[i].R[0][1],ns[i].R[0][2]);
  fprintf(stderr,"% e % e % e\n",ns[i].R[1][0],ns[i].R[1][1],ns[i].R[1][2]);
  fprintf(stderr,"% e % e % e\n",ns[i].R[2][0],ns[i].R[2][1],ns[i].R[2][2]);
}

GeomState &
GeomState::operator=(const GeomState &g2)
{
  // note: unlike the copy constructor, the assignment operator does not copy X0 or refCG
  if(numnodes != g2.numNodes()) {
    numnodes = g2.numNodes();
    ns.resize(g2.numNodes());
    loc.resize(g2.numNodes());
    flag.resize(g2.numNodes());
  }

  haveRot = g2.haveRot;

  // Copy dof locations
  for(int i = 0; i < numnodes; ++i) {
    loc[i].resize(6);
    loc[i][0] = g2.loc[i][0];
    loc[i][1] = g2.loc[i][1];
    loc[i][2] = g2.loc[i][2];
    loc[i][3] = g2.loc[i][3];
    loc[i][4] = g2.loc[i][4];
    loc[i][5] = g2.loc[i][5];
  }

  // Copy node states
  for(int i = 0; i < numnodes; ++i) {
    ns[i]  = g2.ns[i];
    flag[i]= g2.flag[i];
  }
  multiplier_nodes = g2.multiplier_nodes;

  // now deal with element states
  numelems = g2.numelems;
  es.clear();
  es.resize(numelems);
  for(int i = 0; i < numelems; ++i)
    es[i] = g2.es[i];
  emap = g2.emap;

  refCG[0] = g2.refCG[0];
  refCG[1] = g2.refCG[1];
  refCG[2] = g2.refCG[2];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      gRot[i][j] = g2.gRot[i][j];

  numnodesFixed = g2.numnodesFixed;

  return *this;
}

void
GeomState::extract(double *p)
{
  int i;
  for(i=0; i<numnodes; ++i) {

     // Set incremental translational displacements
     if(loc[i][0] >= 0) p[loc[i][0]] = ns[i].x;
     if(loc[i][1] >= 0) p[loc[i][1]] = ns[i].y;
     if(loc[i][2] >= 0) p[loc[i][2]] = ns[i].z;

     double rot[3];
     mat_to_vec(ns[i].R,rot);

     if(loc[i][3] >= 0) p[loc[i][3]] = rot[0];
     if(loc[i][4] >= 0) p[loc[i][4]] = rot[1];
     if(loc[i][5] >= 0) p[loc[i][5]] = rot[2];
  }
}

GeomState::GeomState(const GeomState &g2) : X0(g2.X0), emap(g2.emap), multiplier_nodes(g2.multiplier_nodes)
{
  // Copy number of nodes
  numnodes = g2.numnodes;

  // Allocate memory for node states & dof locations
  ns.resize(numnodes);
  loc.resize(numnodes);

  // flag for node to element connectivity
  flag.resize(numnodes);    

  haveRot = g2.haveRot;

  // Copy dof locations
  for(int i = 0; i < numnodes; ++i) {
    loc[i].resize(6); 
    loc[i][0] = g2.loc[i][0];
    loc[i][1] = g2.loc[i][1];
    loc[i][2] = g2.loc[i][2];
    loc[i][3] = g2.loc[i][3];
    loc[i][4] = g2.loc[i][4];
    loc[i][5] = g2.loc[i][5];
  }

  // Copy node states
  for(int i = 0; i < numnodes; ++i) {
    ns[i]  = g2.ns[i];
    flag[i]= g2.flag[i];
  }

  // now deal with element states
  numelems = g2.numelems;
  es.clear();
  es.resize(numelems);
  for(int i = 0; i < numelems; ++i)
    es[i] = g2.es[i];
 
  // Initialize Global Rotation Matrix & CG position
  refCG[0] = g2.refCG[0];
  refCG[1] = g2.refCG[1];
  refCG[2] = g2.refCG[2];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      gRot[i][j] = g2.gRot[i][j];

  numnodesFixed = g2.numnodesFixed;
}

void
NodeState::operator=(const NodeState &node)
{
  // Assign x, y, and z coordinate values
  x = node.x;
  y = node.y;
  z = node.z;

  // Assign rotation tensor
  R[0][0] = node.R[0][0];
  R[0][1] = node.R[0][1];
  R[0][2] = node.R[0][2];
  R[1][0] = node.R[1][0];
  R[1][1] = node.R[1][1];
  R[1][2] = node.R[1][2];
  R[2][0] = node.R[2][0];
  R[2][1] = node.R[2][1];
  R[2][2] = node.R[2][2];

  // Assign rotation vector
  theta[0] = node.theta[0];
  theta[1] = node.theta[1];
  theta[2] = node.theta[2];

  // Assign velocity and acceleration vectors
  for(int i = 0; i < 6; ++i) {
    v[i] = node.v[i];
    a[i] = node.a[i];
  }

  // Assign temperature
  temp = node.temp;
}

void
ElemState::operator=(const ElemState &elem)
{
  if(numInternalStates != elem.numInternalStates) {
    if(internalStates) delete [] internalStates;
    internalStates = 0;
    numInternalStates = elem.numInternalStates;
  }
  if(internalStates == 0) internalStates = new double[numInternalStates];
  for(int i = 0; i < numInternalStates; ++i) {
    internalStates[i] = elem.internalStates[i];
  }
}

void
GeomState::update(const Vector &v, int SO3param)
{
 // v = incremental displacement vector
 double d[3], dtheta[3];

 int i;
 for(i=0; i<numnodes; ++i) {

     if(flag[i] == -1) {
       double mu = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;
       ns[i].x = mu; // for inequality constraints, we solve for the lagrange multiplier not the increment
       continue;
     }

     // Set incremental translational displacements

     d[0] = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;
     d[1] = (loc[i][1] >= 0) ? v[loc[i][1]] : 0.0;
     d[2] = (loc[i][2] >= 0) ? v[loc[i][2]] : 0.0;

     // Transform from DOF_FRM to basic frame

     NFrameData *cd = X0->dofFrame(i);
     if(cd) cd->invTransformVector3(d);

     // Increment total translational displacements

     ns[i].x += d[0];
     ns[i].y += d[1];
     ns[i].z += d[2];

     if(getNumRotationDof(i) == 2 && SO3param < 2 && !solInfo.getNLInfo().linearelastic) {
#ifdef USE_EIGEN3
       switch(SO3param) {
         case 0: {
           Eigen::Array3d Psi; Psi << ns[i].theta[0], ns[i].theta[1], ns[i].theta[2];
           if(cd) cd->transformVector3(Psi.data());
           Simo::LeftCompositionRule<double> f(Psi,Eigen::Array<int,0,1>::Zero());
           Simo::Jacobian<double,Simo::LeftCompositionRule> dfdq(Psi,Eigen::Array<int,0,1>::Zero());
 
           // solve for constrained (j^th) component of the incremental (left) rotation vector such that
           // the constraint is not violated by the update
           int j;
           Eigen::Vector3d q;
           if(loc[i][3] < 0) { j = 0; q[0] = 0; } else q[0] = v[loc[i][3]];
           if(loc[i][4] < 0) { j = 1; q[1] = 0; } else q[1] = v[loc[i][4]];
           if(loc[i][5] < 0) { j = 2; q[2] = 0; } else q[2] = v[loc[i][5]];

           for(int k=0; k<10; ++k) {
             //std::cerr << "k = " << k << ", residual = " << f(q,0.)[j]-ns[i].theta[j] << std::endl;
             q[j] -= (f(q,0.)[j]-Psi[j])/dfdq(q,0.)(j,j);
           }

           dtheta[0] = q[0];
           dtheta[1] = q[1];
           dtheta[2] = q[2];
           if(cd) cd->invTransformVector3(dtheta);
           inc_rottensor( dtheta, ns[i].R );
           inc_rotvector( ns[i].theta, ns[i].R );
         } break;
         case 1: {
           Eigen::Array3d Psi; Psi << ns[i].theta[0], ns[i].theta[1], ns[i].theta[2];
           if(cd) cd->transformVector3(Psi.data());
           Simo::RightCompositionRule<double> f(Psi,Eigen::Array<int,0,1>::Zero());
           Simo::Jacobian<double,Simo::RightCompositionRule> dfdq(Psi,Eigen::Array<int,0,1>::Zero());

           // solve for constrained (j^th) component of the incremental (right) rotation vector such that
           // the constraint is not violated by the update
           int j;
           Eigen::Vector3d q;
           if(loc[i][3] < 0) { j = 0; q[0] = 0; } else q[0] = v[loc[i][3]];
           if(loc[i][4] < 0) { j = 1; q[1] = 0; } else q[1] = v[loc[i][4]];
           if(loc[i][5] < 0) { j = 2; q[2] = 0; } else q[2] = v[loc[i][5]];

           for(int k=0; k<10; ++k) {
             //std::cerr << "k = " << k << ", residual = " << f(q,0.)[j]-ns[i].theta[j] << std::endl;
             q[j] -= (f(q,0.)[j]-Psi[j])/dfdq(q,0.)(j,j);
           }

           dtheta[0] = q[0];
           dtheta[1] = q[1];
           dtheta[2] = q[2];
           if(cd) cd->invTransformVector3(dtheta);
           inc_rottensor( ns[i].R, dtheta );
           inc_rotvector( ns[i].theta, ns[i].R );
         } break;
       } 
#endif
     }
     else if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0) {

       // Set incremental rotations

       dtheta[0] = (loc[i][3] >= 0) ? v[loc[i][3]] : 0.0;
       dtheta[1] = (loc[i][4] >= 0) ? v[loc[i][4]] : 0.0;
       dtheta[2] = (loc[i][5] >= 0) ? v[loc[i][5]] : 0.0;

       // Transform from DOF_FRM to basic frame

       if(cd) cd->invTransformVector3(dtheta);

       if(solInfo.getNLInfo().linearelastic || SO3param == 2) {
         // Additive update of total rotation vector
         for(int j=0; j<3; ++j) ns[i].theta[j] += dtheta[j];
         vec_to_mat( ns[i].theta, ns[i].R );
       }
       else {
         switch(SO3param) {
           case 0: // Increment rotation tensor from the left R = R(dtheta)*Ra
             inc_rottensor( dtheta, ns[i].R );
             inc_rotvector( ns[i].theta, ns[i].R );
             break;
           case 1: // Increment rotation tensor from the right R = Ra*R(dtheta)
             inc_rottensor( ns[i].R, dtheta );
             inc_rotvector( ns[i].theta, ns[i].R );
             break;
         }
       }
     }
   }
#ifdef COMPUTE_GLOBAL_ROTATION
  computeGlobalRotation();
#endif
}

void
GeomState::update(const Vector &v, const std::vector<int> &weightedNodes, int SO3param)
{
 // v = incremental displacement vector

 double d[3], dtheta[3];

 int i;
  for(std::vector<int>::const_iterator it = weightedNodes.begin(); it != weightedNodes.end(); ++it) {
     i = *it;

     if(flag[i] == -1) {
       double mu = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;
       ns[i].x = mu; // for inequality constraints, we solve for the lagrange multiplier not the increment
       continue;
     }

     // Set incremental translational displacements

     d[0] = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;
     d[1] = (loc[i][1] >= 0) ? v[loc[i][1]] : 0.0;
     d[2] = (loc[i][2] >= 0) ? v[loc[i][2]] : 0.0;

     // Transform from DOF_FRM to basic frame

     NFrameData *cd = X0->dofFrame(i);
     if(cd) cd->invTransformVector3(d);

     // Increment total translational displacements

     ns[i].x += d[0];
     ns[i].y += d[1];
     ns[i].z += d[2];

     if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0) {

       // Set incremental rotations

       dtheta[0] = (loc[i][3] >= 0) ? v[loc[i][3]] : 0.0;
       dtheta[1] = (loc[i][4] >= 0) ? v[loc[i][4]] : 0.0;
       dtheta[2] = (loc[i][5] >= 0) ? v[loc[i][5]] : 0.0;

       // Transform from DOF_FRM to basic frame

       if(cd) cd->invTransformVector3(dtheta);

       if(solInfo.getNLInfo().linearelastic || SO3param == 2) {
         // Additive update of total rotation vector
         for(int j=0; j<3; ++j) ns[i].theta[j] += dtheta[j];
         vec_to_mat( ns[i].theta, ns[i].R );
       }
       else {
         switch(SO3param) {
           case 0: // Increment rotation tensor from the left R = R(dtheta)*Ra
             inc_rottensor( dtheta, ns[i].R );
             inc_rotvector( ns[i].theta, ns[i].R );
             break;
           case 1: // Increment rotation tensor from the right R = Ra*R(dtheta)
             inc_rottensor( ns[i].R, dtheta );
             inc_rotvector( ns[i].theta, ns[i].R );
             break;
         }
       }
     }
   }
#ifdef COMPUTE_GLOBAL_ROTATION
  computeGlobalRotation();
#endif
}

void
GeomState::explicitUpdate(CoordSet &cs, const Vector &v)
{
 // v = total displacement vector (unconstrained dofs only)

 for(int i = 0; i < numnodes; ++i) {
   if(cs[i]) {

     // Set position vector components for non-prescribed dofs
     NFrameData *cd = X0->dofFrame(i);
     if(cd) {
       double d[3] = { ns[i].x - cs[i]->x,
                       ns[i].y - cs[i]->y,
                       ns[i].z - cs[i]->z };
       cd->transformVector3(d);
       if(loc[i][0] >= 0) d[0] = v[loc[i][0]];
       if(loc[i][1] >= 0) d[1] = v[loc[i][1]];
       if(loc[i][2] >= 0) d[2] = v[loc[i][2]];
       cd->invTransformVector3(d);
       ns[i].x = cs[i]->x + d[0];
       ns[i].y = cs[i]->y + d[1];
       ns[i].z = cs[i]->z + d[2];
     }
     else {
       if(loc[i][0] >= 0) ns[i].x = cs[i]->x + v[loc[i][0]];
       if(loc[i][1] >= 0) ns[i].y = cs[i]->y + v[loc[i][1]];
       if(loc[i][2] >= 0) ns[i].z = cs[i]->z + v[loc[i][2]];
     }

     if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0) {

       // Set rotation vector components for non-prescribed dofs
       if(cd) cd->transformVector3(ns[i].theta);
       if(loc[i][3] >= 0) ns[i].theta[0] = v[loc[i][3]];
       if(loc[i][4] >= 0) ns[i].theta[1] = v[loc[i][4]];
       if(loc[i][5] >= 0) ns[i].theta[2] = v[loc[i][5]];
       if(cd) cd->invTransformVector3(ns[i].theta);

       // Set rotation tensor 
       form_rottensor( ns[i].theta, ns[i].R );
     }
   }
 }
#ifdef COMPUTE_GLOBAL_ROTATION
 computeGlobalRotation();
#endif
}

void
GeomState::explicitUpdate(CoordSet &cs, int numNodes, int *nodes, const Vector &v)
{
 // v = total displacement vector

 for(int j = 0; j < numNodes; ++j) {
   int i = nodes[j];
   if(cs[i]) {

     NFrameData *cd = X0->dofFrame(i);
     if(cd) {
       double d[3] = { ns[i].x - cs[i]->x,
                       ns[i].y - cs[i]->y,
                       ns[i].z - cs[i]->z };
       cd->transformVector3(d);
       if(loc[i][0] >= 0) d[0] = v[loc[i][0]];
       if(loc[i][1] >= 0) d[1] = v[loc[i][1]];
       if(loc[i][2] >= 0) d[2] = v[loc[i][2]];
       cd->invTransformVector3(d);
       ns[i].x = cs[i]->x + d[0];
       ns[i].y = cs[i]->y + d[1];
       ns[i].z = cs[i]->z + d[2];
     }
     else {
       // Set position vector components for non-prescribed dofs
       if(loc[i][0] >= 0) ns[i].x = cs[i]->x + v[loc[i][0]];
       if(loc[i][1] >= 0) ns[i].y = cs[i]->y + v[loc[i][1]];
       if(loc[i][2] >= 0) ns[i].z = cs[i]->z + v[loc[i][2]];
     }

     if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0) {

       // Set rotation vector components for non-prescribed dofs
       if(cd) cd->transformVector3(ns[i].theta);
       if(loc[i][3] >= 0) ns[i].theta[0] = v[loc[i][3]];
       if(loc[i][4] >= 0) ns[i].theta[1] = v[loc[i][4]];
       if(loc[i][5] >= 0) ns[i].theta[2] = v[loc[i][5]];
       if(cd) cd->invTransformVector3(ns[i].theta);
 
       // Set rotation tensor 
       form_rottensor( ns[i].theta, ns[i].R );
     }
   }
 }
#ifdef COMPUTE_GLOBAL_ROTATION
 computeGlobalRotation();
#endif
}

void
GeomState::update(GeomState &refState, const Vector &v, int SO3param)
{
 // v = incremental displacement vector w.r.t reference state
 double d[3], dtheta[3];

 int i;
 for(i=0; i<numnodes; ++i) {

     // Set incremental translational displacements

     d[0] = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;
     d[1] = (loc[i][1] >= 0) ? v[loc[i][1]] : 0.0;
     d[2] = (loc[i][2] >= 0) ? v[loc[i][2]] : 0.0;

     // Transform from DOF_FRM to basic frame

     NFrameData *cd = X0->dofFrame(i);
     if(cd) cd->invTransformVector3(d);

     // Increment total translational displacements

     ns[i].x = refState[i].x + d[0];
     ns[i].y = refState[i].y + d[1];
     ns[i].z = refState[i].z + d[2];

     if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0) {

       // Set incremental rotations

       dtheta[0] = (loc[i][3] >= 0) ? v[loc[i][3]] : 0.0;
       dtheta[1] = (loc[i][4] >= 0) ? v[loc[i][4]] : 0.0;
       dtheta[2] = (loc[i][5] >= 0) ? v[loc[i][5]] : 0.0;

       // Transform from DOF_FRM to basic frame

       if(cd) cd->invTransformVector3(dtheta);

       if(solInfo.getNLInfo().linearelastic || SO3param == 2) {
         // Additive update of total rotation vector
         for(int j=0; j<3; ++j) ns[i].theta[j] = refState[i].theta[j] + dtheta[j];
         vec_to_mat( ns[i].theta, ns[i].R );
       }
       else {
         std::cerr << "Error: SO3param = " << SO3param << " case not implemented in GeomState::update\n";
         exit(-1);
       }
     }
   }
}

void
GeomState::setVelocity(const Vector &v, int SO3param)
{
  // set unconstrained components of velocity for all nodes, leaving constrained components unchanged

  for(int i = 0; i < numnodes; ++i) {
    NFrameData *cd = X0->dofFrame(i);
    if(cd) cd->transformVector6(ns[i].v);
    for(int j = 0; j < 6; ++j)
      if(loc[i][j] > -1) {
        ns[i].v[j] = v[loc[i][j]];
      }
    if(cd) cd->invTransformVector6(ns[i].v);
#ifdef USE_EIGEN3
    if(SO3param == 2 && (loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0)
       && !solInfo.getNLInfo().linearelastic) {
      // conversion to convected angular velocity
      Eigen::Vector3d PsiI, PsiIdot, Omega;
      Eigen::Matrix3d R, T;
      PsiI << ns[i].theta[0], ns[i].theta[1], ns[i].theta[2];
      PsiIdot << ns[i].v[3], ns[i].v[4], ns[i].v[5];
      tangential_transf(PsiI, T);
      Omega = T*PsiIdot;
      for(int j=0; j<3; ++j) ns[i].v[3+j] = Omega[j];
    }
#endif
  }
}

void
GeomState::setVelocity(int numNodes, int *nodes, const Vector &v, int SO3param)
{
  // set unconstrained components of velocity for specified nodes, leaving constrained components unchanged

  int i;
  for(int k = 0; k < numNodes; ++k) {
    i = nodes[k];
    NFrameData *cd = X0->dofFrame(i);
    if(cd) cd->transformVector6(ns[i].v);
    for(int j = 0; j < 6; ++j)
      if(loc[i][j] > -1) {
        ns[i].v[j] = v[loc[i][j]];
      }
    if(cd) cd->invTransformVector6(ns[i].v);
#ifdef USE_EIGEN3
    if(SO3param == 2 && (loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0)) {
      // transform from time derivative to total rotation vector to convected angular velocity
      Eigen::Vector3d PsiI, PsiIdot, Omega;
      Eigen::Matrix3d R, T;
      PsiI << ns[i].theta[0], ns[i].theta[1], ns[i].theta[2];
      PsiIdot << ns[i].v[3], ns[i].v[4], ns[i].v[5];
      tangential_transf(PsiI, T);
      Omega = T*PsiIdot;
      for(int j=0; j<3; ++j) ns[i].v[3+j] = Omega[j];
    }
#endif
  }
}

void
GeomState::setAcceleration(const Vector &a, int SO3param)
{
  // set unconstrained components of acceleration for all nodes, leaving constrained components unchanged

  for(int i = 0; i < numnodes; ++i) {
    NFrameData *cd = X0->dofFrame(i);
    if(cd) cd->transformVector6(ns[i].a);
    for(int j = 0; j < 6; ++j)
      if(loc[i][j] > -1) {
        ns[i].a[j] = a[loc[i][j]];
      }
    if(cd) cd->invTransformVector6(ns[i].a);
#ifdef USE_EIGEN3
    if(SO3param == 2 && (loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0)
       && !solInfo.getNLInfo().linearelastic) {
      // conversion from second time derivative of total rotation vector to convected angular acceleration
      Eigen::Vector3d PsiI, PsiIdot, PsiIddot, Omega, Alpha;
      Eigen::Matrix3d R, T, Tdot;
      PsiI << ns[i].theta[0], ns[i].theta[1], ns[i].theta[2];
      Omega << ns[i].v[3], ns[i].v[4], ns[i].v[5];
      PsiIddot << ns[i].a[3], ns[i].a[4], ns[i].a[5];
      tangential_transf(PsiI, T);
      PsiIdot = T.inverse()*Omega;
      tangential_transf_dot(PsiI, PsiIdot, Tdot);
      Alpha = T*PsiIddot + Tdot*PsiIdot;
      for(int j=0; j<3; ++j) ns[i].a[3+j] = Alpha[j];
    }
#endif
  }
}

void
GeomState::setAcceleration(int numNodes, int *nodes, const Vector &a, int SO3param)
{
  // set unconstrained components of acceleration for specified nodes, leaving constrained components unchanged

  int i;
  for(int k = 0; k < numNodes; ++k) {
    i = nodes[k];
    NFrameData *cd = X0->dofFrame(i);
    if(cd) cd->transformVector6(ns[i].a);
    for(int j = 0; j < 6; ++j)
      if(loc[i][j] > -1) {
        ns[i].a[j] = a[loc[i][j]];
      }
    if(cd) cd->invTransformVector6(ns[i].a);
#ifdef USE_EIGEN3
    if(SO3param == 2 && (loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0)) {
      // conversion from second time derivative of total rotation vector to convected angular acceleration
      Eigen::Vector3d PsiI, PsiIdot, PsiIddot, Omega, Alpha;
      Eigen::Matrix3d R, T, Tdot;
      PsiI << ns[i].theta[0], ns[i].theta[1], ns[i].theta[2];
      Omega << ns[i].v[3], ns[i].v[4], ns[i].v[5];
      PsiIddot << ns[i].a[3], ns[i].a[4], ns[i].a[5];
      tangential_transf(PsiI, T);
      PsiIdot = T.inverse()*Omega;
      tangential_transf_dot(PsiI, PsiIdot, Tdot);
      Alpha = T*PsiIddot + Tdot*PsiIdot;
      for(int j=0; j<3; ++j) ns[i].a[3+j] = Alpha[j];
    }
#endif
  }
}

void
GeomState::setVelocityAndAcceleration(const Vector &v, const Vector &a)
{
  // set unconstrained components of velocity and acceleration for all nodes, leaving constrained components unchanged
  // XXX This is not correct for angular velocity/acceleration with one and only one component of the rotation vector constrained,
  // because this constraint implies that the constrained component of the first and second time-derivatives of the rotation vector
  // are zero, not their convected counterparts.
  // Since T^{-1}*Omega = PsiIdot, for example if I know Omega_y and Omega_z and PsiIdot_x then I can compute Omega_x.

  for(int i = 0; i < numnodes; ++i) {
    NFrameData *cd = X0->dofFrame(i);
    if(cd) {
      cd->transformVector6(ns[i].v);
      cd->transformVector6(ns[i].a);
    }
    for(int j = 0; j < 6; ++j)
      if(loc[i][j] > -1) {
        ns[i].v[j] = v[loc[i][j]];
        ns[i].a[j] = a[loc[i][j]];
      }
    if(cd) {
      cd->invTransformVector6(ns[i].v);
      cd->invTransformVector6(ns[i].a);
    }
  }
}

void
GeomState::midpoint_step_update(Vector &vel_n, Vector &acc_n, double delta, GeomState &ss,
                                double beta, double gamma, double alphaf, double alpham, bool zeroRot)
{
  double dt = 2.0*delta;
  double vdcoef, vvcoef, vacoef, avcoef, aacoef;
  vdcoef = (gamma/(dt*beta))/(1-alphaf);
  vvcoef = ((1-(1-alphaf)*gamma/beta)-alphaf)/(1-alphaf);
  vacoef = dt*(2*beta-gamma)/(2*beta);
  avcoef = 1/(dt*gamma);
  aacoef = -(1-gamma)/gamma;

  double d[3];

  int numnodes = std::min(GeomState::numnodes,ss.numnodes);
  for(int i = 0; i < numnodes; ++i) {
    if(flag[i] == -1) continue;

    double d[3] = { ns[i].x - ss[i].x, ns[i].y - ss[i].y, ns[i].z - ss[i].z };

    NFrameData *cd = X0->dofFrame(i);
    if(cd) cd->transformVector3(d);

    // Update translational velocities and accelerations
    if(loc[i][0] >= 0) {
      double v_n = vel_n[loc[i][0]];
      double a_n = acc_n[loc[i][0]];
      vel_n[loc[i][0]] = vdcoef*d[0] + vvcoef*v_n + vacoef*a_n;
      acc_n[loc[i][0]] = avcoef*(vel_n[loc[i][0]] - v_n) + aacoef*a_n;
    }

    if(loc[i][1] >= 0) {
      double v_n = vel_n[loc[i][1]];
      double a_n = acc_n[loc[i][1]];
      vel_n[loc[i][1]] = vdcoef*d[1] + vvcoef*v_n + vacoef*a_n;
      acc_n[loc[i][1]] = avcoef*(vel_n[loc[i][1]] - v_n) + aacoef*a_n;
    }

    if(loc[i][2] >= 0) {
      double v_n = vel_n[loc[i][2]];
      double a_n = acc_n[loc[i][2]];
      vel_n[loc[i][2]] = vdcoef*d[2] + vvcoef*v_n + vacoef*a_n;
      acc_n[loc[i][2]] = avcoef*(vel_n[loc[i][2]] - v_n) + aacoef*a_n;
    }

    // Update angular velocities and accelerations
    if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0) {
      double dtheta[3], dR[3][3];
      if(solInfo.getNLInfo().linearelastic) {
        for(int j=0; j<3; ++j) dtheta[j] = ns[i].theta[j] - ss[i].theta[j];
      }
      else {
        mat_mult_mat(ss[i].R, ns[i].R, dR, 1); // dR = ss[i].R^T * ns[i].R (i.e. ns[i].R = ss[i].R * dR)
        mat_to_vec(dR, dtheta);
      }
      if(cd) cd->transformVector3(dtheta);
      for(int j = 0; j < 3; ++j) {
        if(loc[i][3+j] >= 0) {
          double v_n = vel_n[loc[i][3+j]];
          double a_n = acc_n[loc[i][3+j]];
          vel_n[loc[i][3+j]] = vdcoef*dtheta[j] + vvcoef*v_n + vacoef*a_n;
          acc_n[loc[i][3+j]] = avcoef*(vel_n[loc[i][3+j]] - v_n) + aacoef*a_n;
        }
      }
    }
  }
  setVelocityAndAcceleration(vel_n, acc_n);
  for(int i = 0; i < numnodes; ++i) {
    if(flag[i] == -1) continue;
    for(int j=0; j<6; ++j) {
      ss.ns[i].v[j] = ns[i].v[j];
      ss.ns[i].a[j] = ns[i].a[j];
    }
  }

  // Update step translational displacements
  double tcoef = 1/(1-alphaf);
  for(int i = 0; i < numnodes; ++i) {
    if(flag[i] == -1) continue;
    else if(flag[i] == 0 && ns[i].x == 0) { ss.ns[i].x = 0; continue; }
    ss.ns[i].x = ns[i].x = tcoef*(ns[i].x - alphaf*ss.ns[i].x);
    ss.ns[i].y = ns[i].y = tcoef*(ns[i].y - alphaf*ss.ns[i].y);
    ss.ns[i].z = ns[i].z = tcoef*(ns[i].z - alphaf*ss.ns[i].z);
  }

  // Update step rotational tensor 
  double rcoef  = alphaf/(1-alphaf);
  double result[3][3], result2[3][3], rotVec[3];
  for(int i = 0; i < numnodes; ++i) {
    if(flag[i] == -1 || (loc[i][3] < 0 && loc[i][4] < 0 && loc[i][5] < 0)) continue;
    if(solInfo.getNLInfo().linearelastic) {
      for(int j = 0; j < 3; ++j) 
        ss.ns[i].theta[j] = ns[i].theta[j] = tcoef*(ns[i].theta[j] - alphaf*ss.ns[i].theta[j]);
      vec_to_mat(ns[i].theta, ns[i].R);
      for(int j = 0; j < 3; ++j)
        for(int k = 0; k < 3; ++k)
          ss[i].R[j][k] = ns[i].R[j][k];
    }
    else {
      if(alphaf == 0.0) {
        for(int j = 0; j < 3; ++j) {
          ss[i].theta[j] = ns[i].theta[j];
          for(int k = 0; k < 3; ++k)
            ss[i].R[j][k] = ns[i].R[j][k];
        }
      }
      else {
        mat_mult_mat(ss[i].R, ns[i].R, result2, 1); // result2 = ss[i].R^T * ns[i].R (i.e. ns[i].R = ss[i].R * result2)
        if(alphaf != 0.5) {
          mat_to_vec(result2, rotVec);
          rotVec[0] *= rcoef;
          rotVec[1] *= rcoef;
          rotVec[2] *= rcoef;
          vec_to_mat(rotVec, result2);
        }
        mat_mult_mat(ns[i].R, result2, result, 0); // result = ns[i].R * result2

        for(int j = 0; j < 3; ++j)
          for(int k = 0; k < 3; ++k)
            ss[i].R[j][k] = ns[i].R[j][k] = result[j][k];

        inc_rotvector(ns[i].theta, ns[i].R);
        for(int j=0; j<3; ++j) ss[i].theta[j] = ns[i].theta[j];
      }
    }
  }
#ifdef COMPUTE_GLOBAL_ROTATION
  computeGlobalRotation();
#endif
}

void
GeomState::interp(double alphaf, const GeomState &gs_n, const GeomState &gs_nplus1)
{
  double coef = 1.0 - alphaf;

  // Update displacements: d_{n+1-alphaf} = (1-alphaf)*d_{n+1} + alphaf*d_n
  int inode;
  for(inode=0; inode<numnodes; ++inode) {
    ns[inode].x = alphaf*gs_n[inode].x + coef*gs_nplus1[inode].x;
    ns[inode].y = alphaf*gs_n[inode].y + coef*gs_nplus1[inode].y;
    ns[inode].z = alphaf*gs_n[inode].z + coef*gs_nplus1[inode].z;
  }

  // Update rotations: note this is done using SLERP, interpolating from R_{n+1} to R_n with the incremental rotation defined as a multiplication from the left:
  // R_n = R(l)*R_{n+1} hence R(l) = R_n*R_{n+1}^T and R_{n+1-alphaf} = R(alphaf*l)*R_{n+1}, where l is a rotation "vector" (axis/angle representation)
  // typically SLERP is done using the incremental rotation defined as a multiplication from the right:
  // R_n = R_{n+1}*R(r) hence R(r) = R_{n+1}^T*R_n and R_{n+1-alphaf} = R_{n+1}*R(alphaf*r)
  double rotMat[3][3], rotVec[3];
  for(inode=0; inode<numnodes; ++inode) {
     mat_mult_mat(gs_n[inode].R, gs_nplus1[inode].R, rotMat, 2);
     mat_to_vec(rotMat, rotVec);
     rotVec[0] *= alphaf;
     rotVec[1] *= alphaf;
     rotVec[2] *= alphaf;
     vec_to_mat(rotVec, rotMat);
     mat_mult_mat(rotMat, gs_nplus1[inode].R, ns[inode].R, 0);
  }
}

void
GeomState::diff(const GeomState &un, Vector &vD)
{
  double vec[3], dR[3][3];

  // Loop over all of nodes
  int inode;
  for(inode=0; inode<numnodes; ++inode) {

    // Perform translational dof difference
    if(loc[inode][0] >= 0) vD[loc[inode][0]] = ns[inode].x-un.ns[inode].x;
    if(loc[inode][1] >= 0) vD[loc[inode][1]] = ns[inode].y-un.ns[inode].y;
    if(loc[inode][2] >= 0) vD[loc[inode][2]] = ns[inode].z-un.ns[inode].z;

    // Perform rotational dof difference
    mat_mult_mat(ns[inode].R, un.ns[inode].R, dR, 2);
    mat_to_vec(dR, vec);

    if(loc[inode][3] >= 0) vD[loc[inode][3]] = vec[0];
    if(loc[inode][4] >= 0) vD[loc[inode][4]] = vec[1];
    if(loc[inode][5] >= 0) vD[loc[inode][5]] = vec[2]; 
  }
}

void
GeomState::diff1(const GeomState &un, Vector &vD, int inode)
{
  double vec[3], dR[3][3];

  // Perform translational dof difference
  if(loc[inode][0] >= 0) vD[0] = ns[inode].x-un.ns[inode].x;
  if(loc[inode][1] >= 0) vD[1] = ns[inode].y-un.ns[inode].y;
  if(loc[inode][2] >= 0) vD[2] = ns[inode].z-un.ns[inode].z;

  // Perform rotational dof difference
  mat_mult_mat(ns[inode].R, un.ns[inode].R, dR, 2);
  mat_to_vec(dR, vec);

  if(loc[inode][3] >= 0) vD[3] = vec[0];
  if(loc[inode][4] >= 0) vD[4] = vec[1];
  if(loc[inode][5] >= 0) vD[5] = vec[2];
}

void
GeomState::get_inc_displacement(Vector &incVec, GeomState &ss, bool zeroRot)
{
  // Update incremental translational displacements and rotations
  int inode;
  NFrameData *cd;
  for(inode=0; inode<numnodes; ++inode) {

    if(flag[inode] == -1) continue; // inequality constraint lagrange multiplier dof

    else if(flag[inode] == 0 && inode >= ss.numNodes()) {
      if(loc[inode][0] >= 0) incVec[loc[inode][0]] = ns[inode].x;
      cd = 0;
    }

    else {
      // Update incremental translational displacements
      double d[3] = { ns[inode].x - ss[inode].x,
                      ns[inode].y - ss[inode].y,
                      ns[inode].z - ss[inode].z };
      cd = X0->dofFrame(inode);
      if(cd) cd->transformVector3(d);
      if(loc[inode][0] >= 0) incVec[loc[inode][0]] = d[0];
      if(loc[inode][1] >= 0) incVec[loc[inode][1]] = d[1];
      if(loc[inode][2] >= 0) incVec[loc[inode][2]] = d[2];
    }
    
    if(loc[inode][3] >= 0 || loc[inode][4] >= 0 || loc[inode][5] >= 0) {
      if(zeroRot) {
        // Set rotational displacements equal to zero.
        if(loc[inode][3] >= 0) incVec[loc[inode][3]] = 0.0;
        if(loc[inode][4] >= 0) incVec[loc[inode][4]] = 0.0;
        if(loc[inode][5] >= 0) incVec[loc[inode][5]] = 0.0;
      }
      else {
        double vec[3];
        if(solInfo.getNLInfo().linearelastic) {
          for(int j=0; j<3; ++j) vec[j] = ns[inode].theta[j] - ss[inode].theta[j];
        }
        else {
          double dR[3][3];
          mat_mult_mat( ss[inode].R, ns[inode].R, dR, 1 ); // dR = ss[i].R^T * ns[i].R (i.e. ns[i].R = ss[i].R * dR)
          mat_to_vec( dR, vec );
        }
        if(cd) cd->transformVector3(vec);
        if( loc[inode][3] >= 0 ) incVec[loc[inode][3]] = vec[0];
        if( loc[inode][4] >= 0 ) incVec[loc[inode][4]] = vec[1];
        if( loc[inode][5] >= 0 ) incVec[loc[inode][5]] = vec[2];
      }
    } 
  }

}

void
GeomState::push_forward(Vector &v)
{
#ifdef NEW_PULL_BACK
  // Transform convected quatities (translational only) to spatial frame: v = R*v
  // This new version is required for correct treatment of discrete masses with offsets.
  int inode;
  for(inode=0; inode<numnodes; ++inode) {

    if(flag[inode] == -1) continue;

    if(loc[inode][3] >= 0 || loc[inode][4] >= 0 || loc[inode][5] >= 0) {
      double vec[3], result[3];
      vec[0] = ( loc[inode][0] >= 0 ) ? v[loc[inode][0]] : 0;
      vec[1] = ( loc[inode][1] >= 0 ) ? v[loc[inode][1]] : 0;
      vec[2] = ( loc[inode][2] >= 0 ) ? v[loc[inode][2]] : 0;

      NFrameData *cd = X0->dofFrame(inode);
      if(cd) cd->invTransformVector3(vec);

      mat_mult_vec( ns[inode].R, vec, result, 0 ); // result = R*vec

      if(cd) cd->transformVector3(result);

      if( loc[inode][0] >= 0 ) v[loc[inode][0]] = result[0];
      if( loc[inode][1] >= 0 ) v[loc[inode][1]] = result[1];
      if( loc[inode][2] >= 0 ) v[loc[inode][2]] = result[2];
    }
  }
#endif
}

void
GeomState::pull_back(Vector &v)
{
#ifdef NEW_PULL_BACK
  // Transform spatial quatities (both translational and rotational) to convected frame: v = R^T*v
  // This new version is required for correct treatment of discrete masses with offsets.
  int inode;
  for(inode=0; inode<numnodes; ++inode) {

    if(flag[inode] == -1) continue;

    if(loc[inode][3] >= 0 || loc[inode][4] >= 0 || loc[inode][5] >= 0) {
      double vec[6], result[6];
      vec[0] = ( loc[inode][0] >= 0 ) ? v[loc[inode][0]] : 0;
      vec[1] = ( loc[inode][1] >= 0 ) ? v[loc[inode][1]] : 0;
      vec[2] = ( loc[inode][2] >= 0 ) ? v[loc[inode][2]] : 0;
      vec[3] = ( loc[inode][3] >= 0 ) ? v[loc[inode][3]] : 0;
      vec[4] = ( loc[inode][4] >= 0 ) ? v[loc[inode][4]] : 0;
      vec[5] = ( loc[inode][5] >= 0 ) ? v[loc[inode][5]] : 0;

      NFrameData *cd = X0->dofFrame(inode);
      if(cd) cd->invTransformVector6(vec);

      mat_mult_vec( ns[inode].R, vec, result, 1 ); // result = R^T*vec
      mat_mult_vec( ns[inode].R, vec+3, result+3, 1 );

      if(cd) cd->transformVector6(result);

      if( loc[inode][0] >= 0 ) v[loc[inode][0]] = result[0];
      if( loc[inode][1] >= 0 ) v[loc[inode][1]] = result[1];
      if( loc[inode][2] >= 0 ) v[loc[inode][2]] = result[2];
      if( loc[inode][3] >= 0 ) v[loc[inode][3]] = result[3];
      if( loc[inode][4] >= 0 ) v[loc[inode][4]] = result[4];
      if( loc[inode][5] >= 0 ) v[loc[inode][5]] = result[5];
    }
  }
#else
  // Transform spatial quatities (rotational only) to convected frame: v = R^T*v
  int inode;
  for(inode=0; inode<numnodes; ++inode) {

    if(flag[inode] == -1) continue;

    if(loc[inode][3] >= 0 || loc[inode][4] >= 0 || loc[inode][5] >= 0) {
      double vec[3], result[3];
      vec[0] = ( loc[inode][3] >= 0 ) ? v[loc[inode][3]] : 0;
      vec[1] = ( loc[inode][4] >= 0 ) ? v[loc[inode][4]] : 0;
      vec[2] = ( loc[inode][5] >= 0 ) ? v[loc[inode][5]] : 0;

      NFrameData *cd = X0->dofFrame(inode);
      if(cd) cd->invTransformVector3(vec);

      mat_mult_vec( ns[inode].R, vec, result, 1 ); // result = R^T*vec

      if(cd) cd->transformVector3(result);

      if( loc[inode][3] >= 0 ) v[loc[inode][3]] = result[0];
      if( loc[inode][4] >= 0 ) v[loc[inode][4]] = result[1];
      if( loc[inode][5] >= 0 ) v[loc[inode][5]] = result[2];
    }
  }
#endif
}

void
GeomState::transform(Vector &f, int type, bool unscaled) const
{
  // XXX this function has not been upgraded to support nodal frames yet
#ifdef USE_EIGEN3
  for(int inode = 0; inode < numnodes; ++inode) {

    if(flag[inode] == -1) continue;

    if(loc[inode][3] >= 0 || loc[inode][4] >= 0 || loc[inode][5] >= 0) {
      Eigen::Vector3d vec, result;
      vec[0] = ( loc[inode][3] >= 0 ) ? f[loc[inode][3]] : 0;
      vec[1] = ( loc[inode][4] >= 0 ) ? f[loc[inode][4]] : 0;
      vec[2] = ( loc[inode][5] >= 0 ) ? f[loc[inode][5]] : 0;

      Eigen::Matrix3d R;
      R << ns[inode].R[0][0], ns[inode].R[0][1], ns[inode].R[0][2],
           ns[inode].R[1][0], ns[inode].R[1][1], ns[inode].R[1][2],
           ns[inode].R[2][0], ns[inode].R[2][1], ns[inode].R[2][2];

      Eigen::Vector3d PsiI;
      if(unscaled) {
        PsiI << ns[inode].theta[0], ns[inode].theta[1], ns[inode].theta[2];
      }
      else {
        mat_to_vec<double>(R, PsiI);
      }
      Eigen::Matrix3d T;
      tangential_transf(PsiI, T);

      switch(type) {
        case 0 : // transform time derivative of total rotation vector to convected angular velocity
          result = T*vec;
          break;
        case 1 :   // transform time derivative of total rotation vector to spatial angular velocity
          result = T.transpose()*vec;
          break;
        case 2 :   // transform convected angular velocity to time derivative of total rotation vector
          result = T.inverse()*vec;
          break;
        case 3 :   // transform spatial angular velocity to time derivative of total rotation vector
          result = T.transpose().inverse()*vec;
          break; 
        case 4 : { // transform second time derivative of total rotation vector to convected angular acceleration
                   // assuming ns[inode].v contains the convected angular velocity
          Eigen::Vector3d Omega, PsiIdot;
          Eigen::Matrix3d Tdot;
          Omega << ns[inode].v[3], ns[inode].v[4], ns[inode].v[5];
          PsiIdot = T.inverse()*Omega;
          tangential_transf_dot(PsiI, PsiIdot, Tdot);
          result = T*vec + Tdot*PsiIdot;
        } break;
        case 5 : { // transform second time derivative of total rotation vector to spatial angular acceleration
                   // assuming ns[inode].v contains the convected angular velocity
          Eigen::Vector3d Omega, PsiIdot;
          Eigen::Matrix3d Tdot;
          Omega << ns[inode].v[3], ns[inode].v[4], ns[inode].v[5];
          PsiIdot = T.transpose().inverse()*R*Omega;
          tangential_transf_dot(PsiI, PsiIdot, Tdot);
          result = T.transpose()*vec + Tdot.transpose()*PsiIdot;
        } break;
        case 6 : { // transform convected angular acceleration to second time derivative of total rotation vector
                   // assuming ns[inode].v contains the convected angular velocity
          Eigen::Vector3d Omega, PsiIdot;
          Eigen::Matrix3d Tdot;
          Omega << ns[inode].v[3], ns[inode].v[4], ns[inode].v[5];
          PsiIdot = T.inverse()*Omega;
          tangential_transf_dot(PsiI, PsiIdot, Tdot);
          result = T.inverse()*(vec - Tdot*PsiIdot);
        } break;
        case 7 : { // transform spatial angular acceleration to second time derivative of total rotation vector
                   // assuming ns[inode].v contains the convected angular velocity
          Eigen::Vector3d Omega, PsiIdot;
          Eigen::Matrix3d Tdot;
          Omega << ns[inode].v[3], ns[inode].v[4], ns[inode].v[5];
          PsiIdot = T.transpose().inverse()*R*Omega;
          tangential_transf_dot(PsiI, PsiIdot, Tdot);
          result = T.transpose().inverse()*(vec - Tdot.transpose()*PsiIdot);
        } break;
      }

      if( loc[inode][3] >= 0 ) f[loc[inode][3]] = result[0];
      if( loc[inode][4] >= 0 ) f[loc[inode][4]] = result[1];
      if( loc[inode][5] >= 0 ) f[loc[inode][5]] = result[2];
    }
  }
#endif
}

void
GeomState::transform(Vector &f, const std::vector<int> &weightedNodes, int type, bool unscaled) const
{
  // XXX this function has not been upgraded to support nodal frames yet
#ifdef USE_EIGEN3
  int inode;
  for( std::vector<int>::const_iterator it = weightedNodes.begin(); it != weightedNodes.end(); ++it) {
    inode = *it;

    if(flag[inode] == -1) continue;

    if(loc[inode][3] >= 0 || loc[inode][4] >= 0 || loc[inode][5] >= 0) {
      Eigen::Vector3d vec, result;
      vec[0] = ( loc[inode][3] >= 0 ) ? f[loc[inode][3]] : 0;
      vec[1] = ( loc[inode][4] >= 0 ) ? f[loc[inode][4]] : 0;
      vec[2] = ( loc[inode][5] >= 0 ) ? f[loc[inode][5]] : 0;

      Eigen::Matrix3d R;
      R << ns[inode].R[0][0], ns[inode].R[0][1], ns[inode].R[0][2],
           ns[inode].R[1][0], ns[inode].R[1][1], ns[inode].R[1][2],
           ns[inode].R[2][0], ns[inode].R[2][1], ns[inode].R[2][2];

      Eigen::Vector3d PsiI;
      if(unscaled) {
        PsiI << ns[inode].theta[0], ns[inode].theta[1], ns[inode].theta[2];
      }
      else {
        mat_to_vec<double>(R, PsiI);
      }
      Eigen::Matrix3d T;
      tangential_transf(PsiI, T);

      switch(type) {
        case 0 :
          result = T*vec;
          break;
        case 1 :
          result = T.transpose()*vec;
          break;
        case 2 :
          result = T.inverse()*vec;
          break;
        case 3 :
          result = T.transpose().inverse()*vec;
          break;
      }

      if( loc[inode][3] >= 0 ) f[loc[inode][3]] = result[0];
      if( loc[inode][4] >= 0 ) f[loc[inode][4]] = result[1];
      if( loc[inode][5] >= 0 ) f[loc[inode][5]] = result[2];
    }
  }
#endif
}

void
GeomState::get_tot_displacement(Vector &totVec, bool rescaled)
{
  double x0, y0, z0, vec[3];

  // Loop over all of the nodes
  for(int i = 0; i < numnodes; ++i) {

    // Insert total translational displacements of node i into totVec
    x0 = ((*X0)[i]) ? (*X0)[i]->x : 0;
    vec[0] = ns[i].x - x0;
    y0 = ((*X0)[i]) ? (*X0)[i]->y : 0;
    vec[1] = ns[i].y - y0;
    z0 = ((*X0)[i]) ? (*X0)[i]->z : 0;
    vec[2] = ns[i].z - z0;

    NFrameData *cd = X0->dofFrame(i);
    if(cd) cd->transformVector3(vec);

    if(loc[i][0] >= 0) totVec[loc[i][0]] = vec[0];
    if(loc[i][1] >= 0) totVec[loc[i][1]] = vec[1];
    if(loc[i][2] >= 0) totVec[loc[i][2]] = vec[2];

    // Insert total rotational displacements of node i into totVec
    if(loc[i][3] >= 0 || loc[i][4] >= 0 || loc[i][5] >= 0) {
      if(rescaled) 
        mat_to_vec(ns[i].R, vec);
      else {
        for(int j=0; j<3; ++j) vec[j] = ns[i].theta[j];
      }

      if(cd) cd->transformVector3(vec);

      if(loc[i][3] >= 0) totVec[loc[i][3]] = vec[0];
      if(loc[i][4] >= 0) totVec[loc[i][4]] = vec[1];
      if(loc[i][5] >= 0) totVec[loc[i][5]] = vec[2];
    }
  }
}

void
GeomState::get_temperature(int numNodes, int* nodes, Vector &ndTemps, double Ta)
{
  for(int i=0; i<numNodes; ++i) {
    ndTemps[i] = (ns[nodes[i]].temp == defaultTemp) ? Ta : ns[nodes[i]].temp;
  }
}

void
GeomState::zeroRotDofs(Vector& vec)
{
  for(int inode = 0; inode < numnodes; ++inode) {
    // Set rotational displacements equal to zero.
    if(loc[inode][3] >= 0) vec[loc[inode][3]] = 0.0;
    if(loc[inode][4] >= 0) vec[loc[inode][4]] = 0.0;
    if(loc[inode][5] >= 0) vec[loc[inode][5]] = 0.0;
  }
}

// incremental update of prescribed displacements for nonlinear statics
// i.e. non-zero displacement boundary conditions
void
GeomState::updatePrescribedDisplacement(BCond* dbc, int numDirichlet,
                                        double delta)
{
  // data structure to store rotation vectors of nodes with rotational prescribed dofs
  std::map<int,std::vector<double> > dth;

  // initialize the nodal rotation vectors from the current state
  int i;
  for(i=0; i<numDirichlet; ++i) {
    int nodeNumber = dbc[i].nnum;
    int dofNumber = dbc[i].dofnum;
    if(nodeNumber < ns.size() && dofNumber < 6) {
      std::map<int,std::vector<double> >::iterator it = dth.find(nodeNumber);
      if(it == dth.end()) {
        std::vector<double> v(3);
        v[0] = ns[nodeNumber].theta[0];
        v[1] = ns[nodeNumber].theta[1];
        v[2] = ns[nodeNumber].theta[2];
        dth.insert(it, std::pair<int,std::vector<double> >(nodeNumber,v));
      }
    }
  }

  for(i=0; i<numDirichlet; ++i) {

    int nodeNumber = dbc[i].nnum;
    int dofNumber  = dbc[i].dofnum;
    if(nodeNumber >= ns.size() || dofNumber >= 6) continue;

    // we multiply the total prescribed value by delta which
    // is a parameter prescribed by the user in the input file
    // it is a load control parameter. By default, it is set to 1.0
    // which is applying the full prescribed displacement all at once.

    double prescribedValue = delta*dbc[i].val;

    // if prescribed value is zero, we do nothing.
    if( prescribedValue == 0.0 ) continue;

    NFrameData *cd = X0->dofFrame(nodeNumber);

    switch(dofNumber) {
      case 0:
        if(!cd) {
          ns[nodeNumber].x += prescribedValue;
        }
        else {
          double d[3] = { prescribedValue, 0, 0 };
          cd->invTransformVector3(d);
          ns[nodeNumber].x += d[0];
          ns[nodeNumber].y += d[1];
          ns[nodeNumber].z += d[2];
        }
        break;
      case 1:
        if(!cd) {
          ns[nodeNumber].y += prescribedValue;
        }
        else {
          double d[3] = { 0, prescribedValue, 0 }; 
          cd->invTransformVector3(d);
          ns[nodeNumber].x += d[0];
          ns[nodeNumber].y += d[1];
          ns[nodeNumber].z += d[2];
        }
        break;
      case 2:
        if(!cd) {
          ns[nodeNumber].z += prescribedValue;
        }
        else {
          double d[3] = { 0, 0, prescribedValue };   
          cd->invTransformVector3(d);
          ns[nodeNumber].x += d[0];
          ns[nodeNumber].y += d[1];
          ns[nodeNumber].z += d[2];
        }
        break;
      case 3:
        if(!cd) {
          dth[nodeNumber][0] += prescribedValue;
        }
        else {
          double d[3] = { prescribedValue, 0, 0 };
          cd->invTransformVector3(d);
          dth[nodeNumber][0] += d[0];
          dth[nodeNumber][1] += d[1];
          dth[nodeNumber][2] += d[2];
        }
        break;
      case 4:
        if(!cd) {
          dth[nodeNumber][1] += prescribedValue;
        }
        else {
          double d[3] = { 0, prescribedValue, 0 };
          cd->invTransformVector3(d);
          dth[nodeNumber][0] += d[0];
          dth[nodeNumber][1] += d[1];
          dth[nodeNumber][2] += d[2];
        }
        break;
      case 5:
        if(!cd) {
          dth[nodeNumber][2] += prescribedValue;
        }
        else {
          double d[3] = { 0, 0, prescribedValue };
          cd->invTransformVector3(d);
          dth[nodeNumber][0] += d[0];
          dth[nodeNumber][1] += d[1];
          dth[nodeNumber][2] += d[2];
        }
        break;
      default:
        break;
    }
  }

  // Take care of rotational degrees of freedom for
  // the prescribed displacements
  for(std::map<int,std::vector<double> >::iterator it = dth.begin(); it != dth.end(); ++it) {
    const int &nodeNumber = it->first;
    std::vector<double> &v = it->second;
    form_rottensor( &v[0], ns[nodeNumber].R );
    ns[nodeNumber].theta[0] = v[0];
    ns[nodeNumber].theta[1] = v[1];
    ns[nodeNumber].theta[2] = v[2];
  }
}

// total update of prescribed displacements for nonlinear dynamics
// i.e. displacement boundary and initial conditions prescribed with
// DISP and IDISP
void
GeomState::updatePrescribedDisplacement(BCond* dbc, int numDirichlet,
                                        CoordSet &cs)
{
  // data structure to store rotation and translation vectors of nodes with prescribed dofs
  std::map<int,std::vector<double> > dth, d;

  // initialize the nodal rotation and translation vectors from the current state
  int i;
  for(i=0; i<numDirichlet; ++i) {
    int nodeNumber = dbc[i].nnum;
    int dofNumber  = dbc[i].dofnum;
    if(cs[nodeNumber] && nodeNumber < ns.size()) {
      NFrameData *cd = X0->dofFrame(nodeNumber);
      if(dofNumber == 3 || dofNumber == 4 || dofNumber == 5) {
        std::map<int,std::vector<double> >::iterator it = dth.find(nodeNumber);
        if(it == dth.end()) {
          std::vector<double> v(3);
          v[0] = ns[nodeNumber].theta[0];
          v[1] = ns[nodeNumber].theta[1];
          v[2] = ns[nodeNumber].theta[2];
          if(cd) cd->transformVector3(&v[0]);
          dth.insert(it, std::pair<int,std::vector<double> >(nodeNumber,v));
        }
      }
      else if(cd) {
        std::map<int,std::vector<double> >::iterator it = d.find(nodeNumber);
        if(it == d.end()) {
          std::vector<double> v(3);
          v[0] = ns[nodeNumber].x - cs[nodeNumber]->x;
          v[1] = ns[nodeNumber].y - cs[nodeNumber]->y;
          v[2] = ns[nodeNumber].z - cs[nodeNumber]->z;
          cd->transformVector3(&v[0]);
          d.insert(it, std::pair<int,std::vector<double> >(nodeNumber,v));
        }
      }
    }
  }

  for(i=0; i<numDirichlet; ++i) {

    int nodeNumber = dbc[i].nnum;
    if(!cs[nodeNumber] || nodeNumber >= ns.size()) continue;

    int dofNumber  = dbc[i].dofnum;
    double prescribedValue = dbc[i].val;
    NFrameData *cd = X0->dofFrame(nodeNumber);

    switch(dofNumber) {
      case 0:
        if(!cd) {
          ns[nodeNumber].x = cs[nodeNumber]->x + prescribedValue;
        }
        else {
          d[nodeNumber][0] = prescribedValue;
        }
        break;
      case 1:
        if(!cd) {
          ns[nodeNumber].y = cs[nodeNumber]->y + prescribedValue;
        }
        else {
          d[nodeNumber][1] = prescribedValue;
        }
        break;
      case 2:
        if(!cd) {
          ns[nodeNumber].z = cs[nodeNumber]->z + prescribedValue;
        }
        else {
          d[nodeNumber][2] = prescribedValue;
        }
        break;
      case 3:
        dth[nodeNumber][0] = prescribedValue;
        break;
      case 4:
        dth[nodeNumber][1] = prescribedValue;
        break;
      case 5:
        dth[nodeNumber][2] = prescribedValue;
        break;
      default:
        break;
    }

  }

  for(std::map<int,std::vector<double> >::iterator it = dth.begin(); it != dth.end(); ++it) {
    const int &nodeNumber = it->first;
    std::vector<double> &v = it->second;
    NFrameData *cd = X0->dofFrame(nodeNumber);
    if(cd) cd->invTransformVector3(&v[0]);
    form_rottensor( &v[0], ns[nodeNumber].R );
    ns[nodeNumber].theta[0] = v[0];
    ns[nodeNumber].theta[1] = v[1];
    ns[nodeNumber].theta[2] = v[2];
  }
  for(std::map<int,std::vector<double> >::iterator it = d.begin(); it != d.end(); ++it) {
    const int &nodeNumber = it->first;
    std::vector<double> &v = it->second;
    NFrameData *cd = X0->dofFrame(nodeNumber);
    cd->invTransformVector3(&v[0]);
    ns[nodeNumber].x = cs[nodeNumber]->x + v[0];
    ns[nodeNumber].y = cs[nodeNumber]->y + v[1];
    ns[nodeNumber].z = cs[nodeNumber]->z + v[2];
  }
}

// total update of prescribed displacements and time derivatives for nonlinear dynamics
// i.e. non-zero displacement boundary conditions prescribed with
// USDD which are time dependent user defined displacements
void
GeomState::updatePrescribedDisplacement(double *u, ControlLawInfo *claw,
                                        CoordSet &cs, double *vel, double *acc)
{
  if(claw->numUserDisp == 0) return;

  // data structures to store rotation and translation vectors of nodes with prescribed dofs
  std::map<int,std::vector<double> > dth, d;

  // initialize the nodal rotation and translation vectors from the current state
  int i;
  for(i=0; i<claw->numUserDisp; ++i) {
    int nodeNumber = claw->userDisp[i].nnum;
    int dofNumber = claw->userDisp[i].dofnum;
    if(cs[nodeNumber]) {
      NFrameData *cd = X0->dofFrame(nodeNumber);
      if(dofNumber == 3 || dofNumber == 4 || dofNumber == 5) {
        std::map<int,std::vector<double> >::iterator it = dth.find(nodeNumber);
        if(it == dth.end()) {
          std::vector<double> v(3);
          v[0] = ns[nodeNumber].theta[0];
          v[1] = ns[nodeNumber].theta[1];
          v[2] = ns[nodeNumber].theta[2];
          if(cd) {
            cd->transformVector3(&v[0]);
            cd->transformVector3(&(ns[nodeNumber].v[3]));
            cd->transformVector3(&(ns[nodeNumber].a[3]));
          }
          dth.insert(it, std::pair<int,std::vector<double> >(nodeNumber,v));
        }
      }
      else if(cd) {
        std::map<int,std::vector<double> >::iterator it = d.find(nodeNumber);
        if(it == d.end()) {
          std::vector<double> v(3);
          v[0] = ns[nodeNumber].x - cs[nodeNumber]->x;
          v[1] = ns[nodeNumber].y - cs[nodeNumber]->y;
          v[2] = ns[nodeNumber].z - cs[nodeNumber]->z;
          cd->transformVector3(&v[0]);
          cd->transformVector3(&(ns[nodeNumber].v[0]));
          cd->transformVector3(&(ns[nodeNumber].a[0]));
          d.insert(it, std::pair<int,std::vector<double> >(nodeNumber,v));
        }
      }
    }
  }

  for(i=0; i<claw->numUserDisp; ++i) {
  
    int nodeNumber = claw->userDisp[i].nnum;
    int dofNumber  = claw->userDisp[i].dofnum;
    if(nodeNumber < 0 || nodeNumber >= numnodes) {
      fprintf(stderr, "Bad cntrl law for node number %d(%d) dof(%d)\n", nodeNumber+1,numnodes, claw->userDisp[i].dofnum);
      continue;
    }

    NFrameData *cd = X0->dofFrame(nodeNumber);
    
    switch(dofNumber) {
      case 0: 
        if(!cd) {
          ns[nodeNumber].x = cs[nodeNumber]->x + u[i];
        }
        else {
          d[nodeNumber][0] = u[i];
        }
        ns[nodeNumber].v[0] = vel[i];
        ns[nodeNumber].a[0] = acc[i];
        break;
      case 1:
        if(!cd) {
          ns[nodeNumber].y = cs[nodeNumber]->y + u[i];
        }
        else {
          d[nodeNumber][1] = u[i];
        }
        ns[nodeNumber].v[1] = vel[i];
        ns[nodeNumber].a[1] = acc[i];
        break;
      case 2:
        if(!cd) {
          ns[nodeNumber].z = cs[nodeNumber]->z + u[i];
        }
        else {
          d[nodeNumber][2] = u[i];
        }
        ns[nodeNumber].v[2] = vel[i];
        ns[nodeNumber].a[2] = acc[i];
        break;
      case 3:
        dth[nodeNumber][0] = u[i];
        ns[nodeNumber].v[3] = vel[i];
        ns[nodeNumber].a[3] = acc[i];
        break;
      case 4:
        dth[nodeNumber][1] = u[i];
        ns[nodeNumber].v[4] = vel[i];
        ns[nodeNumber].a[4] = acc[i];
        break;
      case 5:
        dth[nodeNumber][2] = u[i];
        ns[nodeNumber].v[5] = vel[i];
        ns[nodeNumber].a[5] = acc[i];
        break;
      default:
        break;
    }

  }

  for(std::map<int,std::vector<double> >::iterator it = dth.begin(); it != dth.end(); ++it) {
    const int &nodeNumber = it->first;
    std::vector<double> &v = it->second;
    NFrameData *cd = X0->dofFrame(nodeNumber);
    if(cd) {
      cd->invTransformVector3(&v[0]);
      cd->invTransformVector3(&(ns[nodeNumber].v[3]));
      cd->invTransformVector3(&(ns[nodeNumber].a[3]));
    }
    form_rottensor( &v[0], ns[nodeNumber].R );
    ns[nodeNumber].theta[0] = v[0];
    ns[nodeNumber].theta[1] = v[1];
    ns[nodeNumber].theta[2] = v[2];
  }
  for(std::map<int,std::vector<double> >::iterator it = d.begin(); it != d.end(); ++it) {
    const int &nodeNumber = it->first;
    std::vector<double> &v = it->second;
    NFrameData *cd = X0->dofFrame(nodeNumber);
    cd->invTransformVector3(&v[0]);
    cd->invTransformVector3(&(ns[nodeNumber].v[0]));
    cd->invTransformVector3(&(ns[nodeNumber].a[0]));
    ns[nodeNumber].x = cs[nodeNumber]->x + v[0];
    ns[nodeNumber].y = cs[nodeNumber]->y + v[1];
    ns[nodeNumber].z = cs[nodeNumber]->z + v[2];
  }
}

void
GeomState::getPositions(double *positions)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   positions[3*i+0] = ns[i].x;
   positions[3*i+1] = ns[i].y;
   positions[3*i+2] = ns[i].z;
 }
}

void
GeomState::getRotations(double *rotations)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   rotations[12*i+0] = ns[i].R[0][0];
   rotations[12*i+1] = ns[i].R[0][1];
   rotations[12*i+2] = ns[i].R[0][2];
   rotations[12*i+3] = ns[i].R[1][0];
   rotations[12*i+4] = ns[i].R[1][1];
   rotations[12*i+5] = ns[i].R[1][2];
   rotations[12*i+6] = ns[i].R[2][0];
   rotations[12*i+7] = ns[i].R[2][1];
   rotations[12*i+8] = ns[i].R[2][2];
   rotations[12*i+9] = ns[i].theta[0];
   rotations[12*i+10] = ns[i].theta[1];
   rotations[12*i+11] = ns[i].theta[2];
 }
}

void
GeomState::getVelocities(double *velocities)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   velocities[6*i+0] = ns[i].v[0];
   velocities[6*i+1] = ns[i].v[1];
   velocities[6*i+2] = ns[i].v[2];
   velocities[6*i+3] = ns[i].v[3];
   velocities[6*i+4] = ns[i].v[4];
   velocities[6*i+5] = ns[i].v[5];
 }
}

void
GeomState::getAccelerations(double *accelerations)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   accelerations[6*i+0] = ns[i].a[0];
   accelerations[6*i+1] = ns[i].a[1];
   accelerations[6*i+2] = ns[i].a[2];
   accelerations[6*i+3] = ns[i].a[3];
   accelerations[6*i+4] = ns[i].a[4];
   accelerations[6*i+5] = ns[i].a[5];
 }
}

void
GeomState::getElemStates(double *elemStates) const
{
 int i,j,k;
 for(i=0,j=0; i<numelems; ++i) {
   for(k=0; k<es[i].numInternalStates; ++k)
     elemStates[j++] = es[i].internalStates[k];
 }
}

void
GeomState::setPositions(double *positions)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   ns[i].x = positions[3*i+0];
   ns[i].y = positions[3*i+1];
   ns[i].z = positions[3*i+2];
 }
}

void
GeomState::setRotations(double *rotations)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   ns[i].R[0][0] = rotations[12*i+0];
   ns[i].R[0][1] = rotations[12*i+1];
   ns[i].R[0][2] = rotations[12*i+2];
   ns[i].R[1][0] = rotations[12*i+3];
   ns[i].R[1][1] = rotations[12*i+4];
   ns[i].R[1][2] = rotations[12*i+5];
   ns[i].R[2][0] = rotations[12*i+6];
   ns[i].R[2][1] = rotations[12*i+7];
   ns[i].R[2][2] = rotations[12*i+8];
   ns[i].theta[0] = rotations[12*i+9];
   ns[i].theta[1] = rotations[12*i+10];
   ns[i].theta[2] = rotations[12*i+11];
 }
}

void
GeomState::setVelocities(double *velocities)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   ns[i].v[0] = velocities[6*i+0];
   ns[i].v[1] = velocities[6*i+1];
   ns[i].v[2] = velocities[6*i+2];
   ns[i].v[3] = velocities[6*i+3];
   ns[i].v[4] = velocities[6*i+4];
   ns[i].v[5] = velocities[6*i+5];
 }
}

void
GeomState::setAccelerations(double *accelerations)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   ns[i].a[0] = accelerations[6*i+0];
   ns[i].a[1] = accelerations[6*i+1];
   ns[i].a[2] = accelerations[6*i+2];
   ns[i].a[3] = accelerations[6*i+3];
   ns[i].a[4] = accelerations[6*i+4];
   ns[i].a[5] = accelerations[6*i+5];
 }
}

void
GeomState::setNodalTemperatures(double *ndTemps)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   ns[i].temp = (ndTemps) ? ndTemps[i] : defaultTemp;
 }
}

void
GeomState::setElemStates(double *elemStates)
{
 int i,j,k;
 for(i=0,j=0; i<numelems; ++i) {
   for(k=0; k<es[i].numInternalStates; ++k)
     es[i].internalStates[k] = elemStates[j++];
 }
}

int
GeomState::getTotalNumElemStates() const
{
 int n = 0;
 for(int i=0; i<numelems; ++i) {
   n += es[i].numInternalStates;
 }
 return n;
}

void
GeomState::addMultiplierNode(const std::pair<int,int> &lmpc_id, double value)
{
  NodeState n;
  n.x = value;
  ns.push_back(n);
  multiplier_nodes[lmpc_id] = numnodes++;
}

double
GeomState::getMultiplier(const std::pair<int,int> &lmpc_id)
{
  std::map<std::pair<int,int>,int>::iterator it = multiplier_nodes.find(lmpc_id);
  return (it != multiplier_nodes.end()) ? ns[it->second].x : 0;
}

void
GeomState::getMultipliers(std::map<std::pair<int,int>,double> &mu)
{
  for(std::map<std::pair<int,int>,double>::iterator it = mu.begin(); it != mu.end(); it++) {
    it->second = getMultiplier(it->first);
  }
}

void
GeomState::setMultiplier(const std::pair<int,int> &lmpc_id, double mu)
{
  std::map<std::pair<int,int>,int>::iterator it = multiplier_nodes.find(lmpc_id);
  if(it != multiplier_nodes.end()) ns[it->second].x = mu;
}

void
GeomState::setMultipliers(std::map<std::pair<int,int>,double> &mu)
{
  for(std::map<std::pair<int,int>,double>::iterator it = mu.begin(); it != mu.end(); it++) {
    setMultiplier(it->first, it->second);
  }
}

void GeomState::computeGlobalRotation() 
{
  double cg[3];
  computeCG(cg);

  double deltaRot[3][3];  // incremental rotation
  // initialize deltaRot to Identity Matrix
  double zeroRot[3] = {0,0,0};
  computeRotMat(zeroRot, deltaRot);

  int iter = 10;  // maximum number of iterations
  double tol = 1.0e-12; // convergence tolerance
  double jac[3][3], grad[3];  //minimization gradient and jacobian
  double l2;
  bool converged = false;

  for (int it = 0; it < iter; it++)  {
    int i,j;
    // compute minimization gradients and jacobians
    computeRotGradAndJac(cg, grad, jac);

   // check for convergence
    if(it > 0 && sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]) <= tol) {
       converged = true;
       break;
    }

    // rotation results come back in grad
    solve(jac, grad);

    // update deltaRot
    computeRotMat(grad, deltaRot);

    // update global rotation matrix
    double R[3][3];
    mat_mult_mat(deltaRot, gRot, R, 0);
  
    for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
        gRot[i][j] = R[i][j];
  }
/*
#ifndef NDEBUG
  if(!converged) {
    computeRotGradAndJac(cg, grad, jac);
    double l2 = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
    if(l2 > tol)
      std::cerr << "failed to converge in computeGlobalRotation, l2 norm = " << l2 << std::endl;
  }
#endif
*/
} 

void GeomState::computeRotMat(double *angle, double mat[3][3])
{
  // trig functions of angles
  double c1 = cos(angle[0]);
  double s1 = sin(angle[0]);
  double c2 = cos(angle[1]);
  double s2 = sin(angle[1]);
  double c3 = cos(angle[2]);
  double s3 = sin(angle[2]);


  /* Compute rotation matrix
     computed as R1.R2.R3
     where R1 is rotation about z
           R2 is rotation about y
           R3 is rotation about x
  mat[0][0] = c1*c2;
  mat[0][1] = c1*s2*s3 - c3*s1;
  mat[0][2] = c1*c3*s2 + s1*s3;
  
  mat[1][0] = c2*s1;
  mat[1][1] = c1*c3+s1*s2*s3;
  mat[1][2] = c3*s1*s2 - c1*s3;
  
  mat[2][0] = -s2;
  mat[2][1] = c2*s3;
  mat[2][2] = c2*c3; */
  

  /* computed as R1.R2.R3
     where R1 is rotation about x
           R2 is rotation about y
           R3 is rotation about z
  */ 
  mat[0][0] = c2*c3;
  mat[0][1] = -c2*s3;
  mat[0][2] = s2;

  mat[1][0] = c3*s1*s2 + c1*s3;
  mat[1][1] = c1*c3 - s1*s2*s3;
  mat[1][2] = -c2*s1;

  mat[2][0] = s1*s3 - c1*c3*s2;
  mat[2][1] = c3*s1 + c1*s2*s3;
  mat[2][2] = c1*c2;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
void GeomState::solve(double m[3][3], double v[3]) 
{
#ifdef USE_EIGEN3
  Eigen::Matrix3d mat(3,3);
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) mat(i,j) = m[i][j];

  Eigen::Map<Eigen::Vector3d> vec(v);
  mat.selfadjointView<Eigen::Lower>().ldlt().solveInPlace(vec);
#else
  int i,j,k;

  for (i = 0; i < 2; i++)
    for (j = i+1; j < 3; ++j) {
      if (m[i][i] == 0.0)  {
        for (k = 0; k < 3; k++)  m[i][k] = 0;
        m[i][i] = 1.0;
        v[i] = 0;
      }
      double coef = m[j][i]/m[i][i];
      for (k = i+1; k < 3; ++k)
        m[j][k] -= coef*m[i][k];
      v[j] -= coef*v[i];
    }
//NOT VERY CLEAN
  if (m[2][2] == 0.0)  {
    for (k = 0; k < 2; k++)  m[2][k] = 0;
      m[2][2] = 1.0;
      v[2] = 0;
    }

  for (i=2; i >= 0; i--) {
    for (j=2; j > i; j--)
      v[i] -= m[i][j] * v[j];
    v[i] /= m[i][i];

  }
#endif
}

void GeomState::computeRotGradAndJac(double cg[3], double grad[3], double jac[3][3]) 
{
  // init grad and jac
  int i,j;
  for (i = 0; i < 3; i++)  {
    grad[i] = 0.0;
    for (j = 0; j < 3; j++)
      jac[i][j] = 0.0;
  }

  // rotate local vectors using R(n-1)
  for (i = 0; i < numnodes; i++) {
    if (flag[i] <= 0) continue;

    double rd[3];
    rd[0] = (*X0)[i]->x - refCG[0]; 
    rd[1] = (*X0)[i]->y - refCG[1]; 
    rd[2] = (*X0)[i]->z - refCG[2]; 
    rotate(gRot, rd);

    // compute freq. used values
    double eVec[3];  // Deformed X - rotated loc. vec - deformed cg
    eVec[0] = ns[i].x - cg[0] - rd[0];
    eVec[1] = ns[i].y - cg[1] - rd[1];
    eVec[2] = ns[i].z - cg[2] - rd[2];

    // compute gradients 
    // dg = 2*evec dot rotation gradient
    grad[0] += 2*(eVec[2]*rd[1] - eVec[1]*rd[2]);
    grad[1] += 2*(eVec[0]*rd[2] - eVec[2]*rd[0]);
    grad[2] += 2*(eVec[1]*rd[0] - eVec[0]*rd[1]);

    // compute jacobian
    jac[0][0] += 2 * (rd[1]*eVec[1]
                   +  rd[2]*eVec[2]
                   +  rd[2]*rd[2]
 		   +  rd[1]*rd[1]);
  
    jac[0][1] += -2 * (rd[0]*eVec[1]
                    +  rd[1]*rd[0]);

    jac[0][2] += -2 * (rd[0]*eVec[2]
                    +  rd[0]*rd[2]);

    jac[1][1] +=  2 * (rd[0]*eVec[0]
                    +  rd[2]*eVec[2]
                    +  rd[2]*rd[2]
		    +  rd[0]*rd[0]);

    jac[1][2] += -2 * (rd[1]*eVec[2]
                    +  rd[2]*rd[1]);

    jac[2][2] +=  2 * (rd[0]*eVec[0]
                    +  rd[1]*eVec[1]
		    +  rd[1]*rd[1]
  		    +  rd[0]*rd[0]);
  }

  // Fill in symmetric terms
  jac[1][0] = jac[0][1];
  jac[2][0] = jac[0][2];
  jac[2][1] = jac[1][2];
}

void GeomState::computeCG(double cg[3])
{
  // init cg
  cg[0] = 0.0;
  cg[1] = 0.0;
  cg[2] = 0.0;

  int i, numReal = 0;
  for ( i = 0; i < numnodes; ++i)  {
    if (flag[i] == 1)  {
      cg[0] += ns[i].x;
      cg[1] += ns[i].y;
      cg[2] += ns[i].z;
      numReal++;
    }
  }

  double invTotNd = 1.0 / numReal;
  for (i = 0; i < 3; ++i)
    cg[i] *= invTotNd;
}

void GeomState::rotate(double R[3][3], double v[3]) 
{
  double c[3];
  for (int j = 0; j < 3; j++)
    c[j] = R[j][0]*v[0] +
           R[j][1]*v[1] +
           R[j][2]*v[2];
 
  v[0] = c[0];
  v[1] = c[1];
  v[2] = c[2];
}

void GeomState::getGlobalRot(double R[3][3]) 
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      R[i][j] = gRot[i][j];
}

void GeomState::transformCoords(double xScaleFactor, double yScaleFactor, double zScaleFactor)
{
  for(int i = 0; i < numnodes; ++i) {
    if(X0 && (*X0)[i]) {
      (*X0)[i]->x *= xScaleFactor;
      (*X0)[i]->y *= yScaleFactor;
      (*X0)[i]->z *= zScaleFactor;
    }
  }
}

void
GeomState::setNewCoords(const Vector &X)
{
  for(int i = 0; i < numnodes; ++i) {
    if(X0 && (*X0)[i]) {
      (*X0)[i]->x = X[6*i+0];
      (*X0)[i]->y = X[6*i+1];
      (*X0)[i]->z = X[6*i+2];
    }
  }
}

double
NodeState::diff(const Node &un, int dof)
{
  double delta = 0.0;

  // Loop over all of nodes
  switch(dof) {
    case 0:
    {
      delta = x-un.x;
      break;
    }
    case 1:
    {
      delta = y-un.y;
      break;
    }
    case 2:
    {
      delta = z-un.z;
      break;
    }
    case 3:
    case 4:
    case 5:
    {
      double vec[3];
      mat_to_vec(R, vec);
      delta = vec[dof-3];
      break;
    }
  }
  return delta;
}

#include <Corotational.d/TemperatureState.h>
TemperatureState::TemperatureState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs)
 : GeomState(cs)
{
  numnodes = dsa.numNodes();            // Number of nodes
  ns.resize(numnodes);
  loc.resize(numnodes);
  flag.resize(numnodes);

  int i;
  for(i=0; i<numnodes; ++i) {

    loc[i].resize(6);

    // Store location of each degree of freedom
    loc[i][0] = cdsa.locate( i, DofSet::Temp );
    for(int j=1; j<6; ++j) loc[i][j] = -1;

    // Get Node i from the Coordinate (Node) set
    Node *node_i = cs[i];

    if (node_i)  {

      ns[i].x = 0;
      if(dsa[i].list() != 0)  {
        flag[i] = 1;
      }
      else
        flag[i] = 0;

    }
    else  {
      ns[i].x = 0.0;
      int dof;
      if((dof = cdsa.locate( i, DofSet::LagrangeE )) > -1) {
        loc[i][0] = dof;
        flag[i] = 0;
      }
      else if((dof = cdsa.locate( i, DofSet::LagrangeI )) > -1) {
        loc[i][0] = dof;
        flag[i] = -1;
      }
      else 
        flag[i] = 0;
    }

  }

}

TemperatureState::TemperatureState(const TemperatureState &g2) : GeomState(*(g2.X0))
{
  // Copy number of nodes
  numnodes = g2.numnodes;

  // Allocate memory for node states & dof locations
  ns.resize(numnodes);
  loc.resize(numnodes);

  // flag for node to element connectivity
  flag.resize(numnodes);

  // Copy dof locations
  int i;
  for(i = 0; i < numnodes; ++i) {
    loc[i].resize(6);
    loc[i][0] = g2.loc[i][0];
    for(int j=1; j<6; ++j) loc[i][j] = -1;
  }

  // Copy node states
  for(i = 0; i < numnodes; ++i) {
    ns[i]  = g2.ns[i];
    flag[i]= g2.flag[i];
  }
}


void
TemperatureState::update(const Vector &v, int)
{
  // v = incremental displacement vector

  int i;
  for(i=0; i<numnodes; ++i) {

    // Set incremental translational displacements
    double dx = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;

    // Increment total translational displacements
    ns[i].x += dx;

  }
}

void
TemperatureState::explicitUpdate(CoordSet &cs, const Vector &v)
{
  // v = total displacement vector (unconstrained dofs only)

  int i;
  for(i=0; i<numnodes; ++i) {

    // Set total displacements (unconstrained dofs only)
    if (loc[i][0] >= 0) ns[i].x = v[loc[i][0]];

  }
}

void
TemperatureState::updatePrescribedDisplacement(BCond* dbc, int numDirichlet,
                                               double delta)
{
  int i;
  for(int i=0; i<numDirichlet; ++i) {

    int nodeNumber = dbc[i].nnum;
    int dofNumber  = dbc[i].dofnum;


    // we multiply the total prescribed value by delta which
    // is a parameter prescribed by the user in the input file
    // it is a load control parameter. By default, it is set to 1.0
    // which is applying the full prescribed displacement all at once.

    double prescribedValue = delta*dbc[i].val;

    // if prescribed value is zero, we do nothing.
    if( prescribedValue == 0.0 ) continue;

    switch(dofNumber) {
        case 6:
                ns[nodeNumber].x += prescribedValue;
                break;
        default:
                break;
    }

  }
}

// update prescribed displacements for nonlinear dynamics
// i.e. displacement boundary and initial conditions prescribed with
// TEMP and ITEMP
void
TemperatureState::updatePrescribedDisplacement(BCond* dbc, int numDirichlet,
                                               CoordSet &cs)
{
  int i;
  for(i=0; i<numDirichlet; ++i) {

    int nodeNumber = dbc[i].nnum;
    if(!cs[nodeNumber]) continue;

    int dofNumber  = dbc[i].dofnum;

    double prescribedValue = dbc[i].val;

    switch(dofNumber) {
        case 6:
                ns[nodeNumber].x = prescribedValue;
                break;
        default:
                break;
    }

  }
}

void
TemperatureState::get_inc_displacement(Vector &incVec, GeomState &ss, bool)
{
  int inode;
  for(inode=0; inode<numnodes; ++inode) {
    // Update incremental temperature
    if(loc[inode][0] >= 0) incVec[loc[inode][0]] = ns[inode].x - ss[inode].x;
  }
}

void
TemperatureState::midpoint_step_update(Vector &vel_n, Vector &accel_n, double delta, GeomState &ss,
                                       double beta, double gamma, double alphaf, double alpham, bool)
{
 // XXX THIS HASN'T BEEN UPDATED FOR GENERALIZED ALPHA YET
 // Update incremental temperatures at end of step:
 double coef = 2.0/delta;

 // x,y,z velocity step update (velocity_n+1 = 2.0*velocity_n+1/2 - velocity_n)
 // note: we are computing translational velocity_n+1/2 locally
 //       as velocity_n+1/2 = (ns[i] - ss.ns[i])/delta = 2/dt*(u_n+1/2 - u_n)

 int i;
 for(i=0; i<numnodes; ++i) {
   if(loc[i][0] >= 0)
     ns[i].v[0] = vel_n[loc[i][0]] = coef*(ns[i].x - ss[i].x) - vel_n[loc[i][0]];
     ss[i].v[0] = ns[i].v[0];
 }

 // Update step translational displacements
 int inode;
 for(inode=0; inode<numnodes; ++inode) {
   ns[inode].x    = 2.0*ns[inode].x - ss[inode].x;
   ss[inode].x = ns[inode].x;
 }
}

void
TemperatureState::setVelocityAndAcceleration(const Vector &v, const Vector &)
{
  for(int i = 0; i < numnodes; ++i)
    if(loc[i][0] > -1) {
      ns[i].v[0] = v[loc[i][0]];
    }
}

void
TemperatureState::setVelocity(const Vector &v, int)
{
  for(int i = 0; i < numnodes; ++i)
    if(loc[i][0] > -1) {
      ns[i].v[0] = v[loc[i][0]];
    }
}

void
TemperatureState::setAcceleration(const Vector &, int)
{
}

void
TemperatureState::get_tot_displacement(Vector &totVec, bool)
{
  // Loop over all of the nodes
  for(int i = 0; i < numnodes; ++i) {

    if(loc[i][0] >= 0) totVec[loc[i][0]] = ns[i].x;

  }
}

