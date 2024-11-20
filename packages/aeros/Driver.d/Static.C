#include <cstdio>
#include <cmath>
#include <Utils.d/dbg_alloca.h>

#include <Corotational.d/Corotator.h>
#include <Driver.d/Domain.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Solvers.d/Rbm.h>
#include <Solvers.d/Solver.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/Memory.h>
#include <Utils.d/pstress.h>
#include <Timers.d/GetTime.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/DistHelper.h>
#ifdef USE_EIGEN3
#include <Math.d/EiSparseMatrix.h>
#endif
#include <Math.d/CuCSparse.h>
#include <Rom.d/PodProjectionSolver.h>

#include <list>
#include <numeric>
#include <set>
#include <fstream>
#include <iomanip>
#include <string>
#include <limits>

extern int verboseFlag;

void
Domain::buildPrescDisp(Vector &pDisp, double lambda)
{
 int i;
 pDisp.zero();
 // Set the dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
  int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) {
     //fprintf(stderr," *** WARNING: dof does not exist: node %d dof %d ***\n",
     //                dbc[i].nnum+1,dbc[i].dofnum+1);
     continue;
   }

  int cdof = c_dsa->invRCN(dof);
  pDisp[cdof] = lambda*dbc[i].val;
 }
}

void
Domain::buildPrescDisp(Vector &pDisp, double t, double)
{
 int i;
 pDisp.zero();
 // Set the dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
  int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) {
     //fprintf(stderr," *** WARNING: dof does not exist: node %d dof %d ***\n",
     //                dbc[i].nnum+1,dbc[i].dofnum+1);
     continue;
   }

  int cdof = c_dsa->invRCN(dof);
  pDisp[cdof] = t*dbc[i].val;
 }
}

double *
Domain::getNodalTemperatures()
{
  return temprcvd;
}

void
Domain::initNodalTemperatures()
{
  if(sinfo.thermalLoadFlag && sinfo.thermoeFlag == -1) {

    temprcvd = new double[numnodes];

    int i;
    for(i = 0; i < numnodes; ++i)
      temprcvd[i] = defaultTemp;

    for(i = 0; i < numDirichlet; ++i)
      if((1 << dbc[i].dofnum) == DofSet::Temp && dbc[i].nnum < numnodes) {
        temprcvd[dbc[i].nnum] = dbc[i].val;
      }
  }
}

FILE *
Domain::openFile(char *fileName, const char *extension)
{
 // Open decomposition file
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 int len1 = strlen(cinfo->checkfile);
 int len2 = strlen(extension);
 char *file = (char *) dbg_alloca(sizeof(char)*(len1+len2+1));
 strcpy(file, cinfo->checkfile);
 strcat(file, extension);

 FILE *filePtr;
 if((filePtr= fopen(file,"w")) == (FILE *) 0 ) {
   fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n",file );
   exit(0);
 }

 return filePtr;
}

void
Domain::printStatistics(bool domain_decomp)
{
  // Note 1: # Nodes does include nodes which are not attached to any element, but does not
  //         include the so-called "internal nodes" (e.g. Lagrange multiplier nodes).
  // Note 2: # Elements does not include "LMPC elements".
  // Note 3: # Unconstrained dofs does not include dofs of nodes which are not attached to any
  //         elements or LMPC elements.
  // Note 4: # Constrained dofs and Total # dofs do include constraints on dofs of nodes which
  //         (a) are not defined in the input file, and (b) are defined in the input file but
  //         not attached to any elements or LMPC elements.
  if(!nodeTable) {
    exactNumNodes = 0;
    for(int i=0; i<nodes.size(); ++i) {
      if(nodes[i] != 0) exactNumNodes++;
    }
  }
  if(!domain_decomp) setNumDofs(numUncon()+numDirichlet+numComplexDirichlet);

  filePrint(stderr, "\n --------- PROBLEM PARAMETERS ---------");
  filePrint(stderr, "\n ... # Nodes              = %7d ...", exactNumNodes);
  filePrint(stderr, "\n ... # Elements           = %7d ...", numele-geoSource->numMpcElem()-contactSurfElems.size()-geoSource->numContactSurfaceElem());
  filePrint(stderr, "\n ... # Unconstrained dofs = %7d ...", numDofs()-numDirichlet-numComplexDirichlet);
  filePrint(stderr, "\n ... # Constrained dofs   = %7d ...", numDirichlet+numComplexDirichlet);
  filePrint(stderr, "\n ... Total # dofs         = %7d ...", numDofs());
  //filePrint(stderr, "\n ... # Loaded dofs        = %7d ...", numNeuman+numComplexNeuman);
  if(gravityFlag())
    filePrint(stderr, "\n ... Gravity Load is Applied        ...");
  filePrint(stderr, "\n ... # Output Files       = %7d ...", geoSource->getNumOutInfo());
  filePrint(stderr, "\n --------------------------------------\n");
}

double
Domain::computeStructureMass(bool printFlag, int groupId)
{
  // Compute total mass and mass center of gravity
  // Added calculation for moments of inertia

  double totmas = 0.0; // total mass of structural system (including fluid elemets and discrete masses)
  double xc     = 0.0;
  double yc     = 0.0;
  double zc     = 0.0;
  double Mx     = 0.0;
  double My     = 0.0;
  double Mz     = 0.0;
  double Ixx    = 0.0;
  double Iyy    = 0.0;
  double Izz    = 0.0;
  double Ixy    = 0.0;
  double Iyz    = 0.0;
  double Ixz    = 0.0;
  maxNumNodes = 0;

  // compute max number nodes per element
  int iele, i;
  for(iele=0; iele < numele; ++iele) {
    int numNodesPerElement = packedEset[iele]->numNodes();
    maxNumNodes = std::max(maxNumNodes, numNodesPerElement);
  }

  // allocate one array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  for(iele=0; iele < numele; ++iele) {

    int numNodesPerElement = packedEset[iele]->numNodes();
    packedEset[iele]->nodes(nodeNumbers);

    if (groupId > 0) {
      std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
      std::set<int>::iterator it;
      for (int iNode = 0; iNode < numNodesPerElement; ++iNode)
        if((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
      if(it == groupNodes.end()) continue;
    }

    double elementMass = packedEset[iele]->getMass(nodes);
    totmas += elementMass;
    int numRealNodes = 0;
    for(i=0; i<numNodesPerElement; ++i)
      if(nodes[nodeNumbers[i]]) numRealNodes++;

    double massPerNode = elementMass/numRealNodes;

    for(i=0; i<numNodesPerElement; ++i) {
       if(nodes[nodeNumbers[i]] == 0) continue;
       Node &node = nodes.getNode(nodeNumbers[i]);

       double x = node.x;
       double y = node.y;
       double z = node.z;

       xc  += massPerNode*x;
       yc  += massPerNode*y;
       zc  += massPerNode*z;

       Ixx += massPerNode*(y*y + z*z);
       Iyy += massPerNode*(x*x + z*z);
       Izz += massPerNode*(x*x + y*y);

       Ixy -= massPerNode*x*y;
       Ixz -= massPerNode*x*z;
       Iyz -= massPerNode*y*z;
    }
  }

  delete [] nodeNumbers;

  // now add the mass from the fluid elements (these are in a different element set)
  if(geoSource->numElemFluid() > 0 && !(groupId > 0)) {
    for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
      int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
      maxNumNodes = std::max(maxNumNodes, numNodesPerElement);
    }

    nodeNumbers = new int[maxNumNodes];
  
    double fluidmas = 0.0;
    for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
      double elementMass = (*(geoSource->getPackedEsetFluid()))[iele]->getMass(nodes);
      fluidmas += elementMass;
      int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
      (*(geoSource->getPackedEsetFluid()))[iele]->nodes(nodeNumbers);
      int numRealNodes = 0;
      for(i=0; i<numNodesPerElement; ++i)
        if(nodes[nodeNumbers[i]]) numRealNodes++;

      double massPerNode = elementMass/numRealNodes;

      for(i=0; i<numNodesPerElement; ++i) {
         if(nodes[nodeNumbers[i]] == 0) continue;
         Node &node = nodes.getNode(nodeNumbers[i]);
  
         double x = node.x;
         double y = node.y;
         double z = node.z;

         xc  += massPerNode*x;
         yc  += massPerNode*y;
         zc  += massPerNode*z;

         Ixx += massPerNode*(y*y + z*z);
         Iyy += massPerNode*(x*x + z*z);
         Izz += massPerNode*(x*x + y*y);

         Ixy -= massPerNode*x*y;
         Ixz -= massPerNode*x*z;
         Iyz -= massPerNode*y*z;
      }
    }

    delete [] nodeNumbers;

    if(printFlag) {
      filePrint(stderr," Fluid Mass = %e\n",fluidmas);
      filePrint(stderr," --------------------------------------\n");
    }

    totmas += fluidmas;
  
  }

  double Iyx = Ixy;
  double Izx = Ixz;
  double Izy = Iyz;

  // Add discrete masses
  DMassData *current = firstDiMass;
  while(current != 0) {
    int n = current->node;

    if (groupId > 0) {
      std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
      if(groupNodes.find(n) == groupNodes.end()) continue;
    }

    // check if the node is in the model
    //if(nodes.exist(n)) {
      Node &node = nodes.getNode(n);
      switch(current->dof) {
        case 0: {
          Mx += current->diMass;
          xc += current->diMass*node.x;
          Iyy += current->diMass*(node.z*node.z);
          Izz += current->diMass*(node.y*node.y);
          Ixy -= 0.5*current->diMass*(node.x*node.y);
          Iyx -= 0.5*current->diMass*(node.y*node.x);
          Ixz -= 0.5*current->diMass*(node.x*node.z);
          Izx -= 0.5*current->diMass*(node.z*node.x);
        } break;
        case 1: {
          My += current->diMass;
          yc += current->diMass*node.y;
          Ixx += current->diMass*(node.z*node.z);
          Izz += current->diMass*(node.x*node.x);
          Iyx -= 0.5*current->diMass*(node.x*node.y);
          Ixy -= 0.5*current->diMass*(node.y*node.x);
          Iyz -= 0.5*current->diMass*(node.y*node.z);
          Izy -= 0.5*current->diMass*(node.z*node.y);
        } break;
        case 2: {
          Mz += current->diMass;
          zc += current->diMass*node.z;
          Ixx += current->diMass*(node.y*node.y);
          Iyy += current->diMass*(node.x*node.x);
          Izx -= 0.5*current->diMass*(node.x*node.z);
          Ixz -= 0.5*current->diMass*(node.z*node.x);
          Izy -= 0.5*current->diMass*(node.y*node.z);
          Iyz -= 0.5*current->diMass*(node.z*node.y);
        } break;
        case 3: {
          if (current->jdof == -1 || current->jdof == 3) {
            Ixx += current->diMass;
          } else if (current->jdof == 4) {
            Ixy += current->diMass;
            Iyx += current->diMass;
          } else if (current->jdof == 5) {
            Ixz += current->diMass;
            Izx += current->diMass;
          }
        } break;
        case 4: {
          if (current->jdof == -1 || current->jdof == 4) {
            Iyy += current->diMass;
          } else if (current->jdof == 3) {
            Ixy += current->diMass;
            Iyx += current->diMass;
          } else if (current->jdof == 5) {
            Iyz += current->diMass;
            Izy += current->diMass;
          }
        } break;
        case 5: {
          if (current->jdof == -1 || current->jdof == 5) {
            Izz += current->diMass;
          } else if (current->jdof == 3) {
            Ixz += current->diMass;
            Izx += current->diMass;
          } else if (current->jdof == 4) {
            Iyz += current->diMass;
            Izy += current->diMass;
          }
        } break;
      }
    //}
    current = current->next;
  }

  Mx += totmas;
  My += totmas;
  Mz += totmas;
  totmas = (Mx+My+Mz)/3.0;

  if(!(groupId > 0)) { // when groupId is set, only the total mass is computed and printed to a file
  if(Mx != 0.0) xc /= Mx;
  if(My != 0.0) yc /= My;
  if(Mz != 0.0) zc /= Mz;
  // change moments of inertia to centroidal axes using parallel axes theorem: I_z = I_cm + m*d^2
  Ixx -= (My*(yc*yc)+Mz*(zc*zc));
  Iyy -= (Mx*(xc*xc)+Mz*(zc*zc));
  Izz -= (Mx*(xc*xc)+My*(yc*yc));
  Ixy += Mx*xc*yc;
  Iyx += My*yc*xc;
  Ixz += Mx*xc*zc;
  Izx += Mz*zc*xc;
  Iyz += My*yc*zc;
  Izy += Mz*zc*yc;

  if(printFlag) {
    if(Mx != My || Mx != Mz || My != Mz) {
      filePrint(stderr," Directional Mass\n");
      filePrint(stderr," Mx = %e My = %e Mz = %e\n",Mx,My,Mz);
      filePrint(stderr," --------------------------------------\n");
    }

    filePrint(stderr," Moments of Inertia\n");
    filePrint(stderr," Ixx = %e Iyy = %e Izz = %e\n",Ixx,Iyy,Izz);
    filePrint(stderr," --------------------------------------\n");

    filePrint(stderr," Products of Inertia\n");
    filePrint(stderr," Ixy = %e Iyz = %e Ixz = %e\n",Ixy,Iyz,Ixz);
    filePrint(stderr," --------------------------------------\n");

    double tol = 100*std::numeric_limits<double>::epsilon();
    if(std::abs(Iyx-Ixy) > tol || std::abs(Ixz-Izx) > tol || std::abs(Iyz-Izy) > tol)  {
      filePrint(stderr," Warning: Non-Symmetric Products of Inertia\n");
      filePrint(stderr," Iyx = %e Izy = %e Izx = %e\n",Iyx,Izy,Izx);
      filePrint(stderr," --------------------------------------\n");
    }

    filePrint(stderr," Center of Gravity\n");
    filePrint(stderr," x = %e y = %e z = %e\n",xc,yc,zc);
    filePrint(stderr," --------------------------------------\n");
  }

  // Computation to find node closest to the center of gravity
  int nodeMarker = 0;
  double minDistance = 0.0;
  for(i=0; i<numnodes; ++i) {
    Node *node = nodes[i];
    if(node) {
      double x = node->x;
      double y = node->y;
      double z = node->z;

      double dx = xc - x;
      double dy = yc - y;
      double dz = zc - z;

      double distance = sqrt(dx*dx+dy*dy+dz*dz);

      if(i == 0) minDistance = distance;

      if(distance < minDistance) {
        minDistance = distance;
        nodeMarker  = i;
      }
    }
  }

  if(printFlag) {
    filePrint(stderr," Node %d is closest to the center of gravity.\n",nodeMarker+1);
    Node *thisNode = nodes[nodeMarker];
    if(thisNode) {
      filePrint(stderr," Node %d has coordinates: %e %e %e \n",
                nodeMarker + 1, thisNode->x, thisNode->y, thisNode->z);
      filePrint(stderr," It is %e from the center of gravity.\n",minDistance);
      filePrint(stderr," --------------------------------------\n");
    }
  }

  // Compute Geometric center of gravity of structure
  double xmax = -std::numeric_limits<double>::max();
  double ymax = -std::numeric_limits<double>::max();
  double zmax = -std::numeric_limits<double>::max();
  double xmin = std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max();

  int nComponents = renumb.numComp + renumbFluid.numComp;

  if(printFlag) filePrint(stderr," Number of Components = %d\n",renumb.numComp);
  for(int n=0; n<renumb.numComp; ++n) {
    xc = 0.0, yc = 0.0, zc = 0.0;
    int realNodeCnt = 0;
    for(i = renumb.xcomp[n]; i<renumb.xcomp[n+1]; ++i) {
      int inode = renumb.order[i];

      if(dsa->firstdof(inode) == -1 || nodes[inode] == 0) continue;
      realNodeCnt++;
      Node &nd = nodes.getNode(inode);
      xc += nd.x;
      yc += nd.y;
      zc += nd.z;
      if(nd.x > xmax) xmax = nd.x;
      if(nd.y > ymax) ymax = nd.y;
      if(nd.z > zmax) zmax = nd.z;
      if(nd.x < xmin) xmin = nd.x;
      if(nd.y < ymin) ymin = nd.y;
      if(nd.z < zmin) zmin = nd.z;
    }

    double xg = xc/realNodeCnt;
    double yg = yc/realNodeCnt;
    double zg = zc/realNodeCnt;

    if(printFlag) {
      filePrint(stderr," Component %d: Centroid\n", n+1);
      filePrint(stderr," x = %e y = %e z = %e\n",xg,yg,zg);
    }
  }
  for(int n=0; n<renumbFluid.numComp; ++n) {
    xc = 0.0, yc = 0.0, zc = 0.0;
    int realNodeCnt = 0;
    for(i = renumbFluid.xcomp[n]; i<renumbFluid.xcomp[n+1]; ++i) {
      int inode = renumbFluid.order[i];
  
      if(dsa->firstdof(inode) == -1 || nodes[inode] == 0) continue;
      realNodeCnt++;
      Node &nd = nodes.getNode(inode);
      xc += nd.x;
      yc += nd.y;
      zc += nd.z;
      if(nd.x > xmax) xmax = nd.x;
      if(nd.y > ymax) ymax = nd.y;
      if(nd.z > zmax) zmax = nd.z;
      if(nd.x < xmin) xmin = nd.x;
      if(nd.y < ymin) ymin = nd.y;
      if(nd.z < zmin) zmin = nd.z;
    }
  
    double xg = xc/realNodeCnt;
    double yg = yc/realNodeCnt;
    double zg = zc/realNodeCnt;
    
    if(printFlag) {
      filePrint(stderr," Component %d (Fluid): Centroid\n", renumb.numComp+n+1);
      filePrint(stderr," x = %e y = %e z = %e\n",xg,yg,zg);
    }
  }
 
  if(printFlag) {
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," Minimum and maximum x dimension = %e and %e\n",xmin,xmax);
    filePrint(stderr," Minimum and maximum y dimension = %e and %e\n",ymin,ymax);
    filePrint(stderr," Minimum and maximum z dimension = %e and %e\n",zmin,zmax);
    filePrint(stderr," --------------------------------------\n");
  }}

  if(!sinfo.massFile.empty() && (!com || com->cpuNum() == 0)) {
    std::stringstream s;
    s << sinfo.massFile;
    if(groupId > 0) s << '.' << groupId;
    std::ofstream fout(s.str().c_str());
    if(fout.is_open()) {
      fout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << totmas << std::endl;
      fout.close();
    }
  }
  if(!sinfo.massFile.empty() && !(groupId > 0)) {
    std::map<int, std::set<int> > &nodeGroups = geoSource->getNodeGroups();
    for(std::map<int, std::set<int> >::iterator it = nodeGroups.begin(); it != nodeGroups.end(); ++it) {
      if(it->first > 0) {
        double groupMass = computeStructureMass(false, it->first);
        filePrint(stderr, " Mass of group %d = %f\n", it->first, groupMass);
      }
    }
    if(nodeGroups.size() > 0) filePrint(stderr," --------------------------------------\n");
  }

  return totmas;
}

double
Domain::computeFluidMass()
{
  // Compute total mass and mass center of gravity
  // Added calculation for moments of inertia

  double totmas = 0.0; // total mass of fluid
  double xc     = 0.0;
  double yc     = 0.0;
  double zc     = 0.0;
  double Mx     = 0.0;
  double My     = 0.0;
  double Mz     = 0.0;
  double Ixx    = 0.0;
  double Iyy    = 0.0;
  double Izz    = 0.0;
  double Ixy    = 0.0;
  double Iyz    = 0.0;
  double Ixz    = 0.0;
  maxNumNodes = 0;

  // compute max number nodes per element
  int iele, i;
  for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
    int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
    maxNumNodes = std::max(maxNumNodes, numNodesPerElement);
  }

  // allocate one array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
    double elementMass = (*(geoSource->getPackedEsetFluid()))[iele]->getMass(nodes);
    totmas += elementMass;
    int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
    (*(geoSource->getPackedEsetFluid()))[iele]->nodes(nodeNumbers);
    int numRealNodes = 0;
    for(i=0; i<numNodesPerElement; ++i)
      if(nodes[nodeNumbers[i]]) numRealNodes++;

    double massPerNode = elementMass/numRealNodes;

    for(i=0; i<numNodesPerElement; ++i) {
       if(nodes[nodeNumbers[i]] == 0) continue;
       Node &node = nodes.getNode(nodeNumbers[i]);

       double x = node.x;
       double y = node.y;
       double z = node.z;

       xc  += massPerNode*x;
       yc  += massPerNode*y;
       zc  += massPerNode*z;

       Ixx += massPerNode*(y*y + z*z);
       Iyy += massPerNode*(x*x + z*z);
       Izz += massPerNode*(x*x + y*y);

       Ixy -= massPerNode*x*y;
       Ixz -= massPerNode*x*z;
       Iyz -= massPerNode*y*z;
    }
  }

  delete [] nodeNumbers;

  if(totmas != 0.0) xc /= totmas;
  if(totmas != 0.0) yc /= totmas;
  if(totmas != 0.0) zc /= totmas;
  // change moments of inertia to centroidal axes using parallel axes theorem: I_z = I_cm + m*d^2
  Ixx -= (totmas*(yc*yc)+totmas*(zc*zc));
  Iyy -= (totmas*(xc*xc)+totmas*(zc*zc));
  Izz -= (totmas*(xc*xc)+totmas*(yc*yc));

  Ixy += totmas*xc*yc;
  Ixz += totmas*xc*zc;
  Iyz += totmas*yc*zc;

  filePrint(stderr," Moments of Inertia\n");
  filePrint(stderr," Ixx = %e Iyy = %e Izz = %e\n",Ixx,Iyy,Izz);
  filePrint(stderr," --------------------------------------\n");

  filePrint(stderr," Products of Inertia\n");
  filePrint(stderr," Ixy = %e Iyz = %e Ixz = %e\n",Ixy,Iyz,Ixz);
  filePrint(stderr," --------------------------------------\n");

  filePrint(stderr," Center of Gravity\n");
  filePrint(stderr," x = %f y = %f z = %f\n",xc,yc,zc);
  filePrint(stderr," --------------------------------------\n");

  // Computation to find node closest to the center of gravity
  int nodeMarker = 0;
  double minDistance = 0.0;
  for(i=0; i<numnodes; ++i) {
    Node *node = nodes[i];
    if(node) {
      double x = node->x;
      double y = node->y;
      double z = node->z;

      double dx = xc - x;
      double dy = yc - y;
      double dz = zc - z;

      double distance = sqrt(dx*dx+dy*dy+dz*dz);

      if(i == 0) minDistance = distance;

      if(distance < minDistance) {
        minDistance = distance;
        nodeMarker  = i;
      }
    }
  }

  filePrint(stderr," Node %d is closest to the Center of Gravity\n",nodeMarker+1);
  Node *thisNode = nodes[nodeMarker];
  if(thisNode) {
    filePrint(stderr," Node %d has coordinates: %e %e %e \n",
              nodeMarker + 1, thisNode->x, thisNode->y, thisNode->z);
    filePrint(stderr," It is %e from the center of gravity\n",minDistance);
    filePrint(stderr," --------------------------------------\n");
  }

  // Compute Geometric center of gravity of fluid
  double xmax = 0.0;
  double ymax = 0.0;
  double zmax = 0.0;

  int nComponents = renumbFluid.numComp;

  filePrint(stderr," Number of Components = %d\n",renumbFluid.numComp);
  for(int n=0; n<nComponents; ++n) {
    xc = 0.0, yc = 0.0, zc = 0.0;
    int realNodeCnt = 0;
    for(i = renumbFluid.xcomp[n]; i<renumb.xcomp[n+1]; ++i) {
      int inode = renumbFluid.order[i];

      if(dsa->firstdof(inode) == -1 || nodes[inode] == 0) continue;
      realNodeCnt++;
      Node &nd = nodes.getNode(inode);
      xc += nd.x;
      yc += nd.y;
      zc += nd.z;
      if(nd.x > xmax) xmax = nd.x;
      if(nd.y > ymax) ymax = nd.y;
      if(nd.z > zmax) zmax = nd.z;
    }

    double xg = xc/realNodeCnt;
    double yg = yc/realNodeCnt;
    double zg = zc/realNodeCnt;

    filePrint(stderr," Component (%d,%d): Centroid\n",
              renumbFluid.xcomp[n], renumbFluid.xcomp[n+1]);
    filePrint(stderr," x = %f y = %f z = %f\n",xg,yg,zg);
  }
  filePrint(stderr," --------------------------------------\n");
  filePrint(stderr," Maximum x dimension = %f\n",xmax);
  filePrint(stderr," Maximum y dimension = %f\n",ymax);
  filePrint(stderr," Maximum z dimension = %f\n",zmax);
  filePrint(stderr," --------------------------------------\n");

  return totmas;
}

Rbm *
Domain::constructHzem(bool printFlag)
{
  Rbm *rbm = new Rbm(dsa, c_dsa);
  if(printFlag)
    std::cerr << " ... HZEM algorithm detected " << rbm->numRBM() << " zero energy modes ...\n";
  return rbm;
}

// Function to construct zero energy modes
Rbm *
Domain::constructSlzem(bool printFlag)
{
  Rbm *rbm = new Rbm(dsa, c_dsa);
  if(printFlag)
    std::cerr << " ... SZEM algorithm detected " << rbm->numRBM() << " zero energy modes ...\n";
  return rbm;
}

Rbm *
Domain::constructRbm(bool printFlag)
{
  if(printFlag) filePrint(stderr," ... Using Geometric RBM Method     ...\n");
  Rbm *rbm = 0;
  if(numLMPC && sinfo.grbm_use_lmpc) {
    if(renumb_nompc.numComp == 0)
      rbm = new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd, renumb, numLMPC, lmpc);
    else
      rbm = new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd, renumb_nompc, numLMPC, lmpc);
    // delete the LMPCs
    for(int i=0; i<numLMPC; ++i) if(lmpc[i]) delete lmpc[i];
    lmpc.deleteArray();
    lmpc.restartArray();
    numLMPC = 0;
  }
  else 
    rbm = new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd, renumb);
  if(printFlag)
    std::cerr << " ... GRBM algorithm detected " << rbm->numRBM() << " rigid body modes ...\n";
  return rbm;
}

const char* problemTypeMessage[] = {
" ... Linear Static Analysis         ... \n",
" ... Linear Dynamics/Quasistatics   ... \n",
" ... Modal Analysis                 ... \n",
" ... Nonlinear Static Analysis      ... \n",
" ... Nonlinear Dynamics/Quasistatics... \n",
" ... Nonlin. Stat. Anal. Arc Length ... \n",
" ... Condition Number Analysis      ... \n",
" ... Thermal Analysis               ... \n",
" ... Outputing TOPDOMDEC File       ... \n",
" ... Axisymmetric Acoustic Analysis ... \n",
" ... Material Nonlin. Static Analysis ... \n",
" ... Material Nonlin. Dynamic Analysis ... \n",
" ... Acoustic Scattering Analysis   ... \n",
" ... Acoustic Frequency Sweep Analysis   ... \n",
" ... Helmholtz Direction Sweep Analysis ... \n",
" ... HelmholtzMF Analysis           ... \n",
" ... HelmholtzSO Analysis           ... \n",
" ... Computing Decomposition        ... \n",
" ... Nonlin. Themal Dynamic Analysis... \n",
" ... Discontinuous Enrichment Method... \n",
" ... POD ROM Offline Computations   ... \n",
""
};

const char* sensitivityTypeMessage[] = {
" ... Sensitivity Analysis without Static Analysis ...\n",
" ... Sensitivity Analysis with Static Analysis ...\n"
};

const char* solverTypeMessage[] = {
" ... Skyline Solver is Selected     ... \n",
" ... Sparse Solver is Selected      ... \n",
" ... BlockSky Solver is Selected    ... \n",
#ifdef USE_EIGEN3
" ... SimplicalLLT Solver is Selected... \n",
" ... SimplicalLDLT Solver is Selec'd... \n",
#else
" ... Skyline Solver is Selected     ... \n",
" ... Skyline Solver is Selected     ... \n",
#endif
#ifdef EIGEN_CHOLMOD_SUPPORT
" ... Cholmod Solver is Selected     ... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
#ifdef EIGEN_UMFPACK_SUPPORT
" ... UmfPack Solver is Selected     ... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
" ... SuperLU Solver is Selected     ... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
#ifdef USE_SPOOLES
" ... Spooles Solver is Selected     ... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
#ifdef USE_MUMPS
" ... Mumps Solver is Selected       ... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
"",
" ... POD-GN Solver is Selected      ... \n",
" ... POD-Galerkin Solver is Selected... \n",
" ... POD-Galerkin Solver is Selected... \n",
" ... Goldfarb Solver is Selected    ... \n",
#ifdef EIGEN_SPARSELU_SUPPORT
" ... SparseLU Solver is Selected    ... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
#ifdef EIGEN_SPQR_SUPPORT
" ... SuiteSparseQR Solver is Selec'd... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
#ifdef USE_EIGEN3
" ... SparseQR Solver is Selected    ... \n",
#else
" ... Skyline Solver is Selected     ... \n",
#endif
" ... Cholmod Solver is Selected     ... \n",
" ... Pardiso Solver is Selected     ... \n"
};

void
Domain::preProcessing()
{
 // ... CONSTRUCT DOMAIN ELEMENT TO NODE CONNECTIVITY
 matrixTimers->makeConnectivity -= getTime();
 if(elemToNode) delete elemToNode;
 elemToNode = new Connectivity(packedEset.asSet());
 if(sinfo.HEV) {
   if(elemToNodeFluid) delete elemToNodeFluid;
   elemToNodeFluid = new Connectivity(geoSource->getPackedEsetFluid()->asSet());
 }
 matrixTimers->makeConnectivity += getTime();

 // ... RENUMBER IF NECESSARY AND/OR ASKED FOR
 matrixTimers->renumbering -= getTime();
 Renumber rnum = getRenumbering();
 Renumber* rnumFluid = 0;
 if(sinfo.HEV) rnumFluid = getRenumberingFluid();
 matrixTimers->renumbering += getTime();

 // ... CONSTRUCT DOF SET ARRAY
 matrixTimers->createDofs -= getTime();
 if(dsa) delete dsa;
 dsa = new DofSetArray(numnodes, packedEset, rnum.renumb);
 if(sinfo.HEV) {
   dsaFluid = new DofSetArray(numnodesFluid, *(geoSource->getPackedEsetFluid()), rnumFluid->renumb);
 }
 matrixTimers->createDofs += getTime();
}

void
Domain::make_constrainedDSA()
{
 if(c_dsa) delete c_dsa;
 c_dsa = new ConstrainedDSA(*dsa, numDirichlet, dbc, numComplexDirichlet, cdbc);
 if (solInfo().HEV)  {
   c_dsaFluid = new ConstrainedDSA(*dsaFluid, numDirichletFluid, dbcFluid);
 }
}

void
Domain::make_constrainedDSA(const int *bc)
{
 // c_dsa = constrained dof set array
 // dsa   = dof set array
 // dbc   = dirichlet boundary conditions

 // numDirichlet = number of dirichlet boundary conditions
 // cdbc         = complex dirichlet boundary conditions

 // numComplexDirichlet = number of complex
 //                       dirichlet boundary conditions

 // bc = integer array marking dofs that are
 //      constrained or have forces applied

 if(c_dsa) delete c_dsa;
 c_dsa = new ConstrainedDSA(*dsa, dbc, numDirichlet, cdbc,
                            numComplexDirichlet, bc);
}

void
Domain::make_constrainedDSA(int fake)
{ // make a fake constrainedDSA if fake !=0; ie lie to the code by telling
  //   it that there are noconstraints
  if(fake) {
    if(c_dsa) delete c_dsa;
    c_dsa = new ConstrainedDSA(*dsa, 0, dbc);
  }
  else {
    make_constrainedDSA();
  }
}

static const char* topMes[] = {
" ... Writing TOPDOMDEC File         ... \n",
" ... Writing renumbered TOP File    ... \n",
" ... Writing material TOP File      ... \n",
"",
"",
" ... Writing axisymetric TOP File with 2D planes ... \n",
" ... Writing axisymetric 3D TOP File ... \n",
" ... Writing renumbered material TOP File      ... \n"
};

void Domain::writeTopFileElementSets(ControlInfo *cinfo, int * nodeTable, int* nodeNumber, int topFlag)
{
  // ... WRITE ELEMENT CONNECTIVITY
 if(packedEset.last() > 0)  // PJSA: replaced numele with packedEset.last() to include phantoms
   fprintf(cinfo->checkfileptr,"Elements %s using %s\n",
           cinfo->elemSetName, cinfo->nodeSetName);

 int inode, iele;
 int nEls = packedEset.last();

 std::list<int> phantoms;
 std::list<int> constraints;
 std::list<int> dimasses;

 for(iele=0; iele<nEls; ++iele) {
   if(!packedEset[iele]->isPhantomElement() && !packedEset[iele]->isConstraintElement()) {
     int numNodesPerElement = packedEset[iele]->numTopNodes();
     int eletype = packedEset[iele]->getTopNumber();
     if(eletype == 506) { dimasses.push_back(iele); continue; }
     if(numNodesPerElement <= 1) continue;
     int eid = (topFlag == 1 || topFlag == 7) ? iele+1 : packedEset[iele]->getGlNum()+1; // only renumber for -T and -M
     fprintf(cinfo->checkfileptr,"%6d  %4d ",eid,eletype);
     packedEset[iele]->nodes(nodeNumber);
     for(inode=0; inode<numNodesPerElement; ++inode)
       // Avoid to print nodes that are internally created
       if(nodeNumber[inode] < numnodes && nodes[nodeNumber[inode]] != 0)
	 fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber[inode]]);
     fprintf(cinfo->checkfileptr,"\n");
   }
   else {
     if(packedEset[iele]->isPhantomElement())
       phantoms.push_back(iele); //TG create a list of phantom elements
     else
       constraints.push_back(iele);
   }
 }

 //TG output dimasses in a separate element set if there are any
 // as of now (7/5/06) a dimass will be represented by a bar element in xpost between the node the
 // dimass is attached to and itself.
 if(nDimass > 0 || dimasses.size() > 0) {
   fprintf(stderr, " ... Putting %d Dimasses as bars in a separate ElemSet\n", nDimass+int(dimasses.size()));
   fprintf(cinfo->checkfileptr,"Elements %s_dimass using %s\n",
           cinfo->elemSetName, cinfo->nodeSetName);
   int iele = 0;
   int eletype = 506; // TG now 506 to detect dimasses in Xpost
   DMassData* curMass = firstDiMass;
   while(curMass != NULL) {
     int node = curMass->node;
     if(node < numnodes && nodes[node] != 0) {
       fprintf(cinfo->checkfileptr,"%6d  %4d ",iele+1,eletype);

       fprintf(cinfo->checkfileptr,"%6d %6d\n",nodeTable[node], nodeTable[node]);
       curMass = curMass->next;
       iele++;
     }
     else {
       std::cout << " Warning : virtual dimass" << std::endl;
       curMass = curMass->next;
     }
   }

   for(int i=0; i<dimasses.size(); ++i) {
     int jele = dimasses.front();
     dimasses.pop_front();

     // same as main element writing in loop above for non dimass elements
     int eletype = packedEset[jele]->getTopNumber();
     fprintf(cinfo->checkfileptr,"%6d  %4d ",iele+1,eletype);
     packedEset[jele]->nodes(nodeNumber);
     for(inode=0; inode<2; ++inode)
       fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber[0]]);
     fprintf(cinfo->checkfileptr,"\n");
     iele++;
   }
 }

 // output phantom elements in a separate elementset if there are any
 if (phantoms.size() > 0) {
     fprintf(cinfo->checkfileptr,"Elements %s_phantom using %s\n",
	     cinfo->elemSetName, cinfo->nodeSetName);
     int m_phantoms = phantoms.size();
     for(int i=0; i<m_phantoms; ++i){
       iele = phantoms.front();
       phantoms.pop_front();

       // same as main element writing in loop above for non phantom elements
       int numNodesPerElement = packedEset[iele]->numTopNodes();
       if(numNodesPerElement <= 1) continue;
       int eletype = packedEset[iele]->getTopNumber();
       fprintf(cinfo->checkfileptr,"%6d  %4d ",iele+1,eletype);
       packedEset[iele]->nodes(nodeNumber);
       for(inode=0; inode<numNodesPerElement; ++inode)
	 // Avoid printing nodes that are internally created
	 if(nodeNumber[inode] < numnodes && nodes[nodeNumber[inode]] != 0)
	   fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber[inode]]);
       fprintf(cinfo->checkfileptr,"\n");
     }
 }

 // output constraint elements in a separate elementset if there are any
 if (constraints.size() > 0) {
     fprintf(cinfo->checkfileptr,"Elements %s_constraints using %s\n",
	     cinfo->elemSetName, cinfo->nodeSetName);
     int m_constraints = constraints.size();
     for(int i=0; i<m_constraints; ++i){
       iele = constraints.front();
       constraints.pop_front();

       // same as main element writing in loop above for non constraints elements
       int numNodesPerElement = packedEset[iele]->numTopNodes();
       if(numNodesPerElement <= 1) continue;
       int eletype = packedEset[iele]->getTopNumber();
       fprintf(cinfo->checkfileptr,"%6d  %4d ",iele+1,eletype);
       packedEset[iele]->nodes(nodeNumber);
       for(inode=0; inode<numNodesPerElement; ++inode)
	 // Avoid printing nodes that are internally created
	 if(nodeNumber[inode] < numnodes && nodes[nodeNumber[inode]] != 0)
	   fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber[inode]]);
       fprintf(cinfo->checkfileptr,"\n");
     }
   }

 // output surface elements, each in a separate element set if there are any
 // for now, just output the vertices
 for(int iSurf=0; iSurf<nSurfEntity; iSurf++) {
   if(SurfEntities[iSurf]->GetId() == 0) continue;
   fprintf(cinfo->checkfileptr,"Elements surface_%d using %s\n",
           SurfEntities[iSurf]->GetId(), cinfo->nodeSetName);
   FaceElemSet &faceElemSet = SurfEntities[iSurf]->GetFaceElemSet();
   int nFaceEls = faceElemSet.last();
   for(iele=0; iele<nFaceEls; ++iele) {
     int nVertices = faceElemSet[iele]->nVertices();
     int eletype;
     if(nVertices == 3) eletype = 104;
     else if(nVertices == 4) eletype = 2;
     else { std::cerr << "don't know xpost eletype for surface " << SurfEntities[iSurf]->GetId() << " element " << iele << " nVertices = " << nVertices << std::endl; continue; }
     int eleID = iele+1;
     fprintf(cinfo->checkfileptr,"%6d  %4d ",eleID,eletype);
     for(inode=0; inode<nVertices; ++inode) {
       int nodeNumber = faceElemSet[iele]->GetVertex(inode);
       fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber]);
     }
     fprintf(cinfo->checkfileptr,"\n");
   }
 }

 // PJSA 3-24-05 write default pattern for ElemScalar output
  if(nEls > 0)
   fprintf(cinfo->checkfileptr,"Pattern default using %s\n", cinfo->elemSetName);

}

void
Domain::makeNodeTable(int topFlag)
{
 long m4 = - memoryUsed();
 nodeTable = new int[nodes.size()];
 m4 += memoryUsed();
 //fprintf(stderr," ... Node Table %14.3f Mb\n",m4/(1024.0*1024.0));

 int inode;
 for(inode=0; inode<nodes.size(); ++inode)
   nodeTable[inode] = -1;

 // if topFlag = 0, output a topdomdec file with all nodes
 // including nodes that are not defined or used.

 // if topFlag = 1, output a topdomdec file without gaps
 // (nodes renumbered sequentially)

 // if topFlag = 2, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and with all nodes including nodes that are not defined or used.

 // if topFlag = 7, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and without gaps (nodes renumbered sequentially)

 exactNumNodes = 0;
 for(inode=0; inode<nodes.size(); ++inode) {
   if(nodes[inode] == 0) continue;
   exactNumNodes += 1;
   nodeTable[inode] = ((topFlag == 1) || (topFlag == 7)) ? exactNumNodes : inode+1;
 }
}

void
Domain::makeTopFile(int topFlag)
{
 if (solInfo().HEV)  {
   solInfo().HEV = 0;
   packedEset.deleteElems();
   numele = geoSource->getElems(packedEset);
 }

 fprintf(stderr," ... Memory Used so far    %14.3f Mb\n",
                 memoryUsed()/(1024.0*1024.0));

 if(topFlag >= 0)
   fprintf(stderr,"%s",topMes[topFlag]);

 // ... CONSTRUCT DOMAIN ELEMENT TO NODE CONNECTIVITY
 long m1 = - memoryUsed();
 elemToNode = new Connectivity( packedEset.asSet() );
 m1 += memoryUsed();

 // ... CONSTRUCT DOF SET ARRAY
 long m2 = - memoryUsed();
 dsa = new DofSetArray( numnodes, packedEset );
 m2 += memoryUsed();

 // Renumber rnum = getRenumbering();
 // nodeToNode->findMaxDist(rnum.renumb);
 // nodeToNode->findProfileSize(dsa);

 fprintf(stderr," ... Elem. to Node Connectivity %9.3f Mb\n",
         m1/(1024.0*1024.0));
 fprintf(stderr," ... DOF set array              %9.3f Mb\n",
         m2/(1024.0*1024.0));

 // ... WRITE INPUT FILE FOR TOP/DOMDEC
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 cinfo->checkfileptr = openFile(cinfo->checkfile, ".top");

 // Allocate memory for a node table to compact node numbers
 // into sequential format
 //fprintf(stderr," ... Total Memory         = %13.3f Mb\n",
 //               memoryUsed()/(1024.0*1024.0));
 long m4 = - memoryUsed();

 int *nodeTable = new int[numnodes];
 m4 += memoryUsed();
 //fprintf(stderr," ... Node Table %14.3f Mb\n",m4/(1024.0*1024.0));

 int inode;
 for(inode=0; inode<numnodes; ++inode)
   nodeTable[inode] = -1;

 // if topFlag = 0, output a topdomdec file with all nodes
 // including nodes that are not defined or used.

 // if topFlag = 1, output a topdomdec file without gaps
 // (nodes renumbered sequentially)

 // if topFlag = 2, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and with all nodes including nodes that are not defined or used.

 // if topFlag = 7, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and without gaps (nodes renumbered sequentially)

 // ... WRITE NODE COORDINATES
 fprintf(cinfo->checkfileptr,"Nodes %s\n",cinfo->nodeSetName);
 int exactNumNodes = 0;
 for(inode=0; inode<nodes.size(); ++inode) {
   if(nodes[inode] == 0) continue;
   Node &nd = nodes.getNode(inode);
   exactNumNodes += 1;
   nodeTable[inode] = ((topFlag == 1) || (topFlag == 7)) ? exactNumNodes : inode+1;
   double x = nd.x;
   double y = nd.y;
   double z = nd.z;
   fprintf(cinfo->checkfileptr,"%d\t % 14.6f\t% 14.6f\t % 14.6f\n",
                              nodeTable[inode],x,y,z);
 }

 // First find the maximum number of nodes per Element
 // this used to be done in the write element connectivity loop
 maxNumNodes=0;
 int iele;
 int nEls = packedEset.last(); //HB: to avoid calling packedEset.last() at each loop
 for(iele=0; iele<nEls; ++iele) {
   int numNodesPerElement = packedEset[iele]->numNodes();
   maxNumNodes = std::max(numNodesPerElement, maxNumNodes);
 }
 // allocate integer array to store node numbers
 int *nodeNumber = new int[maxNumNodes];

 // ... WRITE ELEMENT CONNECTIVITY
 // TG moved to this function as will need to be outputted in .top file for -m and -M as well.
 writeTopFileElementSets(cinfo, nodeTable, nodeNumber, topFlag);
 //
 // ElemSet creation happens in there
 // ElemSet_phantom, ElemSet_dimass and ElementSet_constraints aswell

 // write structure displacement boundary conditions (dirichlet)
 // if these type of boundary conditions exist in the input file.
 int i;
 if(numDirichlet + numComplexDirichlet > 0)
   fprintf(cinfo->checkfileptr,"SDBoundary %s using %s\n",
           cinfo->bcondSetName,cinfo->nodeSetName);
 for(i=0; i<numDirichlet; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[dbc[i].nnum],dbc[i].dofnum+1,dbc[i].val);
 }

/*
 // ... ADD CONVECTIVE FLUXES
 // KHP: get the correct header for the TOPDOMDEC file!
 if(numConvBC > 0)
   fprintf(cinfo->checkfileptr,"SDTemperature %s using %s\n",
           cinfo->bcondSetName,cinfo->nodeSetName);
 for(i=0; i<numConvBC; ++i) {
  fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[cvbc[i].nnum], cvbc[i].dofnum+1,cvbc[i].val);
  }
*/


 // also print complex dirichlet boundary conditions, for Helmholtz problem
 for(i=0; i<numComplexDirichlet; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[cdbc[i].nnum],cdbc[i].dofnum+1,cdbc[i].reval);
 }

 // write structure force boundary conditions (Neumann)
 // if these type of boundary conditions exist in the input file.
 if(numNeuman + numComplexNeuman > 0)
   fprintf(cinfo->checkfileptr,"SFBoundary %s using %s\n",
           cinfo->bcondSetName,cinfo->nodeSetName);

 for(i=0; i<numNeuman; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[nbc[i].nnum],nbc[i].dofnum+1,nbc[i].val);
 }

 // also write complex Neumann boundary conditions, for Helmholtz problem
 // note we are only outputing the real value of the complex bc as topdomdec
 // will not be able to understand a complex entry (i.e. 2 double values
 // instead of 1)
 for(i=0; i<numComplexNeuman; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[cnbc[i].nnum],cnbc[i].dofnum+1,cnbc[i].reval);
 }

 // If we have more than one component, output a decomposition file
 // specifing the elements in each component. Very useful for debugging
 // a model.

 int *hasAppeared = new int[packedEset.last()];
 for(i=0; i<nEls; ++i)
   hasAppeared[i] = 0;

 // number of components in the structure input file
 //int nComponents = renumb.numComp;
 int nComponents = 1;

 // always output component file if there is more than one component.

 if(nComponents > 1) {

   FILE *componentFile = fopen("component.dec","w");

   int count = 0;

   fprintf(componentFile,"Decomposition COMPONENTS for %s\n%d\n",
           cinfo->elemSetName,nComponents);
   int n;
   for(n = 0; n<nComponents; ++n) {
     for(i = renumb.xcomp[n]; i<renumb.xcomp[n+1]; ++i) {
       inode = renumb.order[i];
       int numElemPerNode = nodeToElem->num(inode);
       int j;
       for(j=0; j<numElemPerNode; ++j) {
         if(hasAppeared[(*nodeToElem)[inode][j]] == 0) {
           hasAppeared[(*nodeToElem)[inode][j]] = (*nodeToElem)[inode][j]+1;
           count++;
         }
       }
     }
     fprintf(componentFile,"%d\n",count);
     for(i=0; i<nEls; ++i) {
       if(hasAppeared[i] != 0) {
         fprintf(componentFile," %d\n",hasAppeared[i]);
         hasAppeared[i] = 0;
       }
     }
     count = 0;
   }
   fflush(componentFile);
 }

 delete [] hasAppeared;

// FOR DANIEL - to output element lists for each attribute

 if((topFlag == 2) || (topFlag == 7)) {
   FILE *matList = fopen("material.top","w");

   // ... WRITE NODE COORDINATES
   fprintf(matList,"Nodes %s\n", cinfo->nodeSetName);
   for(inode=0; inode<numnodes; ++inode) {
     if(nodes[inode] == 0) continue;
     Node &nd = nodes.getNode(inode);
     // if(topFlag && nodes[inode] != 0) continue;
     double x = 0.0;
     double y = 0.0;
     double z = 0.0;
     if(nodes[inode] != 0) {
       x = nd.x;
       y = nd.y;
       z = nd.z;
     }
     fprintf(matList,"%d\t % 14.6f\t% 14.6f\t % 14.6f\n",
             nodeTable[inode],x,y,z);
   }

   int n;
   std::map<int, Attrib> &attrib = geoSource->getAttributes();
   int na = geoSource->getNumAttributes(); // this is actually the number of elements !!
   SPropContainer &sProps = geoSource->getStructProps();
   for(std::map<int, StructProp>::iterator it = sProps.begin(); it != sProps.end(); ++it) {
       n = it->first;
       bool first = true; // PJSA to deal with case of empty EleSet
       for(i=0; i<na; ++i) {
         if(attrib[i].attr == n) {
           int e = attrib[i].nele; 
           Element *elem = geoSource->getElem(e);
           if(elem) {
             if(first) {
               fprintf(matList,"Elements EleSet%d using %s\n", n+1, cinfo->nodeSetName);
               first = false;
             }
             int eletype = elem->getTopNumber();
             fprintf(matList,"%d %d ",e+1,eletype);
             int numNodesPerElement = elem->numTopNodes();
             if(numNodesPerElement <= 1) continue;
             for(inode=0; inode<numNodesPerElement; ++inode) {
               elem->nodes(nodeNumber);
               fprintf(matList,"%d ",nodeTable[nodeNumber[inode]]);
             }
             fprintf(matList,"\n");
           }
         }
       }
       if (output_match_in_top)// Match case for ElemScalar
	 fprintf(matList,"Match EleSet%d with default\n", n+1);
   }

   // PJSA 9-2-03 write boundary conditions to material.top
   if(numDirichlet + numComplexDirichlet > 0)
     fprintf(matList,"SDBoundary %s using %s\n",
             cinfo->bcondSetName, cinfo->nodeSetName);
   for(i=0; i<numDirichlet; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[dbc[i].nnum],dbc[i].dofnum+1,dbc[i].val);
   }

   for(i=0; i<numComplexDirichlet; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[cdbc[i].nnum],cdbc[i].dofnum+1,cdbc[i].reval);
   }
   if(numNeuman + numComplexNeuman > 0)
     fprintf(matList,"SFBoundary %s using %s\n",
             cinfo->bcondSetName,cinfo->nodeSetName);
   for(i=0; i<numNeuman; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[nbc[i].nnum],nbc[i].dofnum+1,nbc[i].val);
   }
   for(i=0; i<numComplexNeuman; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[cnbc[i].nnum],cnbc[i].dofnum+1,cnbc[i].reval);
   }

 }

 // ... CALCULATE STRUCTURE MASS IF REQUESTED
 if(sinfo.massFlag)  {
   Renumber rnum = getRenumbering();
   double mass = computeStructureMass();
   fprintf(stderr," --------------------------------------\n");
   fprintf(stderr," ... Structure mass = %e  ...\n",mass);
   fprintf(stderr," --------------------------------------\n");
 }

 // Compute Total memory requested for constructing the TOP/DOMDEC file
 long m3 = memoryUsed();
 fprintf(stderr," ... Total Memory         = %13.3f Mb\n",m3/(1024.0*1024.0));
}

void
Domain::setsizeSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1)  sizeSfemStress = geoSource->numNode();  // node-based output
  else if(avgnum == 0) {  // element-based output
   sizeSfemStress = 0;
   Connectivity *elemToNode = new Connectivity(domain->getEset()->asSet());
   int numele = geoSource->getNumAttributes();  // number of elements; another option domain->numElements();
   for(int iele=0; iele<numele; ++iele)   {
//     std::cerr << "number of nodes in this element  = " << elemToNode->num(iele) << std::endl;
     sizeSfemStress = sizeSfemStress + elemToNode->num(iele); // add number of nodes for each element
   }
  }
  else {
   std::cerr << "avgnum = " << avgnum << " not implemented in Domain::setsizeSfemStress()" << std::endl;
   sizeSfemStress = 0;
  }
}

double*
Domain::getSfemStress(int fileNumber, double* dummystress)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  if(avgnum == 1)  return stress->data();
  else if(avgnum == 0) return stressAllElems->data();
  else {std::cerr << "avgnum = " << avgnum << " not implemented in Domain::getSfemStress()" << std::endl; return 0;}
}

void
Domain::updateSfemStress(double* str, int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  if(avgnum == 1)  for (int i=0;i<stress->size();++i) (*stress)[i] = str[i];
  else if(avgnum == 0) for (int i=0;i<stressAllElems->size();++i) (*stressAllElems) = str[i];
  else {std::cerr << "avgnum = " << avgnum << " not implemented in Domain::updateSfemStress()" << std::endl;}
}

void 
Domain::computeNormalizedVonMisesStress(Vector &sol, double *bcx, int surface, Vector &stress, bool normalized)
{
  Vector weight(numNodes(),0.0);
  stress.zero();    weight.zero();

  double *nodalTemperatures = 0;
  // Either get the nodal temperatures from the input file or
  // from the thermal model
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >= 0) nodalTemperatures = temprcvd;
  int *nodeNumbers = new int[maxNumNodes];

  for(int iele = 0; iele < numele; ++iele) {

    // Don't do anything if element is a phantom or constraint
    if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;

    packedEset[iele]->nodes(nodeNumbers);
    int NodesPerElement = elemToNode->num(iele);
    Vector elDisp(maxNumDOFs,0.0), elweight(NodesPerElement, 0.0), elstress(NodesPerElement,0.0), elemNodeTemps(maxNumNodes);
    elemNodeTemps.zero();    elDisp.zero();    elstress.zero();    elweight.zero();

    // DETERMINE ELEMENT DISPLACEMENT VECTOR
    for (int k = 0; k < allDOFs->num(iele); ++k) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
      if (cn >= 0)
        elDisp[k] = sol[cn];
      else
        elDisp[k] = bcx[(*allDOFs)[iele][k]];
    }

    for (int iNode = 0; iNode < NodesPerElement; ++iNode) {
      if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp) {
        elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
      } else
        elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
    }

    // transform displacements from DOF_FRM to basic coordinates
    transformVectorInv(elDisp, iele);
    packedEset[iele]->getVonMises(elstress, elweight, nodes, elDisp, 6, surface, elemNodeTemps.data());

    // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
    for(int k = 0; k < NodesPerElement; ++k) {
      int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
      stress[node] += elstress[k];
      weight[node] += elweight[k];
    }
  }

  for(int k = 0; k <numNodes(); ++k)  {
    if(weight[k] == 0.0)
      stress[k] = 0.0;
    else {
      if(normalized) stress[k] /= (*aggregatedStress)*weight[k];
      else stress[k] /= weight[k];
    }
  }

}

void
Domain::scaleToTrueVonMisesStress(Vector &stress)
{
  for(int k = 0; k <numNodes(); ++k) {
    stress[k] *= (*aggregatedStress);
  }
}

void
Domain::getStressStrain(Vector &sol, double *bcx, int fileNumber,
                        int stressIndex, double time, int printFlag)
{
  int numNodes = (outFlag) ? exactNumNodes : geoSource->numNode();

  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS, WHERE
  // ... ONLY THE AXIAL STRAIN (EXX) AND AXIAL STRESS (SXX) ARE CALCULATED

  OutputInfo *oinfo = geoSource->getOutputInfo();

  // avgnum = 2 --> do not include stress/strain of bar/beam element in averaging
  // avgnum = 3 --> only include elements whose nodes are in the ouputGroup in averaging
  int avgnum = oinfo[fileNumber].averageFlg;

  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;
  int surface = oinfo[fileNumber].surface;
    // upper  surface = 1
    // median surface = 2
    // lower  surface = 3
  int str_therm_option = oinfo[fileNumber].str_therm_option;
    // thermomechanical = 0
    // thermal = 1
    // mechanical = 2
  OutputInfo::FrameType oframe = oinfo[fileNumber].oframe;

  int k;
  int iele;
  int iNode;

  double *nodalTemperatures = 0;
  // Either get the nodal temperatures from the input file or
  // from the thermal model
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >= 0) nodalTemperatures = temprcvd;

  if(printFlag != 2) {
    // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
    if(avgnum > 0) {
      if(stress == 0) stress = new Vector(numNodes,0.0);
      if(weight == 0) weight = new Vector(numNodes,0.0);
    }
    else if(printFlag == 1 && stressAllElems == 0) stressAllElems = new Vector(sizeSfemStress,0.0);
    if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

    if((elstress == 0) || (elweight == 0) || (p_elstress == 0 && oframe != OutputInfo::Global)) {
      int NodesPerElement, maxNodesPerElement=0;
      for(iele=0; iele<numele; ++iele) {
        NodesPerElement = elemToNode->num(iele);
        maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
      }
      if(elstress == 0) elstress = new Vector(maxNodesPerElement, 0.0);
      if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
      if(p_elstress == 0 && oframe != OutputInfo::Global) p_elstress = new FullM(maxNodesPerElement,9);
    }

    if(avgnum > 0) {
      // zero the vectors
      stress->zero();
      weight->zero();
    }
    else if(printFlag == 1) stressAllElems->zero();
  }

  int count = 0;
  for(iele = 0; iele < numele; ++iele) {

    int NodesPerElement = elemToNode->num(iele);
    packedEset[iele]->nodes(nodeNumbers);

    if(printFlag != 2) {

      // Don't do anything if element is a phantom or constraint
      if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;

      // Don't include beams or bars in the averaging if nodalpartial (avgnum = 2) is requested
      if ((avgnum == 2 && packedEset[iele]->getElementType() == 6) ||
          (avgnum == 2 && packedEset[iele]->getElementType() == 7) ||
          (avgnum == 2 && packedEset[iele]->getElementType() == 1)) continue;

      // Don't include elements with one or more nodes not in the group if nodalpartialgroup (avgnum = 3) is requested
      if (avgnum == 3) {
        int groupId = oinfo[fileNumber].groupNumber;
        if (groupId > 0) {
          std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
          std::set<int>::iterator it;
          for (iNode = 0; iNode < NodesPerElement; ++iNode)
            if((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
          if(it == groupNodes.end()) continue;
        }
      }

      // Only include specified element in case of single element output
      if ((avgnum == -1 || avgnum == 0) && (oinfo[fileNumber].nodeNumber != -1) &&
          (oinfo[fileNumber].nodeNumber != packedEset[iele]->getGlNum())) continue;

      elDisp->zero();
      elstress->zero();
      elweight->zero();

      // DETERMINE ELEMENT DISPLACEMENT VECTOR
      if(str_therm_option != 1) {
        for (k = 0; k < allDOFs->num(iele); ++k) {
          int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
          if (cn >= 0)
            (*elDisp)[k] = sol[cn];
          else
            (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
        }
      }

      for (iNode = 0; iNode < NodesPerElement; ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp || str_therm_option == 2)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }

      // transform displacements from DOF_FRM to basic coordinates
      transformVectorInv(*elDisp, iele);

      // transform non-invariant stresses/strains from basic frame to DOF_FRM or CFRAME
      if(oframe != OutputInfo::Global && ((stressIndex >= 0 && stressIndex <= 5) || (stressIndex >= 7 && stressIndex <= 12))) {

        // first, calculate stress/strain tensor for each node of the element
        p_elstress->zero();
        int strInd = (stressIndex >= 0 && stressIndex <= 5) ? 0 : 1;
        packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                       *elDisp, strInd, surface,
                                       elemNodeTemps.data());

        // second, transform stress/strain tensor to nodal or material frame coordinates
        transformStressStrain(*p_elstress, iele, oframe);

        // third, extract the requested stress/strain value from the stress/strain tensor
        for (iNode = 0; iNode < NodesPerElement; ++iNode) {
          if(strInd == 0)
            (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex];
          else
            (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex-7]; 
        }

      }
      else {
        int index = (stressIndex == 31) ? 6 : stressIndex;
        // CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
        packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                      *elDisp, index, surface,
                                      elemNodeTemps.data(), ylayer, zlayer, avgnum);
      }

      if(avgnum > 0) {
        // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
        for(k = 0; k < NodesPerElement; ++k) {
          int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
          (*stress)[node] += (*elstress)[k];
          (*weight)[node] += (*elweight)[k];
        }
      }

    } // end of (printFlag != 2)

    // PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
    if(avgnum == 0 || avgnum == -1) {
      std::vector<size_t> offset(2);
      offset[0] = 0;
      offset[1] = packedEset[iele]->numTopNodes(); //HB 06-25-05: avoid the internal nodes for MpcElement

      if(stressIndex != 31) {
        if(printFlag == 0) {
          if(oinfo[fileNumber].nodeNumber == -1) {
            if(iele == 0) geoSource->outputElemStress(fileNumber, (double *) 0, 0, offset, time); // print time
            geoSource->outputElemStress(fileNumber, elstress->data(), 1, offset); // print stresses
          }
          else {
            geoSource->outputElemStress(fileNumber, elstress->data(), 1, offset, time); // print time and stresses
          }
        }
        if(printFlag == 1) {
          for(k = 0; k < NodesPerElement; ++k) {
            stressAllElems[count] = (*elstress)[k];
            count++;
          }
        }
        if(printFlag == 2) {
          if(iele == 0)
            geoSource->outputElemStress(fileNumber, (double *) 0, 0, offset, time); // print time
          count=count+NodesPerElement;
        }
      }
    }


  } // end of the iele loop

  // AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
  // ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.
  if(avgnum > 0) {

    if(printFlag != 2) {
    // assemble stress vector
     for(k = 0; k < numNodes; ++k)  {
       if((*weight)[k] == 0.0)
         (*stress)[k] = 0.0;
       else
         (*stress)[k] /= (*weight)[k];
     }
    }

    if(stressIndex != 31) {
      if(printFlag != 1) {
        if(oinfo[fileNumber].nodeNumber == -1)
          geoSource->outputNodeScalars(fileNumber, stress->data(), numNodes, time);
      else
        geoSource->outputNodeScalars(fileNumber, stress->data()+oinfo[fileNumber].nodeNumber, 1, time);
      }
    }

  }

  if(stressIndex == 31) {
    double stressmax(0);
    for(k = 0; k < numNodes; ++k) {
      if((*stress)[k]>stressmax) stressmax = (*stress)[k];
    }
    *aggregatedStress = 0.0;
    for(k = 0; k < numNodes; ++k) *aggregatedStress += exp(sinfo.ksParameter*((*stress)[k]-stressmax));
    *aggregatedStress = log(*aggregatedStress);
    *aggregatedStress /= sinfo.ksParameter;
    *aggregatedStress += stressmax;
    geoSource->outputEnergy(fileNumber, 0, *aggregatedStress);
  }

  delete [] nodeNumbers;
}

void
Domain::getStressStrain(ComplexVector &sol, DComplex *bcx, int fileNumber,
                        int stressIndex, double time, int printFlag)
{
  int numNodes = (outFlag) ? exactNumNodes : geoSource->numNode();

  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS, WHERE
  // ... ONLY THE AXIAL STRAIN (EXX) AND AXIAL STRESS (SXX) ARE CALCULATED

  OutputInfo *oinfo = geoSource->getOutputInfo();

  // avgnum = 2 --> do not include stress/strain of bar/beam element in averaging
  // avgnum = 3 --> only include elements whose nodes are in the ouputGroup in averaging
  int avgnum = oinfo[fileNumber].averageFlg;

  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;
  int surface = oinfo[fileNumber].surface;
    // upper  surface = 1
    // median surface = 2
    // lower  surface = 3
  OutputInfo::FrameType oframe = oinfo[fileNumber].oframe;

  int k;
  int iele;
  int iNode;

  double *nodalTemperatures = 0;
  // Either get the nodal temperatures from the input file or
  // from the thermal model
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >= 0) nodalTemperatures = temprcvd;

  ComplexVector *stress = 0;
  ComplexVector *stressAllElems = 0;
  ComplexVector *elstress = 0;
  ComplexVector *elDisp = 0;
  FullMC *p_elstress = 0;

  if(printFlag != 2) {
    // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
    if(avgnum != 0) {
      if(stress == 0) stress = new ComplexVector(numNodes,0.0);
      if(weight == 0) weight = new Vector(numNodes,0.0);
    }
    else if(printFlag == 1 && stressAllElems == 0) stressAllElems = new ComplexVector(sizeSfemStress,0.0);
    if(elDisp == 0) elDisp = new ComplexVector(maxNumDOFs,0.0);


    if((elstress == 0) || (elweight == 0) || (p_elstress == 0 && oframe != OutputInfo::Global)) {
      int NodesPerElement, maxNodesPerElement=0;
      for(iele=0; iele<numele; ++iele) {
        NodesPerElement = elemToNode->num(iele);
        maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
      }
      if(elstress == 0) elstress = new ComplexVector(maxNodesPerElement, 0.0);
      if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
      if(p_elstress == 0 && oframe != OutputInfo::Global) p_elstress = new FullMC(maxNodesPerElement,9);
    }

    if(avgnum != 0) {
    // zero the vectors
      stress->zero();
      weight->zero();
    }
    else if (printFlag == 1) stressAllElems->zero();
  }

  int count = 0;
  for(iele = 0; iele < numele; ++iele) {

    int NodesPerElement = elemToNode->num(iele);
    packedEset[iele]->nodes(nodeNumbers);

    if(printFlag != 2) {

      // Don't do anything if element is a phantom or constraint
      if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;

      // Don't include beams or bars in the averaging if nodalpartial (avgnum = 2) is requested
      if ((avgnum == 2 && packedEset[iele]->getElementType() == 6) ||
          (avgnum == 2 && packedEset[iele]->getElementType() == 7) ||
          (avgnum == 2 && packedEset[iele]->getElementType() == 1)) continue;

      // Don't include elements with one or more nodes not in the group if nodalpartialgroup (avgnum = 3) is requested
      if (avgnum == 3) { 
        int groupId = oinfo[fileNumber].groupNumber;
        if (groupId > 0) {
          std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
          std::set<int>::iterator it;
          for (iNode = 0; iNode < NodesPerElement; ++iNode)
            if((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
          if(it == groupNodes.end()) continue;
        } 
      }

      elDisp->zero();
      elstress->zero();
      elweight->zero();

      // DETERMINE ELEMENT DISPLACEMENT VECTOR
      for (k=0; k < allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if (cn >= 0)
          (*elDisp)[k] = sol[cn];
        else
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
      }

      for (iNode = 0; iNode < NodesPerElement; ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }

      // transform displacements from DOF_FRM to basic coordinates
      transformVectorInv(*elDisp, iele);

      // transform non-invariant stresses/strains from basic frame to DOF_FRM or CFRAME
      if(oframe != OutputInfo::Global && ((stressIndex >= 0 && stressIndex <= 5) || (stressIndex >= 7 && stressIndex <= 12))) {

        // first, calculate stress/strain tensor for each node of the element
        p_elstress->zero();
        int strInd = (stressIndex >= 0 && stressIndex <= 5) ? 0 : 1;
        packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                       *elDisp, strInd, surface,
                                       elemNodeTemps.data());

        // second, transform stress/strain tensor to nodal or material frame coordinates
        transformStressStrain(*p_elstress, iele, oframe);

        // third, extract the requested stress/strain value from the stress/strain tensor
        for (iNode = 0; iNode < NodesPerElement; ++iNode) {
          if(strInd == 0)
            (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex];
          else
            (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex-7];
        }

      }
      else {

        // CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
        packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                      *elDisp, stressIndex, surface,
                                      elemNodeTemps.data(), ylayer, zlayer, avgnum);
      }

      if(avgnum > 0) {
        // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
        for(k = 0; k < NodesPerElement; ++k) {
          int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
          (*stress)[node] += (*elstress)[k];
          (*weight)[node] += (*elweight)[k];
        }
      }

    } // end of (printFlag != 2)

    // PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
    if(avgnum == 0) {
      std::vector<size_t> offset(2);
      offset[0] = 0;
      offset[1] = packedEset[iele]->numTopNodes(); //HB 06-25-05: avoid the internal nodes for MpcElement

      if(printFlag == 0) {
        if(iele == 0)
          geoSource->outputElemStress(fileNumber, (DComplex *) 0, 0, offset, time); // print time
        geoSource->outputElemStress(fileNumber, elstress->data(), 1, offset); // print stresses
      }
      if(printFlag == 1) {
        for(k = 0; k < NodesPerElement; ++k) {
          stressAllElems[count] = (*elstress)[k];
          count++;
        }
      }
      if(printFlag == 2) {
        if(iele == 0)
          geoSource->outputElemStress(fileNumber, (DComplex *) 0, 0, offset, time); // print time
        count=count+NodesPerElement;
      }
    }


  } // end of the iele loop

  // AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
  // ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.
  if(avgnum > 0) {

    if(printFlag != 2) {
    // assemble stress vector
     for(k = 0; k < numNodes; ++k)  {
       if((*weight)[k] == 0.0)
         (*stress)[k] = 0.0;
       else
         (*stress)[k] /= (*weight)[k];
     }
    }

    if(printFlag != 1) {
     if(oinfo[fileNumber].nodeNumber == -1)
       geoSource->outputNodeScalars(fileNumber, stress->data(), numNodes, time);
     else
       geoSource->outputNodeScalars(fileNumber, stress->data()+oinfo[fileNumber].nodeNumber, 1, time);
    }

  }

  if(stress != 0) delete stress;
  if(stressAllElems != 0) delete stressAllElems;
  if(elstress != 0) delete elstress;
  if(elDisp != 0) delete elDisp;
  if(p_elstress != 0) delete p_elstress;
  delete [] nodeNumbers;
}

void
Domain::getPrincipalStress(Vector &sol, double *bcx, int fileNumber,
                           int stressIndex, double time)
{
  int numNodes = (outFlag) ? exactNumNodes : geoSource->numNode();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // set stress VS. strain for element subroutines
  int strInd;
  int strDir;
  if ((stressIndex==0)||(stressIndex==1)||(stressIndex==2)) {
    strInd = 0;
    strDir = stressIndex+1;
  }
  else if ((stressIndex==3)||(stressIndex==4)||(stressIndex==5)) {
    strInd = 1;
    strDir = stressIndex-2;
  }
  else {
    fprintf(stderr,"*** ERROR: Bad Principal Stress Direction ***\n");
    exit(-1);
  }

  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS

  int avgnum = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  int str_therm_option = oinfo[fileNumber].str_therm_option;

  // thermomechanical = 0
  // thermal = 1
  // mechanical = 2

  int j,k,iNode;
  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;
  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;
  // ... OUTPUT NODE NUMBER
  int n = oinfo[fileNumber].nodeNumber;

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(p_stress == 0) p_stress = new FullM(numNodes,6);
  if(weight == 0) weight = new Vector(numNodes,0.0);
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  int iele;
  if((p_elstress == 0)||(elweight == 0)) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    if(p_elstress == 0) p_elstress = new FullM(maxNodesPerElement,9);
    if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
  }

  // zero the vectors
  p_stress->zero();
  weight->zero();

  // ... WRITE CURRENT TIME VALUE
  if(avgnum == 0) {
    fprintf(oinfo[fileNumber].filptr,"  % *.*E\n",w,p,time);
  }

  for(iele=0; iele<numele; ++iele) {

    // Don't do anything if element is a phantom or constraint
    if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;

    int NodesPerElement = elemToNode->num(iele);
    packedEset[iele]->nodes(nodeNumbers);

    // Don't include elements with one or more nodes not in the group if nodalpartialgroup (avgnum = 3) is requested
    if (avgnum == 3) {
      int groupId = oinfo[fileNumber].groupNumber;
      if (groupId > 0) {
        std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
        std::set<int>::iterator it;
        for (iNode = 0; iNode < NodesPerElement; ++iNode)
          if((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
        if(it == groupNodes.end()) continue;
      }
    }

    elDisp->zero();
    p_elstress->zero();
    elweight->zero();

// ... DETERMINE ELEMENT DISPLACEMENT VECTOR

    if(str_therm_option != 1) {
      for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] = sol[cn];
        else
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
      }
    }

    for(iNode=0; iNode<NodesPerElement; ++iNode) {
      if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp || str_therm_option == 2)
        elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
      else
        elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
    }

    // transform displacements from DOF_FRM to basic coordinates
    transformVectorInv(*elDisp, iele);

// ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
    if(packedEset[iele]->getProperty())
      packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                     *elDisp, strInd, surface,
                                     elemNodeTemps.data());

// ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

    for(k = 0; k < NodesPerElement; ++k) {
      int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
      for(j = 0; j < 6; ++j) {
        (*p_stress)[node][j] += (*p_elstress)[k][j];
      }
      (*weight)[node] += (*elweight)[k];
    }

// ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
//     THIS WRITES THE CHOSEN PRINCIPAL STRESS FOR EACH ELEMENT
    if(avgnum == 0) {
      for(k=0; k<NodesPerElement; ++k)
         fprintf(oinfo[fileNumber].filptr," % *.*E",w,p,(*p_elstress)[k][5+strDir]);
       fprintf(oinfo[fileNumber].filptr,"\n");
    }
  }

// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

  if(avgnum > 0) {

    if(n == -1) {
      for(k=0; k<numNodes; ++k) {
        if((*weight)[k] == 0.0) {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] = 0.0;
          }
        } else  {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] /= (*weight)[k];
          }
        }
      }
    } else {
      if((*weight)[n] == 0.0) {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] = 0.0;
        }
      } else {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] /= (*weight)[n];
        }
      }
    }

    // CALCULATE PRINCIPALS AT EACH NODE
    double svec[6], pvec[3];
    if(n == -1) {
      double *globalPVec = new double[numNodes];
      for(k=0; k<numNodes; ++k) {
        for (j=0; j<6; ++j)
          svec[j] = (*p_stress)[k][j];

        // Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        globalPVec[k] = pvec[strDir-1];
      }
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodes, time);
    }
    else {
      for (j=0; j<6; ++j)
        svec[j] = (*p_stress)[n][j];

      // Convert Engineering to Tensor Strains
      if(strInd != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      geoSource->outputNodeScalars(fileNumber, pvec+strDir-1, 1, time);
    }

  }

  delete [] nodeNumbers;
}

void
Domain::getElementForces(Vector &sol, double *bcx, int fileNumber,
                         int forceIndex, double time)
{
  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();

  // ...REMARK:  FORCES ARE CALCULATED FOR BARS AND BEAMS CURRENTLY
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  // ... Watch out for this as the max number of nodes may increase
  int NodesPerElement = 2;

  FullM forces(numele, NodesPerElement);

  int i, iele;
  // note, we are reusing elstress here for an elemental force vector
  if(elstress == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    elstress = new Vector(maxNodesPerElement, 0.0);
  }

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  for(iele=0; iele<numele; ++iele) {
    packedEset[iele]->nodes(nodeNumbers);
    NodesPerElement = elemToNode->num(iele);

    // ... DETERMINE ELEMENT DISPLACEMENT VECTOR
    int k;
    for(k=0; k<allDOFs->num(iele); ++k) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
      if(cn >= 0)
        (*elDisp)[k] = sol[cn];
      else
        (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
    }

    // ... CALCULATE INTERNAL FORCE VALUE FOR EACH ELEMENT
    // ... taking into account the new temperature in case of thermal coupling
    if(packedEset[iele]->getProperty()) {
      int iNode;
      for(iNode=0; iNode<NodesPerElement; ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
    }

    // transform displacements from DOF_FRM to basic coordinates
    transformVectorInv(*elDisp, iele);

    if(geoSource->getOutputInfo()[fileNumber].oframe == OutputInfo::Local) {
      Vector fx(NodesPerElement,0.0), fy(NodesPerElement,0.0), fz(NodesPerElement,0.0);
      if(forceIndex == INX || forceIndex == INY || forceIndex == INZ) {
        packedEset[iele]->getIntrnForce(fx,nodes,elDisp->data(),INX,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fy,nodes,elDisp->data(),INY,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fz,nodes,elDisp->data(),INZ,elemNodeTemps.data());
      }
      else {
        packedEset[iele]->getIntrnForce(fx,nodes,elDisp->data(),AXM,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fy,nodes,elDisp->data(),AYM,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fz,nodes,elDisp->data(),AZM,elemNodeTemps.data());
      }
      Vector f(3*NodesPerElement);
      for(int iNode=0; iNode<NodesPerElement; iNode++) {
        double data[3] = { fx[iNode], fy[iNode], fz[iNode] };
        transformVector(data, nodeNumbers[iNode], false);
        switch(forceIndex) {
          case INX: case AXM: elstress[iNode] = data[0]; break;
          case INY: case AYM: elstress[iNode] = data[1]; break;
          case INZ: case AZM: elstress[iNode] = data[2]; break;
        }
      }
    }
    else {
      packedEset[iele]->getIntrnForce(*elstress,nodes,elDisp->data(),
                                      forceIndex,elemNodeTemps.data());
    }

    // ... COPY ELEMENT'S NODAL FORCES INTO A TOTAL FORCE MATRIX
    for(i=0; i<2; ++i) // note: element force output is currently only supported for 2-node elements
      forces[iele][i] = (*elstress)[i];

  }

  // ... PRINT THE ELEMENT FORCES TO A FILE
  geoSource->outputElemVectors(fileNumber, forces.data(), numele, time);

  delete [] nodeNumbers;
}

void
Domain::getKtimesU(Vector &dsp, double *bcx, Vector &ext_f, double eta,
                   FullSquareMatrix *kelArray)
{
  int size = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *karray = (double *) dbg_alloca(size);

  for(int iele=0; iele<numele; ++iele) {

     int numEleDOFs = allDOFs->num(iele);
     Vector elForce(numEleDOFs,0.0);     

     getElemKtimesU(iele,numEleDOFs,dsp,elForce.data(),kelArray,karray); 

     for(int k=0; k<numEleDOFs; ++k) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
       if(cn >= 0) {
         ext_f[cn] += eta * elForce[k]; }
     }
  }
}

void
Domain::getWeightedKtimesU(const std::map<int, double> &weights,
                          Vector &dsp, double *bcx, Vector &ext_f, double eta,
                          FullSquareMatrix *kelArray)
{
  int size = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *karray = (double *) dbg_alloca(size);

  for(std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
     const int iele = it->first;
     const double lumpingWeight = it->second;
 
     int numEleDOFs = allDOFs->num(iele);
     Vector elForce(numEleDOFs,0.0);

     getElemKtimesU(iele,numEleDOFs,dsp,elForce.data(),kelArray,karray);

     elForce *= lumpingWeight;

     for(int k=0; k<numEleDOFs; ++k) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
       if(cn >= 0) {
         ext_f[cn] += eta * elForce[k]; }
     }
  }
}

void 
Domain::getUnassembledKtimesU(const std::map<int, std::vector<int> > &weights,
                              Vector &dsp, double *bcx, Vector &ext_f, double eta,
                              FullSquareMatrix *kelArray)
{
  int size = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *karray = (double *) dbg_alloca(size);

  int uDofCounter = 0;

  for(std::map<int, std::vector<int> >::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
     const int iele = it->first;
     const std::vector<int> DOFvector(it->second);

     int numEleDOFs = allDOFs->num(iele);
     Vector elForce(numEleDOFs,0.0);

     getElemKtimesU(iele,numEleDOFs,dsp,elForce.data(),kelArray,karray);

     for(std::vector<int>::const_iterator DOFit = DOFvector.begin(); DOFit != DOFvector.end(); DOFit++) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][*DOFit]);
       if(cn >= 0) {
         ext_f[uDofCounter] += eta * elForce[*DOFit]; 
         uDofCounter += 1;
       }
     }
  }
}

void
Domain::getElemKtimesU(int iele, int numEleDOFs, Vector &dsp, double *elForce,
                       FullSquareMatrix *kelArray, double *karray)
{
  Vector elDisp(numEleDOFs,0.0);

  FullSquareMatrix kel(numEleDOFs,karray);

  for(int k=0; k<numEleDOFs; ++k) {
     int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
     if(cn >= 0) {
       elDisp[k] = dsp[cn];
     }
     else {
       elDisp[k] = 0.0;
     }
  }

  if(kelArray) {
    kel.copy(kelArray[iele]);
  }
  else {
    kel=packedEset[iele]->stiffness(nodes,karray);
  }

  for(int i=0;i<numEleDOFs;i++) {
    for(int j=0;j<numEleDOFs;j++) {
      elForce[i] += kel[i][j]*elDisp[j] ;
    }
  }
}

double
Domain::getStructureMass()
{
  double totmas = 0.0; // total mass of structure

  int iele;
  for(iele=0; iele < numele; ++iele) {
    double elementMass = packedEset[iele]->getMass(nodes);
    totmas += elementMass;
  }

  double *nodeCount = (double *) dbg_alloca(sizeof(double)*numnodes);
  double *nodeMass  = (double *) dbg_alloca(sizeof(double)*numnodes);

  int n;
  for(n=0; n<numnodes; ++n)
    nodeCount[n] = nodeMass[n] = 0.0;

  // ADD DISCRETE MASS
  DMassData *current = firstDiMass;
  while(current != 0) {
    int n = current->node;
    nodeMass[n]  += current->diMass;
    nodeCount[n] += 1;
    current = current->next;
  }

  for(n=0; n<numnodes; ++n) {
    if(nodeCount[n] == 0.0) continue;
    if(nodeCount[n] > 6) nodeCount[n] = 6;
    totmas += nodeMass[n]/nodeCount[n];
  }

   return totmas;
}

void
Domain::transformMatrix(FullSquareMatrix &kel, int iele)
{
  // transform element matrix from basic to DOF_FRM coordinates 
  if(domain->solInfo().basicDofCoords) return;

  if(packedEset[iele]->isMpcElement()) {
    LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]);
    if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
       lmpcons->getSource() == mpc::TiedSurfaces)) return;
  }

  // TODO: check for thermal/acoustic 
  // TODO: consider how to deal with elements that don't have either 3 translation dof on every node,
  //       or 3 translation + 3 rotation dofs on every node. E.g. torsional spring, ParallelAxesConstraint,
  //       and StraightLinePointFollowerConstraint
  //use: DofSetArray elem_dsa(packedEset[iele]);
  //     int dofs[6]; for(int i=0; i<packedEset[iele]->numNodes(); ++i) elem_dsa->number(i, DofSet::XYZdisp|DofSet::XYZrot, dofs); etc...
#ifdef USE_EIGEN3 
  NFrameData *nfd = geoSource->getNFrames();

  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >
    K(kel.data(),packedEset[iele]->numDofs(),packedEset[iele]->numDofs());

  int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
  int *nn = packedEset[iele]->nodes();

  for(int k = 0; k < numNodes; ++k) {
    int cd = nodes[nn[k]]->cd;
    if(cd == 0) continue;
    
    Eigen::Matrix3d T;
    T << nfd[cd].frame[0][0], nfd[cd].frame[0][1], nfd[cd].frame[0][2],
         nfd[cd].frame[1][0], nfd[cd].frame[1][1], nfd[cd].frame[1][2],
         nfd[cd].frame[2][0], nfd[cd].frame[2][1], nfd[cd].frame[2][2];

    if(packedEset[iele]->hasRot()) {

      Eigen::Matrix<double,6,6> TT;
      TT << T, Eigen::Matrix3d::Zero(),
            Eigen::Matrix3d::Zero(), T;

      for(int l=0; l<numNodes; ++l) {
        K.block<6,6>(6*k,6*l) = (TT*K.block<6,6>(6*k,6*l)).eval();
        K.block<6,6>(6*l,6*k) = (K.block<6,6>(6*l,6*k)*TT.transpose()).eval();
      }
      for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
        K.block<6,1>(6*k,6*numNodes+l) = (TT*K.block<6,1>(6*k,6*numNodes+l)).eval();
        K.block<1,6>(6*numNodes+l,6*k) = (K.block<1,6>(6*numNodes+l,6*k)*TT.transpose()).eval();
      }
    }
    else {
      for(int l=0; l<numNodes; ++l) {
        K.block<3,3>(3*k,3*l) = (T*K.block<3,3>(3*k,3*l)).eval();
        K.block<3,3>(3*l,3*k) = (K.block<3,3>(3*l,3*k)*T.transpose()).eval();
      }
      for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
        K.block<3,1>(3*k,3*numNodes+l) = (T*K.block<3,1>(3*k,3*numNodes+l)).eval();
        K.block<1,3>(3*numNodes+l,3*k) = (K.block<1,3>(3*numNodes+l,3*k)*T.transpose()).eval();
      }
    }
  }
  delete [] nn;
#endif
}

void
Domain::transformMatrixInv(FullSquareMatrix &kel, int iele)
{
  // transform element matrix from DOF_FRM to basic coordinates 
  if(domain->solInfo().basicDofCoords) return;

  if(packedEset[iele]->isMpcElement()) {
    LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]); 
    if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
       lmpcons->getSource() == mpc::TiedSurfaces)) return;
  }

  // TODO: check for thermal/acoustic 
  // TODO: consider how to deal with elements that don't have either 3 translation dof on every node,
  //       or 3 translation + 3 rotation dofs on every node. E.g. torsional spring, ParallelAxesConstraint,
  //       and StraightLinePointFollowerConstraint
  //use: DofSetArray elem_dsa(packedEset[iele]);
  //     int dofs[6]; for(int i=0; i<packedEset[iele]->numNodes(); ++i) elem_dsa->number(i, DofSet::XYZdisp|DofSet::XYZrot, dofs); etc...
#ifdef USE_EIGEN3 
  NFrameData *nfd = geoSource->getNFrames();

  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >
    K(kel.data(),packedEset[iele]->numDofs(),packedEset[iele]->numDofs());

  int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
  int *nn = packedEset[iele]->nodes();

  for(int k = 0; k < numNodes; ++k) {
    int cd = nodes[nn[k]]->cd;
    if(cd == 0) continue;
    
    Eigen::Matrix3d T;
    T << nfd[cd].frame[0][0], nfd[cd].frame[0][1], nfd[cd].frame[0][2],
         nfd[cd].frame[1][0], nfd[cd].frame[1][1], nfd[cd].frame[1][2],
         nfd[cd].frame[2][0], nfd[cd].frame[2][1], nfd[cd].frame[2][2];

    if(packedEset[iele]->hasRot()) {

      Eigen::Matrix<double,6,6> TT;
      TT << T, Eigen::Matrix3d::Zero(),
            Eigen::Matrix3d::Zero(), T;

      for(int l=0; l<numNodes; ++l) {
        K.block<6,6>(6*k,6*l) = (TT.transpose()*K.block<6,6>(6*k,6*l)).eval();
        K.block<6,6>(6*l,6*k) = (K.block<6,6>(6*l,6*k)*TT).eval();
      }
      for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
        K.block<6,1>(6*k,6*numNodes+l) = (TT.transpose()*K.block<6,1>(6*k,6*numNodes+l)).eval();
        K.block<1,6>(6*numNodes+l,6*k) = (K.block<1,6>(6*numNodes+l,6*k)*TT).eval();
      }
    }
    else {
      for(int l=0; l<numNodes; ++l) {
        K.block<3,3>(3*k,3*l) = (T.transpose()*K.block<3,3>(3*k,3*l)).eval();
        K.block<3,3>(3*l,3*k) = (K.block<3,3>(3*l,3*k)*T).eval();
      }
      for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
        K.block<3,1>(3*k,3*numNodes+l) = (T.transpose()*K.block<3,1>(3*k,3*numNodes+l)).eval();
        K.block<1,3>(3*numNodes+l,3*k) = (K.block<1,3>(3*numNodes+l,3*k)*T).eval();
      }
    }
  }
  delete [] nn;
#endif
}

void
Domain::transformVector(Vector &vec, int iele)
{
  // transform element vector from basic to DOF_FRM coordinates
  if(domain->solInfo().basicDofCoords) return;

  if(packedEset[iele]->isMpcElement()) {
    LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]); 
    if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
       lmpcons->getSource() == mpc::TiedSurfaces)) return;
  }

  int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
  int *nn = packedEset[iele]->nodes();
  if(packedEset[iele]->hasRot()) {
    for(int k=0; k<numNodes; ++k) 
      transformVector(vec.data()+6*k, nn[k], true);
  }    
  else {
    for(int k=0; k<numNodes; ++k)
      transformVector(vec.data()+3*k, nn[k], false);
  }    
  delete [] nn;
}

void
Domain::transformNeumVector(Vector &vec, int iele)
{
  // transform neumann b.c. vector from basic to DOF_FRM coordinates
  if(domain->solInfo().basicDofCoords) return;
  int numNodes = neum[iele]->numNodes();
  int *nn = neum[iele]->nodes();
  for(int k=0; k<numNodes; ++k)
    transformVector(vec.data()+3*k, nn[k], false);
  delete [] nn;
}

void
Domain::transformVector(ComplexVector &vec, int iele)
{
  // transform element vector from basic to DOF_FRM coordinates
  if(domain->solInfo().basicDofCoords) return;
  int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
  int *nn = packedEset[iele]->nodes();
  if(packedEset[iele]->hasRot()) {
    for(int k=0; k<numNodes; ++k)
      transformVector(vec.data()+6*k, nn[k], true);
  }
  else {
    for(int k=0; k<numNodes; ++k)
      transformVector(vec.data()+3*k, nn[k], false);
  }
  delete [] nn;
}

void
Domain::transformElementSensitivityInv(GenFullM<double> *dStressdDisp, int iele)
{
  // transform element vector from DOF_FRM to basic coordinates
  if(domain->solInfo().basicDofCoords) return;
  int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
  int *nn = packedEset[iele]->nodes();
  if(packedEset[iele]->hasRot()) {
    for(int k=0; k<numNodes; ++k)
      transformElementSensitivityInv(dStressdDisp->data()+6*k*dStressdDisp->numCol(), nn[k], numNodes, true);
  }
  else {
    for(int k=0; k<numNodes; ++k)
      transformElementSensitivityInv(dStressdDisp->data()+3*k*dStressdDisp->numCol(), nn[k], numNodes, false);
  }
  delete [] nn;
}

void
Domain::transformVectorInv(Vector &vec, int iele)
{
  // transform element vector from DOF_FRM to basic coordinates
  if(domain->solInfo().basicDofCoords) return;

  if(packedEset[iele]->isMpcElement()) {
    LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]); 
    if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
       lmpcons->getSource() == mpc::TiedSurfaces)) return;
  }

  int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
  int *nn = packedEset[iele]->nodes();
  if(packedEset[iele]->hasRot()) {
    for(int k=0; k<numNodes; ++k)            
      transformVectorInv(vec.data()+6*k, nn[k], true);
  }
  else {
    for(int k=0; k<numNodes; ++k)
      transformVectorInv(vec.data()+3*k, nn[k], false);
  }
  delete [] nn;
}

void
Domain::transformVectorInv(ComplexVector &vec, int iele)
{
  // transform element vector from DOF_FRM to basic coordinates
  if(domain->solInfo().basicDofCoords) return;
  int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
  int *nn = packedEset[iele]->nodes();
  if(packedEset[iele]->hasRot()) {
    for(int k=0; k<numNodes; ++k)
      transformVectorInv(vec.data()+6*k, nn[k], true);
  }
  else {
    for(int k=0; k<numNodes; ++k)
      transformVectorInv(vec.data()+3*k, nn[k], false);
  }
  delete [] nn;
}

void
Domain::transformVector(double *data, int inode, bool hasRot)
{
  // transform node vector from basic to DOF_FRM coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(hasRot) cd->transformVector6(data);
    else cd->transformVector3(data);
  }
}

void
Domain::transformVector(complex<double> *data, int inode, bool hasRot)
{
  // transform node vector from basic to DOF_FRM coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(hasRot) cd->transformVector6(data);
    else cd->transformVector3(data);
  }
}

void
Domain::transformVectorInv(double *data, int inode, bool hasRot)
{
  // transform node vector from DOF_FRM to basic coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(hasRot) cd->invTransformVector6(data);
    else cd->invTransformVector3(data);
  }
}

void
Domain::transformVectorInv(complex<double> *data, int inode, bool hasRot)
{
  // transform node vector from DOF_FRM to basic coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(hasRot) cd->invTransformVector6(data);
    else cd->invTransformVector3(data);
  }
}

void 
Domain::transformElementSensitivityInv(double *data, int inode, int numNodes, bool hasRot)
{
  // dStressdDisp->data()+6*k*dStressdDisp->numCol(), nn[k], true);
  // transform node vector from DOF_FRM to basic coordinates
  if(inode >= numnodes || nodes[inode] == NULL) return;
  int cd = nodes[inode]->cd;
  if(cd == 0) return;

  NFrameData *nfd = geoSource->getNFrames();
#ifdef USE_EIGEN3
  Eigen::Matrix3d T;
  T << nfd[cd].frame[0][0], nfd[cd].frame[0][1], nfd[cd].frame[0][2],
       nfd[cd].frame[1][0], nfd[cd].frame[1][1], nfd[cd].frame[1][2],
       nfd[cd].frame[2][0], nfd[cd].frame[2][1], nfd[cd].frame[2][2];

  if(hasRot) {
    for(int row = 0; row < numNodes ; ++row) {
      Eigen::Matrix<double,6,1> v;
      v << data[row], data[row+numNodes], data[row+2*numNodes], data[row+3*numNodes], data[row+4*numNodes], data[row+5*numNodes];
      v.head<3>() = (T*v.head<3>()).eval();
      v.tail<3>() = (T*v.tail<3>()).eval();
      data[row] = v[0];            data[row+numNodes] = v[1];
      data[row+2*numNodes] = v[2]; data[row+3*numNodes] = v[3];
      data[row+4*numNodes] = v[4]; data[row+5*numNodes] = v[5];
    }
  }
  else {
    for(int row = 0; row < numNodes ; ++row) {
      Eigen::Matrix<double,3,1> v;
      v << data[row], data[row+numNodes], data[row+2*numNodes]; 
      v = (T*v).eval();
      data[row] = v[0];      data[row+numNodes] = v[1];     data[row+2*numNodes] = v[2];
    }
  }
#endif
}

void
Domain::transformStressStrain(FullM &mat, int iele, OutputInfo::FrameType oframe)
{
  // transform element stress or strain tensors from basic to DOF_FRM or CFRAME coordinates
  switch(oframe) {
    default:
    case OutputInfo::Local: {
      if(domain->solInfo().basicDofCoords) return;
      int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
      int *nn = packedEset[iele]->nodes();
      for(int k=0; k<numNodes; ++k)
        transformMatrix(mat[k], nn[k]);
      delete [] nn;
    } break;
    case OutputInfo::Global: 
      break;
    case OutputInfo::Material: {
#ifdef USE_EIGEN3
      double cFrame[3][3];
      packedEset[iele]->getCFrame(nodes, cFrame);
      Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&cFrame[0][0]);
      int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
      for(int k=0; k<numNodes; ++k) {
        Eigen::Matrix<double,3,3> M;
        M << mat[k][0], mat[k][3], mat[k][5],
             mat[k][3], mat[k][1], mat[k][4],
             mat[k][5], mat[k][4], mat[k][2];

        M = T*M*T.transpose();

        mat[k][0] = M(0,0);
        mat[k][1] = M(1,1);
        mat[k][2] = M(2,2);
        mat[k][3] = M(0,1);
        mat[k][4] = M(1,2);
        mat[k][5] = M(0,2);
      }
#else
      std::cerr << " *** WARNING: USE_EIGEN3 is not defined in Domain::transformStressStrain\n";
#endif
    } break;
  }
}

void
Domain::transformStressStrain(FullMC &mat, int iele, OutputInfo::FrameType oframe)
{
  // transform element stress or strain tensors from basic to DOF_FRM or CFRAME coordinates
  switch(oframe) {
    default:
    case OutputInfo::Local: {
      if(domain->solInfo().basicDofCoords) return;
      int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
      int *nn = packedEset[iele]->nodes();
      for(int k=0; k<numNodes; ++k)
        transformMatrix(mat[k], nn[k]);
      delete [] nn;
    } break;
    case OutputInfo::Global:
      break;
    case OutputInfo::Material: {
#ifdef USE_EIGEN3
      double cFrame[3][3];
      packedEset[iele]->getCFrame(nodes, cFrame);
      Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&cFrame[0][0]);
      int numNodes = packedEset[iele]->numNodes()-packedEset[iele]->numInternalNodes();
      for(int k=0; k<numNodes; ++k) {
        Eigen::Matrix<complex<double>,3,3> M;
        M << mat[k][0], mat[k][3], mat[k][5],
             mat[k][3], mat[k][1], mat[k][4],
             mat[k][5], mat[k][4], mat[k][2];

        M = T*M*T.transpose();

        mat[k][0] = M(0,0);
        mat[k][1] = M(1,1);
        mat[k][2] = M(2,2);
        mat[k][3] = M(0,1);
        mat[k][4] = M(1,2);
        mat[k][5] = M(0,2);
      }
#else
      std::cerr << " *** WARNING: USE_EIGEN3 is not defined in Domain::transformStressStrain\n";
#endif
    } break;
  }
}

void
Domain::transformMatrix(double *data, int inode, bool sym)
{
  // transform node 3x3 matrix from basic to DOF_FRM coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(sym) cd->transformSymMatrix3(data);
    else cd->transformMatrix3(data);
  }
}

void
Domain::transformMatrix(complex<double> *data, int inode, bool sym)
{
  // transform node 3x3 matrix from basic to DOF_FRM coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(sym) cd->transformSymMatrix3(data);
    else cd->transformMatrix3(data);
  }
}

void
Domain::transformMatrixInv(double *data, int inode, bool sym)
{
  // transform node 3x3 matrix from DOF_FRM to basic coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(sym) cd->invTransformSymMatrix3(data);
    else cd->invTransformMatrix3(data);
  }
}

void
Domain::transformMatrixInv(complex<double> *data, int inode, bool sym)
{
  // transform node 3x3 matrix from DOF_FRM to basic coordinates
  if(NFrameData *cd = nodes.dofFrame(inode)) {
    if(sym) cd->invTransformSymMatrix3(data);
    else cd->invTransformMatrix3(data);
  }
}

void
Domain::computeWeightWRTShapeVariableSensitivity(int sindex, AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
     // ... COMPUTE TOTAL STRUCTURAL WEIGHT AND DERIVATIVE WRT SHAPE VARIABLES
     double weight(0);
     allSens.weightWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numShapeVars);
     allSens.weightWRTshape->setZero();
     std::map<int, Attrib> &attributes = geoSource->getAttributes();
     for(int iele = 0; iele < numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int nnodes = packedEset[iele]->numNodes();
       Vector weightDerivative(3*nnodes);
       StructProp *prop = packedEset[iele]->getProperty();
       if(prop == 0) continue; // phantom element

       weight += packedEset[iele]->weight(nodes, gravityAcceleration);
       packedEset[iele]->getWeightNodalCoordinateSensitivity(weightDerivative, nodes, gravityAcceleration);
       for(int ishap=0; ishap<numShapeVars; ++ishap) {
         for(int i=0; i<nnodes; ++i) {
           int node2 = (outFlag) ? nodeTable[(*elemToNode)[iele][i]]-1 : (*elemToNode)[iele][i];
           int inode = shapeSenData.nodes[node2];
           for(int xyz=0; xyz<3; ++xyz) {
             (*allSens.weightWRTshape)[ishap] += weightDerivative[3*i+xyz] * shapeSenData.sensitivities[ishap][node2][xyz];
           }
         }
       }
     }

     allSens.weight = weight;
#endif
}

void
Domain::computeWeightWRTthicknessSensitivity(int sindex, AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
     // ... COMPUTE TOTAL STRUCTURAL WEIGHT AND DERIVATIVE WRT THICKNESS
     int iele;
     double weight(0);
     std::map<int, Group> &group = geoSource->group;
     std::map<int, AttributeToElement> &atoe = geoSource->atoe;

     allSens.weightWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numThicknessGroups);
     allSens.weightWRTthick->setZero();
     std::map<int, Attrib> &attributes = geoSource->getAttributes();

     for(iele = 0; iele < numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       StructProp *prop = packedEset[iele]->getProperty();
       if(prop == 0) continue; // phantom element

       weight += packedEset[iele]->weight(nodes, gravityAcceleration);
     }

     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       int groupIndex = thicknessGroups[iparam];
       for(int aindex = 0; aindex < group[groupIndex].attributes.size(); ++aindex) {
         for(int eindex =0; eindex < atoe[group[groupIndex].attributes[aindex]].elems.size(); ++eindex) {
           iele = atoe[group[groupIndex].attributes[aindex]].elems[eindex];
           (*allSens.weightWRTthick)[iparam] += packedEset[iele]->getWeightThicknessSensitivity(nodes, gravityAcceleration);
         }
       }
     }
     allSens.weight = weight;
#endif
}

void 
Domain::computeStiffnessWRTthicknessSensitivity(int sindex, AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
     // ... COMPUTE SENSITIVITY OF STIFFNESS MATRIX WRT THICKNESS
     allSens.stiffnessWRTthickSparse = new GenSparseMatrix<double>*[numThicknessGroups];
     allSens.dKucdthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
//     allSens.dKucdthickSparse = new GenSparseMatrix<double>*[numThicknessGroups];
     std::map<int, Group> &group = geoSource->group;
     std::map<int, AttributeToElement> &atoe = geoSource->atoe;
     for(int g=0; g<numThicknessGroups; ++g) {
       allSens.stiffnessWRTthickSparse[g] =
           constructEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >(
               c_dsa, nodeToNode.get(), false);
       allSens.dKucdthick[g] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numUncon(),numDirichlet); 
       allSens.dKucdthick[g]->setZero();
//       allSens.dKucdthickSparse[g] = constructCuCSparse<double>(); 
     } 


     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *stifWRTthic
           = dynamic_cast<GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *>(allSens.stiffnessWRTthickSparse[iparam]);
       int groupIndex = thicknessGroups[iparam];
       for(int aindex = 0; aindex < group[groupIndex].attributes.size(); ++aindex) {
         for(int eindex =0; eindex < atoe[group[groupIndex].attributes[aindex]].elems.size(); ++eindex) {
           int iele = atoe[group[groupIndex].attributes[aindex]].elems[eindex];
           if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
           int DofsPerElement = packedEset[iele]->numDofs();
           FullSquareMatrix dStiffnessdThick(DofsPerElement);
           packedEset[iele]->getStiffnessThicknessSensitivity(nodes, dStiffnessdThick,1);
           auto dofs = (*allDOFs)[iele];
           allSens.stiffnessWRTthickSparse[iparam]->add(dStiffnessdThick, dofs);
           const auto &unconstrNum = c_dsa->getUnconstrNum();
           const auto &constrndNum = c_dsa->getConstrndNum();
           for(int k = 0; k < DofsPerElement; ++k) {
             int dofk = unconstrNum[dofs[k]];
             if(dofs[k] < 0 || dofk < 0) continue;  // Skip undefined/constrained dofs
             for(int j = 0; j < DofsPerElement; ++j) {
               int dofj = constrndNum[dofs[j]];
               if(dofj == -1) continue;
               (*allSens.dKucdthick[iparam])(dofk, dofj) += dStiffnessdThick[k][j];
             }
           }
         }
       }
     }
#endif
}


void 
Domain::computeStiffnessWRTShapeVariableSensitivity(int sindex, AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
     // ... COMPUTE SENSITIVITY OF STIFFNESS MATRIX WRT NODAL COORDINATES
     allSens.stiffnessWRTshapeSparse = new GenSparseMatrix<double>*[numShapeVars];
     allSens.dKucdshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
     allSens.dKucdshapeSparse = new GenSparseMatrix<double>*[numShapeVars]; 
     for(int g=0; g<numShapeVars; ++g) {
       allSens.stiffnessWRTshapeSparse[g] =
           constructEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >(
               c_dsa, nodeToNode.get(), false);
       allSens.dKucdshape[g] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numUncon(),numDirichlet); 
       allSens.dKucdshape[g]->setZero();
       allSens.dKucdshapeSparse[g] = constructCuCSparse<double>(); 
     } 

     for(int iele = 0; iele < numele; iele++) { 
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int DofsPerElement = packedEset[iele]->numDofs();
       int nnodes = packedEset[iele]->numNodes();
       FullSquareMatrix *dStiffnessdCoord = new FullSquareMatrix[3*nnodes];
       for(int i=0; i<3*nnodes; ++i) { dStiffnessdCoord[i].setSize(DofsPerElement); dStiffnessdCoord[i].zero(); }
       packedEset[iele]->getStiffnessNodalCoordinateSensitivity(dStiffnessdCoord, nodes);
       // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
       auto dofs = (*allDOFs)[iele];
       const auto &unconstrNum = c_dsa->getUnconstrNum();
       const auto &constrndNum = c_dsa->getConstrndNum();
       for(int isen = 0; isen < numShapeVars; ++isen) { 
         GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *stifWRTsha = dynamic_cast<GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *>(allSens.stiffnessWRTshapeSparse[isen]);
         for(int i = 0; i < nnodes; ++i) {
           int node2 = (outFlag) ? nodeTable[(*elemToNode)[iele][i]]-1 : (*elemToNode)[iele][i];
           for(int k = 0; k < DofsPerElement; ++k) {
             int dofk = unconstrNum[dofs[k]];
             if(dofs[k] < 0 || dofk < 0) continue;  // Skip undefined/constrained dofs
             for(int j = 0; j < DofsPerElement; ++j) {
               int dofj = unconstrNum[dofs[j]];
               if(dofs[j] < 0 || dofj < 0) continue;  // Skip undefined/constrained dofs
               for(int xyz = 0; xyz < 3; ++xyz) { 
                 stifWRTsha->addCoef(dofk,dofj,dStiffnessdCoord[3*i+xyz][k][j]*shapeSenData.sensitivities[isen][node2][xyz]);
               }
             }
             for(int j = 0; j < DofsPerElement; ++j) {
               int dofj = constrndNum[dofs[j]];
               if(dofj == -1) continue;
               for(int xyz = 0; xyz < 3; ++xyz)
                 (*allSens.dKucdshape[isen])(dofk, dofj) += dStiffnessdCoord[3*i+xyz][k][j]*shapeSenData.sensitivities[isen][node2][xyz];
             }
           }
         }
       }
       delete [] dStiffnessdCoord;
     }
#endif
}


void
Domain::makePreSensitivities(AllSensitivities<double> &allSens, double *bcx)
{
#ifdef USE_EIGEN3
 for(int sindex=0; sindex < numSensitivity; ++sindex) {
  switch(senInfo[sindex].type) {
   case SensitivityInfo::WeightWRTthickness:
   {
     computeWeightWRTthicknessSensitivity(sindex, allSens); 
     break;
   }
   case SensitivityInfo::WeightWRTshape:
   {
     computeWeightWRTShapeVariableSensitivity(sindex, allSens);
     break;
   }
   default:
     break;
  }
 }
#endif
}

void 
Domain::computeLinearStaticWRTthicknessSensitivity(int sindex, 
                                                   AllSensitivities<double> &allSens,
                                                   GenVector<double> *sol,
                                                   GeomState *refState,
                                                   GeomState *geomState,
                                                   Corotator **allCorot)
{
#ifdef USE_EIGEN3
     allSens.linearstaticWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
     std::map<int, Group> &group = geoSource->group;
     std::map<int, AttributeToElement> &atoe = geoSource->atoe;
     for(int g=0; g<numThicknessGroups; ++g) {
       allSens.linearstaticWRTthick[g] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numUncon(),1);
       allSens.linearstaticWRTthick[g]->setZero();   
     }
     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > disp(sol->data(),numUncon(),1);
       if(allSens.stiffnessWRTthickSparse) {
         Eigen::MappedSparseMatrix<double, Eigen::ColMajor, int> dKdthick = 
               dynamic_cast<GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,
                            Eigen::Upper> > *>(allSens.stiffnessWRTthickSparse[iparam])->getEigenSparse();
         *allSens.linearstaticWRTthick[iparam] = dKdthick * disp;
       } else {
         std::cerr << "ERROR! stiffnessWRTthickSparse is not defined yet\n";
         exit(-1);
       }
       if(numDirichlet) {
         Eigen::Matrix<double, Eigen::Dynamic, 1> Vc(numDirichlet);
         Vc.setZero();
         for(int i = 0; i < numDirichlet; ++i) {
           int dof = dsa->locate(dbc[i].nnum, (1 << dbc[i].dofnum));
           if(dof < 0) continue;
           int dof2 = c_dsa->invRCN(dof);
           if(dof2 >= 0) Vc[dof2] = dbc[i].val;
         }
         *allSens.linearstaticWRTthick[iparam] += (*allSens.dKucdthick[iparam]) * Vc;
       }
     }
     if(domain->gravityFlag()) subtractGravityForceSensitivityWRTthickness(sindex,allSens); // TODO-> must consider
                                                                                            // other external forces
                                                                                            // that depend on thickness
#endif
}

void
Domain::computeLinearStaticWRTShapeVariableSensitivity(int sindex,
                                                       AllSensitivities<double> &allSens,
                                                       GenVector<double> *sol)
{
#ifdef USE_EIGEN3
     allSens.linearstaticWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
     if(numShapeVars != shapeSenData.numVars) {
       std::cerr << " *** ERROR: number of shape variables is not equal to the one in Shape Derivative file\n";
       exit(-1);
     }
     for(int s=0; s<numShapeVars; ++s) {
       allSens.linearstaticWRTshape[s] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numUncon(),1);
       allSens.linearstaticWRTshape[s]->setZero();   // this line and next line are for DEBUG purpose
     }
     for(int ishape = 0; ishape < numShapeVars; ++ishape) {
       Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > disp(sol->data(),numUncon(),1);
       if(allSens.stiffnessWRTshapeSparse) {
         Eigen::MappedSparseMatrix<double, Eigen::ColMajor, int> M = 
            dynamic_cast<GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,
                         Eigen::Upper> > *>(allSens.stiffnessWRTshapeSparse[ishape])->getEigenSparse();
         *allSens.linearstaticWRTshape[ishape] = M * disp;
       } else {
         std::cerr << "ERROR! stiffnessWRTshape is not defined yet\n";
         exit(-1);
       }
       if(numDirichlet) {
         Eigen::Matrix<double, Eigen::Dynamic, 1> Vc(numDirichlet);
         Vc.setZero();
         for(int i = 0; i < numDirichlet; ++i) {
           int dof = dsa->locate(dbc[i].nnum, (1 << dbc[i].dofnum));
           if(dof < 0) continue;
           int dof2 = c_dsa->invRCN(dof);
           if(dof2 >= 0) Vc[dof2] = dbc[i].val;
         }
         *allSens.linearstaticWRTshape[ishape] += (*allSens.dKucdshape[ishape]) * Vc;
       }
     }
     if(domain->gravityFlag()) {
       subtractGravityForceSensitivityWRTShapeVariable(sindex,allSens);
      } //TODO-> must consider other external forces
#endif
}

void 
Domain::subtractGravityForceSensitivityWRTthickness(int sindex, AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
     Vector elementGravityForceSen(maxNumDOFs);
     int gravflg;
     std::map<int, Group> &group = geoSource->group;
     std::map<int, AttributeToElement> &atoe = geoSource->atoe;
     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
      int groupIndex = thicknessGroups[iparam];
      for(int aindex = 0; aindex < group[groupIndex].attributes.size(); ++aindex) {  
       for(int eindex = 0; eindex < atoe[group[groupIndex].attributes[aindex]].elems.size(); ++eindex) {
         int iele = atoe[group[groupIndex].attributes[aindex]].elems[eindex];
         if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
         if(packedEset[iele]->getProperty() == 0) continue; // phantom element
         if(geoSource->consistentQFlag() && !(sinfo.isDynam() && packedEset[iele]->getMassType() == 0))
           gravflg = 2;
         else gravflg = geoSource->fixedEndM;
         elementGravityForceSen.zero();
         packedEset[iele]->getGravityForceThicknessSensitivity(nodes, gravityAcceleration, elementGravityForceSen, gravflg);
         // transform vector from basic to DOF_FRM coordinates
         transformVector(elementGravityForceSen, iele);
  
         for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
           int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
           if(cn >= 0) {
             (*allSens.linearstaticWRTthick[iparam])(cn,0) -= elementGravityForceSen[idof];
           } 
         }      
       }
      }
     }
#endif
}

void 
Domain::subtractGravityForceSensitivityWRTShapeVariable(int sindex, AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
     int gravflg;
     if(numShapeVars != shapeSenData.numVars) {
       std::cerr << " *** ERROR: number of shape variables is not equal to the one in Shape Derivative file\n"; 
       exit(-1);
     }

     for(int iele = 0; iele < numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int DofsPerElement = packedEset[iele]->numDofs();
       int NodesPerElement = elemToNode->num(iele);
       GenFullM<double> elementGravityForceSen(3*NodesPerElement,DofsPerElement,double(0.0));
       if(packedEset[iele]->getProperty() == 0) continue; // phantom element
       if(geoSource->consistentQFlag() && !(sinfo.isDynam() && packedEset[iele]->getMassType() == 0))
         gravflg = 2;
       else gravflg = geoSource->fixedEndM;
       elementGravityForceSen.zero();
       packedEset[iele]->getGravityForceNodalCoordinateSensitivity(nodes, gravityAcceleration, elementGravityForceSen, gravflg);

       // transform vector from basic to DOF_FRM coordinates
//       transformVector(elementGravityForceSen, iele); //TODO->commented out for now, but if nodes does not use basic coordinate frame, then shouldn't comment this out
 
       for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
         if(cn >= 0) {
           for(int k = 0; k < NodesPerElement; ++k) {
             int node1 = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
             int inode = shapeSenData.nodes[node1];
             for(int xyz = 0; xyz < 3; ++xyz) {
               for(int ishap = 0; ishap < numShapeVars; ++ishap) {
                 (*allSens.linearstaticWRTshape[ishap])(cn,0) -= elementGravityForceSen[3*k+xyz][idof]*shapeSenData.sensitivities[ishap][node1][xyz];
               }
             }
           } 
         } 
       }      
     }
#endif
}

void
Domain::computeStressVMDualSensitivity(int sindex,
                                     GenSolver<double> *sysSolver,
                                     AllSensitivities<double> &allSens,
                                     GenSparseMatrix<double> *spm,
                                     GenSparseMatrix<double> *K)
{
#ifdef USE_EIGEN3
     std::map<int, AttributeToElement> &atoe = geoSource->atoe;
     allSens.lambdaStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>*[numStressNodes];
     for(int inode = 0; inode < numStressNodes; ++inode) {
       allSens.lambdaStressVM[inode] = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numUncon());
       Vector rhs(numUncon(),0.0), lambdaStress(numUncon(),0.0);
       for(int i=0; i<numUncon(); ++i) rhs[i] = (*allSens.vonMisesWRTdisp)(stressNodes[inode],i);
       sysSolver->solve(rhs,lambdaStress);
       *allSens.lambdaStressVM[inode] = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> >(lambdaStress.data(),numUncon(),1);
     }
#endif
}

void
Domain::computeAggregatedStressVMDualSensitivity(int sindex,
                                                GenSolver<double> *sysSolver,
                                                AllSensitivities<double> &allSens,
                                                GenSparseMatrix<double> *spm,
                                                GenSparseMatrix<double> *K)
{
#ifdef USE_EIGEN3
     allSens.lambdaAggregatedStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numUncon());
     Vector rhs(numUncon(),0.0), lambdaAggregatedStress(numUncon(),0.0);
     for(int i=0; i<numUncon(); ++i) rhs[i] = (*allSens.aggregatedVonMisesWRTdisp)(i);
     sysSolver->solve(rhs, lambdaAggregatedStress);
     *allSens.lambdaAggregatedStressVM = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> >(lambdaAggregatedStress.data(),numUncon(),1);
#endif
}

int 
Domain::returnLocalDofNum(int node, int dof)
{
  switch(dof) {
    case 0:
      return c_dsa->locate( node, DofSet::Xdisp);
    case 1:
      return c_dsa->locate( node, DofSet::Ydisp);
    case 2:
      return c_dsa->locate( node, DofSet::Zdisp);
    case 3:
      return c_dsa->locate( node, DofSet::Xrot);
    case 4:
      return c_dsa->locate( node, DofSet::Yrot);
    case 5:
      return c_dsa->locate( node, DofSet::Zrot);
  }
  return 0; 
}

void
Domain::computeDisplacementDualSensitivity(int sindex,
                                           GenSolver<double> *sysSolver,
                                           AllSensitivities<double> &allSens,
                                           GenSparseMatrix<double> *spm,
                                           GenSparseMatrix<double> *K)
{
#ifdef USE_EIGEN3
     allSens.lambdaDisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>*[numTotalDispDofs];
     int dispDofIndex = 0;
     for(int inode=0; inode<numDispNodes; ++inode) {
       int node = dispNodes[inode].nodeID, loc;
       int numDispDofs = dispNodes[inode].numdofs;
       for(int idof=0; idof<numDispDofs; ++idof) {
         allSens.lambdaDisp[dispDofIndex] = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numUncon());
         allSens.lambdaDisp[dispDofIndex]->setZero();
         Vector rhs(numUncon(),0.0), lambdadisp(numUncon(),0.0);
         int dof = dispNodes[inode].dofs[idof];
         loc = returnLocalDofNum(node, dof);
         if (loc >= 0) rhs[loc] = 1.0;
         else continue;
         sysSolver->solve(rhs,lambdadisp);
         *allSens.lambdaDisp[dispDofIndex] = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> >(lambdadisp.data(),numUncon(),1);
         dispDofIndex++;
       }
     }
#endif
}

void
Domain::computeDisplacementWRTthicknessAdjointSensitivity(int sindex,
                                                          AllSensitivities<double> &allSens,
                                                          GenSparseMatrix<double> *spm)
{
#ifdef USE_EIGEN3
     allSens.dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       allSens.dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numTotalDispDofs,1);
       allSens.dispWRTthick[iparam]->setZero();
       for(int idof=0; idof<numTotalDispDofs; ++idof) {
         (*allSens.dispWRTthick[iparam])(idof,0) -= allSens.lambdaDisp[idof]->adjoint()*(allSens.linearstaticWRTthick[iparam]->col(0));
       }
     }
#endif
}

void
Domain::computeDisplacementWRTShapeVariableAdjointSensitivity(int sindex,
                                                              AllSensitivities<double> &allSens,
                                                              GenSparseMatrix<double> *spm)
{
#ifdef USE_EIGEN3
     allSens.dispWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
     for(int iparam = 0; iparam < numShapeVars; ++iparam) {
       allSens.dispWRTshape[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numTotalDispDofs,1);
       allSens.dispWRTshape[iparam]->setZero();
       for(int idof=0; idof<numTotalDispDofs; ++idof) {
         double a = allSens.lambdaDisp[idof]->adjoint()*(allSens.linearstaticWRTshape[iparam]->col(0));
         (*allSens.dispWRTshape[iparam])(idof,0) -= a;
       }
     }
#endif
}

void
Domain::computeDisplacementWRTthicknessDirectSensitivity(int sindex,
                                                         GenSolver<double> *sysSolver,
                                                         AllSensitivities<double> &allSens,
                                                         GenSparseMatrix<double> *spm,
                                                         GenSparseMatrix<double> *K)
{
#ifdef USE_EIGEN3
     std::map<int, AttributeToElement> &atoe = geoSource->atoe;
     allSens.dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       allSens.dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numUncon(),1);
       Vector rhs(allSens.linearstaticWRTthick[iparam]->data(), numUncon()), sol(numUncon(),0.0);
       rhs *= -1;
       sysSolver->solve(rhs,sol);
       *allSens.dispWRTthick[iparam] = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >(sol.data(),numUncon(),1);
     }
#endif
}

void
Domain::computeDisplacementWRTShapeVariableDirectSensitivity(int sindex,
                                                             GenSolver<double> *sysSolver,
                                                             AllSensitivities<double> &allSens,
                                                             GenSparseMatrix<double> *spm,
                                                             GenSparseMatrix<double> *K)
{
#ifdef USE_EIGEN3
     allSens.dispWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
     for(int ishap = 0; ishap < numShapeVars; ++ishap) {
       allSens.dispWRTshape[ishap] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numUncon(),1);
       Vector rhs(allSens.linearstaticWRTshape[ishap]->data(), numUncon()), sol(numUncon(),0.0);
       rhs *= -1;
       sysSolver->solve(rhs,sol);
       *allSens.dispWRTshape[ishap] = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >(sol.data(),numUncon(),1);
     }
#endif
}
                                      
void 
Domain::computeStressVMWRTMachNumberSensitivity(AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
   // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO Mach Number
   allSens.vonMisesWRTmach = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), 1);
   allSens.vonMisesWRTmach->setZero();
#endif
}

void 
Domain::computeStressVMWRTangleOfAttackSensitivity(AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
   // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO Mach Number
   allSens.vonMisesWRTalpha = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), 1);
   allSens.vonMisesWRTalpha->setZero();
#endif
}

void 
Domain::computeStressVMWRTyawAngleSensitivity(AllSensitivities<double> &allSens)
{
#ifdef USE_EIGEN3
   // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO Mach Number
   allSens.vonMisesWRTbeta = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), 1);
   allSens.vonMisesWRTbeta->setZero();
#endif
}
                                      
void 
Domain::computeStressVMWRTthicknessDirectSensitivity(int sindex, AllSensitivities<double> &allSens, GenVector<double> *sol, double *bcx, 
                                               GeomState *refState, GeomState *geomState, Corotator **allCorot, bool isDynam)
{
#ifdef USE_EIGEN3
     // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO THICKNESS
     allSens.vonMisesWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), numThicknessGroups);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numNodes(), 1);
     allSens.vonMisesWRTthick->setZero();     stressWeight.setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     
     for(int iele = 0; iele< numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenVector<double> dStressdThick(NodesPerElement);
       GenVector<double> weight(NodesPerElement,0.0);
       int surface = senInfo[sindex].surface;
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele);
       if(thgreleFlag[iele]) { 
         packedEset[iele]->getVonMisesThicknessSensitivity(dStressdThick, weight, nodes, *elDisp, 6, surface); 
       } else {
         dStressdThick.zero();
         Vector elstress(NodesPerElement,0.0);
         packedEset[iele]->getVonMises(elstress, weight, nodes, *elDisp, -1, surface, NULL, 0, 0, 0);
       }
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           int iparam = thpaIndex[iele];
           if(iparam >= 0) (*allSens.vonMisesWRTthick)(node, iparam) += dStressdThick[k];
           stressWeight(node, 0) += weight[k];
         }
       }
     } 

     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
      // average vonMisesWRTthick vector
      for(int inode = 0; inode <numNodes(); ++inode)  {
        if(stressWeight(inode, 0) == 0.0)
          (*allSens.vonMisesWRTthick)(inode, iparam) = 0.0;
        else
          (*allSens.vonMisesWRTthick)(inode, iparam) /= stressWeight(inode, 0);
      }
      if(!isDynam) allSens.vonMisesWRTthick->col(iparam) += *allSens.vonMisesWRTdisp * (*allSens.dispWRTthick[iparam]);
     }
#endif
}

void 
Domain::computeAggregatedStressVMWRTthicknessSensitivity(int sindex,
                                                         AllSensitivities<double> &allSens,
                                                         GenVector<double> *sol,
                                                         double *bcx, 
                                                         bool isDynam)
{
#ifdef USE_EIGEN3
     // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO THICKNESS
     int surface = senInfo[sindex].surface;
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressvmWRTthick(numNodes(), numThicknessGroups);
     Eigen::Matrix<double, Eigen::Dynamic, 1> stressWeight(numNodes());
     Vector stress(numNodes(),0.0);
     computeNormalizedVonMisesStress(*sol, bcx, surface, stress);
     allSens.aggregatedVonMisesWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numThicknessGroups);
     stressvmWRTthick.setZero(); stressWeight.setZero(); allSens.aggregatedVonMisesWRTthick->setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele = 0; iele < numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenVector<double> dStressdThick(NodesPerElement);
       GenVector<double> elweight(NodesPerElement,0.0);
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele);         
       if(thgreleFlag[iele]) { 
         packedEset[iele]->getVonMisesThicknessSensitivity(dStressdThick, elweight, nodes, *elDisp, 6, surface); 
       } else {
         dStressdThick.zero();
         Vector elstress(NodesPerElement,0.0);
         packedEset[iele]->getVonMises(elstress, elweight, nodes, *elDisp, -1, surface, NULL, 0, 0, 0);
       }
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           int iparam = thpaIndex[iele];
           if(iparam >= 0) stressvmWRTthick(node, iparam) += dStressdThick[k]; 
           stressWeight(node) += elweight[k]; 
         }
       }
     } 

     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       for(int inode = 0; inode < numNodes(); ++inode)  {
         if(stressWeight(inode) == 0.0)
           stressvmWRTthick(inode, iparam) = 0.0;
         else
           stressvmWRTthick(inode, iparam) /= stressWeight(inode);
       }
     }

     if(!(*aggregatedStressDenom)) computeAggregatedStressDenom(stress);

     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       for(int inode = 0; inode < numNodes(); ++inode)  {
         (*allSens.aggregatedVonMisesWRTthick)(iparam) += stressvmWRTthick(inode,iparam)*exp(sinfo.ksParameter*(stress[inode]));
       }
       (*allSens.aggregatedVonMisesWRTthick)(iparam) /= *aggregatedStressDenom;
     }
     if(!isDynam) {
       for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
         (*allSens.aggregatedVonMisesWRTthick)(iparam) -= 
            allSens.lambdaAggregatedStressVM->adjoint()*(allSens.linearstaticWRTthick[iparam]->col(0));
       }
     }
#endif
}

void
Domain::computeAggregatedStressDenom(Vector &stress)
{
  *aggregatedStressDenom = 0.0;
  for(int k = 0; k < numNodes(); ++k) *aggregatedStressDenom += exp(sinfo.ksParameter*(stress[k]));
} 

void 
Domain::computeStressVMWRTthicknessAdjointSensitivity(int sindex,
                                                      AllSensitivities<double> &allSens,
                                                      GenVector<double> *sol,
                                                      double *bcx, 
                                                      bool *includeStressNodes,
                                                      bool isDynam)
{
#ifdef USE_EIGEN3
     // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO THICKNESS
     allSens.vonMisesWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numStressNodes, numThicknessGroups);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numStressNodes, 1);
     allSens.vonMisesWRTthick->setZero();     stressWeight.setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele=0; iele<numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       if (!includeStressNodes[iele]) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenVector<double> dStressdThick(NodesPerElement);
       GenVector<double> weight(NodesPerElement,0.0);
       int surface = senInfo[sindex].surface;
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele);    
       if(thgreleFlag[iele]) {     
         packedEset[iele]->getVonMisesThicknessSensitivity(dStressdThick, weight, nodes, *elDisp, 6, surface);  
       } else {
         dStressdThick.zero();
         Vector elstress(NodesPerElement,0.0);
         packedEset[iele]->getVonMises(elstress, weight, nodes, *elDisp, -1, surface, NULL, 0, 0, 0);
       }
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           int inode;
           if(checkIsInStressNodes(node,inode)) {
             int iparam = thpaIndex[iele];
             if(iparam >= 0) (*allSens.vonMisesWRTthick)(inode, iparam) += dStressdThick[k]; 
             stressWeight(inode, 0) += weight[k]; 
           }
         }
       }
     } 

     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
      // average vonMisesWRTthick vector
      for(int inode = 0; inode < numStressNodes; ++inode)  {
        if(stressWeight(inode, 0) == 0.0)
          (*allSens.vonMisesWRTthick)(inode, iparam) = 0.0;
        else
          (*allSens.vonMisesWRTthick)(inode, iparam) /= stressWeight(inode, 0);
      }
      if(!isDynam) {
        for(int inode = 0; inode < numStressNodes; ++inode)  {
          double a = allSens.lambdaStressVM[inode]->adjoint()*(allSens.linearstaticWRTthick[iparam]->col(0));
          (*allSens.vonMisesWRTthick)(inode, iparam) -= a;
        }
      }
     }
#endif
}

void 
Domain::computeStressVMWRTdisplacementSensitivity(int sindex,
                                                  AllSensitivities<double> &allSens, 
                                                  GenVector<double> *sol,
                                                  double *bcx)
{
#ifdef USE_EIGEN3
     allSens.vonMisesWRTdisp = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), numUncon());
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numNodes(), 1);
     allSens.vonMisesWRTdisp->setZero();     stressWeight.setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele = 0; iele < numele; iele++) { 
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       int DofsPerElement = packedEset[iele]->numDofs();
       GenFullM<double> dStressdDisp(DofsPerElement,NodesPerElement,double(0.0));
       GenVector<double> weight(NodesPerElement,0.0);
       int surface = senInfo[sindex].surface;
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele);       
       packedEset[iele]->getVonMisesDisplacementSensitivity(dStressdDisp, weight, 0, nodes, *elDisp, 6, surface, 0);
//       transformElementSensitivityInv(&dStressdDisp,iele); //TODO: watch out index of dStressdDisp when implementing
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         const auto &unconstrNum = c_dsa->getUnconstrNum();
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           auto dofs = (*allDOFs)[iele];
           stressWeight(node,0) += weight[k];
           for(int j = 0; j < DofsPerElement; ++j) {
             int dofj = unconstrNum[dofs[j]];
             if(dofs[j] < 0 || dofj < 0) continue;  // Skip undefined/constrained dofs
             if(std::isnan(dStressdDisp[j][k])) std::cerr << "nan occurs in dStressdDisp[" << j << "][" << k << "] with iele of " << iele << "\n";
             (*allSens.vonMisesWRTdisp)(node, dofj) += dStressdDisp[j][k]; 
           }
         }
       }
     }   
  
     for(int inode = 0; inode < numNodes(); ++inode)  {
       if(stressWeight(inode, 0) == 0.0)
         for(int dof = 0; dof < numUncon(); ++dof) 
           (*allSens.vonMisesWRTdisp)(inode,dof) = 0;
       else
         for(int dof = 0; dof < numUncon(); ++dof) 
           (*allSens.vonMisesWRTdisp)(inode,dof) /= stressWeight(inode,0);
     }
#endif
}

void 
Domain::computeAggregatedStressVMWRTdisplacementSensitivity(int sindex,
                                                            AllSensitivities<double> &allSens, 
                                                            GenVector<double> *sol,
                                                            double *bcx)
{
#ifdef USE_EIGEN3
     int surface = senInfo[sindex].surface;
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vonMisesWRTdisp(numNodes(), numUncon());
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numNodes(), 1);
     allSens.aggregatedVonMisesWRTdisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numUncon());
     vonMisesWRTdisp.setZero(); stressWeight.setZero(); allSens.aggregatedVonMisesWRTdisp->setZero();
     Vector stress(numNodes(),0.0);
     computeNormalizedVonMisesStress(*sol, bcx, surface, stress);
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele = 0; iele < numele; iele++) { 
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       int DofsPerElement = packedEset[iele]->numDofs();
       GenFullM<double> dStressdDisp(DofsPerElement,NodesPerElement,double(0.0));
       GenVector<double> elweight(NodesPerElement,0.0);
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele);       
       packedEset[iele]->getVonMisesDisplacementSensitivity(dStressdDisp, elweight, 0, nodes, *elDisp, 6, surface, 0);
//       transformElementSensitivityInv(&dStressdDisp,iele); //TODO: watch out index of dStressdDisp when implementing
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         const auto &unconstrNum = c_dsa->getUnconstrNum();
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           auto dofs = (*allDOFs)[iele];
           stressWeight(node,0) += elweight[k];
           for(int j = 0; j < DofsPerElement; ++j) {
             int dofj = unconstrNum[dofs[j]];
             if(dofs[j] < 0 || dofj < 0) continue;  // Skip undefined/constrained dofs
             if(std::isnan(dStressdDisp[j][k])) std::cerr << "nan occurs in dStressdDisp[" << j << "][" << k << "] with iele of " << iele << "\n";
             vonMisesWRTdisp(node, dofj) += dStressdDisp[j][k]; 
           }
         }
       }
     }   
 
     for(int inode = 0; inode < numNodes(); ++inode)  {
       if(stressWeight(inode, 0) == 0.0)
         for(int dof = 0; dof < numUncon(); ++dof) 
           vonMisesWRTdisp(inode,dof) = 0;
       else
         for(int dof = 0; dof < numUncon(); ++dof) 
           vonMisesWRTdisp(inode,dof) /= stressWeight(inode,0);
     }

     if(!(*aggregatedStressDenom)) computeAggregatedStressDenom(stress);
     for(int dof = 0; dof < numUncon(); ++dof) {
       for(int inode = 0; inode < numNodes(); ++inode)  {
         (*allSens.aggregatedVonMisesWRTdisp)(dof) += vonMisesWRTdisp(inode,dof)*exp(sinfo.ksParameter*(stress[inode]));
       }
       (*allSens.aggregatedVonMisesWRTdisp)(dof) /= *aggregatedStressDenom;
     }

#endif
}

void 
Domain::computeStressVMWRTShapeVariableDirectSensitivity(int sindex,
                                                         AllSensitivities<double> &allSens, 
                                                         GenVector<double> *sol,
                                                         double *bcx,
                                                         bool isDynam)
{
#ifdef USE_EIGEN3
     allSens.vonMisesWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), numShapeVars);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numNodes(), 1);
     allSens.vonMisesWRTshape->setZero();     stressWeight.setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele = 0; iele < numele; iele++) { 
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenFullM<double> dStressdCoord(3*NodesPerElement,NodesPerElement,double(0.0));
       GenVector<double> weight(NodesPerElement,0.0);
       int surface = senInfo[sindex].surface;
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele); 
       packedEset[iele]->getVonMisesNodalCoordinateSensitivity(dStressdCoord, weight, nodes, *elDisp, 6, surface, 0);
//       transformElementSensitivityInv(&dStressdCoord,iele); //TODO: watch out index of dStressdCoord when implementing
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         const auto &unconstrNum = c_dsa->getUnconstrNum();
         for(int k = 0; k < NodesPerElement; ++k) {
           int node1 = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           stressWeight(node1,0) += weight[k];
           for(int j = 0; j < NodesPerElement; ++j) {
             int node2 = (outFlag) ? nodeTable[(*elemToNode)[iele][j]]-1 : (*elemToNode)[iele][j];
             int inode = shapeSenData.nodes[node2];
             for(int xyz = 0; xyz < 3; ++xyz) {
               for(int iShap = 0; iShap < numShapeVars; ++iShap) {
                 (*allSens.vonMisesWRTshape)(node1, iShap) += dStressdCoord[3*j+xyz][k]*shapeSenData.sensitivities[iShap][node2][xyz]; 
               }
             }
           }
         }
       }
     } 
  
     for(int jSen = 0; jSen < numShapeVars; ++jSen) { 
       for(int inode = 0; inode < numNodes(); ++inode)  {
         if(stressWeight(inode, 0) == 0.0)
           (*allSens.vonMisesWRTshape)(inode,jSen) = 0;
         else
           (*allSens.vonMisesWRTshape)(inode,jSen) /= stressWeight(inode,0);
       }
       if(!isDynam) {
         allSens.vonMisesWRTshape->col(jSen) += *allSens.vonMisesWRTdisp * (*allSens.dispWRTshape[jSen]);
       }
     }
#endif
}

void 
Domain::computeAggregatedStressVMWRTShapeVariableSensitivity(int sindex,
                                                             AllSensitivities<double> &allSens, 
                                                             GenVector<double> *sol,
                                                             double *bcx,
                                                             bool isDynam)
{
#ifdef USE_EIGEN3
     int surface = senInfo[sindex].surface;
     allSens.aggregatedVonMisesWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numShapeVars);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vonMisesWRTshape(numNodes(), numShapeVars);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numNodes(), 1);
     vonMisesWRTshape.setZero();     stressWeight.setZero();   allSens.aggregatedVonMisesWRTshape->setZero();
     Vector stress(numNodes(),0.0);
     computeNormalizedVonMisesStress(*sol, bcx, surface, stress);
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele = 0; iele < numele; iele++) { 
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenFullM<double> dStressdCoord(3*NodesPerElement,NodesPerElement,double(0.0));
       GenVector<double> weight(NodesPerElement,0.0);
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele); 
       packedEset[iele]->getVonMisesNodalCoordinateSensitivity(dStressdCoord, weight, nodes, *elDisp, 6, surface, 0);
//       transformElementSensitivityInv(&dStressdCoord,iele); //TODO: watch out index of dStressdCoord when implementing
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         const auto &unconstrNum = c_dsa->getUnconstrNum();
         for(int k = 0; k < NodesPerElement; ++k) {
           int node1 = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           stressWeight(node1,0) += weight[k];
           for(int j = 0; j < NodesPerElement; ++j) {
             int node2 = (outFlag) ? nodeTable[(*elemToNode)[iele][j]]-1 : (*elemToNode)[iele][j];
             int inode = shapeSenData.nodes[node2];
             for(int xyz = 0; xyz < 3; ++xyz) {
               for(int iShap = 0; iShap < numShapeVars; ++iShap) {
                 vonMisesWRTshape(node1, iShap) += dStressdCoord[3*j+xyz][k]*shapeSenData.sensitivities[iShap][node2][xyz]; 
               }
             }
           }
         }
       }
     } 
  
     for(int jSen = 0; jSen < numShapeVars; ++jSen) { 
       for(int inode = 0; inode < numNodes(); ++inode) {
         if(stressWeight(inode, 0) == 0.0)
           vonMisesWRTshape(inode,jSen) = 0;
         else
           vonMisesWRTshape(inode,jSen) /= stressWeight(inode,0);
       }
     }
      
     if(!(*aggregatedStressDenom)) computeAggregatedStressDenom(stress);

     for(int jSen = 0; jSen < numShapeVars; ++jSen) { 
       for(int inode = 0; inode < numNodes(); ++inode)  {
         double a = vonMisesWRTshape(inode,jSen) * exp(sinfo.ksParameter*(stress[inode]));
         (*allSens.aggregatedVonMisesWRTshape)(jSen) += a;
       }
       (*allSens.aggregatedVonMisesWRTshape)(jSen) /= *aggregatedStressDenom;
     }

     if(!isDynam) {
       for(int jSen = 0; jSen < numShapeVars; ++jSen) 
         (*allSens.aggregatedVonMisesWRTshape)(jSen) -= 
            allSens.lambdaAggregatedStressVM->adjoint()*(allSens.linearstaticWRTshape[jSen]->col(0));
     }

#endif
}

void 
Domain::computeStressVMWRTShapeVariableAdjointSensitivity(int sindex,
                                                         AllSensitivities<double> &allSens, 
                                                         GenVector<double> *sol,
                                                         double *bcx,
                                                         bool *includeStressNodes,
                                                         bool isDynam)
{
#ifdef USE_EIGEN3
     allSens.vonMisesWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numStressNodes, numShapeVars);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numStressNodes, 1);
     allSens.vonMisesWRTshape->setZero();     stressWeight.setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele = 0; iele < numele; iele++) { 
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       if (!includeStressNodes[iele]) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenFullM<double> dStressdCoord(3*NodesPerElement,NodesPerElement,double(0.0));
       GenVector<double> weight(NodesPerElement,0.0);
       int surface = senInfo[sindex].surface;
       elDisp->zero();       
       // Determine element displacement vector
       for (int k=0; k < allDOFs->num(iele); ++k) {
         int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
         if (cn >= 0)
           (*elDisp)[k] = (*sol)[cn];
         else
           (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
       }
       transformVectorInv(*elDisp, iele); 
       packedEset[iele]->getVonMisesNodalCoordinateSensitivity(dStressdCoord, weight, nodes, *elDisp, 6, surface, 0);
//       transformElementSensitivityInv(&dStressdCoord,iele); //TODO: watch out index of dStressdCoord when implementing
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         const auto &unconstrNum = c_dsa->getUnconstrNum();
         for(int k = 0; k < NodesPerElement; ++k) {
           int node1 = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           int inode;
           if(checkIsInStressNodes(node1,inode)) {
             stressWeight(inode,0) += weight[k];
             for(int j = 0; j < NodesPerElement; ++j) {
               int node2 = (outFlag) ? nodeTable[(*elemToNode)[iele][j]]-1 : (*elemToNode)[iele][j];
               for(int xyz = 0; xyz < 3; ++xyz) {
                 for(int ishap = 0; ishap < numShapeVars; ++ishap) {
                   (*allSens.vonMisesWRTshape)(inode, ishap) += dStressdCoord[3*j+xyz][k]*shapeSenData.sensitivities[ishap][node2][xyz]; 
                 }
               }
             }
           }
         }
       }
     } 
  
     for(int jSen = 0; jSen < numShapeVars; ++jSen) { 
       for(int inode = 0; inode < numStressNodes; ++inode)  {
         if(stressWeight(inode, 0) == 0.0)
           (*allSens.vonMisesWRTshape)(inode,jSen) = 0;
         else
           (*allSens.vonMisesWRTshape)(inode,jSen) /= stressWeight(inode,0);
       }
       if(!isDynam) {
         for(int inode = 0; inode < numStressNodes; ++inode)  {
           double a = allSens.lambdaStressVM[inode]->adjoint()*(allSens.linearstaticWRTshape[jSen]->col(0));
           (*allSens.vonMisesWRTshape)(inode,jSen) -= a;
         }
       }
     }
#endif
}

void
Domain::makePreSensitivities(AllSensitivities<DComplex> &allSens, DComplex *bcx) {}

void
Domain::makePostSensitivities(GenSolver<DComplex> *sysSolver, 
                              GenSparseMatrix<DComplex> *spm,
                              AllSensitivities<DComplex> &allSens, 
                              GenVector<DComplex> *sol, DComplex *bcx,
                              GenSparseMatrix<DComplex> *K,
                              bool isDynam,
                              GeomState *refState,
                              GeomState *geomState,
                              Corotator **allCorot,
                              bool isNonLin) {}

void
Domain::makePostSensitivities(GenSolver<double> *sysSolver, 
                              GenSparseMatrix<double> *spm,
                              AllSensitivities<double> &allSens, 
                              GenVector<double> *sol, double *bcx,
                              GenSparseMatrix<double> *K,
                              bool isDynam,
                              GeomState *refState,
                              GeomState *geomState,
                              Corotator **allCorot,
                              bool isNonLin)
{
#ifdef USE_EIGEN3
 makeThicknessGroupElementFlag();
 bool *includeStressNodes = new bool[numele];
 if(solInfo().sensitivityMethod == SolverInfo::Adjoint) setIncludeStressNodes(includeStressNodes); //TODO: need to include this function in other sensitivity routine

 for(int sindex=0; sindex < numSensitivity; ++sindex) {
  switch(senInfo[sindex].type) {
   case SensitivityInfo::DisplacementWRTthickness: 
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       if(!allSens.stiffnessWRTthickSparse && !isNonLin) computeStiffnessWRTthicknessSensitivity(sindex, allSens);
       if(!allSens.linearstaticWRTthick) computeLinearStaticWRTthicknessSensitivity(sindex,allSens,sol,refState, geomState, allCorot);
       if(!isDynam) if(!allSens.dispWRTthick) computeDisplacementWRTthicknessDirectSensitivity(sindex, sysSolver, allSens, spm, K);
     }
     else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
       if(!allSens.stiffnessWRTthickSparse && !isNonLin) computeStiffnessWRTthicknessSensitivity(sindex, allSens);
       if(!allSens.linearstaticWRTthick) computeLinearStaticWRTthicknessSensitivity(sindex,allSens,sol, refState, geomState, allCorot);
       if(!isDynam) { 
         if(!allSens.lambdaDisp) { 
           if(!domain->solInfo().readInAdjointROB.empty()) {
             Rom::PodProjectionSolver* podSolver = dynamic_cast<Rom::PodProjectionSolver*>(sysSolver);
             if(podSolver) {
               std::map<OutputInfo::Type,int>::iterator it = domain->solInfo().adjointMap.find(OutputInfo::DispThic);
               if(it != domain->solInfo().adjointMap.end()) {
                 int adjointBasisId = it->second;
                 int blockCols = domain->solInfo().maxSizeAdjointBasis[adjointBasisId];
                 int startCol = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(), 
                                                domain->solInfo().maxSizeAdjointBasis.begin()+adjointBasisId, 0);
                 podSolver->setLocalBasis(startCol, blockCols);
                 podSolver->factor();
               }
               else { std::cerr << "ERROR: adjoint basis is not defined for dispthic quantity of interest\n"; }
             }
           }
           computeDisplacementDualSensitivity(sindex, sysSolver, allSens, spm, K);
         }
         computeDisplacementWRTthicknessAdjointSensitivity(sindex, allSens, spm);
         if (allSens.residual !=0 && allSens.dwrDisp==0) {
           allSens.dwrDisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numTotalDispDofs);
           for (int i=0; i<numTotalDispDofs; ++i) { 
             (*allSens.dwrDisp)[i] = allSens.lambdaDisp[i]->dot(*(allSens.residual));
           }
         }
       }
     }
     break;
   }
   case SensitivityInfo::DisplacementWRTshape: 
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       if(!allSens.stiffnessWRTshapeSparse && !isNonLin) computeStiffnessWRTShapeVariableSensitivity(sindex, allSens);
       if(!allSens.linearstaticWRTshape) computeLinearStaticWRTShapeVariableSensitivity(sindex,allSens,sol);
       if(!isDynam) if(!allSens.dispWRTshape) computeDisplacementWRTShapeVariableDirectSensitivity(sindex, sysSolver, allSens, spm, K);
     }
     else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
       if(!allSens.stiffnessWRTshapeSparse && !isNonLin) computeStiffnessWRTShapeVariableSensitivity(sindex, allSens);
       if(!allSens.linearstaticWRTshape) computeLinearStaticWRTShapeVariableSensitivity(sindex,allSens,sol);
       if(!isDynam) {
         if(!allSens.lambdaDisp) {
           if(!domain->solInfo().readInAdjointROB.empty()) {
             Rom::PodProjectionSolver* podSolver = dynamic_cast<Rom::PodProjectionSolver*>(sysSolver);
             if(podSolver) {  
               std::map<OutputInfo::Type,int>::iterator it = domain->solInfo().adjointMap.find(OutputInfo::DispThic); 
               if(it != domain->solInfo().adjointMap.end()) {
                 int adjointBasisId = it->second;
                 int blockCols = domain->solInfo().maxSizeAdjointBasis[adjointBasisId];
                 int startCol = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(), 
                                                domain->solInfo().maxSizeAdjointBasis.begin()+adjointBasisId, 0);
                 podSolver->setLocalBasis(startCol, blockCols);
                 podSolver->factor();
               }
               else { std::cerr << "ERROR: adjoint basis is not defined for dispthic quantity of interest\n"; }
             }
           }
           computeDisplacementDualSensitivity(sindex, sysSolver, allSens, spm, K);
         }
         computeDisplacementWRTShapeVariableAdjointSensitivity(sindex, allSens, spm);
         if (allSens.residual !=0 && allSens.dwrDisp==0) {
           allSens.dwrDisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numTotalDispDofs);
           for (int i=0; i<numTotalDispDofs; ++i) {
             (*allSens.dwrDisp)[i] = allSens.lambdaDisp[i]->dot(*(allSens.residual));
           }
         }
       }
     }
     break;
   }
   case SensitivityInfo::StressVMWRTthickness: 
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       if(!allSens.vonMisesWRTdisp) computeStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
       if(!allSens.stiffnessWRTthickSparse && !isNonLin) computeStiffnessWRTthicknessSensitivity(sindex, allSens);
       if(!allSens.linearstaticWRTthick) computeLinearStaticWRTthicknessSensitivity(sindex,allSens,sol, refState, geomState, allCorot);
       if(!isDynam) if(!allSens.dispWRTthick) computeDisplacementWRTthicknessDirectSensitivity(sindex, sysSolver, allSens, spm, K);
       computeStressVMWRTthicknessDirectSensitivity(sindex,allSens,sol,bcx,refState,geomState,allCorot,isDynam);
     }
     else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
       if(!allSens.vonMisesWRTdisp) computeStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
       if(!allSens.stiffnessWRTthickSparse && !isNonLin) computeStiffnessWRTthicknessSensitivity(sindex, allSens);
       if(!allSens.linearstaticWRTthick) computeLinearStaticWRTthicknessSensitivity(sindex,allSens,sol, refState, geomState, allCorot);
       if(!isDynam) {
         if(!allSens.lambdaStressVM) {
           if(!domain->solInfo().readInAdjointROB.empty()) {
             Rom::PodProjectionSolver* podSolver = dynamic_cast<Rom::PodProjectionSolver*>(sysSolver);
             if(podSolver) {
               std::map<OutputInfo::Type,int>::iterator it = domain->solInfo().adjointMap.find(OutputInfo::VMstThic);
               if(it != domain->solInfo().adjointMap.end()) {
                 int adjointBasisId = it->second;
                 int blockCols = domain->solInfo().maxSizeAdjointBasis[adjointBasisId];
                 int startCol = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(),
                                                domain->solInfo().maxSizeAdjointBasis.begin()+adjointBasisId, 0);
                 podSolver->setLocalBasis(startCol, blockCols);
                 podSolver->factor();
               }
               else { std::cerr << "ERROR: adjoint basis is not defined for vmstthic quantity of interest\n"; }
             }
           }
           computeStressVMDualSensitivity(sindex, sysSolver, allSens, spm, K);
         }
       }
       computeStressVMWRTthicknessAdjointSensitivity(sindex,allSens,sol,bcx,includeStressNodes,isDynam);
       if (allSens.residual !=0 && allSens.dwrStressVM==0) {
         allSens.dwrStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numStressNodes);
         for (int i=0; i<numStressNodes; i++) {
           (*allSens.dwrStressVM)[i] = allSens.lambdaStressVM[i]->dot(*(allSens.residual));
         }
       }
     }
     break;
   }
   case SensitivityInfo::StressVMWRTshape:
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       if(!allSens.vonMisesWRTdisp) computeStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
       if(!allSens.stiffnessWRTshapeSparse) computeStiffnessWRTShapeVariableSensitivity(sindex, allSens); 
       if(!allSens.linearstaticWRTshape) computeLinearStaticWRTShapeVariableSensitivity(sindex,allSens,sol);
       if(!isDynam) if(!allSens.dispWRTshape) computeDisplacementWRTShapeVariableDirectSensitivity(sindex, sysSolver, allSens, spm, K);
       computeStressVMWRTShapeVariableDirectSensitivity(sindex,allSens,sol,bcx,isDynam);
     }
     else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
       if(!allSens.vonMisesWRTdisp) computeStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
       if(!allSens.stiffnessWRTshapeSparse) computeStiffnessWRTShapeVariableSensitivity(sindex, allSens);
       if(!allSens.linearstaticWRTshape) computeLinearStaticWRTShapeVariableSensitivity(sindex,allSens,sol);
       if(!isDynam) {  
         if(!allSens.lambdaStressVM) {
           if(!domain->solInfo().readInAdjointROB.empty()) {
             Rom::PodProjectionSolver* podSolver = dynamic_cast<Rom::PodProjectionSolver*>(sysSolver);
             if(podSolver) {
               std::map<OutputInfo::Type,int>::iterator it = domain->solInfo().adjointMap.find(OutputInfo::VMstThic);
               if(it != domain->solInfo().adjointMap.end()) {
                 int adjointBasisId = it->second;
                 int blockCols = domain->solInfo().maxSizeAdjointBasis[adjointBasisId];
                 int startCol = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(),
                                                domain->solInfo().maxSizeAdjointBasis.begin()+adjointBasisId, 0);
                 podSolver->setLocalBasis(startCol, blockCols);
                 podSolver->factor();
               }
               else { std::cerr << "ERROR: adjoint basis is not defined for vmstthic quantity of interest\n"; }
             }
           }
           computeStressVMDualSensitivity(sindex, sysSolver, allSens, spm, K);
         }
       }
       computeStressVMWRTShapeVariableAdjointSensitivity(sindex,allSens,sol,bcx,includeStressNodes,isDynam);
       if (allSens.residual !=0 && allSens.dwrStressVM==0) {
         allSens.dwrStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numStressNodes);
         for (int i=0; i<numStressNodes; i++) {
            (*allSens.dwrStressVM)[i] = allSens.lambdaStressVM[i]->dot(*(allSens.residual));
         }
       }  
     }
     break;
   }
   case SensitivityInfo::AggregatedStressVMWRTthickness: 
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       filePrint(stderr, " ... WARNING : Only the adjoint method is available for KS function sensitivities. Switing to the adjoint method.\n");
     }
     if(aggregatedFlag) { getStressStrain(*sol, bcx, aggregatedFileNumber, AGGREGATEDVON, 0.0); aggregatedFlag = false; }
     if(!allSens.aggregatedVonMisesWRTdisp) computeAggregatedStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
     if(!allSens.stiffnessWRTthickSparse && !isNonLin) computeStiffnessWRTthicknessSensitivity(sindex, allSens);
     if(!allSens.linearstaticWRTthick) computeLinearStaticWRTthicknessSensitivity(sindex,allSens,sol, refState, geomState, allCorot);
     if(!isDynam) {
       if(!allSens.lambdaAggregatedStressVM) {
         if(!domain->solInfo().readInAdjointROB.empty()) {
           Rom::PodProjectionSolver* podSolver = dynamic_cast<Rom::PodProjectionSolver*>(sysSolver);
           if(podSolver) {
             std::map<OutputInfo::Type,int>::iterator it = domain->solInfo().adjointMap.find(OutputInfo::AGstThic);  
             if(it != domain->solInfo().adjointMap.end()) {
               int adjointBasisId = it->second;
               int blockCols = domain->solInfo().maxSizeAdjointBasis[adjointBasisId];
               int startCol = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(),
                                              domain->solInfo().maxSizeAdjointBasis.begin()+adjointBasisId, 0);
               podSolver->setLocalBasis(startCol, blockCols);
               podSolver->factor();
             }
             else { std::cerr << "ERROR: adjoint basis is not defined for agstthic quantity of interest\n"; }
           }
         }
         computeAggregatedStressVMDualSensitivity(sindex, sysSolver, allSens, spm, K);
       }
     }
     computeAggregatedStressVMWRTthicknessSensitivity(sindex,allSens,sol,bcx,isDynam);
     if (allSens.residual !=0 && allSens.dwrAggregatedStressVM == 0) {
       allSens.dwrAggregatedStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(1);
       (*allSens.dwrAggregatedStressVM)[0] = allSens.lambdaAggregatedStressVM->dot(*(allSens.residual));
     }
     break;
   }
   case SensitivityInfo::AggregatedStressVMWRTshape:
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       filePrint(stderr, " ... WARNING : Only the adjoint method is available for KS function sensitivities. Switing to the adjoint method.\n");
     }
     if(aggregatedFlag) { getStressStrain(*sol, bcx, aggregatedFileNumber, AGGREGATEDVON, 0.0); aggregatedFlag = false; }
     if(!allSens.aggregatedVonMisesWRTdisp) computeAggregatedStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
     if(!allSens.stiffnessWRTshapeSparse) computeStiffnessWRTShapeVariableSensitivity(sindex, allSens);
     if(!allSens.linearstaticWRTshape) computeLinearStaticWRTShapeVariableSensitivity(sindex,allSens,sol);
     if(!isDynam) {
       if(!allSens.lambdaAggregatedStressVM) {
         if(!domain->solInfo().readInAdjointROB.empty()) {
           Rom::PodProjectionSolver* podSolver = dynamic_cast<Rom::PodProjectionSolver*>(sysSolver);
           if(podSolver) {
             std::map<OutputInfo::Type,int>::iterator it = domain->solInfo().adjointMap.find(OutputInfo::AGstThic);
             if(it != domain->solInfo().adjointMap.end()) {
               int adjointBasisId = it->second;
               int blockCols = domain->solInfo().maxSizeAdjointBasis[adjointBasisId];
               int startCol = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(),
                                              domain->solInfo().maxSizeAdjointBasis.begin()+adjointBasisId, 0);
               podSolver->setLocalBasis(startCol, blockCols);
               podSolver->factor();
             }
             else { std::cerr << "ERROR: adjoint basis is not defined for agstthic quantity of interest\n"; }
           }
         }
         computeAggregatedStressVMDualSensitivity(sindex, sysSolver, allSens, spm, K);
       }
     }
     computeAggregatedStressVMWRTShapeVariableSensitivity(sindex,allSens,sol,bcx,isDynam);
     if (allSens.residual !=0 && allSens.dwrAggregatedStressVM == 0) {
       allSens.dwrAggregatedStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(1);
       (*allSens.dwrAggregatedStressVM)[0] = allSens.lambdaAggregatedStressVM->dot(*(allSens.residual));
    }
    break;
   }
   case SensitivityInfo::StressVMWRTmach:
   {
     if(!allSens.vonMisesWRTdisp) computeStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
     computeStressVMWRTMachNumberSensitivity(allSens);
     break;
   }
   case SensitivityInfo::StressVMWRTalpha:
   {
     if(!allSens.vonMisesWRTdisp) computeStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
     computeStressVMWRTangleOfAttackSensitivity(allSens);
     break;
   }
   case SensitivityInfo::StressVMWRTbeta:
   {
     if(!allSens.vonMisesWRTdisp) computeStressVMWRTdisplacementSensitivity(sindex,allSens,sol,bcx);
     computeStressVMWRTyawAngleSensitivity(allSens);
     break;
   }
   default:
     break;
  }
 }
 delete [] includeStressNodes;
#endif
}
