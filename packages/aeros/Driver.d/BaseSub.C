#include <Utils.d/dofset.h>
#include <Utils.d/Memory.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/IntFullM.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/CornerMaker.h>
#include <Utils.d/GlobalToLocalMap.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <algorithm>
#include <iostream>
#include <Corotational.d/utilities.h>

extern int verboseFlag;
extern int isFeti3;
extern Domain * domain;
extern int salinasFlag;

using std::list;

BaseSub::BaseSub() 
 : Domain()
{ 
}

BaseSub::BaseSub(Domain &dom, int sn, Connectivity &con, Connectivity &nds, int gn) 
 : Domain(dom, con.num(gn), con[gn].data(), nds.num(gn), nds[gn].data())
{
  subNumber = gn;
  localSubNumber = sn;
  sinfo = dom.solInfo();
  initHelm(dom);

  glNumNodes = dom.numNode();
  glElems = con[subNumber].data();
  auto subNodes = nds[subNumber];
  glNums = vec_const_int{subNodes.begin(), subNodes.end() };
  makeGlobalToLocalNodeMap(); 
  makeGlobalToLocalElemMap(); 
}

BaseSub::BaseSub(Domain &dom, int sn, int nNodes, int *nds, int nElems, int *elems, int gn) :
  Domain(dom, nElems, elems, nNodes, nds)
{
  subNumber = gn;
  localSubNumber = sn;
  sinfo = dom.solInfo();
  initHelm(dom);

  glNumNodes = nNodes;
  glNums = vec_const_int{nds, nds+nNodes};
  glElems = elems;
  makeGlobalToLocalNodeMap(); 
  makeGlobalToLocalElemMap(); 
}

BaseSub::BaseSub(Domain &dom, int sn, CoordSet* _nodes, Elemset* _elems, int *glNodeNums, int *glElemNums, int gn) :
  Domain(dom, _elems, _nodes)
{
  subNumber      = gn;
  localSubNumber = sn;  
  sinfo = dom.solInfo();
  initHelm(dom);

  glNumNodes = _nodes->size();
  glNums = vec_const_int{glNodeNums, glNodeNums+_nodes->size()};
  glElems = glElemNums;
  makeGlobalToLocalNodeMap(); 
  makeGlobalToLocalElemMap(); 
}

void
BaseSub::initHelm(Domain &dom)
{
  iRHS = 0;
  iWaveDir = 0;
//  kappa = dom.getWaveNumber();
//  omega = dom.omega;
//  fluidDensity = dom.fluidDensity;
  numWaveDirections = dom.getNumWaveDirections();
  waveDirections = (numWaveDirections) ? new double[3*numWaveDirections] : 0;
  int i;
  for(i=0;i<3*numWaveDirections;i++) {
    waveDirections[i] = dom.waveDirections[i];
  }
  implicitFlag = dom.getImplicitFlag();
  pointSourceFlag = dom.pointSourceFlag;
  sommerfeldType = dom.sommerfeldType;
  curvatureFlag = dom.curvatureFlag;
  curvatureConst1 = dom.curvatureConst1;
  curvatureConst2 = dom.curvatureConst2;
  numFFPDirections = dom.numFFPDirections;


  if(dom.solInfo().isCoupled) { // check if subdomain is mixed fluid & structure
    int nFluid = 0, nStruct = 0;
    for(int i=0; i<numele; ++i) {
      if(isFluid(i)) nFluid++;
      else nStruct++;
    }
    isMixedSub = (nFluid && nStruct) ? true : false;
	  this->isCoupled = true;
  }
  if(solInfo().getFetiInfo().dph_flag && salinasFlag) isMixedSub = true; // assume salinas doesn't know
  //if(isMixedSub) filePrint(stderr," -> sub %d is a mixed fluid/structure subdomain\n",subNumber); //HB
  if(!isMixedSub) {
    isThermalSub = packedEset[0]->getCategory() == Element::Thermal;
    isUndefinedSub = packedEset[0]->getCategory() == Element::Undefined;
    isFluidSub = isFluid(0);
  }
}


void
BaseSub::makeDSA()
{
  elemToNode = new Connectivity(packedEset.asSet()); // used in getRenumbering()
  Renumber rnum = getRenumbering();
  dsa = new DofSetArray(numnodes, packedEset, rnum.renumb);
  delete elemToNode; elemToNode = 0;
}

void
BaseSub::makeCDSA()
{
  int numdofs = dsa->size();

  int *bc = (int *)alloca(sizeof(int)*numdofs);
  //if(solInfo().isAcousticHelm()) { // FETI-DPH Acoustic
  if(implicitFlag || numComplexDirichlet) {
    if(implicitFlag) bcxC = new DComplex[numdofs * numWaveDirections];
    else bcxC = new DComplex[numdofs];
    HData::make_bc(this, bc, bcxC);      // Make bc vectors
    c_dsa = new ConstrainedDSA(*dsa, dbc, numDirichlet, cdbc, numComplexDirichlet, bc);
  }
  else {
    bcx = new double[numdofs];
    make_bc(bc, bcx);

    if(solInfo().isDynam() || solInfo().timeIntegration == 1) { // PJSA: initialize prescribed velocities
      vcx = new double[numdofs];
      acx = new double[numdofs];
      for(int i=0; i<numdofs; ++i) acx[i] = vcx[i] = 0.0;
    }

    // dof set array with boundary conditions applied
    c_dsa = new ConstrainedDSA(*dsa, numDirichlet, dbc, bc);
  }
}

FetiSubCornerHandler *
BaseSub::getCornerHandler()
{
  bool noCorners = (solInfo().solvercntl->fetiInfo.corners == FetiInfo::noCorners) ? true : false; // PJSA 4-26-06
  bool noKcw    = false;  //HB: for testing only ...
  FetiSubCornerHandler *newCH;

  // for coupled_dph remove wet interface nodes & virtual dual mpc nodes from sharedNodes list
  // so they won't be used by corner selection algorithm
  // note: once properly tested this could also be used for dual mpcs and contact in non-coupled feti-dp
  // so corner selection algorithm would not have to check for virtual nodes
  Connectivity &sharedNodes = *(scomm->sharedNodes);
  if((solInfo().isCoupled && solInfo().solvercntl->fetiInfo.fsi_corner == 0) || noCorners) { // JLchange
    int *ptr = new int[sharedNodes.csize()+1];
    int *target = new int[sharedNodes.numConnect()];
    ptr[0] = 0; int ntg = 0;
    for(int i=0; i<sharedNodes.csize(); ++i) {
      ptr[i+1] = ptr[i];
      for(int j=0; j<sharedNodes.num(i); ++j) {
        int nodeij = sharedNodes[i][j];
        bool ok = true;
        if(noKcw)
          for(int k=0; k<domain->getNodeToNode()->num(glNums[nodeij]); ++k)
            if(domain->glWetNodeMap[(*domain->getNodeToNode())[glNums[nodeij]][k]] != -1) {
              ok = false;
              break;
            }
        if(!ok || (nodeij>=numnodes) || noCorners) continue; //HB
        if((wetInterfaceNodeMap.size() > 0) && (wetInterfaceNodeMap[nodeij] == -1)) {
          ptr[i+1]++;
          target[ntg++] = nodeij;
        }
      }
    }
      drySharedNodes = new Connectivity(sharedNodes.csize(), ptr, target);
      newCH = new SubCornerHandler(subNumber, numnodes, nodes, packedEset, *nodeToNode, *dsa, *drySharedNodes,
                                   scomm->subNums,
                                   c_dsa, this);
  }
  else {
      newCH = new SubCornerHandler(subNumber, numnodes, nodes, packedEset, *nodeToNode, *dsa, *scomm->sharedNodes,
                                   scomm->subNums,
                                   c_dsa, this);
  }

  return newCH;
}

void
BaseSub::makeCCDSA()
{
 // dof set array with boundary conditions and corner node boundary conditions
 cc_dsa = std::make_unique<ConstrainedDSA>(*dsa, numDirichlet, dbc, numCRN, cornerNodes, cornerDofs,
                             numComplexDirichlet, cdbc, numWInodes, wetInterfaceNodes.data(),
                             wetInterfaceDofs.data());  // modified for coupled_dph
                                                 // (without the wet interface dofs)

 // Map from cc_dsa numbering to c_dsa numbering and vice versa
 ccToC.assign(cc_dsa->size(), -1);
 cToCC.assign(c_dsa->size(), -1);

 for(int i = 0; i < numnodes; ++i) {
   int c_dofs[DofSet::max_known_dof];
   int cc_dofs[DofSet::max_known_dof];
   int ndofs = cc_dsa->number(i, (*cc_dsa)[i], cc_dofs);
   c_dsa->number(i, (*cc_dsa)[i], c_dofs);
   for(int j = 0; j < ndofs; ++j) {
     if(cc_dofs[j] > -1) ccToC[cc_dofs[j]] = c_dofs[j];
     if(c_dofs[j] > -1) cToCC[c_dofs[j]] = cc_dofs[j];
   }
 }

  // store global corner Node numbers
  glCornerNodes.clear();
  glCornerNodes.reserve(numCRN);
  for(auto localNode: cornerNodes)
    glCornerNodes.push_back(glNums[localNode]);
}

const bool *
BaseSub::getInternalMasterFlag()
{
  if (!internalMasterFlag) {
    computeInternalMasterFlag();
  }
  return internalMasterFlag;
}

//  BClocal. Mathmatically, the size of BClocal should be
//  InterfSize * Total crnDofSize
//  However, for each column, there is only one non-zero entry.
//  Therefore,  we store an array of the size 1 * crnDofSize

IntFullM*
BaseSub::getC(int &crnDofSize, FSCommPattern<int> *sPat)
{
  // Only for FETI-2
  int i, iNode, iSub;

  // numCRN     : physical # of corner nodes
  // crnDofSize : number of lambda
  int numCRNdof = 0, numCRN = 0;
  crnDofSize = 0;
  Connectivity &sharedNodes = *(scomm->sharedNodes);
  int numNeighb = sharedNodes.csize();
  int *nodeweight = (int *) dbg_alloca(sizeof(int)*numnodes);
  int *mask       = (int *) dbg_alloca(sizeof(int)*numnodes);
  for(i=0; i< numnodes; i++) mask[i] = nodeweight[i] = 0;

  for(iSub = 0; iSub < numNeighb; ++iSub) {
     for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode)
        nodeweight[sharedNodes[iSub][iNode]]++;
  }

 if(solInfo().getFetiInfo().corners == FetiInfo::allCorners3 ||
    solInfo().getFetiInfo().corners == FetiInfo::allCorners6) {

   for(i=0; i< numnodes; i++) mask[i] = -1;
   for(iSub = 0; iSub < numNeighb; ++iSub) {
     for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode)
        mask[sharedNodes[iSub][iNode]] = iSub;
     for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
       int count = 0;
       int node = sharedNodes[iSub][iNode];
       for(i = 0; i < nodeToNode->num(node); ++i)
         if(mask[(*nodeToNode)[node][i]] == iSub) count++;
       if(count < 3) nodeweight[node] = 2; 
     }
   }
 }

  numCRN = 0; 

  for(i=0; i< numnodes; i++)  {

    if(solInfo().getFetiInfo().version == FetiInfo::feti2) {

      if ( nodeweight[i] > 1 && ((*c_dsa)[i].contains(DofSet::XYZrot) != 0)) {
        numCRN++;
        if ( c_dsa->locate(i,DofSet::Xdisp) >= 0 )  numCRNdof++;
        if ( c_dsa->locate(i,DofSet::Ydisp) >= 0 )  numCRNdof++;
        if ( c_dsa->locate(i,DofSet::Zdisp) >= 0 )  numCRNdof++;
        if(solInfo().getFetiInfo().corners == FetiInfo::noEndCorners6 ||
           solInfo().getFetiInfo().corners == FetiInfo::allCorners6) {
           if( c_dsa->locate(i,DofSet::Xrot) >= 0) numCRNdof++;
           if( c_dsa->locate(i,DofSet::Yrot) >= 0) numCRNdof++;
           if( c_dsa->locate(i,DofSet::Zrot) >= 0) numCRNdof++;
        }
      }
    } 
    else if(solInfo().getFetiInfo().version == FetiInfo::feti3) {
      if(nodeweight[i] > 1) {
        numCRN++;
        if ( c_dsa->locate(i,DofSet::Xdisp) >= 0 )  numCRNdof++;
        if ( c_dsa->locate(i,DofSet::Ydisp) >= 0 )  numCRNdof++;
        if ( c_dsa->locate(i,DofSet::Zdisp) >= 0 )  numCRNdof++;
      }
    }
  }

  if(verboseFlag)
    fprintf(stderr, " ... Sub %3d has %2d Corner Modes %d   ... \n",
            subNumber+1,numCRN,numCRNdof);

  auto dofToNeighb = scomm->getTypeSpecificList(SComm::std)->alloc_reverse();

  int numDOFs = c_dsa->size();

  int *dofCorner = (int *) dbg_alloca(sizeof(int)*numDOFs);
  int *cornerNum = (int *) dbg_alloca(sizeof(int)*numDOFs);

  for(i = 0; i < numDOFs; ++i)
    cornerNum[i] = dofCorner[i] = -1;

  auto subNums = scomm->neighbsT(SComm::std);
  int offset = 0;
  int cDof;
  for(i=0; i< numnodes; i++)  {

//#ifdef REGULAR
//     if (nodeweight[i] > 1 && ((*c_dsa)[i].contains(DofSet::XYZrot) != 0)) {
//#else
//     if (nodeweight[i] > 1) {
//#endif

// KHP:
    if( nodeweight[i] > 1 ) {
     
      if(solInfo().getFetiInfo().version == FetiInfo::feti2 && 
         ((*c_dsa)[i].contains(DofSet::XYZrot) == 0)) continue;
// KHP:

        int k;
        if ( (cDof = c_dsa->locate(i,DofSet::Xdisp)) >= 0 ) {
           int smallest = subNumber;
           for(k = 0; k < dofToNeighb->num(cDof); ++k)
            if(subNums[(*dofToNeighb)[cDof][k]] < smallest)
               smallest = subNums[(*dofToNeighb)[cDof][k]];
           dofCorner[cDof] = smallest;
         }
        if ( (cDof = c_dsa->locate(i,DofSet::Ydisp)) >= 0 ) {
           int smallest = subNumber;
           for(k = 0; k < dofToNeighb->num(cDof); ++k)
            if(subNums[(*dofToNeighb)[cDof][k]] < smallest)
               smallest = subNums[(*dofToNeighb)[cDof][k]];
           dofCorner[cDof] = smallest;
         }
        if ( (cDof = c_dsa->locate(i,DofSet::Zdisp)) >= 0 ) {
           int smallest = subNumber;
           for(k = 0; k < dofToNeighb->num(cDof); ++k)
            if(subNums[(*dofToNeighb)[cDof][k]] < smallest)
               smallest = subNums[(*dofToNeighb)[cDof][k]];
           dofCorner[cDof] = smallest;
         }
        if(solInfo().getFetiInfo().corners == FetiInfo::noEndCorners6 ||
           solInfo().getFetiInfo().corners == FetiInfo::allCorners6) {
          if ( (cDof = c_dsa->locate(i,DofSet::Xrot)) >= 0 ) {
             int smallest = subNumber;
             for(k = 0; k < dofToNeighb->num(cDof); ++k)
              if(subNums[(*dofToNeighb)[cDof][k]] < smallest)
                 smallest = subNums[(*dofToNeighb)[cDof][k]];
             dofCorner[cDof] = smallest;
           }
          if ( (cDof = c_dsa->locate(i,DofSet::Yrot)) >= 0 ) {
             int smallest = subNumber;
             for(k = 0; k < dofToNeighb->num(cDof); ++k)
              if(subNums[(*dofToNeighb)[cDof][k]] < smallest)
                 smallest = subNums[(*dofToNeighb)[cDof][k]];
             dofCorner[cDof] = smallest;
           }
          if ( (cDof = c_dsa->locate(i,DofSet::Zrot)) >= 0 ) {
             int smallest = subNumber;
             for(k = 0; k < dofToNeighb->num(cDof); ++k)
              if(subNums[(*dofToNeighb)[cDof][k]] < smallest)
                 smallest = subNums[(*dofToNeighb)[cDof][k]];
             dofCorner[cDof] = smallest;
           }
         }
     }
  }

  int j;
  for (i=0; i < scomm->numT(SComm::all); i++) {
    for (j=0; j < scomm->lenT(SComm::all,i); ++ j)
       if(dofCorner[ scomm->boundDofT(SComm::all,i,j) ] == subNumber)
               crnDofSize++;
  }

  crnPerNeighb = new int[scomm->numT(SComm::all)];
  IntFullM *BClocal = new IntFullM(4,crnDofSize);
  BClocal->zero();

  offset = 0;
  for (i=0; i < scomm->numT(SComm::all); i++) {
    crnPerNeighb[i] = 0;
    for (j=0; j < scomm->lenT(SComm::all,i); ++ j) {
      int bdof = scomm->boundDofT(SComm::all,i,j);
      if(dofCorner[bdof] == subNumber) {
        (*BClocal)[0][offset] = scomm->offsetT(SComm::all,i) + j;
        (*BClocal)[2][offset] = bdof;
        // uniquely identify a corner number;
        if(cornerNum[bdof] < 0)
          cornerNum[bdof] = offset;
        (*BClocal)[1][offset] = cornerNum[bdof];
        (*BClocal)[3][offset] = i; // Number of the neighbor
        crnPerNeighb[i]++; // Increment the number of corner with this neighb
        offset++;
      }
    }
  }
  // Sending my crnDofSize to neighbors.
  // Receive by FetiOp::getNumNeighbCRNs
  for(int i = 0; i < scomm->numT(SComm::all); ++i) {
    FSSubRecInfo<int> sInfo = sPat->getSendBuffer(subNumber, scomm->neighbT(SComm::all,i));
    sInfo.data[0] = crnDofSize;
  }

  return BClocal;
}

int
BaseSub::findProfileSize()
{
  return nodeToNode->findProfileSize(c_dsa);
}

void
BaseSub::addCornerPoints(int *glCornerList)
{
  Connectivity &sharedNodes = *(scomm->sharedNodes);
  int numShared = sharedNodes.numConnect();
  if(numShared == 0) return;

  // pick the first node on the interface (first node)
  auto allNodes = sharedNodes[0];
  int n1 = allNodes[0];
  glCornerList[glNums[n1]] = 1;
  double maxDist = 0.0, dist;
  Node &nd1 = *nodes[n1];
  int n2 = n1;
  int n;
  // find the most distant node from node 1 (second node)
  for(n=0; n<numShared; ++n) {
    Node &nd2 = *nodes[allNodes[n]];
    dist = nd2.distance(nd1);
    if(dist > maxDist) { maxDist = dist; n2 = allNodes[n]; }
  }
  glCornerList[glNums[n2]] = 1;
  Node &nd2 = *nodes[n2];
  double dx= nd2.x - nd1.x;
  double dy= nd2.y - nd1.y;
  double dz= nd2.z - nd1.z;

  // finally, find the node that maximizes the area with nodes 1 & 2
  double maxArea=0.0, area;
  int n3 = n1;
  for(n=0; n<numShared; ++n) {
    Node &nd3 = *nodes[allNodes[n]];
    double dx2= nd3.x - nd1.x;
    double dy2= nd3.y - nd1.y;
    double dz2= nd3.z - nd1.z;
    double cross[3] = {dy*dz2 - dz*dy2, dx2*dz-dz2*dx, dx*dy2 - dx2*dy};
    area = (cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
    if(area > maxArea) { maxArea = area; n3 = allNodes[n]; }
  }
  glCornerList[glNums[n3]] = 1;
}

bool
BaseSub::checkForColinearCrossPoints(int numCornerPoints,
                                     int *localCornerPoints)
{
 if(numCornerPoints < 3) return true;

 // find the most distant cross point from the first cross point.
 Node &nd1 = *nodes[localCornerPoints[0]];
 int n, n2;
 double maxDist=0.0, dist;
 for(n=1; n<numCornerPoints; ++n) {
    Node &nd2 = *nodes[localCornerPoints[n]];
    dist = nd2.distance(nd1);
    if(dist > maxDist) { maxDist = dist; n2 = localCornerPoints[n]; }
 }
 Node &nd2 = *nodes[n2];
 if(maxDist==0.0) {
   fprintf(stderr,"Coincident Cross Points! Very Bad!\n");
   fprintf(stderr,"Sub %d has %d corners\n",subNumber, numCornerPoints);
 }
 //fprintf(stderr,"Maximum Distance between cross points = %e\n",dist);

 // loop over all cross points checking for colinearity between 3 points
 double v1[3];
 double v2[3];
 double v3[3];
 v1[0] = nd2.x - nd1.x;
 v1[1] = nd2.y - nd1.y;
 v1[2] = nd2.z - nd1.z;
 double v1Norm = (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

 bool answer = true;
 for(n=1; n<numCornerPoints; ++n) {
   if(n==n2) continue;
   Node &nd3 = *nodes[localCornerPoints[n]];
   v2[0] = nd3.x - nd1.x;
   v2[1] = nd3.y - nd1.y;
   v2[2] = nd3.z - nd1.z;
   crossprod(v1,v2,v3);
   if(magnitude(v3) >= 0.01*v1Norm) {
     answer = false;
     break;
   }
 }

 return answer;
}

void
BaseSub::setCorners(gsl::span<const lc_node_idx> cNum)
{
  int i;

  // numCRN     : physical # of corner nodes
  // crnDofSize : number of degrees of freedom associated with a corner
  numCRNdof     = 0;
  numCRN        = 0;
  crnDofSize    = 0;

  int numdofs = dsa->size();
  cornerMap.assign(numdofs, -1);

  DofSet interestingDofs;
  if(solInfo().getFetiInfo().corners == FetiInfo::allCorners3 ||
     solInfo().getFetiInfo().corners == FetiInfo::noEndCorners3 ||
     solInfo().getFetiInfo().corners == FetiInfo::interface3)
    interestingDofs.mark(DofSet::XYZdisp | DofSet::Temp | DofSet::Helm | DofSet::IntPress);
  else
    interestingDofs.mark(DofSet::XYZdisp | DofSet::XYZrot | DofSet::Temp | DofSet::Helm | DofSet::IntPress);

  DofSet fluidDofs; 
  fluidDofs.mark(DofSet::Helm);
  DofSet structureDofs; 
  structureDofs.mark(DofSet::XYZdisp | DofSet::XYZrot);

  int nInterest = interestingDofs.count();
  int nFluidInterest = fluidDofs.count();
  int nStructureInterest = structureDofs.count();

  for (auto thisNode: cNum) {

    if(solInfo().isCoupled && (solInfo().solvercntl->fetiInfo.fsi_corner == 0))
      if(wetInterfaceNodeMap[thisNode] != -1) continue; // skip wet interface nodes

    if((*c_dsa)[thisNode].contains(interestingDofs.list()) != 0) {
      // no wet interface corner dofs
      if (!(solInfo().isCoupled) || (solInfo().isCoupled && (solInfo().solvercntl->fetiInfo.fsi_corner == 0))) { 
        numCRN++;
        int cdofs[9], dofs[9];//DofSet::max_known_nonL_dof
        c_dsa->number(thisNode, interestingDofs, cdofs);
        dsa->number(thisNode, interestingDofs, dofs);
        for(int iDof = 0; iDof < nInterest; ++iDof)
          if(cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
      }
      // wet interface fluid corner dofs 
      if(solInfo().isCoupled && (solInfo().solvercntl->fetiInfo.fsi_corner == 1)) { 
        if (!(onWetInterface(thisNode))) { 
          numCRN++;
          int cdofs[9], dofs[9];//DofSet::max_known_nonL_dof
          c_dsa->number(thisNode, interestingDofs, cdofs);
          dsa->number(thisNode, interestingDofs, dofs);
          for(int iDof = 0; iDof < nInterest; ++iDof) 
            if (cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
        }
        if (onWetInterfaceFluid(thisNode)) {
          numCRN++;
          int cdofs[9], dofs[9];//DofSet::max_known_nonL_dof
          c_dsa->number(thisNode, fluidDofs, cdofs);
          dsa->number(thisNode, fluidDofs, dofs);
          for(int iDof = 0; iDof < nFluidInterest; ++iDof) 
            if (cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
        }
      }
      // wet interface fluid and structure dofs 
      if(solInfo().isCoupled && (solInfo().solvercntl->fetiInfo.fsi_corner == 2)) { 
        numCRN++;
        int cdofs[9], dofs[9];//DofSet::max_known_nonL_dof
        c_dsa->number(thisNode, interestingDofs, cdofs);
        dsa->number(thisNode, interestingDofs, dofs);
        for(int iDof = 0; iDof < nInterest; ++iDof)
          if(cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
      } 
      // wet interface fluid and a few structure dofs 
      if(solInfo().isCoupled && (solInfo().solvercntl->fetiInfo.fsi_corner == 3)) { 
        if (!(onWetInterface(thisNode))) { 
          numCRN++;
          int cdofs[9], dofs[9];//DofSet::max_known_nonL_dof
          c_dsa->number(thisNode, interestingDofs, cdofs);
          dsa->number(thisNode, interestingDofs, dofs);
          for(int iDof = 0; iDof < nInterest; ++iDof) 
            if (cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
        }
        if (onWetInterfaceFluid(thisNode) || onWetInterfaceStructure(thisNode)) {
          numCRN++;
          int cdofs[9], dofs[9];//DofSet::max_known_nonL_dof
          if (onWetInterfaceFluid(thisNode) && onWetInterfaceStructure(thisNode)) {
            c_dsa->number(thisNode, interestingDofs, cdofs);
            dsa->number(thisNode, interestingDofs, dofs);
            for(int iDof = 0; iDof < nInterest; ++iDof) 
              if (cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
          } 
          else { 
            if (onWetInterfaceFluid(thisNode)) { 
              c_dsa->number(thisNode, fluidDofs, cdofs);
              dsa->number(thisNode, fluidDofs, dofs);
              for(int iDof = 0; iDof < nFluidInterest; ++iDof) 
                if (cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
            } 
            else {  
              c_dsa->number(thisNode, structureDofs, cdofs);
              dsa->number(thisNode, structureDofs, dofs);
              for(int iDof = 0; iDof < nStructureInterest; ++iDof) 
                if (cdofs[iDof] >= 0) cornerMap[dofs[iDof]] = numCRNdof++;
            } 
          } 
        }
      }
    }
  }

	// list of corner nodes and a dofset associated with it
	cornerNodes.resize(numCRN);
	cornerDofs.resize(numCRN);
	// Create a variable at the beginning to be DofSet::XYZrot, XYZrot, Temp

	isCornerNode.assign(numnodes, false);

  numCRN = 0;
  for (auto thisNode: cNum) {
    if(solInfo().isCoupled)
      if(wetInterfaceNodeMap[thisNode] != -1) continue; // skip wet interface nodes
    if((*c_dsa)[thisNode].contains(interestingDofs.list())!=0) {
      if (!(solInfo().isCoupled) || (solInfo().isCoupled && (solInfo().solvercntl->fetiInfo.fsi_corner == 0))) { 
        cornerNodes[numCRN] = thisNode;
        isCornerNode[thisNode] = true;
        cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & interestingDofs );
        numCRN++;
      } 
      if(solInfo().isCoupled && (solInfo().solvercntl->fetiInfo.fsi_corner == 1)) { 
        if (!(onWetInterface(thisNode))) { 
          cornerNodes[numCRN] = thisNode;
          isCornerNode[thisNode] = true;
          cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & interestingDofs );
          numCRN++;
        }
        if (onWetInterfaceFluid(thisNode)) {
          cornerNodes[numCRN] = thisNode;
          isCornerNode[thisNode] = true;
          cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & fluidDofs );
          numCRN++;
        }
      }
      if(solInfo().isCoupled && (getFetiInfo().fsi_corner == 2)) {
        cornerNodes[numCRN] = thisNode;
        isCornerNode[thisNode] = true;
        cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & interestingDofs );
        numCRN++;
      }
      if(solInfo().isCoupled && (getFetiInfo().fsi_corner == 3)) {
        if (!(onWetInterface(thisNode))) { 
          cornerNodes[numCRN] = thisNode;
          isCornerNode[thisNode] = true;
          cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & interestingDofs );
          numCRN++;
        }
        if (onWetInterfaceFluid(thisNode) || onWetInterfaceStructure(thisNode)) {
          cornerNodes[numCRN] = thisNode;
          isCornerNode[thisNode] = true;
          if (onWetInterfaceFluid(thisNode) && onWetInterfaceStructure(thisNode))
            cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & interestingDofs );
          else { 
            if (onWetInterfaceFluid(thisNode)) 
              cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & fluidDofs );
            else 
              cornerDofs[numCRN] = ( (*c_dsa)[thisNode] & structureDofs );
          } 
          numCRN++;
        }
      }
    }
  }
  delete drySharedNodes; drySharedNodes = 0;  // this is only needed for corner selection
}

void
BaseSub::markWetInterface(int nWI, int *wiNum, bool soweredInput)
{
  int i;
  int numnodes = dsa->numNodes();

  wetInterfaceMark.resize(numnodes);
  wetInterfaceFluidMark.resize(numnodes);
  wetInterfaceStructureMark.resize(numnodes);
  wetInterfaceCornerMark.resize(numnodes);
  for(i = 0; i < numnodes; ++i) {
    wetInterfaceMark[i] = false;
    wetInterfaceFluidMark[i] = false;
    wetInterfaceStructureMark[i] = false;
    wetInterfaceCornerMark[i] = false;
  }


  int CornerWeightOne = 1; // shared by more than 1 subdomain
  int CornerWeightTwo = 2; // shared by more than 2 subdomain

  for(i = 0; i < nWI; ++i) {
    DofSet interestingDofs;
    int thisNode = glToLocalNode[wiNum[i]];
    if(soweredInput)
        getInterestingDofs(interestingDofs, thisNode);
    else
        domain->getInterestingDofs(interestingDofs, wiNum[i]);
    if ((*c_dsa)[thisNode].contains(interestingDofs.list())) {
      wetInterfaceMark[thisNode] = true;
      // both fluid and structure wet interface corners 
      if (solInfo().isCoupled && (getFetiInfo().fsi_corner == 2) &&
          ((soweredInput && nodeToSub_sparse->num(wiNum[i]) > CornerWeightOne) ||
          (!soweredInput && nodeToSub->num(wiNum[i]) > CornerWeightOne)))
          wetInterfaceCornerMark[thisNode] = true;
      // only fluid wet interface corners 
      if (solInfo().isCoupled && (getFetiInfo().fsi_corner == 1) &&
          (*c_dsa)[thisNode].contains(DofSet::Helm)) 
        if ((soweredInput && nodeToSub_sparse->num(wiNum[i]) > CornerWeightOne)
            || (!soweredInput && nodeToSub->num(wiNum[i]) > CornerWeightOne))
        {
          wetInterfaceCornerMark[thisNode] = true;
          wetInterfaceFluidMark[thisNode] = true;
        }
      // fluid and a few structure wet interface corners 
      if (solInfo().isCoupled && (getFetiInfo().fsi_corner == 3)) {
        if ((*c_dsa)[thisNode].contains(DofSet::Helm))  
          if ((soweredInput && nodeToSub_sparse->num(wiNum[i]) > CornerWeightOne)
              || (!soweredInput && nodeToSub->num(wiNum[i]) > CornerWeightOne))
          {
            wetInterfaceCornerMark[thisNode] = true;
            wetInterfaceFluidMark[thisNode] = true;
          } 
        if ((*c_dsa)[thisNode].contains(DofSet::XYZdisp) || (*c_dsa)[thisNode].contains(DofSet::XYZrot))  
          if ((soweredInput && nodeToSub_sparse->num(wiNum[i]) > CornerWeightTwo) || (!soweredInput && nodeToSub->num(wiNum[i]) > CornerWeightTwo)) {
            wetInterfaceCornerMark[thisNode] = true;
            wetInterfaceStructureMark[thisNode] = true;
          } 
      }
    }
  }
}

#ifdef HB_DOFBASE_8
#define DOFBASE 8
#define HELMDOF 7
#else
#define DOFBASE 7
#define HELMDOF 6
#endif
//I think we could make it "dof-based" by using the fsi and just making the dof used in the fsi
//as really wet dofs. The other dofs of the wet nodes that are not used in the fsi will then be r-dofs
//(this can make a substancial difference for shell problems ...)
void
BaseSub::setWetInterface(int nWI, int *wiNum, bool soweredInput)
{
  int i;
  numWIdof   = 0;  // number of unconstrained degrees of freedom associated with wet interface
  numWInodes = 0;  // physical # of wet interface nodes in this subdomain

  int numdofs = dsa->size();
  wetInterfaceMap.assign(numdofs, -1);
  int numnodes = dsa->numNodes();
  wetInterfaceNodeMap .assign(numnodes, -1);

  // list of wet interface nodes 
  int* myWetInterfaceNodes = new int[numnodes]; //HB: numdofs is an upper bound (or use a STL vect<int>)
  numWInodes = 0;
  for(i=0; i<nWI; ++i) {
    int thisNode = glToLocalNode[wiNum[i]];
    if(solInfo().isCoupled && getFetiInfo().fsi_corner == 0) {
      DofSet interestingDofs;
      if(soweredInput) getInterestingDofs(interestingDofs, thisNode);
      else domain->getInterestingDofs(interestingDofs, wiNum[i]); // COUPLED DEBUG
      if((*c_dsa)[thisNode].contains(interestingDofs.list())!=0)
        myWetInterfaceNodes[numWInodes++] = thisNode;
    }
  }

  if ((numWInodes > 0) && solInfo().isCoupled && (getFetiInfo().fsi_corner != 0))
    std::cout << " ERROR: no wetINterface nodes should exist. " << std::endl; 

  std::sort(myWetInterfaceNodes, myWetInterfaceNodes+numWInodes);
  //small optimization of the memory 
  int* llToGlobalWImap = new int[DOFBASE*numWInodes]; //HB: 7*numWInodes is an upper bound
  wDofToNode = new int[DOFBASE*numWInodes]; //HB: 7*numWInodes is an upper bound
  for(i=0; i<DOFBASE*numWInodes; ++i) llToGlobalWImap[i] = -1;
  numWIdof = 0;
  DofSet sevenDofs;
  sevenDofs.mark(DofSet::XYZdisp | DofSet::XYZrot | DofSet::Helm);
  wetInterfaceNodes.resize(numWInodes);
  wetInterfaceDofs.resize(numWInodes);
  for(i = 0; i < numWInodes; i++) {
    int thisNode = myWetInterfaceNodes[i];
    DofSet interestingDofs;
    //interestingDofs.mark(DofSet::XYZdisp | DofSet::XYZrot | DofSet::Helm);
    domain->getInterestingDofs(interestingDofs, glNums[thisNode]); // COUPLED DEBUG
    wetInterfaceNodeMap[thisNode] = i;
    int cdofs[7], dofs[7], dmap[7];
    sevenDofs.number(interestingDofs, dmap);
    c_dsa->number(thisNode, interestingDofs, cdofs);
    dsa->number(thisNode, interestingDofs, dofs);
    for(int iDof = 0; iDof < interestingDofs.count(); ++iDof) {
      if(cdofs[iDof] >= 0) {
        llToGlobalWImap[numWIdof] = DOFBASE*localToGlobal(thisNode)+dmap[iDof];
        wetInterfaceMap[dofs[iDof]] = numWIdof;
        wDofToNode[numWIdof] = thisNode; //HB
        numWIdof++;
      }
    }
    wetInterfaceNodes[i] = thisNode;
    wetInterfaceDofs[i] = ((*c_dsa)[thisNode] & interestingDofs);
  }
  if(myWetInterfaceNodes) delete [] myWetInterfaceNodes;

  glToLocalWImap.initialize(numWIdof,llToGlobalWImap);//HB: only the first numWIdof data in the  
  delete [] llToGlobalWImap;                          //llToGlobalWImap array are meaningfull
}

int BaseSub::renumberBC(int *map)  
{
  int i;
  for(i = 0; i < numDirichlet; ++i)  {
    if (map[dbc[i].nnum] < 0)
      fprintf(stderr,"No Dirichlet mapping for %d\n", dbc[i].nnum);
    else
      dbc[i].nnum = map[dbc[i].nnum];
  }

  for(i = 0; i < numNeuman; ++i)
    nbc[i].nnum = map[nbc[i].nnum];

  for(i = 0; i < numIDis; ++i)  {
    if (map[iDis[i].nnum] < 0)  {  
      fprintf(stderr,"No IDisp mapping for %d\n", iDis[i].nnum);
      return -1;
    }
    else
      iDis[i].nnum = map[iDis[i].nnum];
  }

  for(i = 0; i < numIVel; ++i)
    iVel[i].nnum = map[iVel[i].nnum];

  return 0;
}

void BaseSub::applyAuxData()  
{
  // get attributes
  //int nAttr = geoSource->getNumAttributes();
  std::map<int, Attrib> &attributes = geoSource->getAttributes();

  // get struct props
  SPropContainer &sProps = geoSource->getStructProps();

  // get eframe data
  EFrameData *efd = geoSource->getEframes();
  int numEframes = geoSource->getNumEframes();

  // get composite layer info
  LayInfo **layInfo = geoSource->getLayInfo();
  //int numLayInfo = geoSource->getNumLayInfo();

  // get cframe data
  double **cframes = geoSource->getCframes();
  //int numCframes = geoSource->getNumCframes();
  
  // find smallest and largest cluster element number
  int minElemNum = glElems[0];
  int maxElemNum = 0;
  int iElem;
  for (iElem = 0; iElem < numele; iElem++)  {
    if (glElems[iElem] < minElemNum)
      minElemNum = glElems[iElem];

    if (glElems[iElem] > maxElemNum)
      maxElemNum = glElems[iElem];
  }

  // make loc to cluster elem map
  int *cl2LocElemMap = new int[++maxElemNum];

  for (iElem = 0; iElem < maxElemNum; iElem++)
    cl2LocElemMap[iElem] = -1;

  for (iElem = 0; iElem < numele; iElem++)
    cl2LocElemMap[ glElems[iElem] ] = iElem;

  //ResizeArray<double *> cframes = geoSource->getCframes();
  // set material properties for element
  for (int iAttr = minElemNum; iAttr < maxElemNum; iAttr++)  {

    int locElemNum = cl2LocElemMap[ attributes[iAttr].nele ];

    if (locElemNum >= 0)  { 
      //packedEset[locElemNum]->setProp(sProps+attributes[iAttr].attr);
      packedEset[locElemNum]->setProp(&sProps[attributes[iAttr].attr]);
      packedEset[locElemNum]->buildFrame(nodes);
    }
  
    if(attributes[iAttr].cmp_attr >= 0) {
/*
      if (coefData[attributes[iAttr].cmp_attr] != 0) {
        packedEset[locElemNum]->setCompositeData(1, 0, 0,
                  coefData[attributes[iAttr].cmp_attr]->values(),
                  cframes[attributes[iAttr].cmp_frm]);
      }
      else  {
*/
      if (layInfo[attributes[iAttr].cmp_attr] == 0) {
        fprintf(stderr," *** WARNING: Attribute found that refers to"
                         " nonexistant composite data: %d\n",
                           attributes[iAttr].cmp_attr+1);
        continue;
      }

      // type is 3 for LAYC 2 for LAYN
      int type = 3-layInfo[attributes[iAttr].cmp_attr]->getType();
      packedEset[locElemNum]->setCompositeData(type, layInfo[attributes[iAttr].cmp_attr]->nLayers(),
                 layInfo[attributes[iAttr].cmp_attr]->values(), 0,
                 cframes[attributes[iAttr].cmp_frm]);

    }

  }  // end of the iAttr for loop

  // set up element frames
  for (int iFrame = 0; iFrame < numEframes; iFrame++)  {

    // get local element number
    if (efd[iFrame].elnum > maxElemNum)
      continue;
    int locElemNum =  cl2LocElemMap[ efd[iFrame].elnum ];

    // check if eframe is in this subdomain
    if (locElemNum >= 0) {
      packedEset[locElemNum]->setFrame(&(efd[iFrame].frame));
      for(int k=0; k<3; ++k) {
        std::cerr << "\nprint eframe^^";
        for(int l=0; l<3; ++l)
          std::cerr << " " << efd[iFrame].frame[k][l];
      }
      std::cerr << std::endl;  
    }
  }

  // create cluster to local node map for node renumbering

  // find largest cluster number in this subdomain
  int maxClusNodeNum = 0;
  int iNode;
  for (iNode = 0; iNode < glNumNodes; iNode++)
    if ( glNums[iNode] > maxClusNodeNum )
      maxClusNodeNum = glNums[iNode];

  maxClusNodeNum++;  // for easier memory allocation

  // allocate memory for cluster to local node map
  int *cl2LocNodeMap = new int[maxClusNodeNum];

  // initialize map
  for (iNode = 0; iNode < maxClusNodeNum; iNode++)
    cl2LocNodeMap[iNode] = -1;

  // create mapping
  for (iNode = 0; iNode < glNumNodes; iNode++)
    cl2LocNodeMap[glNums[iNode]] = iNode;
 
  // renumber elements in this subdomain
  for (iElem = 0; iElem < numele; iElem++)
    packedEset[iElem]->renum(cl2LocNodeMap);


  // distribute bcs
  distributeBCs(cl2LocNodeMap);

  delete [] cl2LocElemMap;
  delete [] cl2LocElemMap;
}

void BaseSub::distributeBCs(int *cl2LocNodeMap)  {

  // get bc's from geoSource
  BCond *dbc = 0;
  BCond *nbc = 0;
  //BCond *cvbc = 0;
  BCond *iDis = 0;
  BCond *iDis6 = 0;
  BCond *iVel = 0;

  int numDirichlet = geoSource->getDirichletBC(dbc);
  int numNeuman    = geoSource->getNeumanBC(nbc);
  //int numConvBC    = geoSource->getConvBC(cvbc);
  //int numRadBC    = geoSource->getRadBC(cvbc);
  int numIDis      = geoSource->getIDis(iDis);
  int numIDis6     = geoSource->getIDis6(iDis6);
  int numIVel      = geoSource->getIVel(iVel);

  // get bc's for this subdomain
  BCond *subBC;

  int numLocDirichlet = getBC(dbc, numDirichlet, cl2LocNodeMap, subBC);
  this->setDirichlet(numLocDirichlet, subBC);

  int numLocNeuman = getBC(nbc, numNeuman, cl2LocNodeMap, subBC);
  this->setNeuman(numLocNeuman, subBC); 
/*
  int numLocConvBC = getBC(cvbc, numConvBC, cl2LocNodeMap, subBC);
  this->setConvBC(numLocConvBC, subBC);

  int numLocRadBC = getBC(cvbc, numRadBC, cl2LocNodeMap, subBC);
  this->setRadBC(numLocRadBC, subBC);
*/
  int numLocIDis = getBC(iDis, numIDis, cl2LocNodeMap, subBC);
  this->setIDis(numLocIDis, subBC);

  int numLocIDis6 = getBC(iDis6, numIDis6, cl2LocNodeMap, subBC);
  this->setIDis6(numLocIDis6, subBC);

  int numLocIVel = getBC(iVel, numIVel, cl2LocNodeMap, subBC);
  this->setIVel(numLocIVel, subBC);

}

int BaseSub::getBC(BCond *bc, int numBC, int *cl2LocNodeMap, BCond *&subBC)   {

  int numLocBC = 0;

  // check if mapping exists for bc nodes
  int iBC;
  for(iBC = 0; iBC < numBC; iBC++)
    if(cl2LocNodeMap[ bc[iBC].nnum ] >= 0)
      numLocBC++; 

  // set BC
  subBC = new BCond[numLocBC];
  
  int count = 0;
  for(iBC = 0; iBC < numBC; iBC++)  {

    if(cl2LocNodeMap[ bc[iBC].nnum ] >= 0)  {
      subBC[count] = bc[iBC];
  
      // renumber to local node number
      subBC[count].nnum = cl2LocNodeMap[ subBC[count].nnum ];
     
      // increment count
      count++;
    }
  }
  
  return numLocBC;
}

void
BaseSub::setDofPlusCommSize(FSCommStructure *pt) const
{
  for(int iSub = 0; iSub < scomm->numNeighb; ++iSub)
    pt->setLen(subNumber, scomm->subNums[iSub], scomm->sharedDOFsPlus->num(iSub));
}

void BaseSub::setControlData(ControlLawInfo *_claw, int *sensorMap,
                int *actuatorMap, int *userDispMap, int *userForceMap)  {

  claw = _claw;
  locToGlSensorMap = sensorMap;
  locToGlActuatorMap = actuatorMap;
  locToGlUserDispMap = userDispMap;
  locToGlUserForceMap = userForceMap;

}

#ifdef DISTRIBUTED
void BaseSub::setOutputNodes(int nNodes, int *nodeList, int *outList)  {

  numNodalOutput = nNodes;
  outputNodes = nodeList;
  outIndex = outList;
}
#endif

int BaseSub::countElemNodes()  {

  int numElemNodes = 0;
  for (int iElem = 0; iElem < numele; iElem++)
    numElemNodes += packedEset[iElem]->numNodes();

  return numElemNodes;
}

void BaseSub::addNodeXYZ(double *centroid, double* nNodes)
{
  for(int i=0; i<nodes.size(); ++i) {
    Node &nd = nodes.getNode(i);
    nNodes[group] += 1.0;
    centroid[3*group] += nd.x;
    centroid[3*group+1] += nd.y;
    centroid[3*group+2] += nd.z;
  }
}

double BaseSub::getSharedDofCount() 
{
  // PJSA 7-29-03: this function is used to compute the total number of DOFs
  // in the domain (required only for timing file)
  int i;
  double count = 0.0;
  for(i = 0; i < totalInterfSize; ++i) {
    int cc_dof = allBoundDofs[i];
    if(cc_dof >= 0) { // don't count mpc or contact virtual nodes
      count += 1.0/dofWeight(cc_dof);
    }
  }
  return count;
}

int BaseSub::getTotalDofCount()
{ 
  // PJSA 7-29-03: this function is used to compute the total number of DOFs
  // in the domain (required only for timing file)
  return dsa->size();
}

int
BaseSub::globalNumNodes()
{
  return geoSource->getNumGlobNodes();
}

BaseSub::~BaseSub()
{
  if(crnPerNeighb) { delete [] crnPerNeighb; crnPerNeighb = 0; }

  if(dualToBoundary) { delete [] dualToBoundary; dualToBoundary = 0; }
  if(neighbNumGRBMs) { delete [] neighbNumGRBMs; neighbNumGRBMs = 0; }
  if(bcx) { delete [] bcx; bcx = 0; }
  if(bcxC) { delete [] bcxC; bcxC = 0; }
  if(vcx) { delete [] vcx; vcx = 0; }
  if(acx) { delete [] acx; acx = 0; }
  if(faceIsSafe) { delete [] faceIsSafe; faceIsSafe = 0; }

  if(neighbK_p) { delete [] neighbK_p; neighbK_p = 0; }
  if(neighbK_s) { delete [] neighbK_s; neighbK_s = 0; }
  if(neighbK_s2) { delete [] neighbK_s2; neighbK_s2 = 0; }
  if(neighbK_f) { delete [] neighbK_f; neighbK_f = 0; }
  if(neighbYmod) { delete [] neighbYmod; neighbYmod = 0; } 
  if(neighbPrat) { delete [] neighbPrat; neighbPrat = 0; }
  if(neighbDens) { delete [] neighbDens; neighbDens = 0; }
  if(neighbThih) { delete [] neighbThih; neighbThih = 0; }
  if(neighbSspe) { delete [] neighbSspe; neighbSspe = 0; }

  if(localToGlobalMPC) { delete [] localToGlobalMPC; localToGlobalMPC = 0; }
  if(localToGlobalFSI) { delete [] localToGlobalFSI; localToGlobalFSI = 0; }
  if(localToGroupMPC) { delete [] localToGroupMPC; localToGroupMPC = 0; }
  if(mpcToDof) { delete mpcToDof; mpcToDof = 0; }
  if(localMpcToMpc) { delete localMpcToMpc; localMpcToMpc = 0; }
  if(localMpcToGlobalMpc) { delete localMpcToGlobalMpc; localMpcToGlobalMpc = 0; }
  if(localMpcToBlock) { delete localMpcToBlock; localMpcToBlock = 0; }
  if(blockToLocalMpc) { delete blockToLocalMpc; blockToLocalMpc = 0; }
  if(localMpcToBlockMpc) { delete localMpcToBlockMpc; localMpcToBlockMpc = 0; }
  if(mpcToBoundDof) { delete mpcToBoundDof; mpcToBoundDof = 0; }
  if(mpcMaster) { delete [] mpcMaster; mpcMaster = 0; } 
  if(neighbNumGroupGrbm) { delete [] neighbNumGroupGrbm; neighbNumGroupGrbm = 0; }
  if(neighbGroupGrbmOffset) { delete [] neighbGroupGrbmOffset; neighbGroupGrbmOffset = 0; }

  if(localLambda) { delete [] localLambda; localLambda = 0; }

  if(neighbGlToLocalWImap) { delete [] neighbGlToLocalWImap; neighbGlToLocalWImap = 0; }
  /*if(wweight) { delete [] wweight; wweight = 0; }*/
  if(drySharedNodes) { delete drySharedNodes; drySharedNodes = 0; }
  if(mpclast) { delete [] mpclast; mpclast = 0; }

  if(scomm) { delete scomm; scomm = 0; }
  if(internalMasterFlag) { delete [] internalMasterFlag; internalMasterFlag = 0; }
  // don't delete subToSub_fsi;
  // don't delete fsiNeighb;
//  if(glToLocalElem) delete [] glToLocalElem;
#ifdef DISTRIBUTED
  if(outputNodes) delete [] outputNodes;
  if(outIndex) delete [] outIndex;
#endif
  if(claw) {
    if(claw->actuator) { delete [] claw->actuator; claw->actuator = 0; }
    if(claw->userForce) { delete [] claw->userForce; claw->userForce = 0; }
    if(claw->userDisp) { delete [] claw->userDisp; claw->userDisp = 0; }
    delete claw; claw = 0;
  }
  if(locToGlSensorMap) { delete [] locToGlSensorMap; locToGlSensorMap = 0; }
  if(locToGlActuatorMap) { delete [] locToGlActuatorMap; locToGlActuatorMap = 0; }
  if(locToGlUserDispMap) { delete [] locToGlUserDispMap; locToGlUserDispMap = 0; }
  if(locToGlUserForceMap) { delete [] locToGlUserForceMap; locToGlUserForceMap = 0; }
}

void
BaseSub::setWaveNumbers(double *waveNumbers)
{
  // this is used by salinas, not fem
  k_p = waveNumbers[0];
  k_s = waveNumbers[1];
  k_s2 = waveNumbers[2];
  k_f = waveNumbers[3];
}

void
BaseSub::computeWaveNumbers()
{
  double omega2 = geoSource->shiftVal();
  double E, nu, rho, eh, di, beta4, kappaHelm;

  if(getFetiInfo().waveMethod == FetiInfo::uniform) {
    // this should only be used when there is one material property
    StructProp *sProps = packedEset[0]->getProperty();
    E = (*sProps).E;
    nu = (*sProps).nu;
    rho = (*sProps).rho;
    eh  = (*sProps).eh;
    kappaHelm = (*sProps).kappaHelm;
    double lambda = (nu*E)/(1+nu)/(1-2*nu);
    double mu = E/2.0/(1+nu);
    if((solInfo().getFetiInfo().waveType == FetiInfo::fluid) ||
       ((solInfo().getFetiInfo().waveType == FetiInfo::any) && packedEset[0]->isFluidElement())) {
      k_f = kappaHelm;
    }
    else if((solInfo().getFetiInfo().waveType == FetiInfo::solid) || 
       ((solInfo().getFetiInfo().waveType == FetiInfo::any) && !packedEset[0]->isShell())) {
      k_p = sqrt(omega2 * rho)/sqrt(lambda + 2*mu);
      k_s = sqrt(omega2 * rho)/sqrt(mu);
      k_s2 = k_s; // temp fix for isotropic subdomains
    }
    else if((solInfo().getFetiInfo().waveType == FetiInfo::shell) ||
            ((solInfo().getFetiInfo().waveType == FetiInfo::any) && packedEset[0]->isShell())) {
      di = E*eh*eh*eh/(12.0*(1-nu*nu));
      beta4 = omega2*rho*eh/di;
      k_p = sqrt(sqrt(beta4));
      k_s = k_p;
      k_s2 = k_s; // temp fix for isotropic subdomains
    }
  }
  else if(getFetiInfo().waveMethod == FetiInfo::averageK) {
    k_f = k_p = k_s = k_s2 = 0.0;
    int fcount = 0, scount = 0;
    for(int i=0; i<numele; ++i) {
      StructProp *sProps = packedEset[i]->getProperty();
      if(!sProps) continue; // phantom
      E = (*sProps).E;
      nu = (*sProps).nu;
      rho = (*sProps).rho;
      eh  = (*sProps).eh;
      kappaHelm = (*sProps).kappaHelm;
      double lambda = (nu*E)/(1+nu)/(1-2*nu);
      double mu = E/2.0/(1+nu);
      if((solInfo().getFetiInfo().waveType == FetiInfo::fluid) ||
         ((solInfo().getFetiInfo().waveType == FetiInfo::any) && packedEset[0]->isFluidElement())) { // fluid
        k_f += kappaHelm;
        fcount++;
      }
      else if((solInfo().getFetiInfo().waveType == FetiInfo::shell) ||
              ((solInfo().getFetiInfo().waveType == FetiInfo::any) && packedEset[i]->isShell())) { // shell
        di = E*eh*eh*eh/(12.0*(1-nu*nu));
        beta4 = omega2*rho*eh/di;
        double elek_p = sqrt(sqrt(beta4));
        k_p += elek_p;
        k_s += elek_p;
        k_s2 += elek_p; // temp fix for isotropic subdomains
        scount++;
      }
      else if((solInfo().getFetiInfo().waveType == FetiInfo::solid) ||
         ((solInfo().getFetiInfo().waveType == FetiInfo::any) && !packedEset[i]->isConstraintElement() && !packedEset[i]->isSommerElement())) { // solid
        double elek_p = sqrt(omega2 * rho)/sqrt(lambda + 2*mu);
        k_p += elek_p;
        double elek_s = sqrt(omega2 * rho)/sqrt(mu);
        k_s += elek_s;
        k_s2 += elek_s; // temp fix for isotropic subdomains
        scount++;
      }
    }
    if(fcount > 0)
      k_f /= fcount;
    if(scount > 0) {
      k_p /= scount;
      k_s /= scount;
      k_s2 /= scount;
    } 
  }
}

void
BaseSub::averageMatProps()
{
  double sumE = 0.0, sumNu = 0.0, sumRho = 0.0, sumEh = 0.0, sumSs = 0.0;
  int numShellElems = 0;
  int scount = 0, fcount = 0;
  for(int i=0; i<numele; ++i) {
    //if(isFluid(i)) continue;
    StructProp *sProps = packedEset[i]->getProperty();
    if(!sProps) continue; // should be null for phantom
    if(packedEset[i]->isFluidElement()) {
      //cerr << "averageMat not supported for fluid and coupled problems\n"; exit(-1); // soundSpeed is not always set under FLUMAT
      sumSs += real((*sProps).soundSpeed);
      fcount++;
    }
    else if(!packedEset[i]->isConstraintElement() && !packedEset[i]->isSommerElement()) {
      sumE += (*sProps).E;
      sumNu += (*sProps).nu;
      sumRho += (*sProps).rho;
      if(packedEset[i]->isShell()) {
        numShellElems++;
        sumEh  += (*sProps).eh;
      }
      scount++;
    }
  }
  if(fcount > 0) {
    Sspe = sumSs/fcount;
  }
  if(scount > 0) {
    Ymod = sumE/scount;
    Prat = sumNu/scount;
    Dens = sumRho/scount;
    if(numShellElems > 0) Thih = sumEh/numShellElems;
  }
   
}

// *****************************************************************************
// the following functions are for the new sower binary distributed input

void 
BaseSub::setDirichletBC(list<BCond *> *_list)
{
  numDirichlet = 0;
  dbc = new BCond[_list->size()];
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    dbc[numDirichlet++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setNeumanBC(list<BCond *> *_list)
{
  numNeuman = 0;
  nbc = new BCond[_list->size()];
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    nbc[numNeuman++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setInitialDisplacement(list<BCond *> *_list)
{
  numIDis = 0;
  iDis = new BCond[_list->size()];
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    iDis[numIDis++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setInitialDisplacement6(list<BCond *> *_list)
{
  numIDis = 0;
  iDis6 = new BCond[_list->size()];
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    iDis6[numIDis6++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setInitialVelocity(list<BCond *> *_list)
{
  numIVel = 0;
  iVel = new BCond[_list->size()];
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    iVel[numIVel++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setSensor(list<BCond *> *_list)
{
  claw->numSensor = 0;
  claw->sensor = new BCond[_list->size()];
  locToGlSensorMap = new int[_list->size()];
  int temp=0;
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    locToGlSensorMap[temp++] = (int)(*it)->val;//the value of the map is store in val
    (*it)->val = 0.0;
    claw->sensor[claw->numSensor++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setActuator(list<BCond *> *_list)
{
  claw->numActuator = 0;
  claw->actuator = new BCond[_list->size()];
  locToGlActuatorMap = new int[_list->size()];
  int temp=0;
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    locToGlActuatorMap[temp++] = (int)(*it)->val;//the value of the map is store in val
    (*it)->val = 0.0;
    claw->actuator[claw->numActuator++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setUsdd(list<BCond *> *_list)
{
  claw->numUserDisp = 0;
  claw->userDisp = new BCond[_list->size()];
  locToGlUserDispMap = new int[_list->size()];
  int temp=0;
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    locToGlUserDispMap[temp++] = (int)(*it)->val;//the value of the map is store in val
    (*it)->val = 0.0;
    claw->userDisp[claw->numUserDisp++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setUsdf(list<BCond *> *_list)
{
  claw->numUserForce = 0;
  claw->userForce = new BCond[_list->size()];
  locToGlUserForceMap = new int[_list->size()];
  int temp = 0;
  for(list<BCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    locToGlUserForceMap[temp++] = (int)(*it)->val;//the value of the map is store in val
    (*it)->val = 0.0;
    claw->userForce[claw->numUserForce++].setData((*it)->nnum, (*it)->dofnum, (*it)->val);
  }
}

void
BaseSub::setComplexDirichletBC(list<ComplexBCond *> *_list)
{
  numComplexDirichlet = 0;
  cdbc = new ComplexBCond[_list->size()];
  for(list<ComplexBCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    cdbc[numComplexDirichlet++].setData((*it)->nnum, (*it)->dofnum, (*it)->reval, (*it)->imval);
  }
}

void
BaseSub::setComplexNeumanBC(list<ComplexBCond *> *_list)
{
  numComplexNeuman = 0;
  cnbc = new ComplexBCond[_list->size()];
  for(list<ComplexBCond *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    cnbc[numComplexNeuman++].setData((*it)->nnum, (*it)->dofnum, (*it)->reval, (*it)->imval);
  }
}

void
BaseSub::setDnb(list<SommerElement *> *_list)
{
  numNeum = 0;
  for(list<SommerElement *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    addNeum((*it));
  }
}

void
BaseSub::setScat(list<SommerElement *> *_list)
{
  numScatter = 0;
  for(list<SommerElement *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    addScatter((*it));
  }
}

void
BaseSub::setArb(list<SommerElement *> *_list)
{
  numSommer = 0;
  for(list<SommerElement *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    addSommer((*it));
    //add this element in the packedEset
    packedEset.elemadd(packedEset.last(),(*it));
    numele++;
    (*it)->dom = this;
  }
}

void
BaseSub::setWet(list<SommerElement *> *_list)
{
  numWet = 0;
  for(list<SommerElement *>::iterator it = _list->begin(); it != _list->end(); ++it) {
    addWet((*it));
  }
}

void
BaseSub::addFsiElements()
{
  localToGlobalFSI = new int[numFSI];
  for(int i = 0; i < numFSI; ++i) {
    int fluidNode = fsi[i]->lmpcnum;
    localToGlobalFSI[i] = glNums[fluidNode];
    for(int j = 0; j < fsi[i]->nterms; j++) {
      LMPCons *thisFsi = new LMPCons(fluidNode, 0.0);
      LMPCTerm thisLmpcTerm(fsi[i]->terms[j], 1.0);
      thisFsi->addterm(&(thisLmpcTerm));
      addSingleFsi(thisFsi);
    }
  }
}

void
BaseSub::addSingleFsi(LMPCons *localFsi)
{
#ifndef SALINAS
  int nEle = packedEset.last();
  packedEset.fsielemadd(nEle, localFsi);
  numele = packedEset.last();
#endif
}

int
BaseSub::isFluid(int i) 
{
  HelmElement *he = dynamic_cast<HelmElement *>(packedEset[i]);
  if(he==0) return 0; // non-Helmholtz element found
  else return he->isFluid();
}

void
BaseSub::makeGlobalToLocalNodeMap()
{
 globalNMax = 0;
 int i;
 for(i = 0; i < numnodes; ++i) 
   if(glNums[i] > globalNMax) globalNMax = glNums[i];
 globalNMax += 1;

// RT: 030813
// glToLocalNode = new int[globalNMax+1];
// for(i = 0; i <= globalNMax; ++i) glToLocalNode[i] = -1;
// for(i = 0; i < numnodes; ++i) glToLocalNode[glNums[i]] = i;
 glToLocalNode.initialize(numnodes,glNums.data());
}

void
BaseSub::makeGlobalToLocalElemMap()
{
 globalEMax = 0;
 int i;
 int nEls = packedEset.last();
 for(i = 0; i < nEls; ++i)
   if(glElems[i] > globalEMax) globalEMax = glElems[i];
 globalEMax += 1;

// RT: 030813
// glToLocalElem = new int[globalEMax+1];
// for(i = 0; i <= globalEMax; ++i) glToLocalElem[i] = -1;
// for(i = 0; i < packedEset.last(); ++i) glToLocalElem[glElems[i]] = i;
 glToLocalElem.initialize(packedEset.last(),glElems);
}


void
BaseSub::makeMpcInterface(Connectivity *subToMpc, const Connectivity &lmpcToSub,
                          Connectivity *subToSub_mpc)
{
    int i,j,k;

    // Step 0: figure out the Mpc neigbors
    int numMpcNeighb = subToSub_mpc->num(subNumber);
    std::vector<gl_sub_idx> mpcNeighb((*subToSub_mpc)[subNumber].begin() ,
                                      (*subToSub_mpc)[subNumber].begin() + numMpcNeighb);
    for(i = 0; i < numMpcNeighb; ++i)  mpcNeighb[i] = (*subToSub_mpc)[subNumber][i];
    int *mpcNeighbSize = new int[numMpcNeighb];
    for(i = 0; i < numMpcNeighb; ++i) mpcNeighbSize[i] = 0;
    int numtarget = 0;
    for(i = 0; i < subToMpc->num(subNumber); ++i) {
        int iMPC = (*subToMpc)[subNumber][i];
        for(j = 0; j < lmpcToSub.num(iMPC); ++j) {
            int jSub = lmpcToSub[iMPC][j];
            if( (jSub != subNumber) || (lmpcToSub.num(iMPC) == 1) ) {
                int jMpcNeighb = subToSub_mpc->cOffset(subNumber, jSub);
                mpcNeighbSize[jMpcNeighb] += 1;
                numtarget++;
            }
        }
    }

    // Step 1: build the mpcNeighbToMpc connectivity
    int *pointer = new int[numMpcNeighb+1]; pointer[0] = 0;
    for(i=0; i<numMpcNeighb; ++i) pointer[i+1] = pointer[i] + mpcNeighbSize[i];
    int *target = new int[numtarget];
    for(i = 0; i < numMpcNeighb; ++i) mpcNeighbSize[i] = 0;
    for(i = 0; i < subToMpc->num(subNumber); ++i) {
        int iMPC = (*subToMpc)[subNumber][i];
        int locMpc = globalToLocalMPC[iMPC];
        for(j = 0; j < lmpcToSub.num(iMPC); ++j) {
            int jSub = lmpcToSub[iMPC][j];
            if( (jSub != subNumber) || (lmpcToSub.num(iMPC) == 1) ) {
                int jMpcNeighb = subToSub_mpc->cOffset(subNumber, jSub);
                int t = pointer[jMpcNeighb] + mpcNeighbSize[jMpcNeighb];
                target[t] = locMpc;
                mpcNeighbSize[jMpcNeighb] += 1;
            }
        }
    }
    Connectivity *mpcNeighbToMpc = new Connectivity(numMpcNeighb, pointer, target);

    std::vector<size_t> ptr(numMpcNeighb+1);
    ptr[0] = 0;
    for(i = 0; i < numMpcNeighb; ++i)
        ptr[i+1] = ptr[i] + mpcNeighbSize[i];
    std::vector<int> tgt;
    tgt.reserve( ptr[numMpcNeighb] );
    for(j = 0; j < numMpcNeighb; ++j) {
        for(k = 0; k < mpcNeighbSize[j]; ++k)
            tgt.push_back((*mpcNeighbToMpc)[j][k]);
    }
    auto mpcInterfaceDOFs = std::make_unique<Connectivity>( numMpcNeighb, std::move(ptr), std::move(tgt) );
    scomm->setTypeSpecificList( SComm::mpc, std::move(mpcNeighb), std::move(mpcInterfaceDOFs) );

    delete [] mpcNeighbSize;
    delete mpcNeighbToMpc;
}

void
BaseSub::makeFsiInterface(const Connectivity *subToFsi, const Connectivity &fsiToSub,
                          const Connectivity *subToSub_fsi)
{
    int i,j,k;

    // Step 0: figure out the Fsi neigbors
    numFsiNeighb = subToSub_fsi->num(subNumber);
    fsiNeighb = new int[numFsiNeighb];
    for(i = 0; i < numFsiNeighb; ++i)  fsiNeighb[i] = (*subToSub_fsi)[subNumber][i];
    int *fsiNeighbSize = new int[numFsiNeighb];
    for(i = 0; i < numFsiNeighb; ++i) fsiNeighbSize[i] = 0;
    int numtarget = 0;
    for(i = 0; i < subToFsi->num(subNumber); ++i) {
        int glFsi = (*subToFsi)[subNumber][i];
        for(j = 0; j < fsiToSub.num(glFsi); ++j) {
            int jSub = fsiToSub[glFsi][j];
            if( (jSub != subNumber) || (fsiToSub.num(glFsi) == 1) ) {
                int jFsiNeighb = subToSub_fsi->cOffset(subNumber, jSub);
                fsiNeighbSize[jFsiNeighb] += 1;
                numtarget++;
            }
        }
    }

    // Step 1: build the fsiNeighbToFsi connectivity
    int *pointer = new int[numFsiNeighb+1]; pointer[0] = 0;
    for(i=0; i<numFsiNeighb; ++i) pointer[i+1] = pointer[i] + fsiNeighbSize[i];
    int *target = new int[numtarget];
    for(i = 0; i < numFsiNeighb; ++i) fsiNeighbSize[i] = 0;
    for(i = 0; i < subToFsi->num(subNumber); ++i) {
        int glFsi = (*subToFsi)[subNumber][i];
        for(j = 0; j < fsiToSub.num(glFsi); ++j) {
            int jSub = fsiToSub[glFsi][j];
            if( (jSub != subNumber) || (fsiToSub.num(glFsi) == 1) ) {
                int jFsiNeighb = subToSub_fsi->cOffset(subNumber, jSub);
                int t = pointer[jFsiNeighb] + fsiNeighbSize[jFsiNeighb];
                target[t] = glFsi;
                fsiNeighbSize[jFsiNeighb] += 1;
            }
        }
    }
    Connectivity fsiNeighbToFsi(numFsiNeighb, pointer, target);

    std::vector<size_t> ptr(numFsiNeighb+1);
    ptr[0] = 0;
    for(i = 0; i < numFsiNeighb; ++i)
    	ptr[i+1] = ptr[i] + fsiNeighbSize[i];
    std::vector<int> tgt;
    tgt.reserve( ptr[numFsiNeighb] );
    for(j = 0; j < numFsiNeighb; ++j) {
        for(k = 0; k < fsiNeighbSize[j]; ++k)
            tgt.push_back( fsiNeighbToFsi[j][k] );
    }
    auto fsiInterfaceDOFs = std::make_unique<Connectivity>(numFsiNeighb, std::move(ptr), std::move(tgt) );
    scomm->setTypeSpecificList(SComm::fsi,
                               std::vector<gl_sub_idx>(fsiNeighb, fsiNeighb+numFsiNeighb),
                               std::move(fsiInterfaceDOFs) );

    delete [] fsiNeighbSize;
}

