#include <iostream>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <Utils.d/dbg_alloca.h>

#include <map>
#include <list>
using std::map;
using std::list;

#include <unistd.h>
#include <Timers.d/GetTime.h>
#include <Threads.d/PHelper.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>

#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/Domain.h>
#include <Element.d/Element.h>
#include <Utils.d/ModeData.h>
#include <Driver.d/GeoSource.h>
#include <Feti.d/DistrVector.h>
#include <Corotational.d/DistrGeomState.h>
#include <Hetero.d/FlExchange.h>

#include <Element.d/Rigid.d/RigidBeam.h>
#include <Element.d/Rigid.d/RigidThreeNodeShell.h>
#include <Element.d/Rigid.d/RigidFourNodeShell.h>
#include <Element.d/Rigid.d/RigidSolid6Dof.h>

#include "Rom.d/BasisBinaryFile.h"
#include "Rom.d/XPostOutputFile.h"
#include <Rom.d/RobCodec.h>
#include <Sfem.d/Sfem.h>

extern int verboseFlag;
extern Sfem *sfem;
extern GeoSource *geoSource;
extern int totalNewtonIter;

// Global variable for mode data
ModeData modeData;
ModeData modeDataIVel;
ModeData modeDataIDis;
ModeData modeDataMode;

//----------------------------------------------------------------------------------

Domain::Domain(Domain &d, int nele, const int *eles, int nnodes, const int *nnums)
  : nodes(*new CoordSet(nnodes)), lmpc(0), fsi(0), ymtt(0), ctett(0), sdetaft(0),
#ifdef USE_EIGEN3
    rubdaft(0),
#endif
    ss1dt(0), ss2dt(0), ysst(0), yssrt(0), ymst(0), SurfEntities(0), MortarConds(0)
{
 initialize();

 int iele;
 numele = nele; // number of elements
 for(int i=0; i < numele; ++i)
   packedEset.elemadd(i, d.packedEset[eles[i]]);

 numnodes = nnodes; // number of nodes
 for(int i=0; i < numnodes; ++i)
   if(d.nodes[nnums[i]] != NULL)
     nodes.nodeadd(i, *d.nodes[nnums[i]]);

 if(d.gravityAcceleration) {
   gravityAcceleration = new double [3];
   gravityAcceleration[0] = d.gravityAcceleration[0];
   gravityAcceleration[1] = d.gravityAcceleration[1];
   gravityAcceleration[2] = d.gravityAcceleration[2];
 }

 mftval = d.mftval;
 hftval = d.hftval;

 if(verboseFlag == 0) setSilent();
 else setVerbose();

 matrixTimers = new MatrixTimers;
 initializeNumbers();
}

Domain::Domain(Domain &d, Elemset *_elems, CoordSet *_nodes)
  : nodes(*_nodes), lmpc(0), fsi(0), ymtt(0), ctett(0), sdetaft(0),
#ifdef USE_EIGEN3
    rubdaft(0),
#endif
    ss1dt(0), ss2dt(0), ysst(0), yssrt(0), ymst(0), SurfEntities(0), MortarConds(0)
{
 initialize();

 int iele;
 numele = _elems->last();
 for(iele=0; iele < numele; ++iele)
   packedEset.elemadd(iele, (*_elems)[iele]);

 numnodes = _nodes->size();

 if(d.gravityAcceleration) {
   gravityAcceleration = new double [3];
   gravityAcceleration[0] = d.gravityAcceleration[0];
   gravityAcceleration[1] = d.gravityAcceleration[1];
   gravityAcceleration[2] = d.gravityAcceleration[2];
 }

 mftval = d.mftval;
 hftval = d.hftval;

 if(verboseFlag == 0) setSilent();
 else setVerbose();
 initializeNumbers();
 matrixTimers = new MatrixTimers;
}

Domain::Domain(int iniSize) : nodes(*(new CoordSet(iniSize*16))), packedEset(iniSize*16), lmpc(0,iniSize),
   fsi(0,iniSize), ymtt(0,iniSize), ctett(0,iniSize), sdetaft(0,iniSize),
#ifdef USE_EIGEN3
   rubdaft(0,iniSize),
#endif
   ss1dt(0, iniSize), ss2dt(0, iniSize), ysst(0,iniSize), yssrt(0,iniSize), ymst(0,iniSize),
   SurfEntities(0,iniSize), MortarConds(0,iniSize)
{
 initialize();

 if(verboseFlag == 0) setSilent();
 else setVerbose();

 matrixTimers = new MatrixTimers;
 initializeNumbers();
}

void
Domain::makeAllDOFs() // build the dof connectivity
{
 int numele = packedEset.last(); // PJSA 5-2-05: include phantoms here

 int iele;
 int *pointers = new int[numele+1];
 pointers[0] = 0;
 maxNumDOFs  = 0;
 for(iele=0; iele < numele; ++iele) {
   int numDOFs = packedEset[iele]->numDofs();
   if(numDOFs > maxNumDOFs) maxNumDOFs = numDOFs;
   pointers[iele+1] = pointers[iele] + numDOFs;
 }
 int *targets = new int[ pointers[numele] ];
 for(iele=0; iele < numele; ++iele)
   packedEset[iele]->dofs(*dsa, targets + pointers[iele]);
 if(allDOFs) delete allDOFs;
 allDOFs = new Connectivity(numele, pointers, targets);

 maxNumNodes = 0;
 // compute maximum number of nodes per element
 for(iele=0; iele < numele; ++iele) {
   int numNodesPerElement = packedEset[iele]->numNodes();
   maxNumNodes = std::max(maxNumNodes, numNodesPerElement);
 }
}

// build the dof connectivity for fluid
void
Domain::makeAllDOFsFluid()
{
 // Test if allDOFsFluid has been made already
 if(allDOFsFluid) return;
 int numele = (*(geoSource->getPackedEsetFluid())).last(); // include phantoms here

 int iele;
 int *pointers = new int[numele+1];
 pointers[0] = 0;
 maxNumDOFsFluid  = 0;
 for(iele=0; iele < numele; ++iele) {
   int numDOFs = (*(geoSource->getPackedEsetFluid()))[iele]->numDofs();
   if(numDOFs > maxNumDOFsFluid) maxNumDOFsFluid = numDOFs;
   pointers[iele+1] = pointers[iele] + numDOFs;
 }
 int *targets = new int[ pointers[numele] ];
 for(iele=0; iele < numele; ++iele)
   (*(geoSource->getPackedEsetFluid()))[iele]->dofs(*dsaFluid, targets + pointers[iele]);
 allDOFsFluid = new Connectivity(numele, pointers, targets);

 maxNumNodesFluid = 0;
 // compute maximum number of nodes per element
 for(iele=0; iele < numele; ++iele) {
   int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
   maxNumNodesFluid = std::max(maxNumNodesFluid, numNodesPerElement);
 }
}

void
Domain::makeThicknessGroupElementFlag()
{
  std::map<int, Group> &group = geoSource->group;
  std::map<int, AttributeToElement> &atoe = geoSource->atoe;
  thgreleFlag = new bool[numele];
  thpaIndex = new int[numele];

  int iele;
  for(iele=0; iele< numele; iele++) { thgreleFlag[iele] = false; thpaIndex[iele] = -1; }
  for(int iparam = 0; iparam < thicknessGroups.size(); ++iparam) {
    int groupIndex = thicknessGroups[iparam];
    for(int aindex = 0; aindex < group[groupIndex].attributes.size(); ++aindex) {
      for(int eindex = 0; eindex < atoe[group[groupIndex].attributes[aindex]].elems.size(); ++eindex) {
        iele = atoe[group[groupIndex].attributes[aindex]].elems[eindex];
        thgreleFlag[iele] = true;
        thpaIndex[iele] = iparam;
      }
    }
  }
}

// This routine creates the array of boundary condition
// in the global DOF vector
void
Domain::make_bc(int *bc, double *bcx)
{
 // Initialize all boundary conditions to free
 // and all boundary condition values to zero
 int numdof = dsa->size();
 int i;
 for(i=0; i<numdof; ++i) {
   bc[i]  = BCFREE;
   bcx[i] = 0.0;
 }

 // Set the Neuman boundary conditions
 for(i=0; i<numNeuman; ++i) {
   int dof  = dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found FORCE on non-existant dof: node %d dof %d\n",
     //          nbc[i].nnum+1,nbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCLOAD) {
     //filePrint(stderr," *** WARNING: check input, multiple FORCEs defined at node %d"
     //          ", dof %d\n",nbc[i].nnum+1,nbc[i].dofnum+1);
     bcx[dof] += nbc[i].val;
   }
   else {
     bc[dof] = BCLOAD;
     bcx[dof] = nbc[i].val;
   }
 }

 // Set the real part of the Complex Neumann boundary conditions
 for(i=0; i<numComplexNeuman; ++i) {
   int dof  = dsa->locate(cnbc[i].nnum, 1 << cnbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found FORCE on non-existant dof: node %d dof %d\n",
     //          cnbc[i].nnum+1,cnbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCLOAD) {
     //filePrint(stderr," *** WARNING: check input, multiple FORCEs defined at node %d"
     //          ", dof %d\n",cnbc[i].nnum+1,cnbc[i].dofnum+1);
     bcx[dof] += cnbc[i].reval;
   }
   else {
     bc[dof] = BCLOAD;
     bcx[dof] = cnbc[i].reval;
   }
 }

 // Set the dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
   int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found DISP on non-existant dof: node %d dof %d\n",
     //                dbc[i].nnum+1,dbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCFIXED) {
     //filePrint(stderr," *** WARNING: check input, found repeated DISP"
     //               " (node %d, dof %d)\n",dbc[i].nnum+1,dbc[i].dofnum+1);
   }

   bc[dof] = BCFIXED;
   bcx[dof] = dbc[i].val;
 }

 // Set the real part of the Complex Dirichlet boundary conditions
 for(i=0; i<numComplexDirichlet; ++i) {
   int dof  = dsa->locate(cdbc[i].nnum, 1 << cdbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found DISP on non-existant dof: node %d dof %d\n",
     //                cdbc[i].nnum+1,cdbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCFIXED) {
     //filePrint(stderr," *** WARNING: check input, found repeated Complex DISP"
     //               " (node %d, dof %d)\n",cdbc[i].nnum+1,cdbc[i].dofnum+1);
   }
   bc[dof] = BCFIXED;
   bcx[dof] = cdbc[i].reval;
 }
}

void
Domain::make_bc(int *bc, DComplex *bcx)
{
 int numdof = dsa->size();
 int i;
 for(i=0; i<numdof; ++i) {
   bc[i]  = BCFREE;
   bcx[i] = 0.0;
 }

// Set the real Neumann boundary conditions
 for(i=0; i<numNeuman; ++i) {
   int dof  = dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCLOAD) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " FORCE (node %d, dof %d)\n",nbc[i].nnum,nbc[i].dofnum);
   }
   bc[dof] = BCLOAD;
   bcx[dof] = DComplex(nbc[i].val,0.0);
 }

// Set the Complex Neuman boundary conditions
 for(i=0; i<numComplexNeuman; ++i) {
   int dof  = dsa->locate(cnbc[i].nnum, 1 << cnbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCLOAD) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " HFORCE (node %d, dof %d)\n",
        //                 cnbc[i].nnum,cnbc[i].dofnum);
   }
   bc[dof] = BCLOAD;
   bcx[dof] = DComplex(cnbc[i].reval,cnbc[i].imval);
 }

// Set the real Dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
   int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCFIXED) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " DISP (node %d, dof %d)\n",dbc[i].nnum,dbc[i].dofnum);
   }
   bc[dof] = BCFIXED;
   bcx[dof] = DComplex(dbc[i].val,0.0);
 }

// Set the Complex Dirichlet boundary condtions
 for(i=0; i<numComplexDirichlet; ++i) {
   int dof  = dsa->locate(cdbc[i].nnum, 1 << cdbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCFIXED) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " HDISP (node %d, dof %d)\n",
        //               cdbc[i].nnum,cdbc[i].dofnum);
   }
   bc[dof] = BCFIXED;
   bcx[dof] = DComplex(cdbc[i].reval, cdbc[i].imval);
 }
}


int
Domain::addDMass(int n, int d, double mass, int jdof)
{
 DMassData *dmass = new DMassData;
 dmass->next = firstDiMass;
 firstDiMass = dmass;

 firstDiMass->node   = n;
 firstDiMass->dof    = d;
 firstDiMass->diMass = mass;
 firstDiMass->jdof   = jdof;

 if(jdof >= 0 && d != jdof) sinfo.inertiaLumping = 2;

 nDimass++;

 return 0;
}

// set Linear Multipoint Constrains
int
Domain::addLMPC(LMPCons *_MPC, bool checkflag)
{
 // as currently implemented all LMPCs (both real & complex) must use unique id numbers
 // eg if LMPC (real) id numbers 1 to 10, CLMPC (complex) should be 11 to 20
 // however you can use CLMPC to add a complex term to a previously real LMPC with the same id number
 // also inequality constraint id numbers can be the same as equality constraint id numbers
 // NOTE: if checkflag is set to false, then always create new LMPC, used for various internal calls but not Parser
 if(!checkflag) { lmpc[numLMPC++] = _MPC; if(_MPC->type == 1) numCTC++; return numLMPC-1; }

 // Verify if lmpc was already defined
 int i=0;
 while((i < numLMPC) && ((lmpc[i]->lmpcnum != _MPC->lmpcnum) || lmpc[i]->type != _MPC->type)) i++;

 // if LMPC not previously defined create new
 if(i==numLMPC) {
   lmpc[numLMPC++] = _MPC;
   if(_MPC->isComplex) numComplexLMPC++;
   if(_MPC->type == 1) numCTC++;
 }
 // if LMPC already defined overwrite rhs and add terms
 else {
  //filePrint(stderr," *** WARNING: rhs of LMPC number %d overwritten\n", lmpc[i]->lmpcnum);
  //filePrint(stderr," ***          previous terms are added\n");
  int j;
  for(j=0;j<lmpc[i]->nterms;j++)
    _MPC->addterm(&lmpc[i]->terms[j]);
  delete lmpc[i];
  lmpc[i] = _MPC;
 }
 return i; // PJSA: return global MPC id
}

void Domain::printLMPC()
{
 filePrint(stderr," ... LMPC list :                    ...\n");
 int i;
 for(i=0; i<numLMPC; i++) {
   filePrint(stderr," ...    lmpc number %d (i = %d) : ", lmpc[i]->lmpcnum, i);
   if(lmpc[i]->isComplex) {
     filePrint(stderr,"rhs = (%f,%f)\n", lmpc[i]->rhs.c_value.real(), lmpc[i]->rhs.c_value.imag());
     filePrint(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       filePrint(stderr,"        %d        %d        (%f,%f)\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.c_value.real(), lmpc[i]->terms[j].coef.c_value.imag());
   }
   else {
     filePrint(stderr,"rhs = %f\n", lmpc[i]->rhs.r_value);
     filePrint(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       filePrint(stderr,"        %d        %d        %20.16f\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.r_value);
   }
 }
}

void Domain::printLMPC2()
{
 // use fprintf instead of filePrint
 fprintf(stderr," ... LMPC list :                    ...\n");
 int i;
 for(i=0; i<numLMPC; i++) {
   fprintf(stderr," ...    lmpc number %d : ", lmpc[i]->lmpcnum);
   if(lmpc[i]->isComplex) {
     fprintf(stderr,"rhs = (%f,%f)\n", lmpc[i]->rhs.c_value.real(), lmpc[i]->rhs.c_value.imag());
     fprintf(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       fprintf(stderr,"        %d        %d        (%f,%f)\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.c_value.real(), lmpc[i]->terms[j].coef.c_value.imag());
   }
   else {
     fprintf(stderr,"rhs = %f\n", lmpc[i]->rhs.r_value);
     fprintf(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       fprintf(stderr,"        %d        %d        %20.16f\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.r_value);
   }
 }
}


void Domain::normalizeLMPC()
{
 // PJSA 5-24-06
 for(int i=0; i<numLMPC; i++) {
   if(lmpc[i]->isComplex) {
     std::cerr << " *** WARNING: normalizeLMPC not implemented for complex coefficients \n";
   }
   else {
     double cnorm = 0.0;
     for(int j=0; j<lmpc[i]->nterms; j++)
       cnorm += lmpc[i]->terms[j].coef.r_value * lmpc[i]->terms[j].coef.r_value;
     cnorm = sqrt(cnorm);
     for(int j=0; j<lmpc[i]->nterms; j++)
       lmpc[i]->terms[j].coef.r_value /= cnorm;
     lmpc[i]->rhs.r_value /= cnorm;
   }
 }
}

void Domain::setPrimalLMPCs(int& numDual, int &numPrimal)
{
 numDual = numLMPC; numPrimal = 0;
 if(solInfo().solvercntl->fetiInfo.mpcflag == 2) { // convert all dual mpcs (type 0) to primal
   for(int i=0; i<numLMPC; i++)
     if(lmpc[i]->type == 0) {
       lmpc[i]->type = 3;
       numPrimal++;
       numDual--;
     }
 }
/*
 if(solInfo().solvercntl->fetiInfo.cmpc) { // convert corner bmpcs (type 2) to primal
   if(cornerWeight) {
     for(int i=0; i<numLMPC; i++)
       if(lmpc[i]->type == 2) {
         if(cornerWeight[9*lmpc[i]->terms[0].nnum+lmpc[i]->terms[0].dofnum] > 0) {
           lmpc[i]->type = 5;
           numPrimal++;
           numDual--;
         }
       }
   }
 }
*/
 //std::cerr << "numDual = " << numDual << ", numPrimal = " << numPrimal << std::endl;
 // note other strategies may be implemented here
}

Connectivity *
Domain::makeLmpcToNode()
{
 int size = 0;
 // find size of target: total number of coefficients
 // involving a different node
 // Note: done by double loop because assume number of terms is small
 int numtarget = 0;
 int i,j,jj;
 for(i=0; i<numLMPC; i++){
   if(lmpc[i]->isPrimalMPC()) continue;  // skip primal mpcs
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j) numtarget++;
   }
   size++;
 }
 // fill target with coefficient nodes
 int *pointer = new int[size+1];
 int *target  = new int[numtarget];
 size = 0; numtarget = 0;
 for(i=0; i<numLMPC; i++){
   if(lmpc[i]->isPrimalMPC()) continue;  // skip primal mpcs
   pointer[size]=numtarget;
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j){ target[numtarget]=lmpc[i]->terms[j].nnum; numtarget++;}
   }
   size++;
 }
 pointer[size]=numtarget;
 return new Connectivity(size,pointer,target);
}

Connectivity *
Domain::makeLmpcToNode_primal()
{
 int size = 0;
 // find size of target: total number of coefficients
 // involving a different node
 // Note: done by double loop because assume number of terms is small
 int numtarget = 0;
 int i,j,jj;
 for(i=0; i<numLMPC; i++){
   if(!lmpc[i]->isPrimalMPC()) continue;  // skip dual mpcs
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j) numtarget++;
   }
   size++;
 }
 // fill target with coefficient nodes
 int *pointer = new int[size+1];
 int *target  = new int[numtarget];
 size = 0; numtarget = 0;
 for(i=0; i<numLMPC; i++){
   if(!lmpc[i]->isPrimalMPC()) continue;  // skip dual mpcs
   pointer[size]=numtarget;
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j){ target[numtarget]=lmpc[i]->terms[j].nnum; numtarget++;}
   }
   size++;
 }
 pointer[size]=numtarget;
 return new Connectivity(size,pointer,target);
}

void
Domain::makeFsiToNode()
{
 int size = numFSI;
 // find size of target: total number of coefficients
 // involving a different node
 // Note: done by double loop because assume number of terms is small
 int numtarget = 0;
 int i,j,jj;
 for(i=0; i<size; i++){
   numtarget++;  // for fluid node
   for(j=0; j<fsi[i]->nterms; j++){
     if(fsi[i]->terms[j].nnum==fsi[i]->lmpcnum) continue; //HB: case of same fluid & structure node
     for(jj=0; jj<j; jj++) // look if not already found among structure nodes
       if(fsi[i]->terms[j].nnum==fsi[i]->terms[jj].nnum) break;
     if(jj==j) numtarget++;
   }
 }
 // fill target with coefficient nodes
 int *pointer = new int[size+1];
 int *target  = new int[numtarget];
 int count = 0;
 for(i=0; i<numFSI; i++){
   pointer[i]=count;
   target[count++] = fsi[i]->lmpcnum;  // fluid node
   for(j=0; j<fsi[i]->nterms; j++){
     if(fsi[i]->terms[j].nnum==fsi[i]->lmpcnum) continue; //HB: case of same fluid & structure node
     for(jj=0; jj<j; jj++) // look if not already found among structure nodes
       if(fsi[i]->terms[j].nnum==fsi[i]->terms[jj].nnum) break;
     if(jj==j) { target[count]=fsi[i]->terms[j].nnum; count++; }
   }
 }
 pointer[i]=numtarget;
 if (count!=pointer[i]) filePrint(stderr,"*** ERROR in Domain::makeFsiToNode");  // Check
 fsiToNode = new Connectivity(size,pointer,target);
 nodeToFsi = fsiToNode->alloc_reverse();
}

void Domain::getInterestingDofs(DofSet &ret, int glNode)
{
  if(glNode > nodeToFsi->csize()) return;
  for(int i=0; i<nodeToFsi->num(glNode); ++i) {
    int iFSI = (*nodeToFsi)[glNode][i];
    if(fsi[iFSI]->lmpcnum == glNode) ret.mark(DofSet::Helm) ;// check fluid dof
    for(int j=0; j<fsi[iFSI]->nterms; j++) // check structure dofs
      if(fsi[iFSI]->terms[j].nnum == glNode) ret.mark(1 << fsi[iFSI]->terms[j].dofnum);
  }
}

double ** Domain::getCMatrix()
//returns the coupling matrix for hydroelastic vibration problems
{
  if(C_condensed) return C_condensed;
  int np = c_dsaFluid->size();
  int nu = c_dsa->size();
  nuNonZero = 0;
  double ** C = new double * [nu];
  int* umap_inv = new int[nu]; // map from u wet interface numbering to c_dsa numbering
  int* umap_add_temp = new int[nu];
  int** pmap_inv = new int* [nu];
  int* npNonZero_full = new int[nu];

  for(int i=0; i<nu; ++i)   {
    C[i] = 0;
    umap_inv[i] = -1;
    umap_add_temp[i] = -1;
    pmap_inv[i] = new int[np];
    npNonZero_full[i] = 0;
    for (int j=0; j<np; ++j)  {
      pmap_inv[i][j] = -1;
    }
  }

  for(int iFSI=0; iFSI<numFSI; ++iFSI) {

    int pindex = c_dsaFluid->locate(fsi[iFSI]->fluid_node, DofSet::Potential);

    if (pindex < 0) {
      continue;
    }

    for(int j=0; j<fsi[iFSI]->nterms; j++) {
      int uindex = c_dsa->locate(fsi[iFSI]->terms[j].nnum, 1 << fsi[iFSI]->terms[j].dofnum);

      if (uindex < 0)  {
        continue;
      }

      double C_up = fsi[iFSI]->terms[j].coef.r_value;
      if (C_up != 0.) {
        if(C[uindex] == 0) {
          C[uindex] = new double[np];
          for (int jp=0; jp<np; ++jp)
            C[uindex][jp]=0;
          umap_inv[uindex] = nuNonZero;
          umap_add_temp[nuNonZero++] = dsa->locate(fsi[iFSI]->terms[j].nnum, 1 << fsi[iFSI]->terms[j].dofnum);
        }
        C[uindex][pindex] = C_up;
        pmap_inv[uindex][pindex] = (npNonZero_full[uindex])++;
      }
    }
  }

  umap = new int[nuNonZero];
  umap_add = new int[nuNonZero];
  pmap = new int*[nuNonZero];
  npNonZero = new int[nuNonZero];
  C_condensed = new double * [nuNonZero];

  for(int i=0; i<nuNonZero; ++i)  {
    C_condensed[i] = 0;
    umap[i] = -1;
    npNonZero[i] = 0;
  }

  for(int i=0; i<nu; ++i)  {
    int iu = umap_inv[i];
    if (iu >= 0 )  {
      umap[iu] = i;
      pmap[iu] = new int[npNonZero_full[i]];
      for (int j=0; j < npNonZero_full[i]; ++j)  {
        pmap[iu][j] = -1;
      }
      for (int jp=0; jp < np; ++jp)  {
        int K = pmap_inv[i][jp];
        if (K>=0)  {
          pmap[iu][K] = jp;
        }
      }
    }
  }

  for(int i=0; i<nuNonZero; ++i)  {
    umap_add[i] = umap_add_temp[i];
    C_condensed[i] = C[umap[i]];
    npNonZero[i] = npNonZero_full[umap[i]];
  }
  delete [] C;
  return C_condensed;
}

void
Domain::trMultCV(const Vector& x, Vector& y)
{
/*
  // computes y = C^T*x
  double** C_NZrows = getCMatrix();

  int nnp = domain->numUnconFluid();
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  int *unconstrNum = c_dsa->getUnconstrNum();

  y.zero();
  for(int i=0; i < domain->nuNonZero; ++i) {
    for(int jp=0; jp < nnp; ++jp) {
      y[jp] += rhoFluid*C_NZrows[i][jp]*x[unconstrNum[umap_add[i]]];
    }
  }
*/
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  y.zero();
  for(int i=0; i<numFSI; ++i) {
    int pindex = c_dsaFluid->locate(fsi[i]->fluid_node, DofSet::Potential);
    if (pindex < 0) continue;
    for(int j=0; j<fsi[i]->nterms; j++) {
      int uindex = c_dsa->locate(fsi[i]->terms[j].nnum, 1 << fsi[i]->terms[j].dofnum);
      if (uindex < 0) continue;
      double coef = fsi[i]->terms[j].coef.r_value;
      y[pindex] += rhoFluid*coef*x[uindex];
    }
  }
}

void
Domain::multCV(const Vector& x, Vector& y)
{
/*
  // computes y = C*x
  double** C_NZrows = getCMatrix();

  int nnp = domain->numUnconFluid();
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  int *unconstrNum = c_dsa->getUnconstrNum();

  y.zero();
  for(int i=0; i<domain->nuNonZero; ++i) {
    for(int k=0; k<domain->npNonZero[i]; ++k) {
      int K = domain->pmap[i][k];
      y[unconstrNum[umap_add[i]]] += rhoFluid*C_NZrows[i][K]*x[K];
    }
  }
*/
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  y.zero();
  for(int i=0; i<numFSI; ++i) {
    int pindex = c_dsaFluid->locate(fsi[i]->fluid_node, DofSet::Potential);
    if (pindex < 0) continue;
    for(int j=0; j<fsi[i]->nterms; j++) {
      int uindex = c_dsa->locate(fsi[i]->terms[j].nnum, 1 << fsi[i]->terms[j].dofnum);
      if (uindex < 0) continue;
      double coef = fsi[i]->terms[j].coef.r_value;
      y[uindex] += rhoFluid*coef*x[pindex];
    }
  }
}

// set dirichlet boundary conditions
int Domain::setDirichlet(int _numDirichlet, BCond *_dbc)
{
  // dbc = dirichlet boundary conditions
  if (dbc) {
    // Allocate memory for correct number of dbc
    BCond *nd = new BCond[numDirichlet+_numDirichlet];

    // copy old dbc
    int i;
    for (i = 0; i < numDirichlet; ++i)
      nd[i] = dbc[i];

    // copy new dbc
    for (i = 0; i<_numDirichlet; ++i)
      nd[i+numDirichlet] = _dbc[i];

    // set correct number of dbc
    numDirichlet += _numDirichlet;

    // delete old array of dbc
    delete [] dbc;

    // set new pointer to correct number of dbc
    dbc = nd;
  }
  else {
    numDirichlet = _numDirichlet;
    dbc          = _dbc;
  }

  return 0;
}

// set dirichlet boundary conditions in Fluid
int Domain::setDirichletFluid(int _numDirichletFluid, BCond *_dbcFluid)
{
  if (dbcFluid) {
    // Allocate memory for correct number of dbcFluid
    BCond *ndFluid = new BCond[numDirichletFluid+_numDirichletFluid];

    // copy old dbcFluid
    int i;
    for (i = 0; i < numDirichletFluid; ++i)
      ndFluid[i] = dbcFluid[i];

    // copy new dbcFluid
    for (i = 0; i<_numDirichletFluid; ++i)
      ndFluid[i+numDirichletFluid] = _dbcFluid[i];

    // set correct number of dbcFluid
    numDirichletFluid += _numDirichletFluid;

    // delete old array of dbcFluid
    delete [] dbcFluid;

    // set new pointer to correct number of dbcFluid
    dbcFluid = ndFluid;
  }
  else {
    numDirichletFluid = _numDirichletFluid;
    dbcFluid          = _dbcFluid;
  }

  return 0;
}

int Domain::setNeuman(int _numNeuman, BCond *_nbc)
{
 if(nbc) {
   // Allocate memory for correct number of dbc
   BCond *nd = new BCond[numNeuman+_numNeuman];

   // copy old nbc
   int i;
   for(i=0; i < numNeuman; ++i)
      nd[i] = nbc[i];

   // copy new nbc
   for(i=0; i<_numNeuman; ++i)
     nd[i+numNeuman] = _nbc[i];

   // set correct number of dbc
   numNeuman += _numNeuman;

   // delete old array of dbc
   delete [] nbc;

   // set new pointer to correct number of dbc
   nbc = nd;
 }
 else {
   numNeuman = _numNeuman;
   nbc          = _nbc;
 }
 return 0;
}

int Domain::setNeumanModal(int _numNeumanModal, BCond *_nbcModal)
{
 if(nbcModal) {
   // Allocate memory for correct number of dbc
   BCond *nd = new BCond[numNeumanModal+_numNeumanModal];

   // copy old nbcModal
   int i;
   for(i=0; i < numNeumanModal; ++i)
      nd[i] = nbcModal[i];

   // copy new nbcModal
   for(i=0; i<_numNeumanModal; ++i)
     nd[i+numNeumanModal] = _nbcModal[i];

   // set correct number of dbc
   numNeumanModal += _numNeumanModal;

   // delete old array of dbc
   delete [] nbcModal;

   // set new pointer to correct number of dbc
   nbcModal = nd;
 }
 else {
   numNeumanModal = _numNeumanModal;
   nbcModal          = _nbcModal;
 }
 return 0;
}

void
Domain::setLoadFactor(int loadcase_id, int loadset_id, double load_factor)
{
 loadfactor[std::pair<int,int>(loadcase_id,loadset_id)] = load_factor;
}

void
Domain::setLoadFactorMFTT(int loadcase_id, int loadset_id, int table_id)
{
 loadfactor_mftt[std::pair<int,int>(loadcase_id,loadset_id)] = table_id;
}

void
Domain::setLoadFactorHFTT(int loadcase_id, int loadset_id, int table_id)
{
 loadfactor_hftt[std::pair<int,int>(loadcase_id,loadset_id)] = table_id;
}

void
Domain::setLoadFactorTemp(int loadcase_id, bool load_factor)
{
 loadfactor_temp[loadcase_id] = load_factor;
}

void
Domain::setLoadFactorGrav(int loadcase_id, bool load_factor)
{
 loadfactor_grav[loadcase_id] = load_factor;
}

void
Domain::checkCases()
{
  for(std::list<int>::const_iterator it1 = sinfo.loadcases.begin(); it1 != sinfo.loadcases.end(); ++it1) {
    bool found = false;
    for(std::map<std::pair<int,int>,double>::const_iterator it2 = loadfactor.begin(); it2 != loadfactor.end(); ++it2) {
      if(it2->first.first == *it1) {
        found = true;
        continue;
      }
    }
    for(std::map<std::pair<int,int>,int>::const_iterator it2 = loadfactor_mftt.begin(); it2 != loadfactor_mftt.end(); ++it2) {
      if(it2->first.first == *it1) {
        found = true;
        continue;
      }
    }
    for(std::map<std::pair<int,int>,int>::const_iterator it2 = loadfactor_hftt.begin(); it2 != loadfactor_hftt.end(); ++it2) {
      if(it2->first.first == *it1) {
        found = true;
        continue;
      }
    }
    if(!found && *it1 != 0) {
      filePrint(stderr, " *** ERROR: Selected load case %d is not defined\n", *it1);
      exit(-1);
    }
  }
}

double
Domain::getLoadFactor(int loadset_id) const
{
  int loadcase_id = (domain->solInfo().loadcases.size() > 0) ? domain->solInfo().loadcases.front() : 0;
  std::map<std::pair<int,int>,double>::const_iterator it1 = loadfactor.find(std::pair<int,int>(loadcase_id,loadset_id));
  if(it1 != loadfactor.end()) {
    // 1. In this case, a load factor has been defined for the (loadcase_id,loadset_id) pair under LOADCASE.
    return it1->second;
  }
  else if(loadset_id == 0 && loadcase_id == 0) {
    // 2. In this case, no load factor has been defined for the (loadcase_id,loadset_id) pair under LOADCASE.
    //    The policy is to use the default load factor (1.0) for the (0,0) (loadcase_id,loadset_id) pair.
    return 1.0;
  }
  else {
    // 3. In this case, loadset_id is not included in the loadcase.
    return 0.0;
  }
}

int
Domain::gravityFlag()
{
  int loadcase_id = (domain->solInfo().loadcases.size() > 0) ? domain->solInfo().loadcases.front() : 0;
  std::map<int,bool>::const_iterator it = loadfactor_grav.find(loadcase_id);
  if(it != loadfactor_grav.end() && it->second == false) {
    // in this case, gravity has been explicitly deactivated for the current loadcase
    return 0;
  }
  else {
//    return (gravityAcceleration ? 1: 0) || (domain->solInfo().soltyp == 2);
    return (gravityAcceleration ? 1: 0); // 101416 JAT
  }
}

int
Domain::thermalFlag()
{
  int loadcase_id = (domain->solInfo().loadcases.size() > 0) ? domain->solInfo().loadcases.front() : 0;
  std::map<int,bool>::const_iterator it = loadfactor_temp.find(loadcase_id);
  if(it != loadfactor_temp.end() && it->second == false) {
    // in this case, thermal loads have been explicitly deactivated for the current loadcase
    return 0;
  }
  else {
    return sinfo.thermalLoadFlag || sinfo.thermoeFlag >= 0;
  }
}

int
Domain::setMFTT(MFTTData *_mftval, int key)
{
 mftval[key] = _mftval;
 return 0;
}

MFTTData *
Domain::getDefaultMFTT() const
{
  // the default table should be used (when defined) for any load types
  // which do not support the LOADSET_ID construct
  std::map<int,MFTTData*>::const_iterator it = domain->mftval.find(0);
  return (it != domain->mftval.end()) ? it->second : NULL;
}

MFTTData *
Domain::getMFTT(int loadset_id) const 
{
  int loadcase_id = (domain->solInfo().loadcases.size() > 0) ? domain->solInfo().loadcases.front() : 0;
  std::map<std::pair<int,int>,int>::const_iterator it1 = loadfactor_mftt.find(std::pair<int,int>(loadcase_id,loadset_id));
  if(it1 != loadfactor_mftt.end()) {
    // 1. In this case, a table has been defined for the (loadcase_id,loadset_id) pair under LOADCASE.
    int table_id = it1->second;
    std::map<int,MFTTData*>::const_iterator it2 = mftval.find(table_id);
    if(it2 != mftval.end()) return it2->second;
    else {
      filePrint(stderr, " *** WARNING: Table %d assigned to loadset %d for loadcase #%d is not defined\n",
                table_id, loadset_id, loadcase_id); 
      return NULL;
    }
  }
  else if(loadset_id == 0 && loadcase_id == 0) {
    // 2. In this case, no table has been defined for the (loadcase_id,loadset_id) pair under LOADCASE.
    //    The policy is to use the default table for the (0,0) (loadcase_id,loadset_id) pair, if it exists.
    std::map<int,MFTTData*>::const_iterator it2 = mftval.find(0);
    return (it2 != mftval.end()) ? it2->second : NULL;
  }
  else {
    // 3. In this case, loadset_id is not included in the loadcase.
    return NULL;
  }
}

int
Domain::getNumMFTT() const
{
  return mftval.size();
}

int
Domain::setHFTT(MFTTData *_hftval, int key)
{
 hftval[key] = _hftval;
 return 0;
}

MFTTData *
Domain::getDefaultHFTT() const
{
  // the default table should be used (when defined) for any load types
  // which do not support the LOADSET_ID construct
  std::map<int,MFTTData*>::const_iterator it = domain->hftval.find(0);
  return (it != domain->hftval.end()) ? it->second : NULL;
}

MFTTData *
Domain::getHFTT(int loadset_id) const
{
  int loadcase_id = (domain->solInfo().loadcases.size() > 0) ? domain->solInfo().loadcases.front() : 0;
  std::map<std::pair<int,int>,int>::const_iterator it1 = loadfactor_hftt.find(std::pair<int,int>(loadcase_id,loadset_id));
  if(it1 != loadfactor_hftt.end()) {
    // 1. In this case, a table has been defined for the (loadcase_id,loadset_id) pair under LOADCASE.
    int table_id = it1->second;
    std::map<int,MFTTData*>::const_iterator it2 = hftval.find(table_id);
    if(it2 != hftval.end()) return it2->second;
    else {
      filePrint(stderr, " *** WARNING: Table %d assigned to loadset %d for loadcase #%d is not defined\n",
                table_id, loadset_id, loadcase_id);
      return NULL;
    }
  }
  else if(loadset_id == 0 && loadcase_id == 0) {
    // 2. In this case, no table has been defined for the (loadcase_id,loadset_id) pair under LOADCASE.
    //    The policy is to use the default table for the (0,0) (loadcase_id,loadset_id) pair, if it exists.
    std::map<int,MFTTData*>::const_iterator it2 = hftval.find(0);
    return (it2 != hftval.end()) ? it2->second : NULL;
  }
  else {
    // 3. In this case, loadset_id is not included in the loadcase.
    return NULL;
  } 
}

int
Domain::getNumHFTT() const
{
  return hftval.size();
}

int
Domain::addYMTT(MFTTData *_ymtt)
{
 //--- Verify if ymtt was already defined
 int i = 0;
 while(i < numYMTT && ymtt[i]->getID() != _ymtt->getID()) i++;

 // if YMTT not previously defined create new
 if(i == numYMTT) ymtt[numYMTT++] = _ymtt;

 // if YMTT already defined print warning message
 else
   filePrint(stderr," *** WARNING: YMTT %d has already been defined \n", _ymtt->getID());

 return 0;
}

int
Domain::addSDETAFT(MFTTData *_sdetaft)
{
 //--- Verify if sdetat was already defined
 int i = 0;
 while(i < numSDETAFT && sdetaft[i]->getID() != _sdetaft->getID()) i++;

 // if SDETAFT not previously defined create new
 if(i == numSDETAFT) sdetaft[numSDETAFT++] = _sdetaft;

 // if SDETAFT already defined print warning message
 else
   filePrint(stderr," *** WARNING: SDETAFT %d has already been defined \n", _sdetaft->getID());

 return 0;
}

#ifdef USE_EIGEN3
int
Domain::addRUBDAFT(GenMFTTData<Eigen::Vector4d> *_rubdaft)
{
 //--- Verify if sdetat was already defined
 int i = 0;
 while(i < numRUBDAFT && rubdaft[i]->getID() != _rubdaft->getID()) i++;

 // if RUBDAFT not previously defined create new
 if(i == numRUBDAFT) rubdaft[numRUBDAFT++] = _rubdaft;

 // if RUBDAFT already defined print warning message
 else
   filePrint(stderr," *** WARNING: RUBDAFT %d has already been defined \n", _rubdaft->getID());

 return 0;
}
#endif

void Domain::printYMTT()
{
  if(numYMTT > 0) filePrint(stderr," ... YMTT list :                    ...\n");
  int i;
  for(i = 0; i < numYMTT; i++) {
    filePrint(stderr," ...    ymtt %d : \n", ymtt[i]->getID());
    filePrint(stderr," ...    temperature   Young's Modulus \n");
    int j;
    for(j = 0; j < ymtt[i]->getNumPoints(); j++)
      filePrint(stderr,"         %f     %f\n",
              ymtt[i]->getT(j), ymtt[i]->getV(j));
  }
}

int
Domain::addCTETT(MFTTData *_ctett)
{
 //--- Verify if ctett was already defined
 int i = 0;
 while(i < numCTETT && ctett[i]->getID() != _ctett->getID()) i++;

 // if CTETT not previously defined create new
 if(i == numCTETT) ctett[numCTETT++] = _ctett;

 // if CTETT already defined print warning message
 else
   filePrint(stderr," *** WARNING: CTETT %d has already been defined \n", _ctett->getID());

 return 0;
}

void Domain::printCTETT()
{
  if(numCTETT > 0) filePrint(stderr," ... TETT list :                    ...\n");
  int i;
  for(i = 0; i < numCTETT; i++) {
    filePrint(stderr," ...    tett %d : \n", ctett[i]->getID());
    filePrint(stderr," ...    temperature   Coeff of Thermal Expansion \n");
    int j;
    for(j = 0; j < ctett[i]->getNumPoints(); j++)
      filePrint(stderr,"         %f     %f\n",
              ctett[i]->getT(j), ctett[i]->getV(j));
  }
}

int
Domain::addSS1DT(MFTTData *_ss1dt)
{
 //--- Verify if ss1dt was already defined
 int i = 0;
 while(i < numSS1DT && ss1dt[i]->getID() != _ss1dt->getID()) i++;

 // if SS1DT not previously defined create new
 if(i == numSS1DT) ss1dt[numSS1DT++] = _ss1dt;

 // if SS1DT already defined print warning message
 else
   filePrint(stderr," *** WARNING: SS1DT %d has already been defined \n", _ss1dt->getID());

 return 0;
}

int
Domain::addSS2DT(SS2DTData *_ss2dt)
{
 //--- Verify if ss2dt was already defined
 int i = 0;
 while(i < numSS2DT && ss2dt[i]->getID() != _ss2dt->getID()) i++;

 // if SS2DT not previously defined create new
 if(i == numSS2DT) ss2dt[numSS2DT++] = _ss2dt;

 // if SS2DT already defined print warning message
 else
   filePrint(stderr," *** WARNING: SS2DT %d has already been defined \n", _ss2dt->getID());

 return 0;
}

int
Domain::addYSST(MFTTData *_ysst)
{
 //--- Verify if ysst was already defined
 int i = 0;
 while(i < numYSST && ysst[i]->getID() != _ysst->getID()) i++;

 // if YSST not previously defined create new
 if(i == numYSST) ysst[numYSST++] = _ysst;

 // if YSST already defined print warning message
 else
   filePrint(stderr," *** WARNING: YSST %d has already been defined \n", _ysst->getID());

 return 0;
}

int
Domain::addYSSRT(MFTTData *_yssrt)
{
 //--- Verify if yssrt was already defined
 int i = 0;
 while(i < numYSSRT && yssrt[i]->getID() != _yssrt->getID()) i++;

 // if YSSRT not previously defined create new
 if(i == numYSSRT) yssrt[numYSSRT++] = _yssrt;

 // if YSSRT already defined print warning message
 else
   filePrint(stderr," *** WARNING: YSSRT %d has already been defined \n", _yssrt->getID());

 return 0;
}

int
Domain::addYMST(MFTTData *_ymst)
{
 //--- Verify if ymst was already defined
 int i = 0;
 while(i < numYMST && ymst[i]->getID() != _ymst->getID()) i++;

 // if YMST not previously defined create new
 if(i == numYMST) ymst[numYMST++] = _ymst;

 // if YMST already defined print warning message
 else
   filePrint(stderr," *** WARNING: YMST %d has already been defined \n", _ymst->getID());

 return 0;
}

int
Domain::setIDis6(int _numIDis6, BCond *_iDis6)
{
 numIDis6 = _numIDis6;
 iDis6    = _iDis6;
 return 0;
}

int
Domain::setIDis(int _numIDis, BCond *_iDis)
{
 numIDis = _numIDis;
 iDis    = _iDis;
 return 0;
}

int
Domain::setIDisModal(int _numIDisModal, BCond *_iDisModal)
{
 numIDisModal = _numIDisModal;
 iDisModal    = _iDisModal;
 return 0;
}

int
Domain::setIVel(int _numIVel, BCond *_iVel)
{
 numIVel = _numIVel;
 iVel    = _iVel;
 return 0;
}

int
Domain::setIVelModal(int _numIVelModal, BCond *_iVelModal)
{
 numIVelModal = _numIVelModal;
 iVelModal    = _iVelModal;
 return 0;
}

int Domain::setIAcc(int , BCond *)
{
// INITIAL ACCELERATION IS NOT NEEDED CURRENTLY
// SO IT IS READ FROM THE INPUT FILE BUT NOT STORED
// numIAcc = _numIAcc;
// iAcc    = _iAcc;
 return 0;
}

void
Domain::setGravity(double ax, double ay, double az)
{
 gravityAcceleration = new double [3];
 gravityAcceleration[0] = ax;
 gravityAcceleration[1] = ay;
 gravityAcceleration[2] = az;
}

void
Domain::setGravitySloshing(double gg)
{
 gravitySloshing = gg;
}

void
Domain::setUpData(int topFlag)
{
  startTimerMemory(matrixTimers->setUpDataTime, matrixTimers->memorySetUp);

  geoSource->setUpData(topFlag);

  if(!haveNodes) { numnodes = geoSource->getNodes(nodes); haveNodes = true; }
  else numnodes = geoSource->totalNumNodes();
  numele = geoSource->getElems(packedEset);
  numele = packedEset.last();

  if(sinfo.solvercntl->type == SolverSelection::Direct) {
    if(numLMPC && domain->solInfo().rbmflg == 1 && domain->solInfo().grbm_use_lmpc) {
      if(elemToNode == 0) elemToNode = new Connectivity(packedEset.asSet());
      if(nodeToElem == 0) nodeToElem = elemToNode->alloc_reverse();
      Connectivity elemToNode_tmp(packedEset, nodeToElem);
      Connectivity nodeToElem_tmp = elemToNode_tmp.reverse();
      Connectivity nodeToNode_tmp = nodeToElem_tmp.transcon(elemToNode_tmp);
      renumb_nompc = nodeToNode_tmp.renumByComponent(0);
      int numnodes_tmp = nodeToNode_tmp.csize();
      int *order = new int[numnodes_tmp];
      for(int i=0; i<numnodes_tmp; ++i) order[i] = -1;
      for(int i=0; i<numnodes_tmp; ++i)
        if(renumb_nompc.renum[i] >= 0)
          order[renumb_nompc.renum[i]] = i;
      renumb_nompc.order = order;
    }
  }

  // set boundary conditions
  int numBC;
  int numTextBC;
  BCond *bc;
  BCond *textBC;

  // set dirichlet
  numBC = geoSource->getDirichletBC(bc);
  numTextBC = geoSource->getTextDirichletBC(textBC);
  if (numTextBC)
    geoSource->augmentBC(numTextBC, textBC, bc, numBC);
  if (solInfo().HEV) {
    int numBCFluid;
    BCond *bcFluid;
    numBCFluid = geoSource->getDirichletBCFluid(bcFluid);
    setDirichletFluid(numBCFluid, bcFluid);
  }
  setDirichlet(numBC, bc);

  // set neuman
  numBC = geoSource->getNeumanBC(bc);
  setNeuman(numBC, bc);
  numBC = geoSource->getNeumanBCModal(bc);
  setNeumanModal(numBC, bc);

  // set initial displacements
  numBC = geoSource->getIDis(bc);
  setIDis(numBC, bc);
  numBC = geoSource->getIDisModal(bc);
  if(numBC) {
    if(solInfo().keepModalInitialDisplacements()) { // keep the modal idisp and non-modal idisp separate
      setIDisModal(numBC, bc);
    }
    else { // convert the modal idisp into non-modal idisp
      filePrint(stderr, " ... Compute initial displacement from given modal basis ...\n");
      ModeData &modeData = (domain->solInfo().idis_modal_id == -1) ? ::modeData : modeDataIDis;
      int numIDisModal = modeData.numNodes*6;
      BCond *iDis_new = new BCond[numIDis+numIDisModal];
      for(int i=0; i<numIDis; ++i) iDis_new[i] = iDis[i]; 
      modeData.addMultY(numBC, bc, iDis_new+numIDis, 6);
      numIDis += numIDisModal;
      delete [] iDis;
      iDis = iDis_new;
    }
  }

  // set init disp6
  numBC = geoSource->getIDis6(bc);
  setIDis6(numBC, bc);

  // set initial velocities
  numBC = geoSource->getIVel(bc);
  setIVel(numBC, bc);
  numBC = geoSource->getIVelModal(bc);
  if(numBC) {
    if(solInfo().keepModalInitialVelocities()) { // keep the modal ivel and non-modal ivel separate
      setIVelModal(numBC, bc);
    }
    else { // convert the modal ivel into non-modal ivel
      filePrint(stderr, " ... Compute initial velocity from given modal basis ...\n");
      ModeData &modeData = (domain->solInfo().ivel_modal_id == -1) ? ::modeData : modeDataIVel;
      int numIVelModal = modeData.numNodes*6;
      BCond *iVel_new = new BCond[numIVel+numIVelModal];
      for(int i=0; i<numIVel; ++i) iVel_new[i] = iVel[i];
      modeData.addMultY(numBC, bc, iVel_new+numIVel, 6);
      numIVel += numIVelModal;
      delete [] iVel;
      iVel = iVel_new;
    }
  }

  // set Control Law
  claw = geoSource->getControlLaw();

  initNodalTemperatures();

  // PJSA: compute temperature dependent material properties
  computeTDProps();

/*
  -- Include here checks on LMPC --
     check for redundant LMPC by applying a stabilized
     Gram-Schmit on LMPC global constrain matrix (domain->lmpc)
        --> yields the left nullspace L in Lu=g
            => if rhs, i.e. g, not orth to left nullspace: ERROR
            => otherwise conserve only a set of non-redundant MPC
        (D.Rixen 04-28-99)
*/

  stopTimerMemory(matrixTimers->setUpDataTime, matrixTimers->memorySetUp);
}

#ifndef OUTPUTMESSAGE
#define OUTPUTMESSAGE
const char* OutputMessage[] = {
" ... Outputing Displacement         ... \n",
" ... Outputing Temperature          ... \n",
" ... Outputing StressXX             ... \n",
" ... Outputing StressYY             ... \n",
" ... Outputing StressZZ             ... \n",
" ... Outputing StressXY             ... \n",
" ... Outputing StressYZ             ... \n",
" ... Outputing StressXZ             ... \n",
" ... Outputing StrainXX             ... \n",
" ... Outputing StrainYY             ... \n",
" ... Outputing StrainZZ             ... \n",
" ... Outputing StrainYZ             ... \n",
" ... Outputing StrainXY             ... \n",
" ... Outputing StrainXZ             ... \n",
" ... Outputing HeatFlXX             ... \n",
" ... Outputing HeatFlXY             ... \n",
" ... Outputing HeatFlXZ             ... \n",
" ... Outputing GrdTempX             ... \n",
" ... Outputing GrdTempY             ... \n",
" ... Outputing GrdTempZ             ... \n",
" ... Outputing Von Mises Stress     ... \n",
" ... Outputing Principal Stress #1  ... \n",
" ... Outputing Principal Stress #2  ... \n",
" ... Outputing Principal Stress #3  ... \n",
" ... Outputing Principal Strain #1  ... \n",
" ... Outputing Principal Strain #2  ... \n",
" ... Outputing Principal Strain #3  ... \n",
" ... Outputing InXForce             ... \n",
" ... Outputing InYForce             ... \n",
" ... Outputing InZForce             ... \n",
" ... Outputing AXMoment             ... \n",
" ... Outputing AYMoment             ... \n",
" ... Outputing AZMoment             ... \n",
" ... Outputing Energies             ... \n",
" ... Outputing AeroForce            ... \n",
" ... Outputing EigenPair            ... \n",
" ... Outputing Von Mises Strain     ... \n",
" ... Outputing Pressure             ... \n",
" ... Outputing All 6 dofs           ... \n",
" ... Outputing Aero X Force         ... \n",
" ... Outputing Aero Y Force         ... \n",
" ... Outputing Aero Z Force         ... \n",
" ... Outputing Aero X Moment        ... \n",
" ... Outputing Aero Y Moment        ... \n",
" ... Outputing Aero Z Moment        ... \n",
" ... Outputing Velocity             ... \n",
" ... Outputing Acceleration         ... \n",
" ... Outputing Youngs Modulus       ... \n",
" ... Outputing Material Density     ... \n",
" ... Outputing Thickness            ... \n",
" ... Outputing Shape Attributes     ... \n",
" ... Outputing Shape Attributes     ... \n",
" ... Outputing Composite            ... \n",
" ... Outputing Dx                   ... \n",
" ... Outputing Dy                   ... \n",
" ... Outputing Dz                   ... \n",
" ... Outputing Rx                   ... \n",
" ... Outputing Ry                   ... \n",
" ... Outputing Rz                   ... \n",
" ... Outputing DispMod              ... \n",
" ... Outputing RotMod               ... \n",
" ... Outputing TotMod               ... \n",
" ... Outputing Rigid Body Modes     ... \n",
" ... Outputing Elem To Node Table   ... \n",
" ... Outputing Node To Elem Table   ... \n",
" ... Outputing Node To Node Table   ... \n",
" ... Outputing Aero Flux            ... \n",
" ... Outputing Heat Flux            ... \n",
" ... Outputing Grd Temp             ... \n",
" ... Outputing Velocity (6 dofs)    ... \n",
" ... Outputing Acceleration (6 dofs)... \n",
" ... Outputing Alphas Modes         ... \n",
" ... Outputing Error for Modes Decomp. ... \n",
" ... Outputing Sloshing Eigen Modes ... \n",
" ... Outputing Sloshing x-disp.     ... \n"
" ... Outputing Sloshing y-disp.     ... \n"
" ... Outputing Sloshing z-disp.     ... \n"
};
#endif

Connectivity
Domain::makeSommerToNode()
{
    int size = numSommer;
    // Find out the number of targets we will have
    std::vector<size_t> pointer(size+1) ;
    size_t pp = 0;
    for(int i=0; i < size; ++i) {
        pointer[i] = pp;
        pp += sommer[i] ? sommer[i]->numNodes() : 0;
    }
    pointer[size] = pp;
    //int numtarget = pp;
    // Create the target array
    std::vector<int> target( pp );
    // Fill it in
    for(int i=0; i < size; ++i) {
        if(sommer[i])
            sommer[i]->nodes(target.data()+pointer[i]);
    }

    return { size, std::move(pointer), std::move(target) };
}

Connectivity *
Domain::prepDirectMPC()
{
  // get all of the lmpcs from the packedEset
  // TODO don't extract elements which use lagrange multipliers or penalty
  for(int i=0; i<numLMPC; ++i) if(lmpc[i]) delete lmpc[i];
  numLMPC = 0;
  lmpc.deleteArray();
  int nEle = packedEset.last();
  int nRigid = 0;
  Elemset *rigidSet = new Elemset;
  for(int i = 0; i < nEle; ++i) {
    Element *ele = packedEset[i];
    if(ele != 0) {
      if(false /*dynamic_cast<RigidBeam*>(ele) || dynamic_cast<RigidThreeNodeShell*>(ele)
         || dynamic_cast<RigidFourNodeShell*>(ele) || dynamic_cast<RigidSolid6Dof*>(ele)*/) // TODO
        rigidSet->elemadd(nRigid++, ele);
      else {
        int n = ele->getNumMPCs();
        if(n > 0) {
          LMPCons **l = ele->getMPCs();
          for(int j = 0; j < n; ++j) {
            // first check to see if nodal frames are used
            bool use_nframes = false;
            if(l[j]->getSource() != mpc::Lmpc && l[j]->getSource() != mpc::NodalContact && l[j]->getSource() != mpc::TiedSurfaces) {
              for(int k=0; k<l[j]->nterms; ++k) {
                if(nodes.dofFrame(l[j]->terms[k].nnum)) { use_nframes = true; break; }
              }
            }
            if(use_nframes) {
              std::vector<LMPCTerm> terms_copy(l[j]->terms);
              l[j]->terms.clear();
              l[j]->nterms = 0;
              for(std::vector<LMPCTerm>::iterator it = terms_copy.begin(); it != terms_copy.end(); ++it) {
                if(NFrameData *cd = nodes.dofFrame(it->nnum)) {
                  double c[3] = { 0, 0, 0 };
                  c[it->dofnum%3] = it->coef.r_value;
                  cd->transformVector3(c);
                  for(int m=0; m<3; ++m) {
                    if(std::abs(c[m]) > std::numeric_limits<double>::epsilon()) {
                      LMPCTerm t(it->nnum, 3*(it->dofnum/3)+m, c[m]);
                      l[j]->addterm(&t);
                    }
                  }
                }
                else l[j]->addterm(&(*it));
              }
            }
            lmpc[numLMPC++] = l[j];
          }
          delete [] l;
        }
      }
    }
  }
  //printLMPC();
  if(nRigid > 0) {
    //std::cerr << "found " << nRigid << " rigid beam/shell/solid6 elements in the element set\n";
    std::set<int> blockedNodes;
    for(int i = 0; i < numDirichlet; ++i) blockedNodes.insert(dbc[i].nnum);
    rigidSet->collapseRigid6(blockedNodes);
    nRigid = rigidSet->last();
    StructProp p; // dummy property
    for(int i = 0; i < nRigid; ++i) {
      int n = (*rigidSet)[i]->getNumMPCs();
      if(n > 0) {
        if((*rigidSet)[i]->getProperty() == 0) { // this element was instantiated in Elemset::collapseRigid6
          (*rigidSet)[i]->buildFrame(nodes);     // need to call buildFrame and setProp to prep it.
          (*rigidSet)[i]->setProp(&p);
        }
        LMPCons **l = (*rigidSet)[i]->getMPCs();
        for(int j = 0; j < n; ++j)
          lmpc[numLMPC++] = l[j];
        delete [] l;
      }
    }
  }
  //std::cerr << "extracted " << numLMPC << " lmpcs from the element set\n";

  geoSource->makeDirectMPCs(numLMPC, lmpc);
  // MPC Connectivity treatment in direct way.
  std::multimap<int, int> mpcConnect;
  std::multimap<int, int>::iterator it;
  std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> ret;
  int ndMax = nodeToElem->csize();
  for(int i = 0; i < numLMPC; ++i) {
    for(int j = 1; j < lmpc[i]->nterms; ++j)
      if(lmpc[i]->terms[0].nnum != lmpc[i]->terms[j].nnum) {

        ret = mpcConnect.equal_range(lmpc[i]->terms[0].nnum);
        bool found = false;
        for(it = ret.first; it != ret.second; ++it)
          if((*it).second == lmpc[i]->terms[j].nnum) { found = true; break; }

        if(!found)
          mpcConnect.insert(std::pair<int,int>(lmpc[i]->terms[0].nnum, lmpc[i]->terms[j].nnum));
      }
  }

  int *ptr = new int[ndMax+1];
  for(int i = 0; i < ndMax; ++i)
    ptr[i] = 1;
  ptr[ndMax] = 0;
  for(it = mpcConnect.begin(); it != mpcConnect.end(); ++it)
    ptr[it->first]++;
  for(int i = 0; i < ndMax; ++i)
    ptr[i+1] += ptr[i];
  int *tg = new int[ptr[ndMax]];

  for(it = mpcConnect.begin(); it != mpcConnect.end(); ++it)
    tg[--ptr[it->first]] = it->second;
  for(int i = 0; i < ndMax; ++i)
    tg[--ptr[i]] = i;

  Connectivity renumToNode(ndMax, ptr, tg);
  Connectivity *elemToRenum = elemToNode->transcon(&renumToNode);
  Connectivity *renumToElem = elemToRenum->alloc_reverse();
  Connectivity *nodeToNodeDirect = renumToElem->transcon(elemToRenum);
  delete renumToElem;
  delete elemToRenum;
  delete rigidSet;
  return nodeToNodeDirect;
}

Renumber
Domain::getRenumbering()
{
 // create node to element connectivity from element to node connectivity
 if(nodeToElem) delete nodeToElem;
 nodeToElem = elemToNode->alloc_reverse();

 // create node to node connectivity
 nodeToNode = std::make_unique<Connectivity>( nodeToElem->transcon(*elemToNode) );

 if(solInfo().HEV == 1 && solInfo().addedMass == 1) {
   int numWetNodes = 0;
   int *wetIFNodes = getAllWetInterfaceNodes(numWetNodes);
   auto nodeToNodeTemp = std::make_unique<Connectivity>( nodeToNode->combineAll(numWetNodes,wetIFNodes) );
   nodeToNode = std::move(nodeToNodeTemp);
 }

 // get number of nodes (actually, this is the maximum node number when there are gaps in the node numbering)
 numnodes = nodeToNode->csize();

 if(solInfo().solvercntl->type == SolverSelection::Direct || solInfo().solvercntl->type == SolverSelection::Iterative)
     makeNodeToNode_sommer(); // single domain solvers

 // delete any previously allocated memory
 if(renumb.order) { delete [] renumb.order; renumb.order=0; }
 if(renumb.renum) { delete [] renumb.renum; renumb.renum=0; }
 if(renumb.xcomp) { delete [] renumb.xcomp; renumb.xcomp=0; }

 // renumber the nodes
 if((!domain->GetnContactSurfacePairs() || !sinfo.isNonLin() || tdenforceFlag()) && !sinfo.printMatLab) {
   renumb = nodeToNode->renumByComponent(sinfo.renum);
 }
 else {
   // for nonlinear mortar contact, the number of "components" of the nodeToNode connectivity may change due to new interactions
   // however, we don't want the node numbering to change therefore renumByComponent is not appropriate for this application
   renumb.numComp = 1;
   renumb.renum = new int[numnodes];
   for(int i=0; i<numnodes; ++i) renumb.renum[i] = i;
   renumb.xcomp = new int[2];
   renumb.xcomp[0] = 0; renumb.xcomp[1] = numnodes;
 }

 int *order = new int[numnodes];

 int i;
 for(i=0; i<numnodes; ++i)
   order[i] = -1;

 // note: order maps from new index to original index and renumb.renum maps from original index to new index
 int nodecount = 0;
 for(i=0; i<numnodes; ++i)
   if(renumb.renum[i] >= 0) {
     order[renumb.renum[i]] = i;
     nodecount++;
   }

 if(domain->solInfo().rbmflg == 1 && renumb.numComp > 1
    && sinfo.solvercntl->type == SolverSelection::Direct && sinfo.solvercntl->subtype == 1) {
   filePrint(stderr, "\x1B[31m *** WARNING: GRBM with sparse solver is not \n"
                          "     supported for multi-component models.  \x1B[0m\n");
 }

 if(renumb_nompc.order && renumb_nompc.numComp > 1
    && sinfo.solvercntl->type == SolverSelection::Direct && (sinfo.solvercntl->subtype == 0 || sinfo.solvercntl->subtype == 1)) {
   // altering the ordering for skyline or sparse with LMPCs due to requirements of GRBM
   DofSetArray *dsa = new DofSetArray(numnodes, packedEset, renumb.renum);
   ConstrainedDSA *c_dsa = new ConstrainedDSA(*dsa, numDirichlet, dbc);
   int p = nodecount-1;
   int min_defblk = 0;
   for(int n=renumb_nompc.numComp-1; n>=0; --n) {
     int count = 0;
     bool check = true;
     for(int i=renumb_nompc.xcomp[n+1]-1; i>=renumb_nompc.xcomp[n]; --i) {
       int inode = renumb_nompc.order[i];
       int q = renumb.renum[inode];
       int jnode = order[p];
       // swap positions of inode and jnode
       renumb.renum[inode] = p;
       renumb.renum[jnode] = q;
       order[p] = inode;
       order[q] = jnode;
       p--;
       count += c_dsa->weight(inode);
#ifdef USE_EIGEN3
       if(dsa->weight(inode) == 6) check = false;
       if(count > sinfo.solvercntl->sparse_defblk) {
         if(!check) break;
         else {
           // check for co-linearity XXX this could/should be done even for case with no LMPCs
           Eigen::Matrix<double,3,Eigen::Dynamic> M(3,renumb_nompc.xcomp[n+1]-i);
           for(int j=renumb_nompc.xcomp[n+1]-1; j>=i; --j) {
             int jnode = renumb_nompc.order[j];
             M.col(renumb_nompc.xcomp[n+1]-1-j) << nodes[jnode]->x, nodes[jnode]->y, nodes[jnode]->z;
           }
           int rank = M.colPivHouseholderQr().rank();
           if(rank > 1) break;
         }
       }
#else
       if(count > sinfo.solvercntl->sparse_defblk) break;
#endif
     }
     min_defblk += count;
   }
   delete c_dsa; delete dsa;
   sinfo.solvercntl->sparse_defblk = std::max(sinfo.solvercntl->sparse_defblk, min_defblk);
 }

 //if(sinfo.renum > 0 && verboseFlag && renumb.numComp > 1)
 //  filePrint(stdout," ... Number of components =%2d	...\n",renumb.numComp);

 Renumber ret;
 ret.order    = order;
 ret.renumb   = renumb.renum;
 renumb.order = order;

 return ret;
}

Renumber*
Domain::getRenumberingFluid()
{
 // create node to element connectivity from element to node connectivity
 if(nodeToElemFluid) delete nodeToElemFluid;
 nodeToElemFluid = elemToNodeFluid->alloc_reverse();

 // create node to node connectivity
 if(nodeToNodeFluid) delete nodeToNodeFluid;
 nodeToNodeFluid = nodeToElemFluid->transcon(elemToNodeFluid);

 // get number of nodes
 numnodesFluid = nodeToNodeFluid->csize();

 // renumber the nodes
#ifdef TFLOP
 alloca(0);
#endif
 renumbFluid = nodeToNodeFluid->renumByComponent(sinfo.renum);

 int *order = new int[numnodesFluid];

 int i;
 for(i=0; i<numnodesFluid; ++i)
   order[i] = -1;

 for(i=0; i<numnodesFluid; ++i)
   if(renumbFluid.renum[i] >= 0)
     order[renumbFluid.renum[i]] = i;

 Renumber* ret = new Renumber;
 ret->order    = order;
 ret->renumb   = renumbFluid.renum;
 renumbFluid.order = order;

 // Just for Gnu alloca function
#ifdef TFLOP
 alloca(0);
#endif

 return ret;

}
void
Domain::makeNodeToNode_sommer()
{
    if(solInfo().doFreqSweep && (numSommer > 0)) {
        Connectivity sommerToNode = makeSommerToNode();
        Connectivity nodeToSommer = sommerToNode.reverse();
        nodeToNode_sommer = nodeToSommer.transcon(&sommerToNode);
    }

    if( (solInfo().newmarkBeta==0.0)&&(solInfo().isAcoustic()) ) {
        Connectivity sommerToNode = makeSommerToNode();
        Connectivity nodeToSommer = sommerToNode.reverse();
        //Watch nodeToNode_sommer may not be of the right size !
        int numnodes = nodeToNode->csize();
        Connectivity temp0(numnodes);
        Connectivity temp1 = temp0.modify();
        Connectivity temp2 = nodeToSommer.transcon(sommerToNode);
        nodeToNode_sommer = temp1.transcon(&temp2);
    }
}

void
Domain::readInModes(int modal_id, ModeData &modeData)
{
 std::map<int,ModalParams>::iterator it = sinfo.readInModes.find(modal_id);
 if(it == sinfo.readInModes.end()) {
   filePrint(stderr," *** ERROR: A basis with rob_id = %d is specified, but has not been defined via READMODE.\n", modal_id);
   exit(-1);
 }
 const char* modeFileName = it->second.fileName.c_str();

 // Open file containing mode shapes and frequencies
 FILE *f;
 if((f=fopen(modeFileName,"r"))==(FILE *) 0) {
   filePrint(stderr," *** ERROR: Cannot open file %s specified in READMODE with rob_id %d.\n", modeFileName, modal_id);
   exit(-1);
 }
 fflush(f);

 // Read in number of modes and number of nodes
 char buf[80];
 char *str = fgets(buf, sizeof buf, f);
 bool b;
 if(strncmp("Vector", buf, 6) == 0) {
   modeData.numModes = sinfo.readInModes[modal_id].numvec;
   b = false;
 }
 else {
   int count = sscanf(buf, "%d", &modeData.numModes);
   b = (count != 1);
 }
 int count = fscanf(f, "%d", &modeData.numNodes);
 b = b || (count != 1);

 // If the file is not in one of the two valid ascii formats (see manual) then it is assumed to be binary
 if(b) {
   filePrint(stderr," ... Convert binary file to ASCII   ...\n");
   const std::string input_file_name(modeFileName);
   const std::string output_file_name = input_file_name + ".xpost";
   convert_rob<Rom::BasisBinaryInputFile, Rom::XPostOutputFile<6> >(input_file_name, output_file_name);
   sinfo.readInModes[modal_id].fileName = output_file_name;
   readInModes(modal_id, modeData);
   return;
 }

 filePrint(stderr," ... Read in Modes from file: %s ...\n",modeFileName);

 // Allocation of memory for frequencies and mode shapes
 modeData.frequencies  = new double[modeData.numModes];

 typedef double (*T6)[6];
 modeData.modes        = new T6[modeData.numModes];
 modeData.nodes        = new int[modeData.numNodes];

 // Read frequencies and mode shapes
 int iMode, iNode, numsubstring;
 char input[500], *substring;
 double tmpinput[7];

 for(iMode=0; iMode<modeData.numModes; ++iMode) {
   modeData.modes[iMode] = new double[modeData.numNodes][6];

   count = fscanf(f,"%lf\n",&modeData.frequencies[iMode]);

   for(iNode=0; iNode<modeData.numNodes; ++iNode) {

     char *c = fgets(input, 500, f);
     substring = strtok(input, " ");
     numsubstring = 0;
     do {
       tmpinput[numsubstring] = strtod(substring, NULL);
       if (strncmp(substring, "\n", 1) == 0)
         break;
       substring = strtok(NULL, " ");
       ++numsubstring;

     } while (substring != NULL);

     switch(numsubstring){
       case 3:
         modeData.nodes[iNode] = iNode;
         modeData.modes[iMode][iNode][0] = tmpinput[0];
         modeData.modes[iMode][iNode][1] = tmpinput[1];
         modeData.modes[iMode][iNode][2] = tmpinput[2];
         modeData.modes[iMode][iNode][3] = 0.0;
         modeData.modes[iMode][iNode][4] = 0.0;
         modeData.modes[iMode][iNode][5] = 0.0;
         break;
       case 4:
         modeData.nodes[iNode] = int(tmpinput[0]+0.5)-1;
         modeData.modes[iMode][iNode][0] = tmpinput[1];
         modeData.modes[iMode][iNode][1] = tmpinput[2];
         modeData.modes[iMode][iNode][2] = tmpinput[3];
         modeData.modes[iMode][iNode][3] = 0.0;
         modeData.modes[iMode][iNode][4] = 0.0;
         modeData.modes[iMode][iNode][5] = 0.0;
         break;
       case 6:
         modeData.nodes[iNode] = iNode;
         modeData.modes[iMode][iNode][0] = tmpinput[0];
         modeData.modes[iMode][iNode][1] = tmpinput[1];
         modeData.modes[iMode][iNode][2] = tmpinput[2];
         modeData.modes[iMode][iNode][3] = tmpinput[3];
         modeData.modes[iMode][iNode][4] = tmpinput[4];
         modeData.modes[iMode][iNode][5] = tmpinput[5];
         break;
       case 7:
         modeData.nodes[iNode] = int(tmpinput[0]+0.5)-1;
         modeData.modes[iMode][iNode][0] = tmpinput[1];
         modeData.modes[iMode][iNode][1] = tmpinput[2];
         modeData.modes[iMode][iNode][2] = tmpinput[3];
         modeData.modes[iMode][iNode][3] = tmpinput[4];
         modeData.modes[iMode][iNode][4] = tmpinput[5];
         modeData.modes[iMode][iNode][5] = tmpinput[6];
         break;
       default:
         filePrint(stderr, " *** ERROR: Check input for mode data (numstrings = %d) in file %s ... exiting\n", numsubstring, modeFileName);
         exit(0);
     }
   }
 }
}

void
Domain::readInShapeDerivatives(char* shapeDerFileName)
{
 filePrint(stderr," ... Read in Shape Derivatives      ...\n"
                  " ... file: %-24s ...\n",shapeDerFileName);

 // Open file containing mode shapes and frequencies.
 FILE *f;
 if((f=fopen(shapeDerFileName,"r"))==(FILE *) 0 ){
   filePrint(stderr," *** ERROR: Cannot open %s, exiting...",shapeDerFileName);
   exit(0);
 }
 fflush(f);

 char str[80];
 // Read in number of shape variables and number of nodes
 int count = fscanf(f,"%s%s%s%s%s%s",str,str,str,str,str,str);
 count = fscanf(f, "%d\n", &shapeSenData.numNodes);
 shapeSenData.numVars = 0;

 typedef double (*T3)[3];
 shapeSenData.index = new int[100];
 shapeSenData.sensitivities = new T3[100];
 shapeSenData.nodes         = new int[shapeSenData.numNodes*2]; //NOTE: assumes that global node number does not exceed twice of number of nodes

 // Read shape sensitivities
 int iSen(0), iNode, numsubstring;
 char input[500], *substring;
 double tmpinput[7];

 //int jj;
// for(iSen=0; iSen<shapeSenData.numVars; ++iSen) {
 while ( !feof(f) ) {
   shapeSenData.sensitivities[iSen] = new double[shapeSenData.numNodes][3];

   count = fscanf(f,"%d\n",&shapeSenData.index[iSen]);
   if(count < 0) break;
   shapeSenData.numVars++;
   
   //int nodeNum;
   for(iNode=0; iNode<shapeSenData.numNodes; ++iNode) {

     char *c = fgets(input, 500, f);
     substring = strtok(input, " ");
     numsubstring = 0;
     do {
       tmpinput[numsubstring] = strtod(substring, NULL);
       if (strncmp(substring, "\n", 1) == 0)
         break;
       substring = strtok(NULL, " ");
       ++numsubstring;

     } while (substring != NULL);

     if(numsubstring != 3) {
       filePrint(stderr, " *** ERROR: Check input for shape sensitivity data (numstrings = %d) in file %s ... must be 3 ... exiting\n", numsubstring, shapeDerFileName);
       exit(-1);
     } 
 
//     int cnode = (int) tmpinput[0] - 1;
//     if(cnode != iNode) filePrint(stderr, " *** ERROR: cnode differs from iNode, cnode = %d, iNode = %d\n", cnode, iNode);
     shapeSenData.nodes[iNode] = iNode;
     shapeSenData.sensitivities[iSen][iNode][0] = tmpinput[0];
     shapeSenData.sensitivities[iSen][iNode][1] = tmpinput[1];
     shapeSenData.sensitivities[iSen][iNode][2] = tmpinput[2];
   }
   iSen++;
 }
 setNumShapeVars(shapeSenData.numVars);
}

void
Domain::getElementAttr(int fileNumber,int iAttr,double time) {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  double eleattr;

  // ... OUTPUT precision
  int p = oinfo[fileNumber].precision;

  // ... Allocate Vectors for Nodal Attributes
  Vector nodattr(numnodes,0.0);
  Vector ndWeight(numnodes,0.0);

  int iele,k;
  for(iele=0; iele<numele; ++iele) {

     int NodesPerElement = elemToNode->num(iele);

     if (packedEset[iele]->getProperty()!=0) {
       switch (iAttr) {
       case YOUNG:
         eleattr=packedEset[iele]->getProperty()->E;
         break;
       case MDENS:
         eleattr=packedEset[iele]->getProperty()->rho;
         break;
       case THICK:
         eleattr=packedEset[iele]->getProperty()->eh;
         break;
       default:
         assert(0);
       }
     }
     else
       eleattr = 0.0;

     for(k=0; k<NodesPerElement; ++k) {
       nodattr[(*elemToNode)[iele][k]]  += eleattr;
       ndWeight[(*elemToNode)[iele][k]] += 1.0;
     }

     if(avgnum == 0 ) {
       for(k=0; k<NodesPerElement; ++k) {
        filePrint(oinfo[fileNumber].filptr,"% .*E",p,eleattr);
       }
       filePrint(oinfo[fileNumber].filptr,"\n");
       fflush(oinfo[fileNumber].filptr);
     }
  }

  if(avgnum == 1 ) {

    for (k=0; k<numnodes; ++k)  {
      if (ndWeight[k] == 0.0)
	{ nodattr[k] = 0.0; }
      else
	{ nodattr[k] /= ndWeight[k]; }
    }

    geoSource->outputNodeScalars(fileNumber, nodattr.data(), numnodes, time);
  }

}

void
Domain::getCompositeData(int iInfo,double time) {

  OutputInfo *oinfo = geoSource->getOutputInfo();

  int w   = oinfo[iInfo].width;
  int p   = oinfo[iInfo].precision;
  int lay = oinfo[iInfo].nodeNumber;

  // ... For time = 0.0: Output of pseudo-structure

  if ( time == 0.0 ) {

    int iele,icp;
    int numCPoints = numele*4;

    // ... Create Arrays

    if ( ! MidPoint ) {

      MidPoint = new double*[numele];

      int maxLayer=0;
      for(iele=0; iele<numele; ++iele) {
        MidPoint[iele] = packedEset[iele]->getMidPoint(nodes);
        maxLayer = std::max(maxLayer,packedEset[iele]->getCompositeLayer());
      }

      CPoint = new double**[maxLayer];

      int icl;
      for (icl=0;icl<maxLayer;icl++) {
        CPoint[icl] = new double*[numCPoints];
        for (icp=0;icp<numCPoints;icp++)
          CPoint[icl][icp] = new double[3];
      }
    }

    // ... Open top-file for pseudo-structure ".cross"

    char * basename = getBasename(oinfo[iInfo].filename);
    int fnamesize   = strlen(basename)+7;
    char * crossfile  = (char*) malloc(sizeof(char*)*(fnamesize));
    sprintf(crossfile,"%s%s%d",basename,".cross",lay+1);
    FILE * crossout = fopen(crossfile,"w");

    // ... Global Scaling factor for length of cross bars

    double minDist=10000000000.0;
    double dx,dy,dz,dd;
    int jele;
    for (iele=0; iele<numele; ++iele) {
      for (jele=iele+1; jele<numele; ++jele) {
         dx = MidPoint[iele][0]-MidPoint[jele][0];
         dy = MidPoint[iele][1]-MidPoint[jele][1];
         dz = MidPoint[iele][2]-MidPoint[jele][2];
	 dd = sqrt(dx*dx+dy*dy+dz*dz);
	 minDist = std::min(minDist,dd);
      }
    }

    crossScale = minDist * 1.25;

    // ... Determine end-points of cross bars

    for(iele=0; iele<numele; ++iele) {

       // ... get E1, E2, Phi of layer & cframe

       double * layData = packedEset[iele]->getCompositeData(lay);
       double * cFrame  = packedEset[iele]->getCompositeFrame();

       if ( ! layData || ! cFrame ) continue;

       double E1  = layData[0];
       double E2  = layData[1];
       double Phi = layData[8];

       double length1 = crossScale * E1 / std::max(E1,E2) / 2.0;
       double length2 = crossScale * E2 / std::max(E1,E2) / 2.0;

       // .... Adjust Phi

       if ( Phi < 0.0 || Phi > 360.0 ) {
         int irot = int(Phi / 360.0);
         Phi = Phi - double (irot) * 360.0;
         if ( Phi < 0.0 ) Phi = 360.0 - Phi;
       }
       Phi = (3.141593*Phi)/180.0;

       // .... determine cross points in cframe

       double * cP1 = new double[3]; double * cP2 = new double[3];
       double * cP3 = new double[3]; double * cP4 = new double[3];

       cP1[0] =  length1 * cos(Phi);
       cP1[1] =  length1 * sin(Phi);
       cP1[2] =  0.0;

       cP2[0] = -length1 * cos(Phi);
       cP2[1] = -length1 * sin(Phi);
       cP2[2] =  0.0;

       cP3[0] =  length2 * cos(Phi+1.5708);
       cP3[1] =  length2 * sin(Phi+1.5708);
       cP3[2] =  0.0;

       cP4[0] = -length2 * cos(Phi+1.5708);
       cP4[1] = -length2 * sin(Phi+1.5708);
       cP4[2] =  0.0;

       // ... transform into global frame

       int n1 = iele*4; int n2 = n1+1; int n3 = n1+2; int n4 = n1+3;

       CPoint[lay][n1][0] = MidPoint[iele][0] + cFrame[0]*cP1[0] +  cFrame[3]*cP1[1] + cFrame[6]*cP1[2];
       CPoint[lay][n1][1] = MidPoint[iele][1] + cFrame[1]*cP1[0] +  cFrame[4]*cP1[1] + cFrame[7]*cP1[2];
       CPoint[lay][n1][2] = MidPoint[iele][2] + cFrame[2]*cP1[0] +  cFrame[5]*cP1[1] + cFrame[8]*cP1[2];

       CPoint[lay][n2][0] = MidPoint[iele][0] + cFrame[0]*cP2[0] +  cFrame[3]*cP2[1] + cFrame[6]*cP2[2];
       CPoint[lay][n2][1] = MidPoint[iele][1] + cFrame[1]*cP2[0] +  cFrame[4]*cP2[1] + cFrame[7]*cP2[2];
       CPoint[lay][n2][2] = MidPoint[iele][2] + cFrame[2]*cP2[0] +  cFrame[5]*cP2[1] + cFrame[8]*cP2[2];

       CPoint[lay][n3][0] = MidPoint[iele][0] + cFrame[0]*cP3[0] +  cFrame[3]*cP3[1] + cFrame[6]*cP3[2];
       CPoint[lay][n3][1] = MidPoint[iele][1] + cFrame[1]*cP3[0] +  cFrame[4]*cP3[1] + cFrame[7]*cP3[2];
       CPoint[lay][n3][2] = MidPoint[iele][2] + cFrame[2]*cP3[0] +  cFrame[5]*cP3[1] + cFrame[8]*cP3[2];

       CPoint[lay][n4][0] = MidPoint[iele][0] + cFrame[0]*cP4[0] +  cFrame[3]*cP4[1] + cFrame[6]*cP4[2];
       CPoint[lay][n4][1] = MidPoint[iele][1] + cFrame[1]*cP4[0] +  cFrame[4]*cP4[1] + cFrame[7]*cP4[2];
       CPoint[lay][n4][2] = MidPoint[iele][2] + cFrame[2]*cP4[0] +  cFrame[5]*cP4[1] + cFrame[8]*cP4[2];

       delete [] cP1;
       delete [] cP2;
       delete [] cP3;
       delete [] cP4;
    }

    // ... print pseudo structure

    filePrint(crossout,"Nodes cnodes%d\n",lay+1);
    for (icp=0;icp<numCPoints;icp++)
       filePrint(crossout,"%d\t % 14.6f\t% 14.6f\t % 14.6f\n",
                              icp+1,CPoint[lay][icp][0],CPoint[lay][icp][1],CPoint[lay][icp][2]);


    filePrint(crossout,"Elements cele%d using cnodes%d\n",lay+1,lay+1);
    for (iele=0;iele<numele;iele++) {
       filePrint(crossout,"%d 1 %d %d \n",iele*2+1,iele*4+1,iele*4+2);
       filePrint(crossout,"%d 1 %d %d \n",iele*2+2,iele*4+3,iele*4+4);
    }

    // ... close top-file for pseudo structure

    fclose(crossout);

    // ... print header

    filePrint(oinfo[iInfo].filptr,
            "Vector COMPOSIT_LAYER%d under %s for %s%d\n%d\n",
            lay+1,"Attributes","cnodes",lay+1,4*numele);

    // .... print zero deformation

    filePrint(oinfo[iInfo].filptr,"  %f\n",0.0);
    for (icp=0;icp<numCPoints;icp++)
       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                                    ,w,p,0.0,w,p,0.0,w,p,0.0);

    fflush(oinfo[iInfo].filptr);
  }
  else {

    filePrint(oinfo[iInfo].filptr,"  %f\n",time);

    int iele;
    for(iele=0; iele<numele; ++iele) {

       // ... get E1, E2, Phi of layer & cframe

       double * layData = packedEset[iele]->getCompositeData(lay);
       double * cFrame  = packedEset[iele]->getCompositeFrame();

       if ( ! layData || ! cFrame ) continue;

       double E1  = layData[0];
       double E2  = layData[1];
       double Phi = layData[8];

       double length1 = crossScale * E1 / std::max(E1,E2) / 2.0;
       double length2 = crossScale * E2 / std::max(E1,E2) / 2.0;

       // .... Adjust Phi

       if ( Phi < 0.0 || Phi > 360.0 ) {
         int irot = int(Phi / 360.0);
         Phi = Phi - double (irot) * 360.0;
         if ( Phi < 0.0 ) Phi = 360.0 - Phi;
       }
       Phi = (3.141593*Phi)/180.0;

       // .... determine cross points in cframe

       double * cP1 = new double[3]; double * cP2 = new double[3];
       double * cP3 = new double[3]; double * cP4 = new double[3];

       cP1[0] =  length1 * cos(Phi);
       cP1[1] =  length1 * sin(Phi);
       cP1[2] =  0.0;

       cP2[0] = -length1 * cos(Phi);
       cP2[1] = -length1 * sin(Phi);
       cP2[2] =  0.0;

       cP3[0] =  length2 * cos(Phi+1.5708);
       cP3[1] =  length2 * sin(Phi+1.5708);
       cP3[2] =  0.0;

       cP4[0] = -length2 * cos(Phi+1.5708);
       cP4[1] = -length2 * sin(Phi+1.5708);
       cP4[2] =  0.0;

       // ... transform into global frame

       int n1 = iele*4; int n2 = n1+1; int n3 = n1+2; int n4 = n1+3;

       double x,y,z;

       x = MidPoint[iele][0] + cFrame[0]*cP1[0] +  cFrame[3]*cP1[1] + cFrame[6]*cP1[2] - CPoint[lay][n1][0];
       y = MidPoint[iele][1] + cFrame[1]*cP1[0] +  cFrame[4]*cP1[1] + cFrame[7]*cP1[2] - CPoint[lay][n1][1];
       z = MidPoint[iele][2] + cFrame[2]*cP1[0] +  cFrame[5]*cP1[1] + cFrame[8]*cP1[2] - CPoint[lay][n1][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       x = MidPoint[iele][0] + cFrame[0]*cP2[0] +  cFrame[3]*cP2[1] + cFrame[6]*cP2[2] - CPoint[lay][n2][0];
       y = MidPoint[iele][1] + cFrame[1]*cP2[0] +  cFrame[4]*cP2[1] + cFrame[7]*cP2[2] - CPoint[lay][n2][1];
       z = MidPoint[iele][2] + cFrame[2]*cP2[0] +  cFrame[5]*cP2[1] + cFrame[8]*cP2[2] - CPoint[lay][n2][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       x = MidPoint[iele][0] + cFrame[0]*cP3[0] +  cFrame[3]*cP3[1] + cFrame[6]*cP3[2] - CPoint[lay][n3][0];
       y = MidPoint[iele][1] + cFrame[1]*cP3[0] +  cFrame[4]*cP3[1] + cFrame[7]*cP3[2] - CPoint[lay][n3][1];
       z = MidPoint[iele][2] + cFrame[2]*cP3[0] +  cFrame[5]*cP3[1] + cFrame[8]*cP3[2] - CPoint[lay][n3][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       x = MidPoint[iele][0] + cFrame[0]*cP4[0] +  cFrame[3]*cP4[1] + cFrame[6]*cP4[2] - CPoint[lay][n4][0];
       y = MidPoint[iele][1] + cFrame[1]*cP4[0] +  cFrame[4]*cP4[1] + cFrame[7]*cP4[2] - CPoint[lay][n4][1];
       z = MidPoint[iele][2] + cFrame[2]*cP4[0] +  cFrame[5]*cP4[1] + cFrame[8]*cP4[2] - CPoint[lay][n4][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       fflush(oinfo[iInfo].filptr);

       delete [] cP1;
       delete [] cP2;
       delete [] cP3;
       delete [] cP4;
    }
  }
}

//------------------------------------------------------------------------------

char * Domain::getBasename(char * fname)
{
   char * bname;

   int    ssize  = strlen(fname);
   char * fpoint = strstr(fname,".");
   int    psize  = fpoint - fname;

   if ( psize < ssize )
   {
     bname  = new char[psize+1];
     int i;
     for (i=0;i<psize;i++) bname[i]=fname[i];
     bname[psize] = 0;
   }
   else
   {
     bname = fname;
   }

   return bname;
}

// New functions added for Helmholtz-Axisymmetric

const Connectivity *
Domain::getNodeToNode()
{
	if(nodeToNode) return nodeToNode.get();
	if(elemToNode == 0)
		elemToNode = new Connectivity(packedEset.asSet());
	if(nodeToElem == 0)
		nodeToElem = elemToNode->alloc_reverse();
	nodeToNode = std::make_unique<Connectivity>( nodeToElem->transcon(*elemToNode) );
	return nodeToNode.get();
}


ConstrainedDSA *
Domain::makeCDSA(int nbc, BCond *bcs) {

  if(!c_dsa) c_dsa = new ConstrainedDSA(*dsa, nbc, bcs);
  if(solInfo().HEV) {
    fprintf(stderr," *** c_dsaFluid NOT built in this version of makeCDSA ***\n");
    //if(!c_dsaFluid) c_dsaFluid = new ConstrainedDSA(*dsaFluid, nbcFluid, bcs);
  }

  return c_dsa;
}

void Domain::computeTDProps()
{
  if((numYMTT > 0) || (numCTETT > 0)) {
    // used to calculate temperature dependent material properties

    // compute maximum number of nodes per element
    // note this is also done in Domain::makeAllDOFs()
    maxNumNodes = 0;
    int iele;
    for(iele=0; iele < numele; ++iele) {
      int numNodesPerElement = packedEset[iele]->numNodes();
      maxNumNodes = std::max(maxNumNodes, numNodesPerElement);
    }

    int i;
    // create id mapping for ymtt
    int maxid = 0;
    for(i = 0; i < numYMTT; ++i)
      if(ymtt[i]->getID() > maxid) maxid = ymtt[i]->getID();
    int *ymttmap = new int[maxid + 1];
    for(i = 0; i < numYMTT; ++i) ymttmap[ymtt[i]->getID()] = i;

    // create id mapping for ctett
    maxid = 0;
    for(i = 0; i < numCTETT; ++i)
      if(ctett[i]->getID() > maxid) maxid = ctett[i]->getID();
    int *ctettmap = new int[maxid + 1];
    for(i = 0; i < numCTETT; ++i) ctettmap[ctett[i]->getID()] = i;

    // compute average temperatures and material properties
    int *nodeNumbers = new int[maxNumNodes];
    Vector elemNodeTemps(maxNumNodes);
    elemNodeTemps.zero();
    double *nodalTemperatures = getNodalTemperatures();
    for(iele = 0; iele < numele; ++iele) {
      if((packedEset[iele]->numNodes() > 1) && !packedEset[iele]->isSpring() && !packedEset[iele]->isPhantomElement()) {

        if(packedEset[iele]->getProperty()->E < 0) {
          int id = (int) -packedEset[iele]->getProperty()->E;
          packedEset[iele]->getProperty()->ymtt = ymtt[ymttmap[id]];
        }
        if(packedEset[iele]->getProperty()->W < 0) {
          int id = (int) -packedEset[iele]->getProperty()->W;
          packedEset[iele]->getProperty()->ctett = ctett[ctettmap[id]];
        }

        if((packedEset[iele]->getProperty()->ymtt) || (packedEset[iele]->getProperty()->ctett)) { // iele has temp-dependent E or W
          int NodesPerElement = packedEset[iele]->numNodes();
          packedEset[iele]->nodes(nodeNumbers);

          // compute average temperature in element
          double avTemp = 0.0;
          int iNode;
          for(iNode = 0; iNode < NodesPerElement; ++iNode) {
            if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
              elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
            else
              elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
            avTemp += elemNodeTemps[iNode];
          }
          avTemp /= NodesPerElement;

          StructProp *newProp = new StructProp(*packedEset[iele]->getProperty());
          // compute E using interp table
          if(packedEset[iele]->getProperty()->ymtt) {
            newProp->E = packedEset[iele]->getProperty()->ymtt->getValAlt(avTemp);
          }

          // compute coeff of thermal expansion using interp table
          if(packedEset[iele]->getProperty()->ctett) {
            newProp->W = packedEset[iele]->getProperty()->ctett->getValAlt(avTemp);
          }
          packedEset[iele]->setProp(newProp);
        }
      }
    }
    delete [] ymttmap;
    delete [] ctettmap;
    delete [] nodeNumbers;
  }
}

// ************************************************************************
// HB: Temporary methods for testing tied Mortar condition by generating
//     the LMPC from the Mortar tied conditions
// ************************************************************************
//                           EXPERIMENTAL
// ************************************************************************
// HB: Add Surface Entity to the SurfEntities array
int
Domain::AddSurfaceEntity(SurfaceEntity* _SurfEntity)
{ //--- Verify if _SurfEntity was already defined
 int i=0;
 while (i<nSurfEntity &&
        SurfEntities[i]->ID()!=_SurfEntity->ID())
        i++;
 // if _SurfEntity not previously defined create new
 if (i==nSurfEntity)
  SurfEntities[nSurfEntity++] = _SurfEntity;
 else 
  if(_SurfEntity->ID() != 0) filePrint(stderr," *** WARNING: Surface Entity of Id %d has already been defined !!!\n", SurfEntities[i]->ID());

 return 0;
}

//HB: Add Surface Entity at the position isurf in the SurfEntities array
int
Domain::AddSurfaceEntity(SurfaceEntity* _SurfEntity, int isurf)
{
 //--- Verify if _SurfEntity was already defined
 if(SurfEntities[isurf]) {
   if(SurfEntities[isurf]->ID()!=_SurfEntity->ID()) {
     filePrint(stderr," *** ERROR in Domain::AddSurfaceEntity(SurfaceEntity*, int):\n");
     filePrint(stderr," -> try to overwrite the surface stored at index %d in the SurfEntities array by a different surface\n",isurf);
     exit(-1);
    } else {
      filePrint(stderr," *** WARNING in Domain::AddSurfaceEntity(SurfaceEntity*, int):\n");
      filePrint(stderr," -> surface %d already exist in the SurfEntities array at position %d\n",isurf);
    }
 } else {
   SurfEntities[isurf] = _SurfEntity;
 }
 return 0;
}

// HB: Print informations about the Surface Entity
//     in the SurfEntities array
void Domain::PrintSurfaceEntities()
{
 if(nSurfEntity>0) {
   filePrint(stderr," ... Model has %2d Surface Entities  ...\n",nSurfEntity);
#ifdef MORTAR_DEBUG
   for(int i=0; i<nSurfEntity; i++)
      SurfEntities[i]->Print();
#endif
 } else {
   filePrint(stderr," ... No Surface Entity in the pb    ...\n");
 }
}

// Add a wet surface Id.
int 
Domain::AddAeroEmbedSurfaceId(int Id)
{
  aeroEmbeddedSurfaceId.insert(Id);
  return 0;
}
// HB: Add Mortar Conditions to the MortarConds array
// Warning: CURRENTLY NO CHECK FOR MULTIPLE DEFINITION OF THE SAME MORTAR
//          CONDITION
int
Domain::AddMortarCond(MortarHandler* _MortarCond)
{
  // By default, set the  MortarCond Id = the order in which
  //we read them (see no check !!)
  _MortarCond->SetId(nMortarCond);

  // PJSA dirty fix for self contact
  if(_MortarCond->GetMasterEntityId() == _MortarCond->GetSlaveEntityId()) {
    SurfaceEntity *dummy = new SurfaceEntity(0); 
    AddSurfaceEntity(dummy);
    _MortarCond->SetMasterEntityId(0);
    _MortarCond->SetSelfContact(true);
  }

  // Add the MortarCond to the MortarConds array
  MortarConds[nMortarCond] = _MortarCond;

  // Count number of MortarCond
  nMortarCond++;

  if(_MortarCond->GetInteractionType() == MortarHandler::CTC) nContactSurfacePairs++;

  return 0;
}

void
Domain::DeleteMortarConds()
{
  if(nMortarCond) {
    for(int i=0; i<nMortarCond; ++i)
      if(MortarConds[i]) delete MortarConds[i];
  }
  nMortarCond = 0;
}

// HB: Print informations about Mortar Conditions
//     in the MortarConds array
void
Domain::PrintMortarConds()
{
  if(nMortarCond>0){
   filePrint(stderr," ... Model has %2d Mortar conditions ...\n",nMortarCond);
#ifdef MORTAR_DEBUG
  for(int i=0; i<nMortarCond; i++)
    MortarConds[i]->Print();
#endif
  } else
    filePrint(stderr," ... No Mortar condition in the pb  ...\n");

}

// HB: Set Ptr to Surface entities in the MortarCond objects
void
Domain::SetMortarPairing()
{
  if(nMortarCond>0){
#ifdef MORTAR_DEBUG
  filePrint(stderr," ... Set Ptr to Surface entities in the MortarCond objects ...\n");
#endif
  std::map<int, SurfaceEntity*> SurfIdToPtrSurfMap;
  for(int iSurf=0; iSurf<nSurfEntity; iSurf++){
    SurfIdToPtrSurfMap[SurfEntities[iSurf]->ID()] = SurfEntities[iSurf];
  }

  maxContactSurfElems = 0;
  for(int iMortar=0; iMortar<nMortarCond; iMortar++){
    MortarConds[iMortar]->SetPtrSurfaceEntity(SurfIdToPtrSurfMap[MortarConds[iMortar]->GetMasterEntityId()],
                                              SurfIdToPtrSurfMap[MortarConds[iMortar]->GetSlaveEntityId()]);

#ifdef MORTAR_DEBUG
    // debug: check pairing
    filePrint(stderr," -> Pairing master surf. Id %d with slave surf. Id %d in Mortar cond Id %d\n",
              MortarConds[iMortar]->GetMasterEntityId(),MortarConds[iMortar]->GetSlaveEntityId(),MortarConds[iMortar]->ID());
    filePrint(stderr," -> check throught ptr: master surf Id = %d, slave surf Id = %d\n",
              MortarConds[iMortar]->GetPtrMasterEntity()->ID(),
              MortarConds[iMortar]->GetPtrSlaveEntity()->ID());
#endif

    if(MortarConds[iMortar]->GetInteractionType() == MortarHandler::CTC && !tdenforceFlag())
      maxContactSurfElems += MortarConds[iMortar]->GetPtrSlaveEntity()->GetnNodes();
  }
  }
}

// HB: to setup internal data & renumber surfaces
void Domain::SetUpSurfaces(CoordSet* cs)
{
  for(int iSurf=0; iSurf<nSurfEntity; iSurf++){

    // For acme shells it is necessary to include both sides of the element in the face block
    // currently we only use acme shells when we are using the tdenforcement module of acme
    // to compute the contact forces. In other cases we don't use acme shells since acme doesn't
    // support the extraction of interactions for shells
    if(SurfEntities[iSurf]->GetIsShellFace() && tdenforceFlag()) {
      int nFaceElems = SurfEntities[iSurf]->GetnFaceElems();
      for(int i = 0; i < nFaceElems; ++i) {
        int etype = SurfEntities[iSurf]->GetFaceElemSet()[i]->GetFaceElemType();
        if(etype != 1 && etype != 3) {
          std::cerr << " *** ERROR: Surface element type " << etype << " not supported with SHELL_THICKNESS option\n";
          exit(-1);
        }
        int nNodes = SurfEntities[iSurf]->GetFaceElemSet()[i]->nNodes();
        int *nodes = new int[nNodes];
        for(int j=0; j<nNodes; ++j) nodes[nNodes-1-j] = SurfEntities[iSurf]->GetFaceElemSet()[i]->GetNode(j);
        SurfEntities[iSurf]->AddFaceElement(nFaceElems+i, etype, nNodes, nodes);
        delete [] nodes;
      }
    }

#ifdef MORTAR_DEBUG
    //filePrint(stderr," ------------------------------------------------------------------------\n");
    //filePrint(stderr,"  average normal of face element of surface %2d BEFORE local renumbering\n",SurfEntities[iSurf]->ID());
    //filePrint(stderr," ------------------------------------------------------------------------\n");
    //if(cs) SurfEntities[iSurf]->PrintFaceNormal(cs);
#endif
    SurfEntities[iSurf]->SetUpData(cs);
    SurfEntities[iSurf]->Renumber();
#ifdef MORTAR_DEBUG
  #ifdef HB_NODALNORMAL
    SurfEntities[iSurf]->ComputeNodalNormals(SurfEntities[iSurf]->GetNodeSet());
    SurfEntities[iSurf]->PrintNodalNormals();
  #endif
    filePrint(stderr," ------------------------------------------------------------------------\n");
    if(SurfEntities[iSurf]->IsRenumbered())
      filePrint(stderr,"  average normal of face element of surface %2d AFTER local renumbering\n",SurfEntities[iSurf]->ID());
    else
      filePrint(stderr,"  average normal of face element of surface %2d AFTER\n",SurfEntities[iSurf]->ID());
    filePrint(stderr," ------------------------------------------------------------------------\n");
    SurfEntities[iSurf]->PrintFaceNormal(SurfEntities[iSurf]->GetNodeSet());
    SurfEntities[iSurf]->Print();
#endif
  }
}

// **********************************************************************************************************************
// These functions use ACME's dynamic 2-configuation search algorithm and contact enforcement model for explicit dynamics
void Domain::InitializeDynamicContactSearch(int numSub, SubDomain **sd)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    CurrentMortarCond->SetDistAcme(sinfo.dist_acme);
    CurrentMortarCond->build_search(true, numSub, sd);
    CurrentMortarCond->build_td_enforcement();
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::TIED)
      CurrentMortarCond->set_search_data(2);
    else
      CurrentMortarCond->set_search_data(1);
    CurrentMortarCond->SetNoSecondary(solInfo().no_secondary);
    CurrentMortarCond->set_search_options();
    if(numSub == 0 && !sd) CurrentMortarCond->set_node_constraints(numDirichlet, dbc);
    else CurrentMortarCond->set_node_constraints(numSub, sd);
  }
}

void Domain::UpdateSurfaceTopology(int numSub, SubDomain **sd)
{
  if(newDeletedElements.empty()) return;

  if(!nodeToFaceElem) {
    nodeToFaceElem = new Connectivity * [nSurfEntity];
    for(int iSurf=0; iSurf<nSurfEntity; ++iSurf) {
      Connectivity faceElemToNode(SetAccess<FaceElemSet>{*SurfEntities[iSurf]->GetPtrFaceElemSet()});
      nodeToFaceElem[iSurf] = faceElemToNode.alloc_reverse();
    }
  }

  int fnodes[12];
  for(std::set<int>::iterator it = newDeletedElements.begin(); it != newDeletedElements.end(); ++it) {
    Element *ele = packedEset[geoSource->glToPackElem(*it)];
    int *enodes = ele->nodes();
    for(int iSurf=0; iSurf<nSurfEntity; ++iSurf) {
      std::map<int,int> *GlToLlNodeMap = SurfEntities[iSurf]->GetPtrGlToLlNodeMap();
      for(int iNode=0; iNode < ele->numNodes(); ++iNode) {
        std::map<int,int>::iterator it2 = GlToLlNodeMap->find(enodes[iNode]);
        if(it2 == GlToLlNodeMap->end()) continue;
        int *GlNodeIds = SurfEntities[iSurf]->GetPtrGlNodeIds();
        for(int j=0; j<nodeToFaceElem[iSurf]->num(it2->second); j++) { // loop over the face elements connected to the iNode-th node
          int k = (*nodeToFaceElem[iSurf])[it2->second][j];
          FaceElement *faceEl = SurfEntities[iSurf]->GetFaceElemSet()[k];
          if(faceEl && (faceEl->nNodes() <= ele->numNodes())) {
            faceEl->GetNodes(fnodes, GlNodeIds);
#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && HAS_CXX11_ALL_OF && HAS_CXX11_LAMBDA
            if(std::all_of(fnodes, fnodes+faceEl->nNodes(),
                           [&](int i){return (std::find(enodes,enodes+ele->numNodes(),i)!=enodes+ele->numNodes());})) {
              //std::cerr << "removing face element " << k+1 << " from surface " << SurfEntities[iSurf]->GetId() << std::endl;
              SurfEntities[iSurf]->RemoveFaceElement(k);
              break;
            }
#else
            std::cerr << " *** ERROR: C++11 support required in Domain::UpdateSurfaceTopology().\n"; exit(-1); 
#endif
          }
        }
      }
    }
    delete [] enodes;
  }

  for(int iSurf=0; iSurf<nSurfEntity; ++iSurf) {
    SurfEntities[iSurf]->Reset(&(geoSource->GetNodes()));
    delete nodeToFaceElem[iSurf];
  }
  delete [] nodeToFaceElem; nodeToFaceElem = 0;

  domain->InitializeDynamicContactSearch(numSub, sd);
}

void Domain::UpdateSurfaces(GeomState *geomState, int config_type) // config_type = 1 for current, 2 for predicted
{
  for(int iSurf=0; iSurf<nSurfEntity; iSurf++) {
    SurfEntities[iSurf]->UpdateNodeData(geomState);
  }
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    CurrentMortarCond->set_node_configuration(config_type);
  }
}

void Domain::UpdateSurfaces(DistrGeomState *geomState, int config_type, SubDomain **sd) // config_type = 1 for current, 2 for predicted
{
  if(solInfo().dist_acme != 2) {
    for(int iSurf=0; iSurf<nSurfEntity; iSurf++) {
      SurfEntities[iSurf]->UpdateNodeData(geomState, sd);
    }
  }
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(solInfo().dist_acme == 2) {
      CurrentMortarCond->GetPtrMasterEntity()->UpdateNodeData(geomState, sd);
      CurrentMortarCond->GetPtrSlaveEntity()->UpdateNodeData(geomState, sd);
    }
    CurrentMortarCond->set_node_configuration(config_type, geomState->getNumSub(), sd);
  }
}

void Domain::PerformDynamicContactSearch(double dt_old, double dt)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    int search_algorithm = (sinfo.old_dynamic_search) ? 3 : 4;
    double dt_old = dt;
    CurrentMortarCond->perform_search(search_algorithm, dt_old, dt);
  }
}

void Domain::AddContactForces(double dt_old, double dt, Vector &f)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->compute_td_contact_force(dt_old, dt, f);
    }
  }
}

void Domain::AddContactForces(double dt_old, double dt, DistrVector &f)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->compute_td_contact_force(dt_old, dt, f);
    }
  }
}

void Domain::MakeNodalMass(SparseMatrix *M, SparseMatrix *Mcc)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->make_nodal_mass(M, c_dsa, Mcc);
      CurrentMortarCond->make_kinematic_partitioning(packedEset, nodeToElem);
    }
  }
}

#include <Paral.d/SubDOp.h>
void Domain::MakeNodalMass(SubDOp *M, SubDomain **sd)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->make_nodal_mass(M, sd);
      CurrentMortarCond->make_kinematic_partitioning(M->getNumSub(), sd);
    }
  }
}
// **********************************************************************************************************************

// These functions use ACME's static 1-configuation search algorithm and face-face interaction, and Henri's mortar LMPCs for statics and implicit dynamics
void Domain::InitializeStaticContactSearch(MortarHandler::Interaction_Type t, int numSub, SubDomain **sd)
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      CurrentMortarCond->SetDistAcme(sinfo.dist_acme);
      CurrentMortarCond->SetMortarScaling(sinfo.mortar_scaling);
      CurrentMortarCond->SetMortarIntegrationRule(sinfo.mortar_integration_rule);
      CurrentMortarCond->build_search(false, numSub, sd);
      CurrentMortarCond->set_search_data(4); // interaction_type = 4 (FaceFace) 
      CurrentMortarCond->set_node_configuration(1);
      CurrentMortarCond->SetNoSecondary(solInfo().no_secondary);
      CurrentMortarCond->set_search_options();
      if(numSub == 0 && !sd) CurrentMortarCond->set_node_constraints(numDirichlet, dbc);
      else {
        CurrentMortarCond->set_node_constraints(numSub, sd);
        CurrentMortarCond->make_share(numSub, sd);
      }
    }
  }
}

void Domain::ReInitializeStaticContactSearch(MortarHandler::Interaction_Type t, int numSub, SubDomain **sd)
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      if(sd && numSub != 0) CurrentMortarCond->make_share(numSub, sd);
    }
  }
}

void Domain::UpdateSurfaces(MortarHandler::Interaction_Type t, GeomState *geomState)
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      CurrentMortarCond->GetPtrMasterEntity()->UpdateNodeData(geomState);
      CurrentMortarCond->GetPtrSlaveEntity()->UpdateNodeData(geomState);
      CurrentMortarCond->set_node_configuration(1); // 1 --> config_type current
    }
  }
}

void Domain::UpdateSurfaces(MortarHandler::Interaction_Type t, DistrGeomState *geomState, SubDomain **sd) 
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      CurrentMortarCond->GetPtrMasterEntity()->UpdateNodeData(geomState, sd);
      CurrentMortarCond->GetPtrSlaveEntity()->UpdateNodeData(geomState, sd);
      CurrentMortarCond->set_node_configuration(1, geomState->getNumSub(), sd); // 1 --> config_type current
    }
  }
}

void Domain::PerformStaticContactSearch(MortarHandler::Interaction_Type t)
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      int search_algorithm = 1; // static 1-configuration
      CurrentMortarCond->perform_search(search_algorithm);
      int interaction_type = 4;
      CurrentMortarCond->get_interactions(interaction_type);
    }
  }
}

void Domain::ExpComputeMortarLMPC(MortarHandler::Interaction_Type t, int nDofs, int *dofs)
{
  int num_interactions = 0;
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      CurrentMortarCond->CreateFFIPolygon();
      CurrentMortarCond->AddMortarLMPCs(&lmpc, numLMPC, numCTC, nDofs, dofs);
      nMortarLMPCs += CurrentMortarCond->GetnMortarLMPCs();
      num_interactions += CurrentMortarCond->GetnFFI();
      CurrentMortarCond->DeleteFFIData();
    }
  }
  if(verboseFlag && nMortarCond > 0) filePrint(stderr," ... Built %d Mortar Surface/Surface Interactions ...\n", nMortarLMPCs);
#ifdef HB_ACME_FFI_DEBUG
  if(sinfo.ffi_debug && num_interactions > 0) {
    char fname[16];
    sprintf(fname,"FFI.top.%d",totalNewtonIter);
    filePrint(stderr," -> Write FFI top file: %s for %d interactions\n", fname, num_interactions);
    FILE* FFITopFile = fopen(fname,"w");
    WriteFFITopFile(FFITopFile);
    fclose(FFITopFile);
  }
#endif
}

// HB: compute & add to the standard LMPC array (lmpc) the Mortar LMPCs:
void Domain::ComputeMortarLMPC(int nDofs, int *dofs)
{
  if(nMortarCond>0){
  double time00=-getTime();
  for(int iMortar=0; iMortar<nMortarCond; iMortar++){
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
#ifdef MORTAR_TIMINGS
    double time0= -getTime();
    double time = -getTime();
#endif
#ifdef MORTAR_DEBUG
    filePrint(stderr," ... Treat Mortar condition Id %d\n",CurrentMortarCond->ID());
    CurrentMortarCond->Print();
#endif
    switch(int(CurrentMortarCond->GetGeomType())){
     case MortarHandler::EQUIVALENCED:
       CurrentMortarCond->CreateACMEFFIData();
       break;
     case MortarHandler::NON_MATCHING:
     default:
       CurrentMortarCond->PerformACMEFFISearch();
       break;
    }
#ifdef MORTAR_TIMINGS
    time += getTime();
    filePrint(stderr," -> time spent in the ACME search: %e s\n",time/1000);
    time = -getTime();
#endif
    CurrentMortarCond->CreateFFIPolygon();
#ifdef MORTAR_TIMINGS
    time += getTime();
    filePrint(stderr," -> time spent in building the Mortar B matrices: %e s\n",time/1000);
#endif
    if(CurrentMortarCond->GetInteractionType()!=MortarHandler::FSI){
#ifdef MORTAR_TIMINGS
      time = -getTime();
#endif
      CurrentMortarCond->AddMortarLMPCs(&lmpc, numLMPC, numCTC, nDofs, dofs);
#ifdef MORTAR_TIMINGS
      time += getTime();
      time0 += getTime();
      filePrint(stderr," -> time spent in creating the Mortar LMPCs: %e s\n",time/1000);
      filePrint(stderr," -> total time spent in building those Mortar LMPCs: %e s\n",time0/1000);
#endif
      // count the total number of Mortar LMPCs generated
      nMortarLMPCs += CurrentMortarCond->GetnMortarLMPCs();
    } else
      CurrentMortarCond->AddWetFSI(&fsi, numFSI);
  }

  time00 += getTime();
  if(verboseFlag && nMortarCond > 0) filePrint(stderr," ... Built %d Mortar Surface/Surface Interactions ...\n", nMortarLMPCs+numFSI);
#ifdef MORTAR_TIMINGS
  filePrint(stderr," ... CPU time for building mortar surface/surface interactions: %e s\n",time00/1000);
#endif
 }
#ifdef HB_ACME_FFI_DEBUG
 filePrint(stderr," -> Write FFI top file: FFI.top\n");
 FILE* FFITopFile = fopen("FFI.top","w");
 WriteFFITopFile(FFITopFile); //valid only if cs!=0 ...
 fclose(FFITopFile);
#endif
}

#ifdef HB_ACME_FFI_DEBUG
void Domain::WriteFFITopFile(FILE* file)
{
  if(nMortarCond){
    int firstNodeId = 1;
    fprintf(file,"Nodes FFINodes\n");
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* MortarCond = MortarConds[iMortar];
      MortarCond->PrintFFIPolyVertices(file, firstNodeId);
    }
    firstNodeId = 1;
    int firstElemId = 1;
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* MortarCond = MortarConds[iMortar];
      fprintf(file,"Elements SlaveFFI_%d using FFINodes\n",MortarCond->GetSlaveEntityId());
      MortarCond->PrintFFIPolyTopo(file, firstElemId, firstNodeId);
      fprintf(file,"Elements MasterFFI_%d using FFINodes\n",MortarCond->GetMasterEntityId());
      MortarCond->PrintFFIPolyTopo(file, firstElemId, firstNodeId);
    }
  }
}
#endif

// HB: create (Tied) MortarToMPC connectivity
// WARNING: (1) WE ASSUME THAT THE MORTAR LMPCs ARE AFTER THE REGULAR LMPCs
//              IN THE STANDARD LMPC ARRAY (lmpc)
//          (2) WE ASSUME THAT IN EACH MORTAR CONDITION, THE MORTAR LMPCs
//              ARE NUMBERED SEQUENTIALY
// We currently do not separate a Tied Mortar LMPC into 2 or 3 (2D/3D)
// Tied Mortar LMPC (one for each axis: x, y, z)
// this splitting can be done at the CCt level (see component renumbering)
// DONE: changed the order in which the Mortar LMPCs are added in the
// standard LMPC array (lmpc):
// before    [(x,y,z),...,(x,y,z)]
// currently [(x,...,x),(y,...,y),(z,...,z)]
void
Domain::CreateMortarToMPC()
{
  int* pointermap = 0;
  int* target     = 0;

  if(mortarToMPC) { delete mortarToMPC; mortarToMPC = 0; }

  if(nMortarCond>0){
#ifdef MORTAR_DEBUG
    filePrint(stderr," ... Create (Tied) Mortar To LMPCs connectivity ...\n");
#endif

    // 1) Set the map pointer array & count number of
    //    targets = total number of (Tied) Mortar LMPCs
    pointermap = new int[nMortarCond+1];
    int ntargets = 0;
    pointermap[0] = 0;
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* CurrentMortarCond = MortarConds[iMortar];
      int nLMPCs = CurrentMortarCond->GetnMortarLMPCs();
      pointermap[iMortar+1] = pointermap[iMortar] + nLMPCs;
      ntargets += nLMPCs;
    }

    // 2) Fill the target (map) array
    target = new int[ntargets];
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* CurrentMortarCond = MortarConds[iMortar];
      int nLMPCs       = CurrentMortarCond->GetnMortarLMPCs();
      int IdFirstLMPC  = CurrentMortarCond->GetIdFirstMortarLMPC();
      for(int i=0; i<nLMPCs; i++)
        target[pointermap[iMortar]+i] = IdFirstLMPC+i;
    }

    mortarToMPC = new Connectivity(nMortarCond, pointermap, target);
  }
}
// HB: return the Mortar to LMPCs connectivity and the total number of
// Mortar LMPCs
Connectivity*
Domain::GetMortarToMPC() { return mortarToMPC; }
int
Domain::GetnMortarLMPCs() { return nMortarLMPCs; }

// HB: write in a file ALL the generated Mortar LMPCs
// make the assumptions that the Mortar LMPC id is NEGATIVE (<0)
// basicaly, loop over ALL the lmpc and look for the negative lmpc ids
// Another possibility will be to use the mortarToMPC connectivity
// !!! TO DO !!!
/*
void Domain::WriteToFileMortarLMPCs(FILE *file)
{

  int i,ilmpc, nmotarlmpc=0;
  int firstmortarlmpcId = numLMPC-nMortarLMPCs;
  if(firstmortarlmpcId<0) firstmortarlmpcId = 0;
  for(ilmpc=0; ilmpc<numLMPC; ilmpc++){
    int lmpcId = lmpc[ilpmc].lmpcnum;
    if(lmpcId<0){
     nmotarlmpc++;
     double rhs = lmpc[ilpmc].rhs;
     lmpcId = firstmortarlmpcId+nmotarlmpc;
     filePrint(file," %4d  %f\n",lmpcId,rhs);
     for(i=0;i<lmpc[ilpmc].nterms;i++){
       int node  = lmpc[ilpmc].terms[i].nnum;
       int dof   = lmpc[ilpmc].terms[i].dofnum;
       double val= lmpc[ilpmc].terms[i].val;
       filePrint(file,"    %4d  %d   %f\n",node,dof,val);
     }
    }
  }
  filePrint(stderr," -> has written %d motar lmpcs to file\n",nmotarlmpc);
}
*/

// ************************************************************************

// returns the value of the pressure force flag
int
Domain::pressureFlag() { return geoSource->pressureFlag(); }

// function that returns composite layer info
LayInfo *Domain::getLayerInfo(int num) { return geoSource->getLayerInfo(num); }

void
Domain::initializeNumbers()
{
 numThicknessGroups = 0; numShapeVars = 0;  thgreleFlag = 0;  thpaIndex = 0;
 senInfo = new SensitivityInfo[50];  // maximum number of sensitivities are fixed to 50
 aggregatedStress = new double;
 aggregatedStressDenom = new double;
 *aggregatedStress = 0;
 *aggregatedStressDenom = 0;
 numStressNodes = stressNodes.size();
 numDispNodes = dispNodes.size();
}

void
Domain::initialize()
{
 numdofs = 0; numDispDirichlet = 0; numContactPairs = 0;
 numIDis = 0; numIDisModal = 0; numIVel = 0; numIVelModal = 0; numDirichlet = 0; numNeuman = 0; numSommer = 0;
 numComplexDirichlet = 0; numComplexNeuman = 0; numNeumanModal = 0;
 firstDiMass = 0; numIDis6 = 0; gravityAcceleration = 0;
 allDOFs = 0; stress = 0; weight = 0; elstress = 0; elweight = 0; claw = 0; com = 0;
 numLMPC = 0; numYMTT = 0; numCTETT = 0; numSDETAFT = 0; numRUBDAFT = 0; numSS1DT = 0; numSS2DT = 0; numYSST = 0; numYSSRT = 0; numYMST = 0;
 MidPoint = 0; temprcvd = 0;
 heatflux = 0; elheatflux = 0; elTemp = 0; dbc = 0; nbc = 0; nbcModal = 0;
 iDis = 0; iDisModal = 0; iVel = 0; iVelModal = 0; iDis6 = 0; elemToNode = 0; nodeToElem = 0;
 dsa = 0; c_dsa = 0; cdbc = 0; cnbc = 0;
 dsaFluid = 0; c_dsaFluid = 0; allDOFsFluid = 0; dbcFluid = 0;
 elemToNodeFluid = 0; nodeToElemFluid = 0; nodeToNodeFluid = 0;
 nSurfEntity = 0; nMortarCond = 0; nMortarLMPCs= 0; mortarToMPC = 0;
 solver = 0; csolver = 0;
 flExchanger = 0; outFile = 0; elDisp = 0; p_stress = 0; p_elstress = 0; stressAllElems = 0;
 previousExtForce = 0; previousAeroForce = 0; previousDisp = 0; previousCq = 0;
 temprcvd = 0; optinputfile = 0;
 nSurfEntity = 0; nMortarLMPCs = 0; mortarToMPC = 0; matrixTimers = 0;
 allDOFs = 0; haveNodes = false; nWetInterface = 0; wetInterfaces = 0;
 numFSI = 0; firstOutput = true; nodeToNode_sommer = 0; sowering = false; nDimass = 0;
 maxNumDOFs = 0; maxNumDOFsFluid = 0;
 fluidDispSlosh = 0; elPotSlosh = 0; elFluidDispSlosh = 0;
 elFluidDispSloshAll = 0;
 nodeToFsi = 0;
 numCTC = 0;
 output_match_in_top = false;
 C_condensed = 0;
 nContactSurfacePairs = 0; maxContactSurfElems = 0;
 nodeToFaceElem = 0;
 outFlag = 0;
 nodeTable = 0;
 MpcDSA = 0; nodeToNodeDirect = 0;
 g_dsa = 0;
 numSensitivity = 0; senInfo = 0; aggregatedStress = 0; aggregatedStressDenom = 0;
 runSAwAnalysis = false;   
 aggregatedFlag = false;
 aggregatedFileNumber = 0;
 numSensitivityQuantityTypes = 0;
 numTotalDispDofs = 0;
}

Domain::~Domain()
{
 if(nodeToNodeFluid) { delete nodeToNodeFluid; nodeToNodeFluid = 0; }
 if(renumb.order) { delete [] renumb.order; renumb.order = 0; }
 if(renumb.renum) { delete [] renumb.renum; renumb.renum = 0; }
 if(renumb.xcomp) { delete [] renumb.xcomp; renumb.xcomp = 0; }
 if(renumbFluid.order) { delete [] renumbFluid.order; renumbFluid.order = 0; }
 if(renumbFluid.renum) { delete [] renumbFluid.renum; renumbFluid.renum = 0; }
 if(gravityAcceleration) { delete [] gravityAcceleration; gravityAcceleration = 0; }
 if(matrixTimers) { delete matrixTimers; matrixTimers = 0; }
 if(allDOFs) { delete allDOFs; allDOFs = 0; }
 if(dsa) { delete dsa; dsa = 0; }
 if(c_dsa) { delete c_dsa; c_dsa = 0; }
 if(elemToNode) { delete elemToNode; elemToNode = 0; }
 if(nodeToElem) { delete nodeToElem; nodeToElem = 0; }
 if(allDOFsFluid) { delete allDOFsFluid; allDOFsFluid = 0; }
 if(dsaFluid) { delete dsaFluid; dsaFluid = 0; }
 if(c_dsaFluid) { delete c_dsaFluid; c_dsaFluid = 0; }
 if(elemToNodeFluid) { delete elemToNodeFluid; elemToNodeFluid = 0; }
 if(nodeToElemFluid) { delete nodeToElemFluid; nodeToElemFluid = 0; }
 if(mortarToMPC) { delete mortarToMPC; mortarToMPC = 0; }
 if(dbc) { delete [] dbc; dbc = 0; }
 if(nbc) { delete [] nbc; nbc = 0; }
 if(nbcModal) { delete [] nbcModal; nbcModal = 0; }
 if(iDis) { delete [] iDis; iDis = 0; }
 if(iDisModal) { delete [] iDisModal; iDisModal = 0; }
 if(iVel) { delete [] iVel; iVel = 0; }
 if(iVelModal) { delete [] iVelModal; iVelModal = 0; }
 if(iDis6) { delete [] iDis6; iDis6 = 0; }
 if(cdbc) { delete [] cdbc; cdbc = 0; }
 if(cnbc) { delete [] cnbc; cnbc = 0; }
 if(stress) { delete stress; stress = 0; }
 if(weight) { delete weight; weight = 0; }
 if(elDisp) { delete elDisp; elDisp = 0; }
 if(elTemp) { delete elTemp; elTemp = 0; }
 if(fluidDispSlosh) { delete fluidDispSlosh; fluidDispSlosh = 0; }
 if(elPotSlosh) { delete elPotSlosh; elPotSlosh = 0; }
 if(elFluidDispSlosh) { delete elFluidDispSlosh; elFluidDispSlosh = 0; }
 if(elstress) { delete elstress; elstress = 0; }
 if(elweight) { delete elweight; elweight = 0; }
 if(p_stress) { delete p_stress; p_stress = 0; }
 if(p_elstress) { delete p_elstress; p_elstress = 0; }
 if(stressAllElems) { delete stressAllElems; stressAllElems = 0; }
 if(claw) {
   if(claw->actuator) { delete [] claw->actuator; claw->actuator = 0; }
   if(claw->userForce) { delete [] claw->userForce; claw->userForce = 0; }
   if(claw->userDisp) { delete [] claw->userDisp; claw->userDisp = 0; }
   delete claw; claw = 0;
 }
 if(com) { delete com; com = 0; }
 if(firstDiMass) { delete firstDiMass; firstDiMass = 0; }
 if(previousExtForce) { delete previousExtForce; previousExtForce = 0; }
 if(previousAeroForce) { delete previousAeroForce; previousAeroForce = 0; }
 if(previousDisp) { delete previousDisp; previousDisp = 0; }
 if(previousCq) { delete previousCq; previousCq = 0; }
 if(heatflux) { delete heatflux; heatflux = 0; }
 if(elheatflux) { delete elheatflux; elheatflux = 0; }
 delete &nodes;
 if(wetInterfaces) { delete wetInterfaces; wetInterfaces = 0; }
 if(nodeToNode_sommer) { delete nodeToNode_sommer; nodeToNode_sommer = 0; }
 if(nodeToFsi) { delete nodeToFsi; nodeToFsi = 0; }
 for(int i=0; i<SurfEntities.max_size(); i++)
   if(SurfEntities[i]) { SurfEntities[i]->~SurfaceEntity(); SurfEntities[i] = 0; }
 for(int i=0; i<nMortarCond; ++i)
   if(MortarConds[i]) delete MortarConds[i];
 for(int i=0; i<numLMPC; ++i)
   if(lmpc[i]) delete lmpc[i];
 if(nodeTable) delete [] nodeTable;
 if(flExchanger) delete flExchanger;
 if(MpcDSA) delete MpcDSA; if(nodeToNodeDirect) delete nodeToNodeDirect;
 for(int i=0; i<contactSurfElems.size(); ++i)
   packedEset.deleteElem(contactSurfElems[i]);
 if(g_dsa) delete g_dsa;
 if(senInfo) delete [] senInfo;
 if(aggregatedStress) delete aggregatedStress;
 if(aggregatedStressDenom) delete aggregatedStressDenom;
 if(thgreleFlag) delete [] thgreleFlag;
 if(thpaIndex) delete [] thpaIndex;
}

#include <Element.d/Helm.d/HelmElement.h>
#include <Corotational.d/utilities.h>

int
Domain::isFluidElement(int i)
{
 HelmElement *he = dynamic_cast<HelmElement *>(packedEset[i]);
 if(he==0) return 0; // non-Helmholtz element found
 else return he->isFluid();
}

int
Domain::isStructureElement(int i)
{
  return (!isFluidElement(i) && !packedEset[i]->isConstraintElement());
}

bool
Domain::isHomogeneous()
{
  if((geoSource->getNumProps() == 1) && (numYMTT == 0) && (numCTETT == 0)) return true;
  else return false;
}

void
Domain::addWetInterface(int _fluidSurfaceID, int _structureSurfaceID, double normal_tol, double tangential_tol)
{
  if(nWetInterface == 0) wetInterfaces = new ResizeArray<WetInterface *>(0, 1);  // initial size is 1
  WetInterface *newWI = new WetInterface();
  newWI->fluidSurfaceID = _fluidSurfaceID;
  newWI->structureSurfaceID = _structureSurfaceID;
  (*wetInterfaces)[nWetInterface++] = newWI;

  // PJSA: temporary fix, use mortar code to build LMPCs to make interface
  MortarHandler* _MortarCond = new MortarHandler(_structureSurfaceID, _fluidSurfaceID, normal_tol, tangential_tol);
  _MortarCond->SetId(nMortarCond);
  _MortarCond->SetMortarType(MortarHandler::STD); //HB: trigger std mortar basis
  _MortarCond->SetInteractionType(MortarHandler::FSI);
  //filePrint(stderr," _fluidSurfaceID = %d, _structureSurfaceID = %d\n",_fluidSurfaceID,_structureSurfaceID); //HB
  if(_fluidSurfaceID == _structureSurfaceID)
    _MortarCond->SetGeomType(MortarHandler::EQUIVALENCED);
  MortarConds[nMortarCond] = _MortarCond;
  //MortarConds[nMortarCond]->Print();
  nMortarCond++;

  return;
}

//HB: I think it may be more safe to get them from the fsi array in case ACME misses
//some interactions (non-matching wet interfaces) -> in which case only some of the nodes of the
//wet surfaces given in the input file are currently trully "wet interface nodes".
int *
Domain::getAllWetInterfaceNodes(int &count)
{
  count = 0;
#define USE_FSI_NODES
  if(numFSI) {
#ifdef USE_FSI_NODES //HB: assume the FSI has already been computed. to be validated ...
    int *nodeMap = new int[numnodes];
    for(int i=0; i<numnodes; ++i) nodeMap[i] = -1;
    for(int i=0; i<numFSI; i++) {
      if(nodeMap[fsi[i]->lmpcnum]<0) nodeMap[fsi[i]->lmpcnum] = count++; // fluid node
      for(int j=0; j<fsi[i]->nterms; j++) // Loops over structure's nodes
        if(nodeMap[fsi[i]->terms[j].nnum]<0) nodeMap[fsi[i]->terms[j].nnum] = count++;
    }
#else
  if(nWetInterface > 0) {
    std::map<int, SurfaceEntity*> SurfIdToPtrSurfMap;
    for(int i=0; i<nSurfEntity; i++){
      int SurfId = SurfEntities[i]->ID();
      SurfIdToPtrSurfMap[SurfId] = SurfEntities[i];
    }

    int *nodeMap = new int[numnodes];
    for(int i=0; i<numnodes; ++i) nodeMap[i] = -1;
    for(int i=0; i<nWetInterface; i++) {
      int fluidId = (*wetInterfaces)[i]->fluidSurfaceID;
      int structureId = (*wetInterfaces)[i]->structureSurfaceID;
      // fluid
      SurfaceEntity* fluidEntity = SurfIdToPtrSurfMap[fluidId];
      int nFluidElem = fluidEntity->nFaceElements();
      FaceElemSet *fluidElemSet = fluidEntity->GetPtrFaceElemSet();
      int *glFluidNodes = fluidEntity->GetPtrGlNodeIds();
      int nFluidNodes = fluidEntity->GetnNodes();
      for(int j=0; j<nFluidNodes; ++j)
        if(nodeMap[glFluidNodes[j]] == -1) nodeMap[glFluidNodes[j]] = count++;
      if(structureId == fluidId) continue;
      // structure
      SurfaceEntity* structureEntity  = SurfIdToPtrSurfMap[structureId] ;
      int *glStructureNodes = structureEntity->GetPtrGlNodeIds();
      int nStructureNodes = structureEntity->GetnNodes();
      for(int j=0; j<nStructureNodes; ++j)
        if(nodeMap[glStructureNodes[j]] == -1) nodeMap[glStructureNodes[j]] = count++;
      }
    }
#endif
    glWetNodeMap = nodeMap;
    int *allWetInterfaceNodes = new int[count];
    for(int i=0; i<numnodes; ++i)
      if(nodeMap[i] != -1) allWetInterfaceNodes[nodeMap[i]] = i;

    //if(nodeMap) delete [] nodeMap;
    return allWetInterfaceNodes;
  }
  else {
    //std::cerr << " *** WARNING: No Wet Interfaces have been defined \n";
    return 0;
  }
}

void Domain::printFSI(FILE* file)
{
 filePrint(file," ... FSI list :                    ...\n");
 int i;
 for(i=0; i<numFSI; i++) {
   filePrint(file," %4d ... fluid node %6d : \n", i+1,fsi[i]->lmpcnum+1);
   if(fsi[i]->isComplex) {
     filePrint(file," ... structure node   dof       coef\n");
     for(int j=0; j<fsi[i]->nterms; j++)
       filePrint(file,"        %6d            %d        (%6e,%6e)\n",
                 fsi[i]->terms[j].nnum+1, fsi[i]->terms[j].dofnum,
                 fsi[i]->terms[j].coef.c_value.real(), fsi[i]->terms[j].coef.c_value.imag());
   }
   else {
     filePrint(stderr," ... structure node   dof       coef\n");
     for(int j=0; j<fsi[i]->nterms; j++)
       filePrint(file,"        %6d            %d        %6e\n",
                 fsi[i]->terms[j].nnum+1, fsi[i]->terms[j].dofnum,
                 fsi[i]->terms[j].coef.r_value);
   }
 }
}

double
Domain::getFrequencyOrWavenumber()
{
  double ret = 0.0;
  if(geoSource->isShifted()) {
    if(domain->solInfo().doFreqSweep) {
      if(isCoarseGridSolve) {
        ret = domain->coarse_frequencies->front();
        ret = getSavedFreq();
      }
      else ret = domain->frequencies->front();
    }
    else {
      ret = geoSource->omega();
    }
    ret /= (2.0*PI);
  }
  return ret;
}

void
Domain::computeAverageProps(int &structure_element_count, int &fluid_element_count, double &global_average_E,
                            double &global_average_nu, double &global_average_rhof)
{
  int nEle = packedEset.last();
  for(int i = 0; i < nEle; ++i) {
    Element *ele = packedEset[i];
    StructProp *prop = ele->getProperty();
    if(prop != NULL && !ele->isConstraintElement()) {
      if(! dynamic_cast<HelmElement *>(ele)) { // not a fluid element
        global_average_E += prop->E;
        global_average_nu += prop->nu;
        structure_element_count++;
      } else {
        global_average_rhof += prop->rho;
        fluid_element_count++;
      }
    }
  }
}

void
Domain::computeCoupledScaleFactors()
{
  // use global average properties to compute Lame constants & coupled scaling factor
  // (equation proposed by Jan Mandel)
  double ymod = geoSource->global_average_E;
  double prat = geoSource->global_average_nu;
  //double rhofavg = geoSource->global_average_rhof;
  double lame_lambda = prat*ymod/((1.0+prat)*(1.0-2.0*prat));
  double lame_mu_times2 = 2.0*ymod/(2.0*(1.0+prat));
  double lame_max = (lame_lambda > lame_mu_times2) ? lame_lambda : lame_mu_times2;
/* old scaling
  double rhoFomega2 = (geoSource->isShifted()) ? domain->fluidDensity * geoSource->shiftVal() : 1000.0;
  coupledScaling = ((lame_max > 0.0) && (rhoFomega2 > 0.0)) ? 1.0/sqrt(rhoFomega2*lame_max) : 1.0;
  coupledScaling *= domain->solInfo().coupled_scale; // PJSA 10-21-05 can be used to adjust scaling up or down
  cscale_factor = rhoFomega2*coupledScaling;
*/
  double omega2 = (geoSource->isShifted()) ? geoSource->shiftVal() : 1000.0;
  coupledScaling = ((lame_max > 0.0) && (geoSource->global_average_rhof> 0.0)) ?
               1.0/sqrt(geoSource->global_average_rhof*omega2*lame_max) : 1.0; // RADEK
  coupledScaling *= domain->solInfo().coupled_scale; // PJSA 10-21-05 can be used to adjust scaling up or down
  cscale_factor = omega2*coupledScaling; // RADEK

  cscale_factor2 = cscale_factor*coupledScaling;
  //std::cerr << "coupledScaling = " << coupledScaling << ", cscale_factor = " << cscale_factor
  //     << ", cscale_factor2 = " << cscale_factor2 << std::endl;
}

void
Domain::initSfem()
{
#ifndef SALINAS
  if(domain->solInfo().noninpc || domain->solInfo().inpc) {
 //   geoSource->printGroups();
    sfem->computeLP();
  }
#endif
}

// add nodal contact mpcs
int
Domain::addNodalCTC(int n1, int n2, double nx, double ny, double nz,
                    double normalGap, int _mode, int lagrangeMult, double penalty)
{
 // contact note: if normalGapPresent is false perhaps the default should be to compute the geometric gap
 // using the nodal coordinates
 int mode = (_mode > -1) ? _mode : domain->solInfo().contact_mode;  // 0 -> normal tied + tangents free, 1 -> normal contact + tangents free
                                                                    // 2 -> normal+tangents tied, 3 -> normal contact + tied tangents
 int lmpcnum = 0;

 // normal constraint
 LMPCons *_CTC = new LMPCons(lmpcnum, normalGap);
 double norm = sqrt(nx*nx + ny*ny + nz*nz);
 if(nx != 0.0) {
   nx /= norm;
   LMPCTerm *term1 = new LMPCTerm(n1, 0, nx);
   _CTC->addterm(term1);
   LMPCTerm *term2 = new LMPCTerm(n2, 0, -nx);
   _CTC->addterm(term2);
 }
 if(ny != 0.0) {
   ny /= norm;
   LMPCTerm *term1 = new LMPCTerm(n1, 1, ny);
   _CTC->addterm(term1);
   LMPCTerm *term2 = new LMPCTerm(n2, 1, -ny);
   _CTC->addterm(term2);
 }
 if(nz != 0.0) {
   nz /= norm;
   LMPCTerm *term1 = new LMPCTerm(n1, 2, nz);
   _CTC->addterm(term1);
   LMPCTerm *term2 = new LMPCTerm(n2, 2, -nz);
   _CTC->addterm(term2);
 }
 _CTC->type = (mode == 1 || mode == 3); // this is to be phased out
 if(mode == 1 || mode == 3) _CTC->setType(mpc::Inequality);
 _CTC->setSource(mpc::NodalContact);
 addLMPC(_CTC,false);

 // PJSA 7-12-2007 tangential constraints ... note this may lead to redundant constraints (singularity in CCt)
 // for 2D you need to set spacedim 2 in the input file, also it is currently assumed that 2D model is defined in XY plane
 if(mode == 2 || mode == 3) {
   // 1. select direction for initial cross product... use z unless normal is parallel or close to parallel to z axis, in which case use y
   double n[3] = { nx, ny, nz };
   double t1[3], t2[3];
   if(!(nx < 0.1 && ny < 0.1) || solInfo().solvercntl->fetiInfo.spaceDimension == 2) { t2[0] = 0.0; t2[1] = 0.0; t2[2] = 1.0; }
   else if(!(nx < 0.1 && nz < 0.1)) { t2[0] = 0.0; t2[1] = 1.0; t2[2] = 0.0; }
   else { t2[0] = 1.0; t2[1] = 0.0; t2[2] = 0.0; }
   crossprod(n,t2,t1);
   normalize(t1);
   LMPCons *_TGT1 = new LMPCons(0, 0.0);
   if(t1[0] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 0, t1[0]); _TGT1->addterm(term1);
                                    LMPCTerm *term2 = new LMPCTerm(n2, 0, -t1[0]); _TGT1->addterm(term2); }
   if(t1[1] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 1, t1[1]); _TGT1->addterm(term1);
                                    LMPCTerm *term2 = new LMPCTerm(n2, 1, -t1[1]); _TGT1->addterm(term2); }
   if(t1[2] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 2, t1[2]); _TGT1->addterm(term1);
                                    LMPCTerm *term2 = new LMPCTerm(n2, 2, -t1[2]); _TGT1->addterm(term2); }
   _TGT1->setSource(mpc::NodalContact);
   addLMPC(_TGT1,false);
   if(solInfo().solvercntl->fetiInfo.spaceDimension == 3) {
     crossprod(n,t1,t2);
     normalize(t2);
     LMPCons *_TGT2 = new LMPCons(0, 0.0);
     if(t2[0] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 0, t2[0]); _TGT2->addterm(term1);
                                      LMPCTerm *term2 = new LMPCTerm(n2, 0, -t2[0]); _TGT2->addterm(term2); }
     if(t2[1] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 1, t2[1]); _TGT2->addterm(term1);
                                      LMPCTerm *term2 = new LMPCTerm(n2, 1, -t2[1]); _TGT2->addterm(term2); }
     if(t2[2] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 2, t2[2]); _TGT2->addterm(term1);
                                      LMPCTerm *term2 = new LMPCTerm(n2, 2, -t2[2]); _TGT2->addterm(term2); }
     _TGT2->setSource(mpc::NodalContact);
     addLMPC(_TGT2,false);
   }
 }

 return numCTC;
}

int
Domain::getNumCTC()
{
  return numCTC + geoSource->getNumConstraintElementsIeq();
}

void
Domain::addDirichletLMPCs(int _numDirichlet, BCond *_dbc)
{
  for(int i=0; i<_numDirichlet; ++i) {
    LMPCons *dmpc = new LMPCons(0, _dbc[i].val,
         new LMPCTerm(_dbc[i].nnum,_dbc[i].dofnum,1.0));
    addLMPC(dmpc,false);
  }
}

void
Domain::checkLMPCs(Connectivity *nodeToSub)
{
  if(numLMPC > 0) {
    // TODO consider the case where there is an mpc involving a node/dof that is not connected to any other elements. For example this may
    // be used to connect node with a force or a lumped mass to the structure
    if(nodeToSub) {
      for(int i=0; i < numLMPC; ++i) {
        for(int j=0; j < lmpc[i]->nterms; ++j) {
          int node = lmpc[i]->terms[j].nnum;
          if(node > -1 && nodeToSub->num(node) <= 0) // salinas mpcs can have node = -1 (indicates that node is not in any subdomains on this cpu)
            fprintf(stderr," *** WARNING: MPC %d involves bad node %d \n", lmpc[i]->lmpcnum, node+1);
        }
      }
    }
  }
}

void Domain::computeMatchingWetInterfaceLMPC() {

 if(numWet==0) return; // HB
// Create coupling LMPc's
 double tWI = 0.0;

 tWI -= getTime();

 //fprintf(stderr," ... In computeMatchingWetInterfaceLMPC(Domain.C), numWet is %d ...\n", numWet);
 //fprintf(stderr," ... wet[0] is %d ...\n", wet[0]);

 Connectivity wetElemToNode(SetAccess<ResizeArray<SommerElement *> >{wet, numWet});
 Connectivity *nodeToWetElem = wetElemToNode.alloc_reverse();

 Connectivity nodeToNode = nodeToWetElem->transcon(wetElemToNode);

 //fprintf(stderr," ... Printing WetLMPC nodeToNode ...\n");
 //nodeToNode->print();

 int maxElNodes = 0;

 int iele;
 for(iele=0;iele<numWet;iele++) {
   int nn = wet[iele]->numNodes();
   if (nn>maxElNodes) maxElNodes = nn;
 }

// double *marray = (double *) alloca(sizeof(double)*maxElNodes*maxElNodes*3);

 BaseSub *subCast = dynamic_cast<BaseSub*>(this);

 int iNode = 0;
 for (iNode=0; iNode<nodeToWetElem->csize(); iNode++) {
   //fprintf(stderr," ... iNode = %d, nodeToWetElem->num(iNode) is %d ...\n", iNode, nodeToWetElem->num(iNode));
   if (nodeToWetElem->num(iNode) > 0) {
     LMPCons *wetFSI = new LMPCons(iNode,0.0);
     int iEle;
     for(iEle = 0; iEle < nodeToWetElem->num(iNode); iEle++) {
       int jEle = (*nodeToWetElem)[iNode][iEle];
       if(subCast) {
         wet[jEle]->wetInterfaceLMPC(nodes,wetFSI,iNode);
       }
       else {
         // wet[jEle]->wetInterfaceMatrix(geoSource->GetNodes(),marray);
         wet[jEle]->wetInterfaceLMPC(geoSource->GetNodes(),wetFSI,iNode);
       }
     }
     fsi[numFSI++] = wetFSI;
     //wetFSI->print();
   }
 }

// fprintf(stderr," ... In computeMatchingWetInterfaceLMPC(Domain.C), numWet is %d, numFSI is %d ...\n", numWet,numFSI);
 //fprintf(stderr," ... numFSI is %d ...\n", numFSI);

 tWI += getTime();
 //fprintf(stderr,"Time to compute wet interface LMPCs: %f\n",tWI/1000.00);
 if(nodeToWetElem) delete nodeToWetElem; // HB: added to avoid memory leaks
}

int Domain::glToPackElem(int e)
{
 return geoSource->glToPackElem(e);
}

void
Domain::ProcessSurfaceBCs(int topFlag)
{
  BCond *surface_dbc;
  int numSurfaceDirichletBC = geoSource->getSurfaceDirichletBC(surface_dbc);
  for(int i=0; i<numSurfaceDirichletBC; ++i) {
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_dbc[i].nnum) {
        int *glNodes = SurfEntities[j]->GetPtrGlNodeIds();
        int nNodes = SurfEntities[j]->GetnNodes();
        switch(surface_dbc[i].type) {
          default:
          case BCond::Displacements : {
            BCond *bc = new BCond[nNodes];
            for(int k=0; k<nNodes; ++k) { 
              bc[k].nnum = glNodes[k];
              bc[k].dofnum = surface_dbc[i].dofnum;
              bc[k].val = surface_dbc[i].val; 
              bc[k].type = surface_dbc[i].type;
            }
            int numDirichlet_copy = geoSource->getNumDirichlet();
            geoSource->setDirichlet(nNodes, bc);
            if(numDirichlet_copy != 0) delete [] bc;
          } break;
          case BCond::Usdd : {
            BCond *bc = new BCond[nNodes];
            for(int k=0; k<nNodes; ++k) { 
              bc[k].nnum = glNodes[k];
              bc[k].dofnum = surface_dbc[i].dofnum;
              bc[k].val = surface_dbc[i].val; 
              bc[k].type = surface_dbc[i].type;
            }
            int numDirichlet_copy = geoSource->getNumDirichlet();
            geoSource->setUsddLocation(nNodes, bc);
            geoSource->setDirichlet(nNodes, bc);
            if(numDirichlet_copy != 0) delete [] bc;
          } break;
          case BCond::Lmpc : {
            for(int k=0; k<nNodes; ++k) {
              LMPCTerm term(glNodes[k], surface_dbc[i].dofnum, 1);
              LMPCons *lmpc = new LMPCons(-1,surface_dbc[i].val,&term);
              lmpc->setSource(mpc::Lmpc);
              addLMPC(lmpc,false);
            }
          } break;
        }
      }
    }
  }

  BCond *surface_nbc;
  int numSurfaceNeumanBC = geoSource->getSurfaceNeumanBC(surface_nbc);
  for(int i=0; i<numSurfaceNeumanBC; ++i) {
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_nbc[i].nnum) {
        int *glNodes = SurfEntities[j]->GetPtrGlNodeIds();
        int nNodes = SurfEntities[j]->GetnNodes();
        BCond *bc = new BCond[nNodes];
        for(int k=0; k<nNodes; ++k) { 
          bc[k].nnum = glNodes[k];
          bc[k].dofnum = surface_nbc[i].dofnum;
          bc[k].val = surface_nbc[i].val;
          bc[k].type = surface_nbc[i].type;
          bc[k].loadsetid = surface_nbc[i].loadsetid;
          bc[k].mtype = surface_nbc[i].mtype;
        }
        int numNeuman_copy = geoSource->getNumNeuman();
        geoSource->setNeuman(nNodes, bc);
        if(numNeuman_copy != 0) delete [] bc;
      }
    }
  }

  if(topFlag >= 0) return;

  PressureBCond *surface_pres;
  int numSurfacePressure = geoSource->getSurfacePressure(surface_pres);
  int nEle = geoSource->getElemSet()->last();
  for(int i=0; i<numSurfacePressure; ++i) {
    surface_pres[i].mftt = getMFTT(surface_pres[i].loadsetid);
    surface_pres[i].conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_pres[i].surfid) {
        FaceElemSet &faceElemSet = SurfEntities[j]->GetFaceElemSet();
        int last = (SurfEntities[j]->GetIsShellFace() && tdenforceFlag()) ? faceElemSet.last()/2 : faceElemSet.last();
        for(int iele = 0; iele < last; ++iele) {
          int nVertices = faceElemSet[iele]->nVertices();
          if(nVertices == 3 || nVertices == 4 || nVertices == 6 || nVertices == 8 || nVertices == 9
             || nVertices == 12 || nVertices == 10) {
             int *nodes = new int[nVertices];
             for(int inode = 0; inode < nVertices; ++inode)
               nodes[inode] = SurfEntities[j]->GetPtrGlVertexIds()[faceElemSet[iele]->GetVertex(inode)];
             int type;
             switch(nVertices) {
               case 3: type = 15; break;
               case 4: type = 16; break;
               case 6: type = 17; break;
               case 8: type = 18; break;
               case 9: type = 19; break;
               case 12: type = 20; break;
               case 10: type = 21; break;
             }
             PressureBCond *pbc = new PressureBCond(surface_pres[i]);
             pbc->elnum = -1;
             addNeumElem(-1, type, surface_pres[i].val, nVertices, nodes, pbc);
             delete [] nodes;
          }
          else {
            std::cerr << " *** ERROR: can't set pressure for surface " << SurfId << " element " << iele << " nVertices = " << nVertices << std::endl;
            continue;
          }
        }
      }
    }
  }

  BCond *surface_cfe;
  int numSurfaceConstraint = geoSource->getSurfaceConstraint(surface_cfe);
  for(int i=0; i<numSurfaceConstraint; ++i) {
    double frame_data[9];
    EFrame *ef;
    for(int j=0; j<geoSource->getNumCSframes(); ++j) {
      EFrameData &efd = geoSource->getCSframes()[j];
      if(int(surface_cfe[i].val) == efd.elnum) {
        ef = &efd.frame;
        for(int k=0; k<3; ++k) for(int l=0; l<3; ++l) frame_data[3*k+l] = efd.frame[k][l];
        for(int k=0; k<3; ++k) {
          std::cerr << "\nprint eframe==";
          for(int l=0; l<3; ++l)
            std::cerr << " " << efd.frame[k][l];
        }
        std::cerr << std::endl;  
        break;
      }
    }
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_cfe[i].nnum) {
        int *glNodes = SurfEntities[j]->GetPtrGlNodeIds();
        int nNodes = SurfEntities[j]->GetnNodes();
        switch(surface_cfe[i].type) {
          case BCond::PointPointDistance : {
            for(int k=0; k<nNodes; ++k) {
              geoSource->getElemSet()->elemadd(nEle, 77, 1, glNodes+k);
              geoSource->setAttrib(nEle, surface_cfe[i].dofnum);
              geoSource->setFrame(nEle, frame_data);
              nEle++;
            }
          } break;
          case BCond::PointLineDistance : {
            for(int k=0; k<nNodes; ++k) {
              geoSource->getElemSet()->elemadd(nEle, 78, 1, glNodes+k);
              geoSource->setAttrib(nEle, surface_cfe[i].dofnum);
              geoSource->setFrame(nEle, frame_data);
              nEle++;
            }
          } break;
          case BCond::PointPlaneDistance : {
            for(int k=0; k<nNodes; ++k) {
              geoSource->getElemSet()->elemadd(nEle, 79, 1, glNodes+k);
              geoSource->setAttrib(nEle, surface_cfe[i].dofnum);
              geoSource->setFrame(nEle, frame_data);
              nEle++;
            }
          } break;
          default :
            break;
        }
      }
    }
  }

}


void Domain::setNewProperties(int s)
{
  int na = geoSource->getNumAttributes();
  map<int, Attrib> &attr = geoSource->getAttributes();
  if(s==0) {
//   elems_copy.setMyData(true);
   int nEls = packedEset.last();
   for(int j=0; j<nEls; ++j) elems_fullcopy.elemadd(j,packedEset[j]);
   for(map<int, Group >::iterator it = geoSource->group.begin(); it != geoSource->group.end(); ++it) {   // loop over all of the groups
     for(int i = 0; i < int(it->second.attributes.size()); ++i) { // loop over attributes
        int iattr = it->second.attributes[i];

        for(int j=0; j<na; ++j) {
          if(attr[j].attr != iattr) continue;
          int jele = glToPackElem(attr[j].nele);
          if(jele == -1) continue;

          StructProp *newProp = new StructProp(*packedEset[jele]->getProperty());
          switch(it->second.randomProperties[0].rprop)  { // switch on the rprop type
            case 0:
               newProp->A = it->second.randomProperties[0].mean; // assuming one group has one random property
              break;
            case 1:
               newProp->E = it->second.randomProperties[0].mean;
              break;
            case 5:
               newProp->kx = it->second.randomProperties[0].mean;
              break;
            case 6:
               newProp->ky = it->second.randomProperties[0].mean;
              break;
            case 7:
               newProp->kz = it->second.randomProperties[0].mean;
              break;
            default:
               std::cerr << "case " << it->second.randomProperties[0].rprop << " not handled\n";
              break;
          }

          packedEset[jele]->setProp(newProp,true);
          elems_copy.elemadd(jele, packedEset[jele]);
        }
      }
    }
  }
  else {
    packedEset.deleteElems();    // clear elemset
    int count = 0;
    for(int i = 0; i < int(geoSource->group[s-1].attributes.size()); ++i) { // loop over attributes
      int iattr = geoSource->group[s-1].attributes[i];
      for(int j=0; j<na; ++j) {
          if(attr[j].attr != iattr) continue;
          int jele = glToPackElem(attr[j].nele);
          if(jele == -1) continue;
          StructProp *newProp = new StructProp(*elems_copy[jele]->getProperty());
          switch(geoSource->group[s-1].randomProperties[0].rprop)  { // switch on the rprop type
            case 0:
               if(sfem->Gauss) newProp->A = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->A = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 1:
               if(sfem->Gauss) newProp->E = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->E = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 5:
               if(sfem->Gauss) newProp->kx = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->kx = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 6:
               if(sfem->Gauss) newProp->ky = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->ky = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 7:
               if(sfem->Gauss) newProp->kz = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->kz = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            default:
               std::cerr << "case " << geoSource->group[s-1].randomProperties[0].rprop << " not handled\n";
              break;
          }
          elems_copy[jele]->setProp(newProp,true);

          Element* elem = elems_copy[jele];
          packedEset.elemadd(count++,elem);
        }
      }
      setNumElements(packedEset.last());
      deleteAllDOFs();
      makeAllDOFs();
  }
}


void Domain::assignRandMat()  // Equivalent to the non-intrusive version, but used for intrusive
{
  int ndim = sfem->getndim();
  double *xitemp = sfem->getxi();
  int na = geoSource->getNumAttributes(); // Returning the total number of elements, another option domain->numElements();
  map<int, Attrib> &attr = geoSource->getAttributes();
  int jele, count, iattr;
  for(map<int, Group >::iterator it = geoSource->group.begin(); it != geoSource->group.end(); ++it) {   // loop over all of the groups
     for(int i = 0; i < int(it->second.attributes.size()); ++i) { // loop over attributes
       iattr = it->second.attributes[i];
       for(int j=0; j<na; ++j) {
         if(attr[j].attr != iattr) continue;
            jele = glToPackElem(attr[j].nele);
          if(jele == -1)  continue;
         StructProp *newProp = new StructProp(*packedEset[jele]->getProperty());
         switch(it->second.randomProperties[0].rprop)  { // switch on the rprop type
           case 0:
              newProp->A = it->second.randomProperties[0].mean; // assuming one group has one random property
             break;
           case 1:
              newProp->E = it->second.randomProperties[0].mean;
             break;
           case 5:
              newProp->kx = it->second.randomProperties[0].mean;
             break;
           case 6:
              newProp->ky = it->second.randomProperties[0].mean;
             break;
           case 7:
              newProp->kz = it->second.randomProperties[0].mean;
             break;
           default:
              std::cerr << "case " << it->second.randomProperties[0].rprop << " not handled\n";
             break;
         }
         packedEset[jele]->setProp(newProp,true);
       }
     }
   }

   for (int s=1;s<=ndim;s++) {
    count = 0;
    for(int i = 0; i < int(geoSource->group[s-1].attributes.size()); ++i) { // loop over attributes
      iattr = geoSource->group[s-1].attributes[i];
      for(int j=0; j<na; ++j) {
          if(attr[j].attr != iattr) continue;
          jele = glToPackElem(attr[j].nele);
          if(jele == -1) continue;
          StructProp *newProp =  new StructProp(*packedEset[jele]->getProperty());
          switch(geoSource->group[s-1].randomProperties[0].rprop)  { // switch on the rprop type
            case 0:
               if(sfem->Gauss) newProp->A = newProp->A +  geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->A = newProp->A +  geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 1:
               if(sfem->Gauss) newProp->E = newProp->E + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->E = newProp->E + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 5:
               if(sfem->Gauss) newProp->kx = newProp->kx + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->kx = newProp->kx + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 6:
               if(sfem->Gauss) newProp->ky = newProp->ky + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->ky = newProp->ky + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 7:
               if(sfem->Gauss) newProp->kz = newProp->kz + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->kz = newProp->kz + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            default:
               std::cerr << "case " << geoSource->group[s-1].randomProperties[0].rprop << " not handled\n";
              break;
          }
          packedEset[jele]->setProp(newProp,true);
        }
      }
  }

}


void Domain::retrieveElemset()
{
 int na = geoSource->getNumAttributes();
 map<int, Attrib> &attr = geoSource->getAttributes();
 packedEset.deleteElems();    // clear elemset
 int count = 0;
 for(int j=0; j<na; ++j) {
   int jele = glToPackElem(attr[j].nele);
   if(jele == -1) continue;
   Element* elem = elems_fullcopy[jele];
   packedEset.elemadd(count++,elem);
 }
 setNumElements(packedEset.last());
 deleteAllDOFs();
 makeAllDOFs();
}

//CBM: new stuff
void
Domain::setEigenValue(double _lbound, int _nshifts, int _maxArnItr)
{
 // for multiple shifts eigen anaylsis
  if((_lbound < 0.0)) {
    filePrint(stderr, " *** WARNING: lbound is negative \n");
    return;
  } else {
    sinfo.lbound = _lbound;
    sinfo.nshifts = _nshifts;
  }

  sinfo.maxArnItr = _maxArnItr;
  sinfo.doEigSweep = true;
}

void
Domain::setEigenValues(double _lbound , double _ubound, int _neigps, int _maxArnItr)
{
  // for multiple shifts eigen anaylsis
  if(_lbound < 0.0 || _ubound < 0.0) {
    filePrint(stderr, " *** WARNING: lbound or ubound is negative \n");
    return;
  } else if (_lbound > _ubound) {
    filePrint(stderr, " *** WARNING: lbound > ubound, resetting interval to [ubound lbound]\n");
    sinfo.lbound = _ubound;
    sinfo.ubound = _lbound;
  } else {
    sinfo.lbound = _lbound;
    sinfo.ubound = _ubound;
  }
  sinfo.neigps = _neigps;
  sinfo.maxArnItr = _maxArnItr;
  sinfo.doEigSweep = true;
}

void
Domain::deleteAllLMPCs()
{
  lmpc.deleteArray();
  lmpc.restartArray();
  numLMPC = 0;
  nMortarLMPCs = 0;
  numCTC = 0;
  if(mortarToMPC) {
    delete mortarToMPC;
    mortarToMPC = 0;
  }
}

void
Domain::deleteSomeLMPCs(mpc::ConstraintSource s)
{ 
  int j = 0;
  for(int i = 0; i < numLMPC; ++i) {
    if(lmpc[i]->getSource() == s) {
      if(lmpc[i]->getType() == mpc::Inequality) numCTC--;
      if(s == mpc::ContactSurfaces || s == mpc::TiedSurfaces) nMortarLMPCs--;
      delete lmpc[i];
      lmpc[i] = 0;
    }
    else {
      if(i != j) {
        lmpc[j] = lmpc[i];
        lmpc[i] = 0;
      }
      j++; 
    }
  }

  numLMPC = j;
  if(mortarToMPC && (s == mpc::ContactSurfaces || s == mpc::TiedSurfaces)) {
    delete mortarToMPC;
    mortarToMPC = 0;
  }
}

void
Domain::UpdateContactSurfaceElements(GeomState *geomState)
{
  SPropContainer &sProps = geoSource->getStructProps();
  std::map<int,int> &mortar_attrib = geoSource->getMortarAttributes();

  // copy the Lagrange multipliers from geomState
  std::map<std::pair<int,int>,double> mu;
  for(int i = 0; i < numLMPC; ++i) {
    if(lmpc[i]->getSource() == mpc::ContactSurfaces) {
      if(sProps[mortar_attrib[lmpc[i]->id.first]].lagrangeMult)
        mu.insert(std::pair<std::pair<int,int>,double>(lmpc[i]->id, 0.0));
    }
  }
  geomState->getMultipliers(mu);
  geomState->clearMultiplierNodes();

  int count = 0;
  int nEle = packedEset.last();
  int count1 = 0;
  int nNode = geomState->numNodes();
  std::map<std::pair<int,int>,double>::iterator it1 = mu.begin(); 
  for(int i = 0; i < numLMPC; ++i) {
    if(lmpc[i]->getSource() == mpc::ContactSurfaces) {
      if(count < contactSurfElems.size()) { // replace
        //std::cerr << "replacing element " << contactSurfElems[count] << " with lmpc " << i << std::endl;
        packedEset.deleteElem(contactSurfElems[count]);
        packedEset.mpcelemadd(contactSurfElems[count], lmpc[i]); // replace 
        packedEset[contactSurfElems[count]]->setProp(&sProps[mortar_attrib[lmpc[i]->id.first]]);
        if(packedEset[contactSurfElems[count]]->numInternalNodes() == 1) { // i.e. Lagrange multiplier
          int in[1] = { nNode++ };
          packedEset[contactSurfElems[count]]->setInternalNodes(in);
          geomState->addMultiplierNode(it1->first, it1->second);
          it1++;
        }
        count1++;
      }
      else { // new
        //std::cerr << "adding lmpc " << i << " to elemset at index " << nEle << std::endl;
        packedEset.mpcelemadd(nEle, lmpc[i]); // new
        packedEset[nEle]->setProp(&sProps[mortar_attrib[lmpc[i]->id.first]]);
        if(packedEset[nEle]->numInternalNodes() == 1) {
          int in[1] = { nNode++ };
          packedEset[nEle]->setInternalNodes(in);
          geomState->addMultiplierNode(it1->first, it1->second);
          it1++;
        }
        contactSurfElems.push_back(nEle);
        nEle++;
      }
      count++;
    }
  }
  int count2 = 0;
  while(count < contactSurfElems.size()) {
    //std::cerr << "deleting elemset " << contactSurfElems.back() << std::endl;
    packedEset.deleteElem(contactSurfElems.back());
    contactSurfElems.pop_back();
    count2++;
  }
  packedEset.setEmax(nEle-count2); // because element set is packed
  //std::cerr << "replaced " << count1 << " and added " << count-count1 << " new elements while removing " << count2 << std::endl;
  numele = packedEset.last(); 
  numnodes = geomState->numNodes();
}

void Domain::updateSDETAF(StructProp* p, double omega) {
 if (p->etaDamp<0.0 && p->betaDamp==0.0) {
 // First time
   int tid = -int(p->etaDamp);
   int i;
   for(i = 0; i < numSDETAFT; ++i)
      if (sdetaft[i]->id == tid) break;
   if (i==numSDETAFT) { 
     fprintf(stderr,"Structural damping table %d does not exist.\n",tid);
     exit(-1);
   }
   p->etaDampTable = i;
 }
 if (p->etaDampTable>=0) {
   double f = omega/(2.0*M_PI);
   double eta, detadf;
   int tid = p->etaDampTable;
   sdetaft[tid]->getValAndSlopeAlt(f,&eta,&detadf);
   p->etaDamp = eta-f*detadf;
   p->betaDamp = detadf/(2.0*M_PI); 
 }
}

void
Domain::setIncludeStressNodes(bool *includeStressNodes)
{
  for(int iele=0; iele<numele; iele++) {
    includeStressNodes[iele] = false;
    int NodesPerElement = elemToNode->num(iele);
    for(int k=0; k<NodesPerElement; ++k) {
      int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
      int inode(0);
      bool isIn = checkIsInStressNodes(node,inode); 
      if(isIn) { 
        includeStressNodes[iele] = isIn;
        break;
      }
    }
  }
}

bool Domain::checkIsInStressNodes(int node, int &iNode)
{
  for(int inode=0; inode<stressNodes.size(); ++inode) {
    if(stressNodes[inode] == node) { iNode = inode; return true; }
  }
  
  return false;
}


void Domain::updateRUBDAFT(StructProp* p, double omega) {
#ifdef USE_EIGEN3
 if (p->eta_E<0.0 && p->rubDampTable<0) {
 // First time
   int tid = -int(p->eta_E);
   int i;
   for(i = 0; i < numRUBDAFT; ++i)
      if (rubdaft[i]->id == tid) break;
   if (i==numRUBDAFT) { 
     fprintf(stderr,"Rubber damping table %d does not exist.\n",tid);
     exit(-1);
   }
   p->rubDampTable = i;
 }
 if (p->rubDampTable>=0) {
   double f = omega/(2.0*M_PI);
   Eigen::Vector4d v, dv;
   int tid = p->rubDampTable;
   rubdaft[tid]->getValAndSlopeAlt(f, &v, &dv);
   p->E0 = v[0]-f*dv[0];
   p->eta_E = v[1]-f*dv[1];
   p->mu0 = v[2]-f*dv[2];
   p->eta_mu = v[3]-f*dv[3];
   p->dE = dv[0]/(2*M_PI);
   p->deta_E = dv[1]/(2*M_PI);
   p->dmu = dv[2]/(2*M_PI);
   p->deta_mu = dv[3]/(2*M_PI);
 }
#endif
}

void Domain::buildSensitivityInfo()
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  bool weight(false), aggregatedStressVMSensitivity(false), vonMisesStress(false), displacement(false), noAggregatedStressVM(true);
  int numOutInfo = geoSource->getNumOutInfo();
  for(int i=0; i<numOutInfo; ++i) {
    if(oinfo[i].type == OutputInfo::WeigThic) {
      senInfo[numSensitivity].type = SensitivityInfo::WeightWRTthickness;
      addSensitivity(oinfo[i]);  weight = true;
    } else if (oinfo[i].type == OutputInfo::WeigShap) {
      senInfo[numSensitivity].type = SensitivityInfo::WeightWRTshape;
      addSensitivity(oinfo[i]);  weight = true;
    } else if (oinfo[i].type == OutputInfo::AGstThic) {
      senInfo[numSensitivity].type = SensitivityInfo::AggregatedStressVMWRTthickness;
      addSensitivity(oinfo[i]);  aggregatedStressVMSensitivity = true;   aggregatedFileNumber = i;
    } else if (oinfo[i].type == OutputInfo::AGstShap) {
      senInfo[numSensitivity].type = SensitivityInfo::AggregatedStressVMWRTshape;
      addSensitivity(oinfo[i]);  aggregatedStressVMSensitivity = true;   aggregatedFileNumber = i;
    } else if (oinfo[i].type == OutputInfo::VMstThic) {
      senInfo[numSensitivity].type = SensitivityInfo::StressVMWRTthickness;
      addSensitivity(oinfo[i]);  vonMisesStress = true;
    } else if (oinfo[i].type == OutputInfo::VMstShap) {
      senInfo[numSensitivity].type = SensitivityInfo::StressVMWRTshape;
      addSensitivity(oinfo[i]);  vonMisesStress = true;
    } else if (oinfo[i].type == OutputInfo::VMstMach) {
      senInfo[numSensitivity].type = SensitivityInfo::StressVMWRTmach;
      addSensitivity(oinfo[i]);  vonMisesStress = true;
    } else if (oinfo[i].type == OutputInfo::VMstAlpha) {
      senInfo[numSensitivity].type = SensitivityInfo::StressVMWRTalpha;
      addSensitivity(oinfo[i]);  vonMisesStress = true;
    } else if (oinfo[i].type == OutputInfo::VMstBeta) {
      senInfo[numSensitivity].type = SensitivityInfo::StressVMWRTbeta;
      addSensitivity(oinfo[i]);  vonMisesStress = true;
    } else if (oinfo[i].type == OutputInfo::DispThic) {
      senInfo[numSensitivity].type = SensitivityInfo::DisplacementWRTthickness;
      addSensitivity(oinfo[i]);  displacement = true;
    } else if (oinfo[i].type == OutputInfo::DispShap) {
      senInfo[numSensitivity].type = SensitivityInfo::DisplacementWRTshape;
      addSensitivity(oinfo[i]);  displacement = true;
    } else if (oinfo[i].type == OutputInfo::DispMach) {
      senInfo[numSensitivity].type = SensitivityInfo::DisplacementWRTmach;
      addSensitivity(oinfo[i]);  displacement = true;
    } else if (oinfo[i].type == OutputInfo::DispAlph) {
      senInfo[numSensitivity].type = SensitivityInfo::DisplacementWRTalpha;
      addSensitivity(oinfo[i]);  displacement = true;
    } else if (oinfo[i].type == OutputInfo::DispBeta) {
      senInfo[numSensitivity].type = SensitivityInfo::DisplacementWRTbeta;
      addSensitivity(oinfo[i]);  displacement = true;
    } else if (oinfo[i].type == OutputInfo::AggrStVM) noAggregatedStressVM = false;
  }
  if(aggregatedStressVMSensitivity) numSensitivityQuantityTypes++;
  if(vonMisesStress) numSensitivityQuantityTypes++;
  if(displacement) numSensitivityQuantityTypes++;
  if(aggregatedStressVMSensitivity && noAggregatedStressVM) aggregatedFlag = true;

}

void Domain::addSensitivity(OutputInfo &oinfo) {
  oinfo.sentype = 1;
  senInfo[numSensitivity].surface = oinfo.surface;
  if(oinfo.type != OutputInfo::WeigThic && oinfo.type != OutputInfo::WeigShap) {
    runSAwAnalysis = true;
  }
  numSensitivity++;
}
