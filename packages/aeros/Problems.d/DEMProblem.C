#include <cstdio>
#include <alloca.h>

#include <Utils.d/MyComplex.h>
#include <Math.d/Vector.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/dofset.h>
#include <Utils.d/Connectivity.h>
#include <Solvers.d/Spooles.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/Skyline.d/SkyMatrix.h>


#include <Problems.d/DEMProblem.h>
#include <Element.d/DEM.d/DEMElement.h>
#include <Element.d/DEM.d/DEMHelm2d.h>
#include <Element.d/DEM.d/DEMHelm3d.h>
#include <Element.d/DEM.d/DEMLE2d.h>
#include <Element.d/DEM.d/DEMLE3d.h>

using namespace std;

class DEMFaceCompare: public binary_function<int,int,bool> {
 DEMFace *demF;
 public:
 DEMFaceCompare(DEMFace *_demF) { demF = _demF; }
 bool operator() (const int & i1, const int & i2) {
   DEMFace &d1 = demF[i1], &d2 = demF[i2];
   if (d1.nc<d2.nc) return true;
   else if (d1.nc>d2.nc) return false;
   else {
     for(int ci=0;ci<d1.nc;ci++) {
       if (d1.cornernodes[ci] < d2.cornernodes[ci]) return true;
       else if (d1.cornernodes[ci] > d2.cornernodes[ci]) return false;
     }
     return false;
   }
 }
};

DEMLM* DEM::createLM(int type) {
 switch (type) {
   case 1:  return  new DGMHelm2d_1_LM();
   case 2:  return  new DGMHelm2d_2_LM();
   case 3:  return  new DGMHelm2d_4_LM();
   case 4:  return  new DGMHelm2d_8_LM();
   case 5:  return  new DGMHelm2d_Eva1_2_LM();
   case 201:  return  new DGMLE2d_1_LM();
   case 202:  return  new DGMLE2d_4_LM();
   case 51:  return  new DGMHelm3d_1_LM();
   case 52:  return  new DGMHelm3d_4_LM();
   case 53:  return  new DGMHelm3d_8_LM();
   case 54:  return  new DGMHelm3d_12_LM();
   case 251:  return  new DGMLE3d_3_LM();
   case 252:  return  new DGMLE3d_15_LM();
   case 253:  return  new DGMLE3d_28_LM();
   default: fprintf(stderr,"DEM LM type %d does not exist. Exiting.\n",type); exit(-1); return 0;
 }
}

void DEM::run(Domain *d, GeoSource *g) {
 
 int gNodeOffset = d->getNodes().size();
 int nodeOffset = gNodeOffset;

// Find LM faces and artificial boundary faces and Neumann/Robin boundary faces
 fprintf(stderr,"Number of elements: %d %d %d\n",d->numElements(),d->numSommer,d->numNeum);
 int nf = d->numSommer+d->numNeum;
 for(int iele=0; iele < d->numElements(); ++iele) {
   DEMElement* deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) nf += deme->nFaces();
 }

 DEMFace *demF = new DEMFace[nf];

 nf=0;
 for(int iele=0; iele < d->numElements(); ++iele) {
   DEMElement* deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) {
     int enf = deme->nFaces();
     for(int iface=0;iface<enf;iface++) {
       demF[nf].en = iele;
       demF[nf].fi = iface;
       demF[nf].nc = deme->nFaceCorners(iface+1);
       demF[nf].cornernodes = deme->faceCorners(iface+1);
       sort(demF[nf].cornernodes,demF[nf].cornernodes+demF[nf].nc);
       demF[nf].isFluid = (deme->polyDofType()==DofSet::Helm);
       nf++;
     }
   }
 }

 int nef = nf;

 for(int iele=0;iele<d->numSommer;iele++) {
   demF[nf].en = -1;
   demF[nf].fi = -1;
   demF[nf].nc = d->sommer[iele]->nFaceCorners();
   demF[nf].cornernodes = d->sommer[iele]->faceCorners();
   sort(demF[nf].cornernodes,demF[nf].cornernodes+demF[nf].nc);
   nf++;
 } 

 for(int iele=0;iele<d->numNeum;iele++) {
   demF[nf].en = -2;
   demF[nf].fi = -2;
   demF[nf].nc = d->neum[iele]->nFaceCorners();
   demF[nf].cornernodes = d->neum[iele]->faceCorners();
   sort(demF[nf].cornernodes,demF[nf].cornernodes+demF[nf].nc);
   nf++;
 }
 fprintf(stderr,"Number of faces: %d\n",nf);

 int *demFIndex = new int[nf];
 for(int i=0;i<nf;i++) demFIndex[i] = i;
 DEMFaceCompare dfc(demF);
 sort(demFIndex,demFIndex+nf,dfc);

 FaceInfo **ele2FaceI = new FaceInfo*[d->numElements()];
 ele2FaceI[0] = new FaceInfo[nef];
 for(int iele=0;iele<d->numElements()-1;iele++) {
   DEMElement* deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) {
     ele2FaceI[iele+1] = ele2FaceI[iele] + deme->nFaces();
   }
 }
 fprintf(stderr,"Number of faces: %d\n",nf);

 int nLMFaces = 0;
 int nWetFaces = 0;
 for(int i=0;i<nf-1;i++) {
   DEMFace &d1 = demF[demFIndex[i]];
   DEMFace &d2 = demF[demFIndex[i+1]];
   if ( (!dfc(demFIndex[i],demFIndex[i+1])) && (!dfc(demFIndex[i+1],demFIndex[i])) ) {
     if (d1.en==-1 || d2.en==-1) {
       if (d1.en==-1)
          ele2FaceI[d2.en][d2.fi].bcf = 1;
       else
          ele2FaceI[d1.en][d1.fi].bcf = 1;
     } else if (d1.en==-2 || d2.en==-2) {
       if (d1.en==-2)
          ele2FaceI[d2.en][d2.fi].bcf = 2;
       else
          ele2FaceI[d1.en][d1.fi].bcf = 2;
     } else {
       ele2FaceI[d1.en][d1.fi].en2 = d2.en;
       ele2FaceI[d1.en][d1.fi].fi2 = d2.fi;
       ele2FaceI[d2.en][d2.fi].en2 = d1.en;
       ele2FaceI[d2.en][d2.fi].fi2 = d1.fi;
       if (d1.isFluid==d2.isFluid) {
         nLMFaces++;
       } else {
         ele2FaceI[d1.en][d1.fi].bcf = 3;
         ele2FaceI[d2.en][d2.fi].bcf = 3;
         nWetFaces++;
       }
     }
     i++;
   }
 }

 fprintf(stderr,"Number of LM faces: %d\n",nLMFaces);

 for(int i=0;i<nf;i++) if (demF[i].cornernodes!=0) delete[] demF[i].cornernodes;
 delete[] demF;
 delete[] demFIndex;

// Set element boundary conditions
 for(int iele=0; iele < d->numElements(); ++iele) {
   DEMElement* deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) {
      for(int fi=0;fi<deme->nFaces();fi++) deme->bc[fi] = ele2FaceI[iele][fi].bcf;
   }
 }

 double avgE = 0.0;
 double avgRhoF = 0.0;

 int nFE = 0;
 int nSE = 0;

// Initialize enrichment node offsets
 for(int iele=0; iele < d->numElements(); ++iele) {
   DEMElement *deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) {
     nodeOffset += deme->nodesE(nodeOffset);
     if (deme->polyDofType()==DofSet::Helm) {
       avgRhoF += deme->getProperty()->rho;
       nFE++;
     } else {
       avgE += deme->getProperty()->E;
       nSE++;
     }
   }
 }
 if (nFE>0) avgRhoF /= nFE;
 if (nSE>0) avgE /= nSE;
 double coupledScaling = (nSE>0 && nFE>0)?
           1.0/(geoSource->omega()*sqrt(avgE*avgRhoF)):1.0;

fprintf(stderr,"coupled scaling: %e, nFE: %d, nSE: %d\n",coupledScaling,nFE,nSE);

// Create LM multipliers from input and link with elements
 for(int ilmi=0;ilmi<nlmi;ilmi++) {
   int en = lmi[ilmi].en;
   int fi = lmi[ilmi].fi;
   DEMLM *lm;
   lm = createLM(lmi[ilmi].type);
   DEMElement* deme = dynamic_cast<DEMElement*>(d->getElementSet()[en]);
   if (deme) {
     deme->lm[fi] = lm;
     lm->e1 = deme;
   }
   DEMElement* deme2 = dynamic_cast<DEMElement*>(d->getElementSet()[ele2FaceI[en][fi].en2]);
   if (deme2) {
     deme2->lm[ele2FaceI[en][fi].fi2] = lm;
     lm->e2 = deme2;
   }
//   lm->init(d->getNodes());
   lm->init();
   lm->setNodeOffset(nodeOffset);
   nodeOffset += lm->nDofs(); 
 }

// Create implicit LM multipliers from element type and link with elements
// Also create interface DEM elements
 int numEle = d->numElements();
 for(int iele=0; iele < numEle; ++iele) {
   DEMElement* deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) {
     int type = deme->defaultLMType();
     int enf = deme->nFaces();
     DEMLM *lm;
     for(int fi=0;fi<enf;fi++) {
       if (deme->lm[fi]!=0) continue; // LM already exists
       if (ele2FaceI[iele][fi].en2==-1 && ele2FaceI[iele][fi].fi2==-1)
         continue; // no LM here
       DEMElement* deme2 = dynamic_cast<DEMElement*>
                              (d->getElementSet()[ele2FaceI[iele][fi].en2]);
       if (ele2FaceI[iele][fi].bcf==3)  {
         if (deme->polyDofType()!=DofSet::Helm) continue;
         // wet boundary
         DEMInterfaceElement *di = new DEMInterfaceElement(deme,deme2,fi+1);
         d->getElementSet().elemadd(d->numElements(),di);
         d->setNumElements(d->numElements()+1);
       } else {
         // LM here
         lm = createLM(type);
         deme->lm[fi] = lm;
         lm->e1 = deme;
         if (deme2) {
           deme2->lm[ele2FaceI[iele][fi].fi2] = lm;
           lm->e2 = deme2;
         }
//         lm->init(d->getNodes());
         lm->init();
         lm->setNodeOffset(nodeOffset);
         nodeOffset += lm->nDofs(); 
       }
     } 
   }
 }

// Create nodal connectivity underlying the dof set
 Connectivity ele2Node(d->getElementSet().asSet());
 Connectivity *node2Ele = ele2Node.alloc_reverse();
 Connectivity *node2Node = node2Ele->transcon(&ele2Node);

// Create dof set 
 DofSetArray* dsa = new DofSetArray(node2Node->csize(), d->getElementSet());
//for(int i=0;i<nodeOffset;i++) (*dsa)[i].print();

 GenSparseMatrix<complex<double> > *MK;
 GenSolver<complex<double> > *K;

// Create global matrix
 if (d->solInfo().solvercntl->type==SolverSelection::Direct) {
   switch( d->solInfo().solvercntl->subtype ) {
     default:
     case 8: 
       fprintf(stderr,"Using Spooles.\n");
       MK = new GenSpoolesSolver<complex<double> > (node2Node, dsa, *d->solInfo().solvercntl);
       K = (GenSpoolesSolver<complex<double> >*) MK;
       break;
     case 1: 
       fprintf(stderr,"Using Sparse.\n");
       MK = new GenBLKSparseMatrix<complex<double> > (node2Node, dsa, d->solInfo().solvercntl->trbm, *d->solInfo().solvercntl);
       K = (GenBLKSparseMatrix<complex<double> >*) MK;
       break;
     case 0: {
       fprintf(stderr,"Using Skyline.\n");
        int tgsm =0;
       MK = new GenSkyMatrix<complex<double> > (node2Node, dsa, d->solInfo().solvercntl->trbm, tgsm);
       K = (GenSkyMatrix<complex<double> >*) MK;
       break;
     }
   }
 } else {
   fprintf(stderr,"This solver is not supported by DEM. Exiting.\n");
   exit(-1);
 }


// Create allDOFs

 int *pointers = new int[d->numElements()+1];
 pointers[0] = 0;
 int maxNumDOFs  = 0;
 for(int iele=0; iele < d->numElements(); ++iele) {
   int numDOFs = d->getElementSet()[iele]->numDofs();
   if(numDOFs > maxNumDOFs) maxNumDOFs = numDOFs;
   pointers[iele+1] = pointers[iele] + numDOFs;
 }
 int *targets = new int[ pointers[d->numElements()] ];
 for(int iele=0; iele < d->numElements(); ++iele) {
   d->getElementSet()[iele]->dofs(*dsa, targets + pointers[iele]);
 }
 Connectivity *allDOFs = new Connectivity(d->numElements(), pointers, targets);

// Assemble
 complex<double> *karray = new complex<double> [maxNumDOFs * maxNumDOFs];

 int maxGNodes = 0;
 for(int iele=0; iele < d->numElements(); ++iele) {
   if ((iele%int(1.0+double(d->numElements())/10.0)) == 0)
      fprintf(stderr,"Assembled %.1f%%.\n",double(iele)*100.0/d->numElements());
   DEMElement *deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme!=0) {
//     deme->systemMatrix(d->getNodes(), karray);
     deme->systemMatrix(karray);
     if (deme->nGeomNodes()>maxGNodes) maxGNodes = deme->nGeomNodes(); 
//   GenStackFullM<complex<double> > kel(deme->numDofs(), karray);
     FullSquareMatrixC kel(deme->numDofs(), karray);
     if (deme->polyDofType()!=DofSet::Helm)
       kel *= (coupledScaling*coupledScaling);
     auto dofs = (*allDOFs)[iele];

     MK->add(kel,dofs);
   }
// Wet interface
   DEMInterfaceElement *di = dynamic_cast<DEMInterfaceElement*>
                                              (d->getElementSet()[iele]);
   if (di!=0) {
     fprintf(stderr,"adding wet.\n");
     complex<double> *iarray = new complex<double>[di->numDofs()*di->numDofs()];
//     di->systemMatrix(d->getNodes(), iarray);
     di->systemMatrix(iarray);
     FullSquareMatrixC kel(di->numDofs(), iarray);
     kel *= coupledScaling;
     auto dofs = (*allDOFs)[iele];
     MK->add(kel,dofs);
     delete[] iarray;
   }
 }
 delete[] karray;
 fprintf(stderr,"Done assembling.\n");

// Factorize
// K->print();

 K->factor(); 

// Create rhs and solution

 fprintf(stderr,"Problem size: %d\n",dsa->size());
 complex<double> *rhs = new complex<double>[dsa->size()];
 complex<double> *sol = new complex<double>[dsa->size()];

// Assemble rhs

// double *incdir = d->getWaveDirection();

 complex<double> *local = new complex<double>[maxNumDOFs];

 for(int iele=0; iele < d->numElements(); ++iele) {
   DEMElement *deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) {
     int ndofs = deme->numDofs();
//     deme->systemRHS(d->getNodes(), local);
     deme->systemRHS(local);
     auto dofs = (*allDOFs)[iele];

     if (deme->polyDofType()==DofSet::Helm)
       for(int idof=0;idof<ndofs;idof++)  rhs[dofs[idof]] += local[idof];
     else 
       for(int idof=0;idof<ndofs;idof++)  rhs[dofs[idof]] += local[idof]*coupledScaling;
   }
 }
 delete[] local;

// for(int i=0;i<dsa->size();i++)  fprintf(stderr,"rhs %e %e\n",real(rhs[i]),imag(rhs[i]));

// Solve
 K->solve(rhs,sol);

// for(int i=0;i<dsa->size();i++)  fprintf(stderr,"sol %e %e\n",real(sol[i]),imag(sol[i]));

// Post-process


 complex<double> (*nodalSol)[8] = new complex<double>[d->numNode()][8];
 double (*scaling)[8] = new double[d->numNode()][8];
 for(int i=0;i<d->numNode();i++) for(int j=0;j<8;j++) nodalSol[i][j] = 0.0;
 for(int i=0;i<d->numNode();i++) for(int j=0;j<8;j++) scaling[i][j] = 0.0;

 complex<double> *localSol = new complex<double>[maxNumDOFs];
 complex<double> (*localNodalSol)[8] = new complex<double>[maxGNodes][8];
 for(int iele=0; iele < d->numElements(); ++iele) {
   DEMElement *deme = dynamic_cast<DEMElement*>(d->getElementSet()[iele]);
   if (deme) {
     auto dofs = (*allDOFs)[iele];
     for(int idof=0;idof<deme->numDofs();idof++)
       localSol[idof] = sol[dofs[idof]];

     deme->nodalSolution(d->getNodes(),localSol,localNodalSol);

     int *nd = deme->nn;
     for(int i=0;i<deme->nGeomNodes();i++) {
       nodalSol[nd[i]][0] += localNodalSol[i][0];
       if (deme->polyDofType()==DofSet::Helm) scaling[nd[i]][0] += 1.0;
       for(int j=1;j<8;j++) {
          nodalSol[nd[i]][j] += coupledScaling*localNodalSol[i][j];
          if (deme->polyDofType()!=DofSet::Helm) scaling[nd[i]][j] += 1.0;
       }
     }
   }
 }
 delete[] localSol;
 delete[] localNodalSol;

 for(int i=0;i<d->numNode();i++) 
   for(int j=0;j<8;j++) if (scaling[i][j]!=0) nodalSol[i][j] /= scaling[i][j];

 complex<double>(*xyz)[3] = 
    new complex<double>[d->numNode()][3];
 for( int i = 0; i < d->numNode(); ++i) {
   for (int j = 0 ; j < 3; j++)
     xyz[i][j] = nodalSol[i][j+1];
 }

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  geoSource->openOutputFiles();
  double freq = domain->getFrequencyOrWavenumber();
  for(int i = 0; i < numOutInfo; ++i)  {
      int dof=-1;
      switch(oinfo[i].type)  {
        default:
          fprintf(stderr, " *** WARNING: output %d is not supported \n", i);
          break;
        case OutputInfo::Displacement:
          if (oinfo[i].nodeNumber == -1)
            geoSource->outputNodeVectors(i, xyz, d->numNode(), freq);
          else {
            int iNode = oinfo[i].nodeNumber;
            geoSource->outputNodeVectors(i, &(xyz[iNode]), 1, freq);
          }
          break;
        case OutputInfo::Helmholtz:
          if(dof==-1) dof = 0;
        case OutputInfo::DispX:
          if(dof==-1) dof = 1;
        case OutputInfo::DispY:
          if(dof==-1) dof = 2;
        case OutputInfo::DispZ:
          if(dof==-1) dof = 3;
          complex<double> * globVal = new complex<double>[d->numNode()];
          for(int iNode=0; iNode<d->numNode(); ++iNode) 
            globVal[iNode] = nodalSol[iNode][dof];
          if(oinfo[i].nodeNumber == -1)
            geoSource->outputNodeScalars(i, globVal, d->numNode(), freq);
          else {
            int iNode = oinfo[i].nodeNumber;
            geoSource->outputNodeScalars(i, &(globVal[iNode]), 1, freq);
          }
          delete[] globVal;
          break;
      }
  }
  delete[] xyz;
/*
 for(int i=0;i<d->numNode();i++) {
   double xyz[3];
   int nn = i;
   (d->getNodes()).getCoordinates(&nn,1,xyz,xyz+1,xyz+2);
   fprintf(stderr,"%d %e %e %e:\n!%e %e\n",i,xyz[0],xyz[1],xyz[2],real(nodalSol[i][0]),imag(nodalSol[i][0]));
   fprintf(stderr,"#%e %e %e %e %e %e\n",
                                   real(nodalSol[i][1]),imag(nodalSol[i][1]),
                                   real(nodalSol[i][2]),imag(nodalSol[i][2]),
                                   real(nodalSol[i][3]),imag(nodalSol[i][3]));
 }
*/

 delete[] rhs;
 delete[] sol;
 delete[] scaling;
 delete[] nodalSol;
 delete node2Node;
 delete K;
 delete allDOFs;
 delete dsa;
 delete[] ele2FaceI[0];
 delete[] ele2FaceI; 
// for(int iele=0; iele < d->numElements(); ++iele) {
//   delete (d->getElementSet()[iele]);
// }
 fprintf(stderr,"Done.\n");

}
