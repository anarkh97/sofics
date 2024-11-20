#include <cstdio>
#include <cmath>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/HData.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Element.d/Sommerfeld.d/LineSommerBC.h>
#include <Element.d/Sommerfeld.d/Line2SommerBC.h>
#include <Element.d/Sommerfeld.d/CurvedLine2SommerBC.h>
#include <Element.d/Sommerfeld.d/LagLineSommer.h>
#include <Element.d/Sommerfeld.d/TriangleSommerBC.h>
#include <Element.d/Sommerfeld.d/Triangle6SommerBC.h>
#include <Element.d/Sommerfeld.d/QuadSommerBC.h>
#include <Element.d/Sommerfeld.d/IsoParamQuadSommer.h>
#include <Element.d/Sommerfeld.d/SpectralIsoParamQuadSommer.h>
#include <Element.d/Sommerfeld.d/IsoParamTriSommer.h>
#include <Element.d/Sommerfeld.d/IsoParamLineSommer.h>
#include <Element.d/Sommerfeld.d/IsoParamTriLineSommer.h>
#include <Element.d/Sommerfeld.d/TrianglePressureBC.h>
#include <Element.d/Sommerfeld.d/QuadPressureBC.h>
#include <Element.d/Sommerfeld.d/Triangle6PressureBC.h>
#include <Element.d/Sommerfeld.d/Quad8PressureBC.h>
#include <Element.d/Sommerfeld.d/Quad9PressureBC.h>
#include <Element.d/Sommerfeld.d/Quad12PressureBC.h>
#include <Element.d/Sommerfeld.d/Triangle10PressureBC.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Math.d/DBSparseMatrix.h>

double HData::coupledScaling = 1.0, HData::cscale_factor = 1.0, HData::cscale_factor2 = 1.0;

extern int verboseFlag;

using std::complex;

HData::HData() : sommer(0), scatter(0), neum(0), wet(0), sBoundNodes(0)
{
  iRHS = 0;
  iWaveDir = 0;
  implicitFlag = 0;
  numSommer = 0;
  numNeum = 0;
  numScatter = 0;
  numWet = 0;
  numSBoundNodes = 0;
  numComplexDirichlet = 0;
  numComplexNeuman = 0;
  cdbc = 0;
  cnbc = 0;
  PMLFlag = 0;
  pointSourceFlag = 0;

  limitffp = 0;

  numFFPDirections = 0;
  ffpDirections = 0;

  numKirchhoffLocations = 0;
  kirchhoffLocations = 0;

  sommerfeldType = 0;
  curvatureFlag = 0;
  numWaveDirections = 0;

  scaElemToNode = 0;
  scaNodeToElem = 0;
  scaToEl = 0;

  numComplexLMPC = 0;
//  kappa = 0.0;

  frequencies = 0;
  coarse_frequencies = 0;
  isCoarseGridSolve = true;

  somElemToNode = 0;
  somNodeToElem = 0;
  somNodeToNode = 0;

  sommerChecked = false;
  subScaToSca = 0;
}

int
HData::setComplexDirichlet(int _numDirichlet, ComplexBCond *_cdbc)
{
  numComplexDirichlet = _numDirichlet;
  cdbc                = _cdbc;
  return 0;
}

int
HData::setComplexNeuman(int _numComplexNeuman, ComplexBCond * _cnbc)
{
  numComplexNeuman = _numComplexNeuman;
  cnbc             = _cnbc;
  return 0;
}

void
HData::addFSRHS(ComplexVector &force)
{
  fprintf(stderr,"HData::addFSRHS should never be called.\n");
}

void
HData::make_bc(Domain *dom, int *bc, ComplexD *bcxC)
{
 int i,iDir;
 if (implicitFlag && (numComplexDirichlet>0)) {
   ComplexBCond *newcdbc = new ComplexBCond[numComplexDirichlet * numWaveDirections];

   for (iDir=0; iDir< numWaveDirections; iDir++)
     for (i = 0; i< numComplexDirichlet; i++) newcdbc[iDir * numComplexDirichlet + i] = cdbc[i];
   cdbc = newcdbc;
 }

 int nDir = numWaveDirections;
 if (numWaveDirections==0) nDir = 1;

 for(i=0; i<dom->numdof(); ++i) {
   bc[i]  = BCFREE;
 }

 for (iDir=0; iDir< nDir; iDir++) {

   for(i=0; i<dom->numdof(); ++i) {
     bcxC[i] = ComplexD(0.0, 0.0);
   }

   // Set the real Neumann boundary conditions
   for(i=0; i<dom->numNeuman; ++i) {
     int dof  = dom->dsa->locate(dom->nbc[i].nnum, 1 << dom->nbc[i].dofnum);
     if(dof < 0) continue;
     if(bc[dof] == BCLOAD && iDir==0) {
       //fprintf(stderr,"WARNING: check input, found repeated FORCE (node %d, dof %d)\n",dom->nbc[i].nnum,dom->nbc[i].dofnum);
     }
     bc[dof] = BCLOAD;
     bcxC[dof] = ComplexD(dom->nbc[i].val, 0.0);
   }

   // Set the real Dirichlet boundary conditions
   for(i=0; i<dom->numDirichlet; ++i) {
     int dof  = dom->dsa->locate(dom->dbc[i].nnum, 1 << dom->dbc[i].dofnum);
     if(dof < 0) continue;
     if(bc[dof] == BCFIXED && iDir==0) {
       //fprintf(stderr,"WARNING: check input, found repeated DISP (node %d, dof %d)\n",dom->dbc[i].nnum,dom->dbc[i].dofnum);
     }
     bc[dof] = BCFIXED;
     bcxC[dof] = ComplexD(dom->dbc[i].val, 0.0);
   }

   // Set the Complex Neuman boundary conditions
   for(i=0; i<numComplexNeuman; ++i) {
     int dof  = dom->dsa->locate(cnbc[i].nnum, 1 << cnbc[i].dofnum);
     if(dof < 0) continue;
     if(bc[dof] == BCLOAD && iDir==0) {
       //fprintf(stderr,"WARNING: check input, found repeated HFORCE (node %d, dof %d)\n",cnbc[i].nnum,cnbc[i].dofnum);
     }
     bc[dof] = BCLOAD;
     bcxC[dof] = ComplexD(cnbc[i].reval, cnbc[i].imval);
   }

    // Set the Complex Dirichlet boundary condtions
   double kappa = 0.0;
   if (numComplexDirichlet>0 && implicitFlag )
     kappa = geoSource->kappa();
   for(i=0; i<numComplexDirichlet; ++i) {
     int dof  = dom->dsa->locate(cdbc[i].nnum, 1 << cdbc[i].dofnum);
     if(dof < 0) continue;
     if(bc[dof] == BCFIXED && iDir==0) {
       //fprintf(stderr,"WARNING: check input, found repeated HDISP (node %d, dof %d)\n",cdbc[i].nnum,cdbc[i].dofnum);
     }
     if (implicitFlag) {
       Node nd = dom->nodes.getNode(cdbc[i].nnum);
       double x = nd.x;
       double y = nd.y;
       double z = nd.z;

       double e1 = x*waveDirections[iDir*3+0];
       double e2 = y*waveDirections[iDir*3+1];
       double e3 = z*waveDirections[iDir*3+2];

       ComplexD cc = -exp( DComplex(0.0, kappa*(e1+e2+e3)) );

       cdbc[i + numComplexDirichlet * iDir].reval =  real(cc);
       cdbc[i + numComplexDirichlet * iDir].imval =  imag(cc);

       bc[dof] = BCFIXED;
       bcxC[dof] = cc;
     }
     else {
       bc[dof] = BCFIXED;
       bcxC[dof] = ComplexD(cdbc[i].reval, cdbc[i].imval);
     }
   }
   bcxC += dom->numdof();
 }
}

void
HData::addSBoundNodes()
{
  int *sBoundFlag = new int[domain->numNode()];
  for(int i=0;i<domain->numNode();i++) sBoundFlag[i] = 0;
  for(int iScatter=0; iScatter < numScatter; ++iScatter) {
    for(int iNode = 0;iNode<scatter[iScatter]->numNodes();iNode++) {
      int ndNum = (scatter[iScatter]->getNodes())[iNode];
      sBoundFlag[ndNum] = 1;
    }
  }
  for(int i=0; i<domain->numNode(); i++) if(sBoundFlag[i]) addSBoundNode(i);
  delete [] sBoundFlag;
}

void
HData::makeKss(Domain* dom)
{
 int cLen = dom->dsa->size();
 int iDof;

 int *flag = new int[dom->dsa->size()];
 int i;
 for(i=0;i<dom->dsa->size();i++ ) flag[i] = 0;

 int *eleflag = new int[dom->numele];
 for(i=0;i<dom->numele;i++) eleflag[i] = 0;

 for(i=0;i<numSBoundNodes;i++) {
   int ndNum = sBoundNodes[i];
   int iEle;
   for(iEle = 0;iEle<dom->nodeToElem->num(ndNum);iEle++) {
     int eleNum = (*dom->nodeToElem)[ndNum][iEle];
     if(!dom->isFluidElement(eleNum)) continue; // PJSA fix for coupled ... check with Radek
     if(!eleflag[eleNum]) {
       eleflag[eleNum] = 1;
       auto dofs = (*dom->allDOFs)[eleNum];
       int iDof;
       for(iDof=0;iDof<(dom->packedEset[eleNum])->numDofs();iDof++) {
         flag[dofs[iDof]] = 1;
       }
     }
   }
 }

 sBoundMap = new int[cLen];
 int sBoundLen = 0;
 for(iDof = 0; iDof < cLen; ++iDof) {
   sBoundMap[iDof] = flag[iDof] ? sBoundLen++:-1;
 }
 Kss = std::make_unique<DBSparseMatrix>(dom->nodeToNode.get(), dom->dsa, sBoundMap);

 delete[] flag;
 delete[] eleflag;

}


void
HData::getCurvatures(Domain *dom )
{
 if (curvatureFlag == 1 ) {
   if (numSommer > 0) {
     curvatures = new double[numSommer];
     int i;
     for(i=0; i<numSommer; i++)
       curvatures[i] = curvatureConst1;
   }
 }
 else {

   int (*sommerNodeToNode)[2] = new int[dom->numnodes][2];
   int i;
   for(i=0; i<dom->numnodes; i++) {
     sommerNodeToNode[i][0] = -1;
     sommerNodeToNode[i][1] = -1;
   }

   curvatures = new double[numSommer];

   for(i=0; i<numSommer; i++) {
     auto nds = sommer[i]->getNodes();
     int numElNodes = sommer[i]->numNodes();
     int first = nds[0];
     int last = nds[numElNodes-1];
     if (sommerNodeToNode[first][0] == -1)
       sommerNodeToNode[first][0] = last;
     else
       sommerNodeToNode[first][1] = last;
     if (sommerNodeToNode[last][0] == -1)
       sommerNodeToNode[last][0] = first;
     else
       sommerNodeToNode[last][1] = first;
   }

   // Actual calculation of the curvatures
   // first node by node
   double *curvNodes =  new double[dom->numnodes];
   int *lastNodeFlag = new int[dom->numnodes];
   int iNode;
   for(iNode=0; iNode<dom->numnodes; iNode++) {

     int nn[2];
     lastNodeFlag[iNode]=0;
     double curv;
     double x0, x1, x2, y0, y1, y2;

     // Current node
     Node nd0 = dom->nodes.getNode(iNode);
     x0 = nd0.x; y0 = nd0.y;

     // Its neighbors
     nn[0] = sommerNodeToNode[iNode][0];
     nn[1] = sommerNodeToNode[iNode][1];

     // Check for connections accross subdomains
     if (nn[0] != -1) {
       Node nd1 = dom->nodes.getNode(nn[0]);
       x1 = nd1.x; y1 = nd1.y;
     }
     else
       lastNodeFlag[iNode] = 1;
     if (nn[1] != -1) {
       Node nd2 = dom->nodes.getNode(nn[1]);
       x2 = nd2.x; y2 = nd2.y;
     }
     else
       lastNodeFlag[iNode] = 1;

     if(lastNodeFlag[iNode]) {
       curv = 0.0;
     }
     else {
       double l1 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
       double l2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
       double l3 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
       double area = 0.25 *
          sqrt(fabs((l1+l2+l3)*(l1+l2-l3)*(l1+l3-l2)*(l2+l3-l1)));
       curv = 4.0*fabs(area)/(l1*l2*l3);
     }
     curvNodes[iNode] = curv;
   }

   // Now average element wise
   for (i=0; i<numSommer; i++) {
     auto nds = sommer[i]->getNodes();
     int numElNodes = sommer[i]->numNodes();
     int first = nds[0];
     int last = nds[numElNodes-1];
     if ((curvNodes[first] == 0.0) && (lastNodeFlag[first]))
       curvatures[i] = curvNodes[last];
     else if ((curvNodes[last] == 0.0) && (lastNodeFlag[last]))
       curvatures[i] = curvNodes[first];
     else curvatures[i] = (curvNodes[first]+curvNodes[last])/2.0;
   }
 }
}

/*
static int lud4(double *A, int n)
{
  int i,j,k;
  double m;

  for (k=0;k<n;k++) {
    for (i=k+1;i<n;i++) {
      m = A[i*n+k]/A[k*n+k];
      if (fabs(m)<1e-12) {
        fprintf(stderr,"Zero pivot.\n");
        return -1;
      }
      A[i*n+k] = m;
      for (j=k+1;j<n;j++)
        A[i*n+j] = A[i*n+j] - m*A[k*n+j];
    }
  }
  return 0;
}

static void lub4(double *A, int n, double *b)
{
  int i,j;
  double sum;

  // Ly=b
  for (i=0;i<n;i++) {
      sum = b[i];
      for (j=0;j<i;j++)
        sum -= A[i*n+j]*b[j];
      b[i] = sum;
  }

  // Ux=y
  for (i=n-1;i>=0;i--) {
    sum = b[i];
    for (j=i+1;j<n;j++)
      sum -= A[i*n+j]*b[j];
    b[i] = sum/A[i*n+i];
  }
}
*/

static int ludcmp(double *A, int n, int *indx)// LU decomposition
//form "Numerical recipes in C", cambridge
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv = new double[n];
  for (i=0;i<n;i++) vv[i] = 0.0;
  //*d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
	for (j=0;j<n;j++)
	  if ((temp=fabs(A[i*n+j])) > big) big=temp;
	if (big==0.0) fprintf(stderr,"Singular Matrix in ludcmp-HData.C");
	vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
	  sum=A[i*n+j];
	  for (k=0;k<i;k++) sum -= A[i*n+k]*A[k*n+j];
	  A[i*n+j] = sum;
	}
	big=0.0;
	for (i=j;i<n;i++) {
	  sum = A[i*n+j];
	  for (k=0;k<j;k++)
	    sum -= A[i*n+k]*A[k*n+j];
	  A[i*n+j]=sum;
	  if ((dum=vv[i]*fabs(sum))>=big) {
		big=dum;
		imax=i;
	  }
	}
	if (j!=imax) {
	  for (k=0;k<n;k++) {
		dum=A[imax*n+k];
		A[imax*n+k]=A[j*n+k];
		A[j*n+k]=dum;
	  }
	  //*d = -(*d);
	  vv[imax]=vv[j];
	}
	indx[j]=imax;
	if (A[j*n+j]==0.0) A[j*n+j] = 1.0e-10;
	if (j!=(n-1)) {
	  dum = 1.0/(A[j*n+j]);
	  for (i=j+1;i<n;i++) A[i*n+j] *= dum;
	}
  }
  delete [] vv; vv = NULL;
  return 0;
}

static void lubksb (double *A, int n, int *indx, double *b)//LU factorisation
//from "Numerical recipes in C", cambridge
{
  int i,ii=0,ip,j;
  double sum;
  for (i=0;i<n;i++) {
	ip=indx[i];
	sum=b[ip];
	b[ip]=b[i];
	if (ii)
      for (j=ii-1;j<=i-1;j++)
		sum -= A[i*n+j]*b[j];
	else if (sum) ii=i+1;
	b[i]=sum;
  }
  for (i=n-1;i>=0;i--) {
	sum = b[i];
	for (j=i+1;j<n;j++) sum -= A[i*n+j]*b[j];
	b[i]=sum/A[i*n+i];
  }
}

void
HData::makeSommerConnectivities()
{
  if(!somNodeToNode) {
    SommerElement ** sElems = &(sommer[0]); // sbc.yieldPointer();
    somElemToNode = std::make_unique<Connectivity>(
        Connectivity::fromElements(numSommer,
                                   [&sElems](auto idx) {
                                       return std::min(sElems[idx]->numNodes(),4);
                                   },
                                   [&sElems](auto idx) {
                                       auto nodes = sElems[idx]->getNodes();
                                       auto nnd = std::min(sElems[idx]->numNodes(),4);
                                       return gsl::span<const int>{ nodes, nodes + nnd };
                                   }
        )
    );//4 to be used with bricks - JF
    somNodeToElem = std::make_unique<Connectivity>(somElemToNode->reverse());
    somNodeToNode = std::make_unique<Connectivity>(somNodeToElem->transcon(*somElemToNode));
  }
}

// 3DBAYLISS
void
HData::getCurvatures3D(Domain *dom)
{
  makeSommerConnectivities(); // PJSA reusable
  int numNodesSommer = 0;
  int iNode;
  for(iNode = 0; iNode < somNodeToNode->csize(); ++iNode) {
    if(somNodeToNode->num(iNode) > 0) {
      numNodesSommer++;
    }
  }

  if(curvatureFlag == 1 && sommerfeldType == 3) {
    if(numSommer > 0)
      fprintf (stderr, "\n\n\n Explicit curvatures supported.\n\n\n");

    nodeToSommerNodeMap = new int[dom->numnodes];
    int iCurvNode=0;
    int iNode;
    for(iNode=0; iNode<dom->numnodes; iNode++)
      nodeToSommerNodeMap[iNode] = -1;

    for(iNode=0; iNode<somNodeToElem->csize(); iNode++) {
      if(somNodeToElem->num(iNode) > 0) {
        nodeToSommerNodeMap[iNode] = iCurvNode++;
      }
    }

    curvaturesH = new double[numNodesSommer];
    curvaturesK = new double[numNodesSommer];
    double k1 =  curvatureConst1;
    double k2 =  curvatureConst1;
    for(iNode=0; iNode<numNodesSommer; iNode++) {
      curvaturesH[iNode] = 0.5*(k1+k2);
      curvaturesK[iNode] = k1*k2;
    }
  }
  else {
    curvaturesH = new double[numNodesSommer];
    curvaturesK = new double[numNodesSommer];
    curvatures_e = new double[numNodesSommer];
    curvatures_f = new double[numNodesSommer];
    curvatures_g = new double[numNodesSommer];
    curvatures_normal= new double[numNodesSommer][3];
    // This is for experiments with exact curvatures for ellipsoid
    if(curvatureFlag == 2) {
      curvatures_tau1 = new double[numNodesSommer][3];
      curvatures_tau2 = new double[numNodesSommer][3];
    }

    nodeToSommerNodeMap = new int[dom->numnodes];
    int iCurvNode=0;
    int iNode;
    for(iNode=0; iNode<somNodeToElem->csize(); iNode++) {
      nodeToSommerNodeMap[iNode] = -1;
      if(somNodeToElem->num(iNode) > 0) {
        nodeToSommerNodeMap[iNode] = iCurvNode++;
        if(somNodeToElem->num(iNode) == 1) {
          curvaturesH[iCurvNode-1] = 0.0;
          curvaturesK[iCurvNode-1] = 0.0;
          curvatures_e[iCurvNode-1] = 0.0;
          curvatures_f[iCurvNode-1] = 0.0;
          curvatures_g[iCurvNode-1] = 0.0;
          int nn=(*somNodeToElem)[iNode][0];
          double normalBuf[3];
          sommer[nn]->getNormal(dom->nodes,normalBuf);
          curvatures_normal[iCurvNode-1][0] = normalBuf[0];
          curvatures_normal[iCurvNode-1][1] = normalBuf[1];
          curvatures_normal[iCurvNode-1][2] = normalBuf[2];
        }
        else {
          //fprintf (stderr, "\n\n Getting tangent plane for Node %d ", iNode);
          double eleNormal[3]; // normal to each of the elements connected to a node
          // normal defines the tangent plane on a node
          double normal[3]; // vectorial sum of all the normals above on a node
          normal[0] = 0.0; normal[1]=0.0; normal[2]=0.0;
          // current node
          Node nd = dom->nodes.getNode(iNode);
          double x, y, z;
          x = nd.x; y = nd.y; z = nd.z;
          int iSomEle;
          for (iSomEle=0; iSomEle<somNodeToElem->num(iNode); iSomEle++) {
            int iele = (*somNodeToElem)[iNode][iSomEle];
            sommer[iele]->getNormal(dom->nodes, eleNormal);
            bool test=false;//JF
            if (test) {
              double* c=(double*)dbg_alloca(sizeof(double)*3);//can be created before...
              sommer[iele]->center(dom->nodes,c);//get the center of the element (mean of the nodes)
              double cx=c[0];
              double cy=c[1];
              double cz=c[2];
              double dist=pow(((x-cx)*(x-cx)+(y-cy)*(y-cy)+(z-cz)*(z-cz)),0.5);
              normal[0] += eleNormal[0]/dist; normal[1] += eleNormal[1]/dist; normal[2] += eleNormal[2]/dist;
            }
            else {
              normal[0] += eleNormal[0]; normal[1] += eleNormal[1]; normal[2] += eleNormal[2];
            }
          }
          double l = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
          if (l == 0) l=1.0;
          normal[0] /= l; normal[1] /= l; normal[2] /= l;
          curvatures_normal[iCurvNode-1][0] = normal[0];
          curvatures_normal[iCurvNode-1][1] = normal[1];
          curvatures_normal[iCurvNode-1][2] = normal[2];//the normal at a node is the mean of the normals for each adjacent element
          double aj[3];
          double ajBar[3];
          double *delta = (double*)dbg_alloca(sizeof(double)*(somNodeToNode->num(iNode)-1));
          double (*AN)[3] = (double (*)[3])dbg_alloca(sizeof(double)*(somNodeToNode->num(iNode)-1)*(3));
          int iSomNode;
          double tau1[3];
          double tau2[3];
          getTau(normal, tau1, tau2);//gives normalized vectors for plan
          int nodeCounter = 0;
          for (iSomNode=0; iSomNode<somNodeToNode->num(iNode); iSomNode++) {
            int nn1 = (*somNodeToNode)[iNode][iSomNode];
            if (nn1 != iNode) {
              // neighbor node
              Node nd1 = dom->nodes.getNode(nn1);
              double x1, y1, z1;
              x1 = nd1.x; y1 = nd1.y; z1 = nd1.z;
              aj[0] = x1-x;
              aj[1] = y1-y;
              aj[2] = z1-z;
              double dot1 = normal[0]*aj[0]+normal[1]*aj[1]+normal[2]*aj[2];
              double dot2 = normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2];//peut etre calcule en dehors, logiquement=1...
              double alpha = dot1/dot2;
              ajBar[0] = aj[0] - alpha*normal[0];
              ajBar[1] = aj[1] - alpha*normal[1];
              ajBar[2] = aj[2] - alpha*normal[2];
              // Project ajBar to the coord. system defined the tau1 and tau2
              double ajBar1 = tau1[0]*ajBar[0]+tau1[1]*ajBar[1]+tau1[2]*ajBar[2];
              double ajBar2 = tau2[0]*ajBar[0]+tau2[1]*ajBar[1]+tau2[2]*ajBar[2];
              // Construct AN
              AN[nodeCounter][0] = 0.5*ajBar1*ajBar1;
              AN[nodeCounter][1] = ajBar1*ajBar2;
              AN[nodeCounter][2] = 0.5*ajBar2*ajBar2;
              // Now a vector dd is defined as aj-ajBar
              double dd[3];
              dd[0] = aj[0]-ajBar[0];
              dd[1] = aj[1]-ajBar[1];
              dd[2] = aj[2]-ajBar[2];
              delta[nodeCounter] = dd[0]*normal[0]+dd[1]*normal[1]+dd[2]*normal[2];
              nodeCounter++;
            }
          }
          // Now, with AN and delta, we can assemble the local system
          // and get the local e, f and g from the least-squares problem.
          int i,j,k,n;
          double A[3][3];
          double f[3];
          n = somNodeToNode->num(iNode) - 1;
          // A = AN'*AN
          for(i=0;i<3;i++) for(j=0;j<3;j++) A[i][j] = 0.0;
          for(i=0;i<3;i++)
            for(j=0;j<3;j++)
              for(k=0;k<n;k++)
                A[i][j] += AN[k][i]*AN[k][j];

          // f = AN'*delta
          f[0]=0.0; f[1]=0.0; f[2]=0.0;
          for(i=0;i<3;i++)
            for(j=0;j<n;j++)
              f[i] += AN[j][i]*delta[j];
          // Solve the system

/*          int fl = lud4((double*)A, 3); // LU decomposition
          if (fl==0) {
            lub4((double*)A, 3, f); // LU backsubstitution
          }
          else {
            f[0]=0.0; f[1]=0.0; f[2]=0.0;
          }
          double ee = f[0];
          double ff = f[1];
          double gg = f[2];
          */

          double B[3][3];
          for (i=0;i<3;i++)
            for (j=0;j<3;j++) {
              B[i][j]=0;
            }
          double det=A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
          if (det==0)
            fprintf(stderr,"The Matrix is singular");
          else {
            B[0][0]=(A[1][1]*A[2][2]-A[2][1]*A[1][2])/det;
            B[1][0]=-(A[1][0]*A[2][2]-A[2][0]*A[1][2])/det;
            B[2][0]=(A[1][0]*A[2][1]-A[2][0]*A[1][1])/det;
            B[0][1]=-(A[0][1]*A[2][2]-A[2][1]*A[0][2])/det;
            B[1][1]=(A[0][0]*A[2][2]-A[2][0]*A[0][2])/det;
            B[2][1]=-(A[0][0]*A[2][1]-A[2][0]*A[0][1])/det;
            B[0][2]=(A[0][1]*A[1][2]-A[1][1]*A[0][2])/det;
            B[1][2]=-(A[0][0]*A[1][2]-A[1][0]*A[0][2])/det;
            B[2][2]=(A[0][0]*A[1][1]-A[1][0]*A[0][1])/det;
          }

          double ee=0;double ff=0;double gg=0;
          for (i=0;i<3;i++) {
            ee += B[0][i]*f[i];
            ff += B[1][i]*f[i];
            gg += B[2][i]*f[i];
          }

          curvatures_e[iCurvNode-1] = -ee;
          curvatures_f[iCurvNode-1] = -ff;
          curvatures_g[iCurvNode-1] = -gg;
          double b = ee+gg;
          double c = (ee*gg-ff*ff);
          double ddelta = b*b - 4*c;
          double epsilon = 1e-6;
          if (ddelta >= -epsilon  && ddelta <= epsilon) ddelta = 0.0;
          if (ddelta < 0.0)
            fprintf(stderr, "\n\n\n WARNING: 3D 2nd Order Bayliss-Turkel ABC: \n NEGATIVE DELTA for node %d \n\n\n", iNode);
          else {
            double k1 = (-b+sqrt(ddelta))/2;
            double k2 = (-b-sqrt(ddelta))/2;
            curvaturesH[iCurvNode-1]=0.5*(k1+k2);
            curvaturesK[iCurvNode-1]=k1*k2;
          }
        }
      }

      // This is for experiments with exact curvatures for ellipsoid
      if (curvatureFlag == 2) {
        Node nd = dom->nodes.getNode(iNode);
        double x, y, z;
        x = nd.x; y = nd.y; z = nd.z;

        double b=curvatureConst1, a=curvatureConst2;
        double sinasq = x*x/b/b;

        double k1 = b/a/sqrt(a*a*sinasq+b*b-b*b*sinasq);
        double k2 = a*b/pow(a*a*sinasq+b*b-b*b*sinasq,1.5);
        curvaturesH[iCurvNode-1] = 0.5*(k1+k2);
        curvaturesK[iCurvNode-1] = k1*k2;
        curvatures_e[iCurvNode-1] = k1;
        curvatures_f[iCurvNode-1] = 0.0;
        curvatures_g[iCurvNode-1] = k2;
        if (fabs(fabs(x)-b)>1e-4) {
          curvatures_tau1[iCurvNode-1][0] = 0.0;
          curvatures_tau1[iCurvNode-1][1] = -z;
          curvatures_tau1[iCurvNode-1][2] = y;
          double l = sqrt(y*y+z*z);
          curvatures_tau1[iCurvNode-1][0] = curvatures_tau1[iCurvNode-1][0]/l;
          curvatures_tau1[iCurvNode-1][1] = curvatures_tau1[iCurvNode-1][1]/l;
          curvatures_tau1[iCurvNode-1][2] = curvatures_tau1[iCurvNode-1][2]/l;
          curvatures_tau2[iCurvNode-1][0] = sqrt(b*b-x*x);
          curvatures_tau2[iCurvNode-1][1] = -x*y/sqrt(b*b-x*x);
          curvatures_tau2[iCurvNode-1][2] = -x*z/sqrt(b*b-x*x);
          l = sqrt(b*b-x*x+x*y*x*y/(b*b-x*x)+x*z*x*z/(b*b-x*x));
          curvatures_tau2[iCurvNode-1][0] = curvatures_tau2[iCurvNode-1][0]/l;
          curvatures_tau2[iCurvNode-1][1] = curvatures_tau2[iCurvNode-1][1]/l;
          curvatures_tau2[iCurvNode-1][2] = curvatures_tau2[iCurvNode-1][2]/l;
        }
        else {
          curvatures_tau1[iCurvNode-1][0] = 0.0;
          curvatures_tau1[iCurvNode-1][1] = 1.0;
          curvatures_tau1[iCurvNode-1][2] = 0.0;
          curvatures_tau2[iCurvNode-1][0] = 0.0;
          curvatures_tau2[iCurvNode-1][1] = 0.0;
          curvatures_tau2[iCurvNode-1][2] = 1.0;
        }
      }
    }
    if(curvatureFlag == 1) {
      double k1 =  curvatureConst1;
      double k2 =  curvatureConst1;
      for(iNode=0; iNode<numNodesSommer; iNode++) {
        curvaturesH[iNode] = 0.5*(k1+k2);
        curvaturesK[iNode] = k1*k2;
      }
    }
  if (numNodesSommer!=iCurvNode) fprintf(stderr,"Something is wrong here.\n");
  }
}

// 3DBAYLISS
void
HData::getCurvatures3Daccurate( Domain *dom )
{
//NOTE: curvature_e, _f, _g are define in an irrelevant basis!
//      curvaturesH, K and normal are ok.

  somNodeToNode.reset();
  makeSommerConnectivities();

  int numNodesSommer = 0;
  int iNode;
  for(iNode = 0; iNode < somNodeToNode->csize(); ++iNode) {
    if (somNodeToNode->num(iNode) > 0) {
      numNodesSommer++;
    }
  }

  if (curvatureFlag == 1 && sommerfeldType == 3) {
    if (numSommer > 0)
      fprintf (stderr, "\n\n\n Explicit curvatures supported.\n\n\n");
    dom->nodeToSommerNodeMap = new int[dom->numnodes];
    int iCurvNode=0;
    int iNode;
    for (iNode=0; iNode<dom->numnodes; iNode++)
      dom->nodeToSommerNodeMap[iNode] = -1;
    for (iNode=0; iNode<somNodeToElem->csize(); iNode++) {
      if (somNodeToElem->num(iNode) > 0) {
        dom->nodeToSommerNodeMap[iNode] = iCurvNode++;
      }
    }
    curvaturesH = new double[numNodesSommer];
    curvaturesK = new double[numNodesSommer];
    double k1 =  curvatureConst1;
    double k2 =  curvatureConst1;
    for(iNode=0; iNode<numNodesSommer; iNode++) {
      curvaturesH[iNode] = 0.5*(k1+k2);
      curvaturesK[iNode] = k1*k2;
    }
  }
  else {
    curvaturesH = new double[numNodesSommer];
    curvaturesLapH = new double[numNodesSommer];
    curvaturesK = new double[numNodesSommer];
    curvatures_e = new double[numNodesSommer];
    curvatures_f = new double[numNodesSommer];
    curvatures_g = new double[numNodesSommer];
    curvatures_normal= new double[numNodesSommer][3];
    if (curvatureFlag == 2) {// This is for experiments with exact curvatures for ellipsoid
      fprintf(stderr,"  HData.C this option is for RT and has not been checked here");
      curvatures_tau1 = new double[numNodesSommer][3];
      curvatures_tau2 = new double[numNodesSommer][3];
    }
    dom->nodeToSommerNodeMap = new int[dom->numnodes];
    int iCurvNode=0;
    int iNode;
    for (iNode=0; iNode<somNodeToElem->csize(); iNode++) {
      dom->nodeToSommerNodeMap[iNode] = -1;
      if (somNodeToElem->num(iNode) > 0) {
        dom->nodeToSommerNodeMap[iNode] = iCurvNode++;
        if (somNodeToElem->num(iNode) > 2) {//if node belongs to at least 3 elements
          // Get an approxiamtive normal = mean of the normals of adjascent elements
          double eleNormal[3]; // normal for an element adjascent to node iNode
          double normal[3]; // mean of all the normals around a node
          normal[0] = 0.0; normal[1]=0.0; normal[2]=0.0;
          Node nd = dom->nodes.getNode(iNode);// current node
          double x, y, z;
          x = nd.x; y = nd.y; z = nd.z;
          int iSomEle;
          double eleNormalMem[3];//memory for the normal at first adjascent element
          sommer[(*somNodeToElem)[iNode][0]]->getNormal(dom->nodes, eleNormalMem);
          int counter = 0;
          for (iSomEle=0; iSomEle<somNodeToElem->num(iNode); iSomEle++) {
            int iele = (*somNodeToElem)[iNode][iSomEle];
            sommer[iele]->getNormal(dom->nodes, eleNormal);
            normal[0] += eleNormal[0]; normal[1] += eleNormal[1]; normal[2] += eleNormal[2];
            if (fabs(eleNormal[0]-eleNormalMem[0])+fabs(eleNormal[0]-eleNormalMem[0])+fabs(eleNormal[0]-eleNormalMem[0])>1.0e-20)
              counter++;//stays at 0 if all noramls are the same <=> nodes are on the same plan
          }
          if (counter==0) {
            double l = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
            if (l == 0) l=1.0;
            normal[0] /= l; normal[1] /= l; normal[2] /= l;
            curvatures_normal[iCurvNode-1][0] = normal[0];
            curvatures_normal[iCurvNode-1][1] = normal[1];
            curvatures_normal[iCurvNode-1][2] = normal[2];
            curvaturesH[iCurvNode-1] = 0.0;
            curvaturesK[iCurvNode-1] = 0.0;
            curvatures_e[iCurvNode-1] = 0.0;
            curvatures_f[iCurvNode-1] = 0.0;
            curvatures_g[iCurvNode-1] = 0.0;
	  }
          else {
            double l = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
            if (l == 0) l=1.0;
            normal[0] /= l; normal[1] /= l; normal[2] /= l;
            //use this normal to change the coordinates iNode is the reference point (0,0,0)
            // and all neighbours have coordinates in the reference system attached to the plan defined prviously
            double tau1[3]; double tau2[3];
            getTau(normal, tau1, tau2);//gives normalized vectors for plan
            //rotation matrix
            //    normal[1]   normal[2] normal[3]
            //R =  tau1[1]     tau1[2]   tau1[3]
            //     tau2[1]     tau2[2]   tau2[3]
            int numNeighbours = somNodeToNode->num(iNode);
            double (*newCoords)[3] = new double[numNeighbours-1][3];//new coordinates for the neighbours
            int nodeCounter = 0;
            int iSomNode;
            for (iSomNode=0; iSomNode<numNeighbours; iSomNode++) {//change of basis
              int nn1 = (*somNodeToNode)[iNode][iSomNode];
              if (nn1 != iNode) {
                // neighbor node
                Node nd1 = dom->nodes.getNode(nn1);
                double x1, y1, z1;
                x1 = nd1.x-x; y1 = nd1.y-y; z1 = nd1.z-z;
                newCoords[nodeCounter][0] = normal[0]*x1 + normal[1]*y1 + normal[2]*z1;
                newCoords[nodeCounter][1] = tau1[0]*x1 + tau1[1]*y1 + tau1[2]*z1;
                newCoords[nodeCounter][2] = tau2[0]*x1 + tau2[1]*y1 + tau2[2]*z1;
                nodeCounter++;
              }
            }
            int unknown = 0;//number of unknowns for interpolation
            if (numNeighbours==3)
              unknown = 2;//z(x,y) = z0*x + z1*y
            else if (numNeighbours==4)
              unknown = 3;//z(x,y) = z0*x + z1*y + z2*(x*x + y*y)
            else if (numNeighbours==5)
              unknown = 4;//z(x,y) = z0*x + z1*y + z2*x*x + z3*y*y
            else if ((numNeighbours>5)&&(numNeighbours<10))
              unknown = 5;//z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y
            else
              unknown = 9;// z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y + z5*x*x*x + z6*x*x*y + z7*x*y*y + z8*y*y*y
            //double (**AN) = (double (**))dbg_alloca(sizeof(double)*(numNeighbours-1)*(unknown)); can't be safe ?!?
            double **AN = new double * [numNeighbours-1];
            for(int i=0;i<numNeighbours-1;i++) AN[i] = new double[unknown];
            double *delta = (double*)dbg_alloca(sizeof(double)*(numNeighbours-1));
            for (int i=0;i<unknown;i++) {
              for (int j=0;j<numNeighbours-1;j++)
                AN[j][i] = 0.0;
              delta[i] = 0.0;
            }

            const int Nrefi = 3;
            double Rot[3][3][Nrefi+1];
            for (int ii = 0 ; ii < Nrefi ; ii++) {//loop to refine the normal
              //keep the rotation matrices in memory
              Rot[0][0][ii]=normal[0];
              Rot[0][1][ii]=normal[1];
              Rot[0][2][ii]=normal[2];
              Rot[1][0][ii]=tau1[0];
              Rot[1][1][ii]=tau1[1];
              Rot[1][2][ii]=tau1[2];
              Rot[2][0][ii]=tau2[0];
              Rot[2][1][ii]=tau2[1];
              Rot[2][2][ii]=tau2[2];
              //fill AN
              if (unknown==2)  // z(x,y) = z0*x + z1*y - this case should not happen
                for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                  AN[iSomNode][0] = newCoords[iSomNode][1];
                  AN[iSomNode][1] = newCoords[iSomNode][2];
                }
              else if (unknown==3)// z0*x + z1*y + z2*(x*x + y*y)
                for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                  AN[iSomNode][0] = newCoords[iSomNode][1];
                  AN[iSomNode][1] = newCoords[iSomNode][2];
                  AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1] + newCoords[iSomNode][2]*newCoords[iSomNode][2];
                }
              else if (unknown==4)// z0*x + z1*y + z2*x*x + z3*y*y
                for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                  AN[iSomNode][0] = newCoords[iSomNode][1];
                  AN[iSomNode][1] = newCoords[iSomNode][2];
                  AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
                  AN[iSomNode][3] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
                }
              else if (unknown==5)//z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y
                for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                  AN[iSomNode][0] = newCoords[iSomNode][1];
                  AN[iSomNode][1] = newCoords[iSomNode][2];
                  AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
                  AN[iSomNode][3] = newCoords[iSomNode][1]*newCoords[iSomNode][2];
                  AN[iSomNode][4] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
                }
              else// z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y + z5*x*x*x + z6*x*x*y + z7*x*y*y + z8*y*y*y
                for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                  AN[iSomNode][0] = newCoords[iSomNode][1];
                  AN[iSomNode][1] = newCoords[iSomNode][2];
                  AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
                  AN[iSomNode][3] = newCoords[iSomNode][1]*newCoords[iSomNode][2];
                  AN[iSomNode][4] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
                  AN[iSomNode][5] = newCoords[iSomNode][1]*newCoords[iSomNode][1]*newCoords[iSomNode][1];
                  AN[iSomNode][6] = newCoords[iSomNode][1]*newCoords[iSomNode][1]*newCoords[iSomNode][2];
                  AN[iSomNode][7] = newCoords[iSomNode][1]*newCoords[iSomNode][2]*newCoords[iSomNode][2];
                  AN[iSomNode][8] = newCoords[iSomNode][2]*newCoords[iSomNode][2]*newCoords[iSomNode][2];
                }
              for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++)//fill delta
                delta[iSomNode] = newCoords[iSomNode][0];
              int i,j,k;
              //double **A = (double **)new double[unknown*unknown];
              //PJSA double **A = new double * [unknown];
              //PJSA for(int i=0;i<unknown;i++) A[i] = new double[unknown];
              double *A = new double[unknown*unknown]; // PJSA
              double *f = new double[unknown];
              int *vectemp = new int[unknown];
              // A = AN'*AN
              for(i=0;i<unknown;i++) {
                for(j=0;j<unknown;j++)
                  //PJSA A[i][j] = 0.0;
                  A[i*unknown+j] = 0.0; // PJSA
                f[i]=0.0;
                vectemp[i]=0;
              }
              for(i=0;i<unknown;i++)
                for(j=0;j<unknown;j++)
                  for(k=0;k<numNeighbours-1;k++)
                    //PJSA A[i][j] += AN[k][i]*AN[k][j];
                    A[i*unknown+j] += AN[k][i]*AN[k][j]; // PJSA
              for(i=0;i<unknown;i++)// f = AN'*delta
                for(j=0;j<numNeighbours-1;j++)
                  f[i] += AN[j][i]*delta[j];
              // Solve the system
              int fl = ludcmp((double*)A, unknown,vectemp); // LU decomposition
              if (fl==0)
                lubksb((double*)A, unknown, vectemp, f); // LU backsubstitution
              else
                for (i=0;i<unknown;i++)
                  f[i]=0.0;
              //The normal vector is [1,-f[1],-f[0]], note that f[0] and f[1] should be close to 0;
              normal[0] = 1; normal[1] = -f[0]; normal[2] = -f[1];
              l = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
              if (l == 0) l=1.0;
              normal[0] /= l; normal[1] /= l; normal[2] /= l;
              getTau(normal,tau1,tau2);
              for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {//change of basis
                double x1, y1, z1;
                x1 = newCoords[iSomNode][0]; y1 = newCoords[iSomNode][1]; z1 = newCoords[iSomNode][2];
                newCoords[iSomNode][0] = normal[0]*x1 + normal[1]*y1 + normal[2]*z1;
                newCoords[iSomNode][1] = tau1[0]*x1 + tau1[1]*y1 + tau1[2]*z1;
                newCoords[iSomNode][2] = tau2[0]*x1 + tau2[1]*y1 + tau2[2]*z1;
              }
              //PJSA for(int i=0;i<unknown;i++) delete [] A[i];
              delete [] A;
              delete [] f;
              delete [] vectemp;
            }
            Rot[0][0][Nrefi]=normal[0];
            Rot[0][1][Nrefi]=normal[1];
            Rot[0][2][Nrefi]=normal[2];
            Rot[1][0][Nrefi]=tau1[0];
            Rot[1][1][Nrefi]=tau1[1];
            Rot[1][2][Nrefi]=tau1[2];
            Rot[2][0][Nrefi]=tau2[0];
            Rot[2][1][Nrefi]=tau2[1];
            Rot[2][2][Nrefi]=tau2[2];
            //find the normal back: now the normal is [1,0,0]
            normal[0] = 1; normal[1] = 0; normal[2] = 0;
            double tempN[3];
            for (int ii = Nrefi ; ii >= 0 ; ii--) {
              tempN[0]=normal[0]; tempN[1]=normal[1]; tempN[2]=normal[2];
              //that to invert the rotation matrix, using the transpose is enough
              normal[0] = Rot[0][0][ii]*tempN[0]+Rot[1][0][ii]*tempN[1]+Rot[2][0][ii]*tempN[2];
              normal[1] = Rot[0][1][ii]*tempN[0]+Rot[1][1][ii]*tempN[1]+Rot[2][1][ii]*tempN[2];
              normal[2] = Rot[0][2][ii]*tempN[0]+Rot[1][2][ii]*tempN[1]+Rot[2][2][ii]*tempN[2];
            }
            // Quadratic interpolation when the normal is ok , to be modified 	 
            //do a new interpolation without the linear terms and using a^2*(z+1/a)^2+b^2*x^2+c^2*y^2=1
            //look for b^2/a and c^2/a which are the princial curvatures in hte principal coordinates
            // f0*z^2+2*z+f1*x^2+f2*x*y+f3*y^2=0
            if (unknown==2)//this case should not occur
              unknown=1;
            else if (unknown==3)
              unknown=2;
            else if (unknown==4)
              unknown=3;
            else //if (unknown>=5)
              unknown=4;
            double **AN2 = new double * [numNeighbours-1];
            for(int i=0;i<numNeighbours-1;i++) AN2[i] = new double[unknown];
            double *delta2 = (double*)dbg_alloca(sizeof(double)*(numNeighbours-1));
            for (int i=0;i<unknown;i++) {
              for (int j=0;j<numNeighbours-1;j++)
                AN2[j][i] = 0.0;
              delta2[i] = 0.0;
            }
            if (unknown==2)// f0*z^2+2*z+f1*(x^2+y^2)=0
              for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                AN2[iSomNode][0] = newCoords[iSomNode][0]*newCoords[iSomNode][0];
                AN2[iSomNode][1] = newCoords[iSomNode][1]*newCoords[iSomNode][1]+newCoords[iSomNode][2]*newCoords[iSomNode][2];
              }
            else if (unknown==3)// f0*z^2+2*z+f1*x^2+f2*y^2=0
              for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                AN2[iSomNode][0] = newCoords[iSomNode][0]*newCoords[iSomNode][0];
                AN2[iSomNode][1] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
                AN2[iSomNode][2] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
              }
            else //if (unknown==4)// f0*z^2+2*z+f1*x^2+f2*x*y+f3*y^2=0
              for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
                AN2[iSomNode][0] = newCoords[iSomNode][0]*newCoords[iSomNode][0];
                AN2[iSomNode][1] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
                AN2[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][2];
                AN2[iSomNode][3] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
              }
            for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++)//fill delta
              delta2[iSomNode] = -2*newCoords[iSomNode][0];
            int i,j,k;
            //double **A2 = (double **)new double[unknown*unknown];
            //PJSA double **A2 = new double * [unknown];
            //PJSA for(int i=0;i<unknown;i++) A2[i] = new double[unknown];
            double *A2 = new double[unknown*unknown]; // PJSA
            double *f2 = new double[unknown];
            int *vectemp2 = new int[unknown];
            // A = AN'*AN
            for(i=0;i<unknown;i++) {
              for(j=0;j<unknown;j++)
                //PJSA A2[i][j] = 0.0;
                A2[i*unknown+j] = 0.0; // PJSA
              f2[i]=0.0;
              vectemp2[i] = 0;
            }
            for(i=0;i<unknown;i++)
              for(j=0;j<unknown;j++)
                for(k=0;k<numNeighbours-1;k++)
                  //PJSA A2[i][j] += AN2[k][i]*AN2[k][j];
                  A2[i*unknown+j] += AN2[k][i]*AN2[k][j]; // PJSA
            for(i=0;i<unknown;i++)// f = AN'*delta
              for(j=0;j<numNeighbours-1;j++)
                f2[i] += AN2[j][i]*delta2[j];
            // Solve the system
            int fl = ludcmp((double*)A2, unknown, vectemp2); // LU decomposition
            if (fl==0)
              lubksb((double*)A2, unknown, vectemp2, f2); // LU backsubstitution
            else
              for (i=0;i<unknown;i++) f2[i]=0.0;
            //postcalculations
            double k1,k2;
            if (unknown==2) {
              k1 = f2[1];
              k2 = f2[1];
              curvatures_e[iCurvNode-1] = f2[1];
              curvatures_f[iCurvNode-1] = 0.0;
              curvatures_g[iCurvNode-1] = f2[1];
            }
            else if (unknown==3) {
              k1 = f2[1];
              k2 = f2[2];
              curvatures_e[iCurvNode-1] = f2[1];
              curvatures_f[iCurvNode-1] = 0.0;
              curvatures_g[iCurvNode-1] = f2[2];
            }
            else {//if (unknown==4) 
              double ee = f2[1];
              double ff = f2[2]/2;
              double gg = f2[3];
              curvatures_e[iCurvNode-1] = ee;
              curvatures_f[iCurvNode-1] = ff;
              curvatures_g[iCurvNode-1] = gg;
              //diagonalise the bilinear operator and get the eignevalues 2*2 matrix
              double dett2 = (ee+gg)*(ee+gg)-4*(ee*gg-ff*ff);
              k1 = (ee+gg+sqrt(dett2))/2;
              k2 = (ee+gg-sqrt(dett2))/2;
            }
            curvaturesH[iCurvNode-1] = 0.5*(k1+k2);
            curvaturesK[iCurvNode-1] = k1*k2;
            curvatures_normal[iCurvNode-1][0] = normal[0];
            curvatures_normal[iCurvNode-1][1] = normal[1];
            curvatures_normal[iCurvNode-1][2] = normal[2];

            delete [] newCoords;
            //PJSA for(int i=0;i<unknown;i++) delete [] A2[i];
            delete [] A2;
            delete [] f2;
            delete [] vectemp2;
            for(int i=0;i<numNeighbours-1;i++) { delete [] AN[i]; delete [] AN2[i]; }
            delete [] AN;
            delete [] AN2;
          }
        }
    	else {
    	  //if (somNodeToElem->num(iNode) <= 2) //if node belongs to only one or 2 elements on the boundary
          // first step: set curvature to 0.0
   	  curvaturesH[iCurvNode-1] = 0.0;
    	  curvaturesK[iCurvNode-1] = 0.0;
    	  curvatures_e[iCurvNode-1] = 0.0;
    	  curvatures_f[iCurvNode-1] = 0.0;
    	  curvatures_g[iCurvNode-1] = 0.0;
    	  int nn=(*somNodeToElem)[iNode][0];
    	  double normalBuf[3];
    	  sommer[nn]->getNormal(dom->nodes,normalBuf);
    	  curvatures_normal[iCurvNode-1][0] = normalBuf[0];
    	  curvatures_normal[iCurvNode-1][1] = normalBuf[1];
    	  curvatures_normal[iCurvNode-1][2] = normalBuf[2];
    	}
/*
    	// This is for experiments with exact curvatures for ellipsoid
    	if (curvatureFlag == 2) {
    	  Node nd = dom->nodes.getNode(iNode);
    	  double x, y, z;
    	  x = nd.x; y = nd.y; z = nd.z;

    	  double b=curvatureConst1, a=curvatureConst2;
    	  double sinasq = x*x/b/b;

    	  double k1 = b/a/sqrt(a*a*sinasq+b*b-b*b*sinasq);
    	  double k2 = a*b/pow(a*a*sinasq+b*b-b*b*sinasq,1.5);
    	  curvaturesH[iCurvNode-1] = 0.5*(k1+k2);
    	  curvaturesK[iCurvNode-1] = k1*k2;
    	  curvatures_e[iCurvNode-1] = k1;
    	  curvatures_f[iCurvNode-1] = 0.0;
    	  curvatures_g[iCurvNode-1] = k2;
    	  if (fabs(fabs(x)-b)>1e-4) {
            curvatures_tau1[iCurvNode-1][0] = 0.0;
            curvatures_tau1[iCurvNode-1][1] = -z;
            curvatures_tau1[iCurvNode-1][2] = y;
            double l = sqrt(y*y+z*z);
            curvatures_tau1[iCurvNode-1][0] = curvatures_tau1[iCurvNode-1][0]/l;
            curvatures_tau1[iCurvNode-1][1] = curvatures_tau1[iCurvNode-1][1]/l;
            curvatures_tau1[iCurvNode-1][2] = curvatures_tau1[iCurvNode-1][2]/l;
            curvatures_tau2[iCurvNode-1][0] = sqrt(b*b-x*x);
            curvatures_tau2[iCurvNode-1][1] = -x*y/sqrt(b*b-x*x);
            curvatures_tau2[iCurvNode-1][2] = -x*z/sqrt(b*b-x*x);
            l = sqrt(b*b-x*x+x*y*x*y/(b*b-x*x)+x*z*x*z/(b*b-x*x));
            curvatures_tau2[iCurvNode-1][0] = curvatures_tau2[iCurvNode-1][0]/l;
            curvatures_tau2[iCurvNode-1][1] = curvatures_tau2[iCurvNode-1][1]/l;
            curvatures_tau2[iCurvNode-1][2] = curvatures_tau2[iCurvNode-1][2]/l;
          }
          else {
            curvatures_tau1[iCurvNode-1][0] = 0.0;
            curvatures_tau1[iCurvNode-1][1] = 1.0;
            curvatures_tau1[iCurvNode-1][2] = 0.0;
            curvatures_tau2[iCurvNode-1][0] = 0.0;
            curvatures_tau2[iCurvNode-1][1] = 0.0;
            curvatures_tau2[iCurvNode-1][2] = 1.0;
          }
        }*/
      }
    }
    //refine the model for the edges: make a mean of neighbours curvatures instead of 0
    for (iNode=0; iNode<somNodeToElem->csize(); iNode++) {
      if (somNodeToElem->num(iNode) > 0) {
        if (somNodeToElem->num(iNode) <= 2) {//if node belongs to only one or 2 elements on the boundary
          //fprintf(stderr,"  Refinement at the edges. HData.C\n");
          //find the neighbours
          int numNeighbours = somNodeToNode->num(iNode);
          int iSomNode;
          double HHH=0.0;
          double KKK=0.0;
          int compt=0;
          for (iSomNode=0; iSomNode<numNeighbours; iSomNode++) {
            int nn1 = (*somNodeToNode)[iNode][iSomNode];
            if (nn1 != iNode) {
              //get corresponding number for boundary condition
              int nn2 = dom->nodeToSommerNodeMap[nn1];
              double HH=curvaturesH[nn2];
              double KK=curvaturesK[nn2];
              if ((HH!=0.0)&&(KK!=0.0)) {
                HHH += HH;
                KKK += KK;
                compt++;
              }
            }
          }
          if (compt) {
            HHH /= compt;
            KKK /= compt;
            int nn3 = dom->nodeToSommerNodeMap[iNode];
            curvaturesH[nn3]=HHH;
            curvaturesK[nn3]=KKK;
          }
        }
      }
    }
    //compute the surface lagrangian (Beltrami operator) of the mean curvature
    for (iNode=0; iNode<somNodeToElem->csize(); iNode++) {
      if (somNodeToElem->num(iNode) > 0) {
        int nn3 = dom->nodeToSommerNodeMap[iNode];
        if (somNodeToElem->num(iNode) > 2) {//if node belongs to more than 2 elements on the boundary
          //write the surface and the HH with the same parameters
          Node nd = dom->nodes.getNode(iNode);// current node
          double xi, yi, zi;
          xi = nd.x; yi = nd.y; zi = nd.z;
          //read the normal at central node
          double normal[3];
          normal[0] = curvatures_normal[dom->nodeToSommerNodeMap[iNode]][0];
          normal[1] = curvatures_normal[dom->nodeToSommerNodeMap[iNode]][1];
          normal[2] = curvatures_normal[dom->nodeToSommerNodeMap[iNode]][2];
          //use this normal to change the coordinates iNode is the reference point (0,0,0)
          // and all neighbours have coordinates in the reference system attached to the plan defined prviously
          double tau1[3]; double tau2[3];
          getTau(normal, tau1, tau2);//gives normalized vectors for plan
          //rotation matrix
          //    normal[1]   normal[2] normal[3]
          //R =  tau1[1]     tau1[2]   tau1[3]
          //     tau2[1]     tau2[2]   tau2[3]
          //read the nodes
          int numNeighbours = somNodeToNode->num(iNode);
          double (*newCoords)[3] = new double[numNeighbours-1][3];//new coordinates for the neighbours
          double *delta2 = (double*)dbg_alloca(sizeof(double)*(numNeighbours-1));//store the values for the curvature-curvature_at_(0,0)
          for (int i = 0 ; i < numNeighbours-1 ; i++) delta2[i] = 0.0;
          int nodeCounter = 0;
          int iSomNode;
          for (iSomNode=0; iSomNode<numNeighbours; iSomNode++) {//change of basis
            int nn1 = (*somNodeToNode)[iNode][iSomNode];
            if (nn1 != iNode) {
              // neighbor node
              Node nd1 = dom->nodes.getNode(nn1);
              double x1, y1, z1;
              x1 = nd1.x-xi; y1 = nd1.y-yi; z1 = nd1.z-zi;
          //rotate all coordinates to fit the normal to (1,0,0);
              newCoords[nodeCounter][0] = normal[0]*x1 + normal[1]*y1 + normal[2]*z1;
              newCoords[nodeCounter][1] = tau1[0]*x1 + tau1[1]*y1 + tau1[2]*z1;
              newCoords[nodeCounter][2] = tau2[0]*x1 + tau2[1]*y1 + tau2[2]*z1;
              delta2[nodeCounter] = curvaturesH[dom->nodeToSommerNodeMap[nn1]]-curvaturesH[dom->nodeToSommerNodeMap[iNode]];
              nodeCounter++;
            }
          }
          //polynomial interpolation of z and HH
          int unknown = 0;//number of unknowns for interpolation
          if (numNeighbours==3)
            unknown = 2;//z(x,y) = z0*x + z1*y idem for HH
          else if (numNeighbours==4)
            unknown = 3;//z(x,y) = z0*x + z1*y + z2*(x*x + y*y)
          else if (numNeighbours==5)
            unknown = 4;//z(x,y) = z0*x + z1*y + z2*x*x + z3*y*y
          else if ((numNeighbours>5)&&(numNeighbours<10))
            unknown = 5;//z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y
          else
            unknown = 9;// z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y + z5*x*x*x + z6*x*x*y + z7*x*y*y + z8*y*y*y
          //double **AN = (double (**))dbg_alloca(sizeof(double)*(numNeighbours-1)*unknown);
          double **AN = new double * [numNeighbours-1];
          for(int i=0;i<numNeighbours-1;i++) AN[i] = new double[unknown];
          double *delta = (double*)dbg_alloca(sizeof(double)*(numNeighbours-1));
          for (int i=0;i<unknown;i++) {
            for (int j=0;j<numNeighbours-1;j++)
              AN[j][i] = 0.0;
            delta[i] = 0.0;
          }

          //fill AN
          if (unknown==2)  // z(x,y) = z0*x + z1*y - this case should not happen
            for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
              AN[iSomNode][0] = newCoords[iSomNode][1];
              AN[iSomNode][1] = newCoords[iSomNode][2];
            }
          else if (unknown==3)// z0*x + z1*y + z2*(x*x + y*y)
            for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
              AN[iSomNode][0] = newCoords[iSomNode][1];
              AN[iSomNode][1] = newCoords[iSomNode][2];
              AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1] + newCoords[iSomNode][2]*newCoords[iSomNode][2];
            }
          else if (unknown==4)// z0*x + z1*y + z2*x*x + z3*y*y
            for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
              AN[iSomNode][0] = newCoords[iSomNode][1];
              AN[iSomNode][1] = newCoords[iSomNode][2];
              AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
              AN[iSomNode][3] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
            }
          else if (unknown==5)//z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y
            for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
              AN[iSomNode][0] = newCoords[iSomNode][1];
              AN[iSomNode][1] = newCoords[iSomNode][2];
              AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
              AN[iSomNode][3] = newCoords[iSomNode][1]*newCoords[iSomNode][2];
              AN[iSomNode][4] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
            }
          else// z(x,y) = z0*x + z1*y + z2*x*x + z3*x*y + z4*y*y + z5*x*x*x + z6*x*x*y + z7*x*y*y + z8*y*y*y
            for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++) {
              AN[iSomNode][0] = newCoords[iSomNode][1];
              AN[iSomNode][1] = newCoords[iSomNode][2];
              AN[iSomNode][2] = newCoords[iSomNode][1]*newCoords[iSomNode][1];
              AN[iSomNode][3] = newCoords[iSomNode][1]*newCoords[iSomNode][2];
              AN[iSomNode][4] = newCoords[iSomNode][2]*newCoords[iSomNode][2];
              AN[iSomNode][5] = newCoords[iSomNode][1]*newCoords[iSomNode][1]*newCoords[iSomNode][1];
              AN[iSomNode][6] = newCoords[iSomNode][1]*newCoords[iSomNode][1]*newCoords[iSomNode][2];
              AN[iSomNode][7] = newCoords[iSomNode][1]*newCoords[iSomNode][2]*newCoords[iSomNode][2];
              AN[iSomNode][8] = newCoords[iSomNode][2]*newCoords[iSomNode][2]*newCoords[iSomNode][2];
            }
          for (iSomNode=0; iSomNode<numNeighbours-1; iSomNode++)//fill delta
            delta[iSomNode] = newCoords[iSomNode][0];
          int i,j,k;
          //PJSA double **A = new double*[unknown];
          //PJSA for(i=0; i<unknown; i++)
          //PJSA   A[i] = new double[unknown];
          double *A = new double[unknown*unknown]; // PJSA
          //PJSA double *z = new double[unknown];
          double *z = new double[9]; // PJSA
          //PJSA double *h = new double[unknown];
          double *h = new double[9]; // PJSA
          int *vectemp = new int[unknown];
          // A = AN'*AN
          for(i=0;i<unknown;i++) {
            for(j=0;j<unknown;j++)
              //PJSA A[i][j] = 0.0;
              A[i*unknown+j] = 0.0;
            z[i]=0.0;
            h[i]=0.0;
            vectemp[i]=0;
          }
          for(i=0;i<unknown;i++)
            for(j=0;j<unknown;j++)
              for(k=0;k<numNeighbours-1;k++)
                //PJSA A[i][j] += AN[k][i]*AN[k][j];
                A[i*unknown+j] += AN[k][i]*AN[k][j]; // PJSA
          for(i=0;i<unknown;i++)// f = AN'*delta
            for(j=0;j<numNeighbours-1;j++) {
              z[i] += AN[j][i]*delta[j];
              h[i] += AN[j][i]*delta2[j];
            }
          // Solve the system
          int fl = ludcmp((double*)A, unknown, vectemp); // LU decomposition
          if (fl==0) {
            lubksb((double*)A, unknown, vectemp, z); // LU backsubstitution
            lubksb((double*)A, unknown, vectemp, h);
          }
          else
            for (i=0;i<unknown;i++) {
              z[i] = 0.0;
              h[i] = 0.0;
            }
          for(int i=unknown; i<9; ++i) z[i] = h[i] = 0.0; // PJSA
          //compute the lagrangian in (0,0)
          if (unknown==3)
            curvaturesLapH[nn3]=-((4*z[2]*z[0]+2*z[3]*z[1])*((1+z[1]*z[1])*h[0]-z[0]*z[1]*h[1]) +
              (2*z[3]*z[0]+4*z[4]*z[1])*(-z[0]*z[1]*h[0]+(1+z[0]*z[0])*h[1])) /
             (2*(1+z[0]*z[0]+z[1]*z[1])*(1+z[0]*z[0]+z[1]*z[1])) +
            (2*z[3]*z[1]*h[0]+(1+z[1]*z[1])*h[2]-(z[3]*z[0]+2*z[2]*z[1])*h[1] -
              (2*z[4]*z[0]+z[3]*z[1])*h[0]+2*z[3]*z[0]*h[1]+(1+z[0]*z[0])*h[2]) /
             (1+z[0]*z[0]+z[1]*z[1]);
          else if (unknown==4)
            curvaturesLapH[nn3]=-((4*z[2]*z[0]+2*z[3]*z[1])*((1+z[1]*z[1])*h[0]-z[0]*z[1]*h[1]) +
              (2*z[3]*z[0]+4*z[4]*z[1])*(-z[0]*z[1]*h[0]+(1+z[0]*z[0])*h[1])) /
             (2*(1+z[0]*z[0]+z[1]*z[1])*(1+z[0]*z[0]+z[1]*z[1])) +
            (2*z[3]*z[1]*h[0]+(1+z[1]*z[1])*h[2]-(z[3]*z[0]+2*z[2]*z[1])*h[1] -
              (2*z[4]*z[0]+z[3]*z[1])*h[0]+2*z[3]*z[0]*h[1]+(1+z[0]*z[0])*h[3]) /
             (1+z[0]*z[0]+z[1]*z[1]);
          else {
            curvaturesLapH[nn3]=-((4*z[2]*z[0]+2*z[3]*z[1])*((1+z[1]*z[1])*h[0]-z[0]*z[1]*h[1]) +
              (2*z[3]*z[0]+4*z[4]*z[1])*(-z[0]*z[1]*h[0]+(1+z[0]*z[0])*h[1])) /
             (2*(1+z[0]*z[0]+z[1]*z[1])*(1+z[0]*z[0]+z[1]*z[1])) +
            (2*z[3]*z[1]*h[0]+(1+z[1]*z[1])*h[2]-(z[3]*z[0]+2*z[2]*z[1])*h[1]-2*z[0]*z[1]*h[3] -
              (2*z[4]*z[0]+z[3]*z[1])*h[0]+2*z[3]*z[0]*h[1]+(1+z[0]*z[0])*h[4]) /
             (1+z[0]*z[0]+z[1]*z[1]);
          }
          delete [] newCoords;
          //PJSA delete [] delta2; 
          //PJSA for(i=0; i<unknown; i++)
          //PJSA   delete [] A[i];
          delete [] A;
          delete [] z;
          delete [] h;
          delete [] vectemp;
          for(i=0;i<numNeighbours-1;i++) delete [] AN[i];
          delete [] AN;
        }
        else {
          //set the value to 0
          curvaturesLapH[nn3] = 0.0;
        }
      }
    }
    for (iNode=0; iNode<somNodeToElem->csize(); iNode++) {
      if (somNodeToElem->num(iNode) > 0) {
        if (somNodeToElem->num(iNode) <= 2) {//if node belongs to only one or 2 elements on the boundary
          // compute the mean of the non zero neighbour nodes
          int numNeighbours = somNodeToNode->num(iNode);
          int iSomNode;
          double LapH=0.0;
          int compt=0;
          for (iSomNode=0; iSomNode<numNeighbours; iSomNode++) {
            int nn1 = (*somNodeToNode)[iNode][iSomNode];
            if (nn1 != iNode) {
              //get corresponding number for boundary condition
              int nn2 = dom->nodeToSommerNodeMap[nn1];
              double Lap=curvaturesLapH[nn2];
              if (Lap!=0.0) {
                LapH += Lap;
                compt++;
              }
            }
          }
          if (compt) {
            LapH /= compt;
            int nn3 = dom->nodeToSommerNodeMap[iNode];
            curvaturesLapH[nn3]=LapH;
          }
        }
      }

      if (curvatureFlag == 1) {
        fprintf(stderr,"  HData.C: curvatureFlag should not be 1");
        /*double k1 =  curvatureConst1;
        double k2 =  curvatureConst1;
        for(iNode=0; iNode<numNodesSommer; iNode++) {
          curvaturesH[iNode] = 0.5*(k1+k2);
          curvaturesK[iNode] = k1*k2;
        }*/
      }
    }
    if (numNodesSommer!=iCurvNode) fprintf(stderr,"Something is wrong here.\n");
  }
}


void
HData::getTau(double normal[3], double tau1[3], double tau2[3])
{
  if(fabs(normal[0]) <= fabs(normal[1]) && fabs(normal[0]) <= fabs(normal[2]))
  {
    tau1[0] =  0.0;
    tau1[1] = -normal[2];
    tau1[2] =  normal[1];
  }
  if(fabs(normal[1]) <= fabs(normal[0]) && fabs(normal[1]) <= fabs(normal[2]))
  {
    tau1[0] = -normal[2];
    tau1[1] =  0.0;
    tau1[2] =  normal[0];
  }
  if(fabs(normal[2]) <= fabs(normal[0]) && fabs(normal[2]) <= fabs(normal[0]))
  {
    tau1[0] = -normal[1];
    tau1[1] =  normal[0];
    tau1[2] =  0.0;
  }

  tau2[0] = normal[1]*tau1[2]-normal[2]*tau1[1];
  tau2[1] = normal[2]*tau1[0]-normal[0]*tau1[2];
  tau2[2] = normal[0]*tau1[1]-normal[1]*tau1[0];

  double l = sqrt(tau1[0]*tau1[0]+tau1[1]*tau1[1]+tau1[2]*tau1[2]);
  tau1[0] /= l;
  tau1[1] /= l;
  tau1[2] /= l;
  l = sqrt(tau2[0]*tau2[0]+tau2[1]*tau2[1]+tau2[2]*tau2[2]);
  tau2[0] /= l;
  tau2[1] /= l;
  tau2[2] /= l;

}

void
HData::outputFFP(ComplexVector& sol, int iInfo)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();

// RT new style
 if (oinfo[iInfo].type == OutputInfo::Farfield
     && numFFPDirections==0) {

   OutputInfo *oinfo = geoSource->getOutputInfo();
   int dim= scatter[0]->dim();
   int nsint = oinfo[iInfo].interval;
   int numPhi = (dim == 3) ? nsint/2+1 : 1;
   int i, j;

   DComplex ffpCoef;
   if(dim!=3) ffpCoef = exp(DComplex(0.0,M_PI/4.0))/sqrt(8.0*M_PI*geoSource->kappa())*geoSource->global_average_rhof;
   else ffpCoef = DComplex(0.25/M_PI, 0.0)*geoSource->global_average_rhof;

   DComplex *p = new DComplex[numPhi*nsint];
   for(i=0;i<numPhi*nsint;i++) p[i] = complex<double>(0.0,0.0);

   double (*vectorDir)[3] = new double[numPhi*nsint][3];
   int numSample = 0;

   for(i=0; i<numPhi; ++i) {
     double phi = (numPhi==1) ? 0 : M_PI*(-0.5+((double) i)/(numPhi-1.0));
     for (j=0; j<nsint; ++j) {
       vectorDir[numSample][0] = cos(phi)*cos(2*j*M_PI/double(nsint));
       vectorDir[numSample][1] = cos(phi)*sin(2*j*M_PI/double(nsint));
       vectorDir[numSample][2] = sin(phi);
       numSample += 1;
     }
   }

   ffp(domain, numSample, p, vectorDir, sol.data(), true);

   // OUTPUT result in FILE
   // 2D -> theta | real part | imaginary part | logarithmic value
   // 3D -> theta | phi | real part | imaginary | logarithmic value
   if(dim != 3) {
     for(i=0;i<nsint;i++) {
       p[i] *= ffpCoef; // PJSA
       //double y = 10.0 * log(2.0*M_PI*abs(p[i])*abs(p[i]))/log(10.0);
       fprintf(oinfo[iInfo].filptr,"%e  %.10e  %.10e\n", 2*M_PI*i/nsint,
               real(p[i]), imag(p[i]));
     }
   }
   else {
     numSample = 0;
     for(i=0; i<numPhi; ++i) {
       double phi = (numPhi==1) ? 0 : M_PI*(-0.5+((double) i)/(numPhi-1.0));
       for(j=0; j<nsint; ++j) {
         p[numSample] *= ffpCoef; // PJSA
         //double y = 10.0 * log(2.0*M_PI*abs(p[numSample])*abs(p[numSample]))/log(10.0);
         fprintf(oinfo[iInfo].filptr,"%e  %e  %.10e  %.10e\n", 2*M_PI*j/nsint, phi,
                 ScalarTypes::Real(p[numSample]),ScalarTypes::Imag(p[numSample]));
         numSample += 1;
       }
     }
   }
 } else  {
// RT: new style input/output
   complex<double> *c = (complex<double>*)
    alloca((domain->numDirichlet+numComplexDirichlet)*sizeof(complex<double>));
   int iDof;
   for(iDof = 0; iDof < domain->numDirichlet+numComplexDirichlet; ++iDof)
     c[iDof] = 0.0;
   int i;
   for(i=0; i<domain->numDirichlet; ++i) {
     int dof2 = domain->dsa->locate(domain->dbc[i].nnum,(1 << domain->dbc[i].dofnum));
     dof2 = domain->c_dsa->invRCN(dof2);
     if(dof2 >= 0) c[dof2] = complex<double>(domain->dbc[i].val, 0.0);
   }
   ComplexBCond * cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
   for(i=0; i<numComplexDirichlet; ++i) {
      int dof2 = domain->dsa->locate(cdbcMRHS[i].nnum,(1 << cdbcMRHS[i].dofnum));
      dof2 = domain->c_dsa->invRCN(dof2);
      if(dof2 >= 0)
        c[dof2] = complex<double>( cdbcMRHS[i].reval, cdbcMRHS[i].imval);
   }

   ComplexD * uel = (ComplexD*) alloca(domain->maxNumDOFs*sizeof(ComplexD));
   int *nds = (int*) alloca(domain->maxNumDOFs*sizeof(int));
   int *dofs = (int*) alloca(domain->maxNumDOFs*sizeof(int));

   bool direction;
   int numEvalDirLoc;
   double *evalDirLoc;

   switch (oinfo[iInfo].type) {
   case OutputInfo::Farfield:
     direction = true;
     numEvalDirLoc = numFFPDirections;
     evalDirLoc = ffpDirections;
     break;

   case OutputInfo::Kirchhoff:
     direction = false;
     numEvalDirLoc = numKirchhoffLocations;
     evalDirLoc = kirchhoffLocations;
     break;

   default:
     fprintf(stderr,"Error in HData:outputFFP\n");
     exit(-1);
   }

   complex<double> *ffp = new complex<double>[numEvalDirLoc];

   int iele;
   for(iele=0; iele < numScatter; ++iele) {
     int numEleDofs = scatter[iele]->el->numDofs();
     scatter[iele]->el->nodes(nds);
     scatter[iele]->el->dofs(*domain->dsa,dofs);
     int i;
     for(i=0;i<numEleDofs;i++) {
       int hLoc  = domain->c_dsa->locate(nds[i], DofSet::Helm);
       int hLoc1 =   domain->dsa->locate(nds[i], DofSet::Helm);
       if (hLoc >= 0) uel[i] = sol[hLoc];
       else if (hLoc1 >=0 ) uel[i] = c[domain->c_dsa->invRCN(hLoc1)];
       else fprintf(stderr,"Error in HData:outputFFP\n");
     }
     scatter[iele]->ffp(domain->nodes,numEvalDirLoc,evalDirLoc,uel,ffp,direction);
   }
   for(i=0;i<numEvalDirLoc;i++) {
     fprintf(oinfo[iInfo].filptr,"%e %e %e   %e %e\n",
       evalDirLoc[i*3+0],evalDirLoc[i*3+1],evalDirLoc[i*3+2],
       real(ffp[i]),imag(ffp[i]));
   }
   delete[] ffp;
 }
}

void
HData::addSommer(SommerElement *ele)
{
 ele->dom = domain;
 sommer[numSommer++] = ele;
}

void
HData::addWet(SommerElement *ele)
{
 ele->dom = domain;
 wet[numWet++] = ele;
}

void
HData::addScatter(SommerElement *ele)
{
 ele->dom = domain;
 scatter[numScatter++] = ele;
}

void
HData::addNeum(SommerElement *ele)
{
 ele->dom = domain;
 neum[numNeum++] = ele;
}

void HData::setFFP(int _nffp)
{
 nffp = _nffp;
}

void
HData::addSBoundNode(int nd)
{
 sBoundNodes[numSBoundNodes++] = nd;
}

void
HData::setFFP(int _nffp, int _limitffp)
{
 nffp = _nffp;
 limitffp = _limitffp;
}

void
HData::setKirchhoffLocations(double x, double y, double z)
{
  if (numKirchhoffLocations%100 == 0) {
     int i;
     double *p = new double[(numKirchhoffLocations+100)*3];
     for(i=0;i<numKirchhoffLocations*3;i++) p[i] = kirchhoffLocations[i];
     if (numKirchhoffLocations!=0) delete [] kirchhoffLocations;
     kirchhoffLocations = p;
  }
  kirchhoffLocations[numKirchhoffLocations*3  ] = x;
  kirchhoffLocations[numKirchhoffLocations*3+1] = y;
  kirchhoffLocations[numKirchhoffLocations*3+2] = z;
  numKirchhoffLocations++;
}

void
HData::setFFPDirections(double d1, double d2, double d3)
{
  if (numFFPDirections%100 == 0) {
     int i;
     double *p = new double[(numFFPDirections+100)*3];
     for(i=0;i<numFFPDirections*3;i++) p[i] = ffpDirections[i];
     if (numFFPDirections!=0) delete [] ffpDirections;
     ffpDirections = p;
  }
  double d = sqrt(d1*d1 + d2*d2 + d3*d3);
  if (d==0.0) {
    fprintf(stderr,"HData::setFFPDirections: Direction is zero.\n");
    exit(-1);
  }
  ffpDirections[numFFPDirections*3] = d1 / d;
  ffpDirections[numFFPDirections*3+1] = d2 / d ;
  ffpDirections[numFFPDirections*3+2] = d3 / d ;
  numFFPDirections++;
}

void
HData::setWaveDirections(int numDir, double d1, double d2)
{
  if (numDir < 2) numDir = 2;
  int i;
  waveDirections = new double[ numDir * 3];
  for (i=0; i < numDir; i++) {
    double theta = d1 + i * (d2 - d1) / (numDir - 1);
    waveDirections[i*3] = cos(theta);
    waveDirections[i*3 + 1] = sin(theta);
    waveDirections[i*3 + 2] = 0.0;
  }
  numWaveDirections = numDir;
}

void
HData::setWaveDirections(int iDir, double d1, double d2, double d3)
{
  double d = sqrt(d1*d1 + d2*d2 + d3*d3);
  if (pointSourceFlag) d = 1.0;
  if (iDir==0) {
    if (numWaveDirections %100 == 0) {
       int i;
       double *p = new double[(numWaveDirections+100)*3];
       for(i=0;i<numWaveDirections*3;i++) p[i] = waveDirections[i];
       if (numWaveDirections!=0) delete [] waveDirections;
       waveDirections = p;
    }
    if (d==0.0) {
      fprintf(stderr,"HData::setWaveDirection: Direction is zero.\n");
      exit(-1);
    }
    waveDirections[numWaveDirections*3] = d1 / d;
    waveDirections[numWaveDirections*3+1] = d2 / d ;
    waveDirections[numWaveDirections*3+2] = d3 / d ;
    numWaveDirections++;
  }
  else if (d1==0.0 && d2==0.0 && d3==0.0) {
    waveDirections = new double[ iDir * 3];
    numWaveDirections = iDir;
  } else if (iDir>0) {
    if (iDir>numWaveDirections || d==0.0) {
      fprintf(stderr,"HData::setWaveDirection: Wrong direction number.\n");
      exit(-1);
    }
    waveDirections[(iDir-1)*3] = d1 / d;
    waveDirections[(iDir-1)*3+1] = d2 / d ;
    waveDirections[(iDir-1)*3+2] = d3 / d ;
  }
}

void
HData::addSommerElem(int num, int etype, double sommerConst, int nnodes, int*n)
{
   SommerElement *ele;
   switch(etype)
   {
     case 1:
       ele = new LineSommerBC(n[0], n[1]);
       addSommer(ele);
       break;
     case 2:
       ele = new Line2SommerBC(n[0], n[1], n[2]);
       addSommer(ele);
       break;
     case 3:
       ele = new TriangleSommerBC(n[0], n[1], n[2]);
       addSommer(ele);
       break;
     case 4: {
       //(use for debugging) int ntmp[4] = { n[0], n[1], n[3], n[2] };
       //(use for debugging) ele = new IsoParamQuadSommer(nnodes,ntmp);
       ele = new QuadSommerBC(n[0], n[1], n[2], n[3]);
       addSommer(ele);
       } break;
     case 6:
       ele = new Triangle6SommerBC(n[0], n[1], n[2], n[3], n[4], n[5]);
       addSommer(ele);
       break;
     case 7:
       printf("LETriangleSommerBC is not supported \n");
       // ele = new LETriangleSommerBC(n[0], n[1], n[2]);
       // addSommer(ele);
       break;
     case 8:
       ele = new CurvedLine2SommerBC(n[0], n[1], n[2]);
       addSommer(ele);
       break;
     case 9:
       ele = new LagLineSommer(nnodes,n);
       addSommer(ele);
       break;
#ifndef SALINAS
     case 10:
       ele = new IsoParamQuadSommer(nnodes,n);
       addSommer(ele);
       break;
     case 11:
       ele = new IsoParamTriSommer(nnodes,n);
       addSommer(ele);
       break;
     case 12:
       ele = new IsoParamLineSommer(nnodes,n);
       addSommer(ele);
       break;
     case 13:
       ele = new IsoParamTriLineSommer(nnodes,n);
       addSommer(ele);
       break;
     case 14:
       ele = new SpectralIsoParamQuadSommer(nnodes,n);
       addSommer(ele);
       break;
#endif
     default:
       return;
   }
}


void
HData::addWetElem(int num, int etype, double sommerConst, int nnodes, int*n)
{
   SommerElement *ele;
   switch(etype)
   {
     case 1:
       ele = new LineSommerBC(n[0], n[1]);
       addWet(ele);
       break;
     case 2:
       ele = new Line2SommerBC(n[0], n[1], n[2]);
       addWet(ele);
       break;
     case 3:
       ele = new TriangleSommerBC(n[0], n[1], n[2]);
       addWet(ele);
       break;
     case 4:
       ele = new QuadSommerBC(n[0], n[1], n[2], n[3]);
       addWet(ele);
       break;
     case 6:
       ele = new Triangle6SommerBC(n[0], n[1], n[2], n[3], n[4], n[5]);
       addWet(ele);
       break;
     case 7:
       printf("LETriangleSommerBC is not supported \n");
       // ele = new LETriangleSommerBC(n[0], n[1], n[2]);
       // addWet(ele);
       break;
     case 8:
       ele = new CurvedLine2SommerBC(n[0], n[1], n[2]);
       addWet(ele);
       break;
     case 9:
       ele = new LagLineSommer(nnodes,n);
       addWet(ele);
       break;
#ifndef SALINAS
     case 10:
       ele = new IsoParamQuadSommer(nnodes,n);
       addWet(ele);
       break;
     case 11:
       ele = new IsoParamTriSommer(nnodes,n);
       addWet(ele);
       break;
     case 12:
       ele = new IsoParamLineSommer(nnodes,n);
       addWet(ele);
       break;
     case 13:
       ele = new IsoParamTriLineSommer(nnodes,n);
       addWet(ele);
       break;
#endif
     default:
       return;
   }
}

void
HData::addScatterElem(int num, int etype, double sommerConst, int nnodes, int*n)
{
   SommerElement *ele;

   switch(etype)
   {
     case 1:
       ele = new LineSommerBC(n[0], n[1]);
       addScatter(ele);
       break;
     case 2:
       ele = new Line2SommerBC(n[0], n[1], n[2]);
       addScatter(ele);
       break;
     case 3:
       ele = new TriangleSommerBC(n[0], n[1], n[2]);
       addScatter(ele);
       break;
     case 4:
       ele = new QuadSommerBC(n[0], n[1], n[2], n[3]);
       addScatter(ele);
       break;
     case 6:
       ele = new Triangle6SommerBC(n[0], n[1], n[2], n[3], n[4], n[5]);
       addScatter(ele);
       break;
     case 7:
       printf("LETriangleSommerBC is not supported \n");
       // ele = new LETriangleSommerBC(n[0], n[1], n[2]);
       // addScatter(ele);
       break;
     case 8:
       ele = new CurvedLine2SommerBC(n[0], n[1], n[2]);
       addScatter(ele);
       break;
     case 9:
       ele = new LagLineSommer(nnodes,n);
       addScatter(ele);
       break;
#ifndef SALINAS
     case 10:
       ele = new IsoParamQuadSommer(nnodes,n);
       addScatter(ele);
       break;
     case 11:
       ele = new IsoParamTriSommer(nnodes,n);
       addScatter(ele);
       break;
     case 12:
       ele = new IsoParamLineSommer(nnodes,n);
       addScatter(ele);
       break;
     case 13:
       ele = new IsoParamTriLineSommer(nnodes,n);
       addScatter(ele);
       break;
#endif
     default:
       return;
   }
}

void
HData::addNeumElem(int num, int etype, double sommerConst, int nnodes, int *n, PressureBCond *pbc)
{
   SommerElement *ele;

   switch(etype)
   {
     case 1:
       ele = new LineSommerBC(n[0], n[1]);
       addNeum(ele);
       break;
     case 2:
       ele = new Line2SommerBC(n[0], n[1], n[2]);
       addNeum(ele);
       break;
     case 3:
       ele = new TriangleSommerBC(n[0], n[1], n[2]);
       addNeum(ele);
       break;
     case 4:
       //printf(" ... Add a QuadSommerBC el.\n");
       ele = new QuadSommerBC(n[0], n[1], n[2], n[3]);
       addNeum(ele);
       break;
     case 6:
       ele = new Triangle6SommerBC(n[0], n[1], n[2], n[3], n[4], n[5]);
       addNeum(ele);
       break;
     case 7:
       printf("LETriangleSommerBC is not supported \n");
       // ele = new LETriangleSommerBC(n[0], n[1], n[2]);
       // addNeum(ele);
       break;
     case 8:
       ele = new CurvedLine2SommerBC(n[0], n[1], n[2]);
       addNeum(ele);
       break;
     case 9:
       ele = new LagLineSommer(nnodes,n);
       addNeum(ele);
       break;
#ifndef SALINAS
     case 10:
       ele = new IsoParamQuadSommer(nnodes,n);
       addNeum(ele);
       break;
     case 11:
       ele = new IsoParamTriSommer(nnodes,n);
       addNeum(ele);
       break;
     case 12:
       ele = new IsoParamLineSommer(nnodes,n);
       addNeum(ele);
       break;
     case 13:
       ele = new IsoParamTriLineSommer(nnodes,n);
       addNeum(ele);
       break;
#endif
     case 15:
       ele = new TrianglePressureBC(n, pbc);
       addNeum(ele);
       break;
     case 16:
       ele = new QuadPressureBC(n, pbc);
       addNeum(ele);
       break;
     case 17:
       ele = new Triangle6PressureBC(n, pbc);
       addNeum(ele);
       break;
     case 18:
       ele = new Quad8PressureBC(n, pbc);
       addNeum(ele);
       break;
     case 19:
       ele = new Quad9PressureBC(n, pbc);
       addNeum(ele);
       break;
     case 20:
       ele = new Quad12PressureBC(n, pbc);
       addNeum(ele);
       break;
     case 21:
       ele = new Triangle10PressureBC(n, pbc);
       addNeum(ele);
       break;
     default:
       return;
   }
}

void
HData::inithScatter(int m, int n)
{
  ihScatter = 0;
  lenhScatter = m;
  numhScatter = n;
  hScatter = new double[m*n][3];
}

void
HData::addhScatter(int n, double *d)
{
  int i;
  if (n/3 != numhScatter) fprintf(stderr,"Error in HData::addhScatter.\n");
  for(i=0;i<n/3;i++) {
    hScatter[lenhScatter*i+ihScatter][0] = d[3*i];
    hScatter[lenhScatter*i+ihScatter][1] = d[3*i+1];
    hScatter[lenhScatter*i+ihScatter][2] = d[3*i+2];
  }
  ihScatter++;
}

void
HData::ffp(Domain *dom, int ndir, DComplex *ffp, double (*dir)[3], ComplexD *u, bool direction)
{
  if(dom->numFFPDirections > 0 || !direction) {
    complex<double> *c;
    if(dom->numDirichlet+numComplexDirichlet > 0)
      c = (complex<double> *) alloca((dom->numDirichlet+numComplexDirichlet)*sizeof(complex<double>));
    for(int iDof = 0; iDof < dom->numDirichlet+numComplexDirichlet; ++iDof) c[iDof] = 0.0;
    int i;
    for(i = 0; i < dom->numDirichlet; ++i) {
      int dof2 = dom->dsa->locate(dom->dbc[i].nnum,(1 << dom->dbc[i].dofnum));
      dof2 = dom->c_dsa->invRCN(dof2);
      if(dof2 >= 0) c[dof2] = complex<double>(dom->dbc[i].val, 0.0);
    }
    ComplexBCond * cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
    for(i = 0; i < numComplexDirichlet; ++i) {
      int dof2 = dom->dsa->locate(cdbcMRHS[i].nnum,(1 << cdbcMRHS[i].dofnum));
      dof2 = dom->c_dsa->invRCN(dof2);
      if(dof2 >= 0)
        c[dof2] = complex<double>( cdbcMRHS[i].reval, cdbcMRHS[i].imval);
    }

    ComplexD * uel = (ComplexD*) alloca(dom->maxNumDOFs*sizeof(ComplexD));
    int *nds = (int*) alloca(dom->maxNumDOFs*sizeof(int));
    int *dofs = (int*) alloca(dom->maxNumDOFs*sizeof(int));

    for(int iele = 0; iele < numScatter; ++iele) {
      int numEleDofs = scatter[iele]->el->numDofs();
      scatter[iele]->el->nodes(nds);
      scatter[iele]->el->dofs(*dom->dsa,dofs);
      for(i = 0; i < numEleDofs; i++) {
        int hLoc  = dom->c_dsa->locate(nds[i], DofSet::Helm);
        int hLoc1 =   dom->dsa->locate(nds[i], DofSet::Helm);
        if (hLoc >= 0) uel[i] = u[hLoc];
        else if (hLoc1 >=0 ) uel[i] = c[dom->c_dsa->invRCN(hLoc1)];
        else fprintf(stderr,"Error in HData:outputFFP\n");
      }
// RT 03142013: someone messed up the logic above which now requires this
      if (!direction) scatter[iele]->ffp(dom->nodes,ndir,(double*)dir,uel,ffp,direction);
      else
        scatter[iele]->ffp(dom->nodes,dom->numFFPDirections,(double*)dir,uel,ffp,direction);
    }
  }
  else {
    int i;
    for(i = 0; i < ndir; i++) ffp[i] = DComplex(0.0,0.0);

    double *trace = new double[dom->dsa->size()];
    for(i = 0; i < dom->dsa->size(); i++ ) trace[i] = 0.0;

    for(i = 0; i < numSBoundNodes; i++) {
      int ndNum = sBoundNodes[i];
      int dof = dom->dsa->locate(ndNum,DofSet::Helm);
      if(dof > -1) // PJSA
        trace[dof] = 1.0;
    }

    DComplex *c;
    if(dom->numDirichlet+numComplexDirichlet > 0) {
      c = (DComplex *) alloca((dom->numDirichlet+numComplexDirichlet)*sizeof(DComplex));
      for(int iDof = 0; iDof < dom->numDirichlet+numComplexDirichlet; ++iDof)
        c[iDof] = DComplex (0.0, 0.0);

      for(i = 0; i < dom->numDirichlet; ++i) {
        int dof2 = dom->dsa->locate(dom->dbc[i].nnum,(1 << dom->dbc[i].dofnum));
        dof2 = dom->c_dsa->invRCN(dof2);
        if(dof2 >= 0) c[dof2] = dom->dbc[i].val;
      }
      ComplexBCond * cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
      for(i = 0; i < numComplexDirichlet; ++i) {
        int dof2 = dom->dsa->locate(cdbcMRHS[i].nnum,(1 << cdbcMRHS[i].dofnum));
        dof2 = dom->c_dsa->invRCN(dof2);
        if(dof2 >= 0)
          c[dof2] = DComplex(cdbcMRHS[i].reval, cdbcMRHS[i].imval);
      }
    }

    int sBoundLen = 0;
    for(i = dom->dsa->size()-1; i >= 0; i--)
      if(sBoundMap[i] != -1) {
        sBoundLen = sBoundMap[i]+1;
        break;
      }
    if(sBoundLen==0) return;

    int *invSBoundMap = new int[sBoundLen];
    for(i=0; i < dom->dsa->size(); i++) {
      if(sBoundMap[i] != -1) invSBoundMap[sBoundMap[i]] = i;
    }

    DComplex *expPart = (DComplex*) alloca(sBoundLen*sizeof(DComplex));
    DComplex *uPart = (DComplex*) alloca(sBoundLen*sizeof(DComplex));
    DComplex *res1 = (DComplex*) alloca(sBoundLen*sizeof(DComplex));
    DComplex *res2 = (DComplex*) alloca(sBoundLen*sizeof(DComplex));

    for(i = 0; i < sBoundLen; i++) {
      if(dom->c_dsa->getRCN(invSBoundMap[i]) >= 0)
        uPart[i] = coupledScaling*u[dom->c_dsa->getRCN(invSBoundMap[i])];
      else
        uPart[i] = c[dom->c_dsa->invRCN(invSBoundMap[i])];
    }

    int *sBoundNodeMap = new int[sBoundLen];
    for(i = 0; i < dom->numNodes(); i++) {
      int dof = dom->dsa->locate(i,DofSet::Helm);
      int sDof = (dof > -1) ? sBoundMap[dof] : -1; // PJSA
      if (sDof>=0) sBoundNodeMap[sDof] = i;
    }

    for(i = 0; i < sBoundLen; i++) res1[i] = DComplex(0.0,0.0);
    Kss->multcomplex(uPart,res1);

    for(int iDir = 0; iDir < ndir; iDir++) {
      double kDir[3];
      double kappa = geoSource->kappa();
      kDir[0] = -kappa*dir[iDir][0];
      kDir[1] = -kappa*dir[iDir][1];
      kDir[2] = -kappa*dir[iDir][2];
      for(i = 0; i < sBoundLen; i++) {
        Node nd = dom->nodes.getNode(sBoundNodeMap[i]);
        expPart[i] = exp(DComplex(0.0, kDir[0]*nd.x + kDir[1]*nd.y + kDir[2]*nd.z ));
      }
      for(i = 0; i < sBoundLen; i++) res2[i] = DComplex(0.0,0.0);
      Kss->multcomplex(expPart,res2);
      for(i = 0; i < sBoundLen; i++) {
        ffp[iDir] += trace[invSBoundMap[i]]* (-uPart[i]*res2[i]+expPart[i]*res1[i]);
      }
    }
    delete [] trace;
  }
}

template<class Scalar>
void
HData::wError(Domain *dom, double *l2err, double *h1err, double *l2, double *h1, Scalar *u)
{
 Scalar *c = (Scalar *)
          dbg_alloca((dom->numDirichlet+numComplexDirichlet)*sizeof(Scalar));
 int iDof;
 for(iDof = 0; iDof < dom->numDirichlet+numComplexDirichlet; ++iDof)
   ScalarTypes::initScalar(c[iDof], 0.0, 0.0);
 int i;
 for(i=0; i<dom->numDirichlet; ++i) {
   int dof2 = dom->dsa->locate(dom->dbc[i].nnum,(1 << dom->dbc[i].dofnum));
   dof2 = dom->c_dsa->invRCN(dof2);
   if(dof2 >= 0) ScalarTypes::initScalar(c[dof2], dom->dbc[i].val, 0.0);
 }
 ComplexBCond * cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
 for(i=0; i<numComplexDirichlet; ++i) {
    int dof2 = dom->dsa->locate(cdbcMRHS[i].nnum,(1 << cdbcMRHS[i].dofnum));
    dof2 = dom->c_dsa->invRCN(dof2);
    if(dof2 >= 0)
      ScalarTypes::initScalar(c[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
 }

 double *wd = domain->getWaveDirection();

 ComplexD * uel = (ComplexD*) dbg_alloca(dom->maxNumDOFs*sizeof(ComplexD));
 int *heNodes = (int*) dbg_alloca(dom->maxNumDOFs*sizeof(int));
 int *dofs = (int*) dbg_alloca(dom->maxNumDOFs*sizeof(int));
 double * v= (double*) dbg_alloca(dom->maxNumDOFs*dom->maxNumDOFs*sizeof(double));

 int iele;
 for(iele=0; iele < dom->numele; ++iele) {
   auto dofs = (*dom->allDOFs)[iele];
   HelmElement *he = dynamic_cast<HelmElement*>(dom->packedEset[iele]);
   if (he == 0)
     fprintf(stderr,"HData::wErrorreceived non-Helmholtz element.\n");
   else {
     FullSquareMatrix kel = he->acousticm(dom->nodes,v);
     int numEleDofs = dom->packedEset[iele]->numDofs();
     dom->packedEset[iele]->nodes(heNodes);
     dom->packedEset[iele]->dofs(*dom->dsa,dofs.data());
     int i;
     for(i=0;i<numEleDofs;i++) {
       int hLoc  = dom->c_dsa->locate(heNodes[i], DofSet::Helm);
       int hLoc1 =   dom->dsa->locate(heNodes[i], DofSet::Helm);
       if (hLoc >= 0) uel[i] = coupledScaling*u[hLoc];
       else if (hLoc1 >=0 ) uel[i] = c[dom->c_dsa->invRCN(hLoc1)];
       else fprintf(stderr,"Error in HData:wError\n");
     }
     double kappa = dom->packedEset[iele]->getProperty()->kappaHelm; // PJSA 1-15-08
     he->wErrors(dom->nodes,l2err,h1err,l2,h1,uel,kappa,wd);
   }
 }
}


void
HData::checkSommerTypeBC(Domain *dom, const Connectivity *_elemToNode, const Connectivity *_nodeToElem)
{
 if(sommerChecked) return;
 int totEle = dom->numElements();
 int *eleTouch = new int[totEle];
 int *eleCount = new int[totEle];

 auto *packedElemToNode = (_elemToNode) ? _elemToNode : new Connectivity(dom->packedEset.asSet());
 auto *nodeToPackedElem = (_nodeToElem) ? _nodeToElem : packedElemToNode->alloc_reverse();
// int pos = 0, neg = 0;
 int i;
 for (i=0;i<totEle;i++) eleTouch[i] = -1;
 int iSEle;
 for(iSEle =0; iSEle < numSommer; ++iSEle) {
   int s = sommer[iSEle]->findAndSetEle(dom->nodes,dom->packedEset,
                                  nodeToPackedElem, eleTouch, eleCount, iSEle);
   if (s>0) {
//     if(verboseFlag) fprintf(stderr,"Flipping %d of sommer.\n",iSEle);
//     pos++;
     sommer[iSEle]->flipNormal();
   }
//   else { 
//     neg++;
//   }
 }
 // if (numSommer>0) fprintf(stderr,"Sommer bc had %d positive and %d negative normals.\n",pos,neg);
// pos = neg = 0;
 for (i=0;i<totEle;i++) eleTouch[i] = -1;
 for(iSEle =0; iSEle < numNeum; ++iSEle) {
// RT
   int s = neum[iSEle]->findAndSetEle(dom->nodes,dom->packedEset,
//                                      nodeToPackedElem, eleTouch, eleCount, iSEle,2);
                                      nodeToPackedElem, eleTouch, eleCount, iSEle);
   if (s>0) {
//     if(verboseFlag) fprintf(stderr,"Flipping %d of neum.\n",iSEle);
//     pos++;
     neum[iSEle]->flipNormal();
   }
//   else neg++;
 }
// if (numNeum>0) fprintf(stderr,"Neumann bc had %d positive and %d negative normals.\n",pos,neg);

// pos = neg = 0;
 for (i=0;i<totEle;i++) eleTouch[i] = -1;
 for(iSEle =0; iSEle < numScatter; ++iSEle) {
   int s = scatter[iSEle]->findAndSetEle(dom->nodes,dom->packedEset,
                                         nodeToPackedElem, eleTouch, eleCount, iSEle,1);
   if (s>0) {
//     if(verboseFlag) fprintf(stderr,"Flipping %d of sca.\n",iSEle);
//     pos++;
     scatter[iSEle]->flipNormal();
   }
//   else neg++;
 }
 // if (numScatter>0) fprintf(stderr,"Scatterer had %d positive and %d negative normals.\n",pos,neg);

// pos = neg = 0;
 if (!( (domain->solInfo().HEV) && (domain->probType() != SolverInfo::Modal) )) {
 for (i=0;i<totEle;i++) eleTouch[i] = -1;
 for(iSEle =0; iSEle < numWet; ++iSEle) {
   int s = wet[iSEle]->findAndSetBothEle(dom->nodes,dom->packedEset,
                            nodeToPackedElem, eleTouch, eleCount, iSEle);
   if (s>0) {
//     if(verboseFlag) fprintf(stderr,"Flipping %d of wet.\n",iSEle);
//     pos++;
     wet[iSEle]->flipNormal();
   }
//   else neg++;
 }
 }
 // if (numWet>0) fprintf(stderr,"Wet interface had %d positive and %d negative normals.\n",pos,neg);

 delete[] eleTouch;
 delete[] eleCount;
 if(!_nodeToElem) delete nodeToPackedElem;
 if(!_elemToNode) delete packedElemToNode;
 sommerChecked = true;
}

void
HData::addFrequencies1(double w0, double deltaw, int n_deltaw, bool minus_flag)
{
  // for multiple frequency sweep analysis, build list of circular frequencies
  // this function is used for extrapolation eg taylor or 1 point pade
  if(w0 < 0.0) {
    filePrint(stderr, " *** WARNING: w0 is negative \n");
    return;
  }

  if((frequencies == 0) || (coarse_frequencies == 0)) {
    initFreqSweep(w0);
  }
  coarse_frequencies->push_back(w0);

  // add w0 - n_deltaw*deltaw <= w <= w0 + n_deltaw*deltaw
  double w;
  for(int i=1; i<=n_deltaw; ++i) {
    if(minus_flag) {  // sweep in both directions
      w = w0-double(i)*deltaw;
      if(w >= 0.0) frequencies->push_back(w);  // ignore negative freq's
    }
    w = w0+double(i)*deltaw;
    frequencies->push_back(w);
  }
  numFrequencies = frequencies->size() + coarse_frequencies->size();
}

void
HData::addFrequencies2(double w0, double w1, int n_deltaw, bool first_flag)
{
  // for multiple frequency sweep analysis, build list of circular frequencies
  // this function is used for interpolation eg 2 point pade or coarse grid frequency sweep
  if((w0 < 0.0) || (w1 < 0.0)) {
    filePrint(stderr, " *** WARNING: w0 or w1 is negative \n");
    return;
  }

  if((frequencies == 0) || (coarse_frequencies == 0)) {
    initFreqSweep(w0);
  }
  if(first_flag) coarse_frequencies->push_back(w0);
  coarse_frequencies->push_back(w1);

  // add w0 < w < w1
  double w;
  double deltaw = (w1-w0)/n_deltaw;
  for(int i=1; i<n_deltaw; ++i) {
    w = w0+double(i)*deltaw;
    frequencies->push_back(w);
  }
  numFrequencies = frequencies->size() + coarse_frequencies->size();
}

void
HData::addFrequencies(double w_start, double w_end, int n_coarse, int n_fine)
{
  // for multiple LHS frequency sweep analysis, build list of circular frequencies
  // note: n_coarse is the number of coarse points, not the number of coarse increments
  //       n_fine is the number of fine increments per coarse increment
  if((w_start < 0.0) || (w_end < 0.0)) {
    filePrint(stderr, " *** WARNING: w_start or w_end is negative \n");
    return;
  }

  if((frequencies == 0) || (coarse_frequencies == 0)) {
    initFreqSweep(w_start);
  }

  double w;
  double deltaw_coarse = (w_end - w_start)/(n_coarse-1);
  double deltaw_fine = deltaw_coarse/n_fine;
  for(int i = 0; i < (n_coarse-1); ++i) {
    w = w_start+double(i)*deltaw_coarse;
    coarse_frequencies->push_back(w);
    for(int j = 1; j < n_fine; ++j) {
      w += deltaw_fine;
      frequencies->push_back(w);
    }
  }
  coarse_frequencies->push_back(w_end);
  numFrequencies = frequencies->size() + coarse_frequencies->size();
}

void
HData::addFrequencies(double w_i, int n_fine)
{
  // for multiple LHS frequency sweep, add new coarse freq w_i and all the fine freqs 
  // between w_(i-1) and w_i
  double w = coarse_frequencies->back();

  if(w_i < 0.0) {
    filePrint(stderr, " *** WARNING: w_i is negative \n");
    return;
  }

  if((frequencies == 0) || (coarse_frequencies == 0)) {
    initFreqSweep(w);
  }

  double deltaw_coarse = w_i - w;
  double deltaw_fine = deltaw_coarse/n_fine;
  for(int j = 1; j < n_fine; ++j) {
    w += deltaw_fine;
    frequencies->push_back(w);
  }
  coarse_frequencies->push_back(w_i);
  numFrequencies = frequencies->size() + coarse_frequencies->size();
}

void
HData::addFrequency(double _w)
{
  // for multiple frequency sweep analysis, add single freq to list
  if(frequencies == 0) frequencies = new std::list<double>();
  frequencies->push_back(_w);
  numFrequencies = frequencies->size() + coarse_frequencies->size();
}

void
HData::addCoarseFrequency(double _w)
{
  // for multiple frequency sweep analysis, add single coarse freq/wavenumber to list
  if(coarse_frequencies == 0) {
    initFreqSweep(_w);
  }
  coarse_frequencies->push_back(_w);
  numFrequencies = frequencies->size() + coarse_frequencies->size();
}

void
HData::initFreqSweep(double w0)
{
  if(frequencies == 0) frequencies = new std::list<double>();
  if(coarse_frequencies == 0) coarse_frequencies = new std::list<double>();
  domain->solInfo().doFreqSweep = true;
}

HData::~HData()
{
  if(frequencies) delete frequencies;
  if(coarse_frequencies) delete coarse_frequencies;
  if(subScaToSca) delete [] subScaToSca;
}
