#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/linkfc.h>
#include <Math.d/matrix.h>
#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>
#include <Math.d/IntFullM.h>

#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

extern "C"      {
   void _FORTRAN(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);
}

void
Rbm::init()
{
  myMemory = 0;
  Rmat     = 0;
  Amat     = 0;
  allRbm   = 0;
  nRbmPerComp    = 0;
  firstDofOfComp = 0;
  numDofPerComp  = 0;
  numBC          = 0;
  Zstar  = 0;
  Rc     = 0;
  Zmpc   = 0;
  xyzRot = 0;
  ngrbm  = 0;
  nComponents = 0;
  grbm = 0;
  cgrbm = 0;
}

Rbm::Rbm(Vector *_grbm, int _ngrbm, int _numUncon, int _myMemory)
{
  init();
  numUncon = _numUncon;
  ngrbm    = _ngrbm;
  grbm     = _grbm;
  myMemory = _myMemory;
  nComponents = 1;
  firstDofOfComp = new int[nComponents+1];
  firstDofOfComp[0] = 0;
  firstDofOfComp[1] = numUncon;
  nRbmPerComp    = new int[nComponents];
  nRbmPerComp[0] = ngrbm;
  numDofPerComp  = new int[nComponents];
  numDofPerComp[0] = numUncon;
}

Rbm::Rbm(ComplexVector *_grbm, int _ngrbm, int _numUncon, int _myMemory)
{
  init();
  numUncon = _numUncon;
  ngrbm    = _ngrbm;
  cgrbm     = _grbm;
  myMemory = _myMemory;
  nComponents = 1;
  firstDofOfComp = new int[nComponents+1];
  firstDofOfComp[0] = 0;
  firstDofOfComp[1] = numUncon;
  nRbmPerComp    = new int[nComponents];
  nRbmPerComp[0] = ngrbm;
  numDofPerComp  = new int[nComponents];
  numDofPerComp[0] = numUncon;
}

Rbm::Rbm(DofSetArray *_dsa, ConstrainedDSA *_c_dsa)
{
  init();
  myMemory = 1;
  // Compute rbm (=zero energy mode) in case of heat 
  dsa = _dsa;
  c_dsa = _c_dsa;

  int ndofs = dsa->size();       // Number of degrees of freedom
  numUncon = c_dsa->size();      // Number of unconstraint DOF
  int numCon = ndofs - numUncon; 

  // If at least one Dirichlet BC condition, then no need for rbm
  if(numCon > 0) {
    fprintf(stderr, "... %d Dirichlet BCs were found ...\n",numCon);
    ngrbm=0;
    return;
  }

  nComponents = 1;
  allRbm = new Vector*[nComponents];
  allRbm[0]   = new Vector[ndofs];
  int j;
  ngrbm=1; 

  firstDofOfComp = new int[nComponents+1];
  firstDofOfComp[0] = 0;
  firstDofOfComp[1] = numUncon;
  nRbmPerComp    = new int[nComponents];
  nRbmPerComp[0] = ngrbm;
  numDofPerComp  = new int[nComponents];
  numDofPerComp[0] = numUncon;

  for(j=0; j<numUncon; ++j) {
    allRbm[0][j] = 1.;
  }

  if (domain->solInfo().sloshing)  {
    Vector v1(numUncon,0.0);
    grbm = new Vector[ngrbm];
    grbm[0] = v1;
    for(j=0; j<numUncon; ++j)  {
      grbm[0][j] = 1.0;
    }
  }
}

Rbm::Rbm(DofSetArray *_dsa, ConstrainedDSA *_c_dsa, const CoordSet &cs,
         double _tolgrb, compStruct &_comp, IntFullM *_fm)
    : U(6,6)
{
 init();
 myMemory = 1;
 cornerModes = _fm;
 c_dsa = _c_dsa;
 dsa = _dsa;
 comp = &_comp;
 tolgrb = _tolgrb;

 // ... DECLARATIONS
 int i, n, inode;
 numUncon = c_dsa->size();
 if(numUncon == 0) { 
   ngrbm=0;
   return;
 }

 // ... SET NUMBER OF COMPONENTS
 nComponents = comp->numComp;
 allRbm = new Vector*[nComponents];
 int *compSize  = (int *) dbg_alloca(sizeof(int)*nComponents);
 nRbmPerComp    = new int[nComponents];
 firstDofOfComp = new int[nComponents+1];
 numDofPerComp  = new int[nComponents];
 numBC          = new int[nComponents];
 xyzRot         = new double*[nComponents];

 // ... DETERMINE NUMBER OF DIRICHLET BC PER COMPONENT,
 // ... NUMBER OF RBM PER COMPONENT
 // ... AND NUMBER OF DOFS PER COMPONENT
 for(n=0; n<nComponents; ++n) {
   int totalbc = 0, totaldofs = 0;
   for(i=comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
     inode = comp->order[i];
     int dofsPerNode    =   dsa->weight(inode);
     int conDofsPerNode = c_dsa->weight(inode);
     totalbc   += dofsPerNode - conDofsPerNode;
     totaldofs += dofsPerNode;
   }
   numBC[n]    = totalbc; 
   compSize[n] = totaldofs; 
   allRbm[n] = new Vector[6];
   xyzRot[n]   = new double[3];
   numDofPerComp[n] = totaldofs-totalbc;
   firstDofOfComp[n] = (n > 0) ? firstDofOfComp[n-1]+numDofPerComp[n-1] : 0;
   nRbmPerComp[n] = 6;
 }
 firstDofOfComp[nComponents] = numUncon; 

 // ... INITIALIZE ALLRBM:
 initialize(nComponents, nRbmPerComp, numUncon);

 Rmat = new FullM(dsa->size(),6);
 Amat = new FullM(6,6,0.0);

 computeRbms(cs);
}

// Direct Sky/Sparse or Eigen solver with LMPCs
Rbm::Rbm(DofSetArray *_dsa, const ConstrainedDSA *_c_dsa, const CoordSet &cs,
         double _tolgrb, const compStruct &_comp, int numMPC,
         ResizeArray<LMPCons *> &mpc, IntFullM *_fm)
    : U(6*_comp.numComp)
{
 init();
 myMemory = 1;
 cornerModes = _fm;
 c_dsa = _c_dsa;
 dsa = _dsa;
 comp = &_comp;
 tolgrb = _tolgrb;

 // ... DECLARATIONS
 numUncon = c_dsa->size();
 if(numUncon == 0) {
   ngrbm=0;
   return;
 }

 // ... SET NUMBER OF COMPONENTS
 // when there are MPCs, use only one component
 nComponents = 1;
 allRbm = new Vector*[nComponents];
 int *compSize  = (int *) dbg_alloca(sizeof(int)*nComponents);
 nRbmPerComp    = new int[nComponents];
 firstDofOfComp = new int[nComponents+1];
 numDofPerComp  = new int[nComponents];
 numBC          = new int[nComponents];
 xyzRot         = new double*[nComponents];

 // ... DETERMINE NUMBER OF DIRICHLET BC PER COMPONENT,
 // ... NUMBER OF RBM PER COMPONENT
 // ... AND NUMBER OF DOFS PER COMPONENT
 int totaldofs = dsa->size();
 int totalbc = totaldofs - c_dsa->size();
 numBC[0]    = totalbc;
 compSize[0] = totaldofs;
 allRbm[0] = new Vector[6*comp->numComp];
 xyzRot[0]   = new double[3];
 numDofPerComp[0] = totaldofs-totalbc;
 firstDofOfComp[0] = 0;
 nRbmPerComp[0] = 6*comp->numComp;
 firstDofOfComp[nComponents] = numUncon;

 // ... INITIALIZE ALLRBM:
 initialize(nComponents, nRbmPerComp, numUncon);
 
 Rmat = new FullM(dsa->size(),6*comp->numComp);
 Amat = new FullM(6*comp->numComp,6*comp->numComp,0.0);

 computeRbms(cs, numMPC, mpc);
}

// new constructor to be called by SubDomain::makeZstarAndR()
// works for 1 component only
Rbm::Rbm(const DofSetArray *_dsa, const ConstrainedDSA *_c_dsa, const CoordSet &cs,
         double _tolgrb, double *centroid,
         const std::vector<int> &cornerNodes, int numCRN, int numCRNdof, const std::vector<DofSet> &cornerDofs,
         int numMPC, const std::vector<std::unique_ptr<SubLMPCons<double> > > &mpc)
{
 init();
 myMemory = 1;
 c_dsa = _c_dsa;
 dsa = _dsa;
 tolgrb = _tolgrb;

 computeRbms(cs, centroid, cornerNodes, numCRN, numCRNdof, cornerDofs, numMPC, mpc); 
}

void
Rbm::computeRbms(const CoordSet &cs, double *centroid, const std::vector<int> &cornerNodes,
                 int numCRN, int numCRNdof, const std::vector<DofSet> &cornerDofs,
                 int numMPC, const std::vector<std::unique_ptr<SubLMPCons<double> > > &mpc)
{
  int i,j,k;

  // check for active translations
  int rows[6] = {-1, -1, -1, -1, -1, -1};
  int dofs[6] = {DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp, DofSet::Xrot, DofSet::Yrot, DofSet::Zrot};
  for(i=0; i<dsa->numNodes(); ++i)
    for(j=0; j<6; ++j)  // allows for rotational springs with no active translational dofs
      if(rows[j] == -1)
        if(dsa->locate(i, dofs[j]) > -1) rows[j] = j;
  // add required rotations
  if((rows[0] != -1) && (rows[1] != -1)) rows[5] = 5;
  if((rows[0] != -1) && (rows[2] != -1)) rows[4] = 4;
  if((rows[1] != -1) && (rows[2] != -1)) rows[3] = 3;

  // determine nrow and ncol
  int c_rows[6];
  int ncol = 0;
  for(j=0; j<6; ++j)
    if(rows[j] > -1) c_rows[ncol++] = rows[j];
  int nrow = dsa->size() - c_dsa->size();  // number of constrained dofs
  if(ncol == 0) std::cerr << " *** WARNING: found no active dimensions in Rbm::computeRbms(...) \n";
  
  // intialize R and Z
  R.setNewSize(c_dsa->size(), ncol, 0.0);
  FullM Z(nrow, ncol);
  double cgmax = 0.0;

  // build R matrix containing geometric rbms (includes constrained but not inactive dofs)
  // and Z = E^t * R
  for(i=0; i<cs.size(); ++i) {
    auto &nd = cs.getNode(i);
    double x = (nd.x - centroid[0]); 
    double y = (nd.y - centroid[1]); 
    double z = (nd.z - centroid[2]); 
    double rbm[6][6] = {
           { 1.0, 0.0, 0.0, 0.0,   z,  -y },
           { 0.0, 1.0, 0.0,  -z, 0.0,   x },
           { 0.0, 0.0, 1.0,   y,  -x, 0.0 },
           { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
           { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
           { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
    }; 
    if(nd.cd != 0) {
#ifdef USE_EIGEN3
       Eigen::Map<Eigen::Matrix<double,6,6,Eigen::RowMajor> > Ri(&rbm[0][0]);

       // transform Ri from basic to DOF_FRM coordinates
       NFrameData *nfd = geoSource->getNFrames();
       Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&nfd[nd.cd].frame[0][0]);;

       Ri.topLeftCorner<3,3>() = T;
       Ri.topRightCorner<3,3>() = (T*Ri.topRightCorner<3,3>()).eval();
       Ri.bottomRightCorner<3,3>() = T;
#else
       std::cerr << "USE_EIGEN3 is not defined here in Rbm::computeRbms\n";
#endif
    }

    for(j=0; j<6; ++j) {
      int c_dof = c_dsa->locate(i,dofs[j]);
      if(c_dof > -1) 
        for(k=0; k<ncol; ++k) 
          R[c_dof][k] = rbm[j][c_rows[k]];
      int dof = dsa->locate(i,dofs[j]);
      if(dof > -1) {  // add SPC constraint to Z
        int cn = c_dsa->invRCN(dof);
        for(k=0; k<ncol; ++k) 
          if(cn > -1) {
            Z[cn][k] = rbm[j][c_rows[k]];
            if(std::abs(Z[cn][k]) > cgmax) cgmax = std::abs(Z[cn][k]);
          }
      }
    }
  }

  // now add MPC constraints to Zmpc
  if(numMPC > 0) {
    Zmpc = new FullM(numMPC, ncol);
    Zmpc->zero();
    for(i=0; i<numMPC; ++i) {
      for(j=0; j<mpc[i]->nterms; ++j) {
        int c_dof = c_dsa->locate(mpc[i]->terms[j].nnum,
                                 (1 << mpc[i]->terms[j].dofnum));
        if(c_dof < 0) continue;
        for(k=0; k<ncol; ++k) 
          (*Zmpc)[i][k] += R[c_dof][k]*mpc[i]->terms[j].coef;
      }
    }
  }
  
  // get Rc
  int offset = 0;
  Rc = new FullM(numCRNdof, ncol);
  Rc->zero();
  for(i = 0; i < numCRN; ++i) {
    int num[DofSet::max_known_dof];
    int count = c_dsa->number(cornerNodes[i], cornerDofs[i].list(), num);
    for(j = 0; j < count; ++j)
      for(int k = 0; k < ncol; ++k)
        (*Rc)[j+offset][k] = R[num[j]][k];
    offset += cornerDofs[i].count();
  }

  if(nrow == 0) { // no dirichlet boundary conditions
    Zstar = new FullM(0, ncol);  
    ngrbm = ncol;
  }
  else { // do a SVD on Z to form Zstar:
    int rank = 0;
    FullM Utmp(ncol, ncol);
    Utmp.zero();
    singularValueDecomposition(Z, Utmp, ncol, nrow, cgmax, ngrbm, rank);
    Zstar = new FullM(Utmp, rank, 0, ncol, 0);
  }
}

void
Rbm::computeRbms(const CoordSet &cs)
{
  int debug = 0;
  int rank  = 0;
  int i, inode, j, k, cn, ncol, xloc, yloc, zloc, xrot, yrot, zrot;
  int position[6];
  double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  int numdofs = dsa->size();
  R = *Rmat;
  U.zero();

// ... COUNT THE TOTAL NUMBER OF DIRICHLET BCs

  int n, totNumBC = 0;
  for(n=0; n<nComponents; ++n) totNumBC += numBC[n];

// ... BEGIN COMPONENT LOOP:

 for(n=0; n<nComponents; ++n) {

   R.zero();
   ncol     = nRbmPerComp[n];
   int nrow = numBC[n];
   bool setzero = true;

// ... CALCULATE R MATRIX CONTAINING GEOMETRIC RBM
   for(i = comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
     inode = comp->order[i];

     if(dsa->firstdof(inode) == -1) continue;

     if(cs[inode] == 0) continue;

     auto &nd = cs.getNode(inode);

     if(setzero) { // make sure that node0 exists before trying to access its coordinates
       if(domain->solInfo().grbm_ref.empty() /*|| nComponents > 1*/) {
         x0 = xyzRot[n][0] = nd.x;
         y0 = xyzRot[n][1] = nd.y;
         z0 = xyzRot[n][2] = nd.z;
         fprintf(stderr," ... Using first node of component as reference pt for the rotations: %3.2e %3.2e %3.2e\n",xyzRot[0][0],xyzRot[0][1],xyzRot[0][2]);
       }
       else {
         x0 = xyzRot[n][0] = domain->solInfo().grbm_ref[0];
         y0 = xyzRot[n][1] = domain->solInfo().grbm_ref[1];
         z0 = xyzRot[n][2] = domain->solInfo().grbm_ref[2];
         fprintf(stderr," ... Using specified X,Y,Z coordinates as reference pt for the rotations: %3.2e %3.2e %3.2e\n",xyzRot[0][0],xyzRot[0][1],xyzRot[0][2]);
       }
       setzero = false;
     }

     double x = nd.x - x0;
     double y = nd.y - y0;
     double z = nd.z - z0;

     xloc = dsa->locate( inode, DofSet::Xdisp );
     yloc = dsa->locate( inode, DofSet::Ydisp );
     zloc = dsa->locate( inode, DofSet::Zdisp );
     xrot = dsa->locate( inode, DofSet::Xrot  );
     yrot = dsa->locate( inode, DofSet::Yrot  );
     zrot = dsa->locate( inode, DofSet::Zrot  );

     if(nd.cd == 0) {

       if(xloc >= 0) {
         R[xloc][0] = 1.0;
         R[xloc][1] = 0.0;
         R[xloc][2] = 0.0;
         R[xloc][3] = 0.0;
         R[xloc][4] =   z;
         R[xloc][5] =  -y;
       }
       if(yloc >= 0) {
         R[yloc][0] = 0.0;
         R[yloc][1] = 1.0;
         R[yloc][2] = 0.0;
         R[yloc][3] =  -z;
         R[yloc][4] = 0.0;
         R[yloc][5] =   x;
       }
       if(zloc >= 0) {
         R[zloc][0] = 0.0;
         R[zloc][1] = 0.0;
         R[zloc][2] = 1.0;
         R[zloc][3] =   y;
         R[zloc][4] =  -x;
         R[zloc][5] = 0.0;
       }
       if(xrot >= 0) {
         R[xrot][0] = 0.0;
         R[xrot][1] = 0.0;
         R[xrot][2] = 0.0;
         R[xrot][3] = 1.0;
         R[xrot][4] = 0.0;
         R[xrot][5] = 0.0;
       }
       if(yrot >= 0) {
         R[yrot][0] = 0.0;
         R[yrot][1] = 0.0;
         R[yrot][2] = 0.0;
         R[yrot][3] = 0.0;
         R[yrot][4] = 1.0;
         R[yrot][5] = 0.0;
       }
       if(zrot >= 0) {
         R[zrot][0] = 0.0;
         R[zrot][1] = 0.0;
         R[zrot][2] = 0.0;
         R[zrot][3] = 0.0;
         R[zrot][4] = 0.0;
         R[zrot][5] = 1.0;
       }
     }
     else {
       Eigen::Matrix<double,6,6> Ri;
       Ri << 1, 0, 0,  0,  z, -y,
             0, 1, 0, -z,  0,  x,
             0, 0, 1,  y, -x,  0,
             0, 0, 0,  1,  0,  0,
             0, 0, 0,  0,  1,  0,
             0, 0, 0,  0,  0,  1;

       // transform Ri from basic to DOF_FRM coordinates
       NFrameData *nfd = geoSource->getNFrames();
       Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&nfd[nd.cd].frame[0][0]);;

       Ri.topLeftCorner<3,3>() = T;
       Ri.topRightCorner<3,3>() = (T*Ri.topRightCorner<3,3>()).eval();
       Ri.bottomRightCorner<3,3>() = T;

       if(xloc >= 0)
         for(int j=0; j<6; ++j) R[xloc][j] = Ri(0,j);
       if(yloc >= 0)
         for(int j=0; j<6; ++j) R[yloc][j] = Ri(1,j);
       if(zloc >= 0)
         for(int j=0; j<6; ++j) R[zloc][j] = Ri(2,j);
       if(xrot >= 0)
         for(int j=0; j<6; ++j) R[xrot][j] = Ri(3,j);
       if(yrot >= 0)
         for(int j=0; j<6; ++j) R[yrot][j] = Ri(4,j);
       if(zrot >= 0)
         for(int j=0; j<6; ++j) R[zrot][j] = Ri(5,j);
     }
   }

//  eliminate redundent rbms for 1D and 2D cases

    R.transposeMult(R, *Amat);
    FullM &A = *Amat;

    int dimA = A.dim();

    double *diagA = (double *) dbg_alloca(sizeof(double)*dimA);

    for(i = 0; i < dimA; ++i) {
      diagA[i] = A[i][i];
      if(diagA[i] > 1e-12)
        diagA[i] = 1.0/sqrt(diagA[i]);
      else
        diagA[i] = 1.0;
    }
    for(i = 0; i < dimA; ++i)
      for(j = 0; j < dimA; ++j)
        A[i][j] *= diagA[i]*diagA[j];

    double max_value = A.maxAbs();  // max() finds maximum value of matrix A

    if(debug) A.print("--- A = Rt*R matrix ---");

// ... PERFORM SVD ON R^t*R

    singularValueDecomposition(A, U, dimA, dimA, max_value, ngrbm, rank);

    if(rank == dimA && numBC[n] == 0) {
      if(totNumBC == 0) {
        for(i=0; i<ncol; ++i)
          for(j=0; j<numUncon; ++j)
            allRbm[n][i][j] = R[j][i];
      }
      else {
        for(i=0; i<numdofs; ++i) {
          int cn = c_dsa->getRCN(i);
          if(cn >= 0) {
            for(j=0; j<ncol; ++j)
              allRbm[n][j][cn] = R[i][j];
          }
        }
      }
      nRbmPerComp[n] = ncol;
      continue;
    }

// ... COPY NECESSARY ROWS OF U TO USTAR:

    FullM Ustar(U,rank,0,ncol,0);
    FullM Rstar;
    if(rank == 6)
      Rstar = R;
    else
      Rstar = R%Ustar;

// ... FOR A FREE-FREE PROBLEM, RETURN GEOMETRIC RIGID BODY MODES.

    if(numBC[n] == 0) {
      if(totNumBC == 0) {
        for(i=0; i<rank; ++i)
          for(j=0; j<numUncon; ++j)
            allRbm[n][i][j] = Rstar[j][i];
      }
      else {
        for(i=0; i<numdofs; ++i) {
          int cn = c_dsa->getRCN(i);
          if(cn >= 0) {
            for(j=0; j<rank; ++j)
              allRbm[n][j][cn] = Rstar[i][j];
          }
        }
      }
      nRbmPerComp[n] = rank;
      continue;
    }

    ncol = rank;

    if(debug) {
      R.print("--- R matrix ---");
      U.print("--- U matrix ---");
      Rstar.print("--- Rstar ---");
    }

// ... EXTRACT Z MATRIX CONTAINING CONSTRAINED GRBM
// ... AND CALCULATE MAXIMUM ENTRY PER COMPONENT
    FullM Z(nrow,ncol);

    double cgmax = 0.0;

    int offset = 0;
    int ii;
    for(ii=0; ii<n; ++ii) {
      offset += numBC[ii];
    }

    for(i=comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
      inode = comp->order[i];

      position[0] = dsa->locate( inode, DofSet::Xdisp);
      position[1] = dsa->locate( inode, DofSet::Ydisp);
      position[2] = dsa->locate( inode, DofSet::Zdisp);
      position[3] = dsa->locate( inode, DofSet::Xrot);
      position[4] = dsa->locate( inode, DofSet::Yrot);
      position[5] = dsa->locate( inode, DofSet::Zrot);

      for(j=0; j<6; ++j) {
        if(position[j] >= 0) {
          cn = c_dsa->invRCN(position[j]);
          if(cn >= 0) {
            int k;
            for(k=0; k<ncol; ++k) {
              Z[cn-offset][k] = Rstar[position[j]][k];
              if(std::abs(Z[cn-offset][k]) > cgmax) cgmax = std::abs(Z[cn-offset][k]);
            }
          }
        }
      }
    }

// ... PERFORM A SVD OPERATION ON Z:

   if(debug) Z.print("----- Z matrix -----");

   U.setNewSize(ncol,ncol);

   singularValueDecomposition(Z, U, ncol, nrow, cgmax, ngrbm, rank);

// ... CONSTRUCT NECESSARY VECTORS FOR Vstar:

   FullM Vstar(U,ngrbm,rank,ncol,0);  // i think Vstar = Q transpose

// ... MATRIX-MATRIX MULTIPLY [GRBM] = [R][V]:

   FullM result = Rstar%Vstar;   // result is Rstar

// ... CONFORM RIGID BODY MODES WITH RESPECT TO K MATRIX

    for(i=0; i<numdofs; ++i) {
      int cn = c_dsa->getRCN(i);
      if(cn >= 0) {
        for(j=0; j<ngrbm; ++j)
          allRbm[n][j][cn] = result[i][j];
      }
    }
    nRbmPerComp[n] = ngrbm;

// ... DEBUGGING PRINT STATEMENTS ...
    if(debug) {
           U.print("----- U matrix       -----");
       Vstar.print("----- Vstar matrix   -----");
           R.print("----- R matrix       -----");
      result.print("----- Results matrix -----");
    }

 } // ... END OF COMPONENT LOOP

 Vector v1(numUncon,0.0);

 ngrbm = 0;
 for(i=0; i<nComponents; ++i)
   ngrbm += nRbmPerComp[i];

 int addModes = ( (cornerModes) ? cornerModes->numCol() : 0);
 grbm = new Vector[ngrbm + addModes];

 for(i=0; i<ngrbm+addModes; ++i)
   grbm[i] = v1;

// ... COPY COMPONENT RBMS TO ONE MATRIX GRBM
 int cnt = 0;
 for(n=0; n<nComponents; ++n)
   for(j=0; j<nRbmPerComp[n]; ++j) {
     for(k=0; k<numUncon; ++k)
       grbm[cnt][k] = allRbm[n][j][k];
     cnt += 1;
   }

 if(debug) {
   fprintf(stderr, "%d components rigid body modes were found.\n",nComponents);
   print();
 }

 for(i = 0; i < addModes; ++i)
   grbm[ngrbm+i][(*cornerModes)[2][i]] = 1.0;

 ngrbm += addModes;
}

void
Rbm::computeRbms(const CoordSet& cs, int numMPC, ResizeArray<LMPCons *> &mpc)
{
 int debug = 0;
 int rank  = 0;
 int n, i, inode, j, k, cn, ncol, nrow, xloc, yloc, zloc, xrot, yrot, zrot;
 int position[6];
 double x0 = 0.0, y0 = 0.0, z0 = 0.0;
 int numdofs = dsa->size();
 R = *Rmat;
 R.zero();
 U.zero();

 // ... BEGIN COMPONENT LOOP:
 nrow = numBC[0] + numMPC;
 ncol = comp->numComp*6;
//#define USE_GLOBAL_CG //HB: to use the global cg as reference pt for the rotations
#ifdef USE_GLOBAL_CG
 int count = 0;
 for(n=0; n<comp->numComp; ++n) {
   for(i = comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
     inode = comp->order[i];
     if(dsa->firstdof(inode) == -1) continue;
     if(cs[inode] == 0) continue;
     Node &nd = cs.getNode(inode); 
     x0 += nd.x; y0 += nd.y; z0 += nd.z;
     count++;
   } 
 }
 xyzRot[0][0] = x0/count;
 xyzRot[0][1] = y0/count;
 xyzRot[0][2] = z0/count;
 //fprintf(stderr," ... Use global cg as reference pt for the rotations: %3.2e %3.2e %3.2e\n",xyzRot[0][0],xyzRot[0][1],xyzRot[0][2]); 
#else // use first node as reference pt for the rotations
 if(domain->solInfo().grbm_ref.empty()) {
   auto &nd = cs.getNode(comp->order[comp->xcomp[0]]);
   xyzRot[0][0] = nd.x; xyzRot[0][1] = nd.y; xyzRot[0][2] = nd.z;
   fprintf(stderr," ... Using first node of first component as reference pt for the rotations: %3.2e %3.2e %3.2e\n",xyzRot[0][0],xyzRot[0][1],xyzRot[0][2]);
 }
 else {
   xyzRot[0][0] = domain->solInfo().grbm_ref[0];
   xyzRot[0][1] = domain->solInfo().grbm_ref[1];
   xyzRot[0][2] = domain->solInfo().grbm_ref[2];
   fprintf(stderr," ... Using specified X,Y,Z coordinates as reference pt for the rotations: %3.2e %3.2e %3.2e\n",xyzRot[0][0],xyzRot[0][1],xyzRot[0][2]);
 }
#endif

 for(n=0; n<comp->numComp; ++n) {

   // ... CALCULATE R MATRIX CONTAINING GEOMETRIC RBM
   for(i = comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
     inode = comp->order[i];
     if(dsa->firstdof(inode) == -1) continue;
     if(cs[inode] == 0) continue;

     auto &nd = cs.getNode(inode);
     
     double x = nd.x - xyzRot[0][0]; // use same reference pt for all the components
     double y = nd.y - xyzRot[0][1];
     double z = nd.z - xyzRot[0][2];

     xloc = dsa->locate(inode, DofSet::Xdisp);
     yloc = dsa->locate(inode, DofSet::Ydisp);
     zloc = dsa->locate(inode, DofSet::Zdisp);
     xrot = dsa->locate(inode, DofSet::Xrot);
     yrot = dsa->locate(inode, DofSet::Yrot);
     zrot = dsa->locate(inode, DofSet::Zrot);

     if(nd.cd == 0) {

       if(xloc >= 0) {
         R[xloc][6*n+0] = 1.0;
         R[xloc][6*n+1] = 0.0;
         R[xloc][6*n+2] = 0.0;
         R[xloc][6*n+3] = 0.0;
         R[xloc][6*n+4] =   z;
         R[xloc][6*n+5] =  -y;
       }
       if(yloc >= 0) {
         R[yloc][6*n+0] = 0.0;
         R[yloc][6*n+1] = 1.0;
         R[yloc][6*n+2] = 0.0;
         R[yloc][6*n+3] =  -z;
         R[yloc][6*n+4] = 0.0;
         R[yloc][6*n+5] =   x;
       }
       if(zloc >= 0) {
         R[zloc][6*n+0] = 0.0;
         R[zloc][6*n+1] = 0.0;
         R[zloc][6*n+2] = 1.0;
         R[zloc][6*n+3] =   y;
         R[zloc][6*n+4] =  -x;
         R[zloc][6*n+5] = 0.0;
       }
       if(xrot >= 0) {
         R[xrot][6*n+0] = 0.0;
         R[xrot][6*n+1] = 0.0;
         R[xrot][6*n+2] = 0.0;
         R[xrot][6*n+3] = 1.0;
         R[xrot][6*n+4] = 0.0;
         R[xrot][6*n+5] = 0.0;
       }
       if(yrot >= 0) {
         R[yrot][6*n+0] = 0.0;
         R[yrot][6*n+1] = 0.0;
         R[yrot][6*n+2] = 0.0;
         R[yrot][6*n+3] = 0.0;
         R[yrot][6*n+4] = 1.0;
         R[yrot][6*n+5] = 0.0;
       }
       if(zrot >= 0) {
         R[zrot][6*n+0] = 0.0;
         R[zrot][6*n+1] = 0.0;
         R[zrot][6*n+2] = 0.0;
         R[zrot][6*n+3] = 0.0;
         R[zrot][6*n+4] = 0.0;
         R[zrot][6*n+5] = 1.0;
       }
     }
     else {
#ifdef USE_EIGEN3
       Eigen::Matrix<double,6,6> Ri;
       Ri << 1, 0, 0,  0,  z, -y,
             0, 1, 0, -z,  0,  x,
             0, 0, 1,  y, -x,  0,
             0, 0, 0,  1,  0,  0,
             0, 0, 0,  0,  1,  0,
             0, 0, 0,  0,  0,  1;

       // transform Ri from basic to DOF_FRM coordinates
       NFrameData *nfd = geoSource->getNFrames();
       Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&nfd[nd.cd].frame[0][0]);

       Ri.topLeftCorner<3,3>() = T;
       Ri.topRightCorner<3,3>() = (T*Ri.topRightCorner<3,3>()).eval();
       Ri.bottomRightCorner<3,3>() = T;

       if(xloc >= 0) 
         for(int j=0; j<6; ++j) R[xloc][6*n+j] = Ri(0,j);
       if(yloc >= 0) 
         for(int j=0; j<6; ++j) R[yloc][6*n+j] = Ri(1,j);
       if(zloc >= 0) 
         for(int j=0; j<6; ++j) R[zloc][6*n+j] = Ri(2,j);
       if(xrot >= 0) 
         for(int j=0; j<6; ++j) R[xrot][6*n+j] = Ri(3,j);
       if(yrot >= 0) 
         for(int j=0; j<6; ++j) R[yrot][6*n+j] = Ri(4,j);
       if(zrot >= 0) 
         for(int j=0; j<6; ++j) R[zrot][6*n+j] = Ri(5,j);
#else
       std::cerr << "USE_EIGEN3 is not defined here in Rbm::computeRbms\n";
#endif
     }
   }
 }

 //  eliminate redundent rbms for 1D and 2D cases
#ifdef USE_EIGEN3
 typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd;
 Eigen::Map<RowMajorMatrixXd> R_(R.data(),R.numRow(),R.numCol()), A_(Amat->data(),R.numCol(),R.numCol());
 A_.setZero();
 for(n=0; n<comp->numComp; ++n) {
   Eigen::Block<Eigen::Map<RowMajorMatrixXd> > Rn = R_.block(0,6*n,R_.rows(),6);
   A_.block(6*n,6*n,6,6) = Rn.transpose()*Rn;
 }
#else
 R.transposeMult(R, *Amat);
#endif
 FullM &A = *Amat;
 int dimA = A.dim();
 double *diagA = (double *) dbg_alloca(sizeof(double)*dimA);

 for(i = 0; i < dimA; ++i) {
   diagA[i] = A[i][i];
   if(diagA[i] > 1e-12)
     diagA[i] = 1.0/sqrt(diagA[i]);
   else
     diagA[i] = 1.0;
 }
 for(i = 0; i < dimA; ++i)
   for(j = 0; j < dimA; ++j)
     A[i][j] *= diagA[i]*diagA[j];
  
 double max_value = A.maxAbs();
 if(debug) A.print("--- A = Rt*R matrix ---");

 // ... PERFORM SVD ON R^t*R 
 singularValueDecomposition(A, U, dimA, dimA, max_value, ngrbm, rank);

 FullM Rstar;
 if(rank == dimA) {
   Rstar = R;
 }
 else {
   FullM Ustar(U,rank,0,ncol,0);
   Rstar = R%Ustar;
 }
 if(debug) {
   R.print("--- R matrix ---");
   U.print("--- U matrix ---");
   Rstar.print("--- Rstar ---");
 }

 if(nrow == 0) { // FREE-FREE problem with no MPCs
   for(i=0; i<rank; ++i)
     for(j=0; j<numUncon; ++j)
       allRbm[0][i][j] = Rstar[j][i];
   nRbmPerComp[0] = rank;
 }
 else {
   ncol = rank;

   // ... EXTRACT Z MATRIX CONTAINING CONSTRAINED GRBM
   // ... AND CALCULATE MAXIMUM ENTRY PER COMPONENT
   FullM Z(nrow,ncol);  
   Z.zero();

   double cgmax = 0.0;

   if(numBC[0]) {
     for(inode=0; inode<dsa->numNodes(); ++inode) {

       position[0] = dsa->locate(inode, DofSet::Xdisp);
       position[1] = dsa->locate(inode, DofSet::Ydisp);
       position[2] = dsa->locate(inode, DofSet::Zdisp);
       position[3] = dsa->locate(inode, DofSet::Xrot);
       position[4] = dsa->locate(inode, DofSet::Yrot);
       position[5] = dsa->locate(inode, DofSet::Zrot);

       for(j=0; j<6; ++j) {
         if(position[j] >= 0) {
           cn = c_dsa->invRCN(position[j]);
           if(cn >= 0) {
             for(k=0; k<ncol; ++k) { 
               Z[cn][k] = Rstar[position[j]][k];
               if(std::abs(Z[cn][k]) > cgmax) cgmax = std::abs(Z[cn][k]);
             }
           }
         }
       }
     }
   }

   // ... DEAL WITH MPCs
   if(numMPC) { 
     for(i=0; i<numMPC; ++i) {
       for(j=0; j<mpc[i]->nterms; ++j) {
         int inode = mpc[i]->terms[j].nnum;
         int dof = dsa->locate(inode, (1 << mpc[i]->terms[j].dofnum));
         if(dof >= 0) {
           for(k=0; k<ncol; ++k) {
             Z[numBC[0]+i][k] += Rstar[dof][k]*mpc[i]->terms[j].coef.r_value;
             if(std::abs(Z[numBC[0]+i][k]) > cgmax) cgmax = std::abs(Z[numBC[0]+i][k]);
           }
         }
       }
     }
   }
   // ... PERFORM A SVD OPERATION ON Z:
   if(debug) 
     Z.print("----- Z matrix -----");
   U.setNewSize(ncol,ncol);
   singularValueDecomposition(Z, U, ncol, nrow, cgmax, ngrbm, rank);

   // ... CONSTRUCT NECESSARY VECTORS FOR Vstar:
   FullM Vstar(U,ngrbm,rank,ncol,0);

   // ... MATRIX-MATRIX MULTIPLY [GRBM] = [R][V]:
   FullM result = Rstar%Vstar;

   // ... CONFORM RIGID BODY MODES WITH RESPECT TO K MATRIX
   for(i=0; i<numdofs; ++i) {
     int cn = c_dsa->getRCN(i);
     if(cn >= 0 ) {
       for(j=0; j<ngrbm; ++j)
         allRbm[0][j][cn] = result[i][j];
     }
   }
   nRbmPerComp[0] = ngrbm;

   // ... DEBUGGING PRINT STATEMENTS ...
   if(debug) {
     U.print("----- U matrix       -----");
     Vstar.print("----- Vstar matrix   -----");
     R.print("----- R matrix       -----");
     result.print("----- Results matrix -----");
   }

 } 

 Vector v1(numUncon,0.0);

 ngrbm = 0;
 for(i=0; i<nComponents; ++i)
   ngrbm += nRbmPerComp[i];
 
 int addModes = ( (cornerModes) ? cornerModes->numCol() : 0);
 grbm = new Vector[ngrbm + addModes];

 for(i=0; i<ngrbm+addModes; ++i)
   grbm[i] = v1;
   
 // ... COPY COMPONENT RBMS TO ONE MATRIX GRBM
 int cnt = 0;
 for(n=0; n<nComponents; ++n)
   for(j=0; j<nRbmPerComp[n]; ++j) {
     for(k=0; k<numUncon; ++k)
       grbm[cnt][k] = allRbm[n][j][k];
     cnt += 1;
   }

 if(debug) {
   fprintf(stderr, "%d components rigid body modes were found.\n",nComponents);
   print();
  }

 for(i = 0; i < addModes; ++i)
   grbm[ngrbm+i][(*cornerModes)[2][i]] = 1.0;

 ngrbm += addModes;
}

void
Rbm::print()
{
   int j,k;
     for(k=0; k<numUncon; ++k) {
       for(j=0; j<numRBM(); ++j)
         //fprintf(stderr,"%f\t",allRbm[i][j][k]);
         fprintf(stderr," %8.4f ",grbm[j][k]);
       fprintf(stderr,"\n");
     }
}

void
Rbm::initialize(int numComponents, int *numRbmPerComp, int numdofs)
{
  Vector vec(numdofs, 0.0); // sets all elements of v to zero

  int i, n;
  for(n=0; n<numComponents; ++n) {
    for(i=0; i<numRbmPerComp[n]; ++i)
      allRbm[n][i] = vec;
  }
}

void
Rbm::singularValueDecomposition(FullM &A, FullM &Umat, int ncol, int nrow,
                                double max_value, int &ngrbm, int &rank)
{
   int i;
   int info = 0;
   int mindim = std::min(nrow,ncol);
   int maxdim = std::max(nrow,ncol);
   double *w    = (double*) dbg_alloca(sizeof(double)* maxdim);
   double *e    = (double*) dbg_alloca(sizeof(double)* maxdim);
   double *work = (double*) dbg_alloca(sizeof(double)* maxdim);
#ifdef FILERING
  double maxA= A.maxAbs();
  for(int i=0; i<A.numCol()*A.numRow(); i++) //HB
    if(fabs((A.data())[i])<1.E-16*maxA+1.E-20) (A.data())[i] = 0.0;
#endif

   for(int i=0; i<maxdim; i++) { w[i] = 0.0; e[i] = 0.0; work[i] = 0.0; }

   _FORTRAN(dsvdc)(A.data(), ncol, ncol, nrow, w, e,
                   Umat.data(), ncol, Umat.data(), ncol, work, 10, info);

//  ... CHECK RETURN STATUS OF DSVDC:

   if(info <  0)
     fprintf(stderr," *** ERROR: Illegal value in argument #%d \n",info);

   if(info >  0) { 
     fprintf(stderr," *** WARNING: %d diagonals did not converge in Rbm::singularValueDecomposition() \n",info);
     // fprintf(stderr,"%e %e %e %e %e %e\n",w[0],w[1],w[2],w[3],w[4],w[5]);
     // fprintf(stderr," *** ERROR:\n");
   }

//  ... DETERMINE RANK

   rank = 0;
   double tolerance = max_value*tolgrb;
   for(i=0; i<mindim; ++i) 
     if(std::abs(w[i]) > tolerance) rank++;

   ngrbm = ncol - rank;
}

void
Rbm::getRBMs(Vector* rigidBodyModes)
{
 int i;
 for(i=0; i<ngrbm; ++i)
   rigidBodyModes[i] = grbm[i];
}

void
Rbm::getRBMs(VectorSet& rigidBodyModes)
{
 int i;
 for(i=0; i<ngrbm; ++i)
    rigidBodyModes[i] = grbm[i];
}

void
Rbm::reBuildGeometricRbms(GeomState *gs)
{
  int oldNumRbm = ngrbm;

  int rank  = 0;

  int i,inode,j,k,cn, ncol;

  int position[6];

  int numdofs = dsa->size();

  FullM &R = *Rmat;

  U.zero();

 int n;
 for(n=0; n<nComponents; ++n) {

   R.zero();

   ncol     = nRbmPerComp[n];
   int nrow = numBC[n];

/*
   NodeState &ns1 = (*gs)[comp->order[comp->xcomp[n]]];

   x0 = ns1.x;
   y0 = ns1.y;
   z0 = ns1.z;
*/

   double xc = 0.0, yc = 0.0, zc = 0.0;
   for(i = comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
     inode = comp->order[i];
     if(dsa->firstdof(inode) == -1) continue;

     NodeState &nd = (*gs)[inode];

     xc += nd.x;
     yc += nd.y;
     zc += nd.z;
   }
   double x0 = xc/(comp->xcomp[n+1]-comp->xcomp[n]);
   double y0 = yc/(comp->xcomp[n+1]-comp->xcomp[n]);
   double z0 = zc/(comp->xcomp[n+1]-comp->xcomp[n]);
   		  
// ... CALCULATE R MATRIX CONTAINING GEOMETRIC RBM

   for(i = comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
   
     inode = comp->order[i];

     if(dsa->firstdof(inode) == -1) continue;

     NodeState &ns = (*gs)[inode];
     
     double x = ns.x - x0;
     double y = ns.y - y0;
     double z = ns.z - z0;

     int xloc = dsa->locate( inode, DofSet::Xdisp);
     int yloc = dsa->locate( inode, DofSet::Ydisp);
     int zloc = dsa->locate( inode, DofSet::Zdisp);
     int xrot = dsa->locate( inode, DofSet::Xrot);
     int yrot = dsa->locate( inode, DofSet::Yrot);
     int zrot = dsa->locate( inode, DofSet::Zrot);

     if(xloc >= 0) {
         R[xloc][0] = 1.0;
         R[xloc][1] = 0.0;
         R[xloc][2] = 0.0;
         R[xloc][3] = 0.0;
         R[xloc][4] =   z;
         R[xloc][5] =  -y;
     }
     if(yloc >= 0) {
         R[yloc][0] = 0.0;
         R[yloc][1] = 1.0;
         R[yloc][2] = 0.0;
         R[yloc][3] =  -z;
         R[yloc][4] = 0.0;
         R[yloc][5] =   x;
     }
     if(zloc >= 0) {
         R[zloc][0] = 0.0;
         R[zloc][1] = 0.0;
         R[zloc][2] = 1.0;
         R[zloc][3] =   y;
         R[zloc][4] =  -x;
         R[zloc][5] = 0.0;
     }
     if(xrot >= 0) {
         R[xrot][0] = 0.0;
         R[xrot][1] = 0.0;
         R[xrot][2] = 0.0;
         R[xrot][3] = 1.0;
         R[xrot][4] = 0.0;
         R[xrot][5] = 0.0;
     }
     if(yrot >= 0) {
         R[yrot][0] = 0.0;
         R[yrot][1] = 0.0;
         R[yrot][2] = 0.0;
         R[yrot][3] = 0.0;
         R[yrot][4] = 1.0;
         R[yrot][5] = 0.0;
      }
      if(zrot >= 0) {
         R[zrot][0] = 0.0;
         R[zrot][1] = 0.0;
         R[zrot][2] = 0.0;
         R[zrot][3] = 0.0;
         R[zrot][4] = 0.0;
         R[zrot][5] = 1.0;
      }
   }

// ... FOR A FREE-FREE PROBLEM, RETURN GEOMETRIC RIGID BODY MODES.

    R.transposeMult(R,*Amat);
    FullM &A = *Amat;

    int dimA = A.dim();

    double *diagA = (double *) dbg_alloca(sizeof(double)*dimA);

    for(i = 0; i < dimA; ++i) {
      diagA[i] = A[i][i];
      if(diagA[i] > 1e-12)
        diagA[i] = 1.0/sqrt(diagA[i]);
      else
        diagA[i] = 1.0;
    }
    for(i = 0; i < dimA; ++i)
      for(j = 0; j < dimA; ++j)
        A[i][j] *= diagA[i]*diagA[j];

    double max_value = A.maxAbs();  // max() finds maximum value of matrix A

// ... PERFORM SVD ON R^t*R

    singularValueDecomposition(A, U, dimA, dimA, max_value, ngrbm, rank);

    if(rank == dimA && numBC[n] == 0) {
      for(i=0; i<ncol; ++i)
        for(j=0; j<numUncon; ++j)
          allRbm[n][i][j] = R[j][i];
        nRbmPerComp[n] = ncol;
      continue;
    }

// ... COPY NECESSARY ROWS OF U TO USTAR:

    FullM Ustar(U,rank,0,ncol,0);

    FullM Rstar;
    if(rank == 6)
      Rstar = R;
    else
      Rstar = R%Ustar;

    if(numBC[n] == 0) {
      for(i=0; i<rank; ++i)
        for(j=0; j<numUncon; ++j)
          allRbm[n][i][j] = Rstar[j][i];
      nRbmPerComp[n] = rank;
      continue;
    }

    ncol = rank;

// ... EXTRACT Z MATRIX CONTAINING CONSTRAINED GRBM
// ... AND CALCULATE MAXIMUM ENTRY PER COMPONENT

    FullM Z(nrow,ncol);

    double cgmax = 0.0;

    int offset = 0;
    int ii;
    for(ii=0; ii<n; ++ii) {
      offset += numBC[ii];
    }

    for(i=comp->xcomp[n]; i<comp->xcomp[n+1]; ++i) {
      inode = comp->order[i];

      position[0] = dsa->locate( inode, DofSet::Xdisp);
      position[1] = dsa->locate( inode, DofSet::Ydisp);
      position[2] = dsa->locate( inode, DofSet::Zdisp);
      position[3] = dsa->locate( inode, DofSet::Xrot);
      position[4] = dsa->locate( inode, DofSet::Yrot);
      position[5] = dsa->locate( inode, DofSet::Zrot);

      for(j=0; j<6; ++j) {
        if(position[j] >= 0) {
          if((cn = c_dsa->invRCN(position[j])) >= 0) {
            int k;
            for(k=0; k<ncol; ++k) {
              Z[cn-offset][k] = Rstar[position[j]][k];
              if(std::abs(Z[cn-offset][k]) > cgmax) cgmax = std::abs(Z[cn-offset][k]);
            }
          }
        }
      }
    }

// ... PERFORM A SVD OPERATION ON Z:

   U.setNewSize(ncol,ncol);

   singularValueDecomposition(Z, U, ncol, nrow, cgmax, ngrbm, rank);

// ... CONSTRUCT NECESSARY VECTORS FOR Vstar:

   FullM Vstar(U,ngrbm,rank,ncol,0);

// ... MATRIX-MATRIX MULTIPLY [GRBM] = [R][V]:

   FullM result = Rstar%Vstar;

// ... CONFORM RIGID BODY MODES WITH RESPECT TO K MATRIX

    for(i=0; i<numdofs; ++i) {
      int cn = c_dsa->getRCN(i);
      if(cn >= 0 ) {
        for(j=0; j<ngrbm; ++j)
          allRbm[n][j][cn] = result[i][j];
      }
    }
    nRbmPerComp[n] = ngrbm;

 } // ... END OF COMPONENT LOOP

 // Count the number of geometric rigid body modes
 ngrbm = 0;
 for(i=0; i<nComponents; ++i)
   ngrbm += nRbmPerComp[i];

 // MODIFICATION 1: only allocate space if number of rbm changes
 // from one newton iteration to the next
 if(ngrbm != oldNumRbm) {
   grbm = new Vector[ngrbm];
   Vector v1(numUncon, 0.0);
   for(i=0; i<ngrbm; ++i)
     grbm[i] = v1;
 }
  
// ... COPY COMPONENT RBMS TO ONE MATRIX GRBM
 int cnt = 0;
 for(n=0; n<nComponents; ++n)
   for(j=0; j<nRbmPerComp[n]; ++j) {
     for(k=0; k<numUncon; ++k)
       grbm[cnt][k] = allRbm[n][j][k];
     cnt += 1;
   }

 }

void
Rbm::getRBMs(double *rigidBodyModes, int ndof, int *dof, int numGRBM, int offset)
{
 int i,iDof;
 if(numGRBM < 0) numGRBM = ngrbm;
 for(i=0; i<numGRBM; ++i)
   for(iDof=0; iDof<ndof; ++iDof)
     if(dof[iDof] >= 0)
       rigidBodyModes[iDof+i*ndof] = grbm[i+offset][dof[iDof]];
     else
       rigidBodyModes[iDof+i*ndof] = 0.0;
}

void
Rbm::getRBMs(DComplex *rigidBodyModes, int ndof, int *dof, int numGRBM, int offset)
{
 int i,iDof;
 if(numGRBM < 0) numGRBM = ngrbm;
 for(i=0; i<numGRBM; ++i)
   for(iDof=0; iDof<ndof; ++iDof)
     if(dof[iDof] >= 0)
       rigidBodyModes[iDof+i*ndof] = cgrbm[i+offset][dof[iDof]];
     else
       rigidBodyModes[iDof+i*ndof] = 0.0;
}

void
Rbm::getScaledRBMs(double *rigidBodyModes, int ndof, int *dof, double *scaling, int numGRBM, int offset)
{
 int i,iDof;
 if(numGRBM < 0) numGRBM = ngrbm;
 for(i=0; i<numGRBM; ++i)
   for(iDof=0; iDof<ndof; ++iDof)
     if(dof[iDof] >= 0)
       rigidBodyModes[iDof+i*ndof] = grbm[i+offset][dof[iDof]] * scaling[iDof];
     else
       rigidBodyModes[iDof+i*ndof] = 0.0;
}

void
Rbm::clean_up()
{
 if(Rmat) {
  Rmat->clean_up();
  Rmat=0;
 }

 int i;
 if(allRbm) 
   for(i=0; i<nComponents; ++i)
     if(allRbm[i]) {
       allRbm[i]->clean_up();
       allRbm[i]=0;
     }
}

Rbm::~Rbm()
{
 if(grbm && myMemory) { delete [] grbm; grbm=0; }
 if(cgrbm && myMemory) { delete [] cgrbm; cgrbm=0; }
 
 if(nRbmPerComp) { delete [] nRbmPerComp; nRbmPerComp=0; }
 if(firstDofOfComp) { delete [] firstDofOfComp; firstDofOfComp=0; }
 if(numDofPerComp) { delete [] numDofPerComp; numDofPerComp=0; }
 if(numBC) { delete [] numBC; numBC=0; }

 if(xyzRot) {
   for(int i=0; i<nComponents; ++i)
     if(xyzRot[i]) { delete [] xyzRot[i]; xyzRot[i] = 0; }
   delete [] xyzRot; xyzRot = 0;
 }
 if(Zstar) { delete Zstar; Zstar = 0; }
 if(Rc) { delete Rc; Rc = 0; }
 if(Zmpc) { delete Zmpc; Zmpc = 0; }
 if(Rmat) { delete Rmat; Rmat=0; }
 if(Amat) { delete Amat; Amat=0; }
 if(allRbm) {
   for(int i=0; i<nComponents; ++i)
     if(allRbm[i]) { delete [] allRbm[i]; allRbm[i] = 0; }
   delete [] allRbm; allRbm = 0;
 }
}

template<>
void
Rbm::getRBMs(double *rigidBodyModes, bool transpose)
{
  int i, iDof;
  if(transpose) {
    for(i=0; i<ngrbm; ++i)
      for(iDof=0; iDof<numUncon; ++iDof)
         rigidBodyModes[i+iDof*ngrbm] = grbm[i][iDof];
  } 
  else {
    for(i=0; i<ngrbm; ++i)
      for(iDof=0; iDof<numUncon; ++iDof)
         rigidBodyModes[iDof+i*numUncon] = grbm[i][iDof];
  }
}

template<>
void
Rbm::getRBMs(DComplex *rigidBodyModes, bool transpose)
{
  int i, iDof;
  if(transpose) {
    for(i=0; i<ngrbm; ++i)
      for(iDof=0; iDof<numUncon; ++iDof)
         rigidBodyModes[i+iDof*ngrbm] = cgrbm[i][iDof];
  } 
  else {
    for(i=0; i<ngrbm; ++i)
      for(iDof=0; iDof<numUncon; ++iDof)
        rigidBodyModes[iDof+i*numUncon] = cgrbm[i][iDof];
  }
}

void
Rbm::getRBMs(double *rigidBodyModes, std::set<int> &rbmFilters)
{
  int i = 0;
  for(std::set<int>::iterator it = rbmFilters.begin(); it != rbmFilters.end(); ++it) {
    int iMode = *it;
    if(iMode < ngrbm) {
      for(int iDof = 0; iDof < numUncon; ++iDof)
        rigidBodyModes[iDof+i*numUncon] = grbm[iMode][iDof];
      i++;
    }
  }
}
