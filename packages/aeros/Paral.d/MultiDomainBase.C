#include <Paral.d/MultiDomainBase.h>
#include <Math.d/FullMatrix.h>
#include <Paral.d/SubDOp.h>
#include <Solvers.d/MultiDomainRbm.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/DistHelper.h>

MultiDomainBase::MultiDomainBase(SolverInfo &_sinfo)
 : sinfo(_sinfo), X(NULL), R(NULL), numR(0)
{
}

MultiDomainBase::~MultiDomainBase()
{
  if(R) delete R;
  if(X) delete X;
}

void
MultiDomainBase::projector_prep(MultiDomainRbm<double> *rbms, GenSubDOp<double> *M)
{
  if(R) return; // already done it

  numR = rbms->numRBM(); 
  if (!numR) return;

  if(sinfo.isDynam() && !sinfo.rbmFilters.empty()) {
    int count = 0;
    for(std::set<int>::iterator it = sinfo.rbmFilters.begin(); it != sinfo.rbmFilters.end(); ++it) {
      if(*it < numR) count++;
      else filePrint(stderr," *** WARNING: mode %d specified under RBMFILTER does not exist.\n", *it+1);
    }
    numR = count;
  }

  R = new GenDistrVectorSet<double>(numR, rbms->solVecInfo());

  if(sinfo.filterFlags) {
    filePrint(stderr," ... Building the RBM Projector     ...\n");
    if(sinfo.isDynam() && !sinfo.rbmFilters.empty()) {
      filePrint(stderr," ... Number of Filtered Modes = %-4d...\n",numR);   
      rbms->getRBMs(*R, sinfo.rbmFilters);
    }
    else {
      filePrint(stderr," ... Number of RBMs = %-4d          ...\n",numR);
      rbms->getRBMs(*R);
    }
  }

  double *y = (double *) dbg_alloca(numR*sizeof(double));
  double *x = (double *) dbg_alloca(numR*sizeof(double));

  X = new DistrVectorSet(numR, rbms->solVecInfo());

  if(sinfo.isDynam() || sinfo.filterQ == 0) {
    DistrVectorSet MR(numR, rbms->solVecInfo());
    for(int i=0; i<numR; ++i) {
      M->mult((*R)[i], MR[i]);
    }

    FSFullMatrix RtMR(numR,numR);
    for(int i=0; i<numR; ++i)
      for(int j=i; j<numR; ++j)
        RtMR[i][j] = RtMR[j][i] = (*R)[i]*MR[j];

    FSFullMatrix RtMRinverse(numR, numR);
    RtMRinverse = RtMR.invert();

    for(int i=0; i<MR.size()/numR; ++i) {
      for(int j=0; j<numR; ++j) y[j] = MR[j].data()[i];
      RtMRinverse.mult(y, x);
      for(int j=0; j<numR; ++j) (*X)[j].data()[i] = x[j];
    }
  }
  else {
    FSFullMatrix RtR(numR,numR);
    for(int i=0; i<numR; ++i)
      for(int j=i; j<numR; ++j)
        RtR[i][j] = RtR[j][i] = (*R)[i]*(*R)[j];

    FSFullMatrix RtRinverse(numR, numR);
    RtRinverse = RtR.invert();

    for(int i=0; i<R->size()/numR; ++i) {
      for(int j=0; j<numR; ++j) y[j] = (*R)[j].data()[i];
      RtRinverse.mult(y, x);
      for(int j=0; j<numR; ++j) (*X)[j].data()[i] = x[j];
    }
  }
}

void
MultiDomainBase::eigmode_projector_prep(DistrInfo &solVecInfo)
{
  if(R) return; // already done it or requested some other filter

  // Read computed eigenvectors from file EIGENMODES
  // ======================================
#ifdef DISTRIBUTED
  char *filename = new char[40];
  sprintf(filename,"EIGENMODES%d",structCom->myID());
  BinFileHandler modefile(filename, "r");
#else
  BinFileHandler modefile("EIGENMODES", "r");
#endif
  if(modefile.get_fileid() <= 0) { fprintf(stderr, " *** Error: Failed to open EIGENMODES file ***\n"); exit(-1); }

  int numR;
  modefile.read(&numR, 1);
  filePrint(stderr," ... Reading %d modes from EIGENMODES file ...\n", numR);

  int eigsize;
  modefile.read(&eigsize, 1);
  if(eigsize != solVecInfo.totLen()) {
    fprintf(stderr, " *** Error: Bad data in EIGENMODES file %d %d ***\n", eigsize, solVecInfo.totLen()); exit(-1);
  }

  R = new DistrVectorSet(numR, solVecInfo);
  double *data = new double[eigsize];
  for(int i = 0; i < numR; ++i) {
    modefile.read(data, eigsize);
    for(int j=0; j<eigsize; ++j) (*R)[i].data()[j] = data[j];
  }
  delete [] data;

  FSFullMatrix RtR(numR,numR);
  for(int i=0; i<numR; ++i)
    for(int j=i; j<numR; ++j)
      RtR[i][j] = RtR[j][i] = (*R)[i]*(*R)[j];

  FSFullMatrix RtRinverse(numR, numR);
  RtRinverse = RtR.invert();

  double *y = (double *) dbg_alloca(numR*sizeof(double));
  double *x = (double *) dbg_alloca(numR*sizeof(double));

  X = new DistrVectorSet(numR, solVecInfo);
  for(int i=0; i<R->size()/numR; ++i) {
    for(int j=0; j<numR; ++j) y[j] = (*R)[j].data()[i];
    RtRinverse.mult(y, x);
    for(int j=0; j<numR; ++j) (*X)[j].data()[i] = x[j];
  }
}

void
MultiDomainBase::projector_prep(MultiDomainRbm<std::complex<double> > *, GenSubDOp<std::complex<double> > *)
{
}

void
MultiDomainBase::trProject(DistrVector &b)
{
  int numR = (R) ? R->numVec() : 0;
  if(numR == 0) return;

  double *y = (double *) dbg_alloca(numR*sizeof(double));

  // y = Rt*b
  for(int i=0; i<numR; ++i)
    y[i] = (*R)[i]*b;

  // b = b - X*y
  for(int i=0; i<b.size(); ++i)
    for(int j=0; j<numR; ++j) b.data()[i] -= (*X)[j].data()[i]*y[j];
}

void
MultiDomainBase::trProject(ComplexDistrVector &)
{
}

void
MultiDomainBase::project(DistrVector &v)
{
  int numR = (R) ? R->numVec() : 0;
  if(numR == 0) return;

  double *y = (double *) dbg_alloca(numR*sizeof(double));

  // y = Xt*v
  for(int i=0; i<numR; ++i)
    y[i] = (*X)[i]*v;

  // v = v - R*y
  for(int i=0; i<v.size(); ++i)
    for(int j=0; j<numR; ++j) v.data()[i] -= (*R)[j].data()[i]*y[j];
}

void
MultiDomainBase::project(ComplexDistrVector &)
{
}

