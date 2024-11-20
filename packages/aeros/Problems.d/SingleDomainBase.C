#include <Problems.d/SingleDomainBase.h>
#include <Math.d/FullMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/DistHelper.h>

SingleDomainBase::SingleDomainBase(SolverInfo &_sinfo)
 : sinfo(_sinfo), X(NULL), Rmem(NULL), numR(0)
{
}

SingleDomainBase::~SingleDomainBase()
{
  if(Rmem) delete [] Rmem;
  if(X) delete X;
}

void
SingleDomainBase::projector_prep(Rbm *rbms, SparseMatrix *M)
{
  if(Rmem) return; // already done it

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
  // KHP: store this pointer to the RBMs to use in the actual
  //      projection step within the time loop.
  int ndof = rbms->numDof();
  Rmem = new double[numR*ndof];

  if(sinfo.filterFlags) {
    filePrint(stderr," ... Building the RBM Projector     ...\n");
    if(sinfo.isDynam() && !sinfo.rbmFilters.empty()) {
      filePrint(stderr," ... Number of Filtered Modes = %-4d...\n",numR);   
      rbms->getRBMs(Rmem, sinfo.rbmFilters);
    }
    else {
      filePrint(stderr," ... Number of RBMs = %-4d          ...\n",numR);
      rbms->getRBMs(Rmem);
    }
  }
  else if(sinfo.hzemFilterFlag) {
    filePrint(stderr," ... Building the HZEM Projector    ...\n");
    filePrint(stderr," ... Number of HZEMs = %-4d         ...\n",numR);
    for(int n=0; n<ndof; ++n) Rmem[n] = 1.;
  }
  else if(sinfo.slzemFilterFlag) {
    filePrint(stderr," ... Building the SLZEM Projector    ...\n");
    filePrint(stderr," ... Number of SLZEMs = %-4d         ...\n",numR);
    for(int n=0; n<ndof; ++n) Rmem[n] = 1.;
  }

  StackFSFullMatrix Rt(numR, ndof, Rmem);
  X = new FSFullMatrix(ndof,numR);

  if(sinfo.isDynam() || sinfo.filterQ == 0) {
    double *MRmem = new double[numR*ndof];
    StackFSFullMatrix MRt(numR, ndof, MRmem);

    for(int n=0; n<numR; ++n)
      M->mult(Rmem+n*ndof, MRmem+n*ndof);

    FSFullMatrix MR = MRt.transpose();

    FSFullMatrix RtMR(numR,numR);
    Rt.mult(MR,RtMR);

    FSFullMatrix RtMRinverse = RtMR.invert();

    MR.mult(RtMRinverse,(*X));

    delete [] MRmem;
  }
  else {
    FSFullMatrix R = Rt.transpose();

    FSFullMatrix RtR(numR,numR);
    Rt.mult(R,RtR);

    FSFullMatrix RtRinverse = RtR.invert();

    R.mult(RtRinverse,(*X));
  }
}

void
SingleDomainBase::projector_prep(Rbm *, ComplexSparseMatrix *)
{
}

void
SingleDomainBase::eigmode_projector_prep()
{
  if(Rmem) return; // already done it or requested some other filter

  // Read computed eigenvectors from file EIGENMODES
  // ======================================
  BinFileHandler modefile("EIGENMODES" ,"r");
  if(modefile.get_fileid() <= 0) { fprintf(stderr, " *** Error: Failed to open EIGENMODES file ***\n"); exit(-1); }

  modefile.read(&numR, 1);
  fprintf(stderr," ... Reading %d modes from EIGENMODES file ...\n", numR);

  int eigsize;
  modefile.read(&eigsize, 1);

  Rmem = new double[numR*eigsize];
  for(int i = 0; i < numR; ++i)
    modefile.read(Rmem+i*eigsize, eigsize);

/*
  // Check if eigenvectors are M-orthonormal: Phi_i*M*Phi_i = 1 and Phi_i*M*Phi_j = 0
  // ======================================
  GenSparseMatrix<double> *M = (GenSparseMatrix<double> *) allOps.M; // XXX won't work for complex
  double *tPhiM =  new double[numR*eigsize];
  for(int i = 0; i < numR; ++i)
    M->mult(Rmem+i*eigsize, tPhiM+i*eigsize);  // taking advantage of symmetry of M and computing
                                               // M*Phi_i instead of transpose(Phi_i)*M
  for(int i = 0; i < numR; ++i) {
    for(int j = 0; j < numR; ++j) {
      double PhiMPhi = 0;
      for(int k = 0; k < eigsize; ++k) {
        PhiMPhi += Rmem[k+i*eigsize]*tPhiM[k+j*eigsize];
      }
      fprintf(stderr, "Phi_%d*M*Phi_%d = %19.11e\n", i, j, PhiMPhi);
    }
    fprintf(stderr,"\n");
  }
*/
  // Build U_c(U_c^T*U_c)^{-1} for projector
  // ======================================
  StackFSFullMatrix Rt(numR, eigsize, Rmem);
  FSFullMatrix R = Rt.transpose();
  FSFullMatrix RtR(numR, numR);
  Rt.mult(R, RtR);
  FSFullMatrix RtRinverse = RtR.invert();

  X = new FSFullMatrix(eigsize, numR); // X = U(U^t*U)^{-1}
  R.mult(RtRinverse,(*X));
}

void
SingleDomainBase::trProject(Vector &f)
{
  if (!numR) return;

  int ndof = f.size();

  double *yMem = (double *) dbg_alloca(numR*sizeof(double));
  double *zMem = new double[ndof];

  StackVector y(numR,yMem);
  StackVector z(ndof,zMem);

  StackFSFullMatrix Rt(numR, ndof, Rmem);

  // y = Rt*f
  Rt.mult(f,y);

  // z = X*y
  X->mult(y,z);

  // f = f - z;
  f -= z;

  delete [] zMem;
}

void
SingleDomainBase::trProject(ComplexVector &)
{
}

void
SingleDomainBase::project(Vector &v)
{
  if (!numR) return;

  int ndof = v.size();

  double *yMem = (double *) dbg_alloca(numR*sizeof(double));
  double *zMem = new double[ndof];

  StackVector y(numR,yMem);
  StackVector z(ndof,zMem);

  StackFSFullMatrix Rt(numR, ndof, Rmem);

  // y = Xt*v
  X->trMult(v,y);

  // z = R*y
  Rt.trMult(y,z);

  // v = v - z;
  v -= z;

  delete [] zMem;
}

void
SingleDomainBase::project(ComplexVector &)
{
}

