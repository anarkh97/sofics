#ifdef USE_SCALAPACK
#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Timers.d/GetTime.h>
#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Plh.h"

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
splh(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool project, double &dtime,
       int npMax, int scpkMB, int scpkNB, int scpkMP, int scpkNP, int option) {
  // Setup
  int mypid, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Instantiate the solver
  Plh solver(A);

  // Set Parameters prior to initializing
  solver.setRtol(reltol);
  solver.setMaxIterRatio(maxite);
  if (scaling) {
    if (mypid == 0) {
      std::cout << "Using scaled matrix" << std::endl;    
    }
    solver.setColumnScaling();
  }

  if(option) // check if running Polytope Faces Pursuit
    solver.setPolytopeFacesPursuit();

  if (npMax > 0) {
    if (mypid == 0) {
      std::cout << "Reduced mesh size limited to " << npMax << std::endl;    
    }
    solver.setMaxNP(npMax);
  }
  if (scpkMB > 0 && scpkNB > 0) {
    if (mypid == 0) {
      std::cout << "ScaLAPACK block sizes are (" << scpkMB << "," << scpkNB << ")" << std::endl;    
    }
    solver.setABlockSize(scpkMB, scpkNB);
  }
  if (scpkMP > 0 && scpkNP > 0 && scpkMP*scpkNP <= nprocs) {
    if (mypid == 0) {
      std::cout << "ScaLAPACK processor grid is (" << scpkMP << "," << scpkNP << ")" << std::endl;    
    }
    solver.setAProcGrid(scpkMP, scpkNP);
  }

  // Loads the matrix and RHS and print a summary. Can't change some parameters after this.
  solver.init(A, b);
  solver.summary();

  // Solve
  int nfree = solver.solve();

  // Output
  solver.printTimes();
  info = solver.getStatus(); // XXX should be set to 2 if problem size is illegal and 3 if solver did not converge
  if (mypid == 0) {
    std::cout << "splh info = " << info << std::endl;
  }
  rnorm = solver.getResidualNorm();
  dtime = solver.getDownDateTime();
  return solver.getSolution();
}

#endif
#endif
