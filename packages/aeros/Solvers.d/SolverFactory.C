#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/Skyline.d/DistSky.h>
#include <Math.d/Skyline.d/BlockSky.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/DistBLKSparse.h>
#include <Math.d/EiSparseMatrix.h>
#include <Solvers.d/Spooles.h>
#include <Solvers.d/Mumps.h>
#include <Solvers.d/BCGSolver.h>
#include <Solvers.d/CRSolver.h>
#include <Solvers.d/GmresSolver.h>
#include <Solvers.d/PCGSolver.h>
#include <Solvers.d/GoldfarbIdnani.h>
#include <Rom.d/GalerkinProjectionSolver.h>
#include <Rom.d/EiGalerkinProjectionSolver.h>
#include <Math.d/DiagMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Driver.d/Communicator.h>
#include <Utils.d/DistHelper.h>
#include "CholmodSolver.h"
#include "PardisoSolver.h"
#include <Solvers.d/SolverFactory.h>

template<class Scalar>
GenSolver<Scalar> *
GenSolverFactory<Scalar>::createSolver(const Connectivity *con, const EqNumberer *eqnum, const SolverCntl &cntl,
                                       GenSparseMatrix<Scalar> *&sparse,
                                       int ngrbm, FSCommunicator *com, std::string name) const
{
	//std::cerr << "1. creating solver: type = " << cntl.type << ", subtype = " << cntl.subtype << ", name = " << name << std::endl;
	GenSolver<Scalar> *solver = nullptr;
	switch(cntl.type) {
		case SolverSelection::Direct : { // direct solver

			switch(cntl.subtype) {
				default:
				case 0: {
					GenSkyMatrix<Scalar> *s = new GenSkyMatrix<Scalar>(con, eqnum, cntl.trbm, cntl.scaled);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
				case 1: {
					GenBLKSparseMatrix<Scalar> *s = new GenBLKSparseMatrix<Scalar>(con, eqnum, cntl.trbm, cntl, ngrbm);
					s->zeroAll();
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
				case 2: {
					GenBlockSky<Scalar> *s = new GenBlockSky<Scalar>(con, eqnum, cntl.trbm);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#ifdef USE_EIGEN3
				case 3: {
					GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *s =
							new GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(con, eqnum, true);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
				case 4: {
					GenEiSparseMatrix<Scalar,Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *s =
							new GenEiSparseMatrix<Scalar,Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(con, eqnum, true);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#ifdef EIGEN_CHOLMOD_SUPPORT
				case 5: {
					GenEiSparseMatrix<Scalar,Eigen::CholmodDecomposition<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *s =
							new GenEiSparseMatrix<Scalar,Eigen::CholmodDecomposition<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(con, eqnum, true);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
#ifdef EIGEN_UMFPACK_SUPPORT
				case 6: {
					GenEiSparseMatrix<Scalar,Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar> > > *s =
							new GenEiSparseMatrix<Scalar,Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar> > >(con, eqnum, false);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
				case 7: {
          GenEiSparseMatrix<Scalar,Eigen::SuperLU<Eigen::SparseMatrix<Scalar> > > *s =
            new GenEiSparseMatrix<Scalar,Eigen::SuperLU<Eigen::SparseMatrix<Scalar> > >(con, eqnum, false);
          solver = (GenSolver<Scalar> *) s;
          sparse = (GenSparseMatrix<Scalar> *) s;
        } break;
#endif
#endif
#ifdef USE_SPOOLES
				case 8: {
          GenSpoolesSolver<Scalar> *s = new GenSpoolesSolver<Scalar>(con, eqnum, cntl);
          solver = (GenSolver<Scalar> *) s;
          sparse = (GenSparseMatrix<Scalar> *) s;
        } break;
#endif
#ifdef USE_MUMPS
				case 9: {
          GenMumpsSolver<Scalar> *s = new GenMumpsSolver<Scalar>(con, eqnum, cntl, (int *)0, com);
          solver = (GenSolver<Scalar> *) s;
          sparse = (GenSparseMatrix<Scalar> *) s;
        }
        break;
#endif
#ifdef USE_EIGEN3
#ifdef EIGEN_SPARSELU_SUPPORT
				case 15: {
          GenEiSparseMatrix<Scalar,Eigen::SparseLU<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > > *s =
            new GenEiSparseMatrix<Scalar,Eigen::SparseLU<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > >(con, eqnum, false);
          solver = (GenSolver<Scalar> *) s;
          sparse = (GenSparseMatrix<Scalar> *) s;
        } break;
#endif
#ifdef EIGEN_SPQR_SUPPORT
				case 16: {
					GenEiSparseMatrix<Scalar,Eigen::SPQR<Eigen::SparseMatrix<Scalar> > > *s =
							new GenEiSparseMatrix<Scalar,Eigen::SPQR<Eigen::SparseMatrix<Scalar> > >(con, eqnum, false);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
				case 17: {
					auto *s =
						new GenEiSparseMatrix<Scalar,Eigen::SparseQR<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > >(con, eqnum, false);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
				case 18: { // CholMod
					auto p = getCholmod<Scalar>(con, eqnum);
					solver = p.first;
					sparse = p.second;
				}
				break;
				case 19: { // Pardiso in MKL
					auto p = getPardiso<Scalar>(con, eqnum);
					solver = p.first;
					sparse = p.second;
				}
				break;
			}

		} break; // case Direct solver

		default: break;
	}
	//solver->mpicomm = com; // XXXX
	return solver;
}

template<class Scalar>
GenSolver<Scalar> *
GenSolverFactory<Scalar>::createDistSolver(const Connectivity *con, const EqNumberer *eqnum, const SolverCntl &cntl,
                                           GenSparseMatrix<Scalar> *&sparse,
                                           FSCommunicator *com, std::string name) const
{
  //std::cerr << "2. creating solver: type = " << cntl.type << ", subtype = " << cntl.subtype << ", name = " << name << std::endl;
  int neq       = eqnum->size();
  int myCPU = (com) ? com->cpuNum() : 0;
  int numCPUs = (com) ? com->size() : 1;
  int neqPerCPU = neq/numCPUs;
  int remainder = neq%numCPUs;
  int firstAlpha = myCPU * neqPerCPU + ((myCPU < remainder) ? myCPU : remainder);
  int nRowAlpha  = neqPerCPU + ((myCPU < remainder) ? 1 : 0);
  GenSolver<Scalar> *solver = 0;
  switch(cntl.subtype) {
    case 0: {
      GenSkyMatrix<Scalar> *s = new GenDistSky<Scalar>(con, eqnum, cntl.trbm, firstAlpha, nRowAlpha/*, com*/); // TODO: pass com
      solver = (GenSolver<Scalar> *) s;
      sparse = (GenSparseMatrix<Scalar> *) s;
    } break;
    case 1: {
      GenBLKSparseMatrix<Scalar> *s = new GenDistBLKSparse<Scalar>(con, eqnum, cntl.trbm, cntl, firstAlpha, nRowAlpha/*, com*/); // TODO: pass com
      s->zeroAll();
      solver = (GenSolver<Scalar> *) s;
      sparse = (GenSparseMatrix<Scalar> *) s;
    } break;
  }
  //solver->mpicomm = com; // XXXX
  return solver;
}

template<class Scalar>
GenSolver<Scalar> *
GenSolverFactory<Scalar>::createSolver(const Connectivity *con, const DofSetArray *dsa, const ConstrainedDSA *cdsa,
                                       const SolverCntl &cntl, GenSparseMatrix<Scalar> *&sparse,
                                       Rbm *rbm, GenSparseMatrix<Scalar> *&spp, GenSolver<Scalar> *&prec,
                                       FSCommunicator *com, std::string name) const
{
	//std::cerr << "3. creating solver: type = " << cntl.type << ", subtype = " << cntl.subtype << ", name = " << name << std::endl;
	GenSolver<Scalar> *solver = nullptr;
	switch(cntl.type) {
		case SolverSelection::Direct : { // direct solver
			switch(cntl.subtype) {
				default:
				case 0: {
					GenSkyMatrix<Scalar> *s = new GenSkyMatrix<Scalar>(con, cdsa, cntl.trbm, rbm);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
				case 1: {
					GenBLKSparseMatrix<Scalar> *s = new GenBLKSparseMatrix<Scalar>(con, dsa, cdsa, cntl.trbm, cntl, rbm);
					s->zeroAll();
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#ifdef USE_EIGEN3
				case 3: {
					GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *s =
						new GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(con, dsa, cdsa, true);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
				case 4: {
					GenEiSparseMatrix<Scalar,Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *s =
						new GenEiSparseMatrix<Scalar,Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(con, dsa, cdsa, true);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#ifdef EIGEN_CHOLMOD_SUPPORT
				case 5: {
          GenEiSparseMatrix<Scalar,Eigen::CholmodDecomposition<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *s =
            new GenEiSparseMatrix<Scalar,Eigen::CholmodDecomposition<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(con, dsa, cdsa, true);
          solver = (GenSolver<Scalar> *) s;
          sparse = (GenSparseMatrix<Scalar> *) s;
        } break;
#endif
#ifdef EIGEN_UMFPACK_SUPPORT
				case 6: {
				    GenEiSparseMatrix<Scalar,Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar> > > *s =
				      new GenEiSparseMatrix<Scalar,Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar> > >(con, dsa, cdsa, false);
				    solver = (GenSolver<Scalar> *) s;
				    sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
				case 7: {
				    GenEiSparseMatrix<Scalar,Eigen::SuperLU<Eigen::SparseMatrix<Scalar> > > *s =
				      new GenEiSparseMatrix<Scalar,Eigen::SuperLU<Eigen::SparseMatrix<Scalar> > >(con, dsa, cdsa, false);
				    solver = (GenSolver<Scalar> *) s;
				    sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
#endif
#ifdef USE_SPOOLES
				case 8: {
                    GenSpoolesSolver<Scalar> *s = new GenSpoolesSolver<Scalar>(con, dsa, cdsa, cntl);
                    solver = (GenSolver<Scalar> *) s;
                    sparse = (GenSparseMatrix<Scalar> *) s;
                } break;
#endif
#ifdef USE_MUMPS
				case 9: {
                    GenMumpsSolver<Scalar> *s = new GenMumpsSolver<Scalar>(con, dsa, cdsa, cntl, com);
                    solver = (GenSolver<Scalar> *) s;
                    sparse = (GenSparseMatrix<Scalar> *) s;
                }
                break;
#endif
				case 10: {
					GenDiagMatrix<Scalar> *s = new GenDiagMatrix<Scalar>(cdsa);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
				case 12: {
					Rom::GenGalerkinProjectionSolver<Scalar> *s = new Rom::GenGalerkinProjectionSolver<Scalar>(con, dsa, cdsa, cntl.pivot);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
					sparse->zeroAll();
				} break;
#ifdef USE_EIGEN3
				case 13: {
					Rom::GenEiSparseGalerkinProjectionSolver<Scalar> *s = new Rom::GenEiSparseGalerkinProjectionSolver<Scalar>(con, dsa, cdsa);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
					sparse->zeroAll();
				} break;
				case 14: {
					typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> SolverClass;
					typename WrapEiSparseMat<Scalar,SolverClass>::CtorData baseArg(con, dsa, new ConstrainedDSA(*dsa, *cdsa));
					GoldfarbIdnaniQpSolver<WrapEiSparseMat<Scalar,SolverClass>, Scalar> * s = new GoldfarbIdnaniQpSolver<WrapEiSparseMat<Scalar,SolverClass>, Scalar>
						(baseArg, cdsa, cntl.goldfarb_tol);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#ifdef EIGEN_SPARSELU_SUPPORT
				case 15: {
					GenEiSparseMatrix<Scalar,Eigen::SparseLU<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > > *s =
					  new GenEiSparseMatrix<Scalar,Eigen::SparseLU<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > >(con, dsa, cdsa, false);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
#ifdef EIGEN_SPQR_SUPPORT
				case 16: {
					GenEiSparseMatrix<Scalar,Eigen::SPQR<Eigen::SparseMatrix<Scalar> > > *s =
					  new GenEiSparseMatrix<Scalar,Eigen::SPQR<Eigen::SparseMatrix<Scalar> > >(con, dsa, cdsa, false);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
				case 17: {
					GenEiSparseMatrix<Scalar,Eigen::SparseQR<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > > *s =
						new GenEiSparseMatrix<Scalar,Eigen::SparseQR<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > >(con, dsa, cdsa, false);
					solver = (GenSolver<Scalar> *) s;
					sparse = (GenSparseMatrix<Scalar> *) s;
				} break;
#endif
				case 18: { // CholMod
					auto p = getCholmod<Scalar>(con, *dsa, *cdsa);
					solver = p.first;
					sparse = p.second;
				}
				break;
//				case 19: { // Pardiso in MKL
//					auto p = getPardiso<Scalar>(con, eqnum);
//					solver = p.first;
//					sparse = p.second;
//				}
//				break;
			}

		} break;

		case SolverSelection::Iterative : {

			switch(cntl.iterSubtype) {
				case 2: {
					filePrint(stderr," ... Node Based Sparse Matrix       ...\n");
					sparse = new GenNBSparseMatrix<Scalar>(con, cdsa);
				} break;
				default:
				case 3: {
					filePrint(stderr," ... Dof Based Sparse Matrix        ...\n");
					sparse = new GenDBSparseMatrix<Scalar>(con, dsa, cdsa);
				} break;
#ifdef USE_EIGEN3
				case 4: {
					filePrint(stderr," ... Eigen 3 Sparse Matrix          ...\n");
					sparse = new GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(con, dsa, cdsa, !cntl.unsymmetric);
				} break;
#endif
			}

			switch(cntl.precond) {
				case 1: {
					filePrint(stderr," ... Diagonal Preconditioner        ...\n");
					GenDiagMatrix<Scalar> *diag = new GenDiagMatrix<Scalar>(cdsa);
					prec = (GenSolver<Scalar> *) diag;
					spp = (GenSparseMatrix<Scalar> *) diag;
				} break;
#if defined(USE_EIGEN3) && defined(EIGEN_SUPERLU_SUPPORT) && defined(EIGEN_SUPERLU_HAS_ILU)
				case 4: {
					filePrint(stderr," ... Incomplete LU Preconditioner   ...\n");
					GenEiSparseMatrix<Scalar,Eigen::SuperILU<Eigen::SparseMatrix<Scalar> > > *ilu =
					  new GenEiSparseMatrix<Scalar,Eigen::SuperILU<Eigen::SparseMatrix<Scalar> > >(con, dsa, cdsa, false);
					ilu->getEigenSolver().options().ILU_DropTol = cntl.ilu_droptol;
					prec = (GenSolver<Scalar> *) ilu;
					spp = (GenSparseMatrix<Scalar> *) ilu;
				} break;
#endif
			}

			switch(cntl.iterType) {
				default:
				case 0: {
					filePrint(stderr," ... CG Solver is Selected          ...\n");
					solver =  new GenPCGSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar> >
						(sparse, cntl.precond, cntl.maxit, cntl.tol, cntl.maxvecsize, cntl.verbose);
				} break;
				case 1 : {
					filePrint(stderr," ... GMRES Solver is Selected       ...\n");
					GmresSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar>, GenSolver<Scalar> > *gmresSolver
						= new GmresSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar>, GenSolver<Scalar> >
							(cntl.maxit, cntl.tol, sparse, &GenSparseMatrix<Scalar>::matvec, prec, &GenSolver<Scalar>::apply, NULL,
							 &GenSolver<Scalar>::solve, (FSCommunicator *) NULL);
					if(cntl.maxvecsize > 0) gmresSolver->maxortho = cntl.maxvecsize;
					gmresSolver->verbose = cntl.verbose;
					gmresSolver->printNumber = cntl.fetiInfo.printNumber;
					solver = gmresSolver;
				} break;
				case 4: {
					filePrint(stderr," ... Bi-CG Solver is Selected       ...\n");
					solver = new GenBCGSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar> >
						(cntl.maxit, cntl.tol, sparse, prec);
				} break;
				case 5: {
					filePrint(stderr," ... CR Solver is Selected          ...\n");
					solver = new GenCRSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar> >
						(cntl.maxit, cntl.tol, sparse, prec);
				} break;
			}

		} break;
		default:
			throw std::runtime_error("Selected solver is not available.");
	}
	return solver;
}

template<class Scalar>
SolverAndMatrix<Scalar>
GenSolverFactory<Scalar>::createSolver(const Connectivity *con, const DofSetArray *dsa, int *rCN,
                                       const SolverCntl &cntl,
                                       std::string name) const
{
	//std::cerr << "4. creating solver: type = " << cntl.type << ", subtype = " << cntl.subtype << ", name = " << name << std::endl;
	std::unique_ptr<GenSparseMatrix<Scalar>> sparse;
	GenSolver<Scalar> *solver;
	switch(cntl.type) {

		case SolverSelection::Direct: {
			switch(cntl.subtype) {
				case 0: {
					auto s = std::make_unique<GenSkyMatrix<Scalar>>(con, dsa, cntl.trbm, rCN, (int)0);
					solver = static_cast<GenSolver<Scalar> *>(s.get());
					sparse = std::move(s);
				} break;
				case 1: {
					auto s = std::make_unique<GenBLKSparseMatrix<Scalar>>(con, dsa, rCN, cntl.trbm, cntl);
					s->zeroAll();
					solver = static_cast<GenSolver<Scalar> *>(s.get());
					sparse = std::move(s);
				} break;
#ifdef USE_EIGEN3
				case 3: {
					auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper>>>(con, dsa, rCN, true);
					solver = static_cast<GenSolver<Scalar> *>(s.get());
					sparse = std::move(s);
				} break;
				case 4: {
					auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper>>>(con, dsa, rCN, true);
					solver = static_cast<GenSolver<Scalar> *>(s.get());
					sparse = std::move(s);
				} break;
#ifdef EIGEN_CHOLMOD_SUPPORT
				case 5: {
          auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::CholmodDecomposition<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >>(con, dsa, rCN, true);
                    solver = static_cast<GenSolver<Scalar> *>(s.get());
          sparse = std::move(s);
        } break;
#endif
#ifdef EIGEN_UMFPACK_SUPPORT
				case 6: {
          auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>>>>(con, dsa, rCN, false);
                    solver = static_cast<GenSolver<Scalar> *>(s.get());
          sparse = std::move(s);
        } break;
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
				case 7: {
          auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::SuperLU<Eigen::SparseMatrix<Scalar> > >>(con, dsa, rCN, false);
                    solver = static_cast<GenSolver<Scalar> *>(s.get());
          sparse = std::move(s);
        } break;
#endif
#endif
#ifdef USE_SPOOLES
				case 8: {
          auto s = std::make_unique<GenSpoolesSolver<Scalar>>(con, dsa, cntl, rCN);
                    solver = static_cast<GenSolver<Scalar> *>(s.get());
          sparse = std::move(s);
        } break;
#endif
#ifdef USE_MUMPS
				case 9: {
          auto s = std::make_unique<GenMumpsSolver<Scalar>>(con, dsa, cntl, rCN);
                    solver = static_cast<GenSolver<Scalar> *>(s.get());
          sparse = std::move(s);
        }
        break;
#endif
#ifdef USE_EIGEN3
#ifdef EIGEN_SPARSELU_SUPPORT
				case 15: {
          auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::SparseLU<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > >>(con, dsa, rCN, false);
                  solver = static_cast<GenSolver<Scalar> *>(s.get());
          sparse = std::move(s);
          } break;
#endif
#ifdef EIGEN_SPQR_SUPPORT
				case 16: {
          auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::SPQR<Eigen::SparseMatrix<Scalar> > >>(con, dsa, rCN, false);
          solver = static_cast<GenSolver<Scalar> *>(s.get());
          sparse = std::move(s);
        } break;
#endif
				case 17: {
					auto s = std::make_unique<GenEiSparseMatrix<Scalar,Eigen::SparseQR<Eigen::SparseMatrix<Scalar>,Eigen::COLAMDOrdering<int> > >>(con, dsa, rCN, false);
					solver = static_cast<GenSolver<Scalar> *>(s.get());
					sparse = std::move(s);
				} break;
#endif
			}
		} break;
		default:
			throw std::runtime_error("Selected solver is not available.");
	}
	return { solver, std::move(sparse) };
}

#include <complex>

template class GenSolverFactory<double>;
template class GenSolverFactory<std::complex<double>>;