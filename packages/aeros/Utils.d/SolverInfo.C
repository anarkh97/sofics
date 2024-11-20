#include <Utils.d/SolverInfo.h>
#include <Utils.d/DistHelper.h>

int
SolverInfo::classifySolver()
{
	// return value:
	// 1 : if the selected solver is suitable for indefinite systems arising for example
	//     when the "multipliers" constraint method is used for equality constraints
	// 2 : if the selected solver is suitable for equality and/or inequality-constrained
	//     quadratic programs arising for example when the "multipliers" constraint
	//     method is used for inequality constraints
	// 0 : if the selected solver is known to be unsuitable for both indefinite systems
	//     and constrained QPs.
	// -1: if the selected solver is not able to be classified

	switch(solvercntl->type) {
		case SolverSelection::Direct: { // direct solvers
			switch(solvercntl->subtype) {
				case 0: // skyline
					return 0;
				case 1: // sparse
					return 0;
				case 2: // blocksky
					return 0;
				case 3: // llt
					return 0;
				case 4: // ldlt
					return 0;
				case 5: // cholmod
					return 0;
				case 6: // umfpack
#if defined(USE_EIGEN3) && defined(EIGEN_UMFPACK_SUPPORT)
					return 1;
#else
					return 0;
#endif
				case 7: // superlu
#if defined(USE_EIGEN3) && defined(EIGEN_SUPERLU_SUPPORT)
					return 1;
#else
					return 0;
#endif
				case 8: // spooles
#ifdef USE_SPOOLES
					return (solvercntl->pivot) ? 1 : 0;
#else
					return 0;
#endif
				case 9: // mumps
#ifdef USE_MUMPS
					return (solvercntl->pivot) ? 1 : 0;
#else
					return 0;
#endif
				case 10: // diagonal
					return 0;
				case 12: // dbsgal
					return 0;
				case 13: // eisgal
					return 2;
				case 14: // goldfarb
#ifdef USE_EIGEN3
					return 2;
#else
					return 0;
#endif
				case 15: // splu
#ifdef USE_EIGEN3
					return 1;
#else
					return 0;
#endif
				case 16: // ssqr
#if defined(USE_EIGEN3) && defined(EIGEN_SPQR_SUPPORT)
					return 1;
#else
					return 0;
#endif
				case 17: // spqr
#ifdef USE_EIGEN3
					return 1;
#else
					return 0;
#endif
				default: return -1;
			}
		} break;
		case SolverSelection::Iterative : { // iterative solvers
			switch(solvercntl->iterType) {
				case 0: // cg
					return 0;
				case 1: // gmres
					return 1;
				case 4: // bcg
					return 1;
				case 5: // cr
					return 1;
				default: return -1;
			}
		} break;
		case SolverSelection::Feti : { // FETI solvers
			switch(solvercntl->fetiInfo.version) {
				case FetiInfo::feti1 :
					return 0;
				case FetiInfo::feti2 :
					return 0;
				case FetiInfo::feti3 :
					return 0;
				case FetiInfo::fetidp :
					return 2;
				default: return -1;
			}
		} break;
		case SolverSelection::FetiLib :
			return 2;
		default:
			return -1;
	}
}

void
SolverInfo::activatePiecewise()
{
  if(!isNonLin()) {
    if(probType == SolverInfo::Static || probType == SolverInfo::None)
      probType = SolverInfo::NonLinStatic;
    else if(probType == SolverInfo::Dynamic)
      probType = SolverInfo::NonLinDynam;
    else if(probType == SolverInfo::TempDynamic) {
      order = 1;
      probType = SolverInfo::NonLinDynam;
    }
    setNewton(std::numeric_limits<int>::max());
    getNLInfo().stepUpdateK = std::numeric_limits<int>::max();
    getNLInfo().maxiter = 1;
    if(piecewise) {
      getNLInfo().linearelastic = 1;
      getNLInfo().dlambda = piecewise_dlambda;
      getNLInfo().maxLambda = piecewise_maxLambda;
    }
    else if(freeplay) {
      getNLInfo().linearelastic = 2;
    }
  }
}
