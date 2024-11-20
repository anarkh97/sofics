#ifdef DISTRIBUTED
#include <mpi.h>
#endif
#include <Utils.d/dbg_alloca.h>
//#include <algo.h>

#include <Math.d/SparseMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/Vector.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/Communicator.h>
#include <Comm.d/Communicator.h>

#ifdef DISTRIBUTED
extern Communicator *structCom;
//extern SysCom *scom;
#endif

#define DBL_R_NUM

#include <Solvers.d/DSCsolver.h>
#ifdef DISTRIBUTED
#include <Solvers.d/DSCVer1/DSC_LIB/dscmain.h>
// TODO This is horrible! Global variables?
int tglobal_ns, *ta_index = nullptr, *ta_struc, *treplication,
		*ta_nonz_index;
real_number_type *ta_nonz;
#endif

DSCsolver::~DSCsolver()
{
  if(unonz) { delete unonz; unonz = 0; }
  if(adj)   { delete adj;   adj = 0;   }
  if(xadj)  { delete xadj;  xadj = 0;  }

}


DSCsolver::DSCsolver(Connectivity *cn, EqNumberer *eqNums, int s_number)
 : SparseData(eqNums, cn)
{
 unonz = 0; adj = 0; xadj = 0;
#ifdef DISTRIBUTED
 // sparse structure for matrix A
 unonz = new double[xunonz[numUncon]];

 // zero matrix A
 zeroAll();

 // scheme number 1: traditional method of solves
 //               2: selective inversion of solves
 scheme_number = 2;  //s_number;

 // Number of "nodes"
 numNodes = cn->csize();

 // This is what is needed !!!
 // These things come from the connectivity
 tglobal_ns   = numNodes; // Number of "nodes"
 if(ta_index != nullptr)
 	delete [] ta_index;
 ta_index = new int[cn->ptr().size()];
 for(size_t i = 0; i < cn->ptr().size(); ++i)
 	ta_index[i] = cn->ptr()[i];
 ta_struc     = (*cn)[0].data();

 // This comes from the eqNums
 treplication = eqNums->allWeights();

 ta_nonz_index = xunonz.data();
 ta_nonz       =  unonz; 
#endif
}

void
DSCsolver::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
 communicator->globalSum(xunonz[numUncon], unonz);
#endif
}

void
DSCsolver::print()
{
  int i;
  for(i=0; i < xunonz[numUncon]; ++i)
    fprintf(stderr," %d %e\n",rowu[i],unonz[i]);
}


// matrix assembly routine i.e. to fill the matrix A
void
DSCsolver::add(const FullM &knd, int fRow, int fCol)
{
 int i, j, m;

 int nrow  = knd.numRow();
 int nCol  = knd.numCol();

 for(i = 0; i < nrow; ++i) {
    int rowi = fRow + i;
    for(j = 0; j < nCol; ++j) {
       int colj = fCol + j;
       int mstart = xunonz[colj];
       int mstop  = xunonz[colj+1];
       for(m=mstart; m<mstop; ++m) { 
          if( rowu[m] == rowi ) {
            unonz[m] += knd[i][j];
            break;
          }
       }
    }
 }
}

// zero the matrix A, before assembly
void
DSCsolver::zeroAll()
{
  int i;
  for(i=0; i < xunonz[numUncon]; ++i)
    unonz[i] = 0.0;
}

void
DSCsolver::factor()
{
#ifdef DISTRIBUTED
 int my_pid;

 int cpuNum;

 DSC_Remove_Null_Structs(&tglobal_ns, ta_index, ta_struc, treplication);

 maxNum = DSC_Analyze (numNodes, ta_index, ta_struc, treplication );

 // KHP: just to handle the hexplate problem.
 if(numUncon < 30) maxNum = 1;

 MPI_Comm *salinasComm = 0;
 salinasComm = structCom->getCommunicator();
 if(salinasComm == 0) {
   MPI_Comm_rank(MPI_COMM_WORLD, &cpuNum);
   color   = (cpuNum < maxNum) ? 1 : 0;
   int key = (cpuNum < maxNum) ? cpuNum : cpuNum-maxNum;
   MPI_Comm_split(MPI_COMM_WORLD, color, key, &dscComm);
 } else {
   MPI_Comm_rank(*salinasComm, &cpuNum);
   color   = (cpuNum < maxNum) ? 1 : 0;
   int key = (cpuNum < maxNum) ? cpuNum : cpuNum-maxNum;
   MPI_Comm_split(*salinasComm, color, key, &dscComm);
 }
 
 if(color == 0) return;
 DSC_Open0(maxNum, &my_pid, dscComm);
 
 DSC_Order_S_N_Fact (scheme_number, tglobal_ns, ta_index,
		ta_struc, treplication, ta_nonz_index, ta_nonz);
 if (DSC_STATUS.cont_or_stop == DSC_STOP_TYPE) {
    fprintf(stderr, "Failed to setup order\n");
    DSC_Error_Display();
    exit(-1);
 }
#endif
}

void
DSCsolver::reSolve(double *rhs)
{
#ifdef DISTRIBUTED
 if(color == 1) {
   DSC_Input_Rhs_Global_Vec(rhs, neq);
   DSC_N_Solve(scheme_number);
   int *eq_number = (int *) dbg_alloca(sizeof(int)*neq);
   double *tmpSol = (double *)dbg_alloca(sizeof(double)*neq);
   DSC_Solve_Gather(eq_number, tmpSol);

   int i;
   for(i=0; i<neq; ++i)
     rhs[eq_number[i]] = tmpSol[i];
 }
   
 MPI_Comm *salinasComm = 0;
#ifdef DISTRIBUTED
 salinasComm = structCom->getCommunicator();
#endif
 if(salinasComm==0)
   MPI_Bcast(rhs, neq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 else
   MPI_Bcast(rhs, neq, MPI_DOUBLE, 0, *salinasComm);
#endif
}

void
DSCsolver::reSolve(Vector &rhs)
{
 reSolve(rhs.data());
}

void
DSCsolver::reSolve(int nRHS, double **RHS)
{
 int i;
 for(i=0; i<nRHS; ++i)
   reSolve(RHS[i]);
}

void
DSCsolver::reSolve(int nRHS, Vector *RHS)
{
 int i;
 for(i=0; i<nRHS; ++i)
   reSolve(RHS[i].data());
}
