#include <iostream>
#include <Math.d/CuCSparse.h>
#include <Math.d/matrix.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/Vector.h>
#include <Driver.d/Mpc.h>

inline void init(DComplex &d) { d = DComplex(0.0,0.0); }
inline void init(double &d) { d = 0.0; }

template<class Scalar>
void
GenCuCSparse<Scalar>::clean_up()
{
 Kuc.clear();
 Kuc.shrink_to_fit();
}

template<class Scalar>
void
GenCuCSparse<Scalar>::zeroAll()
{
 if(numConstrained > 0) {
   for(int i=0; i<xunonz[numConstrained]; ++i)
     ScalarTypes::initScalar(Kuc[i], 0.0, 0.0);
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::negate()
{
 int i;
 for(i=0; i<xunonz[numConstrained]; ++i)
   Kuc[i] = -Kuc[i];
}

template<class Scalar>
void
GenCuCSparse<Scalar>::print()
{
 int i, mstart, mstop, m;
 std::cerr << "numConstrained = " << numConstrained << std::endl;
 for(i=0; i<numConstrained; ++i) {
   mstart = xunonz[i];
   mstop  = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     std::cerr << "Kuc(" << rowu[m]+1 << "," << i+1 << ") = " << Kuc[m] << std::endl;
 }
}

template<class Scalar>
GenCuCSparse<Scalar>::GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const int *bc) :
SparseData(con,dsa,bc)
{
 // Allocate Kuc and initialize to zero.
 Kuc.assign(xunonz[numConstrained], Scalar(0));

}

template<class Scalar>
GenCuCSparse<Scalar>::GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const DofSetArray *c_dsa) :
SparseData(con,dsa,c_dsa)
{
 // Allocate Kuc
 if(numConstrained>0) {
   Kuc.assign(xunonz[numConstrained], Scalar(0));
 }
}

template<class Scalar>
GenCuCSparse<Scalar>::GenCuCSparse(const Connectivity *con, const DofSetArray *dsa,
                                   const int *glBoundMap, const int *glInternalMap)
: SparseData(con,dsa,glBoundMap,glInternalMap)
{
 // Allocate Kuc and initialize it to zero.
 Kuc.resize(xunonz[numConstrained]);

 // Initialize it to zero
 zeroAll();

}

template<class Scalar>
GenCuCSparse<Scalar>::GenCuCSparse(LMPCons **mpc, int numMPC, const DofSetArray *c_dsa) :
  SparseData(mpc, numMPC, c_dsa)
{
  Kuc.resize(xunonz[numConstrained]);
  int i, mstart, mstop, m;
  for(i=0; i<numConstrained; ++i) {
    mstart = xunonz[i];
    mstop  = xunonz[i+1];
    for(m=mstart; m<mstop; ++m) {
      Kuc[m] = mpc[i]->terms[m].template val<Scalar>();
    }
  }
}

template<class Scalar>
GenCuCSparse<Scalar>::GenCuCSparse(int numInterface, const int *allBoundDofs,
                                   int numModes, Scalar *modes, int ldm)
 : SparseData(numInterface, allBoundDofs, numModes, ldm)
{
 if(ldm < 0) ldm = numInterface;
 // Allocate Kuc and initialize it to zero.
 Kuc.resize(xunonz[numConstrained]);
 neq=numInterface;

 // fill the matrix with the modes
 int i, mstart, mstop, m;
 for(i=0; i<numModes; ++i) {
   mstart = xunonz[i];
   mstop  = xunonz[i+1];
   for(m=mstart; m<mstop; ++m) {
     Kuc[m] = modes[(m-mstart) + i*ldm];
   }
 }
}

template<class Scalar>
GenCuCSparse<Scalar>::GenCuCSparse(int cnt, int neq, gsl::span<int> count, gsl::span<int> list, std::vector<Scalar> val)
 : SparseData(cnt, count.data(), list.data()),  Kuc(std::move(val))
{
	this->neq = neq;
}

template<class Scalar>
void
GenCuCSparse<Scalar>::doWeighting(Scalar *weight)
{
  for(int i=0; i<xunonz[numConstrained]; ++i)
    Kuc[i] *= weight[rowu[i]];
}

template<class Scalar>
void
GenCuCSparse<Scalar>::doWeighting(int *weight)
{
  for(int i=0; i<xunonz[numConstrained]; ++i)
    Kuc[i] *= Scalar(weight[rowu[i]]);
}

template<class Scalar>
double
GenCuCSparse<Scalar>::getMemoryUsed() const
{
 return 8*xunonz[numConstrained]/(1024.0*1024.0);
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const
{
 multSubtract(rhs.data(),result.data());
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multSubtract(const Scalar *rhs, Scalar *result) const
{
 // mult multiplies a bc value vector by Kuc to
 // get an addition to the user-defined rhs vector.

 int i, m, mstart, mstop;

 for(i=0; i<numConstrained; ++i) {
   if(rhs[i] == 0.0) continue;
   mstart = xunonz[i];
   mstop  = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     result[rowu[m]] -= Kuc[m]*rhs[i];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const
{
 mult(rhs.data(),result.data());
}

template<class Scalar>
void
GenCuCSparse<Scalar>::mult(const Scalar *rhs, Scalar *result) const
{
 for(int i=0; i<numRow(); ++i) result[i] = 0;
 int i, m, mstart, mstop;
 for(i=0; i<numConstrained; ++i) {
   if(rhs[i] == 0.0) continue;
   mstart = xunonz[i];
   mstop  = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     result[rowu[m]] += Kuc[m]*rhs[i];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::transposeMult(const Scalar *rhs, Scalar *result) const
{
 int i, m, mstart, mstop;

 for(i = 0; i < numConstrained; i++) {
   mstart = xunonz[i];
   mstop  = xunonz[i+1];
   for( m = mstart; m < mstop; m++)
     result[rowu[m]] = 0.0;
 }

 for(i=0; i<numConstrained; ++i) {
   if( rhs[i] == 0.0) continue;
   mstart = xunonz[i];
   mstop  = xunonz[i+1];

   // zero the result before accumulating the answer
   for(m=mstart; m<mstop; ++m)
     // result[i] += Kuc[m]*rhs[rowu[m]];
     result[rowu[m]] += Kuc[m] * rhs[i];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::add(const FullSquareMatrix &kel, const int *dofs)
{
 int i, j, m, ri;
 int kndof = kel.dim();
 for(i=0; i<kndof; ++i) {
   if((ri = unconstrNum[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(constrndNum[dofs[j]] == -1) continue;
     for(m=xunonz[constrndNum[dofs[j]]]; m<xunonz[constrndNum[dofs[j]]+1]; ++m) {
       if(rowu[m] == ri) {
         Kuc[m] += kel[i][j];
         break;
       }
     }
   }
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multIdentity(Scalar **result) const
{
 int dof,m;
 for(dof=0; dof<numConstrained; ++dof){
   int mstart = xunonz[dof];
   int mstop  = xunonz[dof+1];
   for(m=mstart; m<mstop; ++m)
     result[dof][rowu[m]] += Kuc[m];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multIdentity(Scalar *result) const
{
 int dof,m;
 for(dof=0; dof<numConstrained; ++dof){
   int mstart = xunonz[dof];
   int mstop  = xunonz[dof+1];
   for(m=mstart; m<mstop; ++m)
     (result +dof*numRow())[rowu[m]] += Kuc[m];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multIdentity(Scalar **result, int start, int stop) const
{
 int dof,m;
 stop = (stop<0) ? numConstrained : stop;

 for(dof=start; dof<stop; ++dof){

   for (m=0; m<neq; ++m)
     init(result[dof-start][m]);

   int mstart = xunonz[dof];
   int mstop  = xunonz[dof+1];

   for(m=mstart; m<mstop; ++m)
     result[dof-start][rowu[m]] += Kuc[m];

 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multSub(const Scalar *rhs, Scalar *result) const
{
  int i;
  for(i=0; i<numConstrained; ++i) {
    int mstart = xunonz[i];
    int mstop  = xunonz[i+1];
    int m;
    for(m=mstart; m<mstop; ++m)
      result[i] -= Kuc[m]*rhs[rowu[m]];
  }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::transposeMultSubtract(const Scalar *rhs, Scalar *result) const
{
  int i;
  for(i=0; i<numConstrained; ++i) {
    if(rhs[i] == 0.0) continue;
    int mstart = xunonz[i];
    int mstop  = xunonz[i+1];
    int m;
    for(m=mstart; m<mstop; ++m) {
      if(rowu[m] < 0)
        std::cerr << "**Error in GenCuCSparse<Scalar>::transposeMultSubtract \n";
      result[rowu[m]] -= Kuc[m]*rhs[i];
    }
  }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::transposeMultSubtractClaw(const Scalar *rhs, Scalar *result, int numUserDisp, int *clawdofs) const {

  int i;
  for(i=0; i < numUserDisp; ++i) {
    if(rhs[i] == 0.0) continue;
    int ci = clawdofs[i];
    int mstart = xunonz[ci];
    int mstop  = xunonz[ci+1];
    int m;

    for(m=mstart; m<mstop; ++m)
      result[rowu[m]] -= Kuc[m]*rhs[i];
  }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multSub(int nRHS, const Scalar **rhs, Scalar **result) const
{

  int i,n=0;
  for( ; n < nRHS-3; n += 4) {
    for(i=0; i<numConstrained; ++i) {
      int mstart = xunonz[i];
      int mstop  = xunonz[i+1];
      int m;
      for(m=mstart; m<mstop; ++m) {
        result[n][i]   -= Kuc[m]*rhs[n][rowu[m]];
        result[n+1][i] -= Kuc[m]*rhs[n+1][rowu[m]];
        result[n+2][i] -= Kuc[m]*rhs[n+2][rowu[m]];
        result[n+3][i] -= Kuc[m]*rhs[n+3][rowu[m]];
      }
    }
  }
  for( ; n < nRHS-1; n += 2) {
    for(i=0; i<numConstrained; ++i) {
      int mstart = xunonz[i];
      int mstop  = xunonz[i+1];
      int m;
      for(m=mstart; m<mstop; ++m) {
        result[n][i]   -= Kuc[m]*rhs[n][rowu[m]];
        result[n+1][i] -= Kuc[m]*rhs[n+1][rowu[m]];
      }
    }
  }
  for(; n < nRHS; ++n)
   multSub(rhs[n],result[n]);

}

template<class Scalar>
void
GenCuCSparse<Scalar>::multAdd(const Scalar *rhs, Scalar *result) const
{
  int i, m;
  for(i=0; i<numConstrained; ++i) {
    int mstart = xunonz[i];
    int mstop  = xunonz[i+1];
    result[i] = 0.0;
    for(m=mstart; m<mstop; ++m) {
      result[i] += Kuc[m]*rhs[rowu[m]];
    }
  }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::transposeMultAdd(const Scalar *rhs, Scalar *result) const
{
  int i;
  for(i=0; i<numConstrained; ++i) {
    if(rhs[i] == 0.0) continue;
    int mstart = xunonz[i];
    int mstop  = xunonz[i+1];
    int m;
    for(m=mstart; m<mstop; ++m)
      result[rowu[m]] += Kuc[m]*rhs[i];
  }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::addBoeing(int nl, const int *Kai, const int *Kaj,
                                const double *nz, int *map, Scalar multiplier)
{
 int i, j;
 int dofI, dofJ;
 for(i = 0; i < nl; ++i) {
   if(map[i] == -1) continue;
   dofI = unconstrNum[map[i]];
   if(dofI < 0) {
      if ( (dofI = constrndNum[map[i]]) < 0 ) continue;
      for(j = Kai[i]; j < Kai[i+1]; ++j) {
        if(map[Kaj[j-1]-1] == -1) continue;
        dofJ = unconstrNum[map[Kaj[j-1]-1]];
        if(dofJ < 0) continue;
        int m;
        for(m=xunonz[dofI]; m<xunonz[dofI+1]; ++m)
           if(rowu[m] == dofJ) {
             Kuc[m] += (nz[j-1]*multiplier);
             break;
           }
      }

   } else {
     for(j = Kai[i]; j < Kai[i+1]; ++j) {
       if(map[Kaj[j-1]-1] == -1) continue;
       dofJ = constrndNum[map[Kaj[j-1]-1]];
       if(dofJ < 0) continue;
       int m;
       for(m=xunonz[dofJ]; m<xunonz[dofJ+1]; ++m)
          if(rowu[m] == dofI) {
            Kuc[m] += (nz[j-1]*multiplier);
            break;
          }
     }
   }
 }
}

template<class Scalar>
long
GenCuCSparse<Scalar>::size() const
{
  return ((xunonz.size()) ? xunonz[numConstrained] : 0) ;
}

template<class Scalar>
void
GenCuCSparse<Scalar>::addDiscreteMass(int dof, Scalar mass)
{
  if(dof < 0) return;
  int m, ri;
  if((ri = unconstrNum[dof]) == -1) return;
  if(constrndNum[dof] == -1) return;
  for(m=xunonz[constrndNum[dof]]; m<xunonz[constrndNum[dof]+1]; ++m) {
    if(rowu[m] == ri) {
      Kuc[m] += mass;
      break;
    }
  }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::add(int dofi, int dofj, Scalar s)
{
 if((dofi < 0) || (dofj < 0)) return;
 int m, rowi, colj;
 if(unconstrNum.size() != 0) {
   if((rowi = unconstrNum[dofi]) == -1 || (colj = unconstrNum[dofj]) == -1) return;
 }
 else { rowi = dofi; colj = dofj; }

 // add K[i][j]
 for(m=xunonz[constrndNum[dofj]]; m<xunonz[constrndNum[dofj]+1]; ++m) {
  if(rowu[m] == rowi) {
    Kuc[m] += s;
    break;
  }
 }

 // add K[j][i]
 for(m=xunonz[constrndNum[dofi]]; m<xunonz[constrndNum[dofi]+1]; ++m) {
  if(rowu[m] == colj) {
    Kuc[m] += s;
    break;
  }
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::multSubAdd(const Scalar *rhs, Scalar *result) const
{
 int i, m, mstart, mstop;

 for(i=0; i<numConstrained; ++i) {
   mstart = xunonz[i];
   mstop  = xunonz[i+1];

   for(m=mstart; m<mstop; ++m)
     result[i] += Kuc[m]*rhs[rowu[m]];
 }
}

// PJSA: I think there is a bug in multAdd and transposeMult
// but i don't want to mess with them so I am writing my own versions to test & use
template<class Scalar>
void
GenCuCSparse<Scalar>::multAddNew(const Scalar *rhs, Scalar *result) const
{
 int i, m, mstart, mstop;
 for(i=0; i<numConstrained; ++i) {
   if(rhs[i] == 0.0) continue;
   mstart = xunonz[i];
   mstop = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     result[rowu[m]] += Kuc[m]*rhs[i];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::transposeMultNew(const Scalar *rhs, Scalar *result) const
{
 int i, m, mstart, mstop;
 for(i = 0; i < numConstrained; i++) {
   result[i] = 0.0;
 }
 for(i=0; i<numConstrained; ++i) {
   mstart = xunonz[i];
   mstop = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     result[i] += Kuc[m]*rhs[rowu[m]];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::transposeMultAddNew(const Scalar *rhs, Scalar *result) const
{
 int i, m, mstart, mstop;
 for(i=0; i<numConstrained; ++i) {
   mstart = xunonz[i];
   mstop = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     result[i] += Kuc[m]*rhs[rowu[m]];
 }
}

template<class Scalar>
void
GenCuCSparse<Scalar>::transposeMultSubNew(const Scalar *rhs, Scalar *result) const
{
 int i, m, mstart, mstop;
 for(i=0; i<numConstrained; ++i) {
   mstart = xunonz[i];
   mstop = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     result[i] -= Kuc[m]*rhs[rowu[m]];
 }
}

// y <- alpha.A.x + beta.y
template<class Scalar>
void
GenCuCSparse<Scalar>::mult(const Scalar *rhs, Scalar *result, Scalar alpha, Scalar beta) const
{
 int i, m, mstart, mstop;
 if(beta!=1.0)
   for(i=0; i<numUncon; ++i)
     result[i] = beta*result[i];

 if(alpha!=0.0)
   for(i=0; i<numConstrained; ++i) {
     if(rhs[i] == 0.0) continue;
     mstart = xunonz[i];
     mstop  = xunonz[i+1];
     for(m=mstart; m<mstop; ++m)
       result[rowu[m]] += alpha*Kuc[m]*rhs[i];
   }
}

// y <- alpha.A^t.x + beta.y
template<class Scalar>
void
GenCuCSparse<Scalar>::trMult(const Scalar *rhs, Scalar *result, Scalar alpha, Scalar beta) const
{
 int i, m, mstart, mstop;
 if(beta!=1.0)
   for(i=0; i<numConstrained; ++i)
     result[i] = beta*result[i];

 if(alpha!=0.0)
   for(i=0; i<numConstrained; ++i) {
     mstart = xunonz[i];
     mstop  = xunonz[i+1];
     for(m=mstart; m<mstop; ++m)
       result[i] += alpha*Kuc[m]*rhs[rowu[m]];
   }
}

