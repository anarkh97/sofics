#ifndef _DYNAM_H_
#define _DYNAM_H_

#include <Solvers.d/Solver.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/SparseMatrix.h>

template <class Scalar>
class GenDynamMat {
 public:
   GenSolver<Scalar> *dynMat;           // use to solve (coeM*M + coeC*C + coeK*K)x = b 
   GenSolver<Scalar> *Msolver;          // use to solve Mx=b
   GenSparseMatrix<Scalar> *K;          // Stiffness Matrix
   GenSparseMatrix<Scalar> *refK;       // Stiffness Matrix for eigensolves
   GenSparseMatrix<Scalar> *C;          // Damping Matrix
   GenSparseMatrix<Scalar> *Cuc;        // constrained to unconstrained Damping Matrix
   GenSparseMatrix<Scalar> *M;          // Mass Matrix
   GenSparseMatrix<Scalar> *Muc;        // constrained to unconstrained Mass Matrix
   GenSparseMatrix<Scalar> *Mcc;        // constrained to constrained Mass Matrix
   GenSparseMatrix<Scalar> *Kuc;
   GenSparseMatrix<Scalar> *Kcc;
   GenSparseMatrix<Scalar> *Ccc;
   int  numdofs;                        // number of dof
   Rbm* rigidBodyModes;
   bool myMemory;

   // Constructors
   GenDynamMat(bool _myMemory = false) { 
     dynMat = 0; Msolver = 0; K = 0; C = 0; M = 0; Cuc = 0; Muc = 0; Mcc = 0; refK = 0; Kuc = 0; 
     Kcc = 0; Ccc = 0; rigidBodyModes = 0; myMemory = _myMemory;
   }

   GenDynamMat(GenDynamMat *d) { 
     dynMat = (*d).dynMat; Msolver = (*d).Msolver;
     K = (*d).K; refK = (*d).refK; 
     C = (*d).C; Cuc = (*d).Cuc; 
     M = (*d).M; Muc = (*d).Muc; Mcc = (*d).Mcc;
     Kuc = (*d).Kuc; Kcc = (*d).Kcc; Ccc = (*d).Ccc;
     numdofs = (*d).numdofs;
     rigidBodyModes = (*d).rigidBodyModes;
     myMemory = false;
   }

   // Destructor
   ~GenDynamMat() {
     if(myMemory) {
       if(dynMat) delete dynMat;
       if(Msolver) delete Msolver;
       if(K) delete K;
       if(refK) delete refK;
       if(C) delete C;
       if(Cuc) delete Cuc;
       if(M) delete M;
       if(Muc) delete Muc;
       if(Mcc) delete Mcc;
       if(Kuc) delete Kuc;
       if(Kcc) delete Kcc;
       if(Ccc) delete Ccc;
       if(rigidBodyModes) delete rigidBodyModes;
     }
   }
};

typedef GenDynamMat<double> DynamMat;

#endif
