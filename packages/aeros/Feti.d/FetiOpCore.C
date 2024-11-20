#include <Feti.d/Feti.h>
#include <Feti.d/FetiOp.h>
#include "FetiOpControler.h"
#if defined(WINDOWS) || defined(MACOSX)
 #include <cfloat>
#else
 #include <climits>
#endif
#include <float.h>
#include <Driver.d/SubDomain.h>

template<>
void
GenFetiOp<DComplex>::sendInterfRBM(FSCommPattern<DComplex> *rbmPat)
{
  fprintf(stderr, "WARNING: GenFetiOp<DComplex>::sendInterfRBM() not implemented  \n");
}

template<> 
void
GenFetiOp<double>::sendInterfRBM(FSCommPattern<double> *rbmPat)
{
 locInterfRBMs.resize(numRBM*sd->interfLen());

 // sd->sendRBMs(numRBM, locRBMs, locInterfRBMs);
 sd->extractInterfRBMs(numRBM, locRBMs.data(), locInterfRBMs.data());
 sd->sendInterfRBMs(numRBM, locInterfRBMs.data(), rbmPat);

 GenCoarseSet<double> &thisSet = control->cset[sd->localSubNum()];
 thisSet.numGs = numRBM;
 thisSet.locGs = locInterfRBMs.data();

 if(QGisLocal == 0) { // If there is a Q
   int i;
   thisSet.locQGs = new double[numRBM*sd->interfLen()];
   if(rbm) {  // In the dynamic case, Q is Fi
     // multiple rhs version of multFi
     sd->multMFi(solver, thisSet.locGs, thisSet.locQGs, numRBM);
   } 
   else {
     if(control->nQ == 4) {
       for(i = 0; i< numRBM; ++i)
         sd->multDiagKbb(thisSet.locGs  + i*sd->interfLen(),
                         thisSet.locQGs + i*sd->interfLen());
     } 
     else {
       for(i = 0; i < numRBM; ++i)
          sd->multKbb(thisSet.locGs  + i*sd->interfLen(),
                      thisSet.locQGs + i*sd->interfLen());
     }
   }
 }
 else
   control->cset[sd->localSubNum()].locQGs = locInterfRBMs.data();// no preconditioning

 if (isFeti2 && isDynamic == 0) {
   if (QGisLocal == 0) {
     thisSet.locFGs = thisSet.locQGs;
   } else {
     thisSet.locFGs = new double[numRBM*sd->interfLen()];
     sd->multMFi(solver, thisSet.locGs, thisSet.locFGs, numRBM);
   }
 }
}
