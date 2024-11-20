/*
 * MappedAssembledSolver.cpp
 *
 *  Created on: Apr 15, 2009
 *      Author: michel
 */

#include "MappedAssembledSolver.h"
#include <Driver.d/Domain.h>
#include <Driver.d/Mpc.h>

DOFMap *getDofMaps(int size) {
  DOFMap *res = new DOFMap[size];
  for(int i = 0; i < size; ++i) {
    res[i].ndofs = 1;
    res[i].dofs = new int[1];
    res[i].coefs = new double[1];
    res[i].dofs[0] = i;
    res[i].coefs[0] = 1;
    res[i].rhs = 0;
  }
  return res;
}

ConstrainedDSA *
Domain::makeMaps(DofSetArray *dsa, ConstrainedDSA *cdsa, DOFMap *baseMap, DOFMap *eqMap)
{
  int size = dsa->size();
  int eqSize = cdsa->size();
  ResizeArray<int> constrainedDOFs(0, numLMPC);
  int mpCount = 0;
  for(int i = 0; i < numLMPC; ++i)
    if(!lmpc[i]->terms[0].isNull()) {
      int dof = dsa->locate(lmpc[i]->terms[0].nnum, 1 << lmpc[i]->terms[0].dofnum);
      if(dof >= 0)
        constrainedDOFs[mpCount++] = dof;
    }

  // PJSA also constrain the lagrange multiplier dofs, if any
  for(int i = 0; i < dsa->numNodes(); ++i) {
    int dof;
    if((dof = dsa->locate(i, DofSet::LagrangeE)) >= 0 || (dof = dsa->locate(i, DofSet::LagrangeI)) >= 0)
      constrainedDOFs[mpCount++] = dof;
  }

  ConstrainedDSA *cdsa2 = new ConstrainedDSA(*dsa, *cdsa, mpCount, constrainedDOFs.data());

  for(int i = 0; i < size; ++i)
    baseMap[i].ndofs = -2;
  for(int i = 0; i < eqSize; ++i)
    eqMap[i].ndofs = -2;

  for(int i = 0; i < numLMPC; ++i)
    if(!lmpc[i]->terms[0].isNull()) {
      int dof = dsa->locate(lmpc[i]->terms[0].nnum, 1 << lmpc[i]->terms[0].dofnum);
      int eqDof = cdsa->locate(lmpc[i]->terms[0].nnum, 1 << lmpc[i]->terms[0].dofnum);
      if(dof >= 0) {
        baseMap[dof].dofs = new int[lmpc[i]->nterms-1];
        baseMap[dof].coefs = new double[lmpc[i]->nterms-1];
        baseMap[dof].rhs = lmpc[i]->rhs.r_value/lmpc[i]->terms[0].coef.r_value;
        double c1 = 1.0/lmpc[i]->terms[0].coef.r_value;
        int nCoefs = 0;
        for(int j = 1; j < lmpc[i]->nterms; ++j) {
          int jDof = dsa->locate(lmpc[i]->terms[j].nnum, 1 << lmpc[i]->terms[j].dofnum);
          if(jDof < 0)
            continue;
          baseMap[dof].dofs[nCoefs] = jDof;
          baseMap[dof].coefs[nCoefs] = -c1*lmpc[i]->terms[j].coef.r_value;
          nCoefs++;
        }
        baseMap[dof].ndofs = nCoefs;
      }
      if(eqDof >= 0) {
        eqMap[eqDof].dofs = new int[lmpc[i]->nterms-1];
        eqMap[eqDof].coefs = new double[lmpc[i]->nterms-1];
        eqMap[eqDof].rhs = lmpc[i]->rhs.r_value/lmpc[i]->terms[0].coef.r_value;
        double c1 = 1.0/lmpc[i]->terms[0].coef.r_value;
        int nCoefs = 0;
        for(int j = 1; j < lmpc[i]->nterms; ++j) {
          int jDof = cdsa2->locate(lmpc[i]->terms[j].nnum, 1 << lmpc[i]->terms[j].dofnum);
          if(jDof < 0)
            continue;
          eqMap[eqDof].dofs[nCoefs] = jDof;
          eqMap[eqDof].coefs[nCoefs] = -c1*lmpc[i]->terms[j].coef.r_value;
          nCoefs++;
        }
        eqMap[eqDof].ndofs = nCoefs;
      }
    }

  for(int i = 0; i < dsa->size(); ++i) {
    if(baseMap[i].ndofs < 0) {
      baseMap[i].ndofs = 1;
      baseMap[i].dofs = new int[1];
      baseMap[i].coefs = new double[1];
      baseMap[i].dofs[0] = i;
      baseMap[i].coefs[0] = 1;
      baseMap[i].rhs = 0;
    }
  }

  for(int i = 0; i < dsa->numNodes(); ++i) {
    int dofNums[DofSet::max_known_dof];
    int mappedNums[DofSet::max_known_dof];
    int nd = cdsa->number(i, (*dsa)[i], dofNums);
    cdsa2->number(i, (*dsa)[i], mappedNums);
    for(int j = 0; j < nd; ++j) {
      if(dofNums[j] >= 0 && mappedNums[j] >= 0 && eqMap[dofNums[j]].ndofs < 0) {
        eqMap[dofNums[j]].ndofs = 1;
        eqMap[dofNums[j]].dofs = new int[1];
        eqMap[dofNums[j]].coefs = new double[1];
        eqMap[dofNums[j]].dofs[0] = mappedNums[j];
        eqMap[dofNums[j]].coefs[0] = 1;
        eqMap[dofNums[j]].rhs = 0;
      }
    }
  }

  return cdsa2;
}
