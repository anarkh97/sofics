#include <Paral.d/MDModal.h>
#include <Paral.d/MDOp.h>
#include <Problems.d/ModalBase.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/DecDomain.h>
#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif
#include <Driver.d/SysState.h>
#include <Utils.d/DistHelper.h>
#include <Math.d/EiSparseMatrix.h>
#include "MDDynam.h"

class DenseMatrix;
extern Communicator *structCom;

MultiDomainModal::MultiDomainModal(Domain *d) : modalOps(*(new ModalOps)), MDModalBase(d){

  previousCq = 0;
  previousDisp = 0;
}

//------------------------------------------------------------------------------

MultiDomainModal::~MultiDomainModal(){

  if(previousCq) delete previousCq;
  if(previousDisp) delete previousDisp;
}

//------------------------------------------------------------------------------

void MultiDomainModal::projectForce(DistrVector &fullF, Vector& modalF){
/*PRE: ModeBase::modesFl have been populated by a call to
         ModeBase::populateFlexModes
 POST: return in modalF, fullF projected onto modesFl
 NOTE: this projection is not valid for displacements, velocities or accelerations
*/
  for(int i = 0; i < numModes; ++i)
    modalF[i] = (*modesFl)[i] * fullF;
}

//------------------------------------------------------------------------------

void MultiDomainModal::expand(const Vector &modalV, DistrVector& fullV){
/*PRE: there is at least one modesFl; modesFl are populated
 POST: return in fullV, modalV projected into full space
 NOTE: this expansion is not valid for forces
*/
  fullV.linC((*modesFl)[0], modalV[0]);
  for(int i = 1; i < numModes; ++i)
    fullV.linAdd(modalV[i], (*modesFl)[i]);
}

//------------------------------------------------------------------------------

void MultiDomainModal::processLastOutput(){

  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

//------------------------------------------------------------------------------

void MultiDomainModal::preProcess(){
/*NOTE: call to populateFlexModes gives 1 for 2nd argument to indicate all
          modes should be read, including those with zero frequency
        see also ModalBase::populateFlexModes
*/
  preProcessBase();
  populateFlexModes(1.0, 1);
  numModes = numFlex;
}

//------------------------------------------------------------------------------

void MultiDomainModal::preProcessSA(){

  filePrint(stderr," *** ERROR: MultiDomainModal::preProcessSA is not implemented.\n");
  exit(-1);
}

//------------------------------------------------------------------------------

void MultiDomainModal::postProcessSA(ModalOps *,Vector &sol){

  filePrint(stderr," *** ERROR: MultiDomainModal::postProcessSA is not implemented.\n");
  exit(-1);
}

//------------------------------------------------------------------------------

void MultiDomainModal::getTimes(double &dt, double &tmax){

  dt   = domain->solInfo().getTimeStep();
  tmax = domain->solInfo().tmax;
}

//------------------------------------------------------------------------------

void MultiDomainModal::getInitState(SysState<Vector> &state){

  initStateBase(state.getDisp(), state.getVeloc(),
    state.getAccel(), state.getPrevVeloc());
}

//------------------------------------------------------------------------------

void MultiDomainModal::getConstForce(Vector &constF){

  MultiDomainOp mdop(&MultiDomainOp::getConstForce, decDomain->getAllSubDomains(), fullTmpGrav, (SubDOp*)0);
  threadManager->execParal(decDomain->getNumSub(), &mdop);

  projectForce(*fullTmpGrav, constF);
}

//------------------------------------------------------------------------------

void MultiDomainModal::addConstForceSensitivity(Vector &constF){

  filePrint(stderr," *** ERROR: MultiDomainModal::addConstForceSensitivity is not implemented.\n");
  exit(-1); 
}

//------------------------------------------------------------------------------

void MultiDomainModal::getSensitivityStateParam(double &tol, double &ratioTol) {

  tol = domain->solInfo().sensitivityTol;
  ratioTol = domain->solInfo().ratioSensitivityTol;
}

//------------------------------------------------------------------------------

void MultiDomainModal::getSteadyStateParam(int &flag, int &min, int &max, double &tol){

  flag = domain->solInfo().steadyFlag;
  min  = domain->solInfo().steadyMin;
  max  = domain->solInfo().steadyMax;
  tol  = domain->solInfo().steadyTol;
}

//------------------------------------------------------------------------------

void MultiDomainModal::getNewMarkParameters(double &beta, double &gamma,
  double &alphaf, double &alpham){

  beta  = domain->solInfo().newmarkBeta;
  gamma = domain->solInfo().newmarkGamma;
  alphaf = domain->solInfo().newmarkAlphaF;
  alpham = domain->solInfo().newmarkAlphaM;
}

//------------------------------------------------------------------------------

void MultiDomainModal::getInitialTime(int &tIndex, double &time){

  tIndex = domain->solInfo().initialTimeIndex;
  time   = domain->solInfo().initialTime;
}

//------------------------------------------------------------------------------

double
MultiDomainModal::getInitialForceNorm()
{
  return domain->solInfo().initExtForceNorm;
}

//------------------------------------------------------------------------------

ModalOps* MultiDomainModal::buildOps(double mcoef, double ccoef, double kcoef){
//PRE: ModalProbDescr is instantiated; modeData is populated
// POST: operator for time integration
//

  ModalParams modalParams = domain->solInfo().readInModes[domain->solInfo().modal_id.front()]; // assume only one basis is parsed for modal analysis 

  // check to see if basis satisfies orthogonality property
  bool checkBasis = false;
  if(modalParams.tolerance != 0)
    checkBasis = true;

  switch(modalParams.type) {
    case ModalParams::Eigen:
    case ModalParams::Undefined : { // if eigen basis, all operators are diagonal
      // every process gets an identical copy of modalOps
      modalOps.M       = new DiagonalMatrix(numModes);
      modalOps.Msolver = new DiagonalMatrix(numModes);
      // see below for damping matrix
      modalOps.K      = new DiagonalMatrix(numModes);
      modalOps.dynMat = new DiagonalMatrix(numModes);

      modalOps.M->setDiag(1.0);
      modalOps.Msolver->setDiag(1.0); // Inverse of M

      if(checkBasis) {
        filePrint(stderr," ... Validating orthogonality of modal eigen basis ... \n");
        SparseMatrix **subM = new SparseMatrix * [decDomain->getNumSub()];
        SparseMatrix **subK = new SparseMatrix * [decDomain->getNumSub()];
#if defined(USE_EIGEN3)
        int sizeM = 0, sizeK = 0;
        double normM = 0, normK = 0;
        for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {
          AllOps<double> allOps;
          allOps.M = subM[iSub] = decDomain->getSubDomain(iSub)->constructEiSparse<double>();
          allOps.K = subK[iSub] = decDomain->getSubDomain(iSub)->constructEiSparse<double>();
          decDomain->getSubDomain(iSub)->makeSparseOps<double>(allOps, 0, 0, 0);
          sizeM += static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().size();
          sizeK += static_cast<EiSparseMatrix*>(allOps.K)->getEigenSparse().size();
          normM += static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().norm();
          normK += static_cast<EiSparseMatrix*>(allOps.K)->getEigenSparse().norm();
        }
        sizeM = structCom->globalSum(sizeM);
        sizeK = structCom->globalSum(sizeK);
        if (sizeM == 0 || sizeK == 0) {
          filePrint(stderr, " *** ERROR: Mass matrix and/or stiffness matrix could not be constructed. Exiting...\n\n");
          exit(-1);
        }
        normM = structCom->globalSum(normM);
        if (normM == 0) {
          filePrint(stderr, " ... WARNING: The Frobenius norm of the mass matrix is zero.\n");
        }
        normK = structCom->globalSum(normK);
        if (normK == 0) {
          filePrint(stderr, " ... WARNING: The Frobenius norm of the stiffness matrix is zero.\n");
        }
#else
        for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {
          AllOps<double> allOps;
          allOps.M = subM[iSub] = decDomain->getSubDomain(iSub)->constructDBSparseMatrix<double>();
          allOps.K = subK[iSub] = decDomain->getSubDomain(iSub)->constructDBSparseMatrix<double>();
          decDomain->getSubDomain(iSub)->makeSparseOps<double>(allOps, 0, 0, 0);
        }
#endif
        SubDOp M(decDomain->getNumSub(), subM);
        SubDOp K(decDomain->getNumSub(), subK);
 
        GenDistrVectorSet<double> tPhiM(numModes, decDomain->solVecInfo());
        GenDistrVectorSet<double> tPhiK(numModes, decDomain->solVecInfo());
        for(int i = 0; i<numModes; ++i){
          M.mult((*modesFl)[i], tPhiM[i]);
          K.mult((*modesFl)[i], tPhiK[i]);
        }

#ifdef USE_EIGEN3
        DenseMatrix redMass(numModes);
        DenseMatrix redStif(numModes);
        double Knorm = 0.0;
        for(int row = 0; row < numModes; ++row) { 
          for(int col = 0; col < numModes; ++col) {
            redMass[row+numModes*col] = tPhiM[row]*(*modesFl)[col];
            redStif[row+numModes*col] = tPhiK[row]*(*modesFl)[col];
            if(col == row) {// subtract identity and circular freqs
              redMass[row+numModes*col] -= 1;
              double Kdiag = freqs[col]*freqs[col];
              redStif[row+numModes*col] -= Kdiag;
              Knorm += Kdiag*Kdiag; 
            }
          }
        }
        Knorm = sqrt(Knorm);

        if(redMass.norm()/sqrt(numModes) > modalParams.tolerance) {
          filePrint(stderr," ... Modal Basis does not satisfy mass orthogonality tolerance, exiting. \n");
          exit(-1);
        }
        if(redStif.norm()/Knorm > modalParams.tolerance) {
          filePrint(stderr," ... Modal Basis does not satisfy stiffness decoupling tolerance, exiting. \n");
          exit(-1);
        }
#endif
      }

      int i;
      for(i = 0; i < numModes; ++i){
        (*modalOps.K)[i]      = freqs[i] * freqs[i];
        (*modalOps.dynMat)[i] = kcoef*(*modalOps.K)[i] + mcoef*(*modalOps.M)[i];
      }

      // damping matrix
      double alpha = domain->solInfo().alphaDamp;
      double beta  = domain->solInfo().betaDamp;

      BCond* damping;
      int numDampedModes = geoSource->getModalDamping(damping);

      if(damping){
        modalOps.C = new DiagonalMatrix(numModes);
        modalOps.C->setDiag(0.0);
        for(i = 0; i < numDampedModes; ++i)
          (*modalOps.C)[damping[i].nnum] += 2*damping[i].val*freqs[damping[i].nnum];
        for(i = 0; i < numModes; ++i)
          (*modalOps.dynMat)[i] += ccoef * (*modalOps.C)[i];
      }
      else if(alpha != 0.0 || beta != 0.0){
        modalOps.C = new DiagonalMatrix(numModes);
        for(i = 0; i < numModes; ++i){
          (*modalOps.C)[i] = alpha*(*modalOps.M)[i] + beta*(*modalOps.K)[i];
          (*modalOps.dynMat)[i] += ccoef * (*modalOps.C)[i];
        }
      }
      else{modalOps.C = 0;}

      modalOps.dynMat->invertDiag();
    } break;
    case ModalParams::Mnorm : { // if Mass normal basis, only Mass matrix is diagonal
#ifdef USE_EIGEN3
      // every process gets an identical copy of modalOps
      modalOps.M       = new DiagonalMatrix(numModes);
      modalOps.Msolver = new DiagonalMatrix(numModes);
      // see below for damping matrix
      modalOps.K       = new DenseMatrix(numModes);
      modalOps.dynMat  = new DenseMatrix(numModes);

      modalOps.M->setDiag(1.0);
      modalOps.Msolver->setDiag(1.0); // Inverse of M

      //Construct and assemble full Stiffness matrix in parallel
      SparseMatrix **subK = new SparseMatrix * [decDomain->getNumSub()];
      SparseMatrix **subM;
      if (checkBasis) subM = new SparseMatrix * [decDomain->getNumSub()];

      int sizeM = 0; 
      double normM = 0;
      for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {
        AllOps<double> allOps;
        if (checkBasis) allOps.M = subM[iSub] = decDomain->getSubDomain(iSub)->constructEiSparse<double>();
        allOps.K = subK[iSub] = decDomain->getSubDomain(iSub)->constructDBSparseMatrix<double>();
        decDomain->getSubDomain(iSub)->makeSparseOps<double>(allOps, 0, 0, 0);
        if (checkBasis) {
          sizeM += static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().size();
          normM += static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().norm();
        }
      }
      SubDOp K(decDomain->getNumSub(), subK);
     
      if(checkBasis) {
        filePrint(stderr," ... Validating orthogonality of modal basis ... \n");

        sizeM = structCom->globalSum(sizeM);
        if (sizeM == 0) {
          filePrint(stderr, " *** ERROR: Mass matrix could not be constructed. Exiting...\n\n");
          exit(-1);
        }
        normM = structCom->globalSum(normM);
        if (normM == 0) {
          filePrint(stderr, " ... WARNING: The Frobenius norm of the mass matrix is zero.\n");
        }

        DenseMatrix redMass(numModes);
        SubDOp M(decDomain->getNumSub(), subM);
        GenDistrVectorSet<double> tPhiM(numModes, decDomain->solVecInfo());
        for(int i = 0; i<numModes; ++i)
          M.mult((*modesFl)[i], tPhiM[i]);

        for(int row = 0; row < numModes; ++row) {
          for(int col = 0; col < numModes; ++col) {
            redMass[row+numModes*col] = tPhiM[row]*(*modesFl)[col];
            if(col == row) // subtract identity
              redMass[row+numModes*col] -= 1; 
          }
        }
    
        if(redMass.norm()/sqrt(numModes) > modalParams.tolerance) {
          filePrint(stderr," ... Modal Basis does not satisfy mass orthogonality tolerance, exiting. \n");
          exit(-1); 
        }
      }

      {
        GenDistrVectorSet<double> tPhiK(numModes, decDomain->solVecInfo());
        for(int i = 0; i<numModes; ++i)
          K.mult((*modesFl)[i], tPhiK[i]);

        for(int row = 0; row < numModes; ++row) {
          for(int col = 0; col < numModes; ++col) {
            (*modalOps.K)[row+numModes*col] = tPhiK[row]*(*modesFl)[col]; 
          }
        }
      } // end parallel section of reduced operator construction 
      
      // add mass and stiffness contribution to dynMat
      for(int row = 0; row < numModes; ++row){
        for(int col = 0; col < numModes; ++col){
          if(row == col)
            (*modalOps.dynMat)[row+numModes*col] += kcoef*(*modalOps.K)[row+numModes*col] + mcoef*(*modalOps.M)[col];
          else
            (*modalOps.dynMat)[row+numModes*col] += kcoef*(*modalOps.K)[row+numModes*col];
        }
      }

      // damping matrix
      double alpha = domain->solInfo().alphaDamp;
      double beta  = domain->solInfo().betaDamp;

      // only Rayleigh damping is supported currently
      if(alpha != 0.0 || beta != 0.0){ 
        modalOps.C = new DenseMatrix(numModes);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            if(row == col) { // only add diagonal Mass elements, off diagonal elements don't exist in data structure
              (*modalOps.C)[row+numModes*col] = alpha*(*modalOps.M)[col] + beta*(*modalOps.K)[row+numModes*col];
            } else {
              (*modalOps.C)[row+numModes*col] = beta*(*modalOps.K)[row+numModes*col];
            }
            // add damping matrix contribution to dynMat
            (*modalOps.dynMat)[row+numModes*col] += ccoef*(*modalOps.C)[row+numModes*col];
          }
        }
      }
      else{modalOps.C = 0;}
      modalOps.dynMat->invertDiag(); // this actually forms the LLT factorization for the DenseMatrix class
#endif    
    } break;
    case ModalParams::Inorm : { // if Identity normal basis, all operators are dense
#ifdef USE_EIGEN3
      // every process gets an identcial copy of modalOps
      modalOps.M       = new DenseMatrix(numModes);
      modalOps.Msolver = new DenseMatrix(numModes);
      // see below for damping matrix
      modalOps.K       = new DenseMatrix(numModes);
      modalOps.dynMat  = new DenseMatrix(numModes);

      if(checkBasis) {
        filePrint(stderr," ... Validating orthogonality of modal basis ... \n");
        DenseMatrix redIdentity(numModes);
        GenDistrVectorSet<double> tPhi(numModes, decDomain->solVecInfo());
        for(int i = 0; i<numModes; ++i)
          tPhi[i] = (*modesFl)[i];

        for(int row = 0; row < numModes; ++row) {
          for(int col = 0; col < numModes; ++col) {
            redIdentity[row+numModes*col] = tPhi[row]*(*modesFl)[col];
            if(col == row) // subtract identity
              redIdentity[row+numModes*col] -= 1.0;
          }
        }

        if(redIdentity.norm()/sqrt(numModes) > modalParams.tolerance) {
          filePrint(stderr," ... Modal Basis does not satisfy self orthogonality tolerance, exiting. \n");
          exit(-1);
        }
      }

      //Construct and assemble full mass and Stiffness matrix
      SparseMatrix **subM = new SparseMatrix * [decDomain->getNumSub()];
      SparseMatrix **subK = new SparseMatrix * [decDomain->getNumSub()];
      for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {
        AllOps<double> allOps;
        allOps.K = subK[iSub] = decDomain->getSubDomain(iSub)->constructDBSparseMatrix<double>();
        allOps.M = subM[iSub] = decDomain->getSubDomain(iSub)->constructDBSparseMatrix<double>();
        decDomain->getSubDomain(iSub)->makeSparseOps<double>(allOps, 0, 0, 0);
      }
      SubDOp K(decDomain->getNumSub(), subK);
      SubDOp M(decDomain->getNumSub(), subM);

      {
        GenDistrVectorSet<double> tPhiM(numModes, decDomain->solVecInfo());
        GenDistrVectorSet<double> tPhiK(numModes, decDomain->solVecInfo());
        for(int i = 0; i<numModes; ++i){
          M.mult((*modesFl)[i], tPhiM[i]);
          K.mult((*modesFl)[i], tPhiK[i]);
        }

        for(int row = 0; row < numModes; ++row) {
          for(int col = 0; col < numModes; ++col) {
            (*modalOps.M)[row+numModes*col] = tPhiM[row]*(*modesFl)[col];
            (*modalOps.K)[row+numModes*col] = tPhiK[row]*(*modesFl)[col];
          }
        }
      } // end parallel section of reduced operator construction 

      // add mass and stiffness contribution to dynMat
      for(int row = 0; row < numModes; ++row){
        for(int col = 0; col < numModes; ++col){
          (*modalOps.dynMat)[row+numModes*col] += kcoef*(*modalOps.K)[row+numModes*col] + mcoef*(*modalOps.M)[row+numModes*col];
        }
      }
    
      //damping matrix paraphernalia  
      double alpha = domain->solInfo().alphaDamp;
      double beta  = domain->solInfo().betaDamp;

      // only Rayleigh damping is supported currently
      if(alpha != 0.0 || beta != 0.0){
        modalOps.C = new DenseMatrix(numModes);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            if(row == col) { // only add diagonal Mass elements, off diagonal elements don't exist in data structure
              (*modalOps.C)[row+numModes*col] = alpha*(*modalOps.M)[row+numModes*col] + beta*(*modalOps.K)[row+numModes*col];
            } else {
              (*modalOps.C)[row+numModes*col] = beta*(*modalOps.K)[row+numModes*col];
            }
            // add damping matrix contribution to dynMat
            (*modalOps.dynMat)[row+numModes*col] += ccoef*(*modalOps.C)[row+numModes*col];
          }
        }
      }
      else{modalOps.C = 0;}
      modalOps.dynMat->invertDiag(); // this actually forms the LLT factorization for the DenseMatrix class
#endif
    }
    default : {
      filePrint(stderr," *** ERROR: unsupported modal basis type specified under READMODE command. \n");
      exit(-1);
    }
  }
  
  return (&modalOps);
}

//------------------------------------------------------------------------------

void MultiDomainModal::getQuasiStaticParameters(double &maxVel, double &delta){
  maxVel = domain->solInfo().qsMaxvel;
  delta  = domain->solInfo().delta;
}

//------------------------------------------------------------------------------

void MultiDomainModal::computeExtForce2(SysState<Vector>& state, Vector &extF,
  Vector &constF, int tIndex, double time, Vector *aeroF,
  double gamma, double alphaf){
/*PRE:
 POST: return in extF, the modalized external force
*/
  MultiDomainOp mdop(&MultiDomainOp::computeExtForce,
                     decDomain->getAllSubDomains(), fullTmpF, fullTmpGrav, time, (SubDOp*)0, (ControlInterface*)0, (SubDOp*)0, 0, (SubDOp*)0);
  threadManager->execParal(decDomain->getNumSub(), &mdop);

  if(domain->solInfo().aeroFlag >= 0 && tIndex >= 0) {
    //domain->buildAeroelasticForce(*aeroF, *prevFrc, tIndex, time, gamma, alphaf);
    //(*fullTmpF) += *aeroF;
    filePrint(stderr, " *** ERROR: MultiDomainModal::computeExtForce2 is not implemented for aeroFlag = %d.\n", domain->solInfo().aeroFlag);
    exit(-1);
  }

  projectForce(*fullTmpF, extF);
}

//------------------------------------------------------------------------------

void MultiDomainModal::getInternalForce(Vector &d, Vector &f, double t, int tIndex){
/*PRE: d is the value of the modal coordinates
 POST: return the modal internal force in f
*/

  f.zero();
  modalOps.K->mult(d,f);
}

//------------------------------------------------------------------------------

double MultiDomainModal::getElasticEnergy(Vector &d){
/*PRE: d is the value of the modal coordinates
 POST: return the elastic energy
*/
 Vector Kd(d.size());
 Kd.zero();
 double Wela = 0.0;
 modalOps.K->mult(d,Kd);
 for(int i = 0; i < d.size(); ++i)
    Wela += 0.5*d[i]*Kd[i];
 return Wela;
}

//------------------------------------------------------------------------------

double MultiDomainModal::getKineticEnergy(Vector &v){
/*PRE: v is the value of the first time derivative of the modal coordinates
 POST: return the kinetic energy
*/
  Vector Mv(v.size());
  Mv.zero();
  double Wkin = 0.0;
  modalOps.M->mult(v,Mv);
  for(int i = 0; i < v.size(); ++i)
    Wkin += 0.5*v[i]*Mv[i];

  return Wkin;
}

//------------------------------------------------------------------------------

double MultiDomainModal::getDampingEnergy(Vector &d, Vector &v, double time){
/*PRE: d is the value of the modal coordinates and v is its first time derivative
 POST: return the damping energy
*/
  Vector tmpVec(v.size());

  if(time == domain->solInfo().initialTime) {
    Wdmp  = 0.0;
    if(modalOps.C) {
      modalOps.C->mult(v, tmpVec);
      previousCq = new Vector(tmpVec);
      previousDisp = new Vector(d);
    }
  }
  else {
    if(modalOps.C) {
      double c = domain->solInfo().newmarkGamma;
      modalOps.C->mult(v, tmpVec);
      Wdmp += (c*tmpVec + (1.0-c)*(*previousCq))*(d - (*previousDisp));
      (*previousCq) = tmpVec;
      (*previousDisp) = d;
    }
  }
  return Wdmp;
}

//------------------------------------------------------------------------------

void MultiDomainModal::dynamOutput(int tIndex, double time, ModalOps &ops, Vector &extF,
  Vector *aeroF, SysState<Vector> &state){

  expand(state.getDisp(), *fullDsp);
  expand(state.getVeloc(), *fullVel);
  expand(state.getAccel(), *fullAcc);
  expand(state.getPrevVeloc(), *fullPrevVel);

  MDDynamMat dumDMat; 
  domain->setModalEnergies(getElasticEnergy(state.getDisp()),
                           getKineticEnergy(state.getVeloc()),
                           getDampingEnergy(state.getDisp(),state.getVeloc(),time));
  SysState<DistrVector> distState(*fullDsp, *fullVel, *fullAcc, *fullPrevVel);
  decDomain->postProcessing(*fullDsp, *fullTmpF, time, fullAeroF, tIndex, &dumDMat, &distState);
  outputModal(state, extF, tIndex, ops);
}

//------------------------------------------------------------------------------

int MultiDomainModal::aeroPreProcess(Vector& d_n, Vector& v_n,
  Vector& a_n, Vector& v_p){
/*PRE: arguments are in modal space
 POST: call domain aeroPreProcess with the expanded state vectors
*/
  expand(d_n, *fullDsp);
  expand(v_n, *fullVel);
  expand(a_n, *fullAcc);
  expand(v_p, *fullPrevVel);

  //domain->aeroPreProcess(fullDsp, fullVel, fullAcc, fullPrevVel, bcx, vcx);
  filePrint(stderr, " *** ERROR: MultiDomainModal::aeroPreProcess is not implemented.\n");
  exit(-1);

  return domain->solInfo().aeroFlag;
}

//------------------------------------------------------------------------------

int MultiDomainModal::aeroSensitivityPreProcess(Vector& d_n, Vector& v_n,
  Vector& a_n, Vector& v_p){
/*PRE: arguments are in modal space
 POST: call domain aeroSensitivityPreProcess with the expanded state vectors
*/
  expand(d_n, *fullDsp);
  expand(v_n, *fullVel);
  expand(a_n, *fullAcc);
  expand(v_p, *fullPrevVel);

  //domain->aeroSensitivityPreProcess(fullDsp, fullVel, fullAcc, fullPrevVel, bcx, vcx);
  filePrint(stderr, " *** ERROR: MultiDomainModal::aeroSensitivityPreProcess is not implemented.\n");
  exit(-1);

  return domain->solInfo().aeroFlag;
}

//------------------------------------------------------------------------------

int MultiDomainModal::sendDisplacements(Vector& d_n, Vector& v_n,
                                        Vector& a_n, Vector& v_p){
/*PRE: arguments are in modal space
 POST: call domain sendDisplacements with the expanded state vectors
*/
  expand(d_n, *fullDsp);
  expand(v_n, *fullVel);
  expand(a_n, *fullAcc);
  expand(v_p, *fullPrevVel);

  //domain->sendDisplacements(fullDsp, fullVel, fullAcc, fullPrevVel, bcx, vcx);
  filePrint(stderr, " *** ERROR: MultiDomainModal::sendDisplacements is not implemented.\n");
  exit(-1);

  return domain->solInfo().aeroFlag;
}

//------------------------------------------------------------------------------

void MultiDomainModal::a5TimeLoopCheck(int& parity, double& t, double dt){
/*NOTE: this function copied from SingleDomainDynamic::a5TimeLoopCheck
*/
  if (domain->solInfo().aeroFlag == 5) {
    if (!parity) t -= dt;
    parity = ( parity ? 0 : 1 );
  }
}

//------------------------------------------------------------------------------

void MultiDomainModal::a5StatusRevise(int parity, SysState<Vector>& curState,
  SysState<Vector>& bkState){
/*NOTE: this function copied from SingleDomainDynamic::a5StatusRevise
*/
  if (domain->solInfo().aeroFlag == 5) {
    if (parity) { // restore

      *prevFrc = *prevFrcBackup;
      curState.getDisp()      = bkState.getDisp();
      curState.getVeloc()     = bkState.getVeloc();
      curState.getAccel()     = bkState.getAccel();
      curState.getPrevVeloc() = bkState.getPrevVeloc();

    } else {      // backup

      *prevFrcBackup = *prevFrc;
      bkState.getDisp()      = curState.getDisp();
      bkState.getVeloc()     = curState.getVeloc();
      bkState.getAccel()     = curState.getAccel();
      bkState.getPrevVeloc() = curState.getPrevVeloc();

    }
  }
}

//------------------------------------------------------------------------------

AllSensitivities<double> * MultiDomainModal::getAllSensitivities(){

  filePrint(stderr," *** ERROR: MultiDomainModal::getAllSensitivities is not implemented.\n");
  exit(-1);
}

//------------------------------------------------------------------------------

void MultiDomainModal::getAeroelasticForceSensitivity(int t_index, double t, Vector * aero_f, double gamma, double alphaf){

  filePrint(stderr," *** ERROR: MultiDomainModal::getAeroelasticForceSensitivity is not implemented.\n");
  exit(-1);
}

//------------------------------------------------------------------------------

void MultiDomainModal::computeStabilityTimeStep(double &dt, ModalOps &dMat){

 // ... Compute Stability Time Step
 double sts = domain->computeStabilityTimeStep(dMat);

 if(sts == std::numeric_limits<double>::infinity()) {
   filePrint(stderr," **************************************\n");
   filePrint(stderr," Stability max. timestep could not be  \n");
   filePrint(stderr," determined for this model.            \n");
   filePrint(stderr," Specified time step is selected\n");
   filePrint(stderr," **************************************\n");
   domain->solInfo().stable = 0;
 }
 else {
   filePrint(stderr," **************************************\n");
   if (domain->solInfo().modifiedWaveEquation) {
     sts = 1.73205*sts;
     filePrint(stderr," CONDITIONALLY STABLE MODIFIED WAVE EQUATION\n");
   }
   else
     filePrint(stderr," CONDITIONALLY STABLE NEWMARK ALGORITHM\n");

   filePrint(stderr," --------------------------------------\n");
   filePrint(stderr," Specified time step      = %10.4e\n",dt);
   filePrint(stderr," Stability max. time step = %10.4e\n",sts);
   filePrint(stderr," **************************************\n");
   if( (domain->solInfo().stable == 1 && sts < dt) || domain->solInfo().stable == 2 ) {
     dt = sts;
     filePrint(stderr," Stability max. time step is selected\n");
   } else
     filePrint(stderr," Specified time step is selected\n");
   filePrint(stderr," **************************************\n");
 }

 domain->solInfo().setTimeStep(dt);
}

