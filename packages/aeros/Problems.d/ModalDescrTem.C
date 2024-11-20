#include <Problems.d/ModalDescr.h>
#include <Driver.d/Domain.h>
#include <Driver.d/SysState.h>
#include <Driver.d/Dynam.h>
#include <Math.d/EiSparseMatrix.h>

template <class Scalar>
ModalDescr<Scalar>::ModalDescr(Domain *d) : ModalBase(d), modalOps(*(new ModalOps)){

  flExchanger = domain->getFileExchanger();
  previousCq = 0;
  previousDisp = 0;
}

//------------------------------------------------------------------------------

template <class Scalar>
ModalDescr<Scalar>::~ModalDescr(){

  if(previousCq) delete previousCq;
  if(previousDisp) delete previousDisp;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::projectForce(Vector &fullF, Vector& modalF){
/*PRE: ModeBase::modesFl have been populated by a call to
         ModeBase::populateFlexModes
 POST: return in modalF, fullF projected onto modesFl
 NOTE: this projection is not valid for displacements, velocities or accelerations
*/
  for(int i = 0; i < numModes; ++i)
    modalF[i] = modesFl[i] * fullF;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::expand(const Vector &modalV, Vector& fullV){
/*PRE: there is at least one modesFl; modesFl are populated
 POST: return in fullV, modalV projected into full space
 NOTE: this expansion is not valid for forces
*/
  fullV.linC(modalV[0], modesFl[0]);
  for(int i = 1; i < numModes; ++i)
    fullV.linAdd(modalV[i], modesFl[i]);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::processLastOutput()  {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
  
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::preProcess(){
/*NOTE: call to populateFlexModes gives 1 for 2nd arguement to indicate all
          modes should be read, including those with zero frequency
        see also ModalBase::populateFlexModes
*/
  preProcessBase();
  populateFlexModes(1.0, 1);
  numModes = numFlex;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::preProcessSA(){

  filePrint(stderr," ... ModalDescr::preProcessSA is not implemented\n");
  exit(-1);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::postProcessSA(ModalOps *,Vector &sol){

  filePrint(stderr," ... ModalDescr::postProcessSA is not implemented\n"); 
  exit(-1);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::sensitivityPostProcessing(Vector *sol){

  filePrint(stderr," ... ModalDescr::sensitivityPostProcessing is not implemented\n"); 
  exit(-1);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getTimes(double &dt, double &tmax){

  dt   = domain->solInfo().getTimeStep();
  tmax = domain->solInfo().tmax;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getInitState(SysState<Vector> &state){

  initStateBase(state.getDisp(), state.getVeloc(),
    state.getAccel(), state.getPrevVeloc());

/*
  state.getDisp().zero();
  state.getVeloc().zero();
  state.getAccel().zero();
  state.getPrevVeloc().zero();

  // values for state obtained directly from domain
  for(int jj = 0; jj <  domain->numInitDisp(); ++jj)
    state.getDisp()[domain->getInitDisp()[jj].nnum] +=
      domain->getInitDisp()[jj].val;
*/
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getConstForce(Vector &constF){

  domain->computeConstantForce(fullTmpGrav); 
/*
  fullTmpGrav.zero();

  if( domain->gravityFlag()  ) domain->buildGravityForce(fullTmpGrav);
  if( domain->pressureFlag() ) domain->buildPressureForce(fullTmpGrav);
*/
  projectForce(fullTmpGrav, constF);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::addConstForceSensitivity(Vector &constF){

  filePrint(stderr," ... ModalDescr<Scalar>::addConstForceSensitivity is not implemented\n"); 
  exit(-1); 
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getSensitivityStateParam(double &tol, double &ratioTol) {
  tol = domain->solInfo().sensitivityTol;
  ratioTol = domain->solInfo().ratioSensitivityTol;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getSteadyStateParam(int &flag, int &min, int &max, double &tol){

  flag = domain->solInfo().steadyFlag;
  min  = domain->solInfo().steadyMin;
  max  = domain->solInfo().steadyMax;
  tol  = domain->solInfo().steadyTol;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getNewMarkParameters(double &beta, double &gamma,
  double &alphaf, double &alpham){

  beta  = domain->solInfo().newmarkBeta;
  gamma = domain->solInfo().newmarkGamma;
  alphaf = domain->solInfo().newmarkAlphaF;
  alpham = domain->solInfo().newmarkAlphaM;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getInitialTime(int &tIndex, double &time){

  tIndex = domain->solInfo().initialTimeIndex;
  time   = domain->solInfo().initialTime;
}

//------------------------------------------------------------------------------

template <class Scalar>
double
ModalDescr<Scalar>::getInitialForceNorm()
{
  return domain->solInfo().initExtForceNorm;
}

//------------------------------------------------------------------------------

template <class Scalar>
ModalOps* ModalDescr<Scalar>::buildOps(double mcoef, double ccoef, double kcoef){
//PRE: ModalProbDescr is instantiated; modeData is populated
// POST: operator for time integration
//  need to check what time of modal basis has been parsed (i.e. Eigen, Mnorm, Inorm)

  ModalParams modalParams = domain->solInfo().readInModes[domain->solInfo().modal_id.front()]; // assume only one basis is parsed for modal analysis
  
  // check to see if basis satisfies orthogonality property
  bool checkBasis = false; 
  if(modalParams.tolerance != 0)
    checkBasis = true; 
 

  switch(modalParams.type) { 
    case ModalParams::Eigen :
    case ModalParams::Undefined : { // if eigen basis, all operators are diagonal 
      modalOps.M       = new DiagonalMatrix(numModes);
      modalOps.Msolver = new DiagonalMatrix(numModes);
      // see below for damping matrix
      modalOps.K       = new DiagonalMatrix(numModes);
      modalOps.dynMat  = new DiagonalMatrix(numModes);

      modalOps.M->setDiag(1.0);
      modalOps.Msolver->setDiag(1.0); // Inverse of M

      if (checkBasis) {
        AllOps<double> allOps;
#if defined(USE_EIGEN3)
        allOps.M = domain->constructEiSparse<double>();
        allOps.K = domain->constructEiSparse<double>();
        domain->makeSparseOps(allOps, 0, 0, 0, (SparseMatrix *) NULL, (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL);
        if (static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().size() == 0 ||
            static_cast<EiSparseMatrix*>(allOps.K)->getEigenSparse().size() == 0) {
          std::cerr << " *** ERROR: Mass matrix and/or stiffness matrix could not be constructed. Exiting...\n\n";
          exit(-1);
        }
        if (static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().norm() == 0) {
          std::cerr << " ... WARNING: The Frobenius norm of the mass matrix is zero.\n";
        }
        if (static_cast<EiSparseMatrix*>(allOps.K)->getEigenSparse().norm() == 0) {
          std::cerr << " ... WARNING: The Frobenius norm of the stiffness matrix is zero.\n";
        }
#else
        allOps.M = domain->constructDBSparseMatrix<double>();
        allOps.K = domain->constructDBSparseMatrix<double>();
        domain->makeSparseOps(allOps, 0, 0, 0, (SparseMatrix *) NULL, (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL);
#endif
        fprintf(stderr," ... Validating orthogonality of modal eigen basis ... \n"); 
        double **tPhiM = new double*[numModes];
        double **tPhiK = new double*[numModes];
        for(int i = 0; i < numModes; ++i){
          tPhiM[i] = new double[domain->numdof()];
          tPhiK[i] = new double[domain->numdof()];
        }

        for(int i = 0 ; i<numModes; ++i){
          allOps.M->mult(modesFl[i].data(), tPhiM[i]);
          allOps.K->mult(modesFl[i].data(), tPhiK[i]);
        }

#ifdef USE_EIGEN3
        DenseMatrix redMass(numModes);    
        DenseMatrix redStif(numModes);
        double Knorm = 0.0; 
        Vector fullDsp(domain->numdof(), 0.0);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            double matEle1 = 0.0; 
            double matEle2 = 0.0;
            for(int dof = 0; dof < fullDsp.size(); ++dof){
              matEle1 += tPhiM[row][dof]*modesFl[col][dof];
              matEle2 += tPhiK[row][dof]*modesFl[col][dof]; 
            }
            redMass[row+numModes*col] = matEle1;
            redStif[row+numModes*col] = matEle2;
            if (col == row) {// subtract identity and circular frequencies 
              redMass[row+numModes*col] -= 1;
              double Kdiag = freqs[row]*freqs[row];
              redStif[row+numModes*col] -= Kdiag; 
              Knorm += Kdiag*Kdiag; 
            }
          }
        }
        Knorm = sqrt(Knorm); 

        if (redMass.norm()/sqrt(numModes) > modalParams.tolerance) {
          fprintf(stderr," ... Modal basis does not satisfy mass matrix orthogonality tolerance, exiting. \n");
          exit(-1);
        }
        if (redStif.norm()/Knorm > modalParams.tolerance) {
          fprintf(stderr," ... Modal basis does not satisfy stiffness decoupling tolerance, exiting. \n");
          exit(-1);
        }
#endif
        for(int i = 0 ; i<numModes; ++i){
          delete [] tPhiM[i];
          delete [] tPhiK[i];
        }
        delete [] tPhiM;
        delete [] tPhiK;
        delete allOps.M;
        delete allOps.K;
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
      modalOps.M       = new DiagonalMatrix(numModes);
      modalOps.Msolver = new DiagonalMatrix(numModes);
      // see below for damping matrix
      modalOps.K       = new DenseMatrix(numModes);
      modalOps.dynMat  = new DenseMatrix(numModes);

      modalOps.M->setDiag(1.0);
      modalOps.Msolver->setDiag(1.0); // Inverse of M

      //Construct and assemble full Stiffness matrix
      AllOps<double> allOps;
      if (checkBasis) allOps.M = domain->constructEiSparse<double>();
      allOps.K = domain->constructDBSparseMatrix<double>();
      domain->makeSparseOps(allOps, 0, 0, 0, (SparseMatrix *) NULL, (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL);
      if (checkBasis) {
        fprintf(stderr," ... Validating orthogonality of modal basis. ... \n");

        if (static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().size() == 0) {
          std::cerr << " *** ERROR: Mass matrix could not be constructed. Exiting...\n\n";
          exit(-1);
        }
        if (static_cast<EiSparseMatrix*>(allOps.M)->getEigenSparse().norm() == 0) {
          std::cerr << " ... WARNING: The Frobenius norm of the mass matrix is zero.\n";
        }
     
        double **tPhiM = new double*[numModes];
        for(int i = 0; i < numModes; ++i)
          tPhiM[i] = new double[domain->numdof()];

        for(int i = 0 ; i<numModes; ++i)
          allOps.M->mult(modesFl[i].data(), tPhiM[i]);

        DenseMatrix redMass(numModes); 
        Vector fullDsp(domain->numdof(), 0.0);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            double matEle = 0.0;
            for(int dof = 0; dof < fullDsp.size(); ++dof)
              matEle += tPhiM[row][dof]*modesFl[col][dof];
            redMass[row+numModes*col]       = matEle;
            if (col == row) // subtract identity 
              redMass[row+numModes*col] -= 1; 
          }
        }
        for(int i = 0 ; i<numModes; ++i)
          delete [] tPhiM[i];
        delete [] tPhiM;

        if (redMass.norm()/sqrt(numModes) > modalParams.tolerance) {
          fprintf(stderr," ... Modal basis does not satisfy mass matrix orthogonality tolerance, exiting. \n");
          exit(-1);   
        }
 
      }
 
      {
        // allocate space for intermediate container
        double **tPhiK = new double*[numModes];
        for(int i = 0; i < numModes; ++i)
          tPhiK[i] = new double[domain->numdof()];

        // taking advantage of symmetry of K and computing K*Phi_i instead of transpose(Phi_i)*K
        for(int i = 0 ; i<numModes; ++i)
          allOps.K->mult(modesFl[i].data(), tPhiK[i]);
  
        // now store transpose(Phi_i)*K*Phi_i in dense matrix container modalOps.K
        Vector fullDsp(domain->numdof(), 0.0);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            double matEle = 0.0;
            for(int dof = 0; dof < fullDsp.size(); ++dof)
              matEle += tPhiK[row][dof]*modesFl[col][dof];
            (*modalOps.K)[row+numModes*col] = matEle;
          }
        }
        for(int i = 0 ; i<numModes; ++i)
          delete [] tPhiK[i]; 
        delete [] tPhiK; 
      }
      delete allOps.K;
      if (checkBasis) delete allOps.M; 

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
      modalOps.M       = new DenseMatrix(numModes);
      modalOps.Msolver = new DenseMatrix(numModes);
      // see below for damping matrix
      modalOps.K       = new DenseMatrix(numModes);
      modalOps.dynMat  = new DenseMatrix(numModes);

       if (checkBasis) {
        fprintf(stderr," ... Validating orthogonality of modal basis. ... \n");
        double **tPhi = new double*[numModes];
        for(int i = 0; i < numModes; ++i)
          tPhi[i] = new double[domain->numdof()];

        for(int i = 0 ; i<numModes; ++i)
          tPhi[i] = modesFl[i].data();

        DenseMatrix redIdentity(numModes);
        Vector fullDsp(domain->numdof(), 0.0);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            double matEle = 0.0;
            for(int dof = 0; dof < fullDsp.size(); ++dof)
              matEle += tPhi[row][dof]*modesFl[col][dof];
            redIdentity[row+numModes*col]       = matEle;
            if (col == row) // subtract identity 
              redIdentity[row+numModes*col] -= 1;
          }
        }
        for(int i = 0 ; i<numModes; ++i)
          delete [] tPhi[i];
        delete [] tPhi;

        if (redIdentity.norm()/sqrt(numModes) > modalParams.tolerance) {
          fprintf(stderr," ... Modal basis does not satisfy self orthogonality tolerance, exiting. \n");
          exit(-1);
        }

      }

      //Construct and assemble full Stiffness matrix
      AllOps<double> allOps;
      allOps.M = domain->constructDBSparseMatrix<double>();
      allOps.K = domain->constructDBSparseMatrix<double>();
      domain->makeSparseOps(allOps, 0, 0, 0, (SparseMatrix *) NULL, (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL);

      { // construction of reduced mass matrix
        // allocate space for intermediate container
        double **tPhiM = new double*[numModes];
        for(int i = 0; i < numModes; ++i)
          tPhiM[i] = new double[domain->numdof()];

        // taking advantage of symmetry of M and computing M*Phi_i instead of transpose(Phi_i)*M
        for(int i = 0 ; i<numModes; ++i)
          allOps.M->mult(modesFl[i].data(), tPhiM[i]);

        // now store transpose(Phi_i)*M*Phi_i in dense matrix container modalOps.M
        Vector fullDsp(domain->numdof(), 0.0);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            double matEle = 0.0;
            for(int dof = 0; dof < fullDsp.size(); ++dof)
              matEle += tPhiM[row][dof]*modesFl[col][dof];
            (*modalOps.M)[row+numModes*col]       = matEle;
            (*modalOps.Msolver)[row+numModes*col] = matEle; 
          }
        }
        for(int i = 0 ; i<numModes; ++i)
          delete [] tPhiM[i];
        delete [] tPhiM;
      }

      modalOps.Msolver->invertDiag(); // compute LLT factorization

      { // construction of reduced stiffness matrix
        // allocate space for intermediate container
        double **tPhiK = new double*[numModes];
        for(int i = 0; i < numModes; ++i)
          tPhiK[i] = new double[domain->numdof()];

        // taking advantage of symmetry of K and computing K*Phi_i instead of transpose(Phi_i)*K
        for(int i = 0 ; i<numModes; ++i)
          allOps.K->mult(modesFl[i].data(), tPhiK[i]);
  
        // now store transpose(Phi_i)*K*Phi_i in dense matrix container modalOps.K
        Vector fullDsp(domain->numdof(), 0.0);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            double matEle = 0.0;
            for(int dof = 0; dof < fullDsp.size(); ++dof)
              matEle += tPhiK[row][dof]*modesFl[col][dof];
            (*modalOps.K)[row+numModes*col] = matEle;
          }
        }
        for(int i = 0 ; i<numModes; ++i)
          delete [] tPhiK[i];
        delete [] tPhiK;
      }
      delete allOps.K;
      delete allOps.M;

      // add mass and stiffness contribution to dynMat
      for(int row = 0; row < numModes; ++row){
        for(int col = 0; col < numModes; ++col){
          (*modalOps.dynMat)[row+numModes*col] += kcoef*(*modalOps.K)[row+numModes*col] + mcoef*(*modalOps.M)[row+numModes*col];
        }
      }
 
      // damping matrix
      double alpha = domain->solInfo().alphaDamp;
      double beta  = domain->solInfo().betaDamp;

      // only reighley damping is supported currently
      if(alpha != 0.0 || beta != 0.0){
        modalOps.C = new DenseMatrix(numModes);
        for(int row = 0; row < numModes; ++row){
          for(int col = 0; col < numModes; ++col){
            (*modalOps.C)[row+numModes*col] = beta*(*modalOps.K)[row+numModes*col];
            // add damping matrix contribution to dynMat
            (*modalOps.dynMat)[row+numModes*col] += ccoef*(*modalOps.C)[row+numModes*col];
          }
        }
      } 
      else{modalOps.C = 0;}
      modalOps.dynMat->invertDiag(); // this actually forms the LLT factorization for the DenseMatrix class
      
#endif
    } break;
    default : {
       fprintf(stderr," *** ERROR: unsupported modal basis type specified under READMODE command. \n");
       exit(-1);
    }
  }
  return (&modalOps);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getQuasiStaticParameters(double &maxVel, double &delta){
  maxVel = domain->solInfo().qsMaxvel;
  delta  = domain->solInfo().delta;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::computeExtForce2(SysState<Vector>& state, Vector &extF,
  Vector &constF, int tIndex, double time, Vector *aeroF,
  double gamma, double alphaf){
/*PRE:
 POST: return in extF, the modalized external force
*/
  domain->template computeExtForce4<double>(fullTmpF, fullTmpGrav, time, 0);
  if(domain->solInfo().aeroFlag >= 0 && tIndex >= 0) {
    domain->buildAeroelasticForce(*aeroF, *prevFrc, tIndex, time, gamma, alphaf);
    fullTmpF += *aeroF;
  }
  projectForce(fullTmpF, extF);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getInternalForce(Vector &d, Vector &f, double t, int tIndex){
/*PRE: d is the value of the modal coordinates
 POST: return the modal internal force in f
  Stiffness matrix may not be diagonal so this code is no longer valid
  for(int i = 0; i < d.size(); ++i)
    f[i] = d[i]*freqs[i]*freqs[i];
*/
  f.zero();
  modalOps.K->mult(d,f); 
}

//------------------------------------------------------------------------------

template <class Scalar>
double ModalDescr<Scalar>::getElasticEnergy(Vector &d){
/*PRE: d is the value of the modal coordinates
 POST: return the elastic energy
  Stiffness matrix may not be diagonal so this code is no longer valid
  double Wela = 0;
  for(int i = 0; i < d.size(); ++i)
    Wela += 0.5*d[i]*freqs[i]*freqs[i]*d[i];
  return Wela;
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

template <class Scalar>
double ModalDescr<Scalar>::getKineticEnergy(Vector &v){
/*PRE: v is the value of the first time derivative of the modal coordinates
 POST: return the kinetic energy

  double Wkin = 0;
  for(int i = 0; i < v.size(); ++i)
    Wkin += 0.5*v[i]*v[i];
*/ 
  // mass matrix may not be identity
  Vector Mv(v.size());
  Mv.zero();
  double Wkin = 0.0;
  modalOps.M->mult(v,Mv);
  for(int i = 0; i < v.size(); ++i)
    Wkin += 0.5*v[i]*Mv[i];

  return Wkin;
}

//------------------------------------------------------------------------------

template <class Scalar>
double ModalDescr<Scalar>::getDampingEnergy(Vector &d, Vector &v, double time){
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

template <class Scalar>
void ModalDescr<Scalar>::dynamOutput(int tIndex, double time, ModalOps &ops, Vector &extF,
  Vector *aeroF, SysState<Vector> &state){

  expand(state.getDisp(), fullDsp);
  expand(state.getVeloc(), fullVel);
  expand(state.getAccel(), fullAcc);
  expand(state.getPrevVeloc(), fullPrevVel);

  DynamMat dumDMat;
  domain->setModalEnergies(getElasticEnergy(state.getDisp()),
                           getKineticEnergy(state.getVeloc()),
                           getDampingEnergy(state.getDisp(),state.getVeloc(),time));
  domain->dynamOutput(tIndex, time, bcx, dumDMat, fullTmpF, fullAeroF,
    fullDsp, fullVel, fullAcc, fullPrevVel, vcx);
  outputModal(state, extF, tIndex, ops);
}

//------------------------------------------------------------------------------

template <class Scalar>
int ModalDescr<Scalar>::aeroPreProcess(Vector& d_n, Vector& v_n,
  Vector& a_n, Vector& v_p){
/*PRE: arguements are in modal space
 POST: call domain aeroPreProcess with the expanded state vectors
*/
  expand(d_n, fullDsp);
  expand(v_n, fullVel);
  expand(a_n, fullAcc);
  expand(v_p, fullPrevVel);

  domain->aeroPreProcess(fullDsp, fullVel, fullAcc, fullPrevVel, bcx, vcx);

  return domain->solInfo().aeroFlag;
}

//------------------------------------------------------------------------------

template <class Scalar>
int ModalDescr<Scalar>::aeroSensitivityPreProcess(Vector& d_n, Vector& v_n,
  Vector& a_n, Vector& v_p){
/*PRE: arguements are in modal space
 POST: call domain aeroSensitivityPreProcess with the expanded state vectors
*/
  expand(d_n, fullDsp);
  expand(v_n, fullVel);
  expand(a_n, fullAcc);
  expand(v_p, fullPrevVel);

  domain->aeroSensitivityPreProcess(fullDsp, fullVel, fullAcc, fullPrevVel, bcx, vcx);

  return domain->solInfo().aeroFlag;
}

//------------------------------------------------------------------------------

template <class Scalar>
int ModalDescr<Scalar>::sendDisplacements(Vector& d_n, Vector& v_n,
                                          Vector& a_n, Vector& v_p){
/*PRE: arguements are in modal space
 POST: call domain sendDisplacements with the expanded state vectors
*/
  expand(d_n, fullDsp);
  expand(v_n, fullVel);
  expand(a_n, fullAcc);
  expand(v_p, fullPrevVel);

  domain->sendDisplacements(fullDsp, fullVel, fullAcc, fullPrevVel, bcx, vcx);
  return domain->solInfo().aeroFlag;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::a5TimeLoopCheck(int& parity, double& t, double dt){
/*NOTE: this function copied from SingleDomainDynamic::a5TimeLoopCheck
*/
  if (domain->solInfo().aeroFlag == 5) {
    if (!parity) t -= dt;
    parity = ( parity ? 0 : 1 );
  }
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::a5StatusRevise(int parity, SysState<Vector>& curState,
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

template <class Scalar>
AllSensitivities<Scalar> * ModalDescr<Scalar>::getAllSensitivities(){

  filePrint(stderr," ... ModalDescr::getAllSensitivities is not implemented\n"); exit(-1);
  return 0;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getAeroelasticForceSensitivity(int t_index, double t, Vector * aero_f, double gamma, double alphaf){

  filePrint(stderr," ... ModalDescr::getAeroelasticForceSensitivity is not implemented\n"); exit(-1);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::computeStabilityTimeStep(double &dt, ModalOps &dMat){

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

