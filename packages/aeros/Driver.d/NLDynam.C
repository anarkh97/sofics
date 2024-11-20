#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Corotational.d/Corotator.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>

#include <Driver.d/GeoSource.h>
#include <Corotational.d/MatNLCorotator.h>
#ifdef USE_EIGEN3
#include <Element.d/Function.d/InertialForce.d/InertialType1ForceFunction.h>
#include <Element.d/Function.d/InertialForce.d/InertialType2ForceFunction.h>
#include <Element.d/Function.d/Rotation.d/IncrementalRotationVector.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#endif
#include <algorithm>

void
Domain::getElemInternalForce(const GeomState &geomState, double time,
                             const GeomState *refState, const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
    const_cast<Corotator &>(elemCorot).getInternalForce(
        const_cast<GeomState *>(refState),
        const_cast<GeomState &>(geomState),
        nodes, elemStiff, elemForce, domain->solInfo().getTimeStep(), time);
}

void
Domain::getElemInternalForce(const GeomState &geomState, double time,
                             const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
  const_cast<Corotator &>(elemCorot).getInternalForce(
      const_cast<GeomState &>(geomState),
      nodes, elemStiff, elemForce, domain->solInfo().getTimeStep(), time);
}

void
Domain::getInternalForce(GeomState &geomState, Vector& elementForce,
                         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time,
                         GeomState *refState, Vector *reactions,
                         FullSquareMatrix *mel, FullSquareMatrix *cel)
/*******************************************************************
 *
 * Purpose :
 *
 *  Compute element internal force
 *  and assemble element internal force into global internal force.
 *  Also compute follower external force contribution to
 *  residual.
 *  Also for dynamics, compute inertial and viscous force
 *  corrections (fictitious force) to residual.
 *
 * Input :
 *
 *  corotators : element corotator array
 *  geomState  : current node geometric state
 *  nodes      : undeformed nodal coordinate set
 *
 * Output :
 *
 *  residual   : residual vector = external force - internal force
 *                                 - fictitious force
 *
 *****************************************************************/

{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // mpc needs lambda for nonlinear statics
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  if(elemAdj.empty()) makeElementAdjacencyLists();
  if(time != domain->solInfo().initialTime) newDeletedElements.clear();

  if(matrixTimers) matrixTimers->formTime -= getTime();
  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    // Get updated tangent stiffness matrix and element internal force
    if(corotators[iele] && !solInfo().getNLInfo().linearelastic) {
      getElemInternalForce(geomState, pseudoTime, refState, *corotators[iele], elementForce.data(), kel[iele]);
      if(sinfo.newmarkBeta == 0) handleElementDeletion(iele, geomState, pseudoTime, *corotators[iele], elementForce.data());
    }
    // Or, get linear elastic element internal force
    else {
      Vector disp(packedEset[iele]->numDofs());
      getElementDisp(iele, geomState, disp);
      kel[iele].multiply(disp, elementForce, 1.0);
      if(solInfo().getNLInfo().linearelastic && packedEset[iele]->isFreeplayElement()) {
        Vector f(packedEset[iele]->numDofs());
        getElemInternalForce(geomState, pseudoTime, refState, *corotators[iele], f.data(), kel[iele]);
        for(int idof = 0; idof < f.size(); ++idof) elementForce[idof] += f[idof];
      }
    }

    // Add configuration-dependent external forces and their element stiffness contributions
    getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
                         (corotators[iele]), kel[iele], lambda, time, false, conwep);

    if((domain->solInfo().galerkinPodRom || domain->solInfo().DEIMBasisPod) && packedEset[iele]->hasRot() && !solInfo().getNLInfo().linearelastic) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), kel[iele], iele, false);
    }

    // Transform internal force vector to nodal frame (note: the stiffness matrix is transformed prior to assembly in Domain::makeSparseOps)
    transformVector(elementForce, iele);

    // Assemble element internal force into residual force vector
    for(int idof = 0; idof < kel[iele].dim(); ++idof) {
      int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(uDofNum >= 0)
        residual[uDofNum] -= elementForce[idof];
      else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] += elementForce[idof];
      }
    }
  }
  if(matrixTimers) matrixTimers->formTime += getTime();

  if(sinfo.isDynam() && mel && !solInfo().getNLInfo().linearelastic && !solInfo().quasistatic)
    getFictitiousForce(geomState, elementForce, kel, residual, time, refState, reactions, mel, false, corotators, cel);
}

void
Domain::getWeightedInternalForceOnly(const std::map<int, double> &weights,
                                     GeomState &geomState, Vector& elementForce,
                                     Corotator **corotators, FullSquareMatrix *kel,
                                     Vector &residual, double lambda, double time,
                                     GeomState *refState, FullSquareMatrix *mel, FullSquareMatrix *kelCopy)
{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // MPC needs lambda for nonlinear statics
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  if(elemAdj.empty()) makeElementAdjacencyLists();
 
  Vector LinearElForce(maxNumDOFs,0.0);
  Vector displacement(residual.size(),0.0);
  if(kelCopy)
    geomState.get_tot_displacement(displacement,false);
 
  for(std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iele = it->first;

    elementForce.zero();
    FullSquareMatrix &elementStiff = kel[iele];
 
    //option of remove linear component from nonlinear force, default is false 
    if(kelCopy) {
      LinearElForce.zero();
      getElemKtimesU(iele,elementStiff.dim(),displacement,LinearElForce.data(),kelCopy,(double *) dbg_alloca(sizeof(double)*maxNumDOFs*maxNumDOFs));
    }
  
    // Get updated tangent stiffness matrix and element internal force
    if(corotators[iele] && !solInfo().getNLInfo().linearelastic) {
      getElemInternalForce(geomState, pseudoTime, refState, *corotators[iele], elementForce.data(), elementStiff);
    }
    else {
      Vector disp(packedEset[iele]->numDofs());
      getElementDisp(iele, geomState, disp);
      kel[iele].copy(packedEset[iele]->stiffness(nodes, kel[iele].data())); // XXX
      kel[iele].multiply(disp, elementForce, 1.0);
      if(solInfo().getNLInfo().linearelastic && packedEset[iele]->isFreeplayElement()) {
        Vector f(packedEset[iele]->numDofs());
        getElemInternalForce(geomState, pseudoTime, refState, *corotators[iele], f.data(), kel[iele]);
        for(int idof = 0; idof < f.size(); ++idof) elementForce[idof] += f[idof];
      }
    }

    // Add configuration-dependent external forces and their element stiffness contributions
    if(domain->solInfo().reduceFollower) {
      getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
                           (corotators[iele]), elementStiff, lambda, time, false, conwep);
    }

    if(packedEset[iele]->hasRot() && !solInfo().getNLInfo().linearelastic) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iele, false);
    }
   
    // Apply lumping weight 
    const double lumpingWeight = it->second;
    elementForce  *= lumpingWeight;
    LinearElForce *= lumpingWeight;

    // Assemble element internal force into residual force vector
    const int elemDofCount = elementStiff.dim();
    for(int iDof = 0; iDof < elemDofCount; ++iDof) {
      const int dofId = c_dsa->getRCN((*allDOFs)[iele][iDof]);
      if (dofId >= 0) {
        residual[dofId] -= elementForce[iDof]; 
        if(kelCopy)
          residual[dofId] += LinearElForce[iDof];        
      }
    }
  }

  if(!domain->solInfo().reduceFollower && !domain->solInfo().DEIMPodRom) {
    getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, NULL, false);
  }

  if(sinfo.isDynam() && mel && !solInfo().getNLInfo().linearelastic && !solInfo().quasistatic) {
    getWeightedFictitiousForceOnly(weights, geomState, elementForce, kel, residual, time, refState, NULL, mel, false);
  }
}

void
Domain::getUDEIMInternalForceOnly(const std::map<int, std::vector<int> > &weights,
                                  GeomState &geomState, Vector& elementForce,
                                  Corotator **corotators, FullSquareMatrix *kel,
                                  Vector &residual, int dispSize, double lambda, double time,
                                  GeomState *refState, FullSquareMatrix *mel, FullSquareMatrix *kelCopy)
{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // MPC needs lambda for nonlinear statics
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  if(elemAdj.empty()) makeElementAdjacencyLists();

  int uDofCounter = 0;
  Vector LinearElForce(maxNumDOFs,0.0);  
  Vector displacement(dispSize,0.0);
  if(kelCopy)
    geomState.get_tot_displacement(displacement,false);

  for(std::map<int, std::vector<int> >::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iele = it->first;
    const std::vector<int> DOFvector(it->second);

    elementForce.zero();
    FullSquareMatrix &elementStiff = kel[iele];

    //option of remove linear component from nonlinear force, default is false 
    if(kelCopy) {
      LinearElForce.zero();
      getElemKtimesU(iele,elementStiff.dim(),displacement,LinearElForce.data(),kelCopy,(double *) dbg_alloca(sizeof(double)*maxNumDOFs*maxNumDOFs));
    }

    // Get updated tangent stiffness matrix and element internal force
    if(const Corotator *elementCorot = corotators[iele]) {
      getElemInternalForce(geomState, pseudoTime, refState, *elementCorot, elementForce.data(), elementStiff);
    }

    // Add configuration-dependent external forces and their element stiffness contributions
    if(domain->solInfo().reduceFollower) {
      getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
                           (corotators[iele]), elementStiff, lambda, time, false, conwep);
    }

    if(packedEset[iele]->hasRot()) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iele, false);
    }

    // Assemble element internal force into residual force vector
    const int elemDofCount = elementStiff.dim();
    for(std::vector<int>::const_iterator DOFit = DOFvector.begin(); DOFit != DOFvector.end(); DOFit++) {
      const int dofId = c_dsa->getRCN((*allDOFs)[iele][*DOFit]);
      if (dofId >= 0) {
        residual[uDofCounter] -= elementForce[*DOFit];
        if(kelCopy)
          residual[uDofCounter] += LinearElForce[*DOFit];
        uDofCounter += 1;
      }
    }
  }

//  if(!domain->solInfo().reduceFollower) {
//    getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, NULL, false);
//  }

  if(sinfo.isDynam() && mel) getUDEIMFictitiousForceOnly(weights, geomState, elementForce, kel, residual, time, refState, NULL, mel, false);
}

void
Domain::getUnassembledNonLinearInternalForce(GeomState &geomState, Vector& elementForce,
                                             Corotator **corotators, FullSquareMatrix *kel,
                                             Vector &unassemResidual, 
                                             std::map<int, std::pair<int,int> > &uDOFaDOFmap,
                                             double lambda, double time, int tIndex,
                                             GeomState *refState, Vector *reactions, FullSquareMatrix *mel,
                                             FullSquareMatrix *kelCopy)
{
  //this function creates the assembled and unassembled force snap shots needed to compute the UDEIM ROM basis
  int DOFcounter = 0; //counter variable for navigating unassembled dofs
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // mpc needs lambda for nonlinear statics
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  if(elemAdj.empty()) makeElementAdjacencyLists();

  Vector LinearElForce(maxNumDOFs,0.0);
  Vector displacement(domain->numUncon(),0.0);
  geomState.get_tot_displacement(displacement,false);

  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();
    LinearElForce.zero();

    // Get updated tangent stiffness matrix and element internal force
    getElemInternalForce(geomState, pseudoTime, refState, *corotators[iele], elementForce.data(), kel[iele]);

    getElemKtimesU(iele,kel[iele].dim(),displacement,LinearElForce.data(),kelCopy,(double *) dbg_alloca(sizeof(double)*maxNumDOFs*maxNumDOFs));

    // Add configuration-dependent external forces and their element stiffness contributions
    getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
    (corotators[iele]), kel[iele], lambda, time, false, conwep);

    if((domain->solInfo().galerkinPodRom || domain->solInfo().UDEIMBasisPod) && packedEset[iele]->hasRot() && !solInfo().getNLInfo().linearelastic) {
     // Transform element stiffness and force to solve for the increment in the total rotation vector
     transformElemStiffAndForce(geomState, elementForce.data(), kel[iele], iele, false);
    }

    // Assemble element internal force into residual force vector
    for(int idof = 0; idof < kel[iele].dim(); ++idof) {
      int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(uDofNum >= 0){
        unassemResidual[DOFcounter] -= elementForce[idof];
        unassemResidual[DOFcounter] += LinearElForce[idof];
        if(tIndex == 0){//initialize map from nunassembled, unconstrained DOFS to their elements/elemental dofs
          std::pair<int,int> uDof_key_value(std::make_pair(iele,idof));
          uDOFaDOFmap.insert(std::make_pair(DOFcounter,uDof_key_value));
        }
        DOFcounter += 1;
      }else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] += elementForce[idof];
      }
    }
  }
  if(sinfo.isDynam() && mel) getUnassembledFictitiousForce(geomState, elementForce, kel, unassemResidual, time, refState, reactions, mel, false);
}

void
Domain::getFictitiousForce(GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel, Vector &residual,
                           double time, GeomState *refState, Vector *reactions, FullSquareMatrix *mel,
                           bool compute_tangents, Corotator **corotators, FullSquareMatrix *cel)
{
  if(matrixTimers) matrixTimers->formTime -= getTime();
  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();
    getElemFictitiousForce(iele, geomState, elementForce.data(), kel[iele], time, refState, mel[iele], compute_tangents, corotators[iele], cel);
    transformVector(elementForce, iele);

    // Assemble element force into residual force vector
    for(int idof = 0; idof < kel[iele].dim(); ++idof) {
      int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(uDofNum >= 0)
        residual[uDofNum] -= elementForce[idof];
      else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] += elementForce[idof];
      }
    }
  }
  if(matrixTimers) matrixTimers->formTime += getTime();
}

void
Domain::getUnassembledFictitiousForce(GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel,
                                      Vector &unassemResidual, double time, GeomState *refState, Vector *reactions, 
                                      FullSquareMatrix *mel, bool compute_tangents)
{
  int DOFcounter = 0;
  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    getElemFictitiousForce(iele, geomState, elementForce.data(), kel[iele], time, refState, mel[iele], compute_tangents);

    // Assemble element force into residual force vector
    for(int idof = 0; idof < kel[iele].dim(); ++idof) {
      int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(uDofNum >= 0){
        unassemResidual[DOFcounter] -= elementForce[idof];
        DOFcounter += 1;
      }else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] += elementForce[idof];
      }
    }
  }
}

void
Domain::getElemFictitiousForce(int iele, GeomState &geomState, double *_fel, FullSquareMatrix &_kel,
                               double time, GeomState *refState, FullSquareMatrix &_mel,
                               bool compute_tangents, Corotator *elemCorot, FullSquareMatrix *celArray)
{
#ifdef USE_EIGEN3
  // add the correction to the residual and tangent stiffness due to the inertial effects of 
  // element iele with rotation dofs. Currently implemented for lumped mass matrix only
  // (or more specifically, mass matrices with decoupled rotational and translational diagonal blocks),
  // and elements with either 3 or 6 dofs per node.
  double &beta = domain->solInfo().newmarkBeta,
         &gamma = sinfo.newmarkGamma,
         &alphaf = sinfo.newmarkAlphaF,
         &alpham = sinfo.newmarkAlphaM,
          dt = domain->solInfo().getTimeStep();

  if(elemCorot && !elemCorot->useDefaultInertialStiffAndForce()) {

    // XXX consider nodal frames here
    int numDofs = packedEset[iele]->numDofs();
    FullSquareMatrix _kel2(numDofs);

    elemCorot->getInertialStiffAndForce(refState, geomState, nodes, _kel2, _fel, dt, time,
                                        beta, gamma, alphaf, alpham);
    _kel += _kel2;
  }
  else if(packedEset[iele]->hasRot()) {
    int numDofs = packedEset[iele]->numDofs();
    int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
    int dofsPerNode = (numDofs-packedEset[iele]->getNumMPCs())/numNodes;
    if((dofsPerNode == 3 || dofsPerNode == 6) && (numDofs-packedEset[iele]->getNumMPCs())%numNodes == 0) {
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mel(&_mel[0][0],numDofs,numDofs), kel(&_kel[0][0],numDofs,numDofs);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > fel(_fel,numDofs);
      Eigen::Matrix3d J, K;
      Eigen::Vector3d f;
      int *nodes = packedEset[iele]->nodes();
      for(int i=0; i<numNodes; ++i) {

        int offset = (dofsPerNode == 6) ? 6*i+3 : 3*i;
        Eigen::Matrix3d J = mel.block(offset,offset,3,3);
        if((J.array() == 0).all()) continue;

        Eigen::Vector3d f;
        Eigen::Matrix3d K;
        getNodeFictitiousForce(nodes[i], geomState, time, refState, J, f, K, compute_tangents);

        fel.segment(offset,3) += f;
        if(compute_tangents) kel.block(offset,offset,3,3) += K;
      }

      if(celArray && (compute_tangents || sinfo.mtypeDamp != 0) && !domain->solInfo().galerkinPodRom && !packedEset[iele]->isConstraintElement()) {
        // contribution to the tangent stiffness matrix and/or force vector due to Rayleigh damping for implicit generalized-alpha method
        // note #1: Rayleigh damping is not applied to constraint elements, by convention
        // note #2: Rayleigh damping is currently not supported for ROM
        // note #3: This function is currently not used for explicit dynamics (celArray is NULL). This means that for explicit dynamics
        //          (a) damping moments are always "follower" type, and (b) the contribution of the damping to the tangent stiffness is not included
        Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > cel(&celArray[iele][0][0],numDofs,numDofs);
        FullSquareMatrix cel_basic;
        if(!domain->solInfo().basicDofCoords) {
          cel_basic.copy(celArray[iele]);
          transformMatrixInv(cel_basic, iele);
          new (&cel) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(&cel_basic[0][0],numDofs,numDofs);
        }
        Eigen::Vector3d Psi,F;
        Eigen::Matrix3d Fx,RFx;
        Eigen::Matrix<double,6,1> inc_displacement;
        Eigen::Matrix<double,Eigen::Dynamic,1> V(numDofs);
        V.setZero();
        for(int i=0; i<numNodes; ++i) {
          Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&geomState[nodes[i]].R[0][0],3,3), R_n(&(*refState)[nodes[i]].R[0][0],3,3);
          mat_to_vec<double>(R_n.transpose()*R, Psi);
          if(dofsPerNode == 6) {
            Eigen::Map<Eigen::Matrix<double,6,1> > V_n(&(*refState)[nodes[i]].v[0]), A_n(&(*refState)[nodes[i]].a[0]);
            inc_displacement << geomState[nodes[i]].x-(*refState)[nodes[i]].x,
                                geomState[nodes[i]].y-(*refState)[nodes[i]].y,
                                geomState[nodes[i]].z-(*refState)[nodes[i]].z,
                                Psi[0], Psi[1], Psi[2];
            V.segment(6*i,6) = gamma/(dt*beta)*inc_displacement + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
          }
          else {
            Eigen::Map<Eigen::Matrix<double,3,1> > V_n(&(*refState)[nodes[i]].v[3]), A_n(&(*refState)[nodes[i]].a[3]);
            V.segment(3*i,3) = (gamma/(dt*beta)*Psi + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n);
          }
        }
        for(int i=0; i<numNodes; ++i) {
          Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&geomState[nodes[i]].R[0][0],3,3), R_n(&(*refState)[nodes[i]].R[0][0],3,3);
          Eigen::Array<double,18,1> dconst;
          dconst << R_n(0,0), R_n(0,1), R_n(0,2), R_n(1,0), R_n(1,1), R_n(1,2), R_n(2,0), R_n(2,1), R_n(2,2),
                    R(0,0), R(0,1), R(0,2), R(1,0), R(1,1), R(1,2), R(2,0), R(2,1), R(2,2);
          Simo::Jacobian<double,Simo::IncrementalRotationVector> dPsidq(dconst, Eigen::Array<int,0,1>::Zero());
          int offset = (dofsPerNode == 6) ? 6*i+3 : 3*i;
          F = cel.block(offset,0,3,numDofs)*V;
          if(sinfo.mtypeDamp == 0) { // damping moments are "axial" type
            kel.block(0,offset,numDofs,3) += (gamma/(dt*beta))*cel.block(0,offset,numDofs,3)*(dPsidq(Eigen::Vector3d::Zero(),time)-Eigen::Matrix3d::Identity());
            Fx <<     0, -F[2],  F[1],
                   F[2],     0, -F[0],
                  -F[1],  F[0],     0;
            kel.block(offset,offset,3,3) += 0.5*Fx;
          }
          else { // damping moments are "follower" type (default)
            Eigen::Vector3d RF = R*F;
            fel.segment(offset,3) += (RF-F);
            if(compute_tangents) {
              Eigen::Matrix3d D = dPsidq(Eigen::Vector3d::Zero(),time);
              for(int j=0; j<numNodes; ++j) {
                Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > Rj(&geomState[nodes[j]].R[0][0],3,3);
                if(dofsPerNode == 6) {
                  kel.block(6*j+0,offset,3,3) += (gamma/(dt*beta))*cel.block(6*j+0,offset,3,3)*(D-Eigen::Matrix3d::Identity());
                  kel.block(6*j+3,offset,3,3) += (gamma/(dt*beta))*(Rj*cel.block(6*j+3,offset,3,3)*D-cel.block(6*j+3,offset,3,3));
                  kel.block(6*j+3,offset-3,3,3) += gamma/(dt*beta)*(Rj-Eigen::Matrix3d::Identity())*cel.block(6*j+3,offset-3,3,3);
                } 
                else {
                  kel.block(3*j,offset,3,3) += (gamma/(dt*beta))*(Rj*cel.block(3*j,offset,3,3)*D-cel.block(3*j,offset,3,3));
                }
              }
              RFx <<      0, -RF[2],  RF[1],
                      RF[2],      0, -RF[0],
                     -RF[1],  RF[0],      0;
              kel.block(offset,offset,3,3) -= 0.5*RFx;
            }
          }
        }
      }
      delete [] nodes;
    }
  }

  if(iele >= elemAdj.size()) return;

  // treatment of discrete inertias adjacent to the element
  for(std::vector<std::pair<DMassData*,std::vector<int> > >::iterator it = elemAdj[iele].dimass.begin(); it != elemAdj[iele].dimass.end(); ++it) {

    DMassData *current = it->first;
    std::vector<int> &eledofs = it->second;

    int idof = current->dof;
    int jdof = (current->jdof > -1) ? current->jdof : idof;

    Eigen::Matrix3d J = Eigen::Matrix3d::Zero(); 
    J(idof-3,jdof-3) = current->diMass;
    if(idof != jdof) J(jdof-3,idof-3) = current->diMass;
    if(!domain->solInfo().basicDofCoords) {
      transformMatrixInv(J.data(), current->node, false);
    }

    Eigen::Vector3d f;
    Eigen::Matrix3d K;
    getNodeFictitiousForce(current->node, geomState, time, refState, J, f, K, compute_tangents);

    if(celArray && (compute_tangents || sinfo.mtypeDamp != 0)  && !domain->solInfo().galerkinPodRom) {
      // contribution to the tangent stiffness matrix and/or force vector due to Rayleigh damping for implicit generalized-alpha method
      Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&geomState[current->node].R[0][0],3,3), R_n(&(*refState)[current->node].R[0][0],3,3);
      Eigen::Array<double,18,1> dconst;
      dconst << R_n(0,0), R_n(0,1), R_n(0,2), R_n(1,0), R_n(1,1), R_n(1,2), R_n(2,0), R_n(2,1), R_n(2,2),
                R(0,0), R(0,1), R(0,2), R(1,0), R(1,1), R(1,2), R(2,0), R(2,1), R(2,2);
      Simo::Jacobian<double,Simo::IncrementalRotationVector> dPsidq(dconst, Eigen::Array<int,0,1>::Zero());
      Eigen::Vector3d Psi;
      mat_to_vec<double>(R_n.transpose()*R, Psi);
      Eigen::Map<Eigen::Matrix<double,3,1> > V_n(&(*refState)[current->node].v[3]), A_n(&(*refState)[current->node].a[3]);
      Eigen::Vector3d F = sinfo.alphaDamp*J*(gamma/(dt*beta)*Psi + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n);
      if(sinfo.mtypeDamp == 0) { // damping moments are "axial" type
        Eigen::Matrix3d Fx;
        Fx <<     0, -F[2],  F[1],
               F[2],     0, -F[0],
              -F[1],  F[0],     0;
        K += 0.5*Fx + (gamma/(dt*beta))*sinfo.alphaDamp*J*(dPsidq(Eigen::Vector3d::Zero(),time)-Eigen::Matrix3d::Identity());
      }
      else { // damping moments are "follower" type (default)
        Eigen::Vector3d RF = R*F;
        f += (RF-F);
        if(compute_tangents) {
          Eigen::Matrix3d RFx;
          RFx <<      0, -RF[2],  RF[1],
                  RF[2],      0, -RF[0],
                 -RF[1],  RF[0],      0;
          K += -0.5*RFx + (gamma/(dt*beta))*sinfo.alphaDamp*(R*J*(dPsidq(Eigen::Vector3d::Zero(),time))-J);
        }
      }
    }

    for(int j = 0; j < 3; ++j) {
      _fel[eledofs[j]] += f[j];
      if(compute_tangents) {
        for(int k = 0; k < 3; ++k)
          _kel[eledofs[j]][eledofs[k]] += K(j,k);
      }
    }
  }
#endif
}

#ifdef USE_EIGEN3
void
Domain::getNodeFictitiousForce(int inode, GeomState &geomState, double time, GeomState *refState, Eigen::Matrix3d &J,
                               Eigen::Vector3d &f, Eigen::Matrix3d &K, bool compute_tangents)
{
  double &beta = domain->solInfo().newmarkBeta,
         &gamma = sinfo.newmarkGamma,
         &alphaf = sinfo.newmarkAlphaF,
         &alpham = sinfo.newmarkAlphaM,
          dt = domain->solInfo().getTimeStep();

  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&geomState[inode].R[0][0]), R_n(&(*refState)[inode].R[0][0]);
  Eigen::Matrix3d T, Tdot;
  Eigen::Map<Eigen::Vector3d> Psi(&geomState[inode].theta[0]), Psi_n(&(*refState)[inode].theta[0]),
                              V_n(&(*refState)[inode].v[3]), A_n(&(*refState)[inode].a[3]);
  Eigen::Vector3d V, A;
  V << geomState[inode].v[3], geomState[inode].v[4], geomState[inode].v[5];
  A << geomState[inode].a[3], geomState[inode].a[4], geomState[inode].a[5];

  if(beta == 0) { // compute the tangent stiffness and/or force correction due to rotary inertia for explicit central difference
    // V is either the convected angular velocity at t^{n+1/2} for FOM or ROM model II or model III,
    //   or the convected angular velocity at current snapshot after projection for explicit ROM "training"
    if(domain->solInfo().galerkinPodRom || domain->solInfo().samplingPodRom || domain->solInfo().DEIMBasisPod || domain->solInfo().UDEIMBasisPod) {
      tangential_transf<double>(Psi, T);
      Eigen::Vector3d Psidot = T.inverse()*V;
      tangential_transf_dot<double>(Psi, Psidot, Tdot);
      // f = J'*(T^t*J'*T)^{-1}*T*R*(J*Tdot*Psidot + V.cross(J*V))
      // note: T*R = T^t and J' is the assembled nodal inertia, i.e. J' = Jn[inode]
      // --> f = J'*(J'*T)^{-1}*(J*Tdot*Psidot + V.cross(J*V))
      f = Jn[inode]*(Jn[inode]*T).inverse()*(J*Tdot*Psidot + V.cross(J*V));
    }
    else {
      f = R*V.cross(J*V);
    }

    // TODO: compute tangents (for critical timestep estimate)
    if(compute_tangents) K.setZero();
  }
  else { // compute the tangent stiffness and/or force correction due to rotary inertia for implicit generalized-alpha
    if(domain->solInfo().samplingPodRom || domain->solInfo().ROMPostProcess) {
      // V and A are the convected angular velocity and acceleration at current snapshot after projection
      tangential_transf<double>(Psi, T);
      Eigen::Matrix3d Tinv = T.inverse();
      Eigen::Vector3d Psidot = Tinv*V;
      tangential_transf_dot<double>(Psi, Psidot, Tdot);
      f = T.transpose()*(J*A + V.cross(J*V)) - J*Tinv*(A-Tdot*Psidot);
    }
    else {
      // for HFM, V and A are the first and second time-derivatives of the total rotation vector at t_n
      // and for ROM, V and A are the first and second time-derivatives of the total rotation vector at t_n
      // at t_0, A should be zero!
      if(time != 0) {
        Eigen::Vector3d incd;
        if(domain->solInfo().galerkinPodRom) {
          incd = Psi - Psi_n;
        }
        else {
          Eigen::Matrix3d dR = R_n.transpose()*R;
          mat_to_vec<double>(dR, incd);
        }
        // compute the convected angular velocity at t^{n+1-alphaf} for HFM, or first time-derivative of total rotation vector at t^{n+1-alphaf} for ROM
        V = gamma/(dt*beta)*incd + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
        // compute the convected angular acceleration at t^{n+1-alpham} for HFM, or second first time-derivative of total rotation vector at t^{n+1-alphaf} for ROM
        A = (1-alpham)/(dt*dt*beta*(1-alphaf))*incd - (1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;
      }

      // compute the inertial force
      if(domain->solInfo().galerkinPodRom) {
        tangential_transf<double>(Psi, T);
        tangential_transf_dot<double>(Psi, V, Tdot);
        f = T.transpose()*( J*(T*A + Tdot*V) + (T*V).cross(J*T*V) ); // note T*R = T.transpose()
      }
      else {
        f = R*(J*A + V.cross(J*V));
      }
      // subtract the linear part which is added in probDesc->formRHScorrector
      if(time != 0) f -= J*A;

      if(compute_tangents) { // tangent stiffness contribution of the inertial force

        Eigen::Matrix<double,3,1> q;
        if(domain->solInfo().galerkinPodRom) {
          Eigen::Array<double,23,1> dconst;
          Eigen::Array<int,0,1> iconst;
          dconst << J(0,0), J(0,1), J(0,2), J(1,0), J(1,1), J(1,2), J(2,0), J(2,1), J(2,2),
                    A_n[0], A_n[1], A_n[2],
                    V_n[0], V_n[1], V_n[2],
                    Psi_n[0], Psi_n[1], Psi_n[2],
                    beta, gamma, alphaf, alpham, dt;

          // evaluate the jacobian of the inertial force
          Simo::Jacobian<double,Simo::InertialType2ForceFunction> dFdq(dconst,iconst);
          q << geomState[inode].theta[0], geomState[inode].theta[1], geomState[inode].theta[2];
          K = dFdq(q, time);
        }
        else {
          Eigen::Array<double,38,1> dconst;
          Eigen::Array<int,0,1> iconst;
          dconst << J(0,0), J(0,1), J(0,2), J(1,0), J(1,1), J(1,2), J(2,0), J(2,1), J(2,2),
                    A_n[0], A_n[1], A_n[2],
                    V_n[0], V_n[1], V_n[2],
                    R_n(0,0), R_n(0,1), R_n(0,2), R_n(1,0), R_n(1,1), R_n(1,2), R_n(2,0), R_n(2,1), R_n(2,2),
                    R(0,0), R(0,1), R(0,2), R(1,0), R(1,1), R(1,2), R(2,0), R(2,1), R(2,2),
                    beta, gamma, alphaf, alpham, dt;

          // evaluate the jacobian of the inertial force
          Simo::Jacobian<double,Simo::InertialType1ForceFunction> dFdq(dconst,iconst);
          q = Eigen::Vector3d::Zero();
          K = dFdq(q, time);
        }

        // subtract the linear part which is added to the dynamic tangent stiffness in probDesc->reBuild
        K -= (1-alpham)/((1-alphaf)*(dt*dt*beta))*J;
      }
    }
  }
}
#endif

void
Domain::getWeightedFictitiousForceOnly(const std::map<int, double> &weights, GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel,
                                       Vector &residual, double time, GeomState *refState, Vector *reactions,
                                       FullSquareMatrix *mel, bool compute_tangents)
{
  for (std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iele = it->first;
    const double lumpingWeight = it->second;

    elementForce.zero();
    FullSquareMatrix kel2(kel[iele].dim());
    kel2.zero();

    getElemFictitiousForce(iele, geomState, elementForce.data(), kel2, time, refState, mel[iele], compute_tangents);

    elementForce *= lumpingWeight;

    if(compute_tangents) {
      kel2 *= lumpingWeight;
      kel[iele] += kel2;
    }

    // Assemble element force into residual force vector
    for(int idof = 0; idof < kel[iele].dim(); ++idof) {
      int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(uDofNum >= 0)
        residual[uDofNum] -= elementForce[idof];
      else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] += elementForce[idof];
      }
    }
  }
}

void
Domain::getUDEIMFictitiousForceOnly(const std::map<int, std::vector<int> > &weights, GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel,
                                    Vector &residual, double time, GeomState *refState, Vector *reactions,
                                    FullSquareMatrix *mel, bool compute_tangents)
{

  int uDofCounter = 0;  
  
  for (std::map<int, std::vector<int> >::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iele = it->first;
    const std::vector<int> DOFvector = it->second;

    elementForce.zero();

    getElemFictitiousForce(iele, geomState, elementForce.data(), kel[iele], time, refState, mel[iele], compute_tangents);

    // Assemble element force into residual force vector
   for(std::vector<int>::const_iterator DOFit = DOFvector.begin(); DOFit != DOFvector.end(); DOFit++) {
      int uDofNum = c_dsa->getRCN((*allDOFs)[iele][*DOFit]);
      if(uDofNum >= 0){
        residual[uDofCounter] -= elementForce[*DOFit];
        uDofCounter += 1;
      }
      else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][*DOFit]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] += elementForce[*DOFit];
      }
    }
  }
}

void
Domain::assembleNodalInertiaTensors(FullSquareMatrix *melArray)
{
#ifdef USE_EIGEN3
  if(elemAdj.empty()) makeElementAdjacencyLists();
  Jn.resize(numnodes);
  for(int i = 0; i < numnodes; ++i) Jn[i].setZero();

  for(int iele = 0; iele < numele; ++iele) {

    if(packedEset[iele]->hasRot()) {
      int numDofs = packedEset[iele]->numDofs();
      int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
      int dofsPerNode = (numDofs-packedEset[iele]->getNumMPCs())/numNodes;
      if((dofsPerNode == 3 || dofsPerNode == 6) && (numDofs-packedEset[iele]->getNumMPCs())%numNodes == 0) {
        Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mel(&melArray[iele][0][0],numDofs,numDofs);
        int *nodes = packedEset[iele]->nodes();
        for(int i=0; i<numNodes; ++i) {
          int offset = (dofsPerNode == 6) ? 6*i+3 : 3*i;
          Jn[nodes[i]] += mel.block(offset,offset,3,3);
        }
        delete [] nodes;
      }
    }

    if(iele >= elemAdj.size()) return;

    // add contribution of of discrete inertias adjacent to the element
    for(std::vector<std::pair<DMassData*,std::vector<int> > >::iterator it = elemAdj[iele].dimass.begin(); it != elemAdj[iele].dimass.end(); ++it) {

      DMassData *current = it->first;
      std::vector<int> &eledofs = it->second;

      int idof = current->dof;
      int jdof = (current->jdof > -1) ? current->jdof : idof;

      Jn[current->node](idof-3,jdof-3) += current->diMass;
      if(idof != jdof) Jn[current->node](jdof-3,idof-3) += current->diMass;
    }
  }
#endif
}

double
Domain::getKineticEnergy(double* velocity, FullSquareMatrix *mel)
{
  // Compute Kinetic Energy as 0.5*(vtmp^t M vtmp)
  // note #1: assuming angular velocity is convected
  // note #2: contribution due to prescribed velocity is currently not included
  Vector vtmp(velocity, numUncon());
  Vector tmpVec(numUncon(), 0.0);
  int iele, idof, jdof, dofn1, dofn2;

  // contribution of elements
  for(iele = 0; iele < numele; ++iele) {
    for(idof = 0; idof < mel[iele].dim(); ++idof) {
      dofn1 = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(dofn1 >= 0) {
        for(jdof = 0; jdof < mel[iele].dim(); ++jdof) {
          dofn2 = c_dsa->getRCN((*allDOFs)[iele][jdof]);
          if(dofn2 >= 0)
            tmpVec[dofn1] += vtmp[dofn2]*mel[iele][idof][jdof];
        }
      }
    }
  }
  double T = 0.5*(vtmp*tmpVec);

  // contribution due to discrete masses
  if(firstDiMass != NULL) {
    DMassData *current = firstDiMass;
    while(current != 0) {
      int cdof = c_dsa->locate(current->node, (1 << current->dof));
      if(current->jdof < 0 || current->dof == current->jdof) { // diagonal entry
        if(cdof > -1) {
          T += 0.5*vtmp[cdof]*current->diMass*vtmp[cdof];
        }
      }
      else { // off-diagonal entry
        int jcdof = c_dsa->locate(current->node, (1 << current->jdof));
        if(cdof > -1 && jcdof > -1) {
          T += vtmp[cdof]*current->diMass*vtmp[jcdof]; // note: 0.5 is deliberately omitted here because only the lower
                                                       // triangular part of the inertia tensor should be specified
        }
      }
      current = current->next;
    }
  }

  return T;
}

void
Domain::handleElementDeletion(int iele, GeomState &geomState, double time,
                              Corotator &elemCorot, double *elemForce)
{
  if(domain->solInfo().elementDeletion) {
    bool deleteElem = false;
    if(!domain->solInfo().deleteElements.empty()) {
      std::map<int,double>::iterator it;
#if defined(_OPENMP)
      #pragma omp critical
#endif
      if((it = domain->solInfo().deleteElements.find(packedEset[iele]->getGlNum())) != domain->solInfo().deleteElements.end() && time >= it->second) {
        domain->solInfo().deleteElements.erase(it);
        deleteElem = true;
      }
    }
    if(deleteElem || (packedEset[iele]->getProperty() && elemCorot.checkElementDeletion(geomState))) {
      std::cerr << "\rDeleting element " << std::setw(22) << std::left << packedEset[iele]->getGlNum()+1 << std::endl;
#if defined(_OPENMP)
      #pragma omp critical
#endif
      { newDeletedElements.insert(packedEset[iele]->getGlNum());
        outDeletedElements.push_back(std::pair<double,int>(time,packedEset[iele]->getGlNum())); }
      packedEset[iele]->setProp((StructProp*)NULL);
      packedEset[iele]->setPressure((PressureBCond*)NULL);
      elemAdj[iele].surfp.clear();
      if(elemForce) for(int i=0; i<packedEset[iele]->numDofs(); ++i) elemForce[i] = 0;
    }
  }
}
