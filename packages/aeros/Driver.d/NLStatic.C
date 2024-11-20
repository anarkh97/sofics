#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <set>
#include <Utils.d/dbg_alloca.h>

// New include files for Restart file
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <Corotational.d/Corotator.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/Domain.h>
#include <Driver.d/SubDomain.h>
#include <Utils.d/dofset.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/TemperatureState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/CuCSparse.h>
#include <Control.d/ControlInterface.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Driver.d/GeoSource.h>
#include <Corotational.d/MatNLCorotator.h>
#include <Element.d/Function.d/utilities.hpp>
#include <Utils.d/DistHelper.h>

void
Domain::getElemStiffAndForce(const GeomState &geomState, double time,
                             const GeomState *refState, const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
    // note: delt is the time increment between the refState and current time for dynamics and zero for statics
    double delt = (sinfo.newmarkBeta == 0) ? sinfo.getTimeStep() : (1-sinfo.newmarkAlphaF)*sinfo.getTimeStep();
    const_cast<Corotator &>(elemCorot).getStiffAndForce(
        const_cast<GeomState *>(refState),
        const_cast<GeomState &>(geomState),
        nodes, elemStiff, elemForce, delt, time);
}

void
Domain::getElemStiffAndForce(const GeomState &geomState, double time,
                             const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
  const_cast<Corotator &>(elemCorot).getStiffAndForce(
      const_cast<GeomState &>(geomState),
      nodes, elemStiff, elemForce, domain->solInfo().getTimeStep(), time);
}

void
Domain::getStiffAndForce(GeomState &geomState, Vector& elementForce,
                         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time,
                         GeomState *refState, Vector *reactions,
                         FullSquareMatrix *mel, FullSquareMatrix *cel)
/*******************************************************************
 *
 * Purpose :
 *
 *  Compute element tangential stiffness and element internal force
 *  and assemble element internal force into global internal force.
 *  Also compute configuration-dependent external force contribution
 *  to both tangential stiffness matrix and residual.
 *  Also for dynamics, compute inertial and viscous
 *  corrections to both tangential stiffness matrix and residual.
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
 *  kel        : array of element tangential stiffness matrices
 *               in current configuration
 *
 *****************************************************************/

{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // mpc needs lambda for nonlinear statics
  const bool initialTime = (sinfo.isDynam() && time == 0 && solInfo().newmarkBeta != 0); // when initialTime is true, will
                                                                                         // be solving for initial acceleration
                                                                                         // (implicit dynamics only)
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  bool compute_tangents = !initialTime && !solInfo().getNLInfo().linearelastic;
  if(elemAdj.empty()) makeElementAdjacencyLists();
  if(time != domain->solInfo().initialTime) newDeletedElements.clear();

  if(matrixTimers) matrixTimers->formTime -= getTime();
  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    // Get updated tangent stiffness matrix and element internal force
    if(corotators[iele] && !solInfo().getNLInfo().linearelastic) {
      getElemStiffAndForce(geomState, pseudoTime, refState, *corotators[iele], elementForce.data(), kel[iele]);
      if(sinfo.newmarkBeta == 0) {
        corotators[iele]->updateStates(refState, geomState, nodes, sinfo.getTimeStep());
        handleElementDeletion(iele, geomState, pseudoTime, *corotators[iele], elementForce.data());
      }
      if(initialTime && packedEset[iele]->isConstraintElement() && packedEset[iele]->hasRot()) {
        // transform constraint jacobian and hessian to solve for the initial convected acceleration
        transformElemStiff(geomState, kel[iele], iele);
      }
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
                         (corotators[iele]), kel[iele], lambda, time, compute_tangents, conwep);

    if(domain->solInfo().galerkinPodRom && packedEset[iele]->hasRot() && !solInfo().getNLInfo().linearelastic) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), kel[iele], iele, true);
    }

    // Transform internal force vector to nodal frame
    transformVector(elementForce, iele); 

    // Assemble element force into residual vector
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

  // XXX consider adding the element fictitious forces inside the loop
  if(sinfo.isDynam() && mel && !solInfo().getNLInfo().linearelastic && !solInfo().quasistatic)
    getFictitiousForce(geomState, elementForce, kel, residual, time, refState, reactions, mel, compute_tangents, corotators, cel);

  if(!solInfo().getNLInfo().unsymmetric && solInfo().newmarkBeta != 0) {
    if(matrixTimers) matrixTimers->formTime -= getTime(); 
    for(int iele = 0; iele < numele; ++iele)
      kel[iele].symmetrize();
    if(matrixTimers) matrixTimers->formTime += getTime();
  }
}

void
Domain::getFollowerForce(GeomState &geomState, Vector& elementForce,
                         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time,
                         Vector *reactions, bool compute_tangents)
{
  // Note #1: If compute_tangents is true, then the jacobian of the element follower forces 
  //          will be added to kel. If compute_tangents is false, then kel can be NULL.
  // Note #2: This function requires followedElemList to be pre-computed; ensure makeFollowedElemList
  //          is set to true in Domain::makeElementAdjacencyLists().

  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;

  for(std::vector<int>::iterator it = followedElemList.begin(), it_end = followedElemList.end(); it != it_end; ++it) {
    const int iele = *it;

    elementForce.zero();
    const int elemDofCount = packedEset[iele]->numDofs();
    FullSquareMatrix elementStiff(elemDofCount);
    elementStiff.zero();

    getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
                         (corotators[iele]), elementStiff, lambda, time, compute_tangents, conwep);

    if(domain->solInfo().galerkinPodRom && packedEset[iele]->hasRot())
      transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iele, compute_tangents);

    if(compute_tangents) kel[iele] += elementStiff;

    // Transform internal force vector to nodal frame
    transformVector(elementForce, iele);

    // Assemble element force into residual vector
    for(int iDof = 0; iDof < elemDofCount; ++iDof) {
      const int dofId = c_dsa->getRCN((*allDOFs)[iele][iDof]);
      if (dofId >= 0) {
        residual[dofId] -= elementForce[iDof];
      }
    }
  }
}

void
Domain::getElemFollowerForce(int iele, GeomState &geomState, double *_f, int bufSize,
                             Corotator *corotator, FullSquareMatrix &kel,
                             double lambda, double time, bool compute_tangents,
                             BlastLoading::BlastData *conwep)
{
  if(solInfo().getNLInfo().linearelastic == 2) return;
  Vector elementForceBuf(_f,bufSize,false);

  // 1. Treatment of element pressure
  PressureBCond *pbc;
  if((pbc = packedEset[iele]->getPressure()) != NULL) {
    Vector elementForce(bufSize);
    elementForce.zero();

    // Compute the amplified pressure due to MFTT and/or dlambda, if applicable
    MFTTData *mftt = domain->getMFTT(pbc->loadsetid);
    double loadFactor = (mftt && sinfo.isDynam()) ? mftt->getVal(std::max(time,0.0)) : domain->getLoadFactor(pbc->loadsetid);
    double p0 = pbc->val;
    pbc->val *= (lambda*loadFactor);

    // Compute element pressure force in the local coordinates using the specified blast loading function and/or amplified pressure value
    pbc->conwep = conwep;
    packedEset[iele]->computePressureForce(nodes, elementForce, &geomState, 1, time);
    pbc->val = p0;

    if(corotator) {
      // Include the "load stiffness matrix" in kel
      if(compute_tangents) {

        FullSquareMatrix k(kel.dim());
        k.zero();
        corotator->getDExternalForceDu(geomState, nodes, k, elementForce.data());
        for(int i = 0; i < kel.dim(); ++i)
          for(int j = 0; j < kel.dim(); ++j)
            kel[i][j] += k[i][j];
      }

      // Determine the elemental force for the corotated system
      corotator->getExternalForce(geomState, nodes, elementForce.data());
    }

    elementForceBuf -= elementForce;
  }

  // 2. Treatment of element thermal forces
  //    By convention phantom elements do not have thermal load
  if(domain->thermalFlag() && packedEset[iele]->getProperty() != 0) {

    Vector elementTemp(packedEset[iele]->numNodes());
    Vector elementForce(packedEset[iele]->numDofs());

    // Extract the element nodal temperatures from temprcvd and/or element property ambient temperature
    int *elemnodes = packedEset[iele]->nodes();
    for(int inod = 0; inod < packedEset[iele]->numNodes(); ++inod) {
      double t = temprcvd[elemnodes[inod]];
      elementTemp[inod] = (t == defaultTemp) ? packedEset[iele]->getProperty()->Ta : t;
    }
    delete [] elemnodes;

    // Compute element thermal force in the local coordinates
    elementForce.zero();
    packedEset[iele]->getThermalForce(nodes, elementTemp, elementForce, 1, &geomState);

    if(corotator) {
      // Include the "load stiffness matrix" in kel
      if(compute_tangents) {
        FullSquareMatrix k(kel.dim());
        k.zero();
        corotator->getDExternalForceDu(geomState, nodes, k, elementForce.data());

        for(int i = 0; i < kel.dim(); ++i)
          for(int j = 0; j < kel.dim(); ++j)
            kel[i][j] += k[i][j];
      }

      // Determine the elemental force for the corrotated system
      corotator->getExternalForce(geomState, nodes, elementForce.data());
    }

    elementForceBuf -= elementForce;
  }

  if(iele >= elemAdj.size()) return;

  // 3. Treatment of surface pressure
  for(std::vector<std::pair<int,std::vector<int> > >::iterator it = elemAdj[iele].surfp.begin(); it != elemAdj[iele].surfp.end(); ++it) {

    int i = it->first;
    std::vector<int> &eledofs = it->second;

    Vector f(neum[i]->numDofs());
    f.zero();

    // Compute the amplified pressure due to MFTT and/or dlambda, if applicable
    PressureBCond *pbc = neum[i]->getPressure();
    MFTTData *mftt = domain->getMFTT(pbc->loadsetid);
    double loadFactor = (mftt && sinfo.isDynam()) ? mftt->getVal(std::max(time,0.0)) : domain->getLoadFactor(pbc->loadsetid);
    double p0 = pbc->val;
    pbc->val *= (lambda*loadFactor);

    // Compute surface element pressure force using the specified blast loading function and/or amplified pressure value
    pbc->conwep = conwep;
    neum[i]->neumVector(nodes, f, 0, &geomState, time);

    for(int j = 0; j < neum[i]->numDofs(); ++j) {
      elementForceBuf[eledofs[j]] -= f[j];
    }

    // Include the "load stiffness matrix" in kel
    if(compute_tangents) {
      FullSquareMatrix kTan(neum[i]->numDofs());
      kTan.zero();

      neum[i]->neumVectorJacobian(nodes, kTan, 0, &geomState, time);

      for(int j = 0; j < neum[i]->numDofs(); ++j)
        for(int k = 0; k < neum[i]->numDofs(); ++k)
          kel[eledofs[j]][eledofs[k]] -= kTan[j][k];
    }

    pbc->val = p0;
  }

  // 4. Treatment of configuration dependent nodal moments
  for(std::vector<std::pair<int,std::vector<int> > >::iterator it = elemAdj[iele].cdnm.begin(); it != elemAdj[iele].cdnm.end(); ++it) {

    int i = it->first;
    std::vector<int> &eledofs = it->second;

    // if the element is an LMPC and the node has a DOF_FRM then the element force and or stiffness must be defined in the DOF_FRM
    NFrameData *cd = 0;
    if(packedEset[iele]->isMpcElement() && nodes.dofFrame(nbc[i].nnum)) {
      LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]);
      if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
         lmpcons->getSource() == mpc::TiedSurfaces)) cd = nodes.dofFrame(nbc[i].nnum);
    }

    double m0[3] = { 0, 0, 0 }, m[3], r[3], rotvar[3][3];
    MFTTData *mftt = domain->getMFTT(nbc[i].loadsetid);
    double loadFactor = (mftt && sinfo.isDynam()) ? mftt->getVal(std::max(time,0.0)) : domain->getLoadFactor(nbc[i].loadsetid);
    m0[nbc[i].dofnum-3] = lambda*loadFactor*nbc[i].val;
    transformVectorInv(m0, nbc[i].nnum, false);

    switch(nbc[i].mtype) {
      case BCond::Axial : // axial (constant) moment: m = m0
        for(int j = 0; j < 3; ++j) m[j] = m0[j];
        break;
      case BCond::Rotational : { // rotational moment: m = T^{-1}*m0
        mat_to_vec(geomState[nbc[i].nnum].R, r);
        pseudorot_var(r, rotvar);
        mat_mult_vec(rotvar, m0, m, 1);
      } break;
      case BCond::Follower : // follower moment: m = R*m0
        mat_mult_vec(geomState[nbc[i].nnum].R, m0, m, 0);
        break;
      default :
        std::cerr << " *** WARNING: selected moment type is not supported\n";
    }

    if(cd) {
      double m_copy[3] = { m[0], m[1], m[2] };
      cd->transformVector3(m_copy);
      for(int j = 0; j < 3; ++j) elementForceBuf[eledofs[j]] -= m_copy[j];
    }
    else {
      for(int j = 0; j < 3; ++j) elementForceBuf[eledofs[j]] -= m[j];
    }

    // tangent stiffness contribution: 
    if(compute_tangents) {
      switch(nbc[i].mtype) {
        case BCond::Axial : { // axial (constant) moment
          double skewm0[3][3] = { {     0, -m0[2],  m0[1] },
                                  {  m0[2],     0, -m0[0] },
                                  { -m0[1],  m0[0],    0  } };
          for(int j = 0; j < 3; ++j)
            for(int k = 0; k < 3; ++k) 
              rotvar[j][k] = 0.5*skewm0[j][k];
        } break;
        case BCond::Rotational : { // rotational moment
          double scndvar[3][3];
          pseudorot_2var(r, m0, scndvar);
          for(int j = 0; j < 3; ++j)
            for(int k = 0; k < 3; ++k)
              rotvar[j][k] = 0.5*(scndvar[j][k] + scndvar[k][j]);
        } break;
        case BCond::Follower : { // follower moment
          double skewm[3][3] = { {     0, -m[2],  m[1] },
                                 {  m[2],     0, -m[0] },
                                 { -m[1],  m[0],    0  } };
          for(int j = 0; j < 3; ++j)
            for(int k = 0; k < 3; ++k)
              rotvar[j][k] = -0.5*skewm[j][k];
        } break;
      }

      if(cd) cd->transformMatrix3(&rotvar[0][0]);

      for(int j = 0; j < 3; ++j)
        for(int k = 0; k < 3; ++k)
          kel[eledofs[j]][eledofs[k]] -= rotvar[j][k];
    }
  }

  // 5. Treatment of configuration dependent nodal forces
  for(std::vector<std::pair<int,std::vector<int> > >::iterator it = elemAdj[iele].cdnf.begin(); it != elemAdj[iele].cdnf.end(); ++it) {

    int i = it->first;
    std::vector<int> &eledofs = it->second;

    // if the element is an LMPC and the node has a DOF_FRM then the element force and or stiffness must be defined in the DOF_FRM
    NFrameData *cd = 0;
    if(packedEset[iele]->isMpcElement() && nodes.dofFrame(nbc[i].nnum)) {
      LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]);
      if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
         lmpcons->getSource() == mpc::TiedSurfaces)) cd = nodes.dofFrame(nbc[i].nnum);
    }

    double f0[3] = { 0, 0, 0 }, f[3];
    MFTTData *mftt = domain->getMFTT(nbc[i].loadsetid);
    double loadFactor = (mftt && sinfo.isDynam()) ? mftt->getVal(std::max(time,0.0)) : domain->getLoadFactor(nbc[i].loadsetid);
    f0[nbc[i].dofnum] = lambda*loadFactor*nbc[i].val;
    transformVectorInv(f0, nbc[i].nnum, false);

    mat_mult_vec(geomState[nbc[i].nnum].R, f0, f, 0); // f = R*f0

    if(cd) cd->transformVector3(f);
    for(int j = 0; j < 3; ++j) elementForceBuf[eledofs[j]] -= f[j];

    // tangent stiffness contribution: 
    if(compute_tangents) {
      double skewf[3][3] = { {     0, -f[2],  f[1] },
                             {  f[2],     0, -f[0] },
                             { -f[1],  f[0],    0  } };

      for(int j = 0; j < 3; ++j)
        for(int k = 0; k < 3; ++k)
          kel[eledofs[j]][eledofs[3+k]] += skewf[j][k];
    }
  }

  // 6. Treatment of constrained rotations
  if(elemAdj[iele].crot.size() > 0) {
    transformElemStiffAndForce_S2E(geomState, _f, kel, iele, compute_tangents);
  }
}

void
Domain::makeElementAdjacencyLists()
{
  elemAdj.resize(numele); 
  bool makeFollowedElemList = ( ((!domain->solInfo().reduceFollower || domain->solInfo().getNLInfo().linearelastic)
                                 && domain->solInfo().galerkinPodRom) || geoSource->energiesOutput() || geoSource->romExtForceOutput() );

  // 1. pressure applied to elements
  if(domain->pressureFlag()) {
    if(makeFollowedElemList) {
      for(int iele = 0; iele < numele; ++iele) {
        if(packedEset[iele]->getPressure() != NULL) followedElemList.push_back(iele);
      }
    }
  }

  // 2. pressure using surfacetopo
  //    build the map from kel array index to neum array index
  SubDomain *subCast = (numNeum > 0) ? dynamic_cast<SubDomain*>(this) : NULL;
  if(!subCast && !sommerChecked && numNeum > 0) { checkSommerTypeBC(this); }
  for(int i = 0; i < numNeum; ++i) {
    if((neum[i]->getPressure()) == NULL) continue;
    int iele = (subCast) ? subCast->globalToLocalElem(neum[i]->getAdjElementIndex()) : neum[i]->getAdjElementIndex();

    int *dofs = new int[neum[i]->numDofs()];
    neum[i]->dofs(*dsa, dofs);

    std::vector<int> eledofs(neum[i]->numDofs());
    for(int j = 0; j < neum[i]->numDofs(); ++j) {
      for(int k = 0; k < allDOFs->num(iele); ++k)
        if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
    }

    elemAdj[iele].surfp.push_back(std::pair<int,std::vector<int> >(i,eledofs));
    if(makeFollowedElemList) followedElemList.push_back(iele);

    delete [] dofs;
  }

  // 3. configuration-dependent nodal moments
  //    build map from kel array index to nbc array index
  for(int i = 0; i < numNeuman; ++i) {
    if((nbc[i].type == BCond::Forces || nbc[i].type == BCond::Usdf || nbc[i].type == BCond::Actuators) 
       && (nbc[i].dofnum == 3 || nbc[i].dofnum == 4 || nbc[i].dofnum == 5)) {
      int dofs[3];
      dsa->number(nbc[i].nnum, DofSet::XYZrot, dofs);
      bool foundAdj = false;
      for(int inode = 0; inode < nodeToElem->num(nbc[i].nnum); ++inode) { // loop over the elements attached to the node
                                                                          // at which the nodal moment is applied
        int iele = (*nodeToElem)[nbc[i].nnum][inode];
        std::vector<int> eledofs(3);
        for(int j = 0; j < 3; ++j) {
          eledofs[j] = -1;
          for(int k = 0; k < allDOFs->num(iele); ++k)
            if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
        }
        if(eledofs[0] != -1 && eledofs[1] != -1 && eledofs[2] != -1) {
          // found an element with the 3 rotation dofs of node nbc[i].nnum
          elemAdj[iele].cdnm.push_back(std::pair<int,std::vector<int> >(i,eledofs));
          if(makeFollowedElemList) followedElemList.push_back(iele);
          foundAdj = true;
          break;
        }
      }
      if(!foundAdj) std::cerr << " *** WARNING: could not find adjacent element for configuration-dependent nodal moment\n";
    }
  }

  // 4. configuration-dependent nodal forces
  //    build map from kel array index to nbc array index
  for(int i = 0; i < numNeuman; ++i) {
    if((nbc[i].type == BCond::Forces || nbc[i].type == BCond::Usdf || nbc[i].type == BCond::Actuators) 
       && (nbc[i].dofnum == 0 || nbc[i].dofnum == 1 || nbc[i].dofnum == 2)
       && nbc[i].mtype == BCond::Follower) {
      int dofs[6];
      dsa->number(nbc[i].nnum, DofSet::XYZdisp | DofSet::XYZrot, dofs);
      bool foundAdj = false;
      for(int inode = 0; inode < nodeToElem->num(nbc[i].nnum); ++inode) { // loop over the elements attached to the node
                                                                          // at which the nodal force is applied
        int iele = (*nodeToElem)[nbc[i].nnum][inode];
        std::vector<int> eledofs(6);
        for(int j = 0; j < 6; ++j) {
          eledofs[j] = -1;
          for(int k = 0; k < allDOFs->num(iele); ++k)
            if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
        }
        if(eledofs[0] != -1 && eledofs[1] != -1 && eledofs[2] != -1 &&
           eledofs[3] != -1 && eledofs[4] != -1 && eledofs[5] != -1) {
          // found an element with the 6 translation & rotation dofs of node nbc[i].nnum 
          elemAdj[iele].cdnf.push_back(std::pair<int,std::vector<int> >(i,eledofs));
          if(makeFollowedElemList) followedElemList.push_back(iele);
          foundAdj = true;
          break;
        }
      }
      if(!foundAdj) std::cerr << " *** WARNING: could not find adjacent element for configuration-dependent nodal force\n";
    }
  }

  // 5. thermal forces
  if(domain->thermalFlag()) {
    if(!temprcvd) initNodalTemperatures();
    if(makeFollowedElemList) {
      for(int iele = 0; iele < numele; ++iele) followedElemList.push_back(iele);
    }
  }

  // 6. discrete inertias
  if(firstDiMass != NULL) {
    DMassData *current = firstDiMass;
    while(current != 0) {
      int idof = current->dof;
      int jdof = (current->jdof > -1) ? current->jdof : idof;
      if((idof == 3 || idof == 4 || idof == 5) && (jdof == 3 || jdof == 4 || jdof == 5)) {

        int dofs[3];
        dsa->number(current->node, DofSet::XYZrot, dofs);

        bool foundAdj = false;
        for(int inode = 0; inode < nodeToElem->num(current->node); ++inode) { // loop over the elements attached to the node
                                                                              // at which the discrete mass is located
          int iele = (*nodeToElem)[current->node][inode];

          if(packedEset[iele]->isMpcElement() && nodes.dofFrame(current->node)) {
            LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]);
            if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
               lmpcons->getSource() == mpc::TiedSurfaces)) continue; // XXX if a discrete mass is associated with an LMPC element 
                                                                     // the the inertial force and tangent stiffness contribution 
                                                                     // must be defined in the DOF_FRM.
          }

          std::vector<int> eledofs(3);
          for(int j = 0; j < 3; ++j) {
            eledofs[j] = -1;
            for(int k = 0; k < allDOFs->num(iele); ++k)
              if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
          }

          if(eledofs[0] != -1 && eledofs[1] != -1 && eledofs[2] != -1) {
            // found an element with the 3 rotation dofs of current->node
            elemAdj[iele].dimass.push_back(std::pair<DMassData*,std::vector<int> >(current,eledofs));
            foundAdj = true;
            break;
          }
        }
        if(!foundAdj) std::cerr << " *** WARNING: could not find adjacent element for discrete mass at node " << current->node+1 << std::endl;
      }
      current = current->next;
    }
  }

  // 7. constrained rotations
  if(!solInfo().getNLInfo().linearelastic)
  for(int i = 0; i < numDirichlet; ++i) {
    if(dbc[i].dofnum == 3 || dbc[i].dofnum == 4 || dbc[i].dofnum == 5) {
      int dofs[3];
      c_dsa->number(dbc[i].nnum, DofSet::XYZrot, dofs);
      if((dbc[i].dofnum == 3 && dofs[1] != -1 && dofs[2] != -1) ||
         (dbc[i].dofnum == 4 && dofs[0] != -1 && dofs[2] != -1) ||
         (dbc[i].dofnum == 5 && dofs[0] != -1 && dofs[1] != -1)) {
        dsa->number(dbc[i].nnum, DofSet::XYZrot, dofs);
        for(int inode = 0; inode < nodeToElem->num(dbc[i].nnum); ++inode) { // loop over the elements attached to the node
                                                                            // at which the nodal force is applied
          int iele = (*nodeToElem)[dbc[i].nnum][inode];
          std::vector<int> eledofs(3);
          for(int j = 0; j < 3; ++j) {
            eledofs[j] = -1;
            for(int k = 0; k < allDOFs->num(iele); ++k)
              if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
          }
          if(eledofs[0] != -1 && eledofs[1] != -1 && eledofs[2] != -1) { 
            // element has all 3 rotation dofs of node dbc[i].nnum 
            elemAdj[iele].crot.insert(i);
            if(makeFollowedElemList) followedElemList.push_back(iele);
            break;
          }
        }
      }
    }
  }

  if(makeFollowedElemList) {
    std::sort(followedElemList.begin(), followedElemList.end());
    std::vector<int>::iterator it = std::unique(followedElemList.begin(), followedElemList.end());
    followedElemList.resize(it-followedElemList.begin());
  }
}

void
Domain::getWeightedStiffAndForceOnly(const std::map<int, double> &weights,
                                     GeomState &geomState, Vector& elementForce,
                                     Corotator **corotators, FullSquareMatrix *kel,
                                     Vector &residual, double lambda, double time,
                                     GeomState *refState, FullSquareMatrix *mel,
                                     FullSquareMatrix *kelCopy)
{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // MPC needs lambda for nonlinear statics
  const bool initialTime = (sinfo.isDynam() && time == 0 && solInfo().newmarkBeta != 0);

  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  bool compute_tangents = !initialTime;
  if(elemAdj.empty()) makeElementAdjacencyLists();

  if(!domain->solInfo().reduceFollower) {
    // Zero element stiffness. otherwise, the follower force tangents will accumulate
    for(std::vector<int>::iterator it = followedElemList.begin(), it_end = followedElemList.end(); it != it_end; ++it) kel[*it].zero();
  }
  
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
      getElemStiffAndForce(geomState, pseudoTime, refState, *corotators[iele], elementForce.data(), elementStiff);
      if(sinfo.newmarkBeta == 0) {
        corotators[iele]->updateStates(refState, geomState, nodes, sinfo.getTimeStep());
        // this is where handleElementDeletions would go
      }
      if(initialTime && packedEset[iele]->isConstraintElement() && packedEset[iele]->hasRot()) {
        // transform constraint jacobian and hessian to solve for the initial convected acceleration
        transformElemStiff(geomState, kel[iele], iele);
      }
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
                           (corotators[iele]), elementStiff, lambda, time, compute_tangents, conwep);
    }

    if(packedEset[iele]->hasRot() && !solInfo().getNLInfo().linearelastic) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iele, true);
    }
   
    // Apply lumping weight 
    const double lumpingWeight = it->second;
    elementForce *= lumpingWeight;
    elementStiff *= lumpingWeight;

    if(kelCopy)
     elementStiff -= kelCopy[iele];

    // Assemble element force into residual vector
    const int elemDofCount = elementStiff.dim();
    for(int iDof = 0; iDof < elemDofCount; ++iDof) {
      const int dofId = c_dsa->getRCN((*allDOFs)[iele][iDof]);
      if (dofId >= 0) {
        residual[dofId] -= elementForce[iDof];
        if(kelCopy){
          residual[dofId] += LinearElForce[iDof];
        }
      }
    }
  }

  if(!domain->solInfo().reduceFollower) {
    getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, NULL, compute_tangents);
  }

  if(sinfo.isDynam() && mel && !solInfo().getNLInfo().linearelastic && !solInfo().quasistatic) {
    getWeightedFictitiousForceOnly(weights, geomState, elementForce, kel, residual, time,
                                   refState, NULL, mel, compute_tangents);
  }

  if(mel && domain->solInfo().printMatLab && domain->solInfo().printMatLabExit) {
    filePrint(stderr, " ... Applying Weights to Mass Matrix...\n");
    for(std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
      const int iele = it->first;
      const double lumpingWeight = it->second;
      FullSquareMatrix &elementMass = mel[iele];
      elementMass *= lumpingWeight;
    }
  }

  if(!solInfo().getNLInfo().unsymmetric && solInfo().newmarkBeta != 0) {
    for(std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
      const int iele = it->first;
      kel[iele].symmetrize();
    }
    if(!domain->solInfo().reduceFollower) {
      for(std::vector<int>::iterator it = followedElemList.begin(), it_end = followedElemList.end(); it != it_end; ++it) {
        kel[*it].symmetrize();
      }
    }
  }
}


void
Domain::getUnassembledStiffAndForceOnly(const std::map<int, std::vector<int> > &weights,
                                        GeomState &geomState, Vector& elementForce,
                                        Corotator **corotators, FullSquareMatrix *kel,
                                        Vector &residual, int dispSize, double lambda, double time,
                                        GeomState *refState, FullSquareMatrix *mel,
                                        FullSquareMatrix *kelCopy)
{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // MPC needs lambda for nonlinear statics
  const bool initialTime = (sinfo.isDynam() && time == 0 && solInfo().newmarkBeta != 0);

  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  bool compute_tangents = !initialTime;
  if(elemAdj.empty()) makeElementAdjacencyLists();

  if(!domain->solInfo().reduceFollower) {
    // Zero element stiffness. otherwise, the follower force tangents will accumulate
    for(std::vector<int>::iterator it = followedElemList.begin(), it_end = followedElemList.end(); it != it_end; ++it) kel[*it].zero();
  }

  Vector LinearElForce(maxNumDOFs,0.0);
  Vector displacement(dispSize,0.0);
  if(kelCopy)                        
    geomState.get_tot_displacement(displacement,false);

  int uDofCounter = 0;
  for(std::map<int, std::vector<int> >::const_iterator mapit = weights.begin(), mapit_end = weights.end(); mapit != mapit_end; ++mapit) {
    const int iele = mapit->first;
    const std::vector<int> DOFvector(mapit->second);

    elementForce.zero();

    FullSquareMatrix &elementStiff = kel[iele];

    //option of remove linear component from nonlinear force, default is false 
    if(kelCopy) {
      LinearElForce.zero();
      getElemKtimesU(iele,elementStiff.dim(),displacement,LinearElForce.data(),kelCopy,(double *) dbg_alloca(sizeof(double)*maxNumDOFs*maxNumDOFs));
    }

    // Get updated tangent stiffness matrix and element internal force
    if (const Corotator *elementCorot = corotators[iele]) {
      getElemStiffAndForce(geomState, pseudoTime, refState, *elementCorot, elementForce.data(), elementStiff);
    }

    // Add configuration-dependent external forces and their element stiffness contributions
    if(domain->solInfo().reduceFollower) {
      getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
                           (corotators[iele]), elementStiff, lambda, time, compute_tangents, conwep);
    }


    if(packedEset[iele]->hasRot()) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iele, true);
    }

    if(kelCopy)
     elementStiff -= kelCopy[iele];

    // Assemble element force into residual vector
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
 //   getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, NULL, compute_tangents);
//  }

  if(sinfo.isDynam() && mel) {
    getUDEIMFictitiousForceOnly(weights, geomState, elementForce, kel, residual, time,
                                   refState, NULL, mel, compute_tangents);
  }

  if(!solInfo().getNLInfo().unsymmetric && solInfo().newmarkBeta != 0) {
    for(std::map<int, std::vector<int> >::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
      const int iele = it->first;
      kel[iele].symmetrize();
    }
    if(!domain->solInfo().reduceFollower) {
      for(std::vector<int>::iterator it = followedElemList.begin(), it_end = followedElemList.end(); it != it_end; ++it) {
        kel[*it].symmetrize();
      }
    }
  }
}

void
Domain::applyResidualCorrection(GeomState &geomState, Corotator **corotators, Vector &residual, double rcoef)
{
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) {
      Vector residualCorrection(packedEset[iele]->numDofs());
      residualCorrection.zero();
      corotators[iele]->getResidualCorrection(geomState, residualCorrection.data());
    
      // Assemble element residualCorrection into global residual vector
      for(int idof = 0; idof < packedEset[iele]->numDofs(); ++idof) {
        int dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
        if(dofNum >= 0)
          residual[dofNum] += rcoef*residualCorrection[idof];
      }
    }
  }
}

void
Domain::updateStates(GeomState *refState, GeomState &geomState, Corotator **corotators, double time)
{
  if(matrixTimers) matrixTimers->updateState -= getTime();
  double delt = std::max(sinfo.newmarkAlphaF,std::numeric_limits<double>::epsilon())*sinfo.getTimeStep();
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) {
      corotators[iele]->updateStates(refState, geomState, nodes, delt);
      if(sinfo.newmarkBeta != 0) handleElementDeletion(iele, geomState, time, *corotators[iele]);
    }
  }
  if(matrixTimers) matrixTimers->updateState += getTime();
}

void
Domain::updateWeightedElemStatesOnly(const std::map<int, double> &weights, GeomState *refState,
                                     GeomState &geomState, Corotator **corotators, double time)
{
  if(matrixTimers) matrixTimers->updateState -= getTime();
  double delt = std::max(sinfo.newmarkAlphaF,std::numeric_limits<double>::epsilon())*sinfo.getTimeStep();
  for(std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iele = it->first;
    if(corotators[iele]) corotators[iele]->updateStates(refState, geomState, nodes, delt);
  }
  if(matrixTimers) matrixTimers->updateState += getTime();
}

void
Domain::updateWeightedElemStatesOnly(const std::set<int> &weightedElems, GeomState *refState,
                                     GeomState &geomState, Corotator **corotators, double time)
{
  if(matrixTimers) matrixTimers->updateState -= getTime();
  double delt = std::max(sinfo.newmarkAlphaF,std::numeric_limits<double>::epsilon())*sinfo.getTimeStep();
  for(std::set<int>::const_iterator it = weightedElems.begin(), it_end = weightedElems.end(); it != it_end; ++it) {
    const int iele = *it;
    if(corotators[iele]) corotators[iele]->updateStates(refState, geomState, nodes, delt);
  }
  if(matrixTimers) matrixTimers->updateState += getTime();
}

void
Domain::initializeMultipliers(GeomState &geomState, Corotator **corotators)
{
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) corotators[iele]->initMultipliers(geomState);
  }
}

void
Domain::initializeParameters(bool flag)
{
  if(flag) geoSource->initializeParameters();
}

void
Domain::updateMultipliers(GeomState &geomState, Corotator **corotators)
{
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) corotators[iele]->updateMultipliers(geomState);
  }
}

void
Domain::updateParameters(bool flag)
{
  if(flag) geoSource->updateParameters();
}

double
Domain::getError(Corotator **corotators, GeomState &gs)
{
  // returns the largest constraint violation
  double err = 0;
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) err = std::max(err,corotators[iele]->getError(gs));
  }
  return err;
}

// used in nonlinear statics
void
Domain::createKelArray(FullSquareMatrix *&kArray)
{
  // Allocate array of pointers to FullSquareMatrix to store the
  // element stiffness matrices
  kArray = new FullSquareMatrix[numele];

  // Allocate the correct size for each elements stiffness matrix
  int iele;
  for(iele = 0; iele<numele; ++iele) {
    kArray[iele].setSize(packedEset[iele]->numDofs());
  }

  // Form and store element stiffness matrices into an array
  for(iele=0; iele<numele; ++iele) {
    kArray[iele].copy(packedEset[iele]->stiffness(nodes, kArray[iele].data()));
  }
}

// used in nonlinear dynamics
void
Domain::createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray)
{
 // Allocate array of pointers to FullSquareMatrix to store
 // the element stiffness matrices and element mass matrices
 kArray = new FullSquareMatrix[maxNumElements()];
 mArray = new FullSquareMatrix[maxNumElements()];

 // Allocate the correct size for each element's stiffness & mass matrix
 int iele;
 for(iele = 0; iele<numele; ++iele) {
   int dimension = packedEset[iele]->numDofs();
   kArray[iele].setSize(dimension);
   mArray[iele].setSize(dimension);
 }

 // Form and store element mass and stiffness matrices into an array
 for(iele=0; iele<numele; ++iele) {
   // note: only lumped mass matrix is supported currently for elements with rotation dofs in nonlinear dynamics
   //       (only the euler beam element is affected)
   double mratio = (packedEset[iele]->hasRot() && sinfo.isNonLin() && sinfo.isDynam()) ? 0 : geoSource->getMRatio();
   mArray[iele].copy(packedEset[iele]->massMatrix(nodes, mArray[iele].data(), mratio));
   if(domain->solInfo().elemLumpPodRom) {
     kArray[iele].zero();
   }
   else {
     kArray[iele].copy(packedEset[iele]->stiffness(nodes, kArray[iele].data()));
   }
 }

 // zero rotational degrees of freedom within element mass matrices
 // for nonlinear implicit dynamics
 if(sinfo.zeroRot && sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0) {
   int *dofType = dsa->makeDofTypeArray();

   int i,j;
   for(iele=0; iele<numele; ++iele) {
     for(i=0; i<mArray[iele].dim(); ++i)
       for(j=0; j<mArray[iele].dim(); ++j)
          if( dofType[ (*allDOFs)[iele][i] ] == 1 ||
              dofType[ (*allDOFs)[iele][j] ] == 1) {
             mArray[iele][i][j] = 0.0;
          }
   }
 }

 if(sinfo.isNonLin() && !domain->solInfo().basicDofCoords) {
   for(iele=0; iele<numele; ++iele) transformMatrix(mArray[iele], iele);
 }
}

// used in nonlinear dynamics
void
Domain::createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray, FullSquareMatrix *&cArray)
{
 // Allocate array of pointers to FullSquareMatrix to store
 // the element stiffness matrices, element damping matrices and element mass matrices
 kArray = new FullSquareMatrix[numele];
 mArray = new FullSquareMatrix[numele];
 cArray = new FullSquareMatrix[numele];

 // Allocate the correct size for each elements stiffness & mass matrix
 int iele;
 for(iele = 0; iele<numele; ++iele) {
   int dimension = packedEset[iele]->numDofs();
   kArray[iele].setSize(dimension);
   mArray[iele].setSize(dimension);
   cArray[iele].setSize(dimension);
 }

 // Form and store element damping, mass and stiffness matrices into arrays
 for(iele=0; iele<numele; ++iele) {
   // note: only lumped mass matrix is supported currently for elements with rotation dofs in nonlinear dynamics
   //       (only the euler beam element is affected)
   double mratio = (packedEset[iele]->hasRot() && sinfo.isNonLin() && sinfo.isDynam()) ? 0 : geoSource->getMRatio();
   mArray[iele].copy(packedEset[iele]->massMatrix(nodes, mArray[iele].data(), mratio));
   cArray[iele].copy(packedEset[iele]->dampingMatrix(nodes, cArray[iele].data()));
   if(domain->solInfo().elemLumpPodRom) {
     kArray[iele].zero();
   }
   else {
     kArray[iele].copy(packedEset[iele]->stiffness(nodes, kArray[iele].data()));
   }
 }

 // add Rayleigh damping
 int i,j;
 double alpha, beta;
 for(iele=0; iele<numele; ++iele) {
   if(!packedEset[iele]->getProperty()) continue; // phantom
   if(packedEset[iele]->isConstraintElement()) cArray[iele].zero();
   else {
     alpha = (packedEset[iele]->isDamped()) ? packedEset[iele]->getProperty()->alphaDamp : sinfo.alphaDamp;
     beta  = (packedEset[iele]->isDamped()) ? packedEset[iele]->getProperty()->betaDamp : sinfo.betaDamp;
     for(i=0; i<cArray[iele].dim(); ++i)
       for(j=0; j<cArray[iele].dim(); ++j)
         cArray[iele][i][j] += alpha*mArray[iele][i][j] + beta*kArray[iele][i][j];
   }
 }

 // zero rotational degrees of freedom within element mass matrices and damping matrices
 // for nonlinear implicit dynamics
 if(sinfo.zeroRot && sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0) {
   int *dofType = dsa->makeDofTypeArray();
   for(iele=0; iele<numele; ++iele) {
     for(i=0; i<mArray[iele].dim(); ++i)
       for(j=0; j<mArray[iele].dim(); ++j)
          if( dofType[ (*allDOFs)[iele][i] ] == 1 ||
              dofType[ (*allDOFs)[iele][j] ] == 1) {
             mArray[iele][i][j] = 0.0;
             cArray[iele][i][j] = 0.0;
         }
   }
 }

 if(sinfo.isNonLin() && !domain->solInfo().basicDofCoords) {
   for(iele=0; iele<numele; ++iele) {
     transformMatrix(mArray[iele], iele);
     transformMatrix(cArray[iele], iele);
   }
 }
}

bool
Domain::reactionsReqd(double time, int step)
{
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for(int iInfo = 0; iInfo < numOutInfo; ++iInfo) {
    if((oinfo[iInfo].type == OutputInfo::Reactions || oinfo[iInfo].type == OutputInfo::Reactions6)
       && (step%oinfo[iInfo].interval == 0 || time == 0.0)) return true;
  }
  return false;
}

void
Domain::postProcessing(GeomState *geomState, Vector& force, Vector &aeroForce,
                       double time, int step, double* velocity, double *vcx,
                       Corotator **allCorot, double *acceleration, double *acx,
                       GeomState *refState, Vector *reactions, SparseMatrix *M,
                       SparseMatrix *C)
{
  if(time == sinfo.initialTime) {
    geoSource->openOutputFiles();
    firstOutput = false;
  }

  if(sinfo.nRestart > 0 && velocity != 0) {
    StackVector v_n(velocity, numUncon());
    StackVector a_n(acceleration, numUncon());
    writeRestartFile(time, step, v_n, a_n, geomState);
  }

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for(int iInfo = 0; iInfo < numOutInfo; ++iInfo)
  {
    if(oinfo[iInfo].sentype > 0) continue;
    postProcessingImpl(iInfo, geomState, force, aeroForce, time, step, velocity, vcx,
                       allCorot, acceleration, acx, refState, reactions, M, C);
  }

}

void
Domain::postProcessingImpl(int iInfo, GeomState *geomState, Vector& force, Vector &aeroForce,
                           double time, int step, double* velocity, double *vcx,
                           Corotator **allCorot, double *acceleration, double *acx,
                           GeomState *refState, Vector *reactions, SparseMatrix *M,
                           SparseMatrix *C)
{
 if(outFlag && !nodeTable) makeNodeTable(outFlag);
 int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes

 enum {SXX=0,SYY=1,SZZ=2,SXY= 3,SYZ= 4,SXZ= 5,VON=6,
       EXX=7,EYY=8,EZZ=9,EXY=10,EYZ=11,EXZ=12,STRAINVON=13,
       VONTOP=14,VONBOT=15};

 enum {INX,INY,INZ,AXM,AYM,AZM};

 enum {PSTRESS1=0,PSTRESS2=1,PSTRESS3=2,
       PSTRAIN1=3,PSTRAIN2=4,PSTRAIN3=5};

  int i, nodeI, realNode; 
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // Check output interval
  if(step%oinfo[iInfo].interval != 0 && time != 0.0) return;

  int w = oinfo[iInfo].width;
  int p = oinfo[iInfo].precision;

  int first_node, last_node;
  if (oinfo[iInfo].nodeNumber == -1) {
    first_node=0;
    last_node=numNodes;
  }
  else  {
    first_node=oinfo[iInfo].nodeNumber;
    last_node=first_node+1;
  }
  int nNodes = last_node - first_node;
  int nPrintNodes = (outFlag && oinfo[iInfo].nodeNumber == -1) ? nodes.nnz() : nNodes;
  int angularintype = (domain->solInfo().galerkinPodRom) ? 2 : 1;
  bool rescalein = (domain->solInfo().galerkinPodRom) ? false : true;

  switch(oinfo[iInfo].type) {
    case OutputInfo::Displacement:  {
      double (*data)[3] = new double[nPrintNodes][3];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        data[nodeI][0] = (nodes[iNode] && iNode<geomState->numNodes()) ? (*geomState)[iNode].x-nodes[iNode]->x : 0;
        data[nodeI][1] = (nodes[iNode] && iNode<geomState->numNodes()) ? (*geomState)[iNode].y-nodes[iNode]->y : 0;
        data[nodeI][2] = (nodes[iNode] && iNode<geomState->numNodes()) ? (*geomState)[iNode].z-nodes[iNode]->z : 0;
        if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, false);
      }

      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Temperature:  {
      if(sinfo.order == 1 || sinfo.soltyp == 2) { // thermal analysis
        double *data = new double[nPrintNodes];
        for (i = 0, realNode = -1; i < nNodes; ++i) {
          int iNode = first_node+i;
          if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
          data[nodeI] = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].x : 0;
        }
        geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
        delete [] data;
      }
      else { // structural analysis with thermal load
        double *nodalTemperatures = 0;
        // Either get the nodal temperatures from the input file or from the thermal model
        if(sinfo.thermoeFlag >= 0) nodalTemperatures = temprcvd;
        else if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
        double *data = new double[nPrintNodes];
        for (i = 0, realNode = -1; i < nNodes; ++i) {
          int iNode = first_node+i;
          if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
          data[nodeI] = (nodalTemperatures) ? nodalTemperatures[iNode] : 0;
        }
        geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
        delete [] data;
      }
    } 
      break;
    case OutputInfo::Disp6DOF:  {
      double (*data)[6] = new double[nPrintNodes][6];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          data[nodeI][0] = (*geomState)[iNode].x - nodes[iNode]->x;
          data[nodeI][1] = (*geomState)[iNode].y - nodes[iNode]->y;
          data[nodeI][2] = (*geomState)[iNode].z - nodes[iNode]->z;
          tran_rvec((*geomState)[iNode].R, (*geomState)[iNode].theta, oinfo[iInfo].rescaling,
                    oinfo[iInfo].rotvecouttype, &data[nodeI][3]);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, true);
        }
        else {
          std::fill_n(&data[nodeI][0], 6, 0.0);
        }
      }
      geoSource->outputNodeVectors6(iInfo, data, nPrintNodes, time);
      delete [] data;
    } 
      break;
    case OutputInfo::RotationVector:  {
      double (*data)[3] = new double[nPrintNodes][3];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          tran_rvec((*geomState)[iNode].R, (*geomState)[iNode].theta, oinfo[iInfo].rescaling,
                    oinfo[iInfo].rotvecouttype, &data[nodeI][0]);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, false);
        }
        else {
          std::fill_n(&data[nodeI][0], 3, 0.0);
        }
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotationMatrix:  {
      double (*data)[9] = new double[nPrintNodes][9];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          data[nodeI][0] = (*geomState)[iNode].R[0][0];
          data[nodeI][1] = (*geomState)[iNode].R[0][1];
          data[nodeI][2] = (*geomState)[iNode].R[0][2];
          data[nodeI][3] = (*geomState)[iNode].R[1][0];
          data[nodeI][4] = (*geomState)[iNode].R[1][1];
          data[nodeI][5] = (*geomState)[iNode].R[1][2];
          data[nodeI][6] = (*geomState)[iNode].R[2][0];
          data[nodeI][7] = (*geomState)[iNode].R[2][1];
          data[nodeI][8] = (*geomState)[iNode].R[2][2];
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformMatrix(data[nodeI], iNode, false);
        }
        else {
          std::fill_n(&data[nodeI][0], 9, 0.0);
        }
      }
      geoSource->outputNodeVectors9(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Quaternion:  {
      double (*data)[4] = new double[nPrintNodes][4];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          if(oinfo[iInfo].oframe == OutputInfo::Local) {
            double R[3][3] = { { (*geomState)[iNode].R[0][0], (*geomState)[iNode].R[0][1], (*geomState)[iNode].R[0][2] },
                               { (*geomState)[iNode].R[1][0], (*geomState)[iNode].R[1][1], (*geomState)[iNode].R[1][2] },
                               { (*geomState)[iNode].R[2][0], (*geomState)[iNode].R[2][1], (*geomState)[iNode].R[2][2] } };
            transformMatrix(&(R[0][0]), iNode, false); 
            mat_to_quat(R, data[nodeI]);
          }
          else mat_to_quat((*geomState)[iNode].R, data[nodeI]);
        }
        else {
          std::fill_n(&data[nodeI][0], 4, 0.0);
        }
      }
      geoSource->outputNodeVectors4(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Velocity: {
      double (*data)[3] = new double[nPrintNodes][3];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          for(int j=0; j<3; ++j) data[nodeI][j] = (*geomState)[iNode].v[j];
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, false);
        }
        else {
          std::fill_n(&data[nodeI][0], 3, 0.0);
        }
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::TemperatureFirstTimeDerivative:  {
      if(sinfo.order == 1 || sinfo.soltyp == 2) {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        data[nodeI] = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].v[0] : 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
      } else fprintf(stderr," *** WARNING: Output case %d not implemented for non-linear direct solver \n", iInfo);
    }
      break;
    case OutputInfo::Velocity6: {
      double (*data)[6] = new double[nPrintNodes][6];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          for(int j=0; j<3; ++j) data[nodeI][j] = (*geomState)[iNode].v[j];
          tran_veloc((*geomState)[iNode].R, (*geomState)[iNode].theta, &(*geomState)[iNode].v[3], angularintype,
                     oinfo[iInfo].angularouttype, rescalein, oinfo[iInfo].rescaling, &data[nodeI][3]);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, true);
        }
        else {
          std::fill_n(&data[nodeI][0], 6, 0.0);
        }
      }
      geoSource->outputNodeVectors6(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::AngularVelocity: {
      double (*data)[3] = new double[nPrintNodes][3];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          tran_veloc((*geomState)[iNode].R, (*geomState)[iNode].theta, &(*geomState)[iNode].v[3], angularintype,
                     oinfo[iInfo].angularouttype, rescalein, oinfo[iInfo].rescaling, &data[nodeI][0]);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, false);
        }
        else {
          std::fill_n(&data[nodeI][0], 3, 0.0);
        }
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Acceleration:  {
      double (*data)[3] = new double[nPrintNodes][3];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          for(int j=0; j<3; ++j) data[nodeI][j] = (*geomState)[iNode].a[j];
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, false);
        }
        else {
          std::fill_n(&data[nodeI][0], 3, 0.0);
        }
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Accel6:  {
      double (*data)[6] = new double[nPrintNodes][6];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          for(int j=0; j<3; ++j) data[nodeI][j] = (*geomState)[iNode].a[j];
          tran_accel((*geomState)[iNode].R, (*geomState)[iNode].theta, &(*geomState)[iNode].v[3], &(*geomState)[iNode].a[3],
                     angularintype, oinfo[iInfo].angularouttype, rescalein, oinfo[iInfo].rescaling, &data[nodeI][3]);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, true);
        }
        else {
          std::fill_n(&data[nodeI][0], 6, 0.0);
        }
      }
      geoSource->outputNodeVectors6(iInfo, data, nPrintNodes, time);
      delete [] data;
    } 
      break;
    case OutputInfo::AngularAcceleration:  {
      double (*data)[3] = new double[nPrintNodes][3];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          tran_accel((*geomState)[iNode].R, (*geomState)[iNode].theta, &(*geomState)[iNode].v[3], &(*geomState)[iNode].a[3],
                     angularintype, oinfo[iInfo].angularouttype, rescalein, oinfo[iInfo].rescaling, &data[nodeI][0]);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(data[nodeI], iNode, false);
        }
        else {
          std::fill_n(&data[nodeI][0], 3, 0.0);
        }
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispX:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if(nodes[iNode] && iNode < geomState->numNodes()) {
          if(oinfo[iInfo].oframe == OutputInfo::Local) {
            double d[3] = { (*geomState)[iNode].x - nodes[iNode]->x,
                            (*geomState)[iNode].y - nodes[iNode]->y,
                            (*geomState)[iNode].z - nodes[iNode]->z };
            transformVector(d, iNode, false);
            data[nodeI] = d[0];
          }
          else {
            data[nodeI] = (*geomState)[iNode].x - nodes[iNode]->x;
          }
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispY:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if(nodes[iNode] && iNode < geomState->numNodes()) {
          if(oinfo[iInfo].oframe == OutputInfo::Local) {
            double d[3] = { (*geomState)[iNode].x - nodes[iNode]->x,
                            (*geomState)[iNode].y - nodes[iNode]->y,
                            (*geomState)[iNode].z - nodes[iNode]->z };
            transformVector(d, iNode, false);
            data[nodeI] = d[1];
          }
          else {
            data[nodeI] = (*geomState)[iNode].y - nodes[iNode]->y;
          }
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispZ:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if(nodes[iNode] && iNode < geomState->numNodes()) {
          if(oinfo[iInfo].oframe == OutputInfo::Local) {
            double d[3] = { (*geomState)[iNode].x - nodes[iNode]->x,
                            (*geomState)[iNode].y - nodes[iNode]->y,
                            (*geomState)[iNode].z - nodes[iNode]->z };
            transformVector(d, iNode, false);
            data[nodeI] = d[2];
          }
          else {
            data[nodeI] = (*geomState)[iNode].z - nodes[iNode]->z;
          }
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotX:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          double rot[3];
          tran_rvec((*geomState)[iNode].R, (*geomState)[iNode].theta, oinfo[iInfo].rescaling,
                    oinfo[iInfo].rotvecouttype, rot);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(rot, iNode, false);
          data[nodeI] = rot[0];
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotY:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          double rot[3];
          tran_rvec((*geomState)[iNode].R, (*geomState)[iNode].theta, oinfo[iInfo].rescaling,
                    oinfo[iInfo].rotvecouttype, rot);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(rot, iNode, false);
          data[nodeI] = rot[1];
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotZ:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          double rot[3];
          tran_rvec((*geomState)[iNode].R, (*geomState)[iNode].theta, oinfo[iInfo].rescaling,
                    oinfo[iInfo].rotvecouttype, rot);
          if(oinfo[iInfo].oframe == OutputInfo::Local) transformVector(rot, iNode, false);
          data[nodeI] = rot[2];
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispMod:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        double x = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].x - nodes[iNode]->x : 0;
        double y = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].y - nodes[iNode]->y : 0;
        double z = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].z - nodes[iNode]->z : 0;
        data[nodeI] = sqrt(x*x+y*y+z*z);
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotMod:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          double rot[3];
          tran_rvec((*geomState)[iNode].R, (*geomState)[iNode].theta, oinfo[iInfo].rescaling,
                    oinfo[iInfo].rotvecouttype, rot);
          data[nodeI] = sqrt(rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::TotMod:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          double x = (*geomState)[iNode].x - nodes[iNode]->x;
          double y = (*geomState)[iNode].y - nodes[iNode]->y;
          double z = (*geomState)[iNode].z - nodes[iNode]->z;
          double rot[3];
          tran_rvec((*geomState)[iNode].R, (*geomState)[iNode].theta, oinfo[iInfo].rescaling,
                    oinfo[iInfo].rotvecouttype, rot);
          data[nodeI] = sqrt(x*x+y*y+z*z+rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
        }
        else {
          data[nodeI] = 0;
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::StressXX:
      getStressStrain(*geomState, allCorot,  iInfo, SXX, time, refState);
      break;
    case OutputInfo::StressYY:
      getStressStrain(*geomState, allCorot,  iInfo, SYY, time, refState);
      break;
    case OutputInfo::StressZZ:
      getStressStrain(*geomState, allCorot,  iInfo, SZZ, time, refState);
      break;
    case OutputInfo::StressXY:
      getStressStrain(*geomState, allCorot,  iInfo, SXY, time, refState);
      break;
    case OutputInfo::StressYZ:
      getStressStrain(*geomState, allCorot,  iInfo, SYZ, time, refState);
      break;
    case OutputInfo::StressXZ:
      getStressStrain(*geomState, allCorot,  iInfo, SXZ, time, refState);
      break;
    case OutputInfo::StrainXX:
      getStressStrain(*geomState, allCorot,  iInfo, EXX, time, refState);
      break;
    case OutputInfo::StrainYY:
      getStressStrain(*geomState, allCorot,  iInfo, EYY, time, refState);
      break;
    case OutputInfo::StrainZZ:
      getStressStrain(*geomState, allCorot,  iInfo, EZZ, time, refState);
      break;
    case OutputInfo::StrainXY:
      getStressStrain(*geomState, allCorot,  iInfo, EXY, time, refState);
      break;
    case OutputInfo::StrainYZ:
      getStressStrain(*geomState, allCorot,  iInfo, EYZ, time, refState);
      break;
    case OutputInfo::StrainXZ:
      getStressStrain(*geomState, allCorot,  iInfo, EXZ, time, refState);
      break;
    case OutputInfo::StressVM:
      getStressStrain(*geomState, allCorot,  iInfo, VON, time, refState);
      break;
    case OutputInfo::AggrStVM:
      getStressStrain(*geomState, allCorot,  iInfo, AGGREGATEDVON, time, refState);
      break;
    case OutputInfo::StrainVM:
      getStressStrain(*geomState, allCorot,  iInfo, STRAINVON, time, refState);
      break;
    case OutputInfo::StressPR1:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS1, time);
      break;
    case OutputInfo::StressPR2:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS2, time);
      break;
    case OutputInfo::StressPR3:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS3, time);
      break;
    case OutputInfo::StrainPR1:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN1, time);
      break;
    case OutputInfo::StrainPR2:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN2, time);
      break;
    case OutputInfo::StrainPR3:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN3, time);
      break;
    case OutputInfo::Damage:
      getStressStrain(*geomState, allCorot,  iInfo, DAMAGE, time, refState);
      break;
    case OutputInfo::EquivalentPlasticStrain:
      getStressStrain(*geomState, allCorot,  iInfo, EQPLSTRN, time, refState);
      break;
    case OutputInfo::BackStressXX:
      getStressStrain(*geomState, allCorot,  iInfo, BACKSXX, time, refState);
      break;
    case OutputInfo::BackStressYY:
      getStressStrain(*geomState, allCorot,  iInfo, BACKSYY, time, refState);
      break;
    case OutputInfo::BackStressZZ:
      getStressStrain(*geomState, allCorot,  iInfo, BACKSZZ, time, refState);
      break;
    case OutputInfo::BackStressXY:
      getStressStrain(*geomState, allCorot,  iInfo, BACKSXY, time, refState);
      break;
    case OutputInfo::BackStressYZ:
      getStressStrain(*geomState, allCorot,  iInfo, BACKSYZ, time, refState);
      break;
    case OutputInfo::BackStressXZ:
      getStressStrain(*geomState, allCorot,  iInfo, BACKSXZ, time, refState);
      break;
    case OutputInfo::PlasticStrainXX:
      getStressStrain(*geomState, allCorot,  iInfo, PLASTICEXX, time, refState);
      break;
    case OutputInfo::PlasticStrainYY:
      getStressStrain(*geomState, allCorot,  iInfo, PLASTICEYY, time, refState);
      break;
    case OutputInfo::PlasticStrainZZ:
      getStressStrain(*geomState, allCorot,  iInfo, PLASTICEZZ, time, refState);
      break;
    case OutputInfo::PlasticStrainXY:
      getStressStrain(*geomState, allCorot,  iInfo, PLASTICEXY, time, refState);
      break;
    case OutputInfo::PlasticStrainYZ:
      getStressStrain(*geomState, allCorot,  iInfo, PLASTICEYZ, time, refState);
      break;
    case OutputInfo::PlasticStrainXZ:
      getStressStrain(*geomState, allCorot,  iInfo, PLASTICEXZ, time, refState);
      break;
    case OutputInfo::YModulus:
      getElementAttr(iInfo, YOUNG, time);
      break;
    case OutputInfo::MDensity:
      getElementAttr(iInfo, MDENS, time);
      break;
    case OutputInfo::Thicknes:
      getElementAttr(iInfo, THICK, time);
      break;
    case OutputInfo::InXForce:
      getElementForces(*geomState, allCorot, iInfo, INX, time);
      break;
    case OutputInfo::InYForce:
      getElementForces(*geomState, allCorot, iInfo, INY, time);
      break;
    case OutputInfo::InZForce:
      getElementForces(*geomState, allCorot, iInfo, INZ, time);
      break;
    case OutputInfo::AXMoment:
      getElementForces(*geomState, allCorot, iInfo, AXM, time);
      break;
    case OutputInfo::AYMoment:
      getElementForces(*geomState, allCorot, iInfo, AYM, time);
      break;
    case OutputInfo::AZMoment:
      getElementForces(*geomState, allCorot, iInfo, AZM, time);
      break;
    case OutputInfo::Energies: {
      double Wkin, Wela, Wdis, error;
      computeEnergies(geomState, force, time, &aeroForce, velocity, allCorot, M, C, Wela, Wkin, Wdis, error);
      geoSource->outputEnergies(iInfo, time, Wext, Waero, Wela, Wkin, Wdmp+Wdis, error);
    } break;
    case OutputInfo::DissipatedEnergy: {
      int groupIndex = geoSource->getGroupNumber(iInfo);
      if(groupIndex<0) {
        double D = getDissipatedEnergy(geomState, allCorot);
        geoSource->outputEnergy(iInfo, time, D);
      }
      else {
        std::map<int, Group> &group = geoSource->group;
        int numAttribs = group[groupIndex].attributes.size();
        double* D = getDissipatedEnergy(geomState, allCorot, groupIndex);
        geoSource->outputEnergyPerAttribute(iInfo, time, D, numAttribs);
        delete [] D;
      }
    } break;
    case OutputInfo::AeroForce: break; // this is done in FlExchange.C
    case OutputInfo::AeroXForce:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        if(oinfo[iInfo].oframe == OutputInfo::Local) {
          int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
          data[nodeI]  = (xloc >= 0) ? aeroForce[xloc] : 0.0;
        }
        else {
          int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
          int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
          int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
          double f[3] = { (xloc >= 0) ? aeroForce[xloc] : 0.0,
                          (yloc >= 0) ? aeroForce[yloc] : 0.0,
                          (zloc >= 0) ? aeroForce[zloc] : 0.0 };
          transformVectorInv(f, first_node+iNode, false);
          data[nodeI] = f[0];
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroYForce:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        if(oinfo[iInfo].oframe == OutputInfo::Local) {
          int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
          data[nodeI]  = (yloc >= 0) ? aeroForce[yloc] : 0.0;
        }
        else {
          int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
          int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
          int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
          double f[3] = { (xloc >= 0) ? aeroForce[xloc] : 0.0,
                          (yloc >= 0) ? aeroForce[yloc] : 0.0,
                          (zloc >= 0) ? aeroForce[zloc] : 0.0 };
          transformVectorInv(f, first_node+iNode, false);
          data[nodeI] = f[1];
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroZForce:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        if(oinfo[iInfo].oframe == OutputInfo::Local) {
          int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
          data[nodeI] = (zloc >= 0) ? aeroForce[zloc] : 0.0;
        }
        else {
          int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
          int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
          int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
          double f[3] = { (xloc >= 0) ? aeroForce[xloc] : 0.0,
                          (yloc >= 0) ? aeroForce[yloc] : 0.0,
                          (zloc >= 0) ? aeroForce[zloc] : 0.0 };
          transformVectorInv(f, first_node+iNode, false);
          data[nodeI] = f[2];
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroXMom:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        if(oinfo[iInfo].oframe == OutputInfo::Local) {
          int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
          data[nodeI] = (xrot >= 0) ? aeroForce[xrot] : 0.0;
        }
        else {
          int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
          int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
          int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
          double m[3] = { (xrot >= 0) ? aeroForce[xrot] : 0.0,
                          (yrot >= 0) ? aeroForce[yrot] : 0.0,
                          (zrot >= 0) ? aeroForce[zrot] : 0.0 };
          transformVectorInv(m, first_node+iNode, false);
          data[nodeI] = m[0];
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroYMom:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        if(oinfo[iInfo].oframe == OutputInfo::Local) {
          int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
          data[nodeI] = (yrot >= 0) ? aeroForce[yrot] : 0.0;
        }
        else {
          int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
          int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
          int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
          double m[3] = { (xrot >= 0) ? aeroForce[xrot] : 0.0,
                          (yrot >= 0) ? aeroForce[yrot] : 0.0,
                          (zrot >= 0) ? aeroForce[zrot] : 0.0 };
          transformVectorInv(m, first_node+iNode, false);
          data[nodeI] = m[1];
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroZMom:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        if(oinfo[iInfo].oframe == OutputInfo::Local) {
          int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
          data[nodeI] = (zrot >= 0) ? aeroForce[zrot] : 0.0;
        }
        else {
          int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
          int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
          int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
          double m[3] = { (xrot >= 0) ? aeroForce[xrot] : 0.0,
                          (yrot >= 0) ? aeroForce[yrot] : 0.0,
                          (zrot >= 0) ? aeroForce[zrot] : 0.0 };
          transformVectorInv(m, first_node+iNode, false);
          data[nodeI] = m[2];
        }
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::Reactions: {
      if(!reactions) break;
      double (*rxyz)[3] = new double[nPrintNodes][3];
      DofSet dofs[3] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp };
      for(int iNode = 0, realNode = -1; iNode < nNodes; ++iNode) {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        for(int k = 0; k < 3; ++k) {
          int dof =   dsa->locate(first_node+iNode, dofs[k].list());
          int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
          rxyz[nodeI][k] = (cdof >= 0) ? (*reactions)[cdof] : 0;     // constrained
        }
        if(oinfo[iInfo].oframe == OutputInfo::Global) transformVectorInv(rxyz[nodeI], first_node+iNode, false);
      }
      geoSource->outputNodeVectors(iInfo, rxyz, nPrintNodes, time);
      delete [] rxyz;
    } break;
    case OutputInfo::Reactions6: {
      if(!reactions) break;
      double (*rxyz)[6] = new double[nPrintNodes][6];
      DofSet dofs[6] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
                         DofSet::Xrot, DofSet::Yrot, DofSet::Zrot };
      for(int iNode = 0, realNode = -1; iNode < nNodes; ++iNode) {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        for(int k = 0; k < 6; ++k) {
          int dof =   dsa->locate(first_node+iNode, dofs[k].list());
          int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
          rxyz[nodeI][k] = (cdof >= 0) ? (*reactions)[cdof] : 0;     // constrained
        }
        if(oinfo[iInfo].oframe == OutputInfo::Global) transformVectorInv(rxyz[nodeI], first_node+iNode, true);
      }
      geoSource->outputNodeVectors6(iInfo, rxyz, nPrintNodes, time);
      delete [] rxyz;
    } break;
    case OutputInfo::DeletedElements: {
      for(std::vector<std::pair<double,int> >::iterator it = outDeletedElements.begin(); it != outDeletedElements.end(); ++it) {
        filePrint(oinfo[iInfo].filptr, " %12.6e  %9d          Undetermined\n", it->first, it->second+1);
        fflush(oinfo[iInfo].filptr);
      }
      outDeletedElements.clear();
    } break;
     case OutputInfo::Statevector:
        break;
     case OutputInfo::Velocvector:
        break;
     case OutputInfo::Accelvector:
        break;
     case OutputInfo::InternalStateVar:
        break;
     case OutputInfo::DualStateVar:
        break;
     case OutputInfo::MuStateVar:
        break;
     case OutputInfo::Residual:
        break;
     case OutputInfo::Jacobian:
        break;
     case OutputInfo::RobData:
        break;
     case OutputInfo::SampleMesh:
        break;
     case OutputInfo::DispThic:
     case OutputInfo::WeigThic:
     case OutputInfo::VMstThic:
     case OutputInfo::AGstThic:
        break;

    default:
      fprintf(stderr," *** WARNING: Output case %d not implemented for non-linear direct solver \n", iInfo);
      break;
  }

}

void
Domain::createCorotators(Corotator **allCorot)
{
 // Allocate memory for element stiffness matrix
 double *kmatrix = new double[maxNumDOFs*maxNumDOFs];

 // Loop over elements to get their Corotators
 for(int iele = 0; iele < numele; ++iele)
   allCorot[iele] = packedEset[iele]->getCorotator(nodes, kmatrix,
                                      sinfo.getNLInfo().fitAlgShell,
                                      sinfo.getNLInfo().fitAlgBeam);

 delete [] kmatrix;
}

void
Domain::createContactCorotators(Corotator **allCorot, FullSquareMatrix *kArray, FullSquareMatrix *mArray)
{
 // Allocate memory for element stiffness matrix
 double *kmatrix = new double[maxNumDOFs*maxNumDOFs];

 // Loop over elements to get their Corotators
 for(int iele = numele-contactSurfElems.size(); iele < numele; ++iele) {
   allCorot[iele] = packedEset[iele]->getCorotator(nodes, kmatrix,
                                      sinfo.getNLInfo().fitAlgShell,
                                      sinfo.getNLInfo().fitAlgBeam);

   int dimension = packedEset[iele]->numDofs();
   kArray[iele].setSize(dimension);
   mArray[iele].setSize(dimension);
 }

 delete [] kmatrix;
}

void
Domain::getGeometricStiffness(GeomState &geomState, Vector& elementInternalForce,
                              Corotator **allCorot, FullSquareMatrix *&geomKelArray)
{

   // Get Geometric Stiffness

   int iele;
   for(iele = 0; iele < numele; ++iele) {

      // Get updated tangent stiffness matrix and element internal force
      allCorot[iele]->formGeometricStiffness(geomState, nodes,
                                             geomKelArray[iele],
                                             elementInternalForce.data());
      if(!solInfo().getNLInfo().unsymmetric) geomKelArray[iele].symmetrize();
   }

   if(domain->pressureFlag()) {
     for(iele = 0; iele < numele;  ++iele) {
       // If there is no pressure defined, skip the element
       if(packedEset[iele]->getPressure() == NULL) continue;
 
       // Compute (linear) element pressure force in the local coordinates
       elementInternalForce.zero();
       packedEset[iele]->computePressureForce(nodes, elementInternalForce, &geomState, 1);
       elementInternalForce *= -1.0; // due to IDISP6 -1

       // Include the "load stiffness matrix" in kel[iele]
       allCorot[iele]->getDExternalForceDu(geomState, nodes, geomKelArray[iele],
                                           elementInternalForce.data());
     }
   }

}

void
Domain::computeGeometricPreStress(Corotator **&allCorot, GeomState *&geomState,
                                  FullSquareMatrix *&kelArray, StaticTimers *times,
                                  FullSquareMatrix *&geomKelArray, FullSquareMatrix *&melArray,
                                  bool melFlag)
{
   // ... ALLOCATE MEMORY FOR THE ARRAY OF COROTATORS
   times->corotatorTime -= getTime();
   allCorot = new Corotator *[maxNumElements()];

   // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
   createCorotators(allCorot);
   times->corotatorTime += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Create Element Corotators %19.5f s\n", times->corotatorTime/1000.0);
#endif

   // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
   times->kelArrayTime -= getTime();
   if(melFlag) createKelArray(kelArray, melArray);
   else createKelArray(kelArray);
   times->kelArrayTime += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Create Element Stiffness Array %14.5f s\n", times->kelArrayTime/1000.0);
#endif

   times->timeGeom -= getTime();
   if(domain->solInfo().soltyp == 2)
     geomState =  new TemperatureState(*getDSA(), *getCDSA(), getNodes());
   else 
     geomState = new GeomState(*getDSA(), *getCDSA(), getNodes(), &getElementSet(), getNodalTemperatures());
   times->timeGeom += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Build GeomState %29.5f s\n", times->timeGeom/1000.0);
#endif

   // geometric prestress for linear: update geomState with initial displacements and compute element stiffness matrices
   if(!sinfo.isNonLin() && (numInitDisp6() > 0) && (sinfo.gepsFlg == 1)) {
     times->timePresc -= getTime();
     geomState->updatePrescribedDisplacement(getInitDisp6(), numInitDisp6(), 1.0);
     times->timePresc += getTime();
#ifdef PRINT_NLTIMERS
     fprintf(stderr," ... Update Geometry to initial state %12.5f s\n", times->timePresc/1000.0);
#endif
   }

// if(!domain->solInfo().samplingPodRom && !domain->solInfo().DEIMBasisPod && !domain->solInfo().UDEIMBasisPod) {
 if(domain->solInfo().stable && !domain->solInfo().samplingPodRom && !domain->solInfo().DEIMBasisPod && !domain->solInfo().UDEIMBasisPod) {
     times->buildStiffAndForce -= getTime();
     Vector elementInternalForce(maxNumDOF(), 0.0);
     Vector residual(numUncon(), 0.0);
     getStiffAndForce(*geomState, elementInternalForce, allCorot,
                      kelArray, residual, 1.0, 0.0, geomState, (Vector*) NULL, melArray);
     times->buildStiffAndForce += getTime();
#ifdef PRINT_NLTIMERS
     fprintf(stderr," ... Build Element Tangent Stiffness %13.5f s\n",
             times->buildStiffAndForce/1000.0);
#endif
   }

   // Buckling analysis if requested
   if(geomKelArray) {
     Vector elementInternalForce(maxNumDOF(), 0.0);
     getGeometricStiffness(*geomState, elementInternalForce, allCorot, geomKelArray);
   }
}


void
Domain::getStressStrain(GeomState &geomState, Corotator **allCorot,
                        int fileNumber, int stressIndex, double time,
                        GeomState *refState)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS, WHERE
  // ... ONLY THE AXIAL STRAIN (EXX) AND AXIAL STRESS (SXX) ARE CALCULATED

  int avgnum  = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  OutputInfo::FrameType oframe = oinfo[fileNumber].oframe;

  int k;

  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;

  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;

  // ... WRITE CURRENT TIME VALUE
  if(avgnum == 0 || avgnum == -1) {
    if(oinfo[fileNumber].nodeNumber == -1)
      fprintf(oinfo[fileNumber].filptr,"  % *.*E\n",w,p,time);
    else
      fprintf(oinfo[fileNumber].filptr," % *.*E",w,p,time);
  }

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(stress == 0)
    stress = new Vector(numnodes,0.0);
  else 
    stress->zero();

  if(weight == 0)
    weight = new Vector(numnodes,0.0);
  else 
    weight->zero();

  if(elDisp == 0)
    elDisp = new Vector(maxNumDOFs,0.0);

  Vector *sol = 0;
  if(sinfo.getNLInfo().linearelastic) {
    sol = new Vector(c_dsa->size());
    geomState.get_tot_displacement(*sol);
  }

  int iele;
  if((elstress == 0) || (elweight == 0) || (p_elstress == 0 && oframe != OutputInfo::Global)) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
      if(avgnum == -1) {
        maxNodesPerElement = std::max(maxNodesPerElement, allCorot[iele]->getNumGaussPoints());
      }
    }
    if(elstress == 0) elstress = new Vector(maxNodesPerElement, 0.0);
    if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
    if(p_elstress == 0 && oframe != OutputInfo::Global) p_elstress = new FullM(maxNodesPerElement,9);
  }

  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  // Either get the nodal temperatures from the input file or
  // from the thermal model
  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >= 0) nodalTemperatures = temprcvd;

  int flag;
  for(iele = 0; iele < numElements(); ++iele) {
    if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isMpcElement()) continue;
    if ((avgnum == 0 || avgnum == -1) && (oinfo[fileNumber].nodeNumber != -1) &&
        (oinfo[fileNumber].nodeNumber != packedEset[iele]->getGlNum())) continue;

    int NodesPerElement = packedEset[iele]->numNodes();

    // Don't include elements with one or more nodes not in the group if nodalpartialgroup (avgnum = 3) is requested
    if (avgnum == 3) {
      int groupId = oinfo[fileNumber].groupNumber;
      if (groupId > 0) {
        packedEset[iele]->nodes(nodeNumbers);
        std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
        std::set<int>::iterator it;
        for (int iNode = 0; iNode < NodesPerElement; ++iNode)
          if((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
        if(it == groupNodes.end()) continue;
      }
    }

    elDisp->zero();
    elstress->zero();
    elweight->zero();

// extract deformations from current Geometry State of structure
    if(sinfo.getNLInfo().linearelastic) {
      // DETERMINE ELEMENT DISPLACEMENT VECTOR
      for (k = 0; k < allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if (cn >= 0)
          (*elDisp)[k] = (*sol)[cn];
        else
          (*elDisp)[k] = 0; // XXX bcx[(*allDOFs)[iele][k]];
      } 
      flag = 1;
    } else
    allCorot[iele]->extractDeformations(geomState, nodes,
                                        elDisp->data(), flag);
    int DofsPerElement = packedEset[iele]->numDofs();

    // get element's nodal temperatures
    if(flag == 1) {
      elemNodeTemps.zero();
      if(nodalTemperatures) packedEset[iele]->nodes(nodeNumbers);
      for(int iNode = 0; iNode < NodesPerElement; ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
    }

    // ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
    if(oframe != OutputInfo::Global && ((stressIndex >= 0 && stressIndex <= 5) || (stressIndex >= 7 && stressIndex <= 12))
       && (flag == 1 || flag == 2)) { // transform non-invariant stresses/strains from basic frame to DOF_FRM or CFRAME

      // First, calculate stress/strain tensor for each node of the element
      p_elstress->zero();
      int strInd = (stressIndex >= 0 && stressIndex <= 5) ? 0 : 1;
      if (flag == 1) {
        // USE LINEAR STRESS ROUTINE
        packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                       *elDisp, strInd, surface,
                                       elemNodeTemps.data());
      }
      else {
        // USE NON-LINEAR STRESS ROUTINE
        // note: in this case the element nodal temperatures are extracted from geomState inside the function
        allCorot[iele]->getNLAllStress(*p_elstress, *elweight, geomState,
                                       refState, nodes, strInd, surface);
      }

      // Second, transform stress/strain tensor to nodal or material frame coordinates
      transformStressStrain(*p_elstress, iele, oframe);

      // Third, extract the requested stress/strain value from the stress/strain tensor
      for(int iNode = 0; iNode < NodesPerElement; ++iNode) {
        if(strInd == 0)
          (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex];
        else
          (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex-7];
      }
    }
    else {
      int index = (stressIndex == 31) ? 6 : stressIndex;
      if (flag == 1) {
        // USE LINEAR STRESS ROUTINE
        packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                      *elDisp, index, surface,
                                      elemNodeTemps.data(), ylayer,
                                      zlayer, avgnum);

      }
      else if (flag == 2) {
        // USE NON-LINEAR STRESS ROUTINE
        // note: in this case the element nodal temperatures are extracted from geomState inside the function
        allCorot[iele]->getNLVonMises(*elstress, *elweight, geomState,
                                      refState, nodes, index, surface,
                                      ylayer, zlayer, avgnum);

      } else {
        // NO STRESS RECOVERY
      }
    }

    // ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
    if(avgnum == 0 || avgnum == -1) {
      if(stressIndex != 31) {
        int numPoints = (avgnum == -1) ? allCorot[iele]->getNumGaussPoints() : NodesPerElement;
        for(k=0; k<numPoints; ++k)
          fprintf(oinfo[fileNumber].filptr," % *.*E",w,p,(*elstress)[k]);
        fprintf(oinfo[fileNumber].filptr,"\n");
      }
    }
    // ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
    else {
      for(k=0; k<NodesPerElement; ++k) {
        (*stress)[(*elemToNode)[iele][k]] += (*elstress)[k];
        (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
      }
    }
  }

// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

  if(avgnum > 0) {

    if(oinfo[fileNumber].nodeNumber == -1) {
      int numNodes = geoSource->numNode();
      int numNodesOut = (outFlag) ? exactNumNodes : numNodes;
      double *data = new double[numNodesOut];
      for(k=0; k<numNodesOut; ++k) data[k] = 0;
      for(k=0; k<numNodes; ++k) {
        int l = (outFlag) ? nodeTable[k]-1 : k;
        if(l < 0) continue;
        if(k < numnodes && (*weight)[k] != 0)
          data[l] = (*stress)[k]/=(*weight)[k];
      }
      if(stressIndex != 31) 
        geoSource->outputNodeScalars(fileNumber, data, numNodesOut, time);
      else {
        double stressmax(0);
        for(k = 0; k < numNodes; ++k) {
          if((*stress)[k]>stressmax) stressmax = (*stress)[k];
        }
        *aggregatedStress = 0.0;
        for(k = 0; k < numNodes; ++k) *aggregatedStress += exp(sinfo.ksParameter*((*stress)[k]-stressmax));
        *aggregatedStress = log(*aggregatedStress);
        *aggregatedStress /= sinfo.ksParameter;
        *aggregatedStress += stressmax;
        geoSource->outputEnergy(fileNumber, 0, *aggregatedStress);
      }
      delete [] data;
    }
    else {
      
      if((*weight)[oinfo[fileNumber].nodeNumber] == 0.0)
        fprintf(oinfo[fileNumber].filptr," %*.*E % *.*E\n",w,p,time,w,p,0.0);
      else
        fprintf(oinfo[fileNumber].filptr," %*.*E % *.*E\n",w,p,time,w,p,(*stress)[oinfo[fileNumber].nodeNumber]/=(*weight)[oinfo[fileNumber].nodeNumber]);
    }
  }
  fflush(oinfo[fileNumber].filptr);

  delete [] nodeNumbers;
  if(sol) delete sol;
}

void
Domain::getPrincipalStress(GeomState &geomState, Corotator **allCorot,
                           int fileNumber, int stressIndex, double time,
                           GeomState *refState)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // set stress VS. strain for element subroutines
  int strInd;
  int strDir;
  if ((stressIndex==0)||(stressIndex==1)||(stressIndex==2)) {
    strInd = 0;
    strDir = stressIndex+1;
  }
  else if ((stressIndex==3)||(stressIndex==4)||(stressIndex==5)) {
    strInd = 1;
    strDir = stressIndex-2;
  }
  else {
    fprintf(stderr,"*** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS

  int avgnum  = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  int j,k;

  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;
  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;

  // ... OUTPUT NODE NUMBER
  int n = oinfo[fileNumber].nodeNumber;

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(p_stress == 0) p_stress = new FullM(numnodes,6);
  if(weight == 0) weight = new Vector(numnodes,0.0);
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  int iele;
  if((p_elstress == 0) || (elweight == 0)) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    if(p_elstress == 0) p_elstress = new FullM(maxNodesPerElement,9);
    if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
  }

  // zero the vectors
  p_stress->zero();
  weight->zero();

  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  // Either get the nodal temperatures from the input file or
  // from the thermal model
  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >= 0) nodalTemperatures = temprcvd;

  // ... WRITE CURRENT TIME VALUE
  if(avgnum == 0) {
    fprintf(oinfo[fileNumber].filptr,"  % *.*E\n",w,p,time);
  }

  int flag;

  for(iele = 0; iele < numElements(); ++iele) {
    if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isMpcElement()) continue;

    int NodesPerElement = elemToNode->num(iele);

    // Don't include elements with one or more nodes not in the group if nodalpartialgroup (avgnum = 3) is requested
    if (avgnum == 3) {
      int groupId = oinfo[fileNumber].groupNumber;
      if (groupId > 0) {
        packedEset[iele]->nodes(nodeNumbers);
        std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
        std::set<int>::iterator it;
        for (int iNode = 0; iNode < NodesPerElement; ++iNode)
          if((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
        if(it == groupNodes.end()) continue;
      }
    }

    elDisp->zero();
    p_elstress->zero();
    elweight->zero();

    // extract deformations from current Geometry State of structure
    allCorot[iele]->extractDeformations(geomState, nodes,
                                        elDisp->data(), flag);

    if (flag == 1) {
      // get element's nodal temperatures
      elemNodeTemps.zero();
      if(nodalTemperatures) packedEset[iele]->nodes(nodeNumbers);
      for(int iNode = 0; iNode < packedEset[iele]->numNodes(); ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }

      // USE LINEAR STRESS ROUTINE
      // ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
      packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                     *elDisp, strInd, surface,
                                     elemNodeTemps.data());
    }
    else if (flag == 2) {
      // USE NON-LINEAR STRESS ROUTINE
      // note: in this case the element nodal temperatures are extracted from geomState inside the function
      allCorot[iele]->getNLAllStress(*p_elstress, *elweight, geomState,
                                     refState, nodes, strInd, surface);
    }
    else {
      // NO STRESS RECOVERY
    }

    // ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
    for(k=0; k<NodesPerElement; ++k) {
      for(j=0; j<6; ++j) {
        (*p_stress)[(*elemToNode)[iele][k]][j] += (*p_elstress)[k][j];
      }
      (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
    }

    // ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
    //     THIS WRITES THE CHOSEN PRINCIPAL STRESS FOR EACH ELEMENT
    if(avgnum == 0) {
      for(k=0; k<NodesPerElement; ++k)
         fprintf(oinfo[fileNumber].filptr," % *.*E",w,p,(*p_elstress)[k][5+strDir]);
       fprintf(oinfo[fileNumber].filptr,"\n");
    }
  }

  // ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
  // ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

  if(avgnum > 0) {

    if(n == -1) {
      for(k=0; k<numnodes; ++k) {
        if((*weight)[k] == 0.0) {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] = 0.0;
          }
        } else  {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] /= (*weight)[k];
          }
        }
      }
    } else {
      if((*weight)[n] == 0.0) {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] = 0.0;
        }
      } else {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] /= (*weight)[n];
        }
      }
    }

    // ... CALCULATE PRINCIPALS AT EACH NODE

    double svec[6], pvec[3];
    if(n == -1) {
      int numNodes = geoSource->numNode();
      int numNodesOut = (outFlag) ? exactNumNodes : numNodes;
      double *globalPVec = new double[numnodes];
      for(k=0; k<numNodes; ++k) {
        int l = (outFlag) ? nodeTable[k]-1 : k;
        if(l < 0) continue;
        for (j=0; j<6; ++j) {
          svec[j] = (*p_stress)[k][j];
        }
        // Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        globalPVec[l] = pvec[strDir-1];
      }
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodesOut, time);
    }
    else {
      for (j=0; j<6; ++j) {
        svec[j] = (*p_stress)[n][j];
      }
      // Convert Engineering to Tensor Strains
      if(strInd != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      geoSource->outputNodeScalars(fileNumber, pvec+strDir-1, 1);
    }
  }
  else {
    fflush(oinfo[fileNumber].filptr);
  }

  delete [] nodeNumbers;
}

// Nonlinear version of getElementForces
void
Domain::getElementForces(GeomState &geomState, Corotator **allCorot,
                         int fileNumber, int forceIndex, double time)
{
  // allocate integer array to store node numbers
  // make sure that maxNumNodes is computed
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();
  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag)  nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  // ... FORCES ARE CALCULATED FOR BARS AND BEAMS CURRENTLY

  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  // ... Watch out for this as the max number of nodes may increase
  int NodesPerElement = 2;

  FullM forces(numele, NodesPerElement);

  int i, iele;

  // note, we are reusing elstress here for an elemental force vector
  if(elstress == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    elstress = new Vector(maxNodesPerElement, 0.0);
  }

  int flag;

  for(iele=0; iele<numele; ++iele) {

    packedEset[iele]->nodes(nodeNumbers);

// ... DETERMINE ELEMENT DISPLACEMENT VECTOR

    allCorot[iele]->extractDeformations(geomState, nodes, elDisp->data(), flag);

// ... CALCULATE INTERNAL FORCE VALUE FOR EACH ELEMENT

    elemNodeTemps.zero();
    if(packedEset[iele]->getProperty()) {
      for(int iNode = 0; iNode < packedEset[iele]->numNodes(); ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
    }

    packedEset[iele]->getIntrnForce(*elstress,nodes,elDisp->data(),forceIndex,elemNodeTemps.data());

// ... COPY ELEMENT'S NODAL FORCES INTO A TOTAL FORCE MATRIX

    for(i=0; i<NodesPerElement; ++i)
      forces[iele][i] = (*elstress)[i];

  }

// ... PRINT THE ELEMENT FORCES TO A FILE
   geoSource->outputElemVectors(fileNumber, forces.data(), numele, time);

   delete [] nodeNumbers;
}

// Nonlinear restart file
void
Domain::writeRestartFile(double time, int timeIndex, Vector &v_n, Vector &a_n,
                         GeomState *geomState, const char *ext)
{
// either test for pointer or frequency > 0

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if((timeIndex % sinfo.nRestart == 0) || (time >= sinfo.tmax-0.1*domain->solInfo().getTimeStep()) || domain->solInfo().stop_AeroS) {

   int fn;
   if(strlen(ext) != 0) {
     char *currentRestartFile = new char[strlen(cinfo->currentRestartFile)+strlen(ext)+1];
     strcpy(currentRestartFile, cinfo->currentRestartFile);
     strcat(currentRestartFile, ext);
     fn = open(currentRestartFile, O_WRONLY | O_CREAT, 0666);
     delete [] currentRestartFile;
   } else
   fn = open(cinfo->currentRestartFile, O_WRONLY | O_CREAT, 0666);
   if(fn >= 0) {
     int numElemStates = geomState->getTotalNumElemStates();
     double *buffer = new double[std::max(numElemStates, 12*numnodes)];
     geomState->setVelocityAndAcceleration(v_n, a_n);

     int writeSize;
     writeSize = write(fn, &timeIndex, sizeof(int));
     if(writeSize != sizeof(int))
       fprintf(stderr," *** ERROR: Writing restart file time index\n");

     writeSize = write(fn, &time, sizeof(double));
     if(writeSize != sizeof(double))
       fprintf(stderr," *** ERROR: Writing restart file time\n");

     geomState->getVelocities(buffer);
     writeSize = write(fn, buffer, numnodes*6*sizeof(double));
     if(int(writeSize) != int(numnodes*6*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file velocity\n");

     geomState->getPositions(buffer);
     writeSize = write(fn, buffer, numnodes*3*sizeof(double));
     if(int(writeSize) != int(numnodes*3*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file position\n");

     geomState->getRotations(buffer);
     writeSize = write(fn, buffer, numnodes*12*sizeof(double));
     if(int(writeSize) != int(numnodes*12*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file rotation\n");

     // new method of storing the element states in the GeomState object
     geomState->getElemStates(buffer);
     writeSize = write(fn, buffer, numElemStates*sizeof(double));
     if(int(writeSize) != int(numElemStates*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file element state\n");

     // PJSA 9-17-2010 (note: this idea of the element storing the internal states is deprecated
     // and will eventually be removed
     int numEle = packedEset.last();
     for(int i = 0; i < numEle; ++i)
       packedEset[i]->writeHistory(fn);

     // write accelerations
     geomState->getAccelerations(buffer);
     writeSize = write(fn, buffer, numnodes*6*sizeof(double));
     if(int(writeSize) != int(numnodes*6*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file acceleration\n");

     delete [] buffer;

     close(fn);
   } else {
      perror(" *** ERROR: Restart file could not be opened: ");
      exit(-1);
   }
 }
}

void
Domain::readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                        Vector &v_p, double *bcx, double *vcx,
                        GeomState &geomState, const char *ext)
{
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(cinfo->lastRestartFile) {
   int fn;
   if(strlen(ext) != 0) {
     char *lastRestartFile = new char[strlen(cinfo->lastRestartFile)+strlen(ext)+1];
     strcpy(lastRestartFile, cinfo->lastRestartFile);
     strcat(lastRestartFile, ext);
     fn = open(lastRestartFile, O_RDONLY);
     delete [] lastRestartFile;
   } else
   fn = open(cinfo->lastRestartFile, O_RDONLY);
   if(fn >= 0) {
     int numElemStates = geomState.getTotalNumElemStates();
     double *buffer = new double[std::max(numElemStates, 12*numnodes)];

     int restartTIndex;
     double restartT;
     int readSize = read(fn, &restartTIndex, sizeof(int));
     if(readSize != sizeof(int))
       fprintf(stderr," *** ERROR: Reading restart file time index\n");
     if(strcmp(cinfo->FlagRST,"new") == 0)
       sinfo.initialTimeIndex = 0;
     else
       sinfo.initialTimeIndex = restartTIndex;

     readSize = read(fn, &restartT, sizeof(double));
     if(readSize != sizeof(double))
       fprintf(stderr," *** ERROR: Reading restart file time\n");
     if(strcmp(cinfo->FlagRST,"new") == 0)
       sinfo.initialTime = 0;
     else
       sinfo.initialTime = restartT;

     readSize = read(fn, buffer, numnodes*6*sizeof(double));
     if(int(readSize) != int(numnodes*6*sizeof(double)))
       fprintf(stderr," *** ERROR: Reading restart file velocity\n");
     geomState.setVelocities(buffer);

     readSize = read(fn, buffer, numnodes*3*sizeof(double));
     if(int(readSize) != int(numnodes*3*sizeof(double)))
       fprintf(stderr," *** ERROR: Reading restart file position\n");
     geomState.setPositions(buffer);

     readSize = read(fn, buffer, numnodes*12*sizeof(double));
     if(int(readSize) != int(numnodes*12*sizeof(double)))
       fprintf(stderr," *** ERROR: Reading restart file rotation\n");
     geomState.setRotations(buffer);

     // new method of storing the element states in the GeomState object
     readSize = read(fn, buffer, numElemStates*sizeof(double));
     if(int(readSize) != int(numElemStates*sizeof(double)))
       fprintf(stderr," *** ERROR: Reading restart file element state\n");
     geomState.setElemStates(buffer);

     // old method of storing the element states in the Element objects (e.g. bt shell)
     int numEle = packedEset.last();
     for(int i = 0; i < numEle; ++i) 
       packedEset[i]->readHistory(fn);

     readSize = read(fn, buffer, numnodes*6*sizeof(double));
     if(int(readSize) != int(numnodes*6*sizeof(double)))
       fprintf(stderr," *** ERROR: Reading restart file acceleration\n");
     geomState.setAccelerations(buffer);

     delete [] buffer;
     close(fn);

     d_n.zero();
     v_n.zero();
     a_n.zero();
     v_p.zero();
     for(int i = 0; i < numNodes(); ++i) {

       int xloc  = c_dsa->locate(i, DofSet::Xdisp);
       int xloc1 =   dsa->locate(i, DofSet::Xdisp);

       if(xloc >= 0) {
         d_n[xloc] = geomState[i].x - nodes[i]->x;
         v_n[xloc] = geomState[i].v[0];
         a_n[xloc] = geomState[i].a[0];
       }
       else if(xloc1 >= 0) {
         bcx[xloc1] = geomState[i].x - nodes[i]->x;
         vcx[xloc1] = geomState[i].v[0];
       }

       int yloc  = c_dsa->locate(i, DofSet::Ydisp);
       int yloc1 =   dsa->locate(i, DofSet::Ydisp);

       if(yloc >= 0) {
         d_n[yloc] = geomState[i].y - nodes[i]->y;
         v_n[yloc] = geomState[i].v[1];
         a_n[yloc] = geomState[i].a[1];
       }
       else if(yloc1 >= 0) {
         bcx[yloc1] = geomState[i].y - nodes[i]->y;
         vcx[yloc1] = geomState[i].v[1];
       }

       int zloc  = c_dsa->locate(i, DofSet::Zdisp);
       int zloc1 =   dsa->locate(i, DofSet::Zdisp);

       if(zloc >= 0) {
         d_n[zloc] = geomState[i].z - nodes[i]->z;
         v_n[zloc] = geomState[i].v[2];
         a_n[zloc] = geomState[i].a[2];
       }
       else if(zloc1 >= 0) {
         bcx[zloc1] = geomState[i].z - nodes[i]->z;
         vcx[zloc1] = geomState[i].v[2];
       }

       int xrot  = c_dsa->locate(i, DofSet::Xrot);
       int xrot1 =   dsa->locate(i, DofSet::Xrot);

       if(xrot >= 0) {
         d_n[xrot] = geomState[i].theta[0];
         v_n[xrot] = geomState[i].v[3];
         a_n[xrot] = geomState[i].a[3];
       }
       else if(xrot1 >= 0) {
         bcx[xrot1] = geomState[i].theta[0];
         vcx[xrot1] = geomState[i].v[3];
       }

       int yrot  = c_dsa->locate(i, DofSet::Yrot);
       int yrot1 =   dsa->locate(i, DofSet::Yrot);

       if(yrot >= 0) {
         d_n[yrot] = geomState[i].theta[1];
         v_n[yrot] = geomState[i].v[4];
         a_n[yrot] = geomState[i].a[4];
       }
       else if(yrot1 >= 0) {
         bcx[yrot1] = geomState[i].theta[1];
         vcx[yrot1] = geomState[i].v[4];
       }

       int zrot  = c_dsa->locate(i, DofSet::Zrot);
       int zrot1 =   dsa->locate(i, DofSet::Zrot);

       if(zrot >= 0) {
         d_n[zrot] = geomState[i].theta[2];
         v_n[zrot] = geomState[i].v[5];
         a_n[zrot] = geomState[i].a[5];
       }
       else if(zrot1 >= 0) {
         bcx[zrot1] = geomState[i].theta[2];
         vcx[zrot1] = geomState[i].v[5];
       }
     }
   }
   else {
     perror(" *** ERROR: Restart file could not be opened: ");
   }
 }
}

// nonlinear statics
void
Domain::computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
                             FullSquareMatrix *kel, double lambda, GeomState *refState)
{
  // TODO: include non-follower external forces on constrained dofs
  Vector elementInternalForce(maxNumDOF(), 0.0);
  Vector residual(numUncon(), 0.0);
  fc.zero();
  getInternalForce(*geomState, elementInternalForce, corotators, kel, residual, lambda, 0.0, refState, &fc);
}

// nonlinear dynamics
void
Domain::computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
                             FullSquareMatrix *kel, double time, GeomState *refState,
                             Vector &Vu, Vector &Au, double *vcx, double *acx,
                             SparseMatrix *_cuc, SparseMatrix *_ccc,
                             SparseMatrix *_muc, SparseMatrix *_mcc)
{
  // TODO: include non-follower external forces on constrained dofs
  Vector elementInternalForce(maxNumDOF(), 0.0);
  Vector residual(numUncon(), 0.0);
  fc.zero();

  getInternalForce(*geomState, elementInternalForce, corotators, kel, residual, 1.0, time, refState, &fc);

  CuCSparse *cuc = dynamic_cast<CuCSparse *>(_cuc);
  if(cuc) cuc->transposeMultAddNew(Vu.data(), fc.data()); // fc += Cuc^T * Vu

  CuCSparse *muc = dynamic_cast<CuCSparse *>(_muc);
  if(muc) muc->transposeMultAddNew(Au.data(), fc.data()); // fc += Muc^T * Vu

  CuCSparse *ccc = dynamic_cast<CuCSparse *>(_ccc);
  CuCSparse *mcc = dynamic_cast<CuCSparse *>(_mcc);
  Vector Dc(numDirichlet, 0.0);
  Vector Vc(numDirichlet, 0.0);
  Vector Ac(numDirichlet, 0.0);

  for(int i=0; i<numDirichlet; ++i) {
    int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
    if(dof < 0) continue;
    int cdof = c_dsa->invRCN(dof);
    if(cdof >= 0) {
      Vc[cdof] = vcx[dof];
      Ac[cdof] = (acx) ? acx[dof] : 0;
    }
  }

  if(ccc) ccc->multAddNew(Vc.data(), fc.data()); // fc += Ccc * Vc
  if(mcc) mcc->multAddNew(Ac.data(), fc.data()); // fc += Mcc * Ac
}

void
Domain::transformElemStiffAndForce(const GeomState &geomState, double *elementForce,
                                   FullSquareMatrix &kel, int iele, bool compute_tangents)
{
#ifdef USE_EIGEN3
  // if the element is an LMPC and the node has a DOF_FRM then the element force and or stiffness are be defined in the DOF_FRM
  bool basicDofCoords = true;
  if(packedEset[iele]->isMpcElement() && !domain->solInfo().basicDofCoords) {
    LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]);
    if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
       lmpcons->getSource() == mpc::TiedSurfaces)) basicDofCoords = false;
  }

  // Convert from eulerian spatial to total lagrangian or updated lagrangian spatial
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> >
    G(elementForce, packedEset[iele]->numDofs());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >
    H(kel.data(),packedEset[iele]->numDofs(),packedEset[iele]->numDofs());
  int numDofs = packedEset[iele]->numDofs();
  int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
  int dofsPerNode = (numDofs-packedEset[iele]->getNumMPCs())/numNodes;
  if((dofsPerNode == 6 || dofsPerNode == 3) && (numDofs-packedEset[iele]->getNumMPCs())%numNodes == 0) {
    int *nodes = packedEset[iele]->nodes();
    for(int k = 0; k < numNodes; ++k) {
      Eigen::Vector3d Psi;
      Psi << geomState[nodes[k]].theta[0], geomState[nodes[k]].theta[1], geomState[nodes[k]].theta[2];
      if(!basicDofCoords) { if(NFrameData *cd = getNodes().dofFrame(nodes[k])) cd->transformVector3(&Psi[0]); } // Transform from basic frame to DOF_FRM

      Eigen::Matrix3d T;
      tangential_transf(Psi, T);

      int p = dofsPerNode*k + (dofsPerNode-3);
      Eigen::Vector3d V = G.segment<3>(p);
      G.segment<3>(p) = (domain->solInfo().newmarkBeta == 0) ? (Jn[nodes[k]]*(T.transpose()*Jn[nodes[k]]*T).inverse()*T*V).eval() : (T*V).eval();
 
      if(compute_tangents) {
        Eigen::Matrix3d C1;
        directional_deriv1(Psi, V, C1);

        for(int l=0; l<dofsPerNode/3*numNodes; ++l) {
          int q = 3*l;
          H.block<3,3>(p,q) = T*H.block<3,3>(p,q);
          H.block<3,3>(q,p) = H.block<3,3>(q,p)*T.transpose();
        }
        for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
          int q = dofsPerNode*numNodes+l;
          H.block<3,1>(p,q) = T*H.block<3,1>(p,q);
          H.block<1,3>(q,p) = H.block<1,3>(q,p)*T.transpose();
        }
        H.block<3,3>(p,p) += 0.5*(C1 + C1.transpose());
      }
    }
    delete [] nodes;
  }
  else {
    std::cerr << " *** WARNING: Domain::transformElemStiffAndForce is not implemented for element " << packedEset[iele]->getGlNum()+1
              << " type " << packedEset[iele]->getElementType() << std::endl;
  }
#else
  std::cerr << "USE_EIGEN3 is not defined here in Domain::transformElemStiffAndForce\n";
  exit(-1);
#endif
}

void
Domain::transformElemStiffAndForce_S2E(const GeomState &geomState, double *elementForce,
                                       FullSquareMatrix &kel, int iele, bool compute_tangents)
{
#ifdef USE_EIGEN3
  // if the element is an LMPC and the node has a DOF_FRM then the element force and or stiffness are be defined in the DOF_FRM
  bool basicDofCoords = true;
  if(packedEset[iele]->isMpcElement() && !domain->solInfo().basicDofCoords) {
    LMPCons *lmpcons = dynamic_cast<LMPCons*>(packedEset[iele]);
    if(lmpcons && (lmpcons->getSource() == mpc::Lmpc || lmpcons->getSource() == mpc::NodalContact ||
       lmpcons->getSource() == mpc::TiedSurfaces)) basicDofCoords = false;
  }

  // Convert from eulerian spatial to constrained subset of SO(3), S²E
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> >
    G(elementForce, packedEset[iele]->numDofs());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >
    H(kel.data(),packedEset[iele]->numDofs(),packedEset[iele]->numDofs());
  int numDofs = packedEset[iele]->numDofs();
  int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
  int dofsPerNode = (numDofs-packedEset[iele]->getNumMPCs())/numNodes;
  if((dofsPerNode == 6 || dofsPerNode == 3) && (numDofs-packedEset[iele]->getNumMPCs())%numNodes == 0) {
    int *nodes = packedEset[iele]->nodes();
    for(std::set<int>::iterator it = elemAdj[iele].crot.begin(); it != elemAdj[iele].crot.end(); ++it) {
      int k =0; for(int i = 0; i < numNodes; ++i) if(nodes[i] == dbc[*it].nnum) { k = i; break; }
      Eigen::Vector3d Psi;
      Psi << geomState[nodes[k]].theta[0], geomState[nodes[k]].theta[1], geomState[nodes[k]].theta[2];
      if(!basicDofCoords) { if(NFrameData *cd = getNodes().dofFrame(nodes[k])) cd->transformVector3(&Psi[0]); } // Transform from basic frame to DOF_FRM

      Eigen::Matrix3d T;
      tangential_transf_S2E(Psi, ((basicDofCoords) ? getNodes().dofFrame(nodes[k]) : 0), dbc[*it].dofnum-3, T);

      int p = dofsPerNode*k + (dofsPerNode-3);
      Eigen::Vector3d V = G.segment<3>(p);
      G.segment<3>(p) = (T*V).eval();
 
      if(compute_tangents) {
        for(int l=0; l<dofsPerNode/3*numNodes; ++l) {
          int q = 3*l;
          H.block<3,3>(p,q) = T*H.block<3,3>(p,q);
          H.block<3,3>(q,p) = H.block<3,3>(q,p)*T.transpose();
        }
        for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
          int q = dofsPerNode*numNodes+l;
          H.block<3,1>(p,q) = T*H.block<3,1>(p,q);
          H.block<1,3>(q,p) = H.block<1,3>(q,p)*T.transpose();
        }
      }
    }
    delete [] nodes;
  }
  else {
    std::cerr << " *** WARNING: Domain::transformElemStiffAndForce_S2E is not implemented for element " << packedEset[iele]->getGlNum()+1
              << " type " << packedEset[iele]->getElementType() << std::endl;
  }
#endif
}

void
Domain::transformNodalMoment(const GeomState &geomState, double _G[],
                             double _H[][3], int inode, bool compute_tangents)
{
#ifdef USE_EIGEN3
  // transform from eulerian spatial description to total lagrangian or updated lagrangian spatial description
  Eigen::Map<Eigen::Matrix<double,3,1> > G(&_G[0]);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > H(&_H[0][0]);

  Eigen::Vector3d Psi;
  Psi << geomState[inode].theta[0], geomState[inode].theta[1], geomState[inode].theta[2];

  Eigen::Matrix3d T;
  tangential_transf(Psi, T);

  Eigen::Vector3d V = G;
  G = (domain->solInfo().newmarkBeta == 0) ? (Jn[inode]*(T.transpose()*Jn[inode]*T).inverse()*T*V).eval() : (T*V).eval();
  if(compute_tangents) {
    Eigen::Matrix3d C1;
    directional_deriv1(Psi, V, C1);

    H = (T*H*T.transpose()).eval() + 0.5*(C1 + C1.transpose());
  }
#else
  std::cerr << "USE_EIGEN3 is not defined here in Domain::transformNodalMoment\n";
  exit(-1);
#endif
}

void
Domain::transformElemStiff(const GeomState &geomState, FullSquareMatrix &kel, int iele)
{
#ifdef USE_EIGEN3
  // Convert from eulerian spatial to eulerian convected
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >
    H(kel.data(),packedEset[iele]->numDofs(),packedEset[iele]->numDofs());
  int numDofs = packedEset[iele]->numDofs();
  int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
  int dofsPerNode = (numDofs-packedEset[iele]->getNumMPCs())/numNodes;
  if((dofsPerNode == 3 || dofsPerNode == 6) && (numDofs-packedEset[iele]->getNumMPCs())%numNodes == 0) {
    int *nodes = packedEset[iele]->nodes();
    for(int k = 0; k < numNodes; ++k) {
      Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&geomState[nodes[k]].R[0][0],3,3);
      if(R.isIdentity()) continue;

      if(dofsPerNode == 6) {
        for(int l=0; l<2*numNodes; ++l) {
          H.block<3,3>(6*k+3,3*l) = (R.transpose()*H.block<3,3>(6*k+3,3*l)).eval();
          H.block<3,3>(3*l,6*k+3) = (H.block<3,3>(3*l,6*k+3)*R).eval();
        }
        for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
          H.block<3,1>(6*k+3,6*numNodes+l) = (R.transpose()*H.block<3,1>(6*k+3,6*numNodes+l)).eval();
          H.block<1,3>(6*numNodes+l,6*k+3) = (H.block<1,3>(6*numNodes+l,6*k+3)*R).eval();
        }
      }
      else if(dofsPerNode == 3) {
        for(int l=0; l<numNodes; ++l) {
          H.block<3,3>(3*k,3*l) = (R.transpose()*H.block<3,3>(3*k,3*l)).eval();
          H.block<3,3>(3*l,3*k) = (H.block<3,3>(3*l,3*k)*R).eval();
        }
        for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
          H.block<3,1>(3*k,3*numNodes+l) = (R.transpose()*H.block<3,1>(3*k,3*numNodes+l)).eval();
          H.block<1,3>(3*numNodes+l,3*k) = (H.block<1,3>(3*numNodes+l,3*k)*R).eval();
        }
      }
    }
    delete [] nodes;
  }
  else { 
    // XXX there are a few elements with 3 dofs at one node and 6 at the other, eg. element 201 and 117
    int *nodes = packedEset[iele]->nodes();
    for(int k = 0; k < numNodes; ++k) {
      Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&geomState[nodes[k]].R[0][0],3,3);
      if(!R.isIdentity()) {
        std::cerr << " *** WARNING: Domain::transformElemStiff is not implemented for element " << packedEset[iele]->getGlNum()+1
                  << " type " << packedEset[iele]->getElementType() << std::endl;
        break;
      }
    }
    delete [] nodes;
  }
#else
  std::cerr << "USE_EIGEN3 is not defined here in Domain::transformElemStiff\n";
  exit(-1);
#endif
}

void
Domain::getElementDisp(int iele, GeomState& geomState, Vector& disp)
{
  int *nn = packedEset[iele]->nodes();
  int dofs[DofSet::max_known_nonL_dof];
  double psi[3];

  for(int i=0,l=0; i<packedEset[iele]->numNodes(); ++i) {
    int ndofs = dsa->number(nn[i], DofSet::nonL_dof, dofs);
    for(int j=0; j<ndofs; ++j) {
      if(dofs[j] > -1) {
        for(int k=0; k<packedEset[iele]->numDofs(); ++k) {
          if(dofs[j] == (*allDOFs)[iele][k]) {
            switch(j) {
              case 0 : // x displacement
                disp[l++] = geomState[nn[i]].x - nodes[nn[i]]->x;
                break;
              case 1 : // y displacement 
                disp[l++] = geomState[nn[i]].y - nodes[nn[i]]->y;
                break;
              case 2 : // z displacement
                disp[l++] = geomState[nn[i]].z - nodes[nn[i]]->z;
                break;
              case 3 : case 4 : case 5 : // x,y,z rotations
                disp[l++] = geomState[nn[i]].theta[j-3];
                break;
              case 6 : case 7 : case 8 : // temperature and lagrange multipliers
                disp[l++] = geomState[nn[i]].x;
                break;
              default :
                disp[l++] = 0;
                break;
            }
            break;
          }
        }
      }
    }
  }
  delete [] nn;
}

void
Domain::getElementVelo(int iele, GeomState& geomState, Vector& velo)
{
  int *nn = packedEset[iele]->nodes();
  int dofs[DofSet::max_known_nonL_dof];

  // place holder in velo vector: l
  for(int i=0,l=0; i<packedEset[iele]->numNodes(); ++i) { // loop over nodes: i
    int ndofs = dsa->number(nn[i], DofSet::nonL_dof, dofs);
    for(int j=0; j<ndofs; ++j) {// loop over element dofs: j 
      if(dofs[j] > -1) { // check if dof is fixed
        for(int k=0; k<packedEset[iele]->numDofs(); ++k) { // loop over numdofs for that node: k 
          switch(j) {
              case 0 : // x velocity
                velo[l++] = geomState[nn[i]].v[j];
                break;
              case 1 : // y velocity
                velo[l++] = geomState[nn[i]].v[j];
                break;
              case 2 : // z velocity
                velo[l++] = geomState[nn[i]].v[j];
                break;
              case 3 : case 4 : case 5 : // x,y,z rotational velocity
                velo[l++] = geomState[nn[i]].v[j];
                break;
              case 6 : case 7 : case 8 : // temperature and lagrange multipliers
                velo[l++] = 0.0;
                break;
              default :
                velo[l++] = 0.0;
                break;
            }
            break;
        }
      }
    }  
  }
  delete [] nn;
}

void
Domain::getElementAccel(int iele, GeomState& geomState, Vector& accel)
{
  int *nn = packedEset[iele]->nodes();
  int dofs[DofSet::max_known_nonL_dof];

  // place holder in accel vector: l
  for(int i=0,l=0; i<packedEset[iele]->numNodes(); ++i) { // loop over nodes in element iele: i
    int ndofs = dsa->number(nn[i], DofSet::nonL_dof, dofs);
    for(int j=0; j<ndofs; ++j) {// for node i, loop over dofs: j 
      if(dofs[j] > -1) { // check if dof is fixed
        for(int k=0; k<packedEset[iele]->numDofs(); ++k) { // loop over numdofs for element iele: k 
          switch(j) {
              case 0 : // x acceleration
                accel[l++] = geomState[nn[i]].a[j];
                break;
              case 1 : // y acceleration 
                accel[l++] = geomState[nn[i]].a[j];
                break;
              case 2 : // z acceleration
                accel[l++] = geomState[nn[i]].a[j];
                break;
              case 3 : case 4 : case 5 : // x,y,z rotational accelerations
                accel[l++] = geomState[nn[i]].a[j];
                break;
              case 6 : case 7 : case 8 : // temperature and lagrange multipliers
                accel[l++] = 0.0;
                break;
              default :
                accel[l++] = 0.0;
                break;
            }
            break;         
        }
      } 
    }
  }
  delete [] nn;
}

void
Domain::computeEnergies(GeomState *geomState, Vector &force, double t, Vector *aeroForce, double *velocity,
                        Corotator **allCorot, SparseMatrix *M, SparseMatrix *C, double &Wela, double &Wkin,
                        double &Wdis, double &error)
{
  // Build displacement vector
  Vector disp(numUncon(), 0.0);
  geomState->get_tot_displacement(disp, false);

  double lambda, time;
  if(sinfo.isDynam()) { time = t; lambda = 1.0; }
  else {
    time = 0.0;
    double &dlambda = solInfo().getNLInfo().dlambda;
    double &maxLambda = solInfo().getNLInfo().maxLambda;
    lambda = (dlambda == maxLambda) ? maxLambda : t;
  }

  // Compute follower force
  Vector elemForce(maxNumDOFs);
  Vector folForce(numUncon(), 0.0);
  getFollowerForce(*geomState, elemForce, allCorot, (FullSquareMatrix *) NULL,
                   folForce, lambda, time, (Vector *) NULL, false);

  if(sinfo.isDynam()) {
    double pWext = Wext, pWdmp = Wdmp;
    StackVector vel(velocity, numUncon());
    computeExtAndDmpEnergies(disp, force, time, aeroForce, &vel, C, &folForce);

    Vector tmpVec(numUncon());
    if(M) {
      M->mult(vel, tmpVec);
      Wkin = 0.5 * (vel * tmpVec);
    }
    Wela = getStrainEnergy(geomState, allCorot);
    Wdis = getDissipatedEnergy(geomState, allCorot);

    // XXX consider sign of Wdmp+Wdis in this equation:
    error = (time == sinfo.initialTime) ? 0.0 : (Wela+Wkin+Wdmp+Wdis-Wext)-(pWela+pWkin+pWdmp+pWdis-pWext);

    pWela = Wela;
    pWkin = Wkin;
    pWdis = Wdis;
  }
  else { // nonlinear statics
    Wext  = (lambda*force + folForce) * disp;
    Waero = 0;
    Wdmp  = 0;
    Wkin  = 0;
    Wela  = getStrainEnergy(geomState, allCorot);
    Wdis  = getDissipatedEnergy(geomState, allCorot);
    error = 0;
  }
}

double
Domain::getStrainEnergy(GeomState *geomState, Corotator **allCorot)
{
  // Compute strain energy
  double W = 0;
  for(int iele = 0; iele < numele; ++iele) {
    if(allCorot[iele] != NULL) {
      W += allCorot[iele]->getElementEnergy(*geomState, nodes);
    }
  }
  return W;
}

double
Domain::getDissipatedEnergy(GeomState *geomState, Corotator **allCorot)
{
  // Compute dissipated energy due to plastic deformation
  double D = 0;
  for(int iele = 0; iele < numele; ++iele) {
    if(allCorot[iele] != NULL) {
      D += allCorot[iele]->getDissipatedEnergy(*geomState, nodes);
    }
  }
  return D;
}

double*
Domain::getDissipatedEnergy(GeomState *geomState, Corotator **allCorot, int groupIndex)
{
  // (AN) Note: Here the function assumes that all the attributes assigned to dissipat
  // output card have a corresponding non-linear material (matusage) card.
std::map<int, Group> &group = geoSource->group;

  int numAttributesInGrp = group[groupIndex].attributes.size();
  int numAttributes = geoSource->getMaxAttributeNum();

  if(numAttributesInGrp > numAttributes) {
    fprintf(stderr, "***Error: Number of attributes provided in group %d is greater than max atributes available\n", groupIndex);
    exit(-1);
  }

  double* D = new double[numAttributesInGrp];
  double* temp = new double[numAttributes];
  for(int a=0; a<numAttributes; a++) temp[a] = 0.0;

  for(int iele=0; iele<numele; iele++) {
    if(allCorot[iele]) {
      int aele = packedEset[iele]->getElementAttribute();
      if(aele > numAttributes) {
        fprintf(stderr, "***Error: unknown attribute %d found for element %d.\n", aele, packedEset[iele]->getGlNum());
        exit(-1);
      }
      if(!group[groupIndex].isAttributeInGroup(aele)) continue;
      temp[aele] += allCorot[iele]->getDissipatedEnergy(*geomState, nodes);
    }
  }

  // for(int a=0; a<numAttributes; a++) fprintf(stdout, "temp[%d] = %e\n", a, temp[a]);
  
  for(int i=0; i<numAttributesInGrp; i++) D[i] = temp[group[groupIndex].attributes[i]];
  
  delete [] temp;
  
  return D;

}

void 
Domain::computeNLStaticWRTthicknessSensitivity(int sindex, 
                                               AllSensitivities<double> &allSens,
                                               GeomState *refState,
                                               GeomState *geomState,
                                               Corotator **allCorot)
{
#ifdef USE_EIGEN3
     allSens.linearstaticWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
     std::map<int, Group> &group = geoSource->group;
     std::map<int, AttributeToElement> &atoe = geoSource->atoe;
     for(int g=0; g<numThicknessGroups; ++g) {
       allSens.linearstaticWRTthick[g] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numUncon(),1);
       allSens.linearstaticWRTthick[g]->setZero();   
     }
     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
       int groupIndex = thicknessGroups[iparam];
       for(int aindex = 0; aindex < group[groupIndex].attributes.size(); ++aindex) {
         for(int eindex =0; eindex < atoe[group[groupIndex].attributes[aindex]].elems.size(); ++eindex) {
           int iele = atoe[group[groupIndex].attributes[aindex]].elems[eindex];
           if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
           int DofsPerElement = packedEset[iele]->numDofs();
           GenVector<double> dFintdThick(DofsPerElement,0.0);
           (*allCorot[iele]).getInternalForceThicknessSensitivity(refState, *geomState, nodes, dFintdThick, 0, 0);
           auto dofs = (*allDOFs)[iele];
           const auto &unconstrNum = c_dsa->getUnconstrNum();
           const auto &constrndNum = c_dsa->getConstrndNum();
           for(int k = 0; k < DofsPerElement; ++k) {
             int dofk = unconstrNum[dofs[k]];
             if(dofs[k] < 0 || dofk < 0) continue;  // Skip undefined/constrained dofs
             (*allSens.linearstaticWRTthick[iparam])(dofk,0) += dFintdThick[k];
           }
         }
       }
     }
     if(domain->gravityFlag()) subtractGravityForceSensitivityWRTthickness(sindex,allSens); // TODO-> must consider other
                                                                                            // external forces that depend
                                                                                            // on thickness
#endif
}


void 
Domain::computeNLStressVMWRTthicknessDirectSensitivity(int sindex, AllSensitivities<double> &allSens,
                                                       GeomState *geomState, Corotator **allCorot, bool isDynam)
{
#ifdef USE_EIGEN3
     // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO THICKNESS
     allSens.vonMisesWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), numThicknessGroups);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numNodes(), 1);
     allSens.vonMisesWRTthick->setZero();     stressWeight.setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     Vector *sol = 0;
     if(sinfo.getNLInfo().linearelastic) {
       sol = new Vector(c_dsa->size());
       geomState->get_tot_displacement(*sol);
     }
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     
     for(int iele = 0; iele< numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenVector<double> dStressdThick(NodesPerElement);
       GenVector<double> weight(NodesPerElement,0.0);
       int surface = senInfo[sindex].surface;
       elDisp->zero();
       int flag;
       if(sinfo.getNLInfo().linearelastic) {
         // DETERMINE ELEMENT DISPLACEMENT VECTOR
         for (int k = 0; k < allDOFs->num(iele); ++k) {
           int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
           if (cn >= 0)
             (*elDisp)[k] = (*sol)[cn];
           else
             (*elDisp)[k] = 0; // XXX bcx[(*allDOFs)[iele][k]];
         } 
         flag = 1;
       } else {
         allCorot[iele]->extractDeformations(*geomState, nodes, elDisp->data(), flag);       
         if(flag != 1) {
          filePrint(stderr,
                " ... Error: nonlinear sensitivity for non-corotational elements is not implemented yet. exiting!\n");
          exit(-1);
        }
       }
       // Determine element displacement vector
       if(thgreleFlag[iele]) { 
         packedEset[iele]->getVonMisesThicknessSensitivity(dStressdThick, weight, nodes, *elDisp, 6, surface); 
       } else {
         dStressdThick.zero();
         Vector elstress(NodesPerElement,0.0);
         packedEset[iele]->getVonMises(elstress, weight, nodes, *elDisp, -1, surface, NULL, 0, 0, 0);
       }
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           int iparam = thpaIndex[iele];
           if(iparam >= 0) (*allSens.vonMisesWRTthick)(node, iparam) += dStressdThick[k];
           stressWeight(node, 0) += weight[k];
         }
       }
     } 

     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
      // average vonMisesWRTthick vector
      for(int inode = 0; inode < numNodes(); ++inode)  {
        if(stressWeight(inode, 0) == 0.0)
          (*allSens.vonMisesWRTthick)(inode, iparam) = 0.0;
        else
          (*allSens.vonMisesWRTthick)(inode, iparam) /= stressWeight(inode, 0);
      }
      if(!isDynam) allSens.vonMisesWRTthick->col(iparam) += *allSens.vonMisesWRTdisp * (*allSens.dispWRTthick[iparam]);
     }
#endif
}

void 
Domain::computeAggregatedNLStressVMWRTthicknessSensitivity(int sindex,
                                                           AllSensitivities<double> &allSens,
                                                           GeomState *geomState,
                                                           GeomState *refState, 
                                                           Corotator **allCorot, Vector &stress, Vector &stressWeight,
                                                           Vector &sensCoef, 
                                                           bool isDynam)
{
#ifdef USE_EIGEN3
     int surface = senInfo[sindex].surface;
     allSens.aggregatedVonMisesWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numThicknessGroups);
     allSens.aggregatedVonMisesWRTthick->setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     Vector *sol = 0;
     if(sinfo.getNLInfo().linearelastic) {
       sol = new Vector(c_dsa->size());
       geomState->get_tot_displacement(*sol);
     }
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele = 0; iele < numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenVector<double> dStressdThick(NodesPerElement);
       GenVector<double> elweight(NodesPerElement,0.0);
       elDisp->zero();
       int flag;
       if(sinfo.getNLInfo().linearelastic) {
         // DETERMINE ELEMENT DISPLACEMENT VECTOR
         for (int k = 0; k < allDOFs->num(iele); ++k) {
           int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
           if (cn >= 0)
             (*elDisp)[k] = (*sol)[cn];
           else
             (*elDisp)[k] = 0; // XXX bcx[(*allDOFs)[iele][k]];
         } 
         flag = 1;
       } else {
         allCorot[iele]->extractDeformations(*geomState, nodes, elDisp->data(), flag);
         if(flag != 1) {
           filePrint(stderr,
               " ... Error: nonlinear sensitivity for non-corotational elements is not implemented yet. exiting!\n");
           exit(-1);
          }
       }
       if(thgreleFlag[iele]) { 
         packedEset[iele]->getVonMisesThicknessSensitivity(dStressdThick, elweight, nodes, *elDisp, 6, surface); 
       } else {
         dStressdThick.zero();
         Vector elstress(NodesPerElement,0.0);
         packedEset[iele]->getVonMises(elstress, elweight, nodes, *elDisp, -1, surface, NULL, 0, 0, 0);
       }
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           int iparam = thpaIndex[iele];
           if(iparam >= 0) (*allSens.aggregatedVonMisesWRTthick)(iparam) += (dStressdThick[k]*sensCoef[node]);
//stressvmWRTthick(node, iparam) += dStressdThick[k]; 
         }
       }
     } 

     if(!isDynam) {
       for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
         (*allSens.aggregatedVonMisesWRTthick)(iparam) -=
                allSens.lambdaAggregatedStressVM->adjoint()*(allSens.linearstaticWRTthick[iparam]->col(0));
       }
     }
#endif
}

void 
Domain::computeNLStressVMWRTthicknessAdjointSensitivity(int sindex,
                                                        AllSensitivities<double> &allSens,
                                                        GeomState *geomState, Corotator **allCorot, 
                                                        bool *includeStressNodes, bool isDynam)
{
#ifdef USE_EIGEN3
     // ... COMPUTE DERIVATIVE OF VON MISES STRESS WITH RESPECT TO THICKNESS
     allSens.vonMisesWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numStressNodes, numThicknessGroups);
     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numStressNodes, 1);
     allSens.vonMisesWRTthick->setZero();     stressWeight.setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     Vector *sol = 0;
     if(sinfo.getNLInfo().linearelastic) {
       sol = new Vector(c_dsa->size());
       geomState->get_tot_displacement(*sol);
     }
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     for(int iele=0; iele<numele; ++iele) {
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       if (!includeStressNodes[iele]) continue;
       int NodesPerElement = elemToNode->num(iele);
       GenVector<double> dStressdThick(NodesPerElement);
       GenVector<double> weight(NodesPerElement,0.0);
       int surface = senInfo[sindex].surface;
       elDisp->zero();       
       int flag;
       if(sinfo.getNLInfo().linearelastic) {
         // DETERMINE ELEMENT DISPLACEMENT VECTOR
         for (int k = 0; k < allDOFs->num(iele); ++k) {
           int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
           if (cn >= 0)
             (*elDisp)[k] = (*sol)[cn];
           else
             (*elDisp)[k] = 0; // XXX bcx[(*allDOFs)[iele][k]];
         } 
         flag = 1;
       } else {
         allCorot[iele]->extractDeformations(*geomState, nodes, elDisp->data(), flag);
         if(flag != 1) {
           filePrint(stderr,
                  " ... Error: nonlinear sensitivity for non-corotational elements is not implemented yet. exiting!\n");
           exit(-1);
          }
       }
       if(thgreleFlag[iele]) {     
         packedEset[iele]->getVonMisesThicknessSensitivity(dStressdThick, weight, nodes, *elDisp, 6, surface);  
       } else {
         dStressdThick.zero();
         Vector elstress(NodesPerElement,0.0);
         packedEset[iele]->getVonMises(elstress, weight, nodes, *elDisp, -1, surface, NULL, 0, 0, 0);
       }
       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           int inode;
           if(checkIsInStressNodes(node,inode)) {
             int iparam = thpaIndex[iele];
             if(iparam >= 0) (*allSens.vonMisesWRTthick)(inode, iparam) += dStressdThick[k]; 
             stressWeight(inode, 0) += weight[k]; 
           }
         }
       }
     } 

     for(int iparam = 0; iparam < numThicknessGroups; ++iparam) {
      // average vonMisesWRTthick vector
      for(int inode = 0; inode < numStressNodes; ++inode)  {
        if(stressWeight(inode, 0) == 0.0)
          (*allSens.vonMisesWRTthick)(inode, iparam) = 0.0;
        else
          (*allSens.vonMisesWRTthick)(inode, iparam) /= stressWeight(inode, 0);
      }
      if(!isDynam) {
        for(int inode = 0; inode < numStressNodes; ++inode)  {
          double a = allSens.lambdaStressVM[inode]->adjoint()*(allSens.linearstaticWRTthick[iparam]->col(0));
          (*allSens.vonMisesWRTthick)(inode, iparam) -= a;
        }
      }
     }
#endif
}

void 
Domain::computeNLStressVMWRTdisplacementSensitivity(int sindex,
                                                    AllSensitivities<double> &allSens, 
                                                    GeomState *geomState,
                                                    Corotator **allCorot)
{
#ifdef USE_EIGEN3
	allSens.vonMisesWRTdisp = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numNodes(), numUncon());
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> stressWeight(numNodes(), 1);
	allSens.vonMisesWRTdisp->setZero();     stressWeight.setZero();
	if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
	int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
	Vector *sol = 0;
	if(sinfo.getNLInfo().linearelastic) {
	  sol = new Vector(c_dsa->size());
	  geomState->get_tot_displacement(*sol);
	}
	for(int iele = 0; iele < numele; iele++) {
	  if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
	  int DofsPerElement = packedEset[iele]->numDofs();
	  int NodesPerElement = elemToNode->num(iele);
	  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DeformDispSen(DofsPerElement, DofsPerElement);
	  GenFullM<double> dStressdDisp(DofsPerElement,NodesPerElement,double(0.0));
	  GenVector<double> weight(NodesPerElement,0.0);
	  int surface = senInfo[sindex].surface;
	  elDisp->zero();
	  int flag;
	  if(sinfo.getNLInfo().linearelastic) {
		 // DETERMINE ELEMENT DISPLACEMENT VECTOR
		 for (int k = 0; k < allDOFs->num(iele); ++k) {
		   int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
		   if (cn >= 0)
			 (*elDisp)[k] = (*sol)[cn];
		   else
			 (*elDisp)[k] = 0; // XXX bcx[(*allDOFs)[iele][k]];
		 }
		 flag = 1;
		 packedEset[iele]->getVonMisesDisplacementSensitivity(dStressdDisp, weight, 0, nodes, *elDisp, 6, surface, 0);
	  } else {
		 allCorot[iele]->extractDeformations(*geomState, nodes, elDisp->data(), flag);
		 if(flag !=1 ) {
		   filePrint(stderr,
				   " ... Error: nonlinear sensitivity for non-corotational elements is not implemented yet. exiting!\n");
		   exit(-1);
		  }
		 allCorot[iele]->extractDeformationsDisplacementSensitivity(*geomState, nodes, DeformDispSen.data());
		 GenFullM<double> dDispDisp(DeformDispSen.data(),DofsPerElement,DofsPerElement);
		 packedEset[iele]->getVonMisesDisplacementSensitivity(dStressdDisp, weight, &dDispDisp, nodes, *elDisp, 6, surface, 0);
	  }

	  if(avgnum != 0) {
		 // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
		 const auto &unconstrNum = c_dsa->getUnconstrNum();
		 for(int k = 0; k < NodesPerElement; ++k) {
		   int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
		   auto dofs = (*allDOFs)[iele];
		   stressWeight(node,0) += weight[k];
		   for(int j = 0; j < DofsPerElement; ++j) {
			 int dofj = unconstrNum[dofs[j]];
			 if(dofs[j] < 0 || dofj < 0) continue;  // Skip undefined/constrained dofs
			 if(std::isnan(dStressdDisp[j][k])) std::cerr << "nan occurs in dStressdDisp[" << j << "]["
														  << k << "] with iele of " << iele << "\n";
			 (*allSens.vonMisesWRTdisp)(node, dofj) += dStressdDisp[j][k];
		   }
		 }
	  }
	}
	const auto &unconstrNum = c_dsa->getUnconstrNum();
	for(int inode = 0; inode < numNodes(); ++inode)  {
	  for(int j=0; j<nodeToNode->num(inode); ++j) { // loop over nodes connected to inode
		 int jnode = (*nodeToNode)[inode][j];
		 int jFirstDof = dsa->firstdof(jnode);
		 for(int k=0; k<dsa->weight(jnode); ++k) { // loop over dofs of jnode
		   int dof = unconstrNum[jFirstDof + k];
		   if(inode >= 0) {
			 if (stressWeight(inode, 0) != 0.0)
			   (*allSens.vonMisesWRTdisp)(inode,dof) /= stressWeight(inode,0);
			 else
			   (*allSens.vonMisesWRTdisp)(inode,dof) = 0;
		   }
		 }
	  }
	}
	if(!sol) delete sol;
#endif
}

void 
Domain::computeNormalizedNLVonMisesStress(GeomState &geomState, GeomState *refState, Corotator **allCorot,
                                          int surface, Vector &stress,  Vector &weight, bool normalized)
{

  //Vector weight(numNodes(),0.0);
  stress.zero();    weight.zero();

  int k;

  Vector elDisp(maxNumDOFs,0.0);
  Vector *sol = 0;

  if(sinfo.getNLInfo().linearelastic) {
    sol = new Vector(c_dsa->size());
    geomState.get_tot_displacement(*sol);
  }

  int iele;
  int NodesPerElement, maxNodesPerElement=0;
  for(iele=0; iele<numele; ++iele) {
    NodesPerElement = elemToNode->num(iele);
    maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
  }
  Vector elstress(maxNodesPerElement, 0.0);
  Vector elweight(maxNodesPerElement, 0.0);

  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  // Either get the nodal temperatures from the input file or
  // from the thermal model
  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >= 0) nodalTemperatures = temprcvd;

  int flag;
  for(iele = 0; iele < numElements(); ++iele) {
    if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isMpcElement()) continue;
    elDisp.zero();
    elstress.zero();
    elweight.zero();

    int NodesPerElement = packedEset[iele]->numNodes();

// extract deformations from current Geometry State of structure
    if(sinfo.getNLInfo().linearelastic) {
      // DETERMINE ELEMENT DISPLACEMENT VECTOR
      for (k = 0; k < allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if (cn >= 0)
          elDisp[k] = (*sol)[cn];
        else
          elDisp[k] = 0; // XXX bcx[(*allDOFs)[iele][k]];
      } 
      flag = 1;
    } else
    allCorot[iele]->extractDeformations(geomState, nodes, elDisp.data(), flag);

    // get element's nodal temperatures
    if(flag == 1) {
      elemNodeTemps.zero();
      if(nodalTemperatures) packedEset[iele]->nodes(nodeNumbers);
      for(int iNode = 0; iNode < NodesPerElement; ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
    }

    if (flag == 1) {
      // USE LINEAR STRESS ROUTINE
      packedEset[iele]->getVonMises(elstress, elweight, nodes,
                                    elDisp, 6, surface,
                                    elemNodeTemps.data(), 0, 0, 1);
    }
    else if (flag == 2) {
      // USE NON-LINEAR STRESS ROUTINE
      // note: in this case the element nodal temperatures are extracted from geomState inside the function
      allCorot[iele]->getNLVonMises(elstress, elweight, geomState,
                                    refState, nodes, 6, surface, 0, 0, 1);
    } else {
      // NO STRESS RECOVERY
    }

    // ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
    for(k=0; k<NodesPerElement; ++k) {
      stress[(*elemToNode)[iele][k]] += elstress[k];
      weight[(*elemToNode)[iele][k]] += elweight[k];
    }
  }

// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

  delete [] nodeNumbers;
  if(sol) delete sol;

  for(int k = 0; k < numNodes(); ++k)  {
    if(weight[k] == 0.0)
      stress[k] = 0.0;
    else {
      if(normalized) stress[k] /= (*aggregatedStress)*weight[k];
      else stress[k] /= weight[k];
    }
  }

}

void 
Domain::computeAggregatedNLStressVMWRTdisplacementSensitivity(int sindex,
                                                              AllSensitivities<double> &allSens, 
                                                              GeomState *geomState,
                                                              GeomState *refState,
                                                              Corotator **allCorot, Vector &stress,
                                                              Vector &stressWeight, Vector &sensCoef)
{
#ifdef USE_EIGEN3
     int surface = senInfo[sindex].surface;
     if(!(*aggregatedStressDenom)) computeAggregatedStressDenom(stress);
     allSens.aggregatedVonMisesWRTdisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numUncon());
     allSens.aggregatedVonMisesWRTdisp->setZero();
     if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);
     int avgnum = 1; //TODO: It is hardcoded to be 1, which corresponds to NODALFULL. It needs to be fixed.
     Vector *sol = 0;
     if(sinfo.getNLInfo().linearelastic) {
       sol = new Vector(c_dsa->size());
       geomState->get_tot_displacement(*sol);
     }
     
     for(int iele = 0; iele < numele; iele++) { 
       if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;
       int NodesPerElement = elemToNode->num(iele);
       int DofsPerElement = packedEset[iele]->numDofs();
       Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DeformDispSen(DofsPerElement, DofsPerElement);
       GenFullM<double> dStressdDisp(DofsPerElement,NodesPerElement,double(0.0));
       GenVector<double> elweight(NodesPerElement,0.0);
       elDisp->zero();       
       int flag;
       if(sinfo.getNLInfo().linearelastic) {
         // DETERMINE ELEMENT DISPLACEMENT VECTOR
         for (int k = 0; k < allDOFs->num(iele); ++k) {
           int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
           if (cn >= 0)
             (*elDisp)[k] = (*sol)[cn];
           else
             (*elDisp)[k] = 0; // XXX bcx[(*allDOFs)[iele][k]];
         } 
         flag = 1;
         packedEset[iele]->getVonMisesDisplacementSensitivity(dStressdDisp, elweight, 0, nodes, *elDisp, 6, surface, 0);
       } else {
         allCorot[iele]->extractDeformations(*geomState, nodes, elDisp->data(), flag);
         if(flag != 1) {
           filePrint(stderr,
              " ... Error: nonlinear sensitivity for non-corotational elements is not implemented yet. exiting!\n");
           exit(-1);
          }
         allCorot[iele]->extractDeformationsDisplacementSensitivity(*geomState, nodes, DeformDispSen.data());
         GenFullM<double> dDispDisp(DeformDispSen.data(),DofsPerElement,DofsPerElement);
         packedEset[iele]->getVonMisesDisplacementSensitivity(dStressdDisp, elweight, &dDispDisp, nodes,
                                                              *elDisp, 6, surface, 0);
       }

       if(avgnum != 0) {
         // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
         const auto &unconstrNum = c_dsa->getUnconstrNum();
         for(int k = 0; k < NodesPerElement; ++k) {
           int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
           auto dofs = (*allDOFs)[iele];
           for(int j = 0; j < DofsPerElement; ++j) {
             int dofj = unconstrNum[dofs[j]];
             if(dofs[j] < 0 || dofj < 0) continue;  // Skip undefined/constrained dofs
             if (stressWeight[node] != 0.0)
               (*allSens.aggregatedVonMisesWRTdisp)(dofj) += dStressdDisp[j][k]*sensCoef[node];
           }
         }
       }
     }   
#endif
}



void
Domain::makeNLPostSensitivities(GenSolver<DComplex> *sysSolver, 
                                AllSensitivities<DComplex> &allSens, 
                                GeomState *refState,
                                GeomState *geomState,
                                Corotator **allCorot,
                                bool isDynam) {}

void
Domain::makeNLPostSensitivities(GenSolver<double> *sysSolver, 
                                AllSensitivities<double> &allSens, 
                                GeomState *refState,
                                GeomState *geomState,
                                Corotator **allCorot,
                                bool isDynam)
{
#ifdef USE_EIGEN3
 makeThicknessGroupElementFlag();
 bool *includeStressNodes = new bool[numele];
 if(solInfo().sensitivityMethod == SolverInfo::Adjoint) setIncludeStressNodes(includeStressNodes); // TODO: need to include this
                                                                                                   // function in other sensitivity routine

 for(int sindex=0; sindex < numSensitivity; ++sindex) {
  switch(senInfo[sindex].type) {
   case SensitivityInfo::DisplacementWRTthickness: 
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       if(!allSens.linearstaticWRTthick) computeNLStaticWRTthicknessSensitivity(sindex,allSens,refState, geomState, allCorot);
       if(!isDynam) if(!allSens.dispWRTthick) computeDisplacementWRTthicknessDirectSensitivity(sindex, sysSolver, allSens);
     }
     else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
       if(!allSens.linearstaticWRTthick) computeNLStaticWRTthicknessSensitivity(sindex,allSens, refState, geomState, allCorot);
       if(!isDynam) { 
         if(!allSens.lambdaDisp) { 
           computeDisplacementDualSensitivity(sindex, sysSolver, allSens);
         }
         computeDisplacementWRTthicknessAdjointSensitivity(sindex, allSens);
         if (allSens.residual !=0 && allSens.dwrDisp==0) {
           allSens.dwrDisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numTotalDispDofs);
           for (int i=0; i<numTotalDispDofs; ++i) {
             (*allSens.dwrDisp)[i] = allSens.lambdaDisp[i]->dot(*(allSens.residual));
           }
         }
       }
     }
     break;
   }
   case SensitivityInfo::StressVMWRTthickness: 
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       if(!allSens.vonMisesWRTdisp) computeNLStressVMWRTdisplacementSensitivity(sindex,allSens,geomState,allCorot);
       if(!allSens.linearstaticWRTthick) computeNLStaticWRTthicknessSensitivity(sindex, allSens, refState, geomState, allCorot);
       if(!isDynam) if(!allSens.dispWRTthick) computeDisplacementWRTthicknessDirectSensitivity(sindex, sysSolver, allSens);
       computeNLStressVMWRTthicknessDirectSensitivity(sindex,allSens,geomState,allCorot,isDynam);
     }
     else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
       if(!allSens.vonMisesWRTdisp) computeNLStressVMWRTdisplacementSensitivity(sindex,allSens,geomState,allCorot);
       if(!allSens.linearstaticWRTthick) computeNLStaticWRTthicknessSensitivity(sindex, allSens, refState, geomState, allCorot);
       if(!isDynam) {
         if(!allSens.lambdaStressVM) {
           computeStressVMDualSensitivity(sindex, sysSolver, allSens);
         }
       }
       computeNLStressVMWRTthicknessAdjointSensitivity(sindex,allSens,geomState,allCorot,includeStressNodes,isDynam);
       if (allSens.residual !=0 && allSens.dwrStressVM==0 && !isDynam) {
         allSens.dwrStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numStressNodes);
         for (int i=0; i<numStressNodes; i++) {
           (*allSens.dwrStressVM)[i] = allSens.lambdaStressVM[i]->dot(*(allSens.residual));
         }
       }
     }
     break;
   } 

   case SensitivityInfo::AggregatedStressVMWRTthickness: 
   {
     if(solInfo().sensitivityMethod == SolverInfo::Direct) {
       filePrint(stderr,
         " ... WARNING : Only the adjoint method is available for KS function sensitivities. Switching to the adjoint method.\n");
     }
     if(aggregatedFlag) {
       getStressStrain(*geomState, allCorot,  aggregatedFileNumber, AGGREGATEDVON, 0.0, refState);
       aggregatedFlag = false;
     }
     int surface = senInfo[sindex].surface;
     Vector stress(numNodes(),0.0);
     Vector stressWeight(numNodes(),0.0);
     computeNormalizedNLVonMisesStress(*geomState, refState, allCorot, surface, stress, stressWeight); 
     if(!(*aggregatedStressDenom)) computeAggregatedStressDenom(stress);
     Vector sensCoef(numNodes(),0.0);
     for(int inode = 0; inode < numNodes(); ++inode) { 
       if (stressWeight[inode]!=0.0)
         sensCoef[inode] = exp(sinfo.ksParameter*(stress[inode])) / (stressWeight[inode]* (*aggregatedStressDenom));
     }
     if(!allSens.aggregatedVonMisesWRTdisp) {
       computeAggregatedNLStressVMWRTdisplacementSensitivity(sindex,allSens,geomState,refState,allCorot,stress,
                                                             stressWeight,sensCoef);
     }
     if(!allSens.linearstaticWRTthick) computeNLStaticWRTthicknessSensitivity(sindex, allSens, refState, geomState, allCorot);
     computeAggregatedNLStressVMWRTthicknessSensitivity(sindex,allSens,geomState,refState,allCorot,stress,
                                                        stressWeight,sensCoef,isDynam);
     if (allSens.residual !=0 && allSens.dwrAggregatedStressVM == 0 && !isDynam) {
       allSens.dwrAggregatedStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(1);
       (*allSens.dwrAggregatedStressVM)[0] = allSens.lambdaAggregatedStressVM->dot(*(allSens.residual));
     }
   } 
   default:
     break;
  }
 }
 delete [] includeStressNodes;
#endif
}

