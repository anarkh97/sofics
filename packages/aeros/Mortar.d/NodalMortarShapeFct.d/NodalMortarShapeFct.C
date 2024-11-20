// -----------------------------------------------------------------------------
// HB - 08/10/03
// -----------------------------------------------------------------------------
// HB - Modified 03/07/04 
//      -> add Dual flag in NodalMortarShapeFct::MakeSlaveLink(...)
//         & use it to force M to be a diagonal matrix 
//      -> clean & few optimization of NodalMortarShapeFct::BuildMortarLMPC(...)
// -----------------------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <limits>

// STL
#include <algorithm>
#include <vector>
#include <map>

// FEM headers
#include <Driver.d/Domain.h>
#include <Element.d/Element.h>
#include <Utils.d/Connectivity.h>

#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/NodalMortarShapeFct.d/NodalMortarShapeFct.h>
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS
// -----------------------------------------------------------------------------------------------------
NodalMortarShapeFct::NodalMortarShapeFct()
: NodalData()
, LinkedSlaveNodes()
, LinkedMasterNodes()
, SlaveMPCCoeffs()
, MasterMPCCoeffs() 
#ifdef USE_EIGEN3
, SlaveMPCCoeffDerivs()
, MasterMPCCoeffDerivs()
#endif
{ }

NodalMortarShapeFct::NodalMortarShapeFct(int ref_node, double ref_coeff)
: NodalData(1, node_pair_t(ref_node, ref_coeff))
, LinkedSlaveNodes()
, LinkedMasterNodes()
, SlaveMPCCoeffs()
, MasterMPCCoeffs() 
#ifdef USE_EIGEN3
, SlaveMPCCoeffDerivs()
, MasterMPCCoeffDerivs()
#endif
{ }

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTORS
// -----------------------------------------------------------------------------------------------------
NodalMortarShapeFct::~NodalMortarShapeFct() { }

// -----------------------------------------------------------------------------------------------------
//                                   INITILIZATION & CLEAN/CLEAR METHODS 
// -----------------------------------------------------------------------------------------------------
void 
NodalMortarShapeFct::ClearData()
{
  NodalData.clear();
  LinkedSlaveNodes.clear();
  LinkedMasterNodes.clear();
  SlaveMPCCoeffs.clear();
  MasterMPCCoeffs.clear();
#ifdef USE_EIGEN3
  SlaveMPCCoeffDerivs.resize(0,0);
  MasterMPCCoeffDerivs.resize(0,0);
#endif
}

void 
NodalMortarShapeFct::Reset()
{
  ClearData();
}

// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS
// -----------------------------------------------------------------------------------------------------
void 
NodalMortarShapeFct::SetRefData(int ref_node, double ref_coeff)
{
  NodalData.insert(NodalData.begin(), node_pair_t(ref_node, ref_coeff)); // ref node is first one
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------
//                                            MORTAR LMPC METHODS
// -----------------------------------------------------------------------------------------------------
void
NodalMortarShapeFct::MakeSlaveLink(Connectivity* SlaveNodeToFaces, FaceElemSet* FaceSet, bool Dual)
// ***************************************************************************************
// Look for the SLAVE nodes that belongs to the support of the nodal mortar shape fct
// -> SlaveNodeToFaces: global slave nodes numbering to indices of connected face elements 
//                      in the (ACTIVE) slave face element set (FaceSet)
//    FaceSet         : face element set of ACTIVE slave face elements 
// ***************************************************************************************
{
  // Step 1.: find the slave nodes connected with the nodal motar shape fct's nodes
  if(Dual){ // optimized for the dual mortar case -> FORCE M TO BE DIAGONAL
    LinkedSlaveNodes.push_back(GetRefNode());
    return;
  } else {
    LinkedSlaveNodes.reserve(8);
    for(int inode=0; inode<GetnNodes(); inode++) {
      int cnode = GetNodeId(inode);
      for(int iface=0, nFaces=SlaveNodeToFaces->num(cnode); iface<nFaces; iface++) {
        int jface = (*SlaveNodeToFaces)[cnode][iface]; 
        FaceElement* face = (*FaceSet)[jface];
        size_t old_size(LinkedSlaveNodes.size());
        LinkedSlaveNodes.insert(LinkedSlaveNodes.end(), face->nNodes(), 0);
        face->GetNodes(&LinkedSlaveNodes[old_size]);
      } 
    }
  }

  // eliminate dupplicate nodes
  std::sort(LinkedSlaveNodes.begin(), LinkedSlaveNodes.end());
  std::vector<int>::iterator Ilast(std::unique(LinkedSlaveNodes.begin(), LinkedSlaveNodes.end()));

  // shrink-to-fit (and automatically correctly resized) 
  std::vector<int>(LinkedSlaveNodes.begin(),Ilast).swap(LinkedSlaveNodes);
}

void
NodalMortarShapeFct::MakeMasterLink(Connectivity* SlaveNodeToFaces, FaceElemSet* FaceSet, Connectivity* SlaveToFFI, FFIPolygon** ContactPolygons)
// ***************************************************************************************
// Look for the MASTER nodes that belongs (throught FFIs) to the support of the nodal 
// mortar shape fct
// -> SlaveNodeToFaces: global slave nodes numbering to indices of connected face elements 
//                      in the (ACTIVE) slave face element set (FaceSet)
//    FaceSet         : face element set of ACTIVE slave face elements 
//    SlaveToFFI      : slave face element (in ACTIVE numbering) to indices in the
//                      ContactPolygons array of the supported FFIs 
//    ContactPolygons : array of the FFIs
// ***************************************************************************************
{
  LinkedMasterNodes.reserve(8);

  for(int inode=0; inode<GetnNodes(); inode++) {
    int cnode = GetNodeId(inode);
    for(int iface=0, nFaces=SlaveNodeToFaces->num(cnode); iface<nFaces; iface++){
      int jface = (*SlaveNodeToFaces)[cnode][iface]; 
      for(int iFFI=0, nFFIs=SlaveToFFI->num(jface); iFFI<nFFIs; iFFI++) {
        int jFFI = (*SlaveToFFI)[jface][iFFI];
        FaceElement* face = ContactPolygons[jFFI]->GetPtrMasterFace();
        size_t old_size(LinkedMasterNodes.size());
        LinkedMasterNodes.insert(LinkedMasterNodes.end(), face->nNodes(), 0);
        face->GetNodes(&LinkedMasterNodes[old_size]);
      }
    } 
  }
  // eliminate dupplicate nodes
  std::sort(LinkedMasterNodes.begin(), LinkedMasterNodes.end());
  std::vector<int>::iterator Ilast(std::unique(LinkedMasterNodes.begin(),LinkedMasterNodes.end()));
  
  // shrink-to-fit (and automatically correctly resized) 
  std::vector<int>(LinkedMasterNodes.begin(),Ilast).swap(LinkedMasterNodes);
}

void
NodalMortarShapeFct::BuildMortarLMPC(Connectivity* SlaveNodeToFaces, FaceElemSet* FaceSet, Connectivity* SlaveToFFI, FFIPolygon** ContactPolygons)
// ***************************************************************************************
// Build the mortar LMPC coefficients by "assembling" the contributions of all the FFIs
// belongins to the support of the nodal mortar shape fct 
// -> compute the row of M and N associated with the current nodal mortar shape fct
// -> SlaveNodeToFaces: global slave nodes numbering to indices of connected face elements 
//                      in the (ACTIVE) slave face element set (FaceSet)
//    FaceSet         : face element set of ACTIVE slave face elements 
//    SlaveToFFI      : slave face element (in ACTIVE numbering) to indices in the
//                      ContactPolygons array of the supported FFIs 
//    ContactPolygons : array of the FFIs
// ***************************************************************************************
{
  SlaveMPCCoeffs.assign(LinkedSlaveNodes.size(), 0.0);
  MasterMPCCoeffs.assign(LinkedMasterNodes.size(), 0.0);
  
  for(int inode=0; inode<GetnNodes(); inode++){
    int cnode = GetNodeId(inode);
    double Aq = GetNodalCoeff(inode); 
    for(int iface=0, nFaces=SlaveNodeToFaces->num(cnode); iface<nFaces; iface++){
      int jface = (*SlaveNodeToFaces)[cnode][iface];
      for(int iFFI=0, nFFIs=SlaveToFFI->num(jface); iFFI<nFFIs; iFFI++){
        int jFFI = (*SlaveToFFI)[jface][iFFI];
        FFIPolygon* FFI = ContactPolygons[jFFI];
        FaceElement* slaveface  = FFI->GetPtrSlaveFace();
        FullM* M = FFI->GetPtrM(); 
        int i = slaveface->GetNodeIndex(cnode); // with a trick this could also be moved outside the iFFI loop
        // can be optimized if dual mortar space used by not looping over all the slaveface's node 
        for(int j=0, nNds=slaveface->nNodes(); j<nNds; j++) {
          std::vector<int>::iterator I(std::lower_bound(LinkedSlaveNodes.begin(), 
                                                        LinkedSlaveNodes.end(), slaveface->GetNode(j)));
          if(I != LinkedSlaveNodes.end()) {
            int k(std::distance(LinkedSlaveNodes.begin(), I));
            SlaveMPCCoeffs[k] += Aq*(*M)[i][j];
          }
        }
        FaceElement* masterface = FFI->GetPtrMasterFace();
        FullM* N = FFI->GetPtrN(); 
        for(int j=0, nNds=masterface->nNodes(); j<nNds; j++) {
          int k(std::distance(LinkedMasterNodes.begin(), 
                              std::lower_bound(LinkedMasterNodes.begin(), 
                                               LinkedMasterNodes.end(), masterface->GetNode(j))));
          MasterMPCCoeffs[k] += Aq*(*N)[i][j]; 
        }
      }
    }
  }
}

void
NodalMortarShapeFct::BuildMortarCtcLMPC(Connectivity* SlaveNodeToFaces, FaceElemSet* FaceSet, Connectivity* SlaveToFFI, FFIPolygon** ContactPolygons)
// ***************************************************************************************
// Build the mortar LMPC coefficients by "assembling" the contributions of all the FFIs
// belonging to the support of the nodal mortar shape fct 
// -> compute the row of M and N associated with the current nodal mortar shape fct
// -> SlaveNodeToFaces: global slave nodes numbering to indices of connected face elements 
//                      in the (ACTIVE) slave face element set (FaceSet)
//    FaceSet         : face element set of ACTIVE slave face elements 
//    SlaveToFFI      : slave face element (in ACTIVE numbering) to indices in the
//                      ContactPolygons array of the supported FFIs 
//    ContactPolygons : array of the FFIs
// ***************************************************************************************
{
  SlaveMPCCoeffs.assign(3*LinkedSlaveNodes.size(), 0.0);
  MasterMPCCoeffs.assign(3*LinkedMasterNodes.size(), 0.0); 
  MPCRhs = 0;
#if defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD) || defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
  SlaveMPCCoeffDerivs.resize(3*LinkedSlaveNodes.size(),3*LinkedSlaveNodes.size()+3*LinkedMasterNodes.size());
  SlaveMPCCoeffDerivs.setZero();
  MasterMPCCoeffDerivs.resize(3*LinkedMasterNodes.size(),3*LinkedSlaveNodes.size()+3*LinkedMasterNodes.size());
  MasterMPCCoeffDerivs.setZero();
#endif

  for(int inode=0; inode<GetnNodes(); inode++){
    int cnode = GetNodeId(inode);
    double Aq = GetNodalCoeff(inode); 
    for(int iface=0, nFaces=SlaveNodeToFaces->num(cnode); iface<nFaces; iface++){
      int jface = (*SlaveNodeToFaces)[cnode][iface];
      for(int iFFI=0, nFFIs=SlaveToFFI->num(jface); iFFI<nFFIs; iFFI++){
        int jFFI = (*SlaveToFFI)[jface][iFFI];
        FFIPolygon* FFI = ContactPolygons[jFFI];
        FaceElement* slaveface  = FFI->GetPtrSlaveFace();
        FaceElement* masterface = FFI->GetPtrMasterFace();
        FullM* M = FFI->GetPtrM(); 
        FullM* dM = FFI->GetdM();
        int i = slaveface->GetNodeIndex(cnode); // with a trick this could also be moved outside the iFFI loop
        // can be optimized if dual mortar space used by not looping over all the slaveface's node 
        for(int j=0; j<slaveface->nNodes(); j++) {
          std::vector<int>::iterator I(std::lower_bound(LinkedSlaveNodes.begin(), 
                                                        LinkedSlaveNodes.end(), slaveface->GetNode(j)));
          if(I != LinkedSlaveNodes.end()) {
            int k(std::distance(LinkedSlaveNodes.begin(), I));
            SlaveMPCCoeffs[3*k]   += Aq*(*M)[i][3*j];
            SlaveMPCCoeffs[3*k+1] += Aq*(*M)[i][3*j+1];
            SlaveMPCCoeffs[3*k+2] += Aq*(*M)[i][3*j+2];
#if defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD) || defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
            for(int l=0; l<slaveface->nNodes(); l++) {
              std::vector<int>::iterator J(std::lower_bound(LinkedSlaveNodes.begin(),
                                                            LinkedSlaveNodes.end(), slaveface->GetNode(l)));
              if(J != LinkedSlaveNodes.end()) {
                int m(std::distance(LinkedSlaveNodes.begin(), J));
                for(int p=0; p<3; ++p)
                  for(int q=0; q<3; ++q)
                    SlaveMPCCoeffDerivs(3*k+p,3*m+q) += Aq*dM[i][3*j+p][MAX_MORTAR_DERIVATIVES/2+3*l+q];
              }
            }
            for(int l=0; l<masterface->nNodes(); l++) {
              int m(std::distance(LinkedMasterNodes.begin(),
                                  std::lower_bound(LinkedMasterNodes.begin(),
                                                   LinkedMasterNodes.end(), masterface->GetNode(l))));
              for(int p=0; p<3; ++p)
                for(int q=0; q<3; ++q)
                  SlaveMPCCoeffDerivs(3*k+p,3*LinkedSlaveNodes.size()+3*m+q) += Aq*dM[i][3*j+p][3*l+q];
            }
#endif
          }
        }
        FullM* N = FFI->GetPtrN(); 
        FullM* dN = FFI->GetdN();
        for(int j=0; j<masterface->nNodes(); j++) {
          int k(std::distance(LinkedMasterNodes.begin(), 
                              std::lower_bound(LinkedMasterNodes.begin(), 
                                               LinkedMasterNodes.end(), masterface->GetNode(j))));
          MasterMPCCoeffs[3*k]   += Aq*(*N)[i][3*j];
          MasterMPCCoeffs[3*k+1] += Aq*(*N)[i][3*j+1];
          MasterMPCCoeffs[3*k+2] += Aq*(*N)[i][3*j+2];
#if defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD) || defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
          for(int l=0; l<slaveface->nNodes(); l++) {
            std::vector<int>::iterator J(std::lower_bound(LinkedSlaveNodes.begin(),
                                                          LinkedSlaveNodes.end(), slaveface->GetNode(l)));
            if(J != LinkedSlaveNodes.end()) {
              int m(std::distance(LinkedSlaveNodes.begin(), J));
              for(int p=0; p<3; ++p)
                for(int q=0; q<3; ++q)
                  MasterMPCCoeffDerivs(3*k+p,3*m+q) += Aq*dN[i][3*j+p][MAX_MORTAR_DERIVATIVES/2+3*l+q];
            }
          }
          for(int l=0; l<masterface->nNodes(); l++) {
            int m(std::distance(LinkedMasterNodes.begin(),
                                std::lower_bound(LinkedMasterNodes.begin(),
                                                 LinkedMasterNodes.end(), masterface->GetNode(l))));
            for(int p=0; p<3; ++p)
              for(int q=0; q<3; ++q)
                MasterMPCCoeffDerivs(3*k+p,3*LinkedSlaveNodes.size()+3*m+q) += Aq*dN[i][3*j+p][3*l+q];
          }
#endif
        }
        Vector* g = FFI->GetPtrNormalGeoGaps();
        MPCRhs -= Aq*(*g)[i];
      }
    }
  }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS
// -----------------------------------------------------------------------------------------------------
void
NodalMortarShapeFct::print(int* SlaveLlToGlNodeMap, int* MasterLlToGlNodeMap)
{
   fprintf(stderr,"---------------------------------------\n");
   fprintf(stderr,"NodalMortarShapeFct object:\n");
   fprintf(stderr,"-> nb nodal fcts : %d\n",GetnNodes());
   int node = (SlaveLlToGlNodeMap) ? SlaveLlToGlNodeMap[GetRefNode()] : GetRefNode();
   fprintf(stderr,"-> ref. node     : %d\n",node);
   for(int i=0; i<GetnNodes(); i++) {
     node = (SlaveLlToGlNodeMap) ? SlaveLlToGlNodeMap[GetNodeId(i)] : GetNodeId(i);
     fprintf(stderr,"-> node %d = %d, coeff = %e\n",i+1,node,GetNodalCoeff(i));
   }
   if(!LinkedSlaveNodes.empty()) {
     fprintf(stderr,"-> linked slave nodes: \n");
     for(int i=0; i<int(LinkedSlaveNodes.size()); i++) {
       node = (SlaveLlToGlNodeMap) ? SlaveLlToGlNodeMap[LinkedSlaveNodes[i]] : LinkedSlaveNodes[i];
       std::cerr << node << ' ';
     }
     std::cerr << std::endl;
   }
   if(!LinkedMasterNodes.empty()) {
     fprintf(stderr,"-> linked master nodes: \n");
     for(int i=0; i<int(LinkedMasterNodes.size()); i++) {
       node = (MasterLlToGlNodeMap) ? MasterLlToGlNodeMap[LinkedMasterNodes[i]] : LinkedMasterNodes[i];
       std::cerr << node << ' ';
     }
     std::cerr << std::endl;
   }
   if(!SlaveMPCCoeffs.empty()) {
     fprintf(stderr,"-> slave mpc coeffs: \n");
     double sum = 0.0;
     for(int i=0; i<int(LinkedSlaveNodes.size()); i++) {
       node = (SlaveLlToGlNodeMap) ? SlaveLlToGlNodeMap[LinkedSlaveNodes[i]] : LinkedSlaveNodes[i];
       std::cerr << " slave node " << node << ", coeff = " << SlaveMPCCoeffs[i] << std::endl;
       sum += SlaveMPCCoeffs[i];
     }
     std::cerr <<" sum slave coeffs = "<< sum << std::endl;
   }
   if(!MasterMPCCoeffs.empty()) {
     fprintf(stderr,"-> master mpc coeffs: \n");
     double sum = 0.0;
     for(int i=0; i<int(LinkedMasterNodes.size()); i++) {
       node = (MasterLlToGlNodeMap) ? MasterLlToGlNodeMap[LinkedMasterNodes[i]] : LinkedMasterNodes[i];
       std::cerr << " master node " << node << ", coeff = " << MasterMPCCoeffs[i] << std::endl;
       sum += MasterMPCCoeffs[i];
     }    
     std::cerr << " sum master coeffs = "<< sum << std::endl;
   }
   fprintf(stderr,"---------------------------------------\n");
}

LMPCons*
NodalMortarShapeFct::CreateMortarLMPCons(int lmpcnum, int dof, double rhs, 
                                         int* SlaveLlToGlNodeMap, int* MasterLlToGlNodeMap)
// ***************************************************************************************
// For reusing the standard LMPC code: create the LMPConstrain associated with the
// current nodal mortar shape fct
// ***************************************************************************************
{
  LMPCTerm SlaveTerm;

  LMPCons* MortarLMPC = NULL;
  double tol = 0.0; // we may use a (relative) tolerance to filter the small term
  for(int i = 0; i < int(LinkedSlaveNodes.size()); i++){
    if(fabs(SlaveMPCCoeffs[i]) > tol) {
      SlaveTerm.nnum   = SlaveLlToGlNodeMap ? SlaveLlToGlNodeMap[LinkedSlaveNodes[i]] : LinkedSlaveNodes[i];
      SlaveTerm.dofnum = dof;
      SlaveTerm.coef.r_value = SlaveMPCCoeffs[i];
      if(MortarLMPC == NULL) {
        MortarLMPC = new LMPCons(lmpcnum, rhs, &SlaveTerm);
        MortarLMPC->type = 0;
        MortarLMPC->setType(mpc::Equality);
        MortarLMPC->setSource(mpc::TiedSurfaces);
      }
      else MortarLMPC->addterm(&SlaveTerm);
    }
  }

 LMPCTerm MasterTerm;
 for(int i = 0; i < int(LinkedMasterNodes.size()); i++){
   if(fabs(MasterMPCCoeffs[i]) > tol) {
     MasterTerm.nnum   = MasterLlToGlNodeMap ? MasterLlToGlNodeMap[LinkedMasterNodes[i]] : LinkedMasterNodes[i];
     MasterTerm.dofnum = dof;
     MasterTerm.coef.r_value = -MasterMPCCoeffs[i];
     if(MortarLMPC == NULL) {
        MortarLMPC = new LMPCons(lmpcnum, rhs, &MasterTerm);
        MortarLMPC->type = 0;
        MortarLMPC->setType(mpc::Equality);
        MortarLMPC->setSource(mpc::TiedSurfaces);
      }
      else MortarLMPC->addterm(&MasterTerm);
    }
  }
#ifdef MORTAR_WARNING
  if(!created){
    fprintf(stderr," ### WARNING: in NodalMortarShapeFct::CreateMortarLMPCons(...): no MortarLMPC has been created !!!\n");
    print();
  }
#endif
  return MortarLMPC;
}

LMPCons*
NodalMortarShapeFct::CreateMortarCtcLMPCons(int lmpcnum, int* SlaveLlToGlNodeMap, int* MasterLlToGlNodeMap, int mode)
{
  double rhs = MPCRhs;
  LMPCTerm SlaveTerm;
                                                                                                                                             
  LMPCons* MortarLMPC = NULL;
  double tol = 0.0; // we may use a (relative) tolerance to filter the small term
  int dofs[3] = {0,1,2};
  std::vector<int> indices;
  for(int i = 0; i < int(LinkedSlaveNodes.size()); i++) {
    for(int idof = 0; idof < 3; idof++) {
// filter based on both the mpc coef AND the row/col of the derivative matrix all zero. I think it is possible that the
// coef is zero but the derivative matrix is non-zero. in this case we should keep the zero coef term.
#if defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD) || defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
      if(fabs(SlaveMPCCoeffs[3*i+idof]) > tol || SlaveMPCCoeffDerivs.row(3*i+idof).norm() > tol) {
#else
      if(fabs(SlaveMPCCoeffs[3*i+idof]) > tol) {
#endif
        indices.push_back(3*i+idof);
        SlaveTerm.nnum   = SlaveLlToGlNodeMap ? SlaveLlToGlNodeMap[LinkedSlaveNodes[i]] : LinkedSlaveNodes[i];
        SlaveTerm.dofnum = dofs[idof];
        SlaveTerm.coef.r_value = SlaveMPCCoeffs[3*i+idof];
        if(MortarLMPC == NULL) {
          MortarLMPC = new LMPCons(lmpcnum, rhs, &SlaveTerm);
          if(mode == 0) {
            MortarLMPC->type = 0;
            MortarLMPC->setType(mpc::Equality);
          }
          else {
            MortarLMPC->type = 1;
            MortarLMPC->setType(mpc::Inequality);
          }
          MortarLMPC->setSource(mpc::ContactSurfaces);
        } else
          MortarLMPC->addterm(&SlaveTerm);
      }
    }
  }

 LMPCTerm MasterTerm;
 for(int i = 0; i < int(LinkedMasterNodes.size()); i++) {
   for(int idof = 0; idof < 3; idof++){
#if defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD) || defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
     if(fabs(MasterMPCCoeffs[3*i+idof]) > tol || MasterMPCCoeffDerivs.row(3*i+idof).norm() > tol) {
#else
     if(fabs(MasterMPCCoeffs[3*i+idof]) > tol) {
#endif
       indices.push_back(3*LinkedSlaveNodes.size()+3*i+idof);
       MasterTerm.nnum   = MasterLlToGlNodeMap ? MasterLlToGlNodeMap[LinkedMasterNodes[i]] : LinkedMasterNodes[i];
       MasterTerm.dofnum = dofs[idof];
       MasterTerm.coef.r_value = -MasterMPCCoeffs[3*i+idof];
       if(MortarLMPC == NULL) {
         MortarLMPC = new LMPCons(lmpcnum, rhs, &MasterTerm);
         MortarLMPC->type = 1; // this is to be phased out
         MortarLMPC->setType(mpc::Inequality);
         MortarLMPC->setSource(mpc::ContactSurfaces);
       } else
         MortarLMPC->addterm(&MasterTerm);
     }
   }
 }
 if(MortarLMPC) MortarLMPC->id.second = SlaveLlToGlNodeMap[GetNodeId(0)];

#if defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD) || defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
/*
 int n = 3*LinkedSlaveNodes.size()+3*LinkedMasterNodes.size();
 MortarLMPC->H.resize(n,n);
 MortarLMPC->H << SlaveMPCCoeffDerivs, -MasterMPCCoeffDerivs;
*/
 if(MortarLMPC) {
   int n = 3*LinkedSlaveNodes.size()+3*LinkedMasterNodes.size();
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hfull(n,n);
   Hfull << SlaveMPCCoeffDerivs, -MasterMPCCoeffDerivs;

   int m = indices.size();
   MortarLMPC->H.resize(m,m);
   for(int i=0; i<m; ++i)
     for(int j=0; j<m; ++j)
       MortarLMPC->H(i,j) = Hfull(indices[i],indices[j]);
 }
#endif
 return MortarLMPC;
}

void
NodalMortarShapeFct::BuildWetFSICoupling(Connectivity* SlaveNodeToFaces, FaceElemSet* FaceSet, 
                                         Connectivity* SlaveToFFI, FFIPolygon** ContactPolygons)
// ***************************************************************************************
// Build the wet FSI LMPC coefficients by "assembling" the contributions of all the FFIs
// belongins to the support of the nodal mortar shape fct 
// -> compute the row of M and N associated with the current nodal mortar shape fct
// -> SlaveNodeToFaces: global slave nodes numbering to indices of connected face elements 
//                      in the (ACTIVE) slave face element set (FaceSet)
//    FaceSet         : face element set of ACTIVE slave face elements 
//    SlaveToFFI      : slave face element (in ACTIVE numbering) to indices in the
//                      ContactPolygons array of the supported FFIs 
//    ContactPolygons : array of the FFIs
// ***************************************************************************************
{
  MasterMPCCoeffs.assign(3*LinkedMasterNodes.size(), 0.0);

  for(int inode=0; inode<GetnNodes(); inode++){
    int cnode = GetNodeId(inode);
    double Aq = GetNodalCoeff(inode); 
    for(int iface=0, nFaces=SlaveNodeToFaces->num(cnode); iface<nFaces; iface++){
      int jface = (*SlaveNodeToFaces)[cnode][iface];
      for(int iFFI=0, nFFIs=SlaveToFFI->num(jface); iFFI<nFFIs; iFFI++){
        int jFFI = (*SlaveToFFI)[jface][iFFI];
        FFIPolygon* FFI = ContactPolygons[jFFI];
        FaceElement* slaveface  = FFI->GetPtrSlaveFace();
        FaceElement* masterface = FFI->GetPtrMasterFace();
        FullM* N = FFI->GetPtrN(); 
        int i = slaveface->GetNodeIndex(cnode); // with a trick this could also be moved outside the iFFI loop
        for(int j=0, nNds=masterface->nNodes(); j<nNds; j++){
          int k(std::distance(LinkedMasterNodes.begin(), 
                              std::lower_bound(LinkedMasterNodes.begin(), 
                                               LinkedMasterNodes.end(), masterface->GetNode(j))));
        
          MasterMPCCoeffs[3*k  ] += Aq*(*N)[i][3*j  ]; 
          MasterMPCCoeffs[3*k+1] += Aq*(*N)[i][3*j+1]; 
          MasterMPCCoeffs[3*k+2] += Aq*(*N)[i][3*j+2]; 
        }
      }
    }
  }
}

LMPCons*
NodalMortarShapeFct::CreateWetFSICons(int* SlaveLlToGlNodeMap, int* MasterLlToGlNodeMap)
// ***************************************************************************************
// For reusing the standard LMPC code: create the LMPConstraint associated with the
// current nodal mortar shape fct
// ***************************************************************************************
{
  LMPCons* WetFSICons = NULL;

  int lmpcnum = SlaveLlToGlNodeMap ? SlaveLlToGlNodeMap[GetRefNode()] : GetRefNode();
  int dofs[3] = {0,1,2};
  LMPCTerm MasterTerm;
  bool created = false;
  for(int i=0; i<int(LinkedMasterNodes.size()); i++){
    for(int idof=0; idof<3; idof++){
      if(MasterMPCCoeffs[3*i+idof]!=0.0) { // skip zero coeff. (we may use a (relative) tolerance to filter the small term)
        MasterTerm.nnum   = MasterLlToGlNodeMap ? MasterLlToGlNodeMap[LinkedMasterNodes[i]] : LinkedMasterNodes[i];
        MasterTerm.dofnum = dofs[idof];
        MasterTerm.coef.r_value = -MasterMPCCoeffs[3*i+idof]; // minus because we compute the coeffs. using the normal to the fluid
        if(!created){
          WetFSICons = new LMPCons(lmpcnum, 0.0, &MasterTerm);
          created = true;
        } else
          WetFSICons->addterm(&MasterTerm);
      }
    }
  }
#ifdef MORTAR_WARNING
  if(!created){
    fprintf(stderr," ### WARNING: in NodalMortarShapeFct::CreateWetFSICons(...): no WetFSICons has been created !!!\n");
    print();
  }
#endif
  return WetFSICons;
}

