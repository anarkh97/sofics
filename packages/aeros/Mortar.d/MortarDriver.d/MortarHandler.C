// -----------------------------------------------------------------------
// HB - 07/05/03
// -----------------------------------------------------------------------
// Last modif: 09/03/03 -> IMPLEMENT CASE OF QUADRATIC FACE ELEMENT
//                         USE UNDERLAYING LINEAR FACE ELEMENT FOR
//                         PERFORMING ACME FFI SEARCH
// -----------------------------------------------------------------------
//             09/09/03 -> In AddMortarLMPCs: change the (dof) order 
//                         when the Mortar LMPCs are added to the
//                         standard LMPC array (lmpc)
// -----------------------------------------------------------------------
//             09/29/03 -> implement the USE_ACME
//                         we pass the default MPI_COMM_SELF
//                         as MPI communicator to ACME: every MPI
//                         process do the SAME job
// -----------------------------------------------------------------------
//             09/29/03 -> add the USE_ACME flag: prevent the
//                         include/use/call to ACME class & fcts
// -----------------------------------------------------------------------
//             11/20/03 -> change DISTRIBUTED flag to USE_MPI flag
// -----------------------------------------------------------------------
//             03/02/04 -> add HB_THREAD_FFI_M_N flag: defined in
//                         MortarHandlerDefs.h -> enable the use of
//                         thread level parallelization of the comp-
//                         putation of the FFI M & N contributions
//                         in MortarHandler::CreateFFIPolygons(...)
// -----------------------------------------------------------------------
//             03/02/04 -> add HB_MORTAR_TIMER flag: defined in
//                         MortarHandlerDefs.h -> enable the timings
//                         in MortarHandler::CreateFFIPolygons(...) 
// -----------------------------------------------------------------------
//             03/03/04 -> add HB_THREAD_NODALMORTAR flag: defined in
//                         MortarHandlerDefs.h -> enable the use of
//                         thread level parallelization of the creation 
//                         of the mortar LMPC at each (active) slave nodes
//                         see MortarHandler::CreateFFIPolygons(...)
//                         -> DO NOT SEEM TO BE EFFICIENT ...
// -----------------------------------------------------------------------
//             09/15/04 -> add MORTAR_LOCALNUMBERING flag: defined in
//                         the Makefile -> for the case when the surface
//                         entities have been renumbered in local 
//                         numbering AND have their OWN node set coord. 
// -----------------------------------------------------------------------
//             02/19/07 -> FFI methods taking the nodes'coord. set as input 
//                         are now DEPRECATED. The  nodes'coord. set as now
//                         accessed throught the SurfaceEntity objects. 
//                         This modification allows simplifications: only 
//                         one implementation for the case of global or
//                         local renumbering of the SurfaceEntities objects.
//                         The compiler MORTAR_LOCALNUMBERING flag is still
//                         needed because the selection between global or
//                         local renumbering of the SurfaceEntities objects
//                         is still done at compile time (see SurfaceEntity.h) 
// -----------------------------------------------------------------------
// Std C/C++ headers 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

// STL
#include <algorithm>
#include <vector>
#include <map>

#ifdef USE_MPI 
#include <mpi.h>
#endif

// FEM headers
#include <Driver.d/Domain.h>
#include <Math.d/matrix.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/resize_array.h>
#include <Utils.d/DistHelper.h>
#include <Paral.d/SubDOp.h>

#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/NodalMortarShapeFct.d/NodalMortarShapeFct.h>

// ACME headers
#ifdef USE_ACME
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactEnforcement.h"
#include "Search_Interface.h"
#include "Enforcement_Interface.h"
#endif

// Locally define flags
#include <Mortar.d/MortarDriver.d/MortarHandlerDefs.h>

#ifdef HB_MORTAR_TIMER
#include <Timers.d/GetTime.h>
#endif

extern Communicator *structCom;

extern Domain *domain;

extern MortarElement* CreateMortarElement(FaceElement*, CoordSet&, bool);

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS
// -----------------------------------------------------------------------------------------------------
MortarHandler::MortarHandler()
: MortarEls()
, NodalMortars()
, ActiveSlaveNodes()
, ActiveSlaveFacesToMortarEl()
{
  Initialize();
}

MortarHandler::MortarHandler(int _Id)
: MortarEls()
, NodalMortars()
, ActiveSlaveNodes()
, ActiveSlaveFacesToMortarEl()
{
  Initialize();
  Id = _Id;
}

MortarHandler::MortarHandler(int _Id, int _MasterEntityId, int _SlaveEntityId)
: MortarEls()
, NodalMortars()
, ActiveSlaveNodes()
, ActiveSlaveFacesToMortarEl()
{
  Initialize();
  Id             = _Id;
  MasterEntityId = _MasterEntityId;
  SlaveEntityId  = _SlaveEntityId ;
}

MortarHandler::MortarHandler(int _MasterEntityId, int _SlaveEntityId)
: MortarEls()
, NodalMortars()
, ActiveSlaveNodes()
, ActiveSlaveFacesToMortarEl()
{
  Initialize();
  MasterEntityId = _MasterEntityId;
  SlaveEntityId  = _SlaveEntityId ;
}

MortarHandler::MortarHandler(int _MasterEntityId, int _SlaveEntityId, 
                             double _NormalTol)
: MortarEls()
, NodalMortars()
, ActiveSlaveNodes()
, ActiveSlaveFacesToMortarEl()
{
  Initialize();
  MasterEntityId = _MasterEntityId;
  SlaveEntityId  = _SlaveEntityId ;
  NormalSearchTol= _NormalTol;
}

MortarHandler::MortarHandler(int _MasterEntityId, int _SlaveEntityId, 
                             double _NormalTol, double _TangTol)
: MortarEls()
, NodalMortars()
, ActiveSlaveNodes()
, ActiveSlaveFacesToMortarEl()
{
  Initialize();
  MasterEntityId = _MasterEntityId;
  SlaveEntityId  = _SlaveEntityId ;
  NormalSearchTol= _NormalTol;
  TangSearchTol  = _TangTol;
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTORS
// -----------------------------------------------------------------------------------------------------
MortarHandler::~MortarHandler()
{
  Id             = -1;
  MasterEntityId = -1;
  SlaveEntityId  = -1;

  PtrMasterEntity = 0;
  PtrSlaveEntity  = 0;

  MortarType      = MortarHandler::STD;
  InteractionType = MortarHandler::TIED;
 
  nFFI         = 0;
  nACMEFFIData = 0; 
  if(Slave_face_block_id)        { delete [] Slave_face_block_id       ; Slave_face_block_id       = 0; }
  if(Slave_face_index_in_block)  { delete [] Slave_face_index_in_block ; Slave_face_index_in_block = 0; }
  if(Master_face_block_id)       { delete [] Master_face_block_id      ; Master_face_block_id      = 0; }
  if(Master_face_index_in_block) { delete [] Master_face_index_in_block; Master_face_index_in_block= 0; }
  if(Slave_face_procs)           { delete [] Slave_face_procs          ; Slave_face_procs          = 0; }
  if(Master_face_procs)          { delete [] Master_face_procs         ; Master_face_procs         = 0; }
  if(ACMEFFI_index)              { delete [] ACMEFFI_index             ; ACMEFFI_index             = 0; }
  if(ACMEFFI_data)               { delete [] ACMEFFI_data              ; ACMEFFI_data              = 0; }

  for(size_t i=0; i<MortarEls.size(); ++i) 
    if(MortarEls[i]) delete MortarEls[i];

  for(std::vector<FFIPolygon*>::iterator it = CtcPolygons.begin(); it != CtcPolygons.end(); ++it) {
    delete *it;
  }

  if(ActiveSlaveNodeToElem) { delete ActiveSlaveNodeToElem ; ActiveSlaveNodeToElem  = 0; }
  if(ActiveMasterNodeToElem){ delete ActiveMasterNodeToElem; ActiveMasterNodeToElem = 0; }
  if(SlaveFaceToFFIConnect) { delete SlaveFaceToFFIConnect ; SlaveFaceToFFIConnect  = 0; }

  // this only delete the array of pointers to the (active) master/slave elements
  // but not the element themselves (i.e. they are currently owned by the SurfaceEntity objects) 
  ActiveSlaveElemSet.deleteElems();
  ActiveMasterElemSet.deleteElems();
 
  nMortarLMPCs = 0;
  gIdFirstLMPC = 0;
  gIdLastLMPC  = 0;
#ifdef USE_ACME
  if(search_obj) delete search_obj;
  if(contact_obj) delete contact_obj;
#endif
  if(mass) delete [] mass;
  if(dofmap) delete [] dofmap;
  if(node_global_ids) delete [] node_global_ids;
  if(data) delete [] data;
  if(share) delete [] share;
  if(density) delete [] density;
  if(wavespeed) delete [] wavespeed;
  if(comm_proc_id) delete [] comm_proc_id;
  if(number_nodes_to_partner) delete [] number_nodes_to_partner;
  if(comm_node) delete [] comm_node; 
#ifdef DISTRIBUTED
  if(DIST_ACME == 2) {
    delete PtrMasterEntity;
    delete PtrSlaveEntity;
  }
#endif
  if(ConstraintOptionsData) delete ConstraintOptionsData;
}

// -----------------------------------------------------------------------------------------------------
//                                   INITILIZATION & CLEAN/CLEAR METHODS
// -----------------------------------------------------------------------------------------------------
void
MortarHandler::Initialize()
{
  Id             = -1;
  MasterEntityId = -1;
  SlaveEntityId  = -1;

  PtrMasterEntity = 0;
  PtrSlaveEntity  = 0;

  // Default (absolute) normal & tangential search tolerance
#ifdef SALINAS
  NormalSearchTol= 0.1e1;
  TangSearchTol  = 0.001e1;
#else
  NormalSearchTol= 0.1; 
  TangSearchTol  = 0.001; 
#endif
  // Use standard mortar basis as default
  MortarType = MortarHandler::STD;
  
  // Assume tied interaction as default
  InteractionType = MortarHandler::TIED;
  
  // Use general surface interaction type (i.e. non-conforming & non-matching) as default
  GeomType = MortarHandler::NON_MATCHING;

  // ACME FFI/NFI search output data
  nFFI                      = 0;
  nACMEFFIData              = 0;
  Slave_face_block_id       = 0;
  Slave_face_index_in_block = 0;
  Master_face_block_id      = 0;
  Master_face_index_in_block= 0;
  Slave_face_procs          = 0;
  Master_face_procs         = 0;
  ACMEFFI_index             = 0;
  ACMEFFI_data              = 0;

  ActiveSlaveNodeToElem = 0;
  ActiveMasterNodeToElem= 0;
  SlaveFaceToFFIConnect = 0;

  // mortar lmpc data
  nMortarLMPCs = 0;
  gIdFirstLMPC = 0;
  gIdLastLMPC  = 0;

  //ActiveSlaveFacesToMortarEl.clear();

  search_obj = 0;
  contact_obj = 0;
  mass = 0;
  dofmap = 0;
  node_global_ids = 0;

  TDEnfNumIter = 5;
  TDEnfConvTol = 1.0e-10;
  FrictionModel = 1; // TD_FRICTIONLESS
  for(int i=0; i<4; ++i) FrictionCoef[i] = 0.0;
  data = 0;
  share = 0;
  density = 0;
  wavespeed = 0;

  comm_proc_id = 0;
  number_nodes_to_partner = 0;
  comm_node = 0;

  NoSecondary = false;
  AveragedNodalNormals = false;
  DIST_ACME = 0;

  SelfContact = false;

  MortarScaling = 1.0;
  MortarIntegrationRule = 6;
  CtcMode = 1;

  ConstraintOptionsData = NULL;
}

/*
void 
MortarHandler::InitDefaultSearchTol()
{
  NormalSearchTol= 0.2;
  TangSearchTol  = 0.001;
}
*/

void
MortarHandler::DeleteFFIData()
{
  nFFI         = 0;
  nACMEFFIData = 0;
  if(Slave_face_block_id)        { delete [] Slave_face_block_id       ; Slave_face_block_id       = 0; }
  if(Slave_face_index_in_block)  { delete [] Slave_face_index_in_block ; Slave_face_index_in_block = 0; }
  if(Master_face_block_id)       { delete [] Master_face_block_id      ; Master_face_block_id      = 0; }
  if(Master_face_index_in_block) { delete [] Master_face_index_in_block; Master_face_index_in_block= 0; }
  if(Slave_face_procs)           { delete [] Slave_face_procs          ; Slave_face_procs          = 0; }
  if(Master_face_procs)          { delete [] Master_face_procs         ; Master_face_procs         = 0; }
  if(ACMEFFI_index)              { delete [] ACMEFFI_index             ; ACMEFFI_index             = 0; }
  if(ACMEFFI_data)               { delete [] ACMEFFI_data              ; ACMEFFI_data              = 0; }

  for(size_t i=0; i<MortarEls.size(); ++i) 
    if(MortarEls[i]) delete MortarEls[i];
  MortarEls.clear();

  NodalMortars.clear();

  for(std::vector<FFIPolygon*>::iterator it = CtcPolygons.begin(); it != CtcPolygons.end(); ++it) {
    delete *it;
  }
  CtcPolygons.clear();

  if(ActiveSlaveNodeToElem) { delete ActiveSlaveNodeToElem ; ActiveSlaveNodeToElem  = 0; }
  if(ActiveMasterNodeToElem){ delete ActiveMasterNodeToElem; ActiveMasterNodeToElem = 0; }
  if(SlaveFaceToFFIConnect) { delete SlaveFaceToFFIConnect ; SlaveFaceToFFIConnect  = 0; }

  ActiveSlaveNodes.clear();
  ActiveSlaveFacesToMortarEl.clear();

  // this only delete the array of pointers to the (active) master/slave elements
  // but not the element themselves (i.e. they are currently owned by the SurfaceEntity objects) 
  ActiveSlaveElemSet.deleteElems();
  ActiveMasterElemSet.deleteElems();
}

// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS
// -----------------------------------------------------------------------------------------------------
void 
MortarHandler::SetId(int _Id) { Id = _Id; }

void 
MortarHandler::SetMasterEntityId(int _MasterEntityId) { MasterEntityId = _MasterEntityId; }

void 
MortarHandler::SetSlaveEntityId(int  _SlaveEntityId ) { SlaveEntityId  = _SlaveEntityId ; }

void 
MortarHandler::SetSurfaceEntityId(int _MasterEntityId, int _SlaveEntityId) 
{ 
  SetMasterEntityId(_MasterEntityId);
  SetSlaveEntityId( _SlaveEntityId ); 
}

void
MortarHandler::SetPtrMasterEntity(SurfaceEntity* _PtrMasterEntity)
{
  PtrMasterEntity = _PtrMasterEntity;
}

void
MortarHandler::SetPtrSlaveEntity(SurfaceEntity* _PtrSlaveEntity)
{
  PtrSlaveEntity = _PtrSlaveEntity;
}

void
MortarHandler::SetPtrSurfaceEntity(SurfaceEntity* _PtrMasterEntity, SurfaceEntity* _PtrSlaveEntity)
{
  SetPtrMasterEntity( _PtrMasterEntity);
  SetPtrSlaveEntity(  _PtrSlaveEntity );
}

void 
MortarHandler::SetNormalSearchTol(double _NormalTol)
{
  NormalSearchTol = _NormalTol;
}
void 
MortarHandler::SetTangSearchTol(double _TangTol)
{
  TangSearchTol = _TangTol;
}

void 
MortarHandler::SetNormalAndTangSearchTol(double _NormalTol, double _TangTol)
{
  SetNormalSearchTol(_NormalTol);
  SetTangSearchTol(_TangTol);
}

void 
MortarHandler::SetSearchTol(double _NormalTol, double _TangTol)
{
  SetNormalAndTangSearchTol(_NormalTol, _TangTol);       
}

void
MortarHandler::SetMortarType(int _MortarType)
{
  //MortarType = MortarHandler::Mortar_Type(_MortarType);
  switch(_MortarType) {
    case MortarHandler::STD:
      MortarType = MortarHandler::STD;
      break;
    case MortarHandler::DUAL:
      MortarType = MortarHandler::DUAL;
      break;
    default:
      std::cerr << " In MortarHandler::SetMortarType(...): selected mortar type is NOT supported" << std::endl;
      std::cerr << " In MortarHandler::SetMortarType(...): use standard mortar instead" << std::endl;
      MortarType = MortarHandler::STD;
      break; 
  }
}

void
MortarHandler::SetMortarType(MortarHandler::Mortar_Type _MortarType)
{
  MortarType = _MortarType;
}

void
MortarHandler::SetInteractionType(int _InteractionType)
{
  //InteractionType = MortarHandler::Interaction_Type(_InteractionType);
  switch(_InteractionType) {
    case MortarHandler::TIED:
      InteractionType = MortarHandler::TIED;
      break;
    case MortarHandler::FSI:
      InteractionType = MortarHandler::FSI;
      break;
    case MortarHandler::CTC:
      InteractionType = MortarHandler::CTC;
      break;
    default:
      std::cerr << " In MortarHandler::SetInteractionType(...): selected interaction type is NOT supported" << std::endl;
      std::cerr << " In MortarHandler::SetInteractionType(...): use tied interaction instead" << std::endl;
      InteractionType = MortarHandler::TIED;
      break;
  }
}
                                                                                                                                       
void
MortarHandler::SetInteractionType(MortarHandler::Interaction_Type _InteractionType)
{
  InteractionType = _InteractionType;
}

void
MortarHandler::SetGeomType(int _GeomType)
{
  switch(_GeomType) {
    case MortarHandler::NON_MATCHING:
      GeomType = MortarHandler::NON_MATCHING;
      break;
    case MortarHandler::EQUIVALENCED:
      GeomType = MortarHandler::EQUIVALENCED;
      break;
    default:
      std::cerr << " In MortarHandler::SetGeomType(...): selected geometry type is NOT supported" << std::endl;      std::cerr << " In MortarHandler::SetGeomType(...): use general geometry type instead" << std::endl;
      GeomType = MortarHandler::NON_MATCHING;
      break;
  }
}
                                                                                                               
void
MortarHandler::SetGeomType(MortarHandler::Geom_Type _GeomType)
{
  GeomType = _GeomType;
}

void
MortarHandler::SetTDEnfParams(int _TDEnfNumIter, double _TDEnfConvTol)
{
  TDEnfNumIter = _TDEnfNumIter;
  TDEnfConvTol = _TDEnfConvTol;
}

void
MortarHandler::SetFrictionCoef(double _FrictionCoef)
{
  if(_FrictionCoef != 0.0) { 
    FrictionModel = 2; // TD_CONSTANT_FRICTION
    FrictionCoef[0] = _FrictionCoef;
  }
}

void
MortarHandler::SetFrictionCoef(double _FrictionCoef, double _ReferencePressure, double _OffsetPressure, double _PressureExponent)
{
  if(_FrictionCoef != 0.0 || _ReferencePressure != 0.0 || _OffsetPressure != 0.0 || _PressureExponent != 0.0) {
    FrictionModel = 5; // TD_PRESSURE_DEPENDENT
    FrictionCoef[0] = _FrictionCoef;
    FrictionCoef[1] = _ReferencePressure;
    FrictionCoef[2] = _OffsetPressure;
    FrictionCoef[3] = _PressureExponent;
  }
}

void
MortarHandler::SetFrictionCoef(double _StaticCoef, double _DynamicCoef, double _VelocityDecay)
{
  if(_StaticCoef != 0.0 || _DynamicCoef != 0.0 || _VelocityDecay != 0.0) {
    FrictionModel = 6; // TD_VELOCITY_DEPENDENT
    FrictionCoef[0] = _StaticCoef;
    FrictionCoef[1] = _DynamicCoef;
    FrictionCoef[2] = _VelocityDecay;
  }
}

void 
MortarHandler::SetNoSecondary(bool _NoSecondary)
{
  NoSecondary = _NoSecondary;
}

void
MortarHandler::SetAveragedNodalNormals(bool _AveragedNodalNormals)
{
  AveragedNodalNormals = _AveragedNodalNormals;
}

void
MortarHandler::SetSelfContact(bool _SelfContact)
{
  SelfContact = _SelfContact;
}

void
MortarHandler::SetDistAcme(int _DistAcme)
{
  DIST_ACME = _DistAcme;
}

void
MortarHandler::SetMortarScaling(double _MortarScaling)
{
  MortarScaling = _MortarScaling;
}

void
MortarHandler::SetMortarIntegrationRule(int _MortarIntegrationRule)
{
  MortarIntegrationRule = _MortarIntegrationRule;
}

void
MortarHandler::SetCtcMode(int _CtcMode)
{
  CtcMode = _CtcMode;
}

void
MortarHandler::SetConstraintOptions(ConstraintOptions& _ConstraintOptionsData)
{
  ConstraintOptionsData = new ConstraintOptions(_ConstraintOptionsData);
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int 
MortarHandler::GetId() { return Id; }

int 
MortarHandler::ID() { return Id; }

int 
MortarHandler::GetMasterEntityId() { return MasterEntityId; }

int 
MortarHandler::GetSlaveEntityId() { return SlaveEntityId; }

int
MortarHandler::GetnFFI() { return nFFI; }

SurfaceEntity*
MortarHandler::GetPtrMasterEntity() { return PtrMasterEntity; }

SurfaceEntity*
MortarHandler::GetPtrSlaveEntity() { return PtrSlaveEntity; }

FaceElemSet*
MortarHandler::GetPtrMasterFaceElemSet()
{
   return PtrMasterEntity->GetPtrFaceElemSet();
}
 
FaceElemSet*
MortarHandler::GetPtrSlaveFaceElemSet()
{
   return PtrSlaveEntity->GetPtrFaceElemSet();
}

int
MortarHandler::GetnMortarLMPCs() { return nMortarLMPCs; }

int
MortarHandler::GetIdFirstMortarLMPC() { return gIdFirstLMPC; } 

int
MortarHandler::GetIdLastMortarLMPC() { return gIdLastLMPC; } 

int
MortarHandler::GetCtcMode() { return CtcMode; }

double
MortarHandler::GetNormalTol() { return NormalSearchTol; }

double
MortarHandler::GetTangentialTol() { return TangSearchTol; }

MortarHandler::Mortar_Type
MortarHandler::GetMortarType() { return(MortarType); }

MortarHandler::Interaction_Type
MortarHandler::GetInteractionType() { return(InteractionType); }

MortarHandler::Geom_Type
MortarHandler::GetGeomType() { return(GeomType); }

ConstraintOptions*
MortarHandler::GetConstraintOptions() { return (ConstraintOptionsData); }

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS
// -----------------------------------------------------------------------------------------------------
void
MortarHandler::Print()
{
  std::cout << " -> Mortar condition Id : " << Id << std::endl;
  std::cout << "   * master entity Id     : " << MasterEntityId  << std::endl;
  std::cout << "   * slave entity Id      : " << SlaveEntityId   << std::endl;
  std::cout << "   * normal search tol    : " << NormalSearchTol << std::endl;
  std::cout << "   * tangential search tol: " << TangSearchTol   << std::endl;
  std::cout << "   * type of mortar space : ";
  switch(MortarType){
   case MortarHandler::STD:
     std::cout << " standard" << std::endl;
     break;
   case MortarHandler::DUAL: 
     std::cout << " dual" << std::endl;
     break;
  }
  std::cout << "   * type of interaction  : ";
  switch(InteractionType) {
   case MortarHandler::TIED:
     std::cout << " TIED" << std::endl;
     break;
   case MortarHandler::FSI: 
     std::cout << " FSI" << std::endl;
     break;
   case MortarHandler::CTC:
     std::cout << " contact" << std::endl;
     break;
  }
  std::cout << "   * type of geometry   : ";
  switch(GeomType){
   case MortarHandler::NON_MATCHING:
     std::cout << " general" << std::endl;
     break;
   case MortarHandler::EQUIVALENCED: 
     std::cout << " equivalenced" << std::endl;
     break;
  }
}

#ifdef HB_ACME_FFI_DEBUG
/* TODO these should be called from MortarHandler::CreateFFIPolygon and passed ctcPolygons as an argument
void
MortarHandler::PrintFFIPolyVertices(FILE* file, int& firstVertId)
{
  CoordSet& cs = PtrSlaveEntity->GetNodeSet();
  for(int i=0; i<CtcPolygons.size(); i++) {
    CtcPolygons[i].PrintSlaveVertices(file, cs, firstVertId);
  }

  CoordSet& cs2 = PtrMasterEntity->GetNodeSet();
  for(int i=0; i<CtcPolygons.size(); i++)
    CtcPolygons[i].PrintMasterVertices(file, cs2, firstVertId);
}

void
MortarHandler::PrintFFIPolyTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode)
{
  for(int i=0; i<CtcPolygons.size(); i++)
    CtcPolygons[i].PrintFFIPolygonTopo(file, EdgeOffset, VertOffset, elCode);
}
*/
#endif

// -----------------------------------------------------------------------------------------------------
//                                            PRIVATE METHODS
// -----------------------------------------------------------------------------------------------------
void
MortarHandler::ComputeOneFFIMandN(int iFFI, CoordSet &cs, std::vector<FFIPolygon*> &CtcPolygons)
{
  Connectivity* SlaveACMEBlocksMap = PtrSlaveEntity->GetPtrACMEBlocksMap();
  int SlaveBlkId = Slave_face_block_id[iFFI]-1;
  int slave_face = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
  int iMortar = ActiveSlaveFacesToMortarEl[slave_face];
  MortarElement* MortarEl = MortarEls[iMortar];
  if(InteractionType==MortarHandler::FSI) {
    CtcPolygons[iFFI]->ComputeNormalN(MortarEl, cs, MortarIntegrationRule);
  }
  else if(InteractionType==MortarHandler::CTC) {
    CoordSet& csm = PtrMasterEntity->GetNodeSet();
    double offset = PtrMasterEntity->GetShellThickness()/2 + PtrSlaveEntity->GetShellThickness()/2;
    CtcPolygons[iFFI]->ComputeNormalGeoGap(MortarEl, cs, csm, MortarIntegrationRule, offset);
    CtcPolygons[iFFI]->ComputeGradNormalGeoGap(MortarEl, cs, csm, MortarIntegrationRule, offset);
  }
  else { // TIED
    CtcPolygons[iFFI]->ComputeM(MortarEl, cs, MortarIntegrationRule);
    CtcPolygons[iFFI]->ComputeN(MortarEl, cs, MortarIntegrationRule);
  }
}

void
MortarHandler::MakeOneNodalMortarLMPC(int i, std::vector<FFIPolygon*> &CtcPolygons, bool Dual)
{
  NodalMortars[i].SetRefData(ActiveSlaveNodes[i], MortarScaling);
  NodalMortars[i].MakeSlaveLink(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, Dual);
  NodalMortars[i].MakeMasterLink(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
  if(InteractionType==MortarHandler::FSI)
    NodalMortars[i].BuildWetFSICoupling(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
  else if(InteractionType==MortarHandler::CTC)
    NodalMortars[i].BuildMortarCtcLMPC(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
  else // TIED
    NodalMortars[i].BuildMortarLMPC(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
}

#ifdef USE_ACME 
void
MortarHandler::PerformACMEFFISearch()
{
#ifdef MORTAR_DEBUG
  filePrint(stderr," ### GET IN PerformACMEFFISearch() ###\n");
#endif
  // NEW VERSION
  build_search();
  set_search_data(4); // interaction_type = 4 (FaceFace)
  set_node_configuration(1); // config_type = 1 (current configuration)
  set_search_options();
  perform_search(1); // search_algorithm = 1 (static 1-configuration search)
  get_interactions(4); // interaction_type = 4 (FaceFace)
  //set_node_constraints(domain->nDirichlet(), domain->getDBC()); // YYYY
}
#else
// !!! DUMMY METHOD: SHOULD NOT BE CALLED !!!
void
MortarHandler::PerformACMEFFISearch()
{
  filePrint(stderr," !!! FATAL PB in MortarHandler::PerformACMEFFISearch: USE_ACME FLAG WAS NOT DEFINED !!!\n");
  filePrint(stderr," => IF YOU GET TO THIS POINT, THIS MEANS YOU HAVE DEFINED SOME MORTAR CONDITION(S)\n");
  filePrint(stderr," => BUT IMPOSSIBLE TO DO MORTAR WITHOUT THE ACME LIB !!!\n");
  filePrint(stderr," => STOP EXECUTION \n");
  exit(-1);
}
#endif

void
MortarHandler::CreateFFIPolygon()
{
  // clear data from previous call to this function
  MortarEls.clear();
  NodalMortars.clear();
  ActiveSlaveNodes.clear();
  ActiveSlaveElemSet.deleteElems();
  ActiveMasterElemSet.deleteElems();
  if(ActiveSlaveNodeToElem) { delete ActiveSlaveNodeToElem; ActiveSlaveNodeToElem = 0; }
  if(ActiveMasterNodeToElem) { delete ActiveMasterNodeToElem; ActiveMasterNodeToElem = 0; }
  if(SlaveFaceToFFIConnect) { delete SlaveFaceToFFIConnect; SlaveFaceToFFIConnect = 0; }
  ActiveSlaveFacesToMortarEl.clear();

  // BRUTE FORCE TEST
  CoordSet& cs = PtrSlaveEntity->GetNodeSet();

#ifdef HB_MORTAR_TIMER
  double t0,t00,dt0,dt1,dt2,dt3,dt4;
  t00 = getTime(); t0 = t00;
#endif  
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * Create the FFI polygon \n");
#endif  
  // Allocate CtcPolygons array
  CtcPolygons.resize(nFFI);

  FaceElemSet* MasterElemSet = PtrGlobalMasterEntity->GetPtrFaceElemSet();
  FaceElemSet* SlaveElemSet  = PtrGlobalSlaveEntity->GetPtrFaceElemSet();
  Connectivity* SlaveACMEBlocksMap = PtrGlobalSlaveEntity->GetPtrACMEBlocksMap();
  Connectivity* MasterACMEBlocksMap= PtrGlobalMasterEntity->GetPtrACMEBlocksMap();

  // Create each FFIPolygon
  // !!! ACME indices in Fortran indexing !!!
  std::vector<int> IndexActiveSlaveFaces ; IndexActiveSlaveFaces.reserve(512);
  std::vector<int> IndexActiveMasterFaces; IndexActiveMasterFaces.reserve(512);
  int SlaveBlkId, MasterBlkId;
  for(int iFFI=0; iFFI < nFFI; iFFI++) {
    SlaveBlkId  = Slave_face_block_id[iFFI]-1;
    if(Slave_face_index_in_block[iFFI] <= 0) std::cerr << "error here in MortarHandler::CreateFFIPolygon, iFFI = " << iFFI
                                                  << ", Slave_face_index_in_block[iFFI] = " << Slave_face_index_in_block[iFFI] << std::endl;
    int slave_face  = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
    FaceElement *SlaveFaceEl  = (*SlaveElemSet)[slave_face];

    if(SelfContact) {
      MasterBlkId = Master_face_block_id[iFFI]-1;
      MasterACMEBlocksMap = SlaveACMEBlocksMap;
      MasterElemSet = SlaveElemSet;
      PtrMasterEntity = PtrSlaveEntity;
    }
    else {
      MasterBlkId = Master_face_block_id[iFFI]-1 - 4; // in all the ACME blocks, the master blocks are last
    }
    if(Master_face_index_in_block[iFFI] <= 0) std::cerr << "error here in MortarHandler::CreateFFIPolygon, iFFI = " << iFFI 
                                                        << ", Master_face_index_in_block[iFFI] = " << Master_face_index_in_block[iFFI] << std::endl;
    int master_face = (*MasterACMEBlocksMap)[MasterBlkId][Master_face_index_in_block[iFFI]-1];
    FaceElement *MasterFaceEl = (*MasterElemSet)[master_face];

    int FFIDataOffset = ACMEFFI_index[iFFI];
    int nVertices = (int) ACMEFFI_data[FFIDataOffset];
    double* ACME_FFI_Data = &ACMEFFI_data[FFIDataOffset];
#ifdef USE_FFI_DERIVATIVES
    int ACME_FFI_Derivatives_Order = (InteractionType == MortarHandler::CTC) ? 1 : 0;
#else
    int ACME_FFI_Derivatives_Order = 0;
#endif

    CtcPolygons[iFFI] = new FFIPolygon(MasterFaceEl, SlaveFaceEl, nVertices, ACME_FFI_Data, ACME_FFI_Derivatives_Order);
    
#ifdef FFI_MORTAR_DEBUG
    filePrint(stderr, " * -------------------------------- \n");
    filePrint(stderr, "   -> create FFIPolygon %d\n", iFFI+1);
    std::cerr << "   -> MasterBlkId = " << MasterBlkId << " master_face = " << master_face << std::endl;
    std::cerr << "   -> SlaveBlkId  = " << SlaveBlkId << " slave_face = " << slave_face << std::endl;
    filePrint(stderr, "   -> FFIDataOffset = %d\n", FFIDataOffset);
    std::cerr << "    -> slave  face element:\n"; SlaveFaceEl->print();
    std::cerr << "    -> master face element:\n"; MasterFaceEl->print();
    CtcPolygons[iFFI]->Print();
#endif

    // Get indices of active slave/master face elements 
    IndexActiveSlaveFaces.push_back(slave_face);
    IndexActiveMasterFaces.push_back(master_face);
  }
#ifdef HB_MORTAR_TIMER
  dt0 = (getTime()-t0)/1000.;
#endif
  // Eliminates duplicate face indices
  std::sort(IndexActiveSlaveFaces.begin() , IndexActiveSlaveFaces.end());
  std::sort(IndexActiveMasterFaces.begin(), IndexActiveMasterFaces.end());

  IndexActiveSlaveFaces.erase(unique(IndexActiveSlaveFaces.begin(),IndexActiveSlaveFaces.end()), 
                               IndexActiveSlaveFaces.end());

  IndexActiveMasterFaces.erase(unique(IndexActiveMasterFaces.begin(),IndexActiveMasterFaces.end()),
                               IndexActiveMasterFaces.end());

  int nActiveSlaveFaces = IndexActiveSlaveFaces.size();
  
  // Create active master & slave face element set
  // & associated connectivity object
#ifdef MORTAR_DEBUG
  int nActiveMasterFaces= IndexActiveMasterFaces.size();
  filePrint(stderr,"   * number of active slave face elem. = %6d\n", nActiveSlaveFaces);
  filePrint(stderr,"   * number of active master face elem.= %6d\n", nActiveMasterFaces);
  filePrint(stderr,"   * active slave & master face element sets\n"); 
#endif
#ifdef HB_MORTAR_TIMER
  t0 = getTime();
#endif
  std::map<int,int> IndexActiveSlaveFacesToActiveSlaveElemSetMap; 
  for(int i=0; i<int(IndexActiveSlaveFaces.size()); ++i) {
    ActiveSlaveElemSet.elemadd(i, (*SlaveElemSet)[IndexActiveSlaveFaces[i]]); 
    IndexActiveSlaveFacesToActiveSlaveElemSetMap[IndexActiveSlaveFaces[i]] = i;
  }  
  Connectivity ActiveSlaveElemToNode(SetAccess<FaceElemSet>{ActiveSlaveElemSet});
  ActiveSlaveNodeToElem = ActiveSlaveElemToNode.alloc_reverse();

  for(int i=0; i<int(IndexActiveMasterFaces.size()); ++i) {
    ActiveMasterElemSet.elemadd(i, (*MasterElemSet)[IndexActiveMasterFaces[i]]); 
  }
  Connectivity ActiveMasterElemToNode(SetAccess<FaceElemSet>{ActiveMasterElemSet});
  ActiveMasterNodeToElem = ActiveMasterElemToNode.alloc_reverse();

  ActiveSlaveNodes.clear();
  ActiveSlaveNodes.reserve(128);
  for(int iel=0; iel<int(IndexActiveSlaveFaces.size()); iel++) {
    size_t old_size(ActiveSlaveNodes.size());
    ActiveSlaveNodes.insert(ActiveSlaveNodes.end(), ActiveSlaveElemSet[iel]->nNodes(), 0);
    ActiveSlaveElemSet[iel]->GetNodes(&ActiveSlaveNodes[old_size]);
  }
  std::sort(ActiveSlaveNodes.begin(),ActiveSlaveNodes.end());
  ActiveSlaveNodes.erase(std::unique(ActiveSlaveNodes.begin(),ActiveSlaveNodes.end()),
                         ActiveSlaveNodes.end());

  int nActiveSlaveNodes = ActiveSlaveNodes.size();
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * nActiveSlaveNodes = %d\n", nActiveSlaveNodes);
#endif

  // Create the map ActiveSlaveFaces to FFIs
  // 1) count the number of FFIs per ActiveSlaveFaces
  int* SlaveFaceToFFIpointer = new int[nActiveSlaveFaces+1];
  int* SlaveFaceToFFItarget = new int[nFFI];
  { 
    std::vector<int> nFFIperSlaveFace(nActiveSlaveFaces, 0);
    for(int iFFI=0; iFFI<nFFI; iFFI++) {
      int SlaveBlkId = Slave_face_block_id[iFFI]-1;
      int slave_face = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
      nFFIperSlaveFace[IndexActiveSlaveFacesToActiveSlaveElemSetMap[slave_face]]++;  
    }
    // 2) set the map pointer array
    SlaveFaceToFFIpointer[0] = 0;
    for(int i=0; i<nActiveSlaveFaces; i++)
      SlaveFaceToFFIpointer[i+1] = SlaveFaceToFFIpointer[i]+nFFIperSlaveFace[i];

    // 3) fill the target (map) array
    std::vector<int> Offset(nActiveSlaveFaces, 0);
    for(int iFFI=0; iFFI<nFFI; iFFI++){
      int SlaveBlkId = Slave_face_block_id[iFFI]-1;
      int slave_face = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
      int i = IndexActiveSlaveFacesToActiveSlaveElemSetMap[slave_face];
      SlaveFaceToFFItarget[SlaveFaceToFFIpointer[i]+Offset[i]] = iFFI; 
      Offset[i]++;
    }
  }
  SlaveFaceToFFIConnect = new Connectivity(nActiveSlaveFaces, SlaveFaceToFFIpointer, SlaveFaceToFFItarget);

#ifdef HB_MORTAR_TIMER
  dt1 = (getTime()-t0)/1000.;
  t0 = getTime();
#endif
  // Create Mortar element
  bool Dual = false;
  if(MortarType==MortarHandler::DUAL) { Dual = true; }
#ifdef MORTAR_DEBUG
  if(Dual) filePrint(stderr,"   * Create Dual Mortar elements \n");
  else     filePrint(stderr,"   * Create Std Mortar elements \n"); 
#endif
  MortarEls.assign(nActiveSlaveFaces, (MortarElement*) 0);

  for(int i=0; i<nActiveSlaveFaces; ++i) {
    FaceElement* SlaveFaceEl  = (*SlaveElemSet)[IndexActiveSlaveFaces[i]];
    MortarEls[i] = CreateMortarElement(SlaveFaceEl, cs, Dual); 
    ActiveSlaveFacesToMortarEl[IndexActiveSlaveFaces[i]] = i;
  }

#ifdef HB_MORTAR_TIMER
  dt2 = (getTime()-t0)/1000.;
  t0 = getTime();
#endif
  // Integrate shape fcts product
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * Compute FFI contributions to the shape fct products\n"); 
#endif
#ifdef HB_THREAD_FFI_M_N
  execParal(nFFI, this, &MortarHandler::ComputeOneFFIMandN, cs, CtcPolygons);
#else
  for(int i=0; i<nFFI; ++i) ComputeOneFFIMandN(i, cs, CtcPolygons);
#endif
#ifdef HB_MORTAR_TIMER
  dt3 = (getTime()-t0)/1000.;
  t0 = getTime();
#endif
  // Create nodal Mortar shape fcts (NO BOUNDARY MODIFICATION) 
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * Create nodal Mortar shape fcts\n"); 
  filePrint(stderr,"     -> NO BOUNDARY MODIFICATION\n"); 
  filePrint(stderr,"     -> NodalMortars.size() = %d\n",NodalMortars.size()); 
#endif
  NodalMortars.assign(nActiveSlaveNodes, NodalMortarShapeFct());
#ifdef HB_THREAD_NODALMORTAR
  execParal(nActiveSlaveNodes, this, &MortarHandler::MakeOneNodalMortarLMPC, CtcPolygons, Dual);
#else
  for(int i=0; i<nActiveSlaveNodes; ++i) MakeOneNodalMortarLMPC(i, CtcPolygons, Dual);
#endif

#ifdef HB_MORTAR_TIMER
  dt4 = (getTime()-t0)/1000.;
  double dttot = dt0+dt1+dt2+dt3+dt4;
  filePrint(stderr,"   * CPU time statistic of MortarHandler::CreateFFIPolygon(...):\n");
  filePrint(stderr,"    # making the FFI contact polygons (seq): %e s (%2.4f %%)\n",dt0,100.*dt0/dttot);
  filePrint(stderr,"    # making mapping & connectivities (seq): %e s (%2.4f %%)\n",dt1,100.*dt1/dttot);
  filePrint(stderr,"    # creating the mortar elements    (seq): %e s (%2.4f %%)\n",dt2,100.*dt2/dttot);
#ifdef HB_THREAD_FFI_M_N
  filePrint(stderr,"    # computing FFI M & N contribution(// ): %e s (%2.4f %%)\n",dt3,100.*dt3/dttot);
#else
  filePrint(stderr,"    # computing FFI M & N contribution(seq): %e s (%2.4f %%)\n",dt3,100.*dt3/dttot);
#endif
#ifdef HB_THREAD_NODALMORTAR
  filePrint(stderr,"    # assembling nodal mortar LMPCs   (// ): %e s (%2.4f %%)\n",dt4,100.*dt4/dttot);
#else
  filePrint(stderr,"    # assembling nodal mortar LMPCs   (seq): %e s (%2.4f %%)\n",dt4,100.*dt4/dttot);
#endif
  filePrint(stderr,"    # total CPU time                       : %e s\n",dttot);
#endif
}

// This method to manually create the ACMEFFIData arrays without performing the ACME FFI search
// This is intended for equivalenced geometry & so bypass the ACME FFI search 
// The idea of creating the ACMEFFIData arrays is to be able to reuse the method
// MortarHandler::CreateFFIPolygon() even for equivalenced geometry
void
MortarHandler::CreateACMEFFIData()
{
  nFFI         = PtrSlaveEntity->nFaceElements();
  filePrint(stderr," -> nFFI = %d\n",nFFI);
  if(nFFI==0) return;
  FaceElemSet* FaceElSet  = GetPtrSlaveFaceElemSet();
  int FFI_data_size= 0;
  for(int iFFI=0; iFFI<nFFI; iFFI++) { // one FFI = one face el.
    FFI_data_size += 7*(*FaceElSet)[iFFI]->nVertices(); // Underlying linear geometry 
  }

  Slave_face_block_id       = new int[nFFI];
  Slave_face_index_in_block = new int[nFFI];
  Master_face_block_id      = new int[nFFI];
  Master_face_index_in_block= new int[nFFI];
  ACMEFFI_index             = new int[nFFI];
  ACMEFFI_data              = new double[FFI_data_size];

  int T3index,T6index,Q4index,Q8index;
  T3index = T6index = Q4index = Q8index = 1;
  int offset = 0;
  for(int iFFI=0; iFFI<nFFI; iFFI++){ // one FFI = one face el.
    ACMEFFI_index[iFFI]             = offset;
    //Slave_face_index_in_block[iFFI] = iFFI;
    //Master_face_block_id[iFFI]      = 0;
    //Master_face_index_in_block[iFFI]= iFFI;
    int N = (*FaceElSet)[iFFI]->nVertices(); // Underlying linear geometry
    ACMEFFI_data[offset++] = N; // number of edges/vertices of the FFI polygon
    for(int i=0; i<N; i++){  
       ACMEFFI_data[offset+i]  = i+1; // master edge flag (Fortran numbering)
       ACMEFFI_data[offset+i+N]= i+1; // slave edge flag (Fortran numbering)
    }
    offset += 2*N;
    int etype = (*FaceElSet)[iFFI]->GetFaceElemType();
    switch(etype)
    {
     case FaceElement::QUADFACEL4:
       Slave_face_block_id[iFFI] = 1;
       Master_face_block_id[iFFI]= 5;
       Slave_face_index_in_block[iFFI] = Q4index;
       Master_face_index_in_block[iFFI]= Q4index;
       Q4index++;
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] = -1.0;
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] = -1.0;
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] = -1.0;
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] = -1.0;
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] =  1.0;
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] =  1.0;
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] =  1.0;
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] =  1.0;
       break;
   
     case FaceElement::QUADFACEQ8:
     case FaceElement::QUADFACEQ9:
     case FaceElement::QUADFACEC12:
       Slave_face_block_id[iFFI] = 2;
       Master_face_block_id[iFFI]= 6;
       Slave_face_index_in_block[iFFI] = Q8index;
       Master_face_index_in_block[iFFI]= Q8index;
       Q8index++;
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] = -1.0; 
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] = -1.0; 
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] = -1.0; 
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] = -1.0; 
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] =  1.0; 
       ACMEFFI_data[offset++] =  1.0; ACMEFFI_data[offset++] =  1.0; 
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] =  1.0; 
       ACMEFFI_data[offset++] = -1.0; ACMEFFI_data[offset++] =  1.0; 
       break;
     case FaceElement::TRIFACEL3:
       Slave_face_block_id[iFFI] = 3;
       Master_face_block_id[iFFI]= 7;
       Slave_face_index_in_block[iFFI] = T3index;
       Master_face_index_in_block[iFFI]= T3index;
       T3index++;
       ACMEFFI_data[offset++] = 1.0; ACMEFFI_data[offset++] = 0.0;
       ACMEFFI_data[offset++] = 1.0; ACMEFFI_data[offset++] = 0.0;
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 1.0;
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 1.0;
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 0.0;
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 0.0;
       break;

     case FaceElement::TRIFACEQ6:
     case FaceElement::TRIFACEC10:
       Slave_face_block_id[iFFI] = 4;
       Master_face_block_id[iFFI]= 8;
       Slave_face_index_in_block[iFFI] = T6index;
       Master_face_index_in_block[iFFI]= T6index;
       T6index++;
       ACMEFFI_data[offset++] = 1.0; ACMEFFI_data[offset++] = 0.0; 
       ACMEFFI_data[offset++] = 1.0; ACMEFFI_data[offset++] = 0.0; 
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 1.0; 
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 1.0; 
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 0.0; 
       ACMEFFI_data[offset++] = 0.0; ACMEFFI_data[offset++] = 0.0; 
       break;
     default:
       std::cerr << "Face element Type " << etype << " is NOT supported/implemented." << std::endl;
       exit(-1);
       return;
    }
  }

  PtrGlobalMasterEntity = PtrMasterEntity;
  PtrGlobalSlaveEntity = PtrSlaveEntity;
}

// 09/09/03: Change the order in which the Mortar LMPCs are added to the
// standard LMPC array (lmpc):
// old: [(x,y,z),...,(x,y,z)]
// new: [(x,...,x),(y,...,y),(z,...,z)]
// see order flag "Order"
void
MortarHandler::AddMortarLMPCs(ResizeArray<LMPCons*>* LMPCArray, int& numLMPC, int &numCTC, int nDofs, int* Dofs)
{
 int* SlaveLlToGlNodeMap = PtrSlaveEntity->GetPtrLlToGlNodeMap();
 int* MasterLlToGlNodeMap= PtrMasterEntity->GetPtrLlToGlNodeMap();

  // !!! USE NEGATIVE NUMBERING FOR MY MORTAR LMPC !!!
#ifdef MORTAR_DEBUG
  std::cerr <<"In MortarHandler::AddMortarLMPCs: numLMPC at input: "<<numLMPC<<std::endl;
#endif
  int Order = 1; // order flag 
  int lmpcnum = 0; // !! NEED TO BE CHANGED FOR THE CASE OF SEVERAL MORTAR CONDITIONS
                   // OK FOR NOW BECAUSE WE DO NOT USE THIS ID NUMBER !!

  double rhs  = 0.0;  

  // set the global ID num. of the first Mortar LMPC
  gIdFirstLMPC = numLMPC;
  nMortarLMPCs = 0;
  if(InteractionType == MortarHandler::CTC) {
    for(int i=0; i<int(NodalMortars.size()); i++) {
      lmpcnum--;
      LMPCons* MortarLMPC = NodalMortars[i].CreateMortarCtcLMPCons(lmpcnum, SlaveLlToGlNodeMap, MasterLlToGlNodeMap, CtcMode);
      if(MortarLMPC) { 
        MortarLMPC->id.first = Id; 
        (*LMPCArray)[numLMPC++] = MortarLMPC; 
        nMortarLMPCs++; 
        numCTC++; 
        if(ConstraintOptionsData) { 
          MortarLMPC->lagrangeMult = ConstraintOptionsData->lagrangeMult;
          MortarLMPC->penalty = ConstraintOptionsData->penalty;
        }
      }
    }
  }
  else {
    // default dofs (tie Ux,Uy,Uz)
    int threedofs[3]= {0,1,2};
    int ndofs = (nDofs) ? nDofs : 3;
    int*dofs  = (nDofs) ? Dofs  : threedofs ;
    if(Order==0){ // -> [(x,y,z),...,(x,y,z)] 
      for(int i=0; i<int(NodalMortars.size()); i++){
        for(int j=0; j<ndofs; j++){
          lmpcnum--; 
          LMPCons* MortarLMPC = NodalMortars[i].CreateMortarLMPCons(lmpcnum, dofs[j], rhs, 
                                                                    SlaveLlToGlNodeMap, MasterLlToGlNodeMap);
          if(MortarLMPC) {
            (*LMPCArray)[numLMPC++] = MortarLMPC;
            nMortarLMPCs++;
            if(ConstraintOptionsData) { 
              MortarLMPC->lagrangeMult = ConstraintOptionsData->lagrangeMult;
              MortarLMPC->penalty = ConstraintOptionsData->penalty;
            }
          }
        }  
      }
    }
    if(Order==1){ // -> [(x,...,x),(y,...,y),(z,...,z)]
      for(int j=0; j<ndofs; j++){
        for(int i=0; i<int(NodalMortars.size()); i++){
          lmpcnum--;
          LMPCons* MortarLMPC = NodalMortars[i].CreateMortarLMPCons(lmpcnum, dofs[j], rhs,
                                                                    SlaveLlToGlNodeMap, MasterLlToGlNodeMap);
          if(MortarLMPC) {
            (*LMPCArray)[numLMPC++] = MortarLMPC;
            nMortarLMPCs++;
            if(ConstraintOptionsData) {
              MortarLMPC->lagrangeMult = ConstraintOptionsData->lagrangeMult;
              MortarLMPC->penalty = ConstraintOptionsData->penalty;
            }
          }
        }
      }
    }
  }

#ifdef MORTAR_DEBUG
  std::cerr <<" -> numLMPC at output: "<<numLMPC<<std::endl;
  std::cerr <<" -> adds "<<nMortarLMPCs<<" LMPCs"<<std::endl;
#endif
  // set the global ID num. of the last Mortar LMPC
  gIdLastLMPC = numLMPC-1; 
}

void
MortarHandler::AddWetFSI(ResizeArray<LMPCons*>* FSIArray, int& numFSI)
{
  int* SlaveLlToGlNodeMap = PtrSlaveEntity->GetPtrLlToGlNodeMap();
  int* MasterLlToGlNodeMap= PtrMasterEntity->GetPtrLlToGlNodeMap();

  for(int i=0; i<int(NodalMortars.size()); i++) {
    LMPCons* WetFSI = NodalMortars[i].CreateWetFSICons(SlaveLlToGlNodeMap, MasterLlToGlNodeMap);
    if(WetFSI) { (*FSIArray)[numFSI++] = WetFSI; nMortarLMPCs++; }
  }
}

//***************************************************************************************************
//***************************************************************************************************
// new ACME interface functions
void
MortarHandler::build_search(bool tdenforceFlag, int numSub, SubDomain **sd)
{
#ifdef USE_ACME
#ifdef MORTAR_DEBUG
  filePrint(stderr," ### IN MortarHandler::build_search() ###\n");
#endif
  // Generate ACME input data
  int dimensionality  = 3;
  int num_states      = 1;
  num_entity_keys = 8; // always assume ONLY 4 face blocks per surface entity 
  int num_node_blocks = 1;

  int num_analytical_surfs     = 0; // required by ACME 2.2

  ContactSearch::ContactNode_Type *node_block_types;
  node_block_types    = new ContactSearch::ContactNode_Type[1];
  node_block_types[0] = ContactSearch::NODE;
  
  int num_nodes_per_block[1];

  int number_face_blocks = 8; // always assume ONLY 4 face blocks per surface entity 

  ContactSearch::ContactFace_Type *face_block_types;
  face_block_types   = new ContactSearch::ContactFace_Type[number_face_blocks];

  ContactSearch::ContactElement_Type* element_block_types = 0;
  int num_element_blocks         = 0;
  int* number_elements_per_block = 0;
  int* element_ids               = 0;
  int* element_connectivity      = 0;

  num_comm_partners         = 0;
  comm_proc_id              = 0;
  number_nodes_to_partner   = 0;
  comm_node                 = 0;
#ifdef DISTRIBUTED
  MPI_Comm mpi_Comm = (DIST_ACME == 0) ? MPI_COMM_SELF : *(structCom->getCommunicator());
  /* perhaps we should use only cpus which have local surface nodes
  Communicator* allCom[2];
  int color = (PtrSlaveEntity->GetnVertices() + PtrMasterEntity->GetnVertices() > 0) ? 0 : 1;
  structCom->split(color, 2, allCom);
  MPI_Comm &mpi_Comm = *(allCom[0]->getCommunicator()); */
#elif defined(USE_MPI)
  MPI_Comm mpi_Comm = MPI_COMM_SELF;
#else 
  int mpi_Comm = 0;
#endif

  // !! Assume all the face elements in a FaceElemSet are of 2 types: Tri3/Tri6 & Quad4/Quad8 !!
  FaceElemSet* MasterFaceElemSet = GetPtrMasterFaceElemSet(); 
  FaceElemSet* SlaveFaceElemSet  = GetPtrSlaveFaceElemSet(); 
  double *face_block_shell_thickness = new double[number_face_blocks];
  // -> !!! UNDERLAYING LINEAR FACE ELEMENT TYPE !!!
  if(PtrSlaveEntity->GetIsShellFace() && tdenforceFlag) {
    face_block_types[0] = ContactSearch::SHELLQUADFACEL4;
    face_block_shell_thickness[0] = PtrSlaveEntity->GetShellThickness();
    face_block_types[1] = ContactSearch::SHELLQUADFACEL4; // slave Quad8 face els
    face_block_shell_thickness[1] = PtrSlaveEntity->GetShellThickness();
  }
  else {
    face_block_types[0] = ContactSearch::QUADFACEL4; // slave Quad4 face els
    face_block_types[1] = ContactSearch::QUADFACEL4; // slave Quad8 face els
  }
  if(PtrSlaveEntity->GetIsShellFace() && tdenforceFlag) {
    face_block_types[2] = ContactSearch::SHELLTRIFACEL3;
    face_block_shell_thickness[2] = PtrSlaveEntity->GetShellThickness();
    face_block_types[3] = ContactSearch::SHELLTRIFACEL3;  // slave Tri6 face els
    face_block_shell_thickness[3] = PtrSlaveEntity->GetShellThickness();
  }
  else {
    face_block_types[2] = ContactSearch::TRIFACEL3;  // slave Tri3 face els
    face_block_types[3] = ContactSearch::TRIFACEL3;  // slave Tri6 face els
  }
  if(PtrMasterEntity->GetIsShellFace() && tdenforceFlag) {
    face_block_types[4] = ContactSearch::SHELLQUADFACEL4;
    face_block_shell_thickness[4] = PtrMasterEntity->GetShellThickness();
    face_block_types[5] = ContactSearch::SHELLQUADFACEL4; // master Quad8 face els
    face_block_shell_thickness[5] = PtrMasterEntity->GetShellThickness();
  }
  else {
    face_block_types[4] = ContactSearch::QUADFACEL4; // master Quad4 face els
    face_block_types[5] = ContactSearch::QUADFACEL4; // master Quad8 face els
  }
  if(PtrMasterEntity->GetIsShellFace() && tdenforceFlag) {
    face_block_types[6] = ContactSearch::SHELLTRIFACEL3;
    face_block_shell_thickness[6] = PtrMasterEntity->GetShellThickness();
    face_block_types[7] = ContactSearch::SHELLTRIFACEL3;  // master Tri6 face els
    face_block_shell_thickness[7] = PtrMasterEntity->GetShellThickness();
  }
  else {
    face_block_types[6] = ContactSearch::TRIFACEL3;  // master Tri3 face els
    face_block_types[7] = ContactSearch::TRIFACEL3;  // master Tri6 face els
  }

  std::map<int, int> nVerticesPerFaceType;
  nVerticesPerFaceType[ContactSearch::QUADFACEL4] = 4;
  nVerticesPerFaceType[ContactSearch::QUADFACEQ8] = 4;
  nVerticesPerFaceType[ContactSearch::TRIFACEL3 ] = 3;
  nVerticesPerFaceType[ContactSearch::TRIFACEQ6 ] = 3;
  nVerticesPerFaceType[ContactSearch::SHELLQUADFACEL4] = 4;
  nVerticesPerFaceType[ContactSearch::SHELLTRIFACEL3] = 3;  

#ifdef MORTAR_DEBUG
  std::cerr << " -> In MortarHandler::PerformACMEFFISearch(): " << std::endl;
  std::cerr << "   * slave face_block_types  (true)= " << face_block_types[0]<< std::endl;
  std::cerr << "   * master face_block_types (true)= " << face_block_types[1]<< std::endl;	
#endif

#ifndef SOWER_SURFS
  // Check for common nodes.
  if(InteractionType != MortarHandler::FSI) { // Bypass this check for wet FSI interface
    std::vector<int> CommonNodes; CommonNodes.reserve(128);
    // Warning: assume that the following data have been previousely SORTED
    int* MastergVertexIds = PtrMasterEntity->GetPtrGlVertexIds();
    int* SlavegVertexIds = PtrSlaveEntity->GetPtrGlVertexIds();
    set_intersection(SlavegVertexIds, SlavegVertexIds+PtrSlaveEntity->GetnVertices(),
                     MastergVertexIds, MastergVertexIds+PtrMasterEntity->GetnVertices(),
                     inserter(CommonNodes,CommonNodes.begin()));
    if(!CommonNodes.empty()){
      filePrint(stderr," *** ERROR: FOUND %d COMMON NODES BETWEEN MASTER-SLAVE SURFACES: CASE NOT SUPPORTED\n",
            CommonNodes.size());
      for(int i=0; i<int(CommonNodes.size()); ++i)
        filePrint(stderr,"%6d ",CommonNodes[i]+1);
    }
  }
#endif
 
  // -> !!! UNDERLAYING LINEAR FACE ELEMENT NODES !!!
  int nMasterFaceElem, nSlaveFaceElem;
  PtrGlobalMasterEntity = PtrMasterEntity;
  PtrGlobalSlaveEntity = PtrSlaveEntity;

#ifdef DISTRIBUTED
  Connectivity *nodeToCpu = 0;
  if(DIST_ACME == 1) {
    nMasterFaceElem = (structCom->myID() == 0) ? MasterFaceElemSet->nElems() : 0;
    nSlaveFaceElem  = (structCom->myID() == 0) ? SlaveFaceElemSet->nElems() : 0;
    nMasterNodes    = (structCom->myID() == 0) ? PtrMasterEntity->GetnVertices() : 0; 
    nSlaveNodes     = (structCom->myID() == 0) ? PtrSlaveEntity->GetnVertices() : 0; 
  }
  else if(DIST_ACME == 2) {
    Connectivity *masterFaceElemToNode = new Connectivity(SetAccess<FaceElemSet>{*MasterFaceElemSet});
    SurfaceEntity *PtrLocalMasterEntity = new SurfaceEntity(-1);
    int num = 0;
    for(int i=0; i<masterFaceElemToNode->csize(); ++i) { 
      bool localEle = true;
      int *glNodes = new int[masterFaceElemToNode->num(i)];
      for(int j=0; j<masterFaceElemToNode->num(i); ++j) {
        glNodes[j] = PtrMasterEntity->GetGlNodeId((*masterFaceElemToNode)[i][j]);
        bool localNode = false;
        for(int k=0; k<numSub; ++k) if(sd[k]->globalToLocal(glNodes[j]) > -1) localNode = true;
        if(!localNode) localEle = false; 
      }
      if(localEle) PtrLocalMasterEntity->AddFaceElement(num++, (*MasterFaceElemSet)[i]->GetFaceElemType(), masterFaceElemToNode->csize(), glNodes);
      delete [] glNodes;
    }
    delete masterFaceElemToNode;
    PtrLocalMasterEntity->SetUpData(&geoSource->GetNodes());
    PtrLocalMasterEntity->Renumber();
    nMasterFaceElem = PtrLocalMasterEntity->GetPtrFaceElemSet()->nElems();
    nMasterNodes = PtrLocalMasterEntity->GetnVertices();

    Connectivity *slaveFaceElemToNode = new Connectivity(SetAccess<FaceElemSet>{*SlaveFaceElemSet});
    SurfaceEntity *PtrLocalSlaveEntity = new SurfaceEntity(-2);
    num = 0;
    for(int i=0; i<slaveFaceElemToNode->csize(); ++i) {
      bool localEle = true;
      int *glNodes = new int[slaveFaceElemToNode->num(i)];
      for(int j=0; j<slaveFaceElemToNode->num(i); ++j) {
        glNodes[j] = PtrSlaveEntity->GetGlNodeId((*slaveFaceElemToNode)[i][j]);
        bool localNode = false;
        for(int k=0; k<numSub; ++k) if(sd[k]->globalToLocal(glNodes[j]) > -1) localNode = true; 
        if(!localNode) localEle = false;
      }
      if(localEle) PtrLocalSlaveEntity->AddFaceElement(num++, (*SlaveFaceElemSet)[i]->GetFaceElemType(), slaveFaceElemToNode->csize(), glNodes);
      delete [] glNodes;
    }
    delete slaveFaceElemToNode;
    PtrLocalSlaveEntity->SetUpData(&geoSource->GetNodes());
    PtrLocalSlaveEntity->Renumber();
    nSlaveFaceElem = PtrLocalSlaveEntity->GetPtrFaceElemSet()->nElems();
    nSlaveNodes = PtrLocalSlaveEntity->GetnVertices();
  
    PtrMasterEntity = PtrLocalMasterEntity;
    PtrSlaveEntity = PtrLocalSlaveEntity;
  
    // make cpuToNode
    int *nn_perCPU = new int[structCom->numCPUs()];
    for(int i=0; i<structCom->numCPUs(); ++i) nn_perCPU[i] = 0;
    nn_perCPU[structCom->myID()] = nMasterNodes+nSlaveNodes;
    structCom->globalSum(structCom->numCPUs(), nn_perCPU);
    int *ptr =  new int[structCom->numCPUs()+1];
    ptr[0] = 0;
    for(int i=0; i<structCom->numCPUs(); ++i) ptr[i+1] = ptr[i] + nn_perCPU[i];
    delete [] nn_perCPU;
    int *tgt = new int[ptr[structCom->numCPUs()]];
    for(int i=0; i<ptr[structCom->numCPUs()]; ++i) tgt[i] = 0;
    for(int i=0; i<nMasterNodes; ++i) 
      tgt[ptr[structCom->myID()]+i] = PtrLocalMasterEntity->GetGlVertexId(i); 
    for(int i=0; i<nSlaveNodes; ++i) 
      tgt[ptr[structCom->myID()]+nMasterNodes+i] = PtrLocalSlaveEntity->GetGlVertexId(i);
    structCom->globalSum(ptr[structCom->numCPUs()],tgt);
    Connectivity *cpuToNode = new Connectivity(structCom->numCPUs(), ptr, tgt);
    nodeToCpu = cpuToNode->alloc_reverse();
    Connectivity *cpuToCpu = cpuToNode->transcon(nodeToCpu);
  
    // now make the interface lists
    int iCpu, jCpu, cpuJ, iNode;
    int totConnect;
    int *nConnect = new int[cpuToNode->csize()];
    int *flag = new int[cpuToNode->csize()];
    for(iCpu = 0; iCpu < cpuToNode->csize(); ++iCpu) {
       flag[iCpu] = -1;
       nConnect[iCpu] = 0;
    }
  
    // Count connectivity
    totConnect = 0;
    for(iCpu = 0; iCpu < cpuToNode->csize(); ++iCpu) {
      for(iNode = 0; iNode < cpuToNode->num(iCpu); ++iNode) { // loop over the nodes
        int thisNode = (*cpuToNode)[iCpu][iNode];
        for(jCpu = 0; jCpu < nodeToCpu->num(thisNode); ++jCpu) {
          // loop over the subdomains connected to this node
          cpuJ = (*nodeToCpu)[thisNode][jCpu];
          if(cpuJ > iCpu) {
            // only deal with connection to highered numbered subdomains, guarantees symmetry of lists
            if(flag[cpuJ] != iCpu) {
              flag[cpuJ] = iCpu;
              nConnect[iCpu]++;
              nConnect[cpuJ]++;
              totConnect += 2;
            }
          }
        }
      }
    }
  
    int **nodeCount = new int*[cpuToNode->csize()];
    int **connectedDomain = new int*[cpuToNode->csize()];
    int **remoteID = new int*[cpuToNode->csize()];
  
    // Allocate memory for list of connected subdomains
    for(iCpu = 0; iCpu < cpuToNode->csize(); ++iCpu) {
      int size = nConnect[iCpu];
      connectedDomain[iCpu] = new int[size];
      remoteID[iCpu] = new int[size];
      nodeCount[iCpu] = new int[size];
      flag[iCpu] = -1;
      nConnect[iCpu] = 0;
    }
  
    int *whichLocal  = new int[cpuToNode->csize()];
    int *whichRemote = new int[cpuToNode->csize()];
  
    for(iCpu=0; iCpu < cpuToNode->csize(); ++iCpu) {
      for(iNode = 0; iNode < cpuToNode->num(iCpu); ++iNode) {
        int nd = (*cpuToNode)[iCpu][iNode];
        for(jCpu = 0; jCpu < nodeToCpu->num(nd); ++jCpu) {
          cpuJ = (*nodeToCpu)[nd][jCpu];
          if(cpuJ > iCpu) {
            if(flag[cpuJ] != iCpu) { // attribute location for this sub
              flag[cpuJ] = iCpu;
              connectedDomain[cpuJ][nConnect[cpuJ]] = iCpu;
              connectedDomain[iCpu][nConnect[iCpu]] = cpuJ;
              remoteID[cpuJ][nConnect[cpuJ]] = nConnect[iCpu];
              remoteID[iCpu][nConnect[iCpu]] = nConnect[cpuJ];
              whichLocal[cpuJ] = nConnect[iCpu]++;
              whichRemote[cpuJ] = nConnect[cpuJ]++;
              nodeCount[iCpu][whichLocal[cpuJ]] = 1;
              nodeCount[cpuJ][whichRemote[cpuJ]] = 1;
            }
            else {
              nodeCount[iCpu][whichLocal[cpuJ]]++;
              nodeCount[cpuJ][whichRemote[cpuJ]]++;
            }
          }
        }
      }
    }
  
    // allocate memory for interface node lists
    Connectivity **interfNode = new Connectivity *[cpuToNode->csize()];
    for(iCpu=0; iCpu < cpuToNode->csize(); ++iCpu)
      interfNode[iCpu] = new Connectivity(nConnect[iCpu], nodeCount[iCpu]);
  
    // fill the lists
    for(iCpu = 0; iCpu < cpuToNode->csize(); ++iCpu) {
       flag[iCpu]     = -1;
       nConnect[iCpu] =  0;
    }
  
    for(iCpu = 0; iCpu < cpuToNode->csize(); ++iCpu) {
      for(iNode = 0; iNode < cpuToNode->num(iCpu); ++iNode) {
        int nd = (*cpuToNode)[iCpu][iNode];
        for(jCpu = 0; jCpu < nodeToCpu->num(nd); ++jCpu) {
          cpuJ = (*nodeToCpu)[nd][jCpu];
          if(cpuJ > iCpu) {
            if(flag[cpuJ] != iCpu) { // attribute location for this sub
              flag[cpuJ] = iCpu;
              whichLocal[cpuJ] = nConnect[iCpu]++;
              whichRemote[cpuJ] = nConnect[cpuJ]++;
              (*interfNode[iCpu])[whichLocal[cpuJ]][0] = nd;
              (*interfNode[cpuJ])[whichRemote[cpuJ]][0] = nd;
              nodeCount[iCpu][whichLocal[cpuJ]]=1;
              nodeCount[cpuJ][whichRemote[cpuJ]]=1;
            }
            else {
              int il = nodeCount[iCpu][whichLocal[cpuJ]]++;
              (*interfNode[iCpu])[whichLocal[cpuJ]][il] = nd;
              int jl = nodeCount[cpuJ][whichRemote[cpuJ]]++;
              (*interfNode[cpuJ])[whichRemote[cpuJ]][jl] = nd;
            }
          }
        }
      }
    }
  
    num_comm_partners = nConnect[structCom->myID()];
    comm_proc_id = new int[num_comm_partners]; for(int i=0; i<num_comm_partners; ++i) comm_proc_id[i] = connectedDomain[structCom->myID()][i];
    number_nodes_to_partner = new int[num_comm_partners]; for(int i=0; i<num_comm_partners; ++i) number_nodes_to_partner[i] = interfNode[structCom->myID()]->num(i);
    comm_node = new int[interfNode[structCom->myID()]->numConnect()];
    
    std::map<int, int> *slave_entity = PtrLocalSlaveEntity->GetPtrGlToLlVertexMap();
    std::map<int, int> *master_entity = PtrLocalMasterEntity->GetPtrGlToLlVertexMap();
    std::map<int, int>::iterator it;
    for(int i=0; i<interfNode[structCom->myID()]->numConnect(); ++i) {
      int glnode = interfNode[structCom->myID()]->getTarget()[i];
      int locnode = -1;
      if((it = slave_entity->find(glnode)) == slave_entity->end()) {
        if((it = master_entity->find(glnode)) == master_entity->end()) {
          std::cerr<<glnode<<" is not in the slave or master entity!"<<std::endl;
        }
        else { locnode = it->second + nSlaveNodes + 1; }
      }
      else { locnode = it->second + 1; }
      comm_node[i] = locnode;
    }
  
    delete [] flag; delete [] whichLocal; delete [] whichRemote;
    for(iCpu = 0; iCpu < cpuToNode->csize(); ++iCpu) { delete [] nodeCount[iCpu]; delete [] remoteID[iCpu]; delete [] connectedDomain[iCpu]; delete interfNode[iCpu]; }
    delete [] nodeCount; delete [] remoteID; delete [] connectedDomain; delete [] nConnect;  delete [] interfNode;
    delete cpuToNode; delete cpuToCpu;
    
#ifdef MORTAR_DEBUG
    for(int i=0; i<structCom->numCPUs(); ++i) {
      if(i == structCom->myID()) {
        std::cerr << "cpu#" << i << ", num_comm_partners = " << num_comm_partners << ", comm_proc_id = ";
        for(int j=0; j<num_comm_partners; ++j) std::cerr << comm_proc_id[j] << " "; std::cerr << std::endl;
        std::cerr << "number_nodes_to_partner = "; 
        for(int j=0; j<num_comm_partners; ++j) std::cerr << number_nodes_to_partner[j] << " "; std::cerr << std::endl;
        std::cerr << "comm_node = ";
        for(int j=0; j<interfNode[i]->numConnect(); ++j) std::cerr << comm_node[j] << " "; std::cerr << std::endl;
        std::cerr << "interfNode = \n"; interfNode[structCom->myID()]->print();
      } 
      structCom->sync();
    }
#endif
  } else 
#endif
  {
    nMasterFaceElem = MasterFaceElemSet->nElems();
    nSlaveFaceElem  = SlaveFaceElemSet->nElems();
    nMasterNodes    = PtrMasterEntity->GetnVertices();
    nSlaveNodes     = PtrSlaveEntity->GetnVertices();
  }

#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * nMasterFaces/nSlaveFaces       = %6d / %6d\n",nMasterFaceElem,nSlaveFaceElem);
  filePrint(stderr,"   * nMasterVertices/nSlaveVertices = %6d / %6d\n",nMasterNodes,nSlaveNodes);
#endif

  // Set ACME node coordinates array
  int nACMENodes  = nMasterNodes + nSlaveNodes; // !! Assume NO COMMON nodes !!  
  num_nodes_per_block[0] = nACMENodes;
  double* ACMENodesCoord = new double[3*nACMENodes];

  if(node_global_ids) delete [] node_global_ids;
  node_global_ids = new int[2*nACMENodes];
  int* node_exodus_ids = new int[nACMENodes];

  std::map<int,int> SlaveNdIdsToACMENdIds; 
  std::map<int,int> MasterNdIdsToACMENdIds; 

#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * fill the ACME node coordinates array\n"); 
#endif
  int inode = 0;
  int lnodeId;
  for(int ivertex=0; ivertex<nSlaveNodes; ++ivertex, ++inode) {
    Node &node = PtrSlaveEntity->GetVertex(ivertex); 
    ACMENodesCoord[3*inode  ] = node.x;
    ACMENodesCoord[3*inode+1] = node.y;
    ACMENodesCoord[3*inode+2] = node.z;
    int glnode = PtrSlaveEntity->GetGlVertexId(ivertex);
    node_exodus_ids[inode]    = glnode + 1;
#ifdef DISTRIBUTED
    if(DIST_ACME == 2) {
      int mastercpu = (*nodeToCpu)[glnode][0];
      for(int i=1; i<nodeToCpu->num(glnode); ++i) if((*nodeToCpu)[glnode][i] < mastercpu) mastercpu = (*nodeToCpu)[glnode][i];
      node_global_ids[2*inode]  = mastercpu;
    } else
#endif
    node_global_ids[2*inode]  = 0; 
    node_global_ids[2*inode+1]= glnode + 1;
    lnodeId = PtrSlaveEntity->GetLlVertexInLlNode(ivertex);
    SlaveNdIdsToACMENdIds[lnodeId]= inode+1; // ACME uses Fortran ordering
  }
  for(int ivertex=0; ivertex<nMasterNodes; ++ivertex, ++inode) {
    Node &node = PtrMasterEntity->GetVertex(ivertex);
    ACMENodesCoord[3*inode  ] = node.x;
    ACMENodesCoord[3*inode+1] = node.y;
    ACMENodesCoord[3*inode+2] = node.z;
    int glnode = PtrMasterEntity->GetGlVertexId(ivertex);
    node_exodus_ids[inode]    = glnode + 1;
#ifdef DISTRIBUTED
    if(DIST_ACME == 2) {
      int mastercpu = (*nodeToCpu)[glnode][0];
      for(int i=1; i<nodeToCpu->num(glnode); ++i) if((*nodeToCpu)[glnode][i] < mastercpu) mastercpu = (*nodeToCpu)[glnode][i];
      node_global_ids[2*inode]  = mastercpu;
    }
    else 
#endif
    node_global_ids[2*inode]  = 0;
    node_global_ids[2*inode+1]= glnode + 1;
    lnodeId = PtrMasterEntity->GetLlVertexInLlNode(ivertex);
    MasterNdIdsToACMENdIds[lnodeId]= inode+1; // ACME uses Fortran ordering
  }
 
  // Set ACME face topologies 
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * fill the ACME face topology array\n");
  filePrint(stderr,"   -> block[0] is ACME master <=> OUR slave \n");
  filePrint(stderr,"   -> block[1] is ACME slave  <=> OUR master\n");
#endif
  int* number_faces_per_block = new int[number_face_blocks];
  Connectivity* SlaveACMEBlocksMap = PtrSlaveEntity->GetPtrACMEBlocksMap();
  Connectivity* MasterACMEBlocksMap= PtrMasterEntity->GetPtrACMEBlocksMap();

#ifdef DISTRIBUTED
  if(DIST_ACME == 1) {
    number_faces_per_block[0] = (structCom->myID() == 0) ? SlaveACMEBlocksMap->num(0) : 0; // number of slave Quad4 face els 
    number_faces_per_block[1] = (structCom->myID() == 0) ? SlaveACMEBlocksMap->num(1) : 0; // number of slave Quad8 face els 
    number_faces_per_block[2] = (structCom->myID() == 0) ? SlaveACMEBlocksMap->num(2) : 0; // number of slave Tri3 face els 
    number_faces_per_block[3] = (structCom->myID() == 0) ? SlaveACMEBlocksMap->num(3) : 0; // number of slave Tri6 face els 

    number_faces_per_block[4] = (structCom->myID() == 0) ? MasterACMEBlocksMap->num(0) : 0; // number of master Quad4 face els 
    number_faces_per_block[5] = (structCom->myID() == 0) ? MasterACMEBlocksMap->num(1) : 0; // number of master Quad8 face els 
    number_faces_per_block[6] = (structCom->myID() == 0) ? MasterACMEBlocksMap->num(2) : 0; // number of master Tri3 face els 
    number_faces_per_block[7] = (structCom->myID() == 0) ? MasterACMEBlocksMap->num(3) : 0; // number of master Tri6 face els 
  }
  else 
#endif
  {
    number_faces_per_block[0] = SlaveACMEBlocksMap->num(0); // number of slave Quad4 face els 
    number_faces_per_block[1] = SlaveACMEBlocksMap->num(1); // number of slave Quad8 face els 
    number_faces_per_block[2] = SlaveACMEBlocksMap->num(2); // number of slave Tri3 face els 
    number_faces_per_block[3] = SlaveACMEBlocksMap->num(3); // number of slave Tri6 face els 

    number_faces_per_block[4] = MasterACMEBlocksMap->num(0); // number of master Quad4 face els 
    number_faces_per_block[5] = MasterACMEBlocksMap->num(1); // number of master Quad8 face els 
    number_faces_per_block[6] = MasterACMEBlocksMap->num(2); // number of master Tri3 face els 
    number_faces_per_block[7] = MasterACMEBlocksMap->num(3); // number of master Tri6 face els 
  }
 
  int  size_face_connectivity = 4*(number_faces_per_block[0] + number_faces_per_block[1])
                              + 3*(number_faces_per_block[2] + number_faces_per_block[3])  
                              + 4*(number_faces_per_block[4] + number_faces_per_block[5])
                              + 3*(number_faces_per_block[6] + number_faces_per_block[7]);
  
#ifdef MORTAR_DEBUG
  for(int i=0;i<number_face_blocks;i++)
    filePrint(stderr," * ACME face block %d contains %6d face els.\n",i,number_faces_per_block[i]);
  filePrint(stderr," * Total ACME face_connectivity = %6d\n",size_face_connectivity);
#endif
  int* face_connectivity = new int[size_face_connectivity];

  int offset;
#ifdef DISTRIBUTED
  if(DIST_ACME == 1) {
    offset = (structCom->myID() == 0) ? PtrSlaveEntity->FillACMEFaceBlocks(face_connectivity, SlaveNdIdsToACMENdIds, true) : 0;
    offset += ((structCom->myID() == 0) ? PtrMasterEntity->FillACMEFaceBlocks(&face_connectivity[offset], MasterNdIdsToACMENdIds, true) : 0);
  } else
#endif
  {
    // first: the slave face elements' connectivity (UNDERLAYING LINEAR FACE ELEMENT NODES)
    offset = PtrSlaveEntity->FillACMEFaceBlocks(face_connectivity, SlaveNdIdsToACMENdIds, true);
    // second: the master face elements' connectivity (UNDERLAYING LINEAR FACE ELEMENT NODES)
    offset += PtrMasterEntity->FillACMEFaceBlocks(&face_connectivity[offset], MasterNdIdsToACMENdIds, true);
  }
 
  // map ACME topology to ACME node numbering
#ifdef MORTAR_DEBUG
  filePrint(stderr," * offset = %6d\n",offset);
  for(int i=0;i<size_face_connectivity;i++)
    filePrint(stderr," * ACME face_connectivity[%6d] = %6d\n",i,face_connectivity[i]);
  // create a top file with the ACME blocks
  char fileName[128]; 
#ifdef DISTRIBUTED
  sprintf(fileName,"%s_%d_%d%s","AcmeMortar",GetId(),structCom->myID(),".top");
#else
  sprintf(fileName,"%s_%d%s","AcmeMortar",GetId(),".top");
#endif
  filePrint(stderr," * Write ACME blocks to top file %s\n",fileName);
  FILE* AcmeTop = fopen(fileName,"w");
  fprintf(AcmeTop,"Nodes AcmeNodes\n");
  for(int i=0; i<nACMENodes; i++) 
    fprintf(AcmeTop," %6d %e  %e  %e\n",i+1,ACMENodesCoord[3*i],ACMENodesCoord[3*i+1],ACMENodesCoord[3*i+2]);
  offset = 0;
  for(int iblk=0; iblk<number_face_blocks; iblk++) { 
    if(!number_faces_per_block[iblk]) continue;
    fprintf(AcmeTop,"Elements AcmeBlock_%d using AcmeNodes\n",iblk+1);
    int faceType = face_block_types[iblk];
    int faceTopType;
    switch(faceType) {
      case(ContactSearch::QUADFACEL4):
      case(ContactSearch::QUADFACEQ8):
      case(ContactSearch::SHELLQUADFACEL4):
      case(FaceElement::QUADFACEC12):
        faceTopType = 2;
        break;
      case(ContactSearch::TRIFACEL3):
      case(ContactSearch::TRIFACEQ6):
      case(ContactSearch::SHELLTRIFACEL3):
      case(FaceElement::TRIFACEC10):
        faceTopType = 4;
        break;
    }
    int nVertices = nVerticesPerFaceType[faceType];
    for(int i=0; i<number_faces_per_block[iblk]; i++) {
      fprintf(AcmeTop,"%6d  %3d  ",i+1,faceTopType);
      for(int j=0; j<nVertices; j++) 
        fprintf(AcmeTop," %6d ",face_connectivity[offset++]);   
      fprintf(AcmeTop,"\n");
    } 
  }
  fclose(AcmeTop);
#endif
 
  int* face_global_ids = new int[2*(nSlaveFaceElem+nMasterFaceElem)];
  for(int i=0; i<nSlaveFaceElem+nMasterFaceElem; i++) {
#ifdef DISTRIBUTED
    if(DIST_ACME == 2) {
      face_global_ids[2*i] = structCom->myID();
    } else 
#endif
    face_global_ids[2*i] = 0;
    face_global_ids[2*i+1] = i+1;
  } 

#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * ACME informations\n");
  filePrint(stderr,"    # dimensionality       : %d\n",dimensionality);
  filePrint(stderr,"    # nb of states         : %d\n",num_states);
  filePrint(stderr,"    # nb entity keys       : %d\n",num_entity_keys);
  filePrint(stderr,"    # nb node blck         : %d\n",num_node_blocks);
  filePrint(stderr,"    # node blck type       : %d\n",node_block_types[0]);
  filePrint(stderr,"    # nb node/blck         : %d\n",num_nodes_per_block[0]);
  filePrint(stderr,"    # nb face blck         : %d\n",number_face_blocks);
  filePrint(stderr,"    # nb face/blck         : %d - %d\n",number_faces_per_block[0],number_faces_per_block[1]);
  filePrint(stderr,"    # nb elem blck         : %d\n",num_element_blocks);
  filePrint(stderr,"    # nb comm part.        : %d\n",num_comm_partners);
#endif

  int tot_num_faces = 0; 
  for(int i = 0; i < number_face_blocks; ++i) tot_num_faces += number_faces_per_block[i];
  double *face_lofting_factors = new double[tot_num_faces];
  for(int i = 0; i < tot_num_faces; ++i) face_lofting_factors[i] = 0;

  // Build ACME search object
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * build ACME search object\n");
#endif
  ContactSearch::ContactErrorCode error;
  if(search_obj) delete search_obj;
  search_obj = new ContactSearch( dimensionality,
                            num_states,
                            num_analytical_surfs,
                            num_node_blocks,
                            node_block_types,
                            &(num_nodes_per_block[0]),
                            node_exodus_ids,
                            node_global_ids,
                            ACMENodesCoord,
                            number_face_blocks,
                            face_block_types,
                            &(number_faces_per_block[0]),
                            face_global_ids,
                            face_connectivity,
                            face_lofting_factors,
                            num_element_blocks,
                            element_block_types,
                            number_elements_per_block,
                            element_ids,
                            element_connectivity,
                            num_comm_partners,
                            comm_proc_id,
                            number_nodes_to_partner,
                            comm_node,
                            mpi_Comm,
                            error );
  
  if(error) {
    std::cerr << "Error in ACME ContactSearch::ContactSearch: error code = " << error << std::endl;
    for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
      std::cerr << search_obj->Error_Message(i) << std::endl;
    exit(error);
  }

  // set face block attributes (shell thickness and lofting)
  if(tdenforceFlag) {
    for(int i = 0; i < number_face_blocks; i++) {
      if(face_block_types[i] == ContactSearch::SHELLQUADFACEL4 || face_block_types[i] == ContactSearch::SHELLTRIFACEL3) {
        double *attributes = new double[number_faces_per_block[i]];
        for(int j = 0; j < number_faces_per_block[i]; ++j) attributes[j] = face_block_shell_thickness[i];
        ContactSearch::ContactErrorCode error;
        error = search_obj->Set_Face_Block_Attributes(ContactSearch::SHELL_THICKNESS, i+1, attributes);
        if(error) {
          std::cerr << "Error in ACME ContactSearch::Set_Face_Block_Attributes: error code = " << error << std::endl;
          for(int j = 1; j <= search_obj->Number_of_Errors(); ++j)
            std::cerr << search_obj->Error_Message(j) << std::endl;
          exit(error);
        }

        for(int j = 0; j < number_faces_per_block[i]; ++j) attributes[j] = 0.5;
        error = search_obj->Set_Face_Block_Attributes(ContactSearch::LOFTING_FACTOR, i+1, attributes);
        if(error) {
          std::cerr << "Error in ACME ContactSearch::Set_Face_Block_Attributes: error code = " << error << std::endl;
          for(int j = 1; j <= search_obj->Number_of_Errors(); ++j)
            std::cerr << search_obj->Error_Message(j) << std::endl;
          exit(error);
        }

        delete [] attributes;
      }
    }
  }

  if(node_block_types)       delete [] node_block_types;
  if(ACMENodesCoord)         delete [] ACMENodesCoord;
  if(face_block_types)       delete [] face_block_types;
  if(node_exodus_ids)        delete [] node_exodus_ids;
  if(face_connectivity)      delete [] face_connectivity;
  if(face_global_ids)        delete [] face_global_ids;
  if(face_lofting_factors)   delete [] face_lofting_factors;
  if(number_faces_per_block) delete [] number_faces_per_block;
  if(face_block_shell_thickness) delete [] face_block_shell_thickness;

#endif
}

void
MortarHandler::set_search_data(int interaction_type)
{
#ifdef USE_ACME
  // Set ACME interaction to look for
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * Set ACME active interactions\n");
#endif
  double Normal_Tol = GetNormalTol();
  if(interaction_type == 4) {
    // the tolerance specified in the input file is measured w.r.t the lofted surfaces
    // this makes the interpretation of the normal tolerance consistent from a user's perspective 
    if(SelfContact) Normal_Tol += PtrSlaveEntity->GetShellThickness();
    else Normal_Tol += (PtrMasterEntity->GetShellThickness()/2 + PtrSlaveEntity->GetShellThickness()/2);
  }
  double Tangential_Tol = GetTangentialTol();
  double Interaction_Typ;
  switch(interaction_type) {
    case 1 : Interaction_Typ = (double)(ContactSearch::SLIDING_INTERACTION); break;
    case 2 : Interaction_Typ = (double)(ContactSearch::TIED_INTERACTION); break;
    case 3 : case 4 : case 6 : Interaction_Typ = (double)(ContactSearch::GENERIC_INTERACTION); break;
    case 5 : Interaction_Typ = (double)(ContactSearch::COVERAGE_INTERACTION); break;
  }

  int step_size = ContactSearch::NSIZSD;

  double *Search_Data = new double[step_size*num_entity_keys*num_entity_keys];
  memset(Search_Data, 0, sizeof(double)*step_size*num_entity_keys*num_entity_keys);

  for(int i=0; i<num_entity_keys*num_entity_keys; i++) {
    Search_Data[step_size*i  ] = (double)(ContactSearch::NO_INTERACTION);
    Search_Data[step_size*i+1] = Normal_Tol;
    Search_Data[step_size*i+2] = Tangential_Tol;
  }
  for(int i=0; i<4; i++) {
    if(SelfContact) { 
      Search_Data[step_size*( 0+i)] = Interaction_Typ;
      Search_Data[step_size*( 8+i)] = Interaction_Typ;
      Search_Data[step_size*(16+i)] = Interaction_Typ;
      Search_Data[step_size*(24+i)] = Interaction_Typ;
    }
    else {
      Search_Data[step_size*( 4+i)] = Interaction_Typ;
      Search_Data[step_size*(12+i)] = Interaction_Typ;
      Search_Data[step_size*(20+i)] = Interaction_Typ;
      Search_Data[step_size*(28+i)] = Interaction_Typ;
    }
  }

#ifdef MORTAR_DEBUG
  filePrint(stderr,"    # normal tol.    : %e\n",Normal_Tol);
  filePrint(stderr,"    # tangential tol.: %e\n",Tangential_Tol);
  filePrint(stderr,"    # Search_Data array:\n");
  for(int i=0; i<8; i++){
    filePrint(stderr,"    --------- \n");
    for(int j=0; j<8; j++){
      filePrint(stderr,"     -> FB %d - FB %d",i,j);
      if(Search_Data[step_size*(j+8*i)]==(double)(ContactSearch::NO_INTERACTION))
        filePrint(stderr,", NO_INTERACTION             ");
      if(Search_Data[step_size*(j+8*i)]==(double)(ContactSearch::SLIDING_INTERACTION))
        filePrint(stderr,", SLIDING_INTERACTION        ");
      if(Search_Data[step_size*(j+8*i)]==(double)(ContactSearch::GENERIC_INTERACTION))
        filePrint(stderr,", GENERIC_INTERACTION        ");
      if(Search_Data[step_size*(j+8*i)]==(double)(ContactSearch::COVERAGE_INTERACTION))
        filePrint(stderr,", COVERAGE_INTERACTION       ");
      filePrint(stderr,", normal tol = %e, tang. tol = %e\n", Search_Data[step_size*(j+8*i)+1], Search_Data[step_size*(j+8*i)+2]);
    }
  }
#endif

  search_obj->Set_Search_Data(Search_Data);

  ContactSearch::ContactErrorCode error;
  error = search_obj->Check_Search_Data_Size(step_size,num_entity_keys);
  if(error) {
    std::cerr << "Error in ACME ContactSearch::Check_Search_Data_Size: error code = " << error << std::endl;
    for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
      std::cerr << search_obj->Error_Message(i) << std::endl;
    exit(error);
  }

  if(Search_Data) delete [] Search_Data;
#endif
}
 
void
MortarHandler::set_node_configuration(int config_type)
{
#ifdef USE_ACME
  int num_nodes = PtrMasterEntity->GetnVertices() + PtrSlaveEntity->GetnVertices(); 
  double* positions = new double[3*num_nodes];
  int inode = 0;

  for(int ivertex=0; ivertex<PtrSlaveEntity->GetnVertices(); ++ivertex, ++inode) {
    Node &node = PtrSlaveEntity->GetVertex(ivertex);
    positions[3*inode+0] = node.x;
    positions[3*inode+1] = node.y;
    positions[3*inode+2] = node.z;
  }
  for(int ivertex=0; ivertex<PtrMasterEntity->GetnVertices(); ++ivertex, ++inode) {
    Node &node = PtrMasterEntity->GetVertex(ivertex);
    positions[3*inode+0] = node.x;
    positions[3*inode+1] = node.y;
    positions[3*inode+2] = node.z;
  }

  ContactSearch::ContactNode_Configuration config = (config_type == 1) ? ContactSearch::CURRENT_CONFIG : ContactSearch::PREDICTED_CONFIG;
  int node_block_id = 1;

  ContactSearch::ContactErrorCode error;
  error = search_obj->Set_Node_Block_Configuration(config, node_block_id, positions);
  if(error) {
    std::cerr << "Error in ACME ContactSearch::Set_Node_Block_Configuration: error code = " << error << std::endl;
    for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
      std::cerr << search_obj->Error_Message(i) << std::endl;
    exit(error);
  }

  delete [] positions;
#endif
}

void
MortarHandler::set_node_configuration(int config_type, int numSub, SubDomain **sd)
{   
  // multiple domain version
#ifdef USE_ACME
  int num_nodes = PtrMasterEntity->GetnVertices() + PtrSlaveEntity->GetnVertices();
  double* positions = new double[3*num_nodes];
  for(int i=0; i<3*num_nodes; ++i) positions[i] = 0.0;

  for(int j=0; j<numSub; ++j) {
    int inode = 0;
    for(int ivertex=0; ivertex< PtrSlaveEntity->GetnVertices(); ++ivertex, ++inode) {
      int glNode = PtrSlaveEntity->GetGlVertexId(ivertex);
      int locNode = sd[j]->globalToLocal(glNode);
      if(locNode > -1) {
        Node &node = PtrSlaveEntity->GetVertex(ivertex);
#ifdef DISTRIBUTED
        if(DIST_ACME != 2) {
          positions[3*inode+0] += (node.x/double(share[3*inode+0]));
          positions[3*inode+1] += (node.y/double(share[3*inode+1]));
          positions[3*inode+2] += (node.z/double(share[3*inode+2]));
        } else 
#endif
        {
          positions[3*inode+0] = node.x;
          positions[3*inode+1] = node.y;
          positions[3*inode+2] = node.z;
        }
      }
    }
    if(PtrMasterEntity != PtrSlaveEntity) { // i.e. not mortar SelfContact
      for(int ivertex=0; ivertex<PtrMasterEntity->GetnVertices(); ++ivertex, ++inode) {
        int glNode = PtrMasterEntity->GetGlVertexId(ivertex);
        int locNode = sd[j]->globalToLocal(glNode);
        if(locNode > -1) {
          Node &node = PtrMasterEntity->GetVertex(ivertex);
#ifdef DISTRIBUTED
          if(DIST_ACME != 2) {
            positions[3*inode+0] += (node.x/double(share[3*inode+0]));
            positions[3*inode+1] += (node.y/double(share[3*inode+1]));
            positions[3*inode+2] += (node.z/double(share[3*inode+2]));
          } else
#endif
          {
            positions[3*inode+0] = node.x;
            positions[3*inode+1] = node.y;
            positions[3*inode+2] = node.z;
          }
        }
      }
    }
  }
#ifdef DISTRIBUTED
  if(DIST_ACME == 0) {
    structCom->globalSum(3*num_nodes, positions);
    int inode = 0;
    for(int ivertex=0; ivertex< PtrSlaveEntity->GetnVertices(); ++ivertex, ++inode) {
      Node &node = PtrSlaveEntity->GetVertex(ivertex);
      node.x = positions[3*inode+0];
      node.y = positions[3*inode+1];
      node.z = positions[3*inode+2];
    }
    if(PtrMasterEntity != PtrSlaveEntity) { // i.e. not mortar SelfContact
      for(int ivertex=0; ivertex<PtrMasterEntity->GetnVertices(); ++ivertex, ++inode) {
        Node &node = PtrMasterEntity->GetVertex(ivertex);
        node.x = positions[3*inode+0];
        node.y = positions[3*inode+1];
        node.z = positions[3*inode+2];
      }
    }
  }
  else if(DIST_ACME == 1) {
    structCom->reduce(3*num_nodes, positions, 0);
    if(structCom->myID() != 0) { delete [] positions; positions = 0; num_nodes = 0; }
  }
#endif

  ContactSearch::ContactNode_Configuration config = (config_type == 1) ? ContactSearch::CURRENT_CONFIG : ContactSearch::PREDICTED_CONFIG;
  int node_block_id = 1;

  ContactSearch::ContactErrorCode error;
  error = search_obj->Set_Node_Block_Configuration(config, node_block_id, positions);
  if(error) {
    std::cerr << "Error in ACME ContactSearch::Set_Node_Block_Configuration: error code = " << error << std::endl;
    for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
      std::cerr << search_obj->Error_Message(i) << std::endl;
    exit(error);
  }

  if(positions != 0) delete [] positions;
#endif
}

void
MortarHandler::set_node_constraints(int numDBC, BCond *dbc)
{
#ifdef USE_ACME
  if(numDBC == 0) return;

  int num_nodes = PtrMasterEntity->GetnVertices() + PtrSlaveEntity->GetnVertices();
  int *constraints_per_node = new int[num_nodes];
  // constraint_vector is a vector for each node that describes the constraint direc-
  // tion. If constraints_per_node is 0 or 3, this vector should be set to 0.      If
  // constraints_per_node is 1, this vector should be the constrained direction.   If
  // constraints_per_node is 2, this vector should be the unconstrained direction
  double* constraint_vector = new double[3*num_nodes];
  for(int i=0; i<num_nodes; ++i) constraints_per_node[i] = 0;
  for(int i=0; i<3*num_nodes; ++i) constraint_vector[i] = 0.0;

  std::map<int, int> *slave_entity = PtrSlaveEntity->GetPtrGlToLlVertexMap();
  std::map<int, int> *master_entity = PtrMasterEntity->GetPtrGlToLlVertexMap();
  std::map<int, int>::iterator it;
  for(int i=0; i<numDBC; ++i) {
    int glnode = dbc[i].nnum;
    int locnode = -1;
    if((it = slave_entity->find(glnode)) != slave_entity->end()) {
      locnode = it->second;
    }
    else if((it = master_entity->find(glnode)) != master_entity->end()) {
      locnode = it->second + nSlaveNodes;
      //cerr << " *** WARNING: Constrained nodes on master surface\n";
      //If you see this message, try switching master and slave surface ids under TIEDSURFACES/CONTACTSURFACES
    }
    if(locnode != -1 && (dbc[i].dofnum == 0 || dbc[i].dofnum == 1 || dbc[i].dofnum == 2)) {
      switch(constraints_per_node[locnode]) {
        case 0 : 
          constraint_vector[3*locnode+dbc[i].dofnum] = 1.0; 
          constraints_per_node[locnode] = 1;
          break;
        case 1 :
          if(constraint_vector[3*locnode+dbc[i].dofnum] == 0.0) {
            for(int j=0; j<3; ++j) 
              constraint_vector[3*locnode+j] = (j == dbc[i].dofnum || constraint_vector[3*locnode+j] == 1.0) ? 0.0 : 1.0;
            constraints_per_node[locnode] = 2;
          }
          break;
        case 2 :
          if(constraint_vector[3*locnode+dbc[i].dofnum] == 1.0) {
            for(int j=0; j<3; ++j) constraint_vector[3*locnode+j] = 0.0; 
            constraints_per_node[locnode] = 3;
          }
          break;
        case 3 :
          /* do nothing */
          break;
      }
    }
  }
  
  int node_block_id = 1;
  ContactSearch::ContactErrorCode error;
  error = search_obj->Set_Node_Block_Kinematic_Constraints(node_block_id, constraints_per_node, constraint_vector);
  if(error) {
    std::cerr << "Error in ACME ContactSearch::Set_Node_Block_Configuration: error code = " << error << std::endl;
    for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
      std::cerr << search_obj->Error_Message(i) << std::endl;
    exit(error);
  } 
  
  delete [] constraints_per_node;
  delete [] constraint_vector;
#endif
}

void
MortarHandler::set_node_constraints(int numSub, SubDomain **sd)
{
#ifdef USE_ACME
  int num_nodes = PtrMasterEntity->GetnVertices() + PtrSlaveEntity->GetnVertices();
  int *constraints_per_node = new int[num_nodes];
  // constraint_vector is a vector for each node that describes the constraint direc-
  // tion. If constraints_per_node is 0 or 3, this vector should be set to 0.      If
  // constraints_per_node is 1, this vector should be the constrained direction.   If
  // constraints_per_node is 2, this vector should be the unconstrained direction
  double* constraint_vector = new double[3*num_nodes];
  for(int i=0; i<num_nodes; ++i) constraints_per_node[i] = 0;
  for(int i=0; i<3*num_nodes; ++i) constraint_vector[i] = 0.0;

  std::map<int, int> *slave_entity = PtrSlaveEntity->GetPtrGlToLlVertexMap();
  std::map<int, int> *master_entity = PtrMasterEntity->GetPtrGlToLlVertexMap();
  std::map<int, int>::iterator it;
  for(int s=0; s<numSub; ++s) {
    for(int i=0; i<sd[s]->nDirichlet(); ++i) {
      BCond *dbc = sd[s]->getDBC();
      int glnode = sd[s]->localToGlobal(dbc[i].nnum);
      int locnode = -1;
      if((it = slave_entity->find(glnode)) != slave_entity->end()) {
        locnode = it->second;
      }
      else if((it = master_entity->find(glnode)) != master_entity->end()) {
        locnode = it->second + nSlaveNodes; // XXX this may not be correct for DIST_ACME > 0
        //cerr << " *** WARNING: ContactSearch::Set_Node_Block_Kinematic_Constraints doesn't seem to work when constrained nodes are on master surface.\n"
        //     << " ***          Try switching master and slave under TIEDSURFACES/CONTACTSURFACES\n";
      }
      if(locnode != -1 && (dbc[i].dofnum == 0 || dbc[i].dofnum == 1 || dbc[i].dofnum == 2)) {
        switch(constraints_per_node[locnode]) {
          case 0 :
            constraint_vector[3*locnode+dbc[i].dofnum] = 1.0;
            constraints_per_node[locnode] = 1;
            break;
          case 1 :
            if(constraint_vector[3*locnode+dbc[i].dofnum] == 0.0) {
              for(int j=0; j<3; ++j)
                constraint_vector[3*locnode+j] = (j == dbc[i].dofnum || constraint_vector[3*locnode+j] == 1.0) ? 0.0 : 1.0;
              constraints_per_node[locnode] = 2;
            }
            break;
          case 2 :
            if(constraint_vector[3*locnode+dbc[i].dofnum] == 1.0) {
              for(int j=0; j<3; ++j) constraint_vector[3*locnode+j] = 0.0;
              constraints_per_node[locnode] = 3;
            }
            break;
          case 3 :
            /* do nothing */
            break;
        }
      }
    }
  }
#ifdef DISTRIBUTED
  if(DIST_ACME == 0) {
    structCom->globalMax(num_nodes, constraints_per_node);
    structCom->globalMax(3*num_nodes, constraint_vector);
  }
  else if(DIST_ACME == 1) {
    structCom->reduce(num_nodes, constraints_per_node, 0, MPI_MAX);
    structCom->reduce(3*num_nodes, constraint_vector, 0, MPI_MAX);
  }
#endif

  int node_block_id = 1;
  ContactSearch::ContactErrorCode error;
  error = search_obj->Set_Node_Block_Kinematic_Constraints(node_block_id, constraints_per_node, constraint_vector);
  if(error) {
    std::cerr << "Error in ACME ContactSearch::Set_Node_Block_Configuration: error code = " << error << std::endl;
    for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
      std::cerr << search_obj->Error_Message(i) << std::endl;
    exit(error);
  }

  delete [] constraints_per_node;
  delete [] constraint_vector;
#endif
}


void
MortarHandler::set_search_options()
{
#ifdef USE_ACME
  ContactSearch::ContactErrorCode error;
  double data[3];
  
  // Deactivate Secondary Decomposition
  if(NoSecondary) {
    error = search_obj->Set_Search_Option(ContactSearch::NO_SECONDARY,
                                          ContactSearch::ACTIVE,
                                          (double *)0);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate multiple interations
  if(!domain->solInfo().no_multiple_interactions) {
    data[0] = domain->solInfo().sharp_non_sharp_angle; // default is 30.
    error = search_obj->Set_Search_Option(ContactSearch::MULTIPLE_INTERACTIONS,
                                          ContactSearch::ACTIVE,
                                          &data[0]);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate normal smoothing
  if(domain->solInfo().normal_smoothing) {
    data[0] = domain->solInfo().sharp_non_sharp_angle; // default is 30;
    data[1] = domain->solInfo().normal_smoothing_distance; // default is 0.5
    data[2] = (domain->solInfo().resolution_method == 0) ? ContactSearch::USE_NODE_NORMAL : ContactSearch::USE_EDGE_BASED_NORMAL; // default is 1
    error = search_obj->Set_Search_Option(ContactSearch::NORMAL_SMOOTHING,
                                          ContactSearch::ACTIVE,
                                          &data[0]);

    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;   
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i) 
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate simple lofting
  if(domain->solInfo().shell_simple_lofting) {
    error = search_obj->Set_Search_Option(ContactSearch::SHELL_SIMPLE_LOFTING,
                                          ContactSearch::ACTIVE,
                                          (double *)0);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate gap partitioning
  if(domain->solInfo().partition_gap) {
    error = search_obj->Set_Search_Option(ContactSearch::PARTITION_GAP,
                                          ContactSearch::ACTIVE,
                                          (double *)0);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate global search cull
  if(domain->solInfo().global_search_cull) {
    data[0] = ContactSearch::SLAVE_CULL;
    error = search_obj->Set_Search_Option(ContactSearch::GLOBAL_SEARCH_CULL,
                                          ContactSearch::ACTIVE,
                                          &data[0]);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate no warped volume
  if(domain->solInfo().no_warped_volume) {
    error = search_obj->Set_Search_Option(ContactSearch::NO_WARPED_VOLUME,
                                          ContactSearch::ACTIVE,
                                          (double *)0);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate auto tol
  if(domain->solInfo().auto_tol) {
    error = search_obj->Set_Search_Option(ContactSearch::AUTO_TOL,
                                          ContactSearch::ACTIVE,
                                          (double*)0);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate aggressive tolerances
  if(domain->solInfo().agressive_tolerances) {
    data[0] = 0.1;
    data[1] = 0.1;
    error = search_obj->Set_Search_Option(ContactSearch::AGGRESSIVE_TOLERANCES,
                                          ContactSearch::INACTIVE,
                                          &data[0]);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate skip physical faces
  if(domain->solInfo().skip_physical_faces) {
    data[0] = ContactSearch::PF_FACE_BASED;
    error = search_obj->Set_Search_Option(ContactSearch::SKIP_PHYSICAL_FACES,
                                          ContactSearch::INACTIVE,
                                          &data[0]);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Deactivate ghosting
  if(domain->solInfo().no_ghosting) {
    error = search_obj->Set_Search_Option(ContactSearch::NO_GHOSTING,
                                          ContactSearch::ACTIVE,
                                          (double*)0);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }

  // Activate ffi partials 
  if(InteractionType == MortarHandler::CTC) {
#ifdef USE_FFI_DERIVATIVES
    data[0] = 1;
    error = search_obj->Set_Search_Option(ContactSearch::COMPUTE_PARTIALS,
                                          ContactSearch::ACTIVE,
                                          &data[0]);
    if(error) {
      std::cerr << "Error in ACME ContactSearch::Set_Search_Option: error code = " << error << std::endl;
      for(int i=1; i<=search_obj->Number_of_Errors(); ++i)
        std::cerr << search_obj->Error_Message(i) << std::endl;
      exit(error);
    }
#endif
  }
#endif
}

void
MortarHandler::perform_search(int search_algorithm, double dt_old, double dt)
{
#ifdef PJSA_ACME_TIMER
  double tt = -getTime();
#endif
#ifdef USE_ACME
  ContactSearch::ContactErrorCode error;
  switch(search_algorithm) {
    default:
    case 1:
      search_obj->Delete_All_Interactions();
      error = search_obj->Static_Search_1_Configuration();
      break;
    case 2:
      error = search_obj->Static_Search_2_Configuration();
      break;
    case 3:
      search_obj->Set_Search_Option(ContactSearch::OLD_DYNAMIC_SEARCH, ContactSearch::ACTIVE, (double *)0);
      error = search_obj->Dynamic_Search_2_Configuration(dt_old, dt);
      break;
    case 4:
      error = search_obj->Dynamic_Search_2_Configuration(dt_old, dt);
      break;
  }  
  if(error) {
    std::cerr << "Error in ACME ContactSearch: error code = " << error << std::endl; 
    for(int i=1; i<=search_obj->Number_of_Errors(); ++i) 
      std::cerr << search_obj->Error_Message(i) << std::endl;
    exit(error);
  }
#else
  filePrint(stderr," !!! FATAL PB in MortarHandler::perform_search: USE_ACME FLAG WAS NOT DEFINED !!!\n");
  filePrint(stderr," => IF YOU GET TO THIS POINT, THIS MEANS YOU HAVE DEFINED SOME CONTACT SURFACE(S)\n");
  filePrint(stderr," => BUT IMPOSSIBLE TO DO CONTACT WITHOUT THE ACME LIB !!!\n");
  filePrint(stderr," => STOP EXECUTION \n");
  exit(-1);
#endif
#ifdef PJSA_ACME_TIMER
  std::cerr << "elapsed time in MortarHandler::perform_search = " << (tt+getTime())/1000 << " seconds\n";
#endif
}

void
MortarHandler::get_interactions(int interaction_type)
{
#ifdef PJSA_ACME_TIMER
  double tt = -getTime();
#endif
#ifdef USE_ACME
/*
  if(PtrMasterEntity->GetIsShellFace() || PtrSlaveEntity->GetIsShellFace()) {
    std::cerr << " *** ERROR: ACME does not currently support the extraction of interactions involving shells.\n";
    exit(-1);
  }
*/
  switch(interaction_type) {

    case 1 :  // NodeFace_Interactions
      std::cerr << "get NodeFace_Interactions not supported\n";
      break;

    case 2 :  // NodeSurface_Interactions
      std::cerr << "get NodeSurface_Interactions not supported\n";
      break;

    case 3 :  // NodeNode_Interactions
      std::cerr << "get NodeNode_Interactions not supported\n";
      break;

    default:
    case 4 : // FaceFace_Interactions
    {
      if(Slave_face_block_id) { delete [] Slave_face_block_id; Slave_face_block_id = 0; }
      if(Slave_face_index_in_block) { delete [] Slave_face_index_in_block; Slave_face_index_in_block = 0; }
      if(Master_face_block_id) { delete [] Master_face_block_id; Master_face_block_id = 0; }
      if(Master_face_index_in_block) { delete [] Master_face_index_in_block; Master_face_index_in_block = 0; }
      if(Slave_face_procs) { delete [] Slave_face_procs; Slave_face_procs = 0; }
      if(Master_face_procs) { delete [] Master_face_procs; Master_face_procs = 0; }
      if(ACMEFFI_index) { delete [] ACMEFFI_index; ACMEFFI_index = 0; }
      if(ACMEFFI_data) { delete [] ACMEFFI_data; ACMEFFI_data = 0; }

      int num_FFI = 0;
      int FFI_data_size = 0; 
#ifdef DISTRIBUTED
      if(structCom->myID() == 0) {
#endif
      search_obj->Size_FaceFace_Interactions(num_FFI,FFI_data_size);
#ifdef MORTAR_DEBUG
      filePrint(stderr,"   * Extract the FFI if any\n");
      filePrint(stderr,"    * nb of FFI     : %d\n",num_FFI);
      filePrint(stderr,"    * nb FFI data   : %d\n",FFI_data_size);
#endif
#ifdef DISTRIBUTED
      }
      structCom->broadcast(1, &num_FFI);
      structCom->broadcast(1, &FFI_data_size);
#endif
      nFFI         = num_FFI;
      nACMEFFIData = FFI_data_size;

      if(num_FFI){
	Slave_face_block_id       = new int[num_FFI];
	Slave_face_index_in_block = new int[num_FFI];
	Master_face_block_id      = new int[num_FFI];
	Master_face_index_in_block= new int[num_FFI];
	Slave_face_procs          = new int[num_FFI];
	Master_face_procs         = new int[num_FFI];
	ACMEFFI_index             = new int[num_FFI];
	ACMEFFI_data              = new double[FFI_data_size];

        // !!! HERE ACME MASTER <=> OUR SLAVE !!! 
        //          ACME SLAVE  <=> OUR MASTER  
#ifdef DISTRIBUTED
        if(structCom->myID() == 0)
#endif
	search_obj->Get_FaceFace_Interactions(Master_face_block_id,
                                              Master_face_index_in_block,
                                              Master_face_procs,
                                              Slave_face_block_id,
                                              Slave_face_index_in_block,
                                              Slave_face_procs,
                                              ACMEFFI_index,
                                              ACMEFFI_data);
#ifdef DISTRIBUTED
        structCom->broadcast(num_FFI, Master_face_block_id); 
        structCom->broadcast(num_FFI, Master_face_index_in_block);
        structCom->broadcast(num_FFI, Master_face_procs);
        structCom->broadcast(num_FFI, Slave_face_block_id);
        structCom->broadcast(num_FFI, Slave_face_index_in_block);
        structCom->broadcast(num_FFI, Slave_face_procs);
        structCom->broadcast(num_FFI, ACMEFFI_index);
        structCom->broadcast(FFI_data_size, ACMEFFI_data);
#endif
        //TODO: for self contact we should eliminate the redundant constraints.
      } else {
#ifdef MORTAR_DEBUG
        filePrint(stderr," *** WARNING: no Face-Face interaction found between master surf %2d & slave surf %2d\n",
                  MasterEntityId, SlaveEntityId);
#endif
      } 
    } break;

    case 5 :  // FaceCoverage_Interactions
      std::cerr << "get FaceCoverage_Interactions not supported\n";
      break;

    case 6 :  // ElementElement_Interactions
      std::cerr << "get ElementElement_Interactions not supported\n";
      break;
  }
#endif
#ifdef PJSA_ACME_TIMER
  std::cerr << "elapsed time in MortarHandler::get_interactions = " << (tt+getTime())/1000 << " seconds\n";
#endif
}

void
MortarHandler::build_td_enforcement()
{
#ifdef USE_ACME
  double *Enforcement_Data = new double[2*num_entity_keys*num_entity_keys];
  // this array stores the kinematic partition factor and friction model ID
  // a kinematic partition factor between 0 and 1 speciﬁes a ﬁxed kinematic partitioning, 
  // 2 tells acme to compute this for each interaction pair based on the physical data (density, wavespeed)
  // (see acme manual section 1.7)

  if(SelfContact) {
    for(int i=0; i<num_entity_keys*num_entity_keys; ++i) {
      Enforcement_Data[2*i+0] = 0.5;
      Enforcement_Data[2*i+1] = 1.0;
    }
  }
  else {
    for(int i=0; i<num_entity_keys*num_entity_keys; ++i) {
      Enforcement_Data[2*i+0] = 0.0;
      Enforcement_Data[2*i+1] = 1.0;
    } 

    double fixed_kpf = 1.0; // kinematic partitioning factor: this can be varied 0 < fixed_kpf <= 1.0
    for(int f_key=0; f_key<4; ++f_key) {
      for(int n_key = 4; n_key < 8; ++n_key) {
        Enforcement_Data[(f_key*num_entity_keys+n_key)*2 + 0] = (MortarType == STD) ? fixed_kpf : 2.0;
      }
    }
  }

  bool get_cvars = true;
  bool calc_plot_force = false;
  ContactSearch::ContactErrorCode error;
  if(contact_obj) delete contact_obj;
  if((ConstraintOptionsData && ConstraintOptionsData->lagrangeMult == 0 && ConstraintOptionsData->penalty != 0) ||
     (ConstraintOptionsData == NULL && domain->solInfo().lagrangeMult == 0 && domain->solInfo().penalty != 0)) {
    contact_obj = new ContactTDEnfPenalty(Enforcement_Data, search_obj, error, get_cvars, calc_plot_force);
  }
  else {
    contact_obj = new ContactTDEnforcement(Enforcement_Data, search_obj, error, get_cvars, calc_plot_force);
  }
  if(error) { 
    std::cerr << "Error in ACME ContactTDEnforcement::ContactTDEnforcement: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  } 

  ContactEnforcement::Enforcement_Model_Types type = (InteractionType == MortarHandler::TIED) ? ContactEnforcement::TD_TIED : (ContactEnforcement::Enforcement_Model_Types) FrictionModel;
#ifdef MORTAR_DEBUG
  switch(type) {
    case ContactEnforcement::TD_FRICTIONLESS: std::cerr << "Enforcement Model: FRICTIONLESS\n"; break;
    case ContactEnforcement::TD_TIED: std::cerr << "Enforcement Model: TIED\n"; break;
    case ContactEnforcement::TD_CONSTANT_FRICTION: std::cerr << "Enforcement Model: CONSTANT FRICTION (friction coeff = " << FrictionCoef[0] << ")\n"; break;
    case ContactEnforcement::TD_VELOCITY_DEPENDENT: std::cerr << "Enforcement Model: VELOCITY DEPENDENT (static coeff = " << FrictionCoef[0] 
           << ", dynamic coeff = " << FrictionCoef[1] << ", velocity decay = " << FrictionCoef[2] << ")\n"; break;
    case ContactEnforcement::TD_PRESSURE_DEPENDENT: std::cerr << "Enforcement Model: PRESSURE DEPENDENT (friction coeff = " << FrictionCoef[0] 
           << ", reference pressure = " << FrictionCoef[1] << ", offset pressure = " << FrictionCoef[2] << ", pressure exponent = " << FrictionCoef[3] << ")\n"; break;
    default: std::cerr << "Error in MortarHandler::build_td_enforcement: Enforcement model " << type << " not supported\n"; exit(-1); break;
  }
#endif

  int ID = 1;
  int* integer_data = 0; // not used
  double* real_data = new double[4]; 
  for(int i=0; i<4; ++i) real_data[i] = FrictionCoef[i];
  contact_obj->Add_Enforcement_Model(type, ID, integer_data, real_data);
  
  int number_iterations;
  if(((ConstraintOptionsData && ConstraintOptionsData->lagrangeMult == 0 && ConstraintOptionsData->penalty > 0) ||
     (ConstraintOptionsData == NULL && domain->solInfo().lagrangeMult == 0 && domain->solInfo().penalty > 0)) &&
     !domain->solInfo().default_penalty) {
    number_iterations = 1;
  }
  else {
    number_iterations = TDEnfNumIter;
  }
  error = contact_obj->Set_Number_of_Iterations(number_iterations);
  if(error) {
    std::cerr << "Error in ACME ContactTDEnforcement::Set_Number_of_Iterations: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  } 
  
  double tol = TDEnfConvTol;
  contact_obj->Set_Convergence_Tolerance(tol);
  if(error) {
    std::cerr << "Error in ACME ContactTDEnforcement::Set_Convergence_Tolerance: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  }

  delete [] Enforcement_Data;
  if(integer_data) delete [] integer_data;
  if(real_data) delete [] real_data;
#endif
}

void
MortarHandler::remove_gap(Vector &d)
{
#ifdef USE_ACME
  double *Enforcement_Data = new double[num_entity_keys*num_entity_keys]; // this array stores the kinematic partition factor
  for(int i=0; i<num_entity_keys*num_entity_keys; ++i) Enforcement_Data[i] = 0.0;
  for(int f_key=0; f_key<4; ++f_key) {
    for(int n_key = 4; n_key < 8; ++n_key) {
     Enforcement_Data[f_key*num_entity_keys+n_key] = 1.0; // kinematic partitioning factor
    }
  }

  ContactSearch::ContactErrorCode error;
  ContactGapRemoval *gap_removal_obj = new ContactGapRemoval(Enforcement_Data, search_obj, error);
  if(error) {
    std::cerr << "Error in ACME ContactGapRemoval::ContactGapRemoval: error code = " << error << std::endl;
    for(int i=1; i<=gap_removal_obj->Number_of_Errors(); ++i)
      std::cerr << gap_removal_obj->Error_Message(i) << std::endl;
    exit(error);
  }
  delete [] Enforcement_Data;

  set_node_configuration(1);
  perform_search(1);

  int num_nodes = PtrSlaveEntity->GetnVertices() + PtrMasterEntity->GetnVertices();
  int max_iterations = 100;
  double trivial_gap = 1.0e-10;
  double* displ_cor = new double[3*num_nodes];
  error = gap_removal_obj->Compute_Gap_Removal(max_iterations, trivial_gap, displ_cor);
  if(error) {
    std::cerr << "Error in ACME ContactGapRemoval::Compute_Gap_Removal: error code = " << error << std::endl;
    for(int i=1; i<=gap_removal_obj->Number_of_Errors(); ++i)
      std::cerr << gap_removal_obj->Error_Message(i) << std::endl;
    exit(error);
  }
  // assemble contact force
  for(int i=0; i<3*num_nodes; ++i) if(dofmap[i] > -1) d[dofmap[i]] += displ_cor[i];
  delete [] displ_cor;
  delete gap_removal_obj;
#endif
}

void
MortarHandler::make_nodal_mass(SparseMatrix *M, ConstrainedDSA *c_dsa, SparseMatrix *Mcc)
{
  int num_nodes = PtrSlaveEntity->GetnVertices() + PtrMasterEntity->GetnVertices();
  mass = new double[num_nodes];
  for(int i=0; i<num_nodes; ++i) mass[i] = 0.0;
  dofmap = new int[3*num_nodes];
  for(int i=0; i<3*num_nodes; ++i) dofmap[i] = -1;
  int inode = 0; int dofs[3]; int count; double m;
  for(int i = 0; i<PtrSlaveEntity->GetnVertices(); ++i, ++inode) {
    int glNode = PtrSlaveEntity->GetGlVertexId(i);
    count = 0; m = 0.0;
    c_dsa->number(glNode,DofSet::XYZdisp,dofs);
    for(int j=0; j<3; ++j) {
      if(dofs[j] > -1) {
        m += M->diag(dofs[j]);
        dofmap[3*inode+j] = dofs[j];
        count++;
      }
    }
    if(count > 0) mass[inode] += m/double(count);
  }  
  for(int i = 0; i<PtrMasterEntity->GetnVertices(); ++i, ++inode) {
    int glNode = PtrMasterEntity->GetGlVertexId(i);
    count = 0; m = 0.0;
    c_dsa->number(glNode,DofSet::XYZdisp,dofs);
    for(int j=0; j<3; ++j) {
      if(dofs[j] > -1) {
        m += M->diag(dofs[j]);
        dofmap[3*inode+j] = dofs[j];
        count++;
      }
    }
    if(count > 0) mass[inode] += m/double(count);
  }
  for(int i=0; i<num_nodes; ++i) if(mass[i] == 0.0) mass[i] = 1.0; // TODO: need Mcc for for constrained nodes' mass
}

void
MortarHandler::make_nodal_mass(SubDOp *M, SubDomain **sd)
{
  // multiple domain version
  int numSub = M->getNumSub();
  int num_nodes = PtrSlaveEntity->GetnVertices() + PtrMasterEntity->GetnVertices();
  mass = new double[num_nodes];
  for(int i=0; i<num_nodes; ++i) mass[i] = 0.0;
  dofmap = new int[numSub*3*num_nodes];
  for(int i=0; i<numSub*3*num_nodes; ++i) dofmap[i] = -1;
  share = new int[3*num_nodes];
  for(int i=0; i<3*num_nodes; ++i) share[i] = 0;

  int k = 0; int dofs[3]; int count; double m;
  for(int s = 0; s < numSub; ++s) {
    ConstrainedDSA *c_dsa = sd[s]->getCDSA();
    int inode = 0;
    for(int i = 0; i<PtrSlaveEntity->GetnVertices(); ++i, ++inode) {
      int glNode = PtrSlaveEntity->GetGlVertexId(i);
      int locNode = sd[s]->globalToLocal(glNode);
      if(locNode > -1) {
        count = 0; m = 0.0;
        c_dsa->number(locNode,DofSet::XYZdisp,dofs);
        for(int j=0; j<3; ++j) {
          share[3*inode+j] += 1;
          if(dofs[j] > -1) {
            m += (*M)[s]->diag(dofs[j]);
            dofmap[3*(s*num_nodes+inode)+j] = dofs[j];
            count++;
          }
        }
        if(count > 0) mass[inode] += m/double(count);
      }
    }
    for(int i = 0; i<PtrMasterEntity->GetnVertices(); ++i, ++inode) {
      int glNode = PtrMasterEntity->GetGlVertexId(i);
      int locNode = sd[s]->globalToLocal(glNode);
      if(locNode > -1) {
        count = 0; m = 0.0;
        c_dsa->number(locNode,DofSet::XYZdisp,dofs);
        for(int j=0; j<3; ++j) {
          share[3*inode+j] += 1;
          if(dofs[j] > -1) { 
            m += (*M)[s]->diag(dofs[j]); 
            dofmap[3*(s*num_nodes+inode)+j] = dofs[j];
            count++;
          }
        }
        if(count > 0) mass[inode] += m/double(count);
      }
    }
  }
#ifdef DISTRIBUTED
  if(DIST_ACME == 0) {
    structCom->globalSum(num_nodes, mass);
    structCom->globalSum(3*num_nodes, share);
  }
  else if(DIST_ACME == 1) {
    structCom->reduce(num_nodes, mass, 0);
    structCom->globalSum(3*num_nodes, share);
    if(structCom->myID() != 0) { delete [] mass; mass = 0; }
  }
  else if(DIST_ACME == 2) {
    FSCommunicator *communicator = new FSCommunicator(structCom);
    int myCPU = communicator->cpuNum();
    int numCPUs = communicator->size();
    int *_ptr = new int[numCPUs+1]; for(int i=0;i<numCPUs+1;++i) _ptr[i] = i;
    int *_tgt = new int[numCPUs]; for(int i=0;i<numCPUs;++i) _tgt[i] = i;
    Connectivity *cpuToSelf = new Connectivity(numCPUs,_ptr,_tgt);
    // comm pattern for shared nodal mass
    FSCommPattern<double> *pat1 = new FSCommPattern<double>(communicator, cpuToSelf, myCPU, FSCommPattern<double>::CopyOnSend);
    for(int i = 0; i < num_comm_partners; ++i) pat1->setLen(myCPU, comm_proc_id[i], number_nodes_to_partner[i]);
    pat1->finalize();
    // comm pattern for shared dof weighting
    FSCommPattern<int> *pat2 = new FSCommPattern<int>(communicator, cpuToSelf, myCPU, FSCommPattern<int>::CopyOnSend);
    for(int i = 0; i < num_comm_partners; ++i) pat2->setLen(myCPU, comm_proc_id[i], number_nodes_to_partner[i]*3);
    pat2->finalize();

    // send and receive the shared nodal mass
    k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<double> info = pat1->getSendBuffer(communicator->cpuNum(), comm_proc_id[i]);
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        info.data[j] = mass[comm_node[k]-1];
      }
    }
    pat1->exchange();
    k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<double> info = pat1->recData(comm_proc_id[i],communicator->cpuNum());
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        mass[comm_node[k]-1] += info.data[j];
      }
    }

    // send and receive the shared dof weighting
    k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<int> info = pat2->getSendBuffer(communicator->cpuNum(), comm_proc_id[i]);
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        info.data[3*j+0] = share[3*(comm_node[k]-1)+0];
        info.data[3*j+1] = share[3*(comm_node[k]-1)+1];
        info.data[3*j+2] = share[3*(comm_node[k]-1)+2];
      }
    }
    pat2->exchange();
    k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<int> info = pat2->recData(comm_proc_id[i],communicator->cpuNum());
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        share[3*(comm_node[k]-1)+0] += info.data[3*j+0];
        share[3*(comm_node[k]-1)+1] += info.data[3*j+1];
        share[3*(comm_node[k]-1)+2] += info.data[3*j+2];
      }
    }
    
    delete cpuToSelf; delete pat1; delete pat2;
  }

  if(!(DIST_ACME == 1 && structCom->myID() != 0))
#endif
  for(int i=0; i<num_nodes; ++i) if(mass[i] == 0.0) mass[i] = 1.0; // TODO: need Mcc for constrained nodes' mass
  //this should not be necessary: for(int i=0; i<3*num_nodes; ++i) if(share[i] == 0) share[i] = 1;
}

void
MortarHandler::make_share(int numSub, SubDomain **sd)
{
  // multiple domain version
  int num_nodes = PtrSlaveEntity->GetnVertices() + PtrMasterEntity->GetnVertices();
  if(share) delete [] share;
  share = new int[3*num_nodes];
  for(int i=0; i<3*num_nodes; ++i) share[i] = 0;

  int k = 0; int dofs[3]; int count; double m;
  for(int s = 0; s < numSub; ++s) {
    //DofSetArray *dsa = sd[s]->getDSA();
    int inode = 0;
    for(int i = 0; i<PtrSlaveEntity->GetnVertices(); ++i, ++inode) {
      int glNode = PtrSlaveEntity->GetGlVertexId(i);
      int locNode = sd[s]->globalToLocal(glNode);
      if(locNode > -1) {
        //dsa->number(locNode,DofSet::XYZdisp,dofs);
        for(int j=0; j<3; ++j) {
          //if(dofs[j] > -1)
            share[3*inode+j] += 1;
        }
      }
    }
    for(int i = 0; i<PtrMasterEntity->GetnVertices(); ++i, ++inode) {
      int glNode = PtrMasterEntity->GetGlVertexId(i);
      int locNode = sd[s]->globalToLocal(glNode);
      if(locNode > -1) {
        //dsa->number(locNode,DofSet::XYZdisp,dofs);
        for(int j=0; j<3; ++j) {
          //if(dofs[j] > -1)
            share[3*inode+j] += 1;
        }
      }
    }
  }
#ifdef DISTRIBUTED
  if(DIST_ACME != 2) {
    structCom->globalSum(3*num_nodes, share);
  }
  else {
    FSCommunicator *communicator = new FSCommunicator(structCom);
    int myCPU = communicator->cpuNum();
    int numCPUs = communicator->size();
    int *_ptr = new int[numCPUs+1]; for(int i=0;i<numCPUs+1;++i) _ptr[i] = i;
    int *_tgt = new int[numCPUs]; for(int i=0;i<numCPUs;++i) _tgt[i] = i;
    Connectivity *cpuToSelf = new Connectivity(numCPUs,_ptr,_tgt);
    // comm pattern for shared dof weighting
    FSCommPattern<int> *pat2 = new FSCommPattern<int>(communicator, cpuToSelf, myCPU, FSCommPattern<int>::CopyOnSend);
    for(int i = 0; i < num_comm_partners; ++i) pat2->setLen(myCPU, comm_proc_id[i], number_nodes_to_partner[i]*3);
    pat2->finalize();

    // send and receive the shared dof weighting
    k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<int> info = pat2->getSendBuffer(communicator->cpuNum(), comm_proc_id[i]);
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        info.data[3*j+0] = share[3*(comm_node[k]-1)+0];
        info.data[3*j+1] = share[3*(comm_node[k]-1)+1];
        info.data[3*j+2] = share[3*(comm_node[k]-1)+2];
      }
    }
    pat2->exchange();
    k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<int> info = pat2->recData(comm_proc_id[i],communicator->cpuNum());
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        share[3*(comm_node[k]-1)+0] += info.data[3*j+0];
        share[3*(comm_node[k]-1)+1] += info.data[3*j+1];
        share[3*(comm_node[k]-1)+2] += info.data[3*j+2];
      }
    }

    delete cpuToSelf; delete pat2;
  }
#endif
  //this should not be necessary: for(int i=0; i<3*num_nodes; ++i) if(share[i] == 0) share[i] = 1;
}

void
MortarHandler::make_kinematic_partitioning(Elemset &packedEset, Connectivity *nodeToElem)
{
  if(MortarType == STD) return;
  int num_nodes = PtrSlaveEntity->GetnVertices() + PtrMasterEntity->GetnVertices();
  density = new double[num_nodes];
  wavespeed = new double[num_nodes];
  for(int i=0; i<num_nodes; ++i) { density[i] = 0.0; wavespeed[i] = 0.0; }

  int inode = 0;
  for(int i = 0; i<PtrSlaveEntity->GetnVertices(); ++i, ++inode) {
    int glNode = PtrSlaveEntity->GetGlVertexId(i);
    for(int j = 0; j < nodeToElem->num(glNode); ++j) {
      int k = (*nodeToElem)[glNode][j];
      Element *ele = packedEset[k];
      StructProp *sprop = ele->getProperty();
      density[inode] += sprop->rho/(nodeToElem->num(glNode));
      wavespeed[inode] += sqrt((sprop->E)/(sprop->rho))/(nodeToElem->num(glNode));
    }
  }
  for(int i = 0; i<PtrMasterEntity->GetnVertices(); ++i, ++inode) {
    int glNode = PtrMasterEntity->GetGlVertexId(i);
    for(int j = 0; j < nodeToElem->num(glNode); ++j) {
      int k = (*nodeToElem)[glNode][j];
      Element *ele = packedEset[k];
      StructProp *sprop = ele->getProperty();
      density[inode] += sprop->rho/(nodeToElem->num(glNode));
      wavespeed[inode] += sqrt((sprop->E)/(sprop->rho))/(nodeToElem->num(glNode));
    }
  }
}

void
MortarHandler::make_kinematic_partitioning(int numSub, SubDomain **sd)
{
  // multiple domain version
  if(MortarType == STD) return;
  int num_nodes = PtrSlaveEntity->GetnVertices() + PtrMasterEntity->GetnVertices();
  density = new double[num_nodes];
  wavespeed = new double[num_nodes];
  int *count = new int[num_nodes];
  for(int i=0; i<num_nodes; ++i) { density[i] = 0.0; wavespeed[i] = 0.0; count[i] = 0; }

  for(int s = 0; s < numSub; ++s) {
    Elemset &packedEset =  sd[s]->getElementSet();
    Connectivity elemToNode(packedEset.asSet());
    Connectivity *nodeToElem = elemToNode.alloc_reverse();
    int inode = 0;
    for(int i = 0; i<PtrSlaveEntity->GetnVertices(); ++i, ++inode) {
      int glNode = PtrSlaveEntity->GetGlVertexId(i);
      int locNode = sd[s]->globalToLocal(glNode);
      if(locNode > -1) {
        for(int j = 0; j < nodeToElem->num(locNode); ++j) {
          int k = (*nodeToElem)[locNode][j];
          Element *ele = packedEset[k];
          StructProp *sprop = ele->getProperty();
          count[inode] += 1;
          density[inode] += sprop->rho;
          wavespeed[inode] += sqrt((sprop->E)/(sprop->rho));
        }
      }
    }
    for(int i = 0; i<PtrMasterEntity->GetnVertices(); ++i, ++inode) {
      int glNode = PtrMasterEntity->GetGlVertexId(i);
      int locNode = sd[s]->globalToLocal(glNode);
      if(locNode > -1) {
        for(int j = 0; j < nodeToElem->num(locNode); ++j) {
          int k = (*nodeToElem)[locNode][j];
          Element *ele = packedEset[k];
          StructProp *sprop = ele->getProperty();
          count[inode] += 1;
          density[inode] += sprop->rho;
          wavespeed[inode] += sqrt((sprop->E)/(sprop->rho));
        }
      }
    }
    delete nodeToElem; 
  }

#ifdef DISTRIBUTED
  if(DIST_ACME == 0) {
    structCom->globalSum(num_nodes, density);
    structCom->globalSum(num_nodes, wavespeed);
    structCom->globalSum(num_nodes, count);
  }
  else if(DIST_ACME == 1) {
    structCom->reduce(num_nodes, density, 0);
    structCom->reduce(num_nodes, wavespeed, 0);
    structCom->reduce(num_nodes, count, 0);
    if(structCom->myID() != 0) { delete [] density; density = 0; delete [] wavespeed; wavespeed = 0; }
  }
  else if(DIST_ACME == 2) {
    FSCommunicator *communicator = new FSCommunicator(structCom);
    int myCPU = communicator->cpuNum();
    int numCPUs = communicator->size();
    int *_ptr = new int[numCPUs+1]; for(int i=0;i<numCPUs+1;++i) _ptr[i] = i;
    int *_tgt = new int[numCPUs]; for(int i=0;i<numCPUs;++i) _tgt[i] = i;
    Connectivity *cpuToSelf = new Connectivity(numCPUs,_ptr,_tgt);
    // comm pattern for shared interface force
    FSCommPattern<double> *pat3 = new FSCommPattern<double>(communicator, cpuToSelf, myCPU, FSCommPattern<double>::CopyOnSend);
    for(int i = 0; i < num_comm_partners; ++i) pat3->setLen(myCPU, comm_proc_id[i], number_nodes_to_partner[i]*3);
    pat3->finalize();

    int k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<double> info = pat3->getSendBuffer(communicator->cpuNum(), comm_proc_id[i]);
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        info.data[3*j+0] = density[comm_node[k]-1];
        info.data[3*j+1] = wavespeed[comm_node[k]-1];
        info.data[3*j+2] = double(count[comm_node[k]-1]);
      }
    }
    pat3->exchange();
    k = 0;
    for(int i = 0; i < num_comm_partners; ++i) {
      FSSubRecInfo<double> info = pat3->recData(comm_proc_id[i],communicator->cpuNum());
      for(int j = 0; j < number_nodes_to_partner[i]; ++j,++k) {
        density[comm_node[k]-1] += info.data[3*j+0];
        wavespeed[comm_node[k]-1] += info.data[3*j+1];
        count[comm_node[k]-1] += int(info.data[3*j+2]);
      }
    }
  
    delete cpuToSelf; delete pat3;
  }

  if(!(DIST_ACME == 1 && structCom->myID() != 0))
#endif
  for(int i=0; i<num_nodes; ++i) { density[i] /= double(count[i]); wavespeed[i] /= double(count[i]); }
  delete [] count;
}

void
MortarHandler::compute_td_contact_force(double dt_old, double dt, Vector &f)
{
#ifdef USE_ACME
  ContactSearch::ContactErrorCode error;

  if(!domain->solInfo().default_penalty) {
  // override the default ACME penalty parameter with the value set in the AERO-S input file
  if((ConstraintOptionsData && ConstraintOptionsData->lagrangeMult == 0 && ConstraintOptionsData->penalty != 0) ||
     (ConstraintOptionsData == NULL && domain->solInfo().lagrangeMult == 0 && domain->solInfo().penalty != 0)) {
    double dt2 = 1.0/(0.5*(dt+dt_old)*dt);
    double penalty = (ConstraintOptionsData) ? ConstraintOptionsData->penalty : domain->solInfo().penalty;
    // if the penalty parameter is less than zero, then by convention it's absolute value is the scale factor
    double penalty_scale = (penalty > 0) ? 2*penalty/dt2 : -penalty;
    error = static_cast<ContactTDEnfPenalty*>(contact_obj)->Set_Penalty_Scale(penalty_scale);
    if(error) {
      std::cerr << "Error in ACME ContactTDEnfPenalty::Set_Penalty_Scale: error code = " << error << std::endl;
      for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
        std::cerr << contact_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }}

  int nACMENodes  = nMasterNodes + nSlaveNodes; // !! Assume NO COMMON nodes !!
  double *force = new double[nACMENodes*3]; // force vectors
  error = contact_obj->Compute_Contact_Force(dt_old, dt, mass, density, wavespeed, force);
  if(error) {
    std::cerr << "Error in ACME ContactTDEnforcement::Compute_Contact_Force: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  }
  // transform to DOF_FRM 
  if(!domain->solInfo().basicDofCoords) {
    for(int i=0; i<nACMENodes; ++i) domain->transformVector(&force[3*i], node_global_ids[2*i+1]-1, false);
  }
  // assemble contact force
  for(int i=0; i<3*nACMENodes; ++i) if(dofmap[i] > -1) f[dofmap[i]] += force[i];
  delete [] force;
#endif
}

void
MortarHandler::compute_td_contact_force(double dt_old, double dt, DistrVector &f)
{
  // multiple domain version
#ifdef USE_ACME
  ContactSearch::ContactErrorCode error;

  // override the ACME default penalty parameter with the value set in the AERO-S input file
  if(!domain->solInfo().default_penalty) {
  if((ConstraintOptionsData && ConstraintOptionsData->lagrangeMult == 0 && ConstraintOptionsData->penalty != 0) ||
     (ConstraintOptionsData == NULL && domain->solInfo().lagrangeMult == 0 && domain->solInfo().penalty != 0)) {
    double dt2 = 1.0/(0.5*(dt+dt_old)*dt);
    double penalty = (ConstraintOptionsData) ? ConstraintOptionsData->penalty : domain->solInfo().penalty;
    // if the penalty parameter is less than zero, then by convention it's absolute value is the scale factor
    double penalty_scale = (penalty > 0) ? 2*penalty/dt2 : -penalty;
    error = static_cast<ContactTDEnfPenalty*>(contact_obj)->Set_Penalty_Scale(penalty_scale);
    if(error) {
      std::cerr << "Error in ACME ContactTDEnfPenalty::Set_Penalty_Scale: error code = " << error << std::endl;
      for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
        std::cerr << contact_obj->Error_Message(i) << std::endl;
      exit(error);
    }
  }}

  int nACMENodes = nMasterNodes + nSlaveNodes; 
  double *force = new double[nACMENodes*3]; 
  error = contact_obj->Compute_Contact_Force(dt_old, dt, mass, density, wavespeed, force);
  if(error) {
    std::cerr << "Error in ACME ContactTDEnforcement::Compute_Contact_Force: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  }

  // transform to DOF_FRM 
  if(!domain->solInfo().basicDofCoords) {
    for(int i=0; i<nACMENodes; ++i) domain->transformVector(&force[3*i], node_global_ids[2*i+1]-1, false);
  }

#ifdef DISTRIBUTED
  if(DIST_ACME == 1) {
    if(structCom->myID() != 0) {
      delete [] force;
      nACMENodes = PtrMasterEntity->GetnVertices() + PtrSlaveEntity->GetnVertices();
      force = new double[nACMENodes*3];
    }
    structCom->broadcast(3*nACMENodes, force, 0);
  }
#endif

  // assemble contact force
  for(int j = 0; j < f.num(); ++j) {
    for(int i=0; i<3*nACMENodes; ++i) {
      int k = j*3*nACMENodes+i;
      if(dofmap[k] > -1) f.subData(j)[dofmap[k]] += (force[i]/double(share[i]));
    }
  }
  delete [] force;
#endif
}

void
MortarHandler::get_plot_variable(int plot_var, double *all_data)
{
#ifdef USE_ACME
  if(!data) { data = new double[nMasterNodes + nSlaveNodes]; return; } // first output there is no td enforcement 

  ContactSearch::ContactErrorCode error = contact_obj->Get_Plot_Variable((ContactTDEnforcement::Contact_TDEnf_Plot_Vars) plot_var, data);
  if(error) {
    std::cerr << "Error in ACME ContactTDEnforcement::Get_Plot_Variable: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  }
  for(int i=0; i<nMasterNodes+nSlaveNodes; ++i) all_data[node_global_ids[2*i+1]-1] = data[i];
#endif
}

void
MortarHandler::get_plot_variable(int plot_var, double **plot_data, int numSub, gsl::span<SubDomain *> sd)
{
#ifdef USE_ACME
  if(!data) { data = new double[nMasterNodes + nSlaveNodes]; return; } // first output there is no td enforcement 

  ContactSearch::ContactErrorCode error = contact_obj->Get_Plot_Variable((ContactTDEnforcement::Contact_TDEnf_Plot_Vars) plot_var, data);
  if(error) {
    std::cerr << "Error in ACME ContactTDEnforcement::Get_Plot_Variable: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  }
  for(int i=0; i<nMasterNodes+nSlaveNodes; ++i) {
    for(int s=0; s<numSub; ++s) {
      int n;
      if((n = sd[s]->globalToLocal(node_global_ids[2*i+1]-1)) > -1) plot_data[s][n] = data[i];
    }
  }
#endif
}


void
MortarHandler::get_global_variable(int var, double &value)
{
#ifdef USE_ACME
  ContactSearch::ContactErrorCode error = contact_obj->Get_Global_Variable((ContactTDEnforcement::Contact_TDEnf_Global_Vars) var, value);
  if(error) {
    std::cerr << "Error in ACME ContactTDEnforcement::Get_Global_Variable: error code = " << error << std::endl;
    for(int i=1; i<=contact_obj->Number_of_Errors(); ++i)
      std::cerr << contact_obj->Error_Message(i) << std::endl;
    exit(error);
  }
#endif
}

