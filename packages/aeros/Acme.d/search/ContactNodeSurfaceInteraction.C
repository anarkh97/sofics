// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#include "ContactTopologyEntityList.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactTopology.h"
#include <cstring>
#include <new>

ContactNodeSurfaceInteraction::ContactNodeSurfaceInteraction( 
   InteractionSource source,
   ContactNode<Real>* Node,
   ContactAnalyticSurface* Surface,
   int node_entity_key,
   Real* physical_face_normal,
   Real  Gap,
   Real* cpoint,
   Real* snorm,
   Real* pb_dir,
   Real  ttc,
   int   Location,
   bool  IS_TIED,
   bool  IS_INFSLIP,
   bool  IS_GLUED,
   bool  IS_TRACKED ):
     ContactNodeEntityInteraction(IS_TIED, IS_INFSLIP, IS_GLUED, IS_TRACKED, Node, Surface, Location, NODE_SURFACE_INTERACTION, CT_NSI),
     surface_id(Surface->ProcArrayIndex())
{
  PRECONDITION( Node );
  PRECONDITION( Surface );
  PRECONDITION( physical_face_normal );
  PRECONDITION( cpoint );
  PRECONDITION( snorm );
  PRECONDITION( pb_dir );
  
  Set_NodeEntityData();
  Set_Entity_Data();
  
  // Currently, we are assuming all the gap occurred this step
  Scalar_Var(GAP_CUR) = Gap;
  Scalar_Var(GAP_OLD) = 0.0;
  Scalar_Var(ContactNodeSurfaceInteraction::SOURCE)  = source;
  Real* coordinates = Vector_Var(COORDINATES);
  coordinates[0] = cpoint[0];
  coordinates[1] = cpoint[1];
  coordinates[2] = cpoint[2];
  Real* contact_point = Vector_Var(CONTACT_POINT);
  contact_point[0] = cpoint[0];
  contact_point[1] = cpoint[1];
  contact_point[2] = cpoint[2];
  Real* surface_normal = Vector_Var(NORMAL_DIR);
  surface_normal[0] = snorm[0];
  surface_normal[1] = snorm[1];
  surface_normal[2] = snorm[2];
  Real* pf_normal = Vector_Var(PHYSICAL_FACE_NORMAL);
  pf_normal[0] = physical_face_normal[0];
  pf_normal[1] = physical_face_normal[1];
  pf_normal[2] = physical_face_normal[2];
  Real* pushback_dir = Vector_Var(PUSHBACK_DIR);
  pushback_dir[0] = pb_dir[0];
  pushback_dir[1] = pb_dir[1];
  pushback_dir[2] = pb_dir[2];
  Scalar_Var(TIME_TO_CONTACT) = ttc;
  Scalar_Var(NODE_ENTITY_KEY) = node_entity_key;
  Scalar_Var(NODE_AREA)       = -1.0;
}

ContactNodeSurfaceInteraction::ContactNodeSurfaceInteraction(const ContactNodeSurfaceInteraction &cnsi):
    ContactNodeEntityInteraction(cnsi.is_tied, cnsi.is_infSlip, cnsi.is_glued, cnsi.is_tracked, cnsi.node, cnsi.entity, cnsi.location, NODE_SURFACE_INTERACTION, CT_NSI),
    surface_id(cnsi.surface_id)
{
  //copy the data
  std::memcpy(DataArray, cnsi.DataArray, DataArray_Length()*sizeof(Real));
}

ContactNodeSurfaceInteraction::ContactNodeSurfaceInteraction():
    ContactNodeEntityInteraction(false, false, false, false, NULL, NULL, 0, NODE_SURFACE_INTERACTION, CT_NSI),
    surface_id(-1)
{
  Scalar_Var(ContactNodeSurfaceInteraction::SOURCE) = UNKNOWN_SOURCE;
  Scalar_Var(GAP_CUR) = 0.0;
  Scalar_Var(GAP_OLD) = 0.0;
  Real* coordinates = Vector_Var(COORDINATES);
  coordinates[0] = 0.0;
  coordinates[1] = 0.0;
  coordinates[2] = 0.0;
  Real* contact_point = Vector_Var(CONTACT_POINT);
  contact_point[0] = 0.0;
  contact_point[1] = 0.0;
  contact_point[2] = 0.0;
  Real* surface_normal = Vector_Var(NORMAL_DIR);
  surface_normal[0] = 0.0;
  surface_normal[1] = 0.0;
  surface_normal[2] = 0.0;
  Real* pf_normal = Vector_Var(PHYSICAL_FACE_NORMAL);
  pf_normal[0] = 0.0;
  pf_normal[1] = 0.0;
  pf_normal[2] = 0.0;
  Real* pushback_dir = Vector_Var(PUSHBACK_DIR);
  pushback_dir[0] = 0.0;
  pushback_dir[1] = 0.0;
  pushback_dir[2] = 0.0;
  Scalar_Var(TIME_TO_CONTACT) = -BIGNUM;
}

ContactNodeSurfaceInteraction*
ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction(
                                     ContactFixedSizeAllocator& alloc,
                                     InteractionSource interaction_source,
                                     ContactNode<Real>* Node,
                                     ContactAnalyticSurface* Surface,
				     int   node_entity_key,
                                     Real* physical_face_normal,
                                     Real  Gap,
				     Real* cpoint,
				     Real* snorm,
                                     Real* pbdir,
                                     Real  ttc,
                                     int   loc,
				     bool  IS_TIED,
                                     bool  IS_INFSLIP,
                                     bool  IS_GLUED,
				     bool  IS_TRACKED )
{
  return new (alloc.New_Frag())
    ContactNodeSurfaceInteraction( interaction_source, Node, Surface, node_entity_key,
				   physical_face_normal, Gap, cpoint, snorm, pbdir,
                                   ttc, loc, IS_TIED, IS_INFSLIP, IS_GLUED, IS_TRACKED );
}                                                             

ContactNodeSurfaceInteraction*
ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction(
				        ContactFixedSizeAllocator& alloc )
{
  return new (alloc.New_Frag())
    ContactNodeSurfaceInteraction( );
}
 
ContactNodeSurfaceInteraction* 
ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactNodeSurfaceInteraction& cnsi )
{
  return new (alloc.New_Frag())ContactNodeSurfaceInteraction( cnsi );
}

ContactNodeEntityInteraction *ContactNodeSurfaceInteraction::New_Instance(ContactFixedSizeAllocator* allocators){
  return new_ContactNodeSurfaceInteraction(allocators[ContactSearch::ALLOC_ContactNodeSurfaceInteraction], *this);
}

void 
ContactNodeSurfaceInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactNodeSurfaceInteraction),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactQNodeSurfaceInteraction allocator" );
}


ContactNodeSurfaceInteraction::~ContactNodeSurfaceInteraction() {}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------

int ContactNodeSurfaceInteraction::Size()
{
  return(ContactInteractionEntity<Real>::Size()+
         2*sizeof(entity_data)+
         6*sizeof(int)+
         (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real));
}

char* ContactNodeSurfaceInteraction::Pack( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Pack( buffer );
  buff += ContactInteractionEntity<Real>::Size();
  
  int  cnt     = 0;
  int* i_buf   = reinterpret_cast<int*>(buff);
  cnt         += PackEntityData(&node_entity_data,   &i_buf[cnt]);
  cnt         += PackEntityData(&entity_entity_data, &i_buf[cnt]);
  i_buf[cnt++] = Surface()->ProcArrayIndex();
  i_buf[cnt++] = is_tied;
  i_buf[cnt++] = is_infSlip;
  i_buf[cnt++] = is_glued;
  i_buf[cnt++] = is_tracked;
  i_buf[cnt++] = location;
  
  buff += cnt*sizeof(int);
  cnt   = (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real);
  std::memcpy( buff,DataArray,cnt);
  buff += cnt;
  return buff;
}

char* ContactNodeSurfaceInteraction::Unpack( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Unpack( buffer );
  buff += ContactInteractionEntity<Real>::Size();
  
  int cnt    = 0;
  int* i_buf = reinterpret_cast<int*>(buff);
  cnt       += UnPackEntityData(&node_entity_data,   &i_buf[cnt]);
  cnt       += UnPackEntityData(&entity_entity_data, &i_buf[cnt]);
  surface_id = i_buf[cnt++];
  is_tied    = i_buf[cnt++];
  is_infSlip = i_buf[cnt++];
  is_glued   = i_buf[cnt++];
  is_tracked = i_buf[cnt++];
  location   = i_buf[cnt++];
  
  buff += cnt*sizeof(int);
  cnt   = (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real);
  std::memcpy( DataArray,buff,cnt);
  buff += cnt;
  return buff;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------

  #define DEFINE_DATA_TO_SEND int vector_data   = 1;\
                              int scalar_data   = 0;\
                              int vector_offset = PHYSICAL_FACE_NORMAL;\
                              int scalar_offset = GAP_CUR;\
                              if (is_tied||is_infSlip||is_glued) {\
                                scalar_data   = 7;\
                                scalar_offset = GAP_CUR;\
                                vector_data   = 5;\
                                vector_offset = NORMAL_DIR;\
                              } else if (is_tracked) {\
                                scalar_data   = 2;\
                                scalar_offset = SOURCE;\
                                vector_data   = 4;\
                                vector_offset = NORMAL_DIR;\
                              }

inline int ContactNodeSurfaceInteraction::Size_ForSecondary()
{
  DEFINE_DATA_TO_SEND;
  return(ContactInteractionEntity<Real>::Size_ForSecondary()+
         2*sizeof(int)+(scalar_data+3*vector_data)*sizeof(Real));
}

inline char* ContactNodeSurfaceInteraction::Pack_ForSecondary( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Pack_ForSecondary( buffer );
  buff += ContactInteractionEntity<Real>::Size_ForSecondary();
  
  int  cnt     = 0;
  int* i_buf   = reinterpret_cast<int*> (buff);
  
  i_buf[cnt++] = Surface()->ProcArrayIndex();
  i_buf[cnt++] = is_tied<<2 | is_infSlip<<3 | is_glued<<4 | is_tracked<<5;
  
  DEFINE_DATA_TO_SEND;
  
  buff += cnt*sizeof(int);
  if (vector_data>0) {
    cnt   = 3*vector_data*sizeof(Real);
    std::memcpy( buff,&DataArray[NUMBER_SCALAR_VARS+3*vector_offset],cnt);
    buff += cnt;
  }
  if (scalar_data>0) {
    cnt   = scalar_data*sizeof(Real);
    std::memcpy( buff,&DataArray[scalar_offset],cnt);
    buff += cnt;
  }
  return buff;
}

inline char* ContactNodeSurfaceInteraction::Unpack_ForSecondary( char* buffer )
{
  char* buff = buffer;
  buff += ContactInteractionEntity<Real>::Size_ForSecondary();
  
  int cnt    = 0;
  int* i_buf = reinterpret_cast<int*> (buff);
  
  surface_id =  i_buf[cnt++];
  is_tied    = (i_buf[cnt] & 0x4 )>>2;
  is_infSlip = (i_buf[cnt] & 0x8 )>>3;
  is_glued   = (i_buf[cnt] & 0x10)>>4;
  is_tracked = (i_buf[cnt] & 0x20)>>5;
  cnt++;
  
  DEFINE_DATA_TO_SEND;
  
  buff += cnt*sizeof(int);
  if (vector_data>0) {
    cnt   = 3*vector_data*sizeof(Real);
    std::memcpy( &DataArray[NUMBER_SCALAR_VARS+3*vector_offset],buff,cnt);
    buff += cnt;
  }
  if (scalar_data>0) {
    cnt   = scalar_data*sizeof(Real);
    std::memcpy( &DataArray[scalar_offset],buff,cnt);
    buff += cnt;
  }
  return buff;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for restart
//--------------------------------------------------------------------

int ContactNodeSurfaceInteraction::Restart_Size()
{
  return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS+
         2*sizeof(entity_data)/sizeof(int)+6;
}

void ContactNodeSurfaceInteraction::Restart_Pack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += PackEntityData(&node_entity_data,   &buf_loc[cnt]);
  cnt += PackEntityData(&entity_entity_data, &buf_loc[cnt]);
  buf_loc[cnt++] = Surface()->ProcArrayIndex();
  buf_loc[cnt++] = is_tied;
  buf_loc[cnt++] = is_infSlip;
  buf_loc[cnt++] = is_glued;
  buf_loc[cnt++] = is_tracked;
  buf_loc[cnt++] = location;

  // Pack the data
  std::memcpy( &buf_loc[cnt], &DataArray, 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );

}

void ContactNodeSurfaceInteraction::Restart_Unpack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += UnPackEntityData(&node_entity_data,   &buf_loc[cnt]);
  cnt += UnPackEntityData(&entity_entity_data, &buf_loc[cnt]);
  surface_id = (int)  buf_loc[cnt++];
  is_tied    = (bool) buf_loc[cnt++];
  is_infSlip = (bool) buf_loc[cnt++];
  is_glued   = (bool) buf_loc[cnt++];
  is_tracked = (bool) buf_loc[cnt++];
  location   = (int ) buf_loc[cnt++];

  // Unpack the Data Array
  std::memcpy( &DataArray, &buf_loc[cnt], 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );
}

void ContactNodeSurfaceInteraction::Get_Face_Normal(Real *return_normal,
                                                    ContactTopology* search_topology) {
  Real* surface_normal = Vector_Var(NORMAL_DIR);
  return_normal[0] = surface_normal[0];
  return_normal[1] = surface_normal[1];
  return_normal[2] = surface_normal[2];
}

void ContactNodeSurfaceInteraction::Get_Avg_Face_Normal(Real *return_normal, 
                                                        ContactTopology* search_topology,
                                                        VariableHandle NODE_COORDS) {
  Get_Face_Normal(return_normal, search_topology);
}

void ContactNodeSurfaceInteraction::Delete_Frag(ContactFixedSizeAllocator* allocators) {
  this->~ContactNodeSurfaceInteraction();
  allocators[ContactSearch::ALLOC_ContactNodeSurfaceInteraction].Delete_Frag(this);
}

int ContactNodeSurfaceInteraction::Get_Entity_Key() {
  return Surface()->ProcArrayIndex();
}

/**
 *  This function ignores the last five variable handles as they are not needed, however, they are part of the
 *  combined interface for a node entity interaction. 
 */
void ContactNodeSurfaceInteraction::Update_Tied_Interaction( VariableHandle   NODE_POSITION,
                                                             VariableHandle   FACE_NORMAL,
                                                             Real   gap0,
                                                             ScratchVariable &NUM_PF,
                                                             ScratchVariable &PF_1,
                                                             ScratchVariable &PF_2,
                                                             ScratchVariable &PF_3){

  PRECONDITION(is_tied || is_infSlip || is_glued);
  
  Real *position = node->Variable(NODE_POSITION);
  Real cpoint[3];
  Real snorm[3];
  Real pbdir[3];
  Real ttc;
  Real gap;
  int  loc;

  bool valid_interaction = false;
  valid_interaction = Surface()->Process(position, gap, cpoint, snorm, pbdir, ttc, loc);
  
  //asserting out here because this interaction is no longer valid. Since this is a tied interaction
  //I don't know what could cause that, but it can't be good.
  PRECONDITION(valid_interaction == true);

  Real *coordinates = Vector_Var(COORDINATES);
  coordinates[0] = cpoint[0];
  coordinates[1] = cpoint[1];
  coordinates[2] = cpoint[2];

  Real *contact_point = Vector_Var(CONTACT_POINT);
  contact_point[0] = cpoint[0];
  contact_point[1] = cpoint[1];
  contact_point[2] = cpoint[2];
  
  Real* pf_normal = Vector_Var(PHYSICAL_FACE_NORMAL);
  pf_normal[0] = 0.0;
  pf_normal[1] = 0.0;
  pf_normal[2] = 0.0;
 
  Real *surface_normal = Vector_Var(NORMAL_DIR);
  surface_normal[0] = snorm[0];
  surface_normal[1] = snorm[1];
  surface_normal[2] = snorm[2];
  
  Real *pushback_dir = Vector_Var(PUSHBACK_DIR);
  pushback_dir[0] = pbdir[0];
  pushback_dir[1] = pbdir[1];
  pushback_dir[2] = pbdir[2];
  
  Scalar_Var(TIME_TO_CONTACT) = ttc;
  
  location = loc;
  
  Scalar_Var(GAP_CUR) = gap;
  Scalar_Var(GAP_OLD) = 0.0;
}

void ContactNodeSurfaceInteraction::Update_Tied_Interaction( VariableHandle   NODE_POSITION,
                                                             VariableHandle   FACE_NORMAL,
                                                             Real   gap0,
                                                             VariableHandle   NODE_NORMAL){

  PRECONDITION(is_tied || is_infSlip || is_glued);
  Real *position = node->Variable(NODE_POSITION);
  Real cpoint[3];
  Real snorm[3];
  Real pbdir[3];
  Real ttc;
  Real gap;
  int  loc;

  bool valid_interaction = false;
  valid_interaction = Surface()->Process(position, gap, cpoint, snorm, pbdir, ttc, loc);
  
  //asserting out here because this interaction is no longer valid. Since this is a tied interaction
  //I don't know what could cause that, but it can't be good.
  PRECONDITION(valid_interaction == true);

  Real *coordinates = Vector_Var(COORDINATES);
  coordinates[0] = cpoint[0];
  coordinates[1] = cpoint[1];
  coordinates[2] = cpoint[2];

  Real *contact_point = Vector_Var(CONTACT_POINT);
  contact_point[0] = cpoint[0];
  contact_point[1] = cpoint[1];
  contact_point[2] = cpoint[2];
  
  Real* pf_normal = Vector_Var(PHYSICAL_FACE_NORMAL);
  pf_normal[0] = 0.0;
  pf_normal[1] = 0.0;
  pf_normal[2] = 0.0;
 
  Real *surface_normal = Vector_Var(NORMAL_DIR);
  surface_normal[0] = snorm[0];
  surface_normal[1] = snorm[1];
  surface_normal[2] = snorm[2];
  
  Real *pushback_dir = Vector_Var(PUSHBACK_DIR);
  pushback_dir[0] = pbdir[0];
  pushback_dir[1] = pbdir[1];
  pushback_dir[2] = pbdir[2];
  
  Scalar_Var(TIME_TO_CONTACT) = ttc;
  
  location = loc;
  
  Scalar_Var(GAP_CUR) = gap;
  Scalar_Var(GAP_OLD) = 0.0;
}

void ContactNodeSurfaceInteraction::Update_Tied_Interaction( VariableHandle   NODE_POSITION,
                                                             VariableHandle   FACE_NORMAL,
                                                             Real   gap0,
                                                             Real* pfn){

  PRECONDITION(is_tied || is_infSlip || is_glued);
  Real *position = node->Variable(NODE_POSITION);
  Real cpoint[3];
  Real snorm[3];
  Real pbdir[3];
  Real ttc;
  Real gap;

  bool valid_interaction = false;
  valid_interaction = Surface()->Process(position, gap, cpoint, snorm, pbdir, ttc, location);
  
  //asserting out here because this interaction is no longer valid. Since this is a tied interaction
  //I don't know what could cause that, but it can't be good.
  PRECONDITION(valid_interaction == true);

  Real *coordinates = Vector_Var(COORDINATES);
  coordinates[0] = cpoint[0];
  coordinates[1] = cpoint[1];
  coordinates[2] = cpoint[2];

  Real *contact_point = Vector_Var(CONTACT_POINT);
  contact_point[0] = cpoint[0];
  contact_point[1] = cpoint[1];
  contact_point[2] = cpoint[2];
  
  Real* pf_normal = Vector_Var(PHYSICAL_FACE_NORMAL);
  pf_normal[0] = pfn[0];
  pf_normal[1] = pfn[1];
  pf_normal[2] = pfn[2];
 
  Real *surface_normal = Vector_Var(NORMAL_DIR);
  surface_normal[0] = snorm[0];
  surface_normal[1] = snorm[1];
  surface_normal[2] = snorm[2];
  
  Real *pushback_dir = Vector_Var(PUSHBACK_DIR);
  pushback_dir[0] = pbdir[0];
  pushback_dir[1] = pbdir[1];
  pushback_dir[2] = pbdir[2];
  
  Scalar_Var(TIME_TO_CONTACT) = ttc;
  
  Scalar_Var(GAP_CUR) = gap;
  Scalar_Var(GAP_OLD) = 0.0;
}

void ContactNodeSurfaceInteraction::Update_Tied_Interaction( VariableHandle POSITION,
                                                             Real* buffer,
                                                             Real  gap0,
                                                             Real* pfn){

  PRECONDITION(is_tied || is_infSlip || is_glued);
  PRECONDITION(0); //should never get called for node/surface interactions
}

void ContactNodeSurfaceInteraction::Connect_Entity( ContactTopology* topology)
{
  Connect_Surface(topology);
}

/*
void ContactNodeFaceInteraction::Get_Entity_Sum(VariableHandle source_var, Real *target_data) {
  //
  //  Don't know how to do the entity sum on a surface, so just zero out the
  //  target data array
  //
  for(j=0 ; j<NDIM ; ++j){ target_data[j] = 0.0; }
}
*/
