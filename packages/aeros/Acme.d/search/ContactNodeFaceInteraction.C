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


#include "allocators.h"
#include "ContactUtilities.h"
#include "CString.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactTopologyEntityHash.h"
#include "ContactTopologyEntityList.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactTopology.h"
#include "ContactNodeBlock.h"
#include "ContactFaceBlock.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactQuadFaceQ9.h"
#include "ContactTriFaceL3.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineFaceQ3.h"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <new>

using acme::Dot;
using acme::Normalize;
using acme::Magnitude;
using acme::Copy;

ContactNodeFaceInteraction::ContactNodeFaceInteraction( ) 
  : ContactNodeEntityInteraction(false, false, false, false, NULL, NULL, 0, NODE_FACE_INTERACTION, CT_NFI)
{
  Scalar_Var(SOURCE)          = UNKNOWN_SOURCE;
  Scalar_Var(TIME_TO_CONTACT) = -BIGNUM;
  location                    = 0;
  number_shared_faces         = 0;
}

ContactNodeFaceInteraction::ContactNodeFaceInteraction( ContactNode<Real>* Node ) 
  : ContactNodeEntityInteraction(false, false, false, false, Node, NULL, 0, NODE_FACE_INTERACTION, CT_NFI)
{
  Set_NodeEntityData();
  Scalar_Var(SOURCE)          = UNKNOWN_SOURCE;
  Scalar_Var(TIME_TO_CONTACT) = -BIGNUM;
  number_shared_faces         = 0;
}

ContactNodeFaceInteraction::ContactNodeFaceInteraction( 
                                                  InteractionSource Source,
                                                  ContactNode<Real>* Node,
                                                  ContactFace<Real>* Face_,
                                                  Real* Data,
                                                  int node_data_key,
                                                  Real* Physical_Face_Normal,
                                                  bool IS_TIED,
                                                  bool IS_INFSLIP,
						  bool IS_TRACKED,
                                                  bool IS_GLUED,
                                                  int  include_neighbors,
                                                  const VariableHandle &POSITION )
  : ContactNodeEntityInteraction(IS_TIED, IS_INFSLIP, IS_TRACKED, IS_GLUED, Node, Face_, (int)Data[11], NODE_FACE_INTERACTION, CT_NFI)
{
  PRECONDITION( Node && Face_ && Data );

  Set_Source(Source);

  Set_NodeEntityData();
  Set_FaceEntityData();
 
  Real* coordinates = Vector_Var(COORDINATES);
  coordinates[0]    = Data[1];
  coordinates[1]    = Data[2];
  coordinates[2]    = Data[3];
  
  Real* contact_point = Vector_Var(CONTACT_POINT);
  contact_point[0]    = Data[1];
  contact_point[1]    = Data[2];
  contact_point[2]    = Data[3];

  Get_Gap_Cur() = Data[4];
  
  Get_Gap_Old() = 0.0;
  
  Real* pushback_dir = Get_Pushback_Dir();
  pushback_dir[0] = Data[5];
  pushback_dir[1] = Data[6];
  pushback_dir[2] = Data[7];

  Real* normal_dir = Get_Normal();
  normal_dir[0]    = Data[8];
  normal_dir[1]    = Data[9];
  normal_dir[2]    = Data[10];

  Get_Time_To_Contact() = Data[12];

  Scalar_Var(NODE_ENTITY_KEY) = node_data_key;

  Real* physical_face_normal = Get_Physical_Face_Normal();
  physical_face_normal[0]    = Physical_Face_Normal[0];
  physical_face_normal[1]    = Physical_Face_Normal[1];
  physical_face_normal[2]    = Physical_Face_Normal[2];

  Get_Node_Area() = -1.0;

  number_shared_faces = 0;

  if (IS_TRACKED && !(IS_TIED || IS_INFSLIP || IS_GLUED) && include_neighbors) {
#ifdef INCLUDE_ALL_FACES
    ContactTopologyEntity<Real>::connection_data* shared_face_info = Face()->NeighborInfo();
    for (int i=0; i<Face()->Edges_Per_Face(); ++i) {
      if (shared_face_info[i].owner>=0) {
        shared_face_data[number_shared_faces++] = shared_face_info[i];
      }
    }
#else
    int nedges, edge1, edge2;
    Real local_coords[3];
    Face()->Compute_Local_Coordinates( POSITION,
                                     coordinates,
                                     local_coords );
    Face()->Get_Close_Edges(local_coords, nedges, edge1, edge2);
    if (Face()->NumberOfNeighbors()>0 && nedges>0) {
      PRECONDITION(number_shared_faces<Face()->NumberOfNeighbors());
      ContactTopologyEntity<Real>::connection_data* shared_face_info = Face()->NeighborInfo();
      if (shared_face_info[edge1].owner>=0) {
        shared_face_data[number_shared_faces++] = shared_face_info[edge1];
      }
      if (nedges>1 && shared_face_info[edge2].owner>=0) {
        shared_face_data[number_shared_faces++] = shared_face_info[edge2];
      }
    }
#endif
  }
 
}

ContactNodeFaceInteraction::ContactNodeFaceInteraction( 
                                       ContactNodeFaceInteraction& nfi )
  : ContactNodeEntityInteraction(nfi.is_tied, nfi.is_infSlip, nfi.is_glued, nfi.is_tracked, nfi.node, nfi.entity, nfi.location, NODE_FACE_INTERACTION, CT_NFI )
{
  std::memcpy( &(DataArray[0]), &(nfi.DataArray[0]), 
	  (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real) );
  node_entity_data      = nfi.node_entity_data;
  entity_entity_data    = nfi.entity_entity_data;
  number_shared_faces   = nfi.number_shared_faces;
  for (int i=0; i<number_shared_faces; ++i) {
    shared_face_data[i] = nfi.shared_face_data[i];
  }
}

ContactNodeFaceInteraction* 
ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     InteractionSource Source,
				     ContactNode<Real>* Node,
				     ContactFace<Real>* Face,
				     Real* Data,
				     int node_entity_key,
				     Real* physical_face_normal,
				     bool is_tied,
                                     bool is_infSlip,
                                     bool is_glued,
				     bool is_tracked,
                                     int  include_neighbors,
                                     const VariableHandle &POSITION )
{
  return new (alloc.New_Frag())
    ContactNodeFaceInteraction( Source, Node, Face, Data, node_entity_key,
				physical_face_normal, is_tied, is_infSlip, is_glued, 
                                is_tracked, include_neighbors, POSITION );
}

ContactNodeFaceInteraction* 
ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactNode<Real>* node )
{
  PRECONDITION( node );
  return new (alloc.New_Frag())
    ContactNodeFaceInteraction( node );
}

ContactNodeFaceInteraction* 
ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
				     ContactFixedSizeAllocator& alloc )
{
  return new (alloc.New_Frag())
    ContactNodeFaceInteraction( );
}

ContactNodeFaceInteraction* 
ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactNodeFaceInteraction& cnfi )
{
  return new (alloc.New_Frag())
    ContactNodeFaceInteraction( cnfi );
}

ContactNodeEntityInteraction *ContactNodeFaceInteraction::New_Instance(ContactFixedSizeAllocator* allocators){
  return new_ContactNodeFaceInteraction(allocators[ContactSearch::ALLOC_ContactNodeFaceInteraction], *this);
}

void ContactNodeFaceInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactNodeFaceInteraction),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactQNodeFaceInteraction allocator" );
}


ContactNodeFaceInteraction::~ContactNodeFaceInteraction()
{}


void BlendUnitVectors(const Real* vec_a, const Real* vec_b, const Real a_ratio, Real* output_vec) {
  Real b_ratio = 1.0-a_ratio;
  output_vec[0] = a_ratio*vec_a[0] + b_ratio*vec_b[0];
  output_vec[1] = a_ratio*vec_a[1] + b_ratio*vec_b[1];
  output_vec[2] = a_ratio*vec_a[2] + b_ratio*vec_b[2];
  Real mag = Normalize(output_vec);
  if(mag == 0.0) {
    output_vec[0] = vec_b[0];
    output_vec[1] = vec_b[1];
    output_vec[2] = vec_b[2];
  }
}

void ComputeZoneNormal(VariableHandle &POSITION, 
                       VariableHandle &FACE_NORMAL, 
                       ContactFace<Real>* face, 
                       const int edge_id, 
                       const Real *pushback_dir,
                       const Real *face_normal,
                       Real *zone_normal) {
  Real edge_normal[3];
  face->Compute_Edge_Normal( POSITION, FACE_NORMAL, edge_id, edge_normal);
  Real pdpmn1 = Dot(pushback_dir,edge_normal);
  zone_normal[0] = pushback_dir[0] - pdpmn1*edge_normal[0];
  zone_normal[1] = pushback_dir[1] - pdpmn1*edge_normal[1];
  zone_normal[2] = pushback_dir[2] - pdpmn1*edge_normal[2];
  Real pbdnmag = Normalize(zone_normal);
  if( pbdnmag == 0.0 ){
    acme::Copy(face_normal, zone_normal);
  }
}

void 
ContactNodeFaceInteraction::Modify_for_Curvature( VariableHandle POSITION,
						  VariableHandle FACE_NORMAL,
						  VariableHandle CURVATURE,
                                                  Real sharp_smooth_curvature,
                                                  ContactSearch::Search_Option_Status multiple_interaction_status )

{ 
  // get the edges the interaction is close to (0, 1 or 2)
  int number, edge_id_1,edge_id_2;
  Real* coordinates = Vector_Var(COORDINATES);
  Face()->Get_Close_Edges( coordinates, number, edge_id_1, edge_id_2 );

  POSTCONDITION( number>= 0 && number<=2 );

  // If we aren't close to any edge return
  if( number == 0 ) return;


  // Initialize the data to its original value in case it doesn't get modified
  Real* pushback_dir = Vector_Var(PUSHBACK_DIR);

  Real* face_normal = Face()->Variable(FACE_NORMAL);
  Real& gap = Scalar_Var(GAP_CUR);

  // Coresponds to 2.0 degrees for nearly flat faces
  static Real alpha_zero = ContactSearch::ComputeCurvatureFromAngle(2.0 ) - ContactSearch::ComputeCurvatureFromAngle( 0.0);
  // Coresponds to 2.0 degress for faces near the sharp/smooth angle (assumed to be ~60 degrees).
  static Real alpha_crit = ContactSearch::ComputeCurvatureFromAngle(62.0) - ContactSearch::ComputeCurvatureFromAngle(60.0); 

  Real Curvature1 = Face()->GetEdgeCurvature(edge_id_1);
  if( number == 1 ){
    //
    //  If multiple interactions then:
    //
    //  If curv < 0.0                               edge is convex, set push_back to face_normal
    //  if 0.0 < curv < sharp_smooth_curvature, edge is moderatly concave leave push_back as is
    //  if sharp_smooth_curvature < curv,       edge is very concave, again use face_normal instead of push_back
    //
    //  Near 0.0 and sharp_smooth_curvature smooth out the transition between push back and normal for faces within
    //  1.0 degrees to reduce solution sensitivity to small tweaks.
    //
    //  NKC, this is strange, is it right?  Why not always use push back when the edge is concave, why switch back to face
    //  normal for the very concave face ????????
    //
    //  If not multiple interactions then:
    //  If curv < 0.0                               edge is convex, set push_back to face_normal
    //  if 0.0 < curv                               edge is concave leave push_back as is
    //
    if(Curvature1 <= 0.0) {
      gap *= Dot(face_normal,pushback_dir);
      acme::Copy(face_normal, pushback_dir);
    } else if(Curvature1 <= alpha_zero) {
      Real pb_dir_new[3];
      BlendUnitVectors(pushback_dir, face_normal, Curvature1/alpha_zero, pb_dir_new);
      gap *= Dot(pushback_dir, pb_dir_new);
      acme::Copy(pb_dir_new, pushback_dir);
    } else if(Curvature1 < sharp_smooth_curvature || (multiple_interaction_status == ContactSearch::INACTIVE)) {
      //  Concave edge, leave pushback_dir and gap as is
    } else if(Curvature1 < sharp_smooth_curvature + alpha_crit) {
      Real pb_dir_new[3];
      BlendUnitVectors(face_normal, pushback_dir, (Curvature1-sharp_smooth_curvature)/alpha_crit, pb_dir_new);
      gap *= Dot(pushback_dir, pb_dir_new);
      acme::Copy(pb_dir_new, pushback_dir);
    } else {
      gap *= Dot(face_normal,pushback_dir);
      acme::Copy(face_normal, pushback_dir);
    }
  } else if( number == 2 ) {
    Real Curvature2 = Face()->GetEdgeCurvature(edge_id_2);
    //
    //     C     |      |       |
    //     A     |      |       |
    //     V     |  B   | B+D   |   D
    //           |      |       |
    //  C        |      |       |
    //  U  alpha +------+-------+-----------
    //  R        |      | A+B+  |
    //  V        | A+B  | C+D   |  C+D
    //           |      |       |
    //  2    0.0 +------+-------+-----------
    //           |      |       |
    //     V     |  A   |  A+C  |   C
    //     E     |      |       |
    //     X     +------+-------+------------
    //            CVEX  0     alpha  CCAVE
    //
    //                  CURVATURE 1
    //
    //  Blend results from each zone with the appropriate weight to get a smoothly varying normal
    //
    //  Note, in Zone_A  new_normal = face_normal
    //        in Zone_D  new_normal = pushback_dir
    //        in Zone B, C, calculate normal as needed
    //
    Real pb_dir_new[3];
    if(Curvature1 <= 0.0) {
      if(Curvature2 <= 0.0) {
        acme::Copy(face_normal, pb_dir_new);
      } else if (Curvature2 < alpha_zero) {
        Real zone_B_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_1, pushback_dir, face_normal, zone_B_normal);
        BlendUnitVectors(zone_B_normal, face_normal, Curvature2/alpha_zero, pb_dir_new);
      } else {
        Real zone_B_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_1, pushback_dir, face_normal, zone_B_normal);
        acme::Copy(zone_B_normal, pb_dir_new);
      }
    } else if(Curvature1 < alpha_zero) {
      if(Curvature2 <= 0.0) {
        Real zone_C_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_2, pushback_dir, face_normal, zone_C_normal);
        BlendUnitVectors(zone_C_normal, face_normal, Curvature1/alpha_zero, pb_dir_new);
      } else if (Curvature2 < alpha_zero) {
        Real curv1_ratio = Curvature1/alpha_zero;
        Real curv2_ratio = Curvature2/alpha_zero;
	Real A_ratio = (1.0-curv1_ratio)*(1.0-curv2_ratio);
        Real B_ratio = (1.0-curv1_ratio)*(    curv2_ratio);
        Real C_ratio = (curv1_ratio    )*(1.0-curv2_ratio);
        Real D_ratio = (curv1_ratio    )*(    curv2_ratio);

        Real zone_B_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_1, pushback_dir, face_normal, zone_B_normal);

        Real zone_C_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_2, pushback_dir, face_normal, zone_C_normal);
        pb_dir_new[0] = A_ratio*face_normal[0] + B_ratio*zone_B_normal[0] + C_ratio*zone_C_normal[0] + D_ratio*pushback_dir[0];
        pb_dir_new[1] = A_ratio*face_normal[1] + B_ratio*zone_B_normal[1] + C_ratio*zone_C_normal[1] + D_ratio*pushback_dir[1];
        pb_dir_new[2] = A_ratio*face_normal[2] + B_ratio*zone_B_normal[2] + C_ratio*zone_C_normal[2] + D_ratio*pushback_dir[2];
        Real mag = Normalize(pb_dir_new);
        if(mag == 0.0) {
          acme::Copy(face_normal, pb_dir_new);
	}
      } else {
        Real zone_B_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_1, pushback_dir, face_normal, zone_B_normal);
        BlendUnitVectors(pushback_dir, zone_B_normal, Curvature1/alpha_zero, pb_dir_new);
      }
    } else {
      if(Curvature2 <= 0.0) {
        Real zone_C_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_2, pushback_dir, face_normal, zone_C_normal);
        acme::Copy(zone_C_normal, pb_dir_new);
      } else if (Curvature2 < alpha_zero) {
        Real zone_C_normal[3];
        ComputeZoneNormal(POSITION, FACE_NORMAL, Face(), edge_id_2, pushback_dir, face_normal, zone_C_normal);
        BlendUnitVectors(pushback_dir, zone_C_normal, Curvature2/alpha_zero, pb_dir_new);
      } else {
        acme::Copy(pushback_dir, pb_dir_new);
      }
    }
    gap *= Dot(pushback_dir, pb_dir_new);
    acme::Copy(pb_dir_new, pushback_dir);
  }
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for DLB and for
// sending interactions from the secondary to the primary decomposition
//----------------------------------------------------------------------
int ContactNodeFaceInteraction::Size()
{
  return(ContactInteractionEntity<Real>::Size()+
         sizeof(entity_data)+sizeof(int)+
         number_shared_faces*sizeof(ContactTopologyEntity<Real>::connection_data)+
         (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real));
}

char* ContactNodeFaceInteraction::Pack( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Pack( buffer );
  buff += ContactInteractionEntity<Real>::Size();
  
  int  cnt     = 0;
  int* i_buf   = reinterpret_cast<int*> (buff);
  cnt         += PackEntityData(&entity_entity_data, &i_buf[cnt]);
  i_buf[cnt++] = location | is_tied<<2 | is_infSlip<<3 | is_glued<<4 | is_tracked<<5 | number_shared_faces<<6;
  for (int i=0; i<number_shared_faces; ++i) {
    cnt += PackConnection(&shared_face_data[i], &i_buf[cnt]);
  }
  
  buff += cnt*sizeof(int);
  cnt   = (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real);
  std::memcpy( buff,DataArray,cnt);
  buff += cnt;
  return buff;
}

char* ContactNodeFaceInteraction::Unpack( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Unpack( buff );
  buff += ContactInteractionEntity<Real>::Size();
  
  int cnt    = 0;
  int* i_buf = reinterpret_cast<int*> (buff);
  cnt       += UnPackEntityData(&entity_entity_data, &i_buf[cnt]);
  location            = (i_buf[cnt]   ) & 0x3;
  is_tied             = (i_buf[cnt]>>2) & 0x1;
  is_infSlip          = (i_buf[cnt]>>3) & 0x1;
  is_glued            = (i_buf[cnt]>>4) & 0x1;
  is_tracked          = (i_buf[cnt]>>5) & 0x1;
  number_shared_faces = (i_buf[cnt]>>6);
  ++cnt;
  for (int i=0; i<number_shared_faces; ++i) {
    cnt += UnPackConnection(&shared_face_data[i], &i_buf[cnt]);
  }
  
  buff += cnt*sizeof(int);
  cnt   = (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real);
  std::memcpy( DataArray,buff,cnt);
  buff += cnt;
  return buff;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for sending
// interactions from the primary to the secondary decomposition
//-----------------------------------------------------------------
  #define DEFINE_DATA_TO_SEND int vector_data   = 1;\
                              int scalar_data   = 0;\
                              int vector_offset = PHYSICAL_FACE_NORMAL;\
                              int scalar_offset = GAP_CUR;\
                              if (is_tied||is_infSlip||is_glued) {\
                                scalar_data   = 7;\
                                scalar_offset = GAP_CUR;\
                                vector_data   = 5;\
                                vector_offset = NORMAL_DIR;\
                              } else if (is_tracked){ \
                                scalar_data   = 2;\
                                scalar_offset = SOURCE;\
                                vector_data   = 4;\
                                vector_offset = NORMAL_DIR;\
                              }

int ContactNodeFaceInteraction::Size_ForSecondary()
{
  DEFINE_DATA_TO_SEND;
  return(ContactInteractionEntity<Real>::Size_ForSecondary()+
         3*sizeof(int)+(scalar_data+3*vector_data)*sizeof(Real));
}

char* ContactNodeFaceInteraction::Pack_ForSecondary( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Pack_ForSecondary( buffer );
  buff += ContactInteractionEntity<Real>::Size_ForSecondary();
  
  int  cnt     = 0;
  int* i_buf   = reinterpret_cast<int*> (buff);

  i_buf[cnt++] = entity_entity_data.host_gid[0];
  i_buf[cnt++] = entity_entity_data.host_gid[1];

  i_buf[cnt++] = location | is_tied<<2 | is_infSlip<<3 | is_glued<<4 | is_tracked<<5;
  
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

char* ContactNodeFaceInteraction::Unpack_ForSecondary( char* buffer )
{
  char* buff = buffer;
  buff += ContactInteractionEntity<Real>::Size_ForSecondary();
  
  int cnt    = 0;
  int* i_buf = reinterpret_cast<int*> (buff);

  entity_entity_data.host_gid[0] = i_buf[cnt++];
  entity_entity_data.host_gid[1] = i_buf[cnt++];

  location            = (i_buf[cnt]   ) & 0x3;
  is_tied             = (i_buf[cnt]>>2) & 0x1;
  is_infSlip          = (i_buf[cnt]>>3) & 0x1;
  is_glued            = (i_buf[cnt]>>4) & 0x1;
  is_tracked          = (i_buf[cnt]>>5) & 0x1;
  number_shared_faces = 0;
  ++cnt;
  
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

int ContactNodeFaceInteraction::PackConnection( ContactTopologyEntity<Real>::connection_data* data,
                                                       int* buffer )
{
  int* buf = buffer;
  *buf++   = data->owner;
  *buf++   = data->owner_proc_array_index;
  *buf++   = data->host_gid[0];
  *buf++   = data->host_gid[1]; 
  return 4;
}

int ContactNodeFaceInteraction::UnPackConnection( ContactTopologyEntity<Real>::connection_data* data, 
                                                         int* buffer)
{
  int* buf                     = buffer;
  data->owner                  = *buf++;
  data->owner_proc_array_index = *buf++;
  data->host_gid[0]            = *buf++;
  data->host_gid[1]            = *buf++;
  return 4;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions to be used for restart
//--------------------------------------------------------------------
int ContactNodeFaceInteraction::Restart_Size()
{
  return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS+
         2*sizeof(entity_data)/sizeof(int)+4;
}

void ContactNodeFaceInteraction::Restart_Pack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += PackEntityData(&node_entity_data,   &buf_loc[cnt]);
  cnt += PackEntityData(&entity_entity_data, &buf_loc[cnt]);
  buf_loc[cnt++] = is_tied;
  buf_loc[cnt++] = is_infSlip;
  buf_loc[cnt++] = is_glued;
  buf_loc[cnt++] = is_tracked;

  // Pack the data
  std::memcpy( &buf_loc[cnt], &DataArray, 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );

}

void ContactNodeFaceInteraction::Restart_Unpack( Real* buffer )
{
  int cnt = 0;
  Real* buf_loc = buffer;

  cnt += UnPackEntityData(&node_entity_data,   &buf_loc[cnt]);
  cnt += UnPackEntityData(&entity_entity_data, &buf_loc[cnt]);
  is_tied    = (int) buf_loc[cnt++];
  is_infSlip = (int) buf_loc[cnt++];
  is_glued   = (int) buf_loc[cnt++];
  is_tracked = (int) buf_loc[cnt++];

  // Unpack the Data Array
  std::memcpy( &DataArray, &buf_loc[cnt], 
	  sizeof(Real)*(NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS ) );
}


void ContactNodeFaceInteraction::Modify_for_Normal_Smoothing(
                                              VariableHandle CURRENT_POSITION,
					      VariableHandle NODE_NORMAL,
					      VariableHandle FACE_NORMAL,
					      VariableHandle CURVATURE,
		               ContactSearch::Smoothing_Resolution resolution,
					      Real smoothing_distance,
                                              Real sharp_smooth_curvature )
{
  Real smoothed_normal[3];

  Real* coordinates = Vector_Var(COORDINATES);

  Face()->Smooth_Normal( CURRENT_POSITION, NODE_NORMAL, FACE_NORMAL,
		       CURVATURE, resolution, smoothing_distance, coordinates, 
		       &(smoothed_normal[0]), sharp_smooth_curvature );
  Real* normal_dir = Vector_Var(NORMAL_DIR);
  normal_dir[0] = smoothed_normal[0];
  normal_dir[1] = smoothed_normal[1];
  normal_dir[2] = smoothed_normal[2];
  
}

void 
ContactNodeFaceInteraction::Update_Tied_Interaction(VariableHandle POSITION,
						    VariableHandle FACE_NORMAL,
						    Real gap0,
						    ScratchVariable &NUM_PF,
						    ScratchVariable &PF_1,
						    ScratchVariable &PF_2,
						    ScratchVariable &PF_3 )
{
  PRECONDITION( is_tied == true || is_infSlip == true || is_glued == true );
  PRECONDITION( Face() );

  Set_FaceEntityData();
  
  // Set the normal direction to be the face normal
  Real* normal_dir  = Vector_Var(NORMAL_DIR);
  Real* face_normal = Face()->Variable(FACE_NORMAL);

  //
  // BUG_KHP: We need to update the normal_direction here instead of using
  // BUG_KHP: the face normal.
  //

  normal_dir[0] = face_normal[0];
  normal_dir[1] = face_normal[1];
  normal_dir[2] = face_normal[2];

  // Set the initial gap
  Real& gap_initial = Scalar_Var(GAP_INT);
  gap_initial = gap0;

  // Compute the pushback direction to be the difference between the node 
  // and the interpolated face position.  The pushback dir is the unit
  // vector along this direction and the gap is magnitude of this vector.
  Real& gap            = Scalar_Var(GAP_CUR);
  Real* pushback_dir   = Vector_Var(PUSHBACK_DIR);
  Real* node_position  = node->Variable(POSITION);
  Real contact_point[] = {0.0, 0.0, 0.0};
  Real shape_functions[MAX_NODES_PER_FACE];
  Real* local_coords = Vector_Var(COORDINATES);
  PRECONDITION( Face() != NULL );
  PRECONDITION( Face()->Nodes_Per_Face() <= MAX_NODES_PER_FACE );
  Face()->Evaluate_Shape_Functions( local_coords, shape_functions );

  for( int i=0 ; i<Face()->Nodes_Per_Face() ; ++i ){
    ContactNode<Real>* nd = Face()->Node(i);
    Real* pp = nd->Variable(POSITION);
    contact_point[0] += shape_functions[i]*pp[0];
    contact_point[1] += shape_functions[i]*pp[1];
    contact_point[2] += shape_functions[i]*pp[2];
  }
  pushback_dir[0] = contact_point[0] - node_position[0];
  pushback_dir[1] = contact_point[1] - node_position[1];
  pushback_dir[2] = contact_point[2] - node_position[2];

  Real mag = Magnitude(pushback_dir);

  if( mag > 0.0){
    Real tmp_mag     = 1.0/mag;
    pushback_dir[0] *= tmp_mag;
    pushback_dir[1] *= tmp_mag;
    pushback_dir[2] *= tmp_mag;
    Real dot = Dot(pushback_dir, face_normal);

    int sign = dot >= 0.0 ? -1 : 1;
    gap = mag*sign;
  } else {
    gap = 0.0;
    pushback_dir[0] = -face_normal[0];
    pushback_dir[1] = -face_normal[1];
    pushback_dir[2] = -face_normal[2];
  }
  
  // Update the Physical Face Normal
  // Assume rotations are small so we can just take the one with the
  // smallest angle between itself and the old physical face normal
  int num_physical_faces = (int) *(NUM_PF.Get_Scratch(node->ProcArrayIndex()));
  if( num_physical_faces ){
    Real* physical_face_normal = Vector_Var( PHYSICAL_FACE_NORMAL );
    Real* pfn[3];
    pfn[0] = PF_1.Get_Scratch(node->ProcArrayIndex());
    pfn[1] = PF_2.Get_Scratch(node->ProcArrayIndex());
    pfn[2] = PF_3.Get_Scratch(node->ProcArrayIndex());
    int closest = 0;
    Real most_opposed = 2.0;
    if( num_physical_faces == 1 ){
      physical_face_normal[0] = pfn[0][0];
      physical_face_normal[1] = pfn[0][1];
      physical_face_normal[2] = pfn[0][2];
    } else {
      Find_Physical_Face(num_physical_faces, pfn, normal_dir, closest, most_opposed);
      physical_face_normal[0] = pfn[closest][0];
      physical_face_normal[1] = pfn[closest][1];
      physical_face_normal[2] = pfn[closest][2];
    }
  }
}

void 
ContactNodeFaceInteraction::Update_Tied_Interaction(VariableHandle POSITION,
						    VariableHandle FACE_NORMAL,
						    Real gap0,
						    VariableHandle NODE_NORMAL )
{
  PRECONDITION( is_tied == true || is_infSlip == true || is_glued == true );
  PRECONDITION( Face() );

  // Set the normal direction to be the face normal
  Real* normal_dir = Vector_Var(NORMAL_DIR);
  Real* face_normal = Face()->Variable(FACE_NORMAL);
  normal_dir[0] = face_normal[0];
  normal_dir[1] = face_normal[1];
  normal_dir[2] = face_normal[2];

  // Set the initial gap
  Real& gap_initial = Scalar_Var(GAP_INT);
  gap_initial = gap0;

  // Compute the pushback direction to be the difference between the node 
  // and the interpolated face position.  The pushback dir is the unit
  // vector along this direction and the gap is magnitude of this vector.
  Real& gap = Scalar_Var(GAP_CUR);
  Real* pushback_dir = Vector_Var(PUSHBACK_DIR);
  Real* node_position = node->Variable(POSITION);
  Real  contact_point[] = {0.0, 0.0, 0.0};
  Real  shape_functions[MAX_NODES_PER_FACE];
  Real* local_coords = Vector_Var(COORDINATES);
  ContactFace<Real>* face = Face();
  PRECONDITION( face != NULL );
  PRECONDITION( face->Nodes_Per_Face() <= MAX_NODES_PER_FACE );
  face->Evaluate_Shape_Functions( local_coords, shape_functions );

  for( int i=0 ; i<face->Nodes_Per_Face() ; ++i ){
    ContactNode<Real>* nd = face->Node(i);
    Real* pp = nd->Variable(POSITION);
    contact_point[0] += shape_functions[i]*pp[0];
    contact_point[1] += shape_functions[i]*pp[1];
    contact_point[2] += shape_functions[i]*pp[2];
  }
  pushback_dir[0] = contact_point[0] - node_position[0];
  pushback_dir[1] = contact_point[1] - node_position[1];
  pushback_dir[2] = contact_point[2] - node_position[2];
  Real mag = Magnitude(pushback_dir);
  if( mag > 0.0){
    Real tmp_mag     = 1.0/mag;
    pushback_dir[0] *= tmp_mag;
    pushback_dir[1] *= tmp_mag;
    pushback_dir[2] *= tmp_mag;
    Real dot = Dot(pushback_dir, face_normal);
    int sign = dot >= 0.0 ? -1 : 1;
    gap = mag*sign;
  } else {
    gap = 0.0;
    pushback_dir[0] = -face_normal[0];
    pushback_dir[1] = -face_normal[1];
    pushback_dir[2] = -face_normal[2];
  }
  
  // Update the Physical Face Normal
  Real* physical_face_normal = Vector_Var( PHYSICAL_FACE_NORMAL );
  Real* node_normal = node->Variable(NODE_NORMAL);
  physical_face_normal[0] = node_normal[0];
  physical_face_normal[1] = node_normal[1];
  physical_face_normal[2] = node_normal[2];
}

void 
ContactNodeFaceInteraction::Update_Tied_Interaction(VariableHandle POSITION,
						    VariableHandle FACE_NORMAL,
						    Real gap0,
						    Real* pfn )
{
  PRECONDITION( is_tied == true || is_infSlip == true || is_glued == true );
  PRECONDITION( Face() );

  // Set the normal direction to be the face normal
  Real* normal_dir = Vector_Var(NORMAL_DIR);
  Real* face_normal = Face()->Variable(FACE_NORMAL);
  normal_dir[0] = face_normal[0];
  normal_dir[1] = face_normal[1];
  normal_dir[2] = face_normal[2];

  // Set the initial gap
  Real& gap_initial = Scalar_Var(GAP_INT);
  gap_initial = gap0;

  // Compute the pushback direction to be the difference between the node 
  // and the interpolated face position.  The pushback dir is the unit
  // vector along this direction and the gap is magnitude of this vector.
  Real& gap = Scalar_Var(GAP_CUR);
  Real* pushback_dir = Vector_Var(PUSHBACK_DIR);
  Real* node_position = node->Variable(POSITION);
  Real contact_point[] = {0.0, 0.0, 0.0};
  Real shape_functions[MAX_NODES_PER_FACE];
  Real* local_coords = Vector_Var(COORDINATES);
  ContactFace<Real>* face = Face();
  PRECONDITION( face != NULL );
  PRECONDITION( face->Nodes_Per_Face() <= MAX_NODES_PER_FACE );
  face->Evaluate_Shape_Functions( local_coords, shape_functions );

  for( int i=0 ; i<face->Nodes_Per_Face() ; ++i ){
    ContactNode<Real>* nd = face->Node(i);
    Real* pp = nd->Variable(POSITION);
    contact_point[0] += shape_functions[i]*pp[0];
    contact_point[1] += shape_functions[i]*pp[1];
    contact_point[2] += shape_functions[i]*pp[2];
  }
  pushback_dir[0] = contact_point[0] - node_position[0];
  pushback_dir[1] = contact_point[1] - node_position[1];
  pushback_dir[2] = contact_point[2] - node_position[2];
  Real mag = Magnitude(pushback_dir);
  if( mag > 0.0){
    Real tmp_mag     = 1.0/mag;
    pushback_dir[0] *= tmp_mag;
    pushback_dir[1] *= tmp_mag;
    pushback_dir[2] *= tmp_mag;
    Real dot = Dot(pushback_dir, face_normal);
    int sign = dot >= 0.0 ? -1 : 1;
    gap = mag*sign;
  } else {
    gap = 0.0;
    pushback_dir[0] = -face_normal[0];
    pushback_dir[1] = -face_normal[1];
    pushback_dir[2] = -face_normal[2];
  }
  
  // Update the Physical Face Normal
  Real* physical_face_normal = Vector_Var( PHYSICAL_FACE_NORMAL );
  physical_face_normal[0] = pfn[0];
  physical_face_normal[1] = pfn[1];
  physical_face_normal[2] = pfn[2];
}

void 
ContactNodeFaceInteraction::Update_Tied_Interaction(VariableHandle POSITION,
						    Real* buffer,
						    Real  gap0,
						    Real* pfn )
{
  PRECONDITION( is_tied == true || is_infSlip == true || is_glued == true );

  ContactHostGlobalID Face_GID( FaceEntityData()->host_gid[0], 
                                FaceEntityData()->host_gid[1] );
  
  // Set the normal direction to be the face normal
  int j = 0;
  Real  face_normal[3];
  Real* normal_dir = Vector_Var(NORMAL_DIR);
  ContactSearch::ContactFace_Type face_type = (ContactSearch::ContactFace_Type)buffer[j++];
  face_normal[0] = buffer[j++];
  face_normal[1] = buffer[j++];
  face_normal[2] = buffer[j++];
  normal_dir[0]  = face_normal[0];
  normal_dir[1]  = face_normal[1];
  normal_dir[2]  = face_normal[2];

  // Set the initial gap
  Real& gap_initial = Scalar_Var(GAP_INT);
  gap_initial       = gap0;

  // Compute the pushback direction to be the difference between the node 
  // and the interpolated face position.  The pushback dir is the unit
  // vector along this direction and the gap is magnitude of this vector.
  Real shape_functions[MAX_NODES_PER_FACE];
  Real& gap            = Scalar_Var(GAP_CUR);
  Real* pushback_dir   = Vector_Var(PUSHBACK_DIR);
  Real* node_position  = node->Variable(POSITION);
  Real* local_coords   = Vector_Var(COORDINATES);
  
  int nodes_per_face = ContactSearch::Number_Nodes_Per_Face(face_type);
  switch (face_type) {
  case ContactSearch::QUADFACEL4:
    ContactQuadFaceL4<Real>::Compute_Shape_Functions( local_coords, 
                                                shape_functions );
    break;
  case ContactSearch::QUADFACEQ8:
    ContactQuadFaceQ8::Compute_Shape_Functions( local_coords, 
                                                shape_functions );
    break;
  case ContactSearch::QUADFACEQ9:
    ContactQuadFaceQ9::Compute_Shape_Functions( local_coords, 
                                                shape_functions );
    break;
  case ContactSearch::TRIFACEL3:
    ContactTriFaceL3<Real>::Compute_Shape_Functions( local_coords, 
                                               shape_functions );
    break;
  case ContactSearch::TRIFACEQ6:
    ContactTriFaceQ6::Compute_Shape_Functions( local_coords, 
                                               shape_functions );
    break;
  case ContactSearch::SHELLQUADFACEL4:
    ContactShellQuadFaceL4<Real>::Compute_Shape_Functions( local_coords, 
                                                     shape_functions );
    break;
  case ContactSearch::SHELLTRIFACEL3:
    ContactShellTriFaceL3<Real>::Compute_Shape_Functions( local_coords, 
                                                    shape_functions );
    break;
  case ContactSearch::LINEFACEL2:
    ContactLineFaceL2::Compute_Shape_Functions( local_coords, 
                                                shape_functions );
    break;
  case ContactSearch::LINEFACEQ3:
    ContactLineFaceQ3::Compute_Shape_Functions( local_coords, 
                                                shape_functions );
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  
  Real contact_point[] = {0.0, 0.0, 0.0};
  for( int i=0 ; i<nodes_per_face ; ++i ){
    Real xyz[3];
    xyz[0]  = buffer[0];
    xyz[1]  = buffer[1];
    xyz[2]  = buffer[2];
    contact_point[0] += shape_functions[i]*buffer[j++];
    contact_point[1] += shape_functions[i]*buffer[j++];
    contact_point[2] += shape_functions[i]*buffer[j++];
  }
  
  pushback_dir[0] = contact_point[0] - node_position[0];
  pushback_dir[1] = contact_point[1] - node_position[1];
  pushback_dir[2] = contact_point[2] - node_position[2];
  Real mag = Magnitude( pushback_dir );
  if( mag > 0.0){
    Real tmp_mag     = 1.0/mag;
    pushback_dir[0] *= tmp_mag;
    pushback_dir[1] *= tmp_mag;
    pushback_dir[2] *= tmp_mag;
    Real dot = Dot(pushback_dir,face_normal);
    int sign = dot >= 0.0 ? -1 : 1;
    gap = mag*sign;
  } else {
    gap = 0.0;
    pushback_dir[0] = -face_normal[0];
    pushback_dir[1] = -face_normal[1];
    pushback_dir[2] = -face_normal[2];
  }
  
  // Update the Physical Face Normal
  Real* physical_face_normal = Vector_Var( PHYSICAL_FACE_NORMAL );
  physical_face_normal[0] = pfn[0];
  physical_face_normal[1] = pfn[1];
  physical_face_normal[2] = pfn[2];
}

void ContactNodeFaceInteraction::Connect_Entity( ContactTopology* topology){
  Connect_Face(topology);
}

void ContactNodeFaceInteraction::Get_Face_Normal(Real *return_normal, 
                                                 ContactTopology* search_topology) {
  VariableHandle FACE_NORMAL = search_topology->Variable_Handle( ContactTopology::Face_Normal );  
  Real* normal      = Face()->Variable(FACE_NORMAL);
  return_normal[0] = normal[0];
  return_normal[1] = normal[1];
  return_normal[2] = normal[2];
}

void ContactNodeFaceInteraction::Get_Avg_Face_Normal(Real *return_normal, 
                                                     ContactTopology* search_topology,
                                                     VariableHandle NODE_COORDS) {
  VariableHandle NODE_NORMAL = search_topology->Variable_Handle( ContactTopology::Node_Normal );

  Real  shape_functions[MAX_NODES_PER_FACE];

  PRECONDITION( Face()->Nodes_Per_Face() <= MAX_NODES_PER_FACE );
  // Note: the coordinates at the point in the algorithm where this is
  // called are the global coordinates.  We must first convert from the
  // global coordinates to the local coordinates.
  Real local_coords[3];
  Real* coordinates = Vector_Var(ContactNodeFaceInteraction::COORDINATES);
  Face()->Compute_Local_Coordinates( NODE_COORDS, coordinates, local_coords );
  Face()->Evaluate_Shape_Functions( local_coords, shape_functions );
  for( int i=0 ; i<Face()->Nodes_Per_Face() ; ++i ){
    ContactNode<Real>* nd = Face()->Node(i);
    Real* ndnrml = nd->Variable(NODE_NORMAL); // assembled node normal
    return_normal[0] += shape_functions[i]*ndnrml[0];
    return_normal[1] += shape_functions[i]*ndnrml[1];
    return_normal[2] += shape_functions[i]*ndnrml[2];
  }
  Normalize(return_normal);
}

void ContactNodeFaceInteraction::Delete_Frag(ContactFixedSizeAllocator* allocators) {
  this->~ContactNodeFaceInteraction();
  allocators[ContactSearch::ALLOC_ContactNodeFaceInteraction].Delete_Frag(this);
}
