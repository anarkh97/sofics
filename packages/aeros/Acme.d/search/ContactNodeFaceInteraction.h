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


#ifndef ContactNodeFaceInteraction_h_
#define ContactNodeFaceInteraction_h_

#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactNode.h"
#include "ContactFace.h"
#include "ContactNodeEntityInteraction.h"
#include "ContactTopologyEntity.h"
#include "ContactTopologyEntityList.h"
#include "ContactTopology.h"

class CString;
class ContactTopologyEntityHash;
class ContactTopologyEntityList;
class ContactFixedSizeAllocator;

class ContactNodeFaceInteraction : public ContactNodeEntityInteraction {
  
 public:

  ContactNodeFaceInteraction( );
  ContactNodeFaceInteraction( ContactNode<Real>* );
  ContactNodeFaceInteraction( InteractionSource, ContactNode<Real>*, ContactFace<Real>*, 
			      Real*, int, Real*, bool, bool, bool, bool,
                              int, const VariableHandle& );
  ContactNodeFaceInteraction( ContactNodeFaceInteraction& );
  static ContactNodeFaceInteraction* new_ContactNodeFaceInteraction(
            ContactFixedSizeAllocator&, InteractionSource,  
	    ContactNode<Real>*, ContactFace<Real>*, Real*, int, Real*, 
            bool, bool, bool, bool, int, const VariableHandle& );
  static ContactNodeFaceInteraction* new_ContactNodeFaceInteraction(
	     ContactFixedSizeAllocator&  );
  static ContactNodeFaceInteraction* new_ContactNodeFaceInteraction(
	     ContactFixedSizeAllocator& , ContactNode<Real>* node );
  static ContactNodeFaceInteraction* new_ContactNodeFaceInteraction(
	     ContactFixedSizeAllocator&, ContactNodeFaceInteraction& );
  ~ContactNodeFaceInteraction();
  
  virtual ContactNodeEntityInteraction *New_Instance(ContactFixedSizeAllocator* allocators);

#ifndef CONTACT_NO_MPI
  inline ContactZoltanLID& Zoltan_LID();
  inline ContactZoltanGID& Zoltan_GID();
  inline void ZoltanFaceLID(LB_ID_PTR lid, int flag=0); 
  inline void ZoltanFaceGID(LB_ID_PTR gid);
#endif

  // Parallel packing/unpacking functions
  virtual int   Size();
  virtual char* Pack( char* buffer );
  virtual char* Unpack( char* buffer );
  
  virtual int   Size_ForSecondary();
  virtual char* Pack_ForSecondary( char* buffer );
  virtual char* Unpack_ForSecondary( char* buffer );

  inline ContactFace<Real>* Face();
  inline entity_data* FaceEntityData();
  inline int Set_FaceEntityData();

  int Is_Better(ContactNodeFaceInteraction*, VariableHandle, VariableHandle );
  void Modify_for_Curvature( VariableHandle, VariableHandle, VariableHandle,
                             double cos_sharp_smooth_angle, 
                             ContactSearch::Search_Option_Status multiple_interaction_status  );
  void Modify_for_Normal_Smoothing( VariableHandle, VariableHandle, 
				    VariableHandle, VariableHandle, 
				    ContactSearch::Smoothing_Resolution, 
				    Real,
                                    Real cos_sharp_smooth_angle );

  virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real,
                                        ScratchVariable&, ScratchVariable&,
                                        ScratchVariable&, ScratchVariable& );
  virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real, VariableHandle );
  virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real, Real* );
  virtual void Update_Tied_Interaction( VariableHandle, Real*, Real, Real* );

  virtual void Connect_Entity( ContactTopology* );

  inline void Connect_Face( ContactTopologyEntityList& );
  inline void Connect_Face( ContactTopologyEntityHash& );
  inline void Connect_Face( ContactTopology* );
  inline void Connect_Face( ContactFace<Real>* );

  // Restart Pack/Unpack functions
  virtual int  Restart_Size();
  virtual void Restart_Pack( Real* buffer );
  virtual void Restart_Unpack( Real* buffer );

  virtual void Get_Face_Normal(Real *return_normal, 
                       ContactTopology* search_topology);
  virtual void Get_Avg_Face_Normal(Real *return_normal, 
                           ContactTopology* search_topology,
                           VariableHandle NODE_COORDS);

  virtual void Delete_Frag(ContactFixedSizeAllocator* allocators);

  inline int NumSharedFaces();
  inline ContactTopologyEntity<Real>::connection_data* SharedFaceData(int i);
    
 protected:

 private:
  inline int  PackConnection(ContactTopologyEntity<Real>::connection_data*, int*);
  inline int  UnPackConnection(ContactTopologyEntity<Real>::connection_data*, int*);
  
  ContactTopologyEntity<Real>::connection_data shared_face_data[4];
  int number_shared_faces;
#ifndef CONTACT_NO_MPI
  ContactZoltanLID zoltan_lid;
  ContactZoltanGID zoltan_gid;
#endif
};

//-------------------------------------------------------------------------------

#ifndef CONTACT_NO_MPI
inline ContactZoltanLID& ContactNodeFaceInteraction::Zoltan_LID()
{
  return zoltan_lid;
}

inline ContactZoltanGID& ContactNodeFaceInteraction::Zoltan_GID()
{
  return zoltan_gid;
}

inline void ContactNodeFaceInteraction::ZoltanFaceLID(LB_ID_PTR lid, int flag) 
{
  //if (flag)
  //{
  //  zoltan_lid.ZoltanLID(CT_FACE, 
  //                       entity_entity_data.index_in_owner_proc_array, 
  //                       lid);
  //}
  //else
  //{
  //  zoltan_lid.ZoltanLID(CT_FACE, 
  //                       entity_entity_data.index_in_proc_array, 
  //                       lid);
  //}
  zoltan_lid.ZoltanLID(CT_FACE, 
                       entity_entity_data.index_in_owner_proc_array, 
                       lid);
}

inline void ContactNodeFaceInteraction::ZoltanFaceGID(LB_ID_PTR gid) 
{
  zoltan_gid.ZoltanGID(CT_FACE, 
                       entity_entity_data.host_gid[0],
                       entity_entity_data.host_gid[1],
                       gid);
}
#endif

inline ContactFace<Real>* ContactNodeFaceInteraction::Face()
{
  return static_cast<ContactFace<Real> *>(entity);
}

inline ContactInteractionEntity<Real>::entity_data* ContactNodeFaceInteraction::FaceEntityData()
{
  return &entity_entity_data;
}

inline int ContactNodeFaceInteraction::Set_FaceEntityData() 
{
  int temp = Set_Entity_Data();
  return temp;
}

inline void ContactNodeFaceInteraction::Connect_Face( ContactTopologyEntityList& hash_table )
{
  entity = static_cast<ContactFace<Real>*>(hash_table.Find( &entity_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( face );
}

inline void ContactNodeFaceInteraction::Connect_Face( ContactTopologyEntityHash& hash_table )
{
  ContactHostGlobalID gid( entity_entity_data.host_gid[0], 
                           entity_entity_data.host_gid[1] );
  entity = static_cast<ContactFace<Real>*>(hash_table.find(gid));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( face );
}

inline void ContactNodeFaceInteraction::Connect_Face( ContactTopology* topology )
{ 
  entity = topology->FaceList()->Find( &entity_entity_data );
}

inline void ContactNodeFaceInteraction::Connect_Face( ContactFace<Real>* Face )
{
  entity = Face;
}

inline int ContactNodeFaceInteraction::NumSharedFaces()
{
  return number_shared_faces;
}

inline ContactTopologyEntity<Real>::connection_data* ContactNodeFaceInteraction::SharedFaceData(int i) 
{
  return &shared_face_data[i];
}
#endif // ContactNodeFaceInteraction_h_
