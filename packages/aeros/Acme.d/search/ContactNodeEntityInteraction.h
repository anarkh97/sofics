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


#ifndef ContactNodeEntityInteraction_h_
#define ContactNodeEntityInteraction_h_

#include "Contact_Defines.h"
#include "CString.h"
#include "ContactInteractionEntity.h"
#include "ContactNode.h"
#include "contact_assert.h"

class CString;
class ContactTopologyEntityHash;
class ContactTopologyEntityList;
class ContactFixedSizeAllocator;
class ContactTopology;
class ContactFixedSizeAllocator;

class ContactNodeEntityInteraction : public ContactInteractionEntity<Real> {

 public:
 
   enum InteractionType{NODE_FACE_INTERACTION    = 0,
                        NODE_SURFACE_INTERACTION = 1};

   enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  NEI_SCALAR_VAR
#define NEI_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
                      NUMBER_SCALAR_VARS };

   enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  NEI_VECTOR_VAR
#define NEI_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		      NUMBER_VECTOR_VARS };

   //Values to be stored in the scratch variable "SOURCE" in the node face interactions and
   //the node surface interactions. Put here to avoid name conflicts as well as to prevent
   //duplicate code.
   enum InteractionSource { UNKNOWN_SOURCE=-1,CLOSEST_POINT_PROJECTION_1=1, 
                            CLOSEST_POINT_PROJECTION_2, 
			    MOVING_INTERSECTION, 
			    CLOSEST_POINT_PROJECTION_2_FROM_MOVING,
			    RETRIEVED_TIED,
                            RETRIEVED_GLUED };

   ContactNodeEntityInteraction(bool is_tied_,
                                bool is_infSlip_,
                                bool is_tracked_,
                                bool is_glued_,
                                //bool track_initial_gap_,
                                ContactNode<Real>* node_,
                                ContactTopologyEntity<Real>* entity_,
                                int Location,
                                InteractionType type_,
                                ContactType base_type_);
   virtual ~ContactNodeEntityInteraction();

   /**
    *  This method will create a near identical copy of the current instance. The interaction
    *  source is not copied as this is a new instance and the source should have changed. This 
    *  makes creating a new node face or node surface interaction in the retrieve tied or
    *  retrieve tracked interactions functions generic enough to keep from having lots of if
    *  statements just for creating the right interaction type. The actual work is done by the
    *  the sub class so you are guaranteed to get the right kind of interaction. 
    */
   virtual ContactNodeEntityInteraction* New_Instance(ContactFixedSizeAllocator* allocators) = 0;

   inline int  DataArray_Length();  
   inline void Initialize_Memory();

   inline ContactTopologyEntity<Real>*                 Entity();
   inline const ContactInteractionEntity<Real>::entity_data* Get_Entity_Data();
   inline int                                          Set_Entity_Data();
   inline int                                          Get_Entity_Owner();

   inline ContactNode<Real>*                                 Node();
   inline const ContactInteractionEntity<Real>::entity_data* NodeEntityData();
   inline int                                          Set_NodeEntityData();

   inline int Get_Entity_Key();
   inline int Get_Node_Key();

   void Modify_for_Kinematic_Constraints(VariableHandle NUM_KIN_CONSTR, VariableHandle KIN_CONSTR_VECTOR, bool& Would_Violate_Constraints);

   inline bool Is_Tied();
   inline void Is_Tied( bool IS_TIED );

   inline bool Is_InfSlip();
   inline void Is_InfSlip(bool IS_INFSLIP);

   inline bool Is_Tracked();
   inline void Is_Tracked( bool IS_TRACKED );

   inline bool Is_Glued();
   inline void Is_Glued( bool IS_GLUED );

   inline int Location();

   void Connect_Node( ContactTopologyEntityList& );
   void Connect_Node( ContactTopologyEntityHash& );
   void Connect_Node( ContactTopology* );
   void Connect_Node( ContactNode<Real>* );

   virtual void Connect_Entity( ContactTopology* ) = 0;

   virtual void Get_Face_Normal(Real *return_normal, 
                                ContactTopology* search_topology) = 0;
   virtual void Get_Avg_Face_Normal(Real *return_normal, 
                                    ContactTopology* search_topology, 
                                    VariableHandle NODE_COORDS) = 0;

   inline Real& Scalar_Var( VariableHandle vh );
   inline Real* Vector_Var( VariableHandle vh );

   virtual void Delete_Frag(ContactFixedSizeAllocator* allocators) = 0;

   inline Real*             Get_Coordinates();
   inline Real*             Get_Contact_Point();
   inline Real*             Get_Physical_Face_Normal();
   inline Real*             Get_Normal();
   inline Real*             Get_Pushback_Dir();
   inline Real&             Get_Gap_Cur();
   inline Real&             Get_Gap_Old();
   inline Real&             Get_Gap_Initial();
   inline Real&             Get_Node_Area();
   inline Real&             Get_Time_To_Contact();
   inline InteractionSource Get_Source();
   inline void              Set_Source(InteractionSource);
   
   CString Source_Name();

   inline void Set_Initial_Gap();

   virtual void Update_Tied_Interaction( VariableHandle, VariableHandle,  Real,
                                         ScratchVariable&, ScratchVariable&,
                                         ScratchVariable&, ScratchVariable& ) = 0;
   virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real, VariableHandle ) = 0;
   virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real, Real* ) = 0;
   virtual void Update_Tied_Interaction( VariableHandle, Real*, Real, Real* ) = 0;

   // Parallel pack/unpack functions
   virtual int   Size() = 0;
   virtual char* Pack( char* buffer ) = 0;
   virtual char* Unpack( char* buffer ) = 0;
   virtual int   Size_ForSecondary() = 0;
   virtual char* Pack_ForSecondary( char* buffer ) = 0;
   virtual char* Unpack_ForSecondary( char* buffer ) = 0;

   // Restart Pack/Unpack functions
   virtual int  Restart_Size() = 0;
   virtual void Restart_Pack( Real* buffer ) = 0;
   virtual void Restart_Unpack( Real* buffer ) = 0;

   InteractionType Get_Type(){return type;};

   //virtual void Add_Contact_Force_To_Entity(Real m, Real *a_sn) = 0;

 protected:
   bool is_tied;
   bool is_infSlip;
   bool is_glued;
   bool is_tracked;
   //bool track_initial_gap;
   int location;
   Real DataArray[NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS];
   ContactNode<Real>* node;
   ContactTopologyEntity<Real>* entity;
   ContactInteractionEntity<Real>::entity_data node_entity_data;
   ContactInteractionEntity<Real>::entity_data entity_entity_data;
   InteractionType type;
};

//------------------------------------------------------------------------------

inline int ContactNodeEntityInteraction::DataArray_Length()
{
  return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;
}

inline void ContactNodeEntityInteraction::Initialize_Memory()
{
  std::memset(DataArray_Buffer(), 0, DataArray_Length()*sizeof(Real));
}

inline ContactTopologyEntity<Real>* ContactNodeEntityInteraction::Entity(){
  return entity;
}

inline int ContactNodeEntityInteraction::Get_Node_Key() {
  return (int)Scalar_Var(NODE_ENTITY_KEY);
}

int ContactNodeEntityInteraction::Get_Entity_Key() {
  return entity->Entity_Key();
}

inline const ContactInteractionEntity<Real>::entity_data* ContactNodeEntityInteraction::Get_Entity_Data()
{
  return &entity_entity_data;
}

inline int ContactNodeEntityInteraction::Set_Entity_Data() 
{
  if (entity) {
    ContactInteractionEntity<Real>::SetEntityData(&entity_entity_data, entity);
    return 1;
  }else{      
    return 0;
  }
}

inline int ContactNodeEntityInteraction::Get_Entity_Owner(){
  return entity_entity_data.owner;
}

inline ContactNode<Real>* ContactNodeEntityInteraction::Node()
{
  return node;
}

inline const ContactInteractionEntity<Real>::entity_data* ContactNodeEntityInteraction::NodeEntityData()
{
  return &node_entity_data;
}

inline int ContactNodeEntityInteraction::Set_NodeEntityData() 
{
  if (node) {
    ContactInteractionEntity<Real>::SetEntityData(&node_entity_data, node);
    return 1;
  }else{      
    return 0;
  }
}

inline bool ContactNodeEntityInteraction::Is_Tied()
{
  return is_tied; 
}

inline void ContactNodeEntityInteraction::Is_Tied( bool IS_TIED )
{
  REMEMBER(if (IS_TIED) PRECONDITION(!(is_glued||is_infSlip)));
  is_tied = IS_TIED; 
}

inline bool ContactNodeEntityInteraction::Is_InfSlip()
{
  return is_infSlip; 
}

inline void ContactNodeEntityInteraction::Is_InfSlip(bool IS_INFSLIP)
{
  REMEMBER(if (IS_INFSLIP) PRECONDITION(!(is_glued||is_infSlip)));
  is_infSlip = IS_INFSLIP;
}

inline bool ContactNodeEntityInteraction::Is_Tracked()
{
  return is_tracked;
}

inline void ContactNodeEntityInteraction::Is_Tracked( bool IS_TRACKED )
{
  is_tracked = IS_TRACKED;
}

inline bool ContactNodeEntityInteraction::Is_Glued()
{
  return is_glued; 
}

inline void ContactNodeEntityInteraction::Is_Glued( bool IS_GLUED )
{
  REMEMBER(if (IS_GLUED) PRECONDITION(!(is_tied||is_infSlip)));
  is_glued = IS_GLUED; 
}

inline int ContactNodeEntityInteraction::Location()
{
  return location;
}

inline Real& ContactNodeEntityInteraction::Scalar_Var( VariableHandle vh )
{
  return DataArray[vh];
}

inline Real* ContactNodeEntityInteraction::Vector_Var( VariableHandle vh )
{
  return (DataArray+NUMBER_SCALAR_VARS+3*vh);
}

inline Real* ContactNodeEntityInteraction::Get_Coordinates()
{
  return Vector_Var(COORDINATES);
}

inline Real* ContactNodeEntityInteraction::Get_Contact_Point()
{
  return Vector_Var(CONTACT_POINT);
}

inline Real* ContactNodeEntityInteraction::Get_Physical_Face_Normal()
{
  return Vector_Var(PHYSICAL_FACE_NORMAL);
}

inline Real* ContactNodeEntityInteraction::Get_Normal()
{
  return Vector_Var(NORMAL_DIR);
}

inline Real* ContactNodeEntityInteraction::Get_Pushback_Dir()
{
  return Vector_Var(PUSHBACK_DIR);
}

inline Real& ContactNodeEntityInteraction::Get_Gap_Cur()
{
  return Scalar_Var(GAP_CUR);
}

inline Real& ContactNodeEntityInteraction::Get_Gap_Old()
{
  return Scalar_Var(GAP_OLD);
}

inline Real& ContactNodeEntityInteraction::Get_Gap_Initial()
{
  return Scalar_Var(GAP_INT);
}

inline Real& ContactNodeEntityInteraction::Get_Node_Area()
{
  return Scalar_Var(NODE_AREA);
}

inline Real& ContactNodeEntityInteraction::Get_Time_To_Contact()
{
  return Scalar_Var(TIME_TO_CONTACT);
}

inline ContactNodeEntityInteraction::InteractionSource ContactNodeEntityInteraction::Get_Source()
{
  return (ContactNodeEntityInteraction::InteractionSource)Scalar_Var(SOURCE);
}

inline void ContactNodeEntityInteraction::Set_Source(ContactNodeEntityInteraction::InteractionSource new_source)
{
  Scalar_Var(SOURCE) = new_source;
}

inline void ContactNodeEntityInteraction::Set_Initial_Gap(void)
{
  Scalar_Var(GAP_INT) = Scalar_Var(GAP_CUR);
}

#endif
