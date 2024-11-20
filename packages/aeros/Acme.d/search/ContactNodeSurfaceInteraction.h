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


#ifndef ContactNodeSurfaceInteraction_h_
#define ContactNodeSurfaceInteraction_h_

#include "Contact_Defines.h"
#include "ContactInteractionEntity.h"
#include "ContactNode.h"
#include "ContactAnalyticSurface.h"
#include "ContactNodeEntityInteraction.h"
#include "ContactTopology.h"

class CString;
class ContactTopologyEntityList;
class ContactFixedSizeAllocator;
class ContactTopology;

class ContactNodeSurfaceInteraction : public ContactNodeEntityInteraction {

 public:

  ContactNodeSurfaceInteraction( InteractionSource, ContactNode<Real>*, 
				 ContactAnalyticSurface*,
                                 int, Real*, Real, Real*, Real*, Real *, Real, int, bool, bool, bool, bool);
  ContactNodeSurfaceInteraction();
  ContactNodeSurfaceInteraction(const ContactNodeSurfaceInteraction &cnsi);
  static ContactNodeSurfaceInteraction* new_ContactNodeSurfaceInteraction(
            ContactFixedSizeAllocator&, InteractionSource, ContactNode<Real>*,
            ContactAnalyticSurface*, int, Real*, Real, Real*, Real*, Real*, Real, int, bool, bool, bool, bool );
  static ContactNodeSurfaceInteraction* new_ContactNodeSurfaceInteraction(
             ContactFixedSizeAllocator& );     
  static ContactNodeSurfaceInteraction* new_ContactNodeSurfaceInteraction(
             ContactFixedSizeAllocator&, ContactNodeSurfaceInteraction& );     
  ~ContactNodeSurfaceInteraction();
  
  virtual ContactNodeEntityInteraction *New_Instance(ContactFixedSizeAllocator* allocators);
  
  inline ContactAnalyticSurface* Surface();
  inline int SurfaceID();

  // Parallel pack/unpack functions
  virtual int   Size();
  virtual char* Pack( char* buffer );
  virtual char* Unpack( char* buffer );
  
  virtual int   Size_ForSecondary();
  virtual char* Pack_ForSecondary( char* buffer );
  virtual char* Unpack_ForSecondary( char* buffer );

  // Restart Pack/Unpack functions
  virtual int  Restart_Size();
  virtual void Restart_Pack( Real* buffer );
  virtual void Restart_Unpack( Real* buffer );

  /**
   *  This function ignores the last five variable handles as they are not needed, however, they are part of the
   *  combined interface for a node entity interaction. 
   */
  virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real,
                                        ScratchVariable&, ScratchVariable&,
                                        ScratchVariable&, ScratchVariable& );
  
  virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real, VariableHandle );
  virtual void Update_Tied_Interaction( VariableHandle, VariableHandle, Real, Real* );
  virtual void Update_Tied_Interaction( VariableHandle, Real*, Real, Real* );

  virtual void Connect_Entity( ContactTopology* );

  inline void Connect_Surface( int, ContactAnalyticSurface** );
  inline void Connect_Surface( ContactTopology* );

  virtual void Get_Face_Normal(Real *return_normal, 
                       ContactTopology* search_topology);
  virtual void Get_Avg_Face_Normal(Real *return_normal, 
                           ContactTopology* search_topology,
                           VariableHandle NODE_COORDS);

  virtual void Delete_Frag(ContactFixedSizeAllocator* allocators);

  virtual int Get_Entity_Key();

  //void Add_Contact_Force_To_Entity(Real m, Real *a_sn) {};

 private:
//  ContactAnalyticSurface* surface;
  int surface_id;
  
};

//-------------------------------------------------------------------------------

inline ContactAnalyticSurface* ContactNodeSurfaceInteraction::Surface()
{
  return static_cast<ContactAnalyticSurface*>(entity);
}

inline int ContactNodeSurfaceInteraction::SurfaceID()
{
  return surface_id;
}

inline void ContactNodeSurfaceInteraction::Connect_Surface(int offset_for_id, 
					       ContactAnalyticSurface** surfs )
{
  entity = surfs[surface_id];
  POSTCONDITION( entity );

  REMEMBER(int surface_data_set = )
  Set_Entity_Data();
  POSTCONDITION( surface_data_set );
}

inline void ContactNodeSurfaceInteraction::Connect_Surface(ContactTopology* topology )
{
  entity = topology->Analytic_Surface(surface_id);
  POSTCONDITION( entity );

  REMEMBER(int surface_data_set = )
  Set_Entity_Data();
  POSTCONDITION( surface_data_set );
}

#endif
