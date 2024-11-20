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


#ifndef ContactElementElementInteraction_h_
#define ContactElementElementInteraction_h_

#include "Contact_Defines.h"
#include "ContactInteractionEntity.h"
#include "ContactElement.h"
#include "ContactSearch.h"

class CString;
class ContactTopologyEntityList;
class ContactTopologyEntityHash;
class ContactGlobalID;
class ContactHostGlobalID;
class ContactFixedSizeAllocator;

class ContactElementElementInteraction : public ContactInteractionEntity<Real> {
  
 public:
  
  enum InteractionSource { UNKNOWN_SOURCE=-1,CLOSEST_POINT_PROJECTION_1=1, 
                           CLOSEST_POINT_PROJECTION_2, MOVING_INTERSECTION };

  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  EEI_SCALAR_VAR
#define EEI_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
                     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactElementElementInteraction();
  ContactElementElementInteraction( ContactElement*, ContactElement*, Real );
  ContactElementElementInteraction( ContactElementElementInteraction& );
  static ContactElementElementInteraction* new_ContactElementElementInteraction(
            ContactFixedSizeAllocator&, ContactElement*, ContactElement*, Real );
  static ContactElementElementInteraction* new_ContactElementElementInteraction(
	     ContactFixedSizeAllocator& );
  static ContactElementElementInteraction* new_ContactElementElementInteraction(
	     ContactFixedSizeAllocator&, ContactElementElementInteraction& );
  ~ContactElementElementInteraction();
  
#ifndef CONTACT_NO_MPI
  inline ContactZoltanLID& Zoltan_LID() { return zoltan_lid; };
  inline ContactZoltanGID& Zoltan_GID() { return zoltan_gid; };
  inline void ZoltanElementLID(LB_ID_PTR lid, int flag=0) 
    {//if (flag) zoltan_lid.ZoltanLID(CT_ELEMENT, 
     //                               master_element_entity_data.index_in_owner_proc_array, 
     //                               lid);
     //else      zoltan_lid.ZoltanLID(CT_ELEMENT, 
     //                               master_element_entity_data.index_in_proc_array, 
     //                               lid);
     zoltan_lid.ZoltanLID(CT_ELEMENT, 
                                    master_element_entity_data.index_in_owner_proc_array, 
                                    lid);};
  inline void ZoltanElementGID(LB_ID_PTR gid) 
    {zoltan_gid.ZoltanGID(CT_ELEMENT, 
                          master_element_entity_data.host_gid[0],
                          master_element_entity_data.host_gid[1],
                          gid);};
#endif

  inline ContactElement* SlaveElement() {return slave_element;};
  inline entity_data* SlaveElementEntityData() {return &slave_element_entity_data;};
  inline ContactElement* MasterElement() {return master_element;};
  inline entity_data* MasterElementEntityData() {return &master_element_entity_data;};
  int Set_SlaveElementEntityData() ;
  int Set_MasterElementEntityData() ;

  inline int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline void Initialize_Memory() {std::memset(DataArray_Buffer(), 0, DataArray_Length()*sizeof(Real));};
  inline Real& Scalar_Var( VariableHandle vh ) {return DataArray[vh];};
  inline Real* Vector_Var( VariableHandle vh ) 
    { return (DataArray+NUMBER_SCALAR_VARS+3*vh); };

  void Connect_SlaveElement ( ContactTopologyEntityList& );
  void Connect_MasterElement( ContactTopologyEntityList& );
  void Connect_SlaveElement ( ContactTopologyEntityHash& );
  void Connect_MasterElement( ContactTopologyEntityHash& );
  void Connect_SlaveElement ( ContactTopology* );
  void Connect_MasterElement( ContactTopology* );
  void Connect_SlaveElement ( ContactElement* );
  void Connect_MasterElement( ContactElement* );

  // Parallel packing/unpacking functions
  int  Size();
  void Pack( char* buffer );
  void Unpack( char* buffer );
  void Copy( ContactElementElementInteraction* src );

  // Restart Pack/Unpack functions
  int  Restart_Size();
  void Restart_Pack( Real* buffer );
  void Restart_Unpack( Real* buffer );

  int Data_Size();
  
 protected:

 private:
  Real DataArray[NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS+1];
  ContactElement* slave_element;
  entity_data slave_element_entity_data;
  ContactElement* master_element;
  entity_data master_element_entity_data;
#ifndef CONTACT_NO_MPI
  ContactZoltanLID zoltan_lid;
  ContactZoltanGID zoltan_gid;
#endif
};

#endif // _ContactElementElementInteraction_h_
