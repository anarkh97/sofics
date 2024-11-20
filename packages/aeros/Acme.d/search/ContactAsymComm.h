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


#ifndef ContactAsymComm_h_
#define ContactAsymComm_h_

#ifndef CONTACT_NO_MPI
#include "zoltan.h"
#include "Contact_Defines.h"
#include "contact_assert.h"
#include "ContactTopologyEntity.h"
#include "ContactCommList.h"

class ContactTopologyEntityHash;
class ContactTopologyEntityList;
class ContactZoltanComm;
class ContactZoltanCommUtils;
class ContactSymComm;
class ContactParOStream;

class ContactAsymComm {
  
 public:
  ContactAsymComm( ContactZoltanComm&, ContactTopologyEntityHash& );
  ContactAsymComm( ContactZoltanComm&, ContactTopologyEntityHash&, ContactTopologyEntityList&);
  ContactAsymComm( ContactZoltanComm&, ContactTopologyEntityList& );
  ContactAsymComm( ContactZoltanCommUtils&, ContactTopologyEntityList*, ContactType );
  ContactAsymComm( ContactSymComm& );
  ContactAsymComm( ContactAsymComm &com1, ContactAsymComm &com2);
  ContactAsymComm( ContactZoltanComm& comm, 
                   Zoltan_Struct* z_ptr,
                   ContactTopologyEntityList& lookup_list );


  ContactAsymComm( const int Num_Export_Comm_Partners,
		   const int Num_Import_Comm_Partners,
		   const int* Num_Export_to_Proc,
		   const int* Num_Import_to_Proc,
		   const int* Export_Comm_Proc_IDs,
		   const int* Import_Comm_Proc_IDs,
		   ContactTopologyEntity<Real>** Export_Entity_List,
		   ContactTopologyEntity<Real>** Import_Entity_List );
  ContactAsymComm( const int Num_Export_Comm_Partners,
		   const int Num_Import_Comm_Partners,
		   const int* Num_Export_to_Proc,
		   const int* Num_Import_to_Proc,
		   const int* Export_Comm_Proc_IDs,
		   const int* Import_Comm_Proc_IDs );
  ~ContactAsymComm();

  void Sort_Comm_Lists();

  void Print();
  void Print(ContactParOStream&);
  void Print(char* s);
  void Print(char* s, ContactParOStream&);

  inline int Size_Export() { return export_list.num_entities; };
  inline int Size_Import() { return import_list.num_entities; };
  inline int Num_Export_Comm_Partners() { return export_list.num_comm_partners; };
  inline int Num_Import_Comm_Partners() { return import_list.num_comm_partners; };
  inline int* Export_Comm_Proc_IDs() { return export_list.comm_proc_ids; };
  inline int* Import_Comm_Proc_IDs() { return import_list.comm_proc_ids; };
  inline int* Num_Export_to_Procs() { return export_list.num_to_proc; };
  inline int* Num_Import_from_Procs() { return import_list.num_to_proc; };
  inline ContactTopologyEntity<Real>** Export_Entity_Lists() { return export_list.entity_list; };
  inline ContactTopologyEntity<Real>** Import_Entity_Lists() { return import_list.entity_list; };

  inline int* Export_Index_List() {return export_list.entity_index_list; };
  inline int* Import_Index_List() {return import_list.entity_index_list; };

  inline int Export_Comm_Proc_ID( int i ){
    PRECONDITION( i>=0 && i<export_list.num_comm_partners );
    return export_list.comm_proc_ids[i];
  }
  inline int Import_Comm_Proc_ID( int i ){
    PRECONDITION( i>=0 && i<import_list.num_comm_partners );
    return import_list.comm_proc_ids[i];
  }
  inline int Num_Export_to_Proc( int i ){
    PRECONDITION( i>=0 && i<export_list.num_comm_partners );
    return export_list.num_to_proc[i];
  }
  inline int Num_Import_from_Proc( int i ){
    PRECONDITION( i>=0 && i<import_list.num_comm_partners );
    return import_list.num_to_proc[i];
  }

  inline ContactTopologyEntity<Real>** Export_Entity_List( int i ){
    PRECONDITION( i>=0 && i<export_list.num_comm_partners );
    return export_list.entity_list+export_list.offset[i];
  }
  inline ContactTopologyEntity<Real>** Import_Entity_List( int i ){
    PRECONDITION( i>=0 && i<import_list.num_comm_partners );
    return import_list.entity_list+import_list.offset[i];
  }

  inline int* Export_Index_List( int i ){
    PRECONDITION( i>=0 && i<export_list.num_comm_partners );
    return export_list.entity_index_list+export_list.offset[i];
  }
  inline int* Import_Index_List( int i ){
    PRECONDITION( i>=0 && i<import_list.num_comm_partners );
    return import_list.entity_index_list+import_list.offset[i];
  }

  void Set_Index_From_EnfArrayIndex() {
    import_list.Set_Index_From_EnfArrayIndex();
    export_list.Set_Index_From_EnfArrayIndex();
  }

 private:
  ContactCommList import_list;
  ContactCommList export_list;
};
#endif
#endif
