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


#ifndef ContactCommList_h_
#define ContactCommList_h_

#ifndef CONTACT_NO_MPI
#include "zoltan.h"
#include "Contact_Defines.h"
#include "contact_assert.h"
#include "ContactTopologyEntity.h"

class ContactTopologyEntityHash;
class ContactTopologyEntityList;
class ContactZoltanComm;

class ContactCommList {
 public:
  ContactCommList();
  ContactCommList(const int Num_Comm_Partners,
                  const int* Num_to_Proc,
                  const int* Export_Proc_IDs,
                  ContactTopologyEntity<Real>** Entity_List);
  ContactCommList(const int  Num_Comm_Partners,
	       	  const int* Num_to_Proc,
		  const int* Comm_Proc_IDs);
  
  ~ContactCommList();

  void Sort();
  void CreateFromListAddition(const ContactCommList &list1, 
                              const ContactCommList &list2);

  void Allocate();

  void Print(const char*);

  void Print(const char*, ContactParOStream&);

  void CalculateOffset();

  void CreateFromIDs(const int num,
                     const LB_ID_PTR gids,
                     const LB_ID_PTR lids,
                     const int* procs,
                     ContactTopologyEntityHash& entity_hash);

  void CreateFromIDs(const int num,
                     const LB_ID_PTR gids,
                     const LB_ID_PTR lids,
                     const int* procs,
                     ContactTopologyEntityList& entity_list);

  void CreateFromIDs(const int num,
                     const LB_ID_PTR gids,
                     const LB_ID_PTR lids,
                     const int* procs,
                     ContactTopologyEntityHash& entity_hash,
                     ContactTopologyEntityList& entity_list);

  void CreateFromIDs(const int num,
                     const LB_ID_PTR gids,
                     const LB_ID_PTR lids,
                     const int* procs,
                     ContactTopologyEntityList* entity_list,
                     ContactType base_type);

  void Set_Index_From_EnfArrayIndex();

  friend class ContactAsymComm;
  friend class ContactSymComm;

 private:
 
  ContactCommList(ContactCommList&);
  ContactCommList& operator=(ContactCommList&);
  
  int  num_comm_partners;
  int  num_entities;
  int* comm_proc_ids;
  int* num_to_proc;
  int* offset;
  ContactTopologyEntity<Real>** entity_list;
  int* entity_index_list;
};

#endif
#endif
