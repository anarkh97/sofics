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


#ifndef ContactSymComm_h_
#define ContactSymComm_h_

#include "Contact_Defines.h"
#include "contact_assert.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
template <typename DataType> class ContactTopologyEntity;
class ContactTopologyEntityHash;
class ContactTopologyEntityList;
class ContactHostGlobalID;

class ContactSymComm {
  
 public:
  ContactSymComm( int, const int*, const int*, ContactTopologyEntity<Real>** );
  ContactSymComm();
  ~ContactSymComm();
  
  inline int Size() {return num_entities;};
  inline int Num_Comm_Partners() {return num_comm_partners;};
  inline int Comm_Proc_ID( int i ){
    PRECONDITION( i<num_comm_partners || num_comm_partners==0 );
    return comm_proc_ids[i];
  };

  inline int Num_to_Proc( int i ){
    PRECONDITION( i<num_comm_partners || num_comm_partners==0 );
    return num_to_proc[i];
  }
  inline ContactTopologyEntity<Real>** Entity_List( int i ){
    PRECONDITION( i<num_comm_partners || num_comm_partners==0 );
    if (num_comm_partners==0)
      return NULL;
    else
      return entity_list+offset[i];
  }

  void Build_Comm_Plan( ContactSymComm*, int, int, int*, ContactHostGlobalID*,
			ContactTopologyEntityList*, MPI_Comm& );
  void Build_Comm_Plan( ContactSymComm*, int, int, int*, ContactHostGlobalID*,
			ContactTopologyEntityHash*, MPI_Comm& );
  void Build_Subset_Comm_Plan_Using_Temp_Tag( ContactSymComm& );

  void Check(MPI_Comm Comm);
  void Sort();

 private:
  
  ContactSymComm(ContactSymComm&);
  ContactSymComm& operator=(ContactSymComm&);

  int num_comm_partners;    // How many procs I communicate to
  int num_entities;         // Number of entities in entity_list
  int* comm_proc_ids;       // Processor ids (num_com_partners long)
  int* num_to_proc;         // Number to each proc (num_com_partners long)
  int* offset;              // Offset in entity_list for comm_partner i 
  ContactTopologyEntity<Real>** entity_list;

  int allocated_procs;
  int allocated_entities;

};

#endif
#endif
