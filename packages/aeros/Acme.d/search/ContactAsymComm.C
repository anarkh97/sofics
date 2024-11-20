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


#ifndef CONTACT_NO_MPI

#include "ContactAsymComm.h"
#include "ContactSymComm.h"
#include "ContactZoltanComm.h"
#include "ContactZoltanCommUtils.h"
#include "ContactZoltanID.h"
#include "contact_assert.h"
#include "ContactTopologyEntity.h"
#include "ContactTopologyEntityHash.h"
#include "ContactTopologyEntityList.h"
#include "ContactParOStream.h"

#include <cstring>

using namespace std;

//*******************************************************************************************************
//
//  Sort all export and import lists by Global_ID.  The sorting allows merging of lists and should
//  improve cache coherence when looping over the lists.
//
//*******************************************************************************************************
void ContactAsymComm::Sort_Comm_Lists() {
  export_list.Sort();
  import_list.Sort();
}

//*******************************************************************************************************
//
//  Build an asymmetric communication specification by composing together
//  two other communication specifications
//
//  Assumptions:
//    If a node is a ghost node in one communication specification, then it is also a ghost node
//    in the other communication specification (or does not exist in the other.)
//
//*******************************************************************************************************
ContactAsymComm::ContactAsymComm( ContactAsymComm &com1,
                                  ContactAsymComm &com2) {
  //
  //  Sort the import and export lists by global_id
  //
  com1.Sort_Comm_Lists();
  com2.Sort_Comm_Lists();
  //
  //  Merge the communcation lists to create the new list
  //
  import_list.CreateFromListAddition(com1.import_list, com2.import_list);
  export_list.CreateFromListAddition(com1.export_list, com2.export_list);
}

ContactAsymComm::ContactAsymComm( ContactZoltanComm& comm, 
                                  ContactTopologyEntityList& lookup_list )
{
  // Extract the information from the ContactZoltanComm object
  int num_import        = comm.Num_Import();
  int num_export        = comm.Num_Export();
  LB_ID_PTR import_gids = comm.Import_GIDS();
  LB_ID_PTR export_gids = comm.Export_GIDS();
  LB_ID_PTR import_lids = comm.Import_LIDS();
  LB_ID_PTR export_lids = comm.Export_LIDS();
  int* import_procs     = comm.Import_Procs();
  int* export_procs     = comm.Export_Procs();

  import_list.CreateFromIDs(num_import, import_gids, import_lids, import_procs, lookup_list);
  export_list.CreateFromIDs(num_export, export_gids, export_lids, export_procs, lookup_list);
}

ContactAsymComm::ContactAsymComm( ContactZoltanComm& comm, 
				  ContactTopologyEntityHash& entity_hash )
{
  // Extract the information from the ContactZoltanComm object
  int num_import = comm.Num_Import();
  int num_export = comm.Num_Export();
  LB_ID_PTR import_gids = comm.Import_GIDS();
  LB_ID_PTR export_gids = comm.Export_GIDS();
  LB_ID_PTR import_lids = comm.Import_LIDS();
  LB_ID_PTR export_lids = comm.Export_LIDS();
  int* import_procs = comm.Import_Procs();
  int* export_procs = comm.Export_Procs();

  import_list.CreateFromIDs(num_import, import_gids, import_lids, import_procs, entity_hash);
  export_list.CreateFromIDs(num_export, export_gids, export_lids, export_procs, entity_hash);
}

ContactAsymComm::ContactAsymComm( ContactZoltanComm& comm, 
                                  Zoltan_Struct* zz,
                                  ContactTopologyEntityList& lookup_list )
{
  // Extract the information from the ContactZoltanComm object
  int dir               = -1;
  int num_import        = comm.Num_Import();
  int num_export        = comm.Num_Export();
  LB_ID_PTR import_gids = comm.Import_GIDS();
  LB_ID_PTR export_gids = comm.Export_GIDS();
  LB_ID_PTR import_lids = comm.Import_LIDS();
  LB_ID_PTR export_lids = comm.Export_LIDS();
  int* import_procs     = comm.Import_Procs();
  int* export_procs     = comm.Export_Procs();
  int gid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_GID_SIZE);
  //int lid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_LID_SIZE);
  if (num_export==-1) {
    dir          = 0;
    num_export   = 0;
    export_gids  = NULL;
    export_lids  = NULL;
    export_procs = NULL;
    int msgtag   = 15;
    ZOLTAN_COMM_OBJ *comm_plan=NULL;
    //  Create zoltan comm plan based on import lists.
    Zoltan_Comm_Create(&comm_plan, 
                       num_import, 
                       import_procs, 
                       zz->Communicator, 
                       msgtag, 
                       &num_export);
    if (num_export>0) {
      export_gids  = new ZOLTAN_ID_TYPE [num_export * ZOLTAN_GID_SIZE];
      export_lids  = new ZOLTAN_ID_TYPE [num_export * ZOLTAN_LID_SIZE];
      export_procs = new int [num_export];
    }
    Zoltan_Comm_Do(comm_plan, 
                   msgtag+1, 
                   (char *) import_gids,
                   gid_size,
                   (char *) export_gids);
    //Zoltan_Comm_Do(comm_plan, 
    //               msgtag+2, 
    //               (char *) import_lids,
    //               lid_size,
    //               (char *) export_lids);
    Zoltan_Comm_Info(comm_plan, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, export_procs, NULL);
    if (comm_plan != NULL) Zoltan_Comm_Destroy(&comm_plan);
  
  } else if (num_import==-1) {
    dir          = 1;
    num_import   = 0;
    import_gids  = NULL;
    import_lids  = NULL;
    import_procs = NULL;
    int msgtag   = 15;
    ZOLTAN_COMM_OBJ *comm_plan=NULL;
    //  Create zoltan comm plan based on export lists.
    Zoltan_Comm_Create(&comm_plan, 
                       num_export, 
                       export_procs, 
                       zz->Communicator, 
                       msgtag, 
                       &num_import);
    if (num_import>0) {
      import_gids  = new ZOLTAN_ID_TYPE [num_import * ZOLTAN_GID_SIZE];
      import_lids  = new ZOLTAN_ID_TYPE [num_import * ZOLTAN_LID_SIZE];
      import_procs = new int [num_export];
    }
    Zoltan_Comm_Do(comm_plan, 
                   msgtag+1, 
                   (char *) export_gids,
                   gid_size,
                   (char *) import_gids);
    //Zoltan_Comm_Do(comm_plan, 
    //               msgtag+2, 
    //               (char *) export_lids,
    //               lid_size,
    //               (char *) import_lids);
    Zoltan_Comm_Info(comm_plan, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, import_procs, NULL);
    if (comm_plan != NULL) Zoltan_Comm_Destroy(&comm_plan);
  }
  import_list.CreateFromIDs(num_import, import_gids, import_lids, import_procs, lookup_list);
  export_list.CreateFromIDs(num_export, export_gids, export_lids, export_procs, lookup_list);
  if (dir==0) {
    delete [] export_gids;
    delete [] export_lids;
    delete [] export_procs;
  }
  if (dir==1) {
    delete [] import_gids;
    delete [] import_lids;
    delete [] import_procs;
  }
}

ContactAsymComm::ContactAsymComm( ContactZoltanComm& comm, 
				  ContactTopologyEntityHash& entity_hash,
                                  ContactTopologyEntityList& lookup_list )
{
  // Extract the information from the ContactZoltanComm object
  int num_import        = comm.Num_Import();
  int num_export        = comm.Num_Export();
  LB_ID_PTR import_gids = comm.Import_GIDS();
  LB_ID_PTR export_gids = comm.Export_GIDS();
  LB_ID_PTR import_lids = comm.Import_LIDS();
  LB_ID_PTR export_lids = comm.Export_LIDS();
  int* import_procs     = comm.Import_Procs();
  int* export_procs     = comm.Export_Procs();
  
  if (num_import>0) {
    POSTCONDITION(import_gids);
    POSTCONDITION(import_lids);
    POSTCONDITION(import_procs);
  }
  if (num_export>0) {
    POSTCONDITION(export_gids);
    POSTCONDITION(export_lids);
    POSTCONDITION(export_procs);
  }

  import_list.CreateFromIDs(num_import, import_gids, import_lids, import_procs, entity_hash, lookup_list);
  export_list.CreateFromIDs(num_export, export_gids, export_lids, export_procs, entity_hash, lookup_list);
}

ContactAsymComm::ContactAsymComm( ContactZoltanCommUtils& comm, 
                                  ContactTopologyEntityList* lookup_list,
                                  ContactType base_type )
{
  // Extract the information from the ContactZoltanComm object
  int num_import = comm.Num_Import();
  int num_export = comm.Num_Export();
  LB_ID_PTR import_gids = comm.Import_GIDS();
  LB_ID_PTR export_gids = comm.Export_GIDS();
  LB_ID_PTR import_lids = comm.Import_LIDS();
  LB_ID_PTR export_lids = comm.Export_LIDS();
  int* import_procs = comm.Import_Procs();
  int* export_procs = comm.Export_Procs();

  import_list.CreateFromIDs(num_import, import_gids, import_lids, import_procs, lookup_list, base_type);
  export_list.CreateFromIDs(num_export, export_gids, export_lids, export_procs, lookup_list, base_type);
}

ContactAsymComm::ContactAsymComm( ContactSymComm& symcomm)
{
  //
  // This constructor builds a comm plan for pushing from the owner of an
  // entity to all the ghosts of this entity on other processors
  //
  
  int i,j;
  ContactTopologyEntity<Real>** entity_list;

  // Find sizing information
  for( i=0 ; i<symcomm.Num_Comm_Partners() ; ++i ){
    bool need_import = false;
    bool need_export = false;
    int proc_id = symcomm.Comm_Proc_ID( i );
    int num_to_proc = symcomm.Num_to_Proc( i );
    entity_list = symcomm.Entity_List( i );
    for( j=0 ; j<num_to_proc ; ++j ){
      if( entity_list[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ){
	need_export = true;
	export_list.num_entities += 1;
      } else if( entity_list[j]->Owner() == proc_id ){
	need_import = true;
	import_list.num_entities += 1;
      }
    }
    if( need_import ) import_list.num_comm_partners += 1;
    if( need_export ) export_list.num_comm_partners += 1;
  }

  // Allocate Memory
  export_list.Allocate();
  import_list.Allocate();
  
  // Fill Data Structures
  int import_entity_offset = 0;
  int export_entity_offset = 0;
  int import_proc_index = 0;
  int export_proc_index = 0;
  for( i=0 ; i<symcomm.Num_Comm_Partners() ; ++i ){
    int num_to_proc = symcomm.Num_to_Proc( i );
    int proc_id = symcomm.Comm_Proc_ID( i );
    entity_list = symcomm.Entity_List( i );
    int num_to_export = 0;
    int num_to_import = 0;
    for( j=0 ; j<num_to_proc ; ++j ){
      if( entity_list[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ){
	num_to_export += 1;
	export_list.entity_list[export_entity_offset++] = entity_list[j];
      } else if( entity_list[j]->Owner() == proc_id ){
	num_to_import += 1;
	import_list.entity_list[import_entity_offset++] = entity_list[j];
      }
    }
    if( num_to_import ){
      import_list.comm_proc_ids[import_proc_index] = proc_id;
      import_list.num_to_proc[import_proc_index] = num_to_import;
      ++import_proc_index;
    }
    if( num_to_export ){
      export_list.comm_proc_ids[export_proc_index] = proc_id;
      export_list.num_to_proc[export_proc_index] = num_to_export;
      ++export_proc_index;
    }
  }
 
  import_list.CalculateOffset();
  export_list.CalculateOffset();
}

ContactAsymComm::ContactAsymComm( const int Num_Export_Comm_Partners,
				  const int Num_Import_Comm_Partners,
				  const int* Num_Export_to_Proc,
				  const int* Num_Import_to_Proc,
				  const int* Export_Comm_Proc_IDs,
				  const int* Import_Comm_Proc_IDs,
				  ContactTopologyEntity<Real>** Export_Entity_List,
				  ContactTopologyEntity<Real>** Import_Entity_List ) :
  import_list(Num_Import_Comm_Partners, Num_Import_to_Proc, Import_Comm_Proc_IDs, Import_Entity_List),
  export_list(Num_Export_Comm_Partners, Num_Export_to_Proc, Export_Comm_Proc_IDs, Export_Entity_List)
{}

ContactAsymComm::ContactAsymComm( const int Num_Export_Comm_Partners,
				  const int Num_Import_Comm_Partners,
				  const int* Num_Export_to_Proc,
				  const int* Num_Import_to_Proc,
				  const int* Export_Comm_Proc_IDs,
				  const int* Import_Comm_Proc_IDs ) :
  import_list(Num_Import_Comm_Partners, Num_Import_to_Proc, Import_Comm_Proc_IDs),
  export_list(Num_Export_Comm_Partners, Num_Export_to_Proc, Export_Comm_Proc_IDs)
{}



ContactAsymComm::~ContactAsymComm()
{}

//
//  Print the current communication specification to the terminal
//
void ContactAsymComm::Print() {
  cout<<"ContactAsymComm:"<<endl;
  cout<<"  Export Comm List"<<endl;
  export_list.Print("    ");
  cout<<"  Import Comm List"<<endl;
  import_list.Print("    ");
}

void ContactAsymComm::Print(ContactParOStream& os) {
  os<<"ContactAsymComm\n";
  os<<"  Export Comm List\n";
  export_list.Print("    ",os);
  os<<"  Import Comm List\n";
  import_list.Print("    ",os);
}
void ContactAsymComm::Print(char* s) {
  cout<<"ContactAsymComm: "<<s<<endl;
  cout<<"  Export Comm List"<<endl;
  export_list.Print("    ");
  cout<<"  Import Comm List"<<endl;
  import_list.Print("    ");
}

void ContactAsymComm::Print(char* s, ContactParOStream& os) {
  os<<"ContactAsymComm: "<<s<<"\n";
  os<<"  Export Comm List\n";
  export_list.Print("    ",os);
  os<<"  Import Comm List\n";
  import_list.Print("    ",os);
}


#endif
