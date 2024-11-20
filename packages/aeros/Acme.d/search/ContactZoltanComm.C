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


#include "ContactZoltanComm.h"
#include "ContactZoltanID.h"
#include "contact_assert.h"
#include "ContactEntity.h"
#include "ContactHostGlobalID.h"
#include "ContactParOStream.h"

#include <algorithm>
#include <cstring>
#include <iostream>

#ifndef CONTACT_NO_MPI

using namespace std;
	

ContactZoltanComm::ContactZoltanComm()
{
  direction = ZOLTAN_UNKNOWN;
  import_gids  = NULL;
  import_lids  = NULL;
  import_procs = NULL;
  export_gids  = NULL;
  export_lids  = NULL;
  export_procs = NULL;
  num_import   = 0;
  num_export   = 0;
  cur_import_capacity = 0;
  cur_export_capacity = 0;
}
   
ContactZoltanComm::ContactZoltanComm( Zoltan_Dir dir)
{
  direction = dir;
  import_gids  = NULL;
  import_lids  = NULL;
  import_procs = NULL;
  export_gids  = NULL;
  export_lids  = NULL;
  export_procs = NULL;
  num_import   = 0;
  num_export   = 0;
  cur_import_capacity = 0;
  cur_export_capacity = 0;
}


ContactZoltanComm::~ContactZoltanComm() {
  CleanUp();
}

void ContactZoltanComm::Initialize(Zoltan_Dir dir) {
  CleanUp();
  direction = dir;
}


void ContactZoltanComm::CleanUp() 
{
  switch (direction) {
  case ZOLTAN_IMPORT: {
    if(import_gids)  delete [] import_gids;
    if(import_gids)  delete [] import_lids;
    if(import_procs) delete [] import_procs;

    LB_ID_TYPE** export_gids_ptr  = (export_gids == NULL)  ? NULL : &export_gids ;
    LB_ID_TYPE** export_lids_ptr  = (export_lids == NULL)  ? NULL : &export_lids ;
    int**        export_procs_ptr = (export_procs == NULL) ? NULL : &export_procs;

    Zoltan_LB_Free_Data( NULL, NULL, NULL,
                  export_gids_ptr, 
                  export_lids_ptr, 
                  export_procs_ptr );
    break;
  }
  case ZOLTAN_EXPORT: {
    if(export_gids)  delete [] export_gids;
    if(export_lids)  delete [] export_lids;
    if(export_procs) delete [] export_procs;

    LB_ID_TYPE** import_gids_ptr  = (import_gids == NULL)  ? NULL : &import_gids ;
    LB_ID_TYPE** import_lids_ptr  = (import_lids == NULL)  ? NULL : &import_lids ;
    int**        import_procs_ptr = (import_procs == NULL) ? NULL : &import_procs;

    Zoltan_LB_Free_Data( import_gids_ptr, 
                  import_lids_ptr, 
                  import_procs_ptr, 
                  NULL, NULL, NULL );
    break;
  }
  default:
    break;
  }
  direction = ZOLTAN_UNKNOWN;
  import_gids  = NULL;
  import_lids  = NULL;
  import_procs = NULL;
  export_gids  = NULL;
  export_lids  = NULL;
  export_procs = NULL;
  num_import   = 0;
  num_export   = 0;
  cur_import_capacity = 0;
  cur_export_capacity = 0;
  lookup.clear();
}

void ContactZoltanComm::Add_Import(LB_ID_PTR lid, LB_ID_PTR gid, int proc, int check)
{
  // If the node is in the list, just return

  if(check) {
    if (find_hash_entry(gid, proc)) return;
    add_hash_entry(gid, proc);
  }
  static int nbytes_lid  = ZOLTAN_LID_SIZE*sizeof(LB_ID_TYPE);
  static int nbytes_gid  = ZOLTAN_GID_SIZE*sizeof(LB_ID_TYPE);
  static int nbytes_proc = sizeof(int);

  if(num_import + 1 > cur_import_capacity) {
    cur_import_capacity = max((cur_import_capacity*3)/2, 16);

    LB_ID_PTR new_import_gids  = new LB_ID_TYPE[cur_import_capacity*ZOLTAN_GID_SIZE];
    LB_ID_PTR new_import_lids  = new LB_ID_TYPE[cur_import_capacity*ZOLTAN_LID_SIZE];
    int*      new_import_procs = new int[cur_import_capacity];

    std::memcpy(new_import_gids,  import_gids,  nbytes_gid*num_import);
    std::memcpy(new_import_lids,  import_lids,  nbytes_lid*num_import);
    std::memcpy(new_import_procs, import_procs, nbytes_proc*num_import);
    if(import_gids) delete [] import_gids;
    if(import_lids) delete [] import_lids;
    if(import_procs) delete [] import_procs;
    import_gids = new_import_gids;
    import_lids = new_import_lids;
    import_procs = new_import_procs;
  }

  std::memcpy(&import_lids[num_import*ZOLTAN_LID_SIZE],lid,nbytes_lid);
  std::memcpy(&import_gids[num_import*ZOLTAN_GID_SIZE],gid,nbytes_gid);
  import_procs[num_import] = proc;
  ++num_import;
}

void ContactZoltanComm::Set_Export(int num, LB_ID_PTR gids, LB_ID_PTR lids, int* procs) {
  cur_export_capacity = num;
  num_export = num;
  if(export_gids) delete [] export_gids;
  if(export_lids) delete [] export_lids;
  if(export_procs) delete [] export_procs;
  export_gids = gids;
  export_lids = lids;
  export_procs = procs;
}


void ContactZoltanComm::Add_Export(LB_ID_PTR lid, LB_ID_PTR gid, 
                                   int proc, int check)
{  
  if (check) {
    if (find_hash_entry(gid, proc)) return;
    add_hash_entry(gid, proc);
  }

  static int nbytes_lid  = ZOLTAN_LID_SIZE*sizeof(LB_ID_TYPE);
  static int nbytes_gid  = ZOLTAN_GID_SIZE*sizeof(LB_ID_TYPE);
  static int nbytes_proc = sizeof(int);

  if(num_export + 1 > cur_export_capacity) {
    cur_export_capacity = max((cur_export_capacity*3)/2, 16);

    LB_ID_PTR new_export_gids  = new LB_ID_TYPE[cur_export_capacity*ZOLTAN_GID_SIZE];
    LB_ID_PTR new_export_lids  = new LB_ID_TYPE[cur_export_capacity*ZOLTAN_LID_SIZE];
    int*      new_export_procs = new int[cur_export_capacity];

    std::memcpy(new_export_gids,  export_gids,  nbytes_gid*num_export);
    std::memcpy(new_export_lids,  export_lids,  nbytes_lid*num_export);
    std::memcpy(new_export_procs, export_procs, nbytes_proc*num_export);
    if(export_gids) delete [] export_gids;
    if(export_lids) delete [] export_lids;
    if(export_procs) delete [] export_procs;
    export_gids = new_export_gids;
    export_lids = new_export_lids;
    export_procs = new_export_procs;
  }

  std::memcpy(&export_lids[num_export*ZOLTAN_LID_SIZE],lid,nbytes_lid);
  std::memcpy(&export_gids[num_export*ZOLTAN_GID_SIZE],gid,nbytes_gid);
  export_procs[num_export] = proc;
  ++num_export;
}

void ContactZoltanComm::add_hash_entry(  LB_ID_PTR gid, int proc )
{
  std::pair<ContactZoltanGID, int> key = pair<ContactZoltanGID, int>(ContactZoltanGID(gid), proc);
  lookup[key] = 1;
}




int ContactZoltanComm::find_hash_entry( LB_ID_PTR gid, int proc )
{
  std::pair<ContactZoltanGID, int> key = pair<ContactZoltanGID, int>(ContactZoltanGID(gid), proc);
  if(lookup.find(key) != lookup.end()) {
    return 1;
  } else {
    return 0;
  }
}

void ContactZoltanComm::Print(ContactParOStream& os)
{
  if (export_gids!=NULL && export_lids!=NULL && export_procs!=NULL) {
    os << "ContactZoltanComm EXPORTS ("<<num_export<<"):\n";
    for (int i=0; i< num_export; ++i){
      ContactHostGlobalID GlobalID(&export_gids[i*ZOLTAN_GID_SIZE]);
      int entity_type = ContactZoltanGID::Type(&export_gids[i*ZOLTAN_GID_SIZE]);
      int LocalID = ContactZoltanLID::Index(&export_lids[i*ZOLTAN_LID_SIZE]);
      POSTCONDITION(entity_type==ContactZoltanLID::Type(&export_lids[i*ZOLTAN_LID_SIZE]));
      os << "  " << i << ": ";
      switch (entity_type) {
      case CT_NODE:
        os<<"NODE";
        break;
      case CT_SHELL_NODE:
        os<<"SHELL_NODE";
        break;
      case CT_EDGE:
        os<<"EDGE";
        break;
      case CT_FACE:
        os<<"FACE";
        break;
      case CT_ELEMENT:
        os<<"ELEMENT";
        break;
      default:
        POSTCONDITION(0);
        break;
      }
      os<<",  gid = "<<GlobalID<<",  index = "<<LocalID<<",  proc = "<<export_procs[i]<<"\n";
    }
  }
  if (import_gids!=NULL && import_lids!=NULL && import_procs!=NULL) {
    os << "ContactZoltanComm IMPORTS ("<<num_import<<"):\n";
    for (int i=0; i< num_import; ++i){
      ContactHostGlobalID GlobalID(&import_gids[i*ZOLTAN_GID_SIZE]);
      int entity_type = ContactZoltanGID::Type(&import_gids[i*ZOLTAN_GID_SIZE]);
      int LocalID = ContactZoltanLID::Index(&import_lids[i*ZOLTAN_LID_SIZE]);
      POSTCONDITION(entity_type==ContactZoltanLID::Type(&import_lids[i*ZOLTAN_LID_SIZE]));
      os << "  " << i << ": ";
      switch (entity_type) {
      case CT_NODE:
        os<<"NODE";
        break;
      case CT_SHELL_NODE:
        os<<"SHELL_NODE";
        break;
      case CT_EDGE:
        os<<"EDGE";
        break;
      case CT_FACE:
        os<<"FACE";
        break;
      case CT_ELEMENT:
        os<<"ELEMENT";
        break;
      default:
        POSTCONDITION(0);
        break;
      }
      os<<",  gid = "<<GlobalID<<",  index = "<<LocalID<<",  proc = "<<import_procs[i]<<"\n";
    }
  }
}



#endif

