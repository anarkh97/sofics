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

#include "ContactZoltanCommUtils.h"
#include "ContactHostGlobalID.h"
#include "ContactZoltanID.h"
#include <cstdlib>
#include <cstring>

ContactZoltanCommUtils::ContactZoltanCommUtils(CommDirection dir, 
                                               Zoltan_Struct* z_ptr,
                                               int num, 
                                               LB_ID_TYPE* gids, 
                                               LB_ID_TYPE* lids, 
                                               int* pids,
                                               ContactParOStream&  os)
{
  type         = 0;
  zz           = z_ptr;
  direction    = dir;
  import_gids  = NULL;
  import_lids  = NULL;
  import_pids  = NULL;
  export_gids  = NULL;
  export_lids  = NULL;
  export_pids  = NULL;
  send_buffer  = NULL;
  recv_buffer  = NULL;
  send_sizes   = NULL;
  recv_sizes   = NULL;
  send_indices = NULL;
  recv_indices = NULL;
  import_plan  = NULL;
  export_plan  = NULL;
  int msgtag   = 15;
  int err      = 0;
  int gid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_GID_SIZE);
  int lid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_LID_SIZE);
  switch (dir) {
  case IMPORT:
    num_import  = num;
    num_export  = 0;
    //  Create zoltan comm plan based on import lists.
    err = Zoltan_Comm_Create(&import_plan, 
                             num_import, 
                             pids, 
                             zz->Communicator, 
                             msgtag, 
                             &num_export);
    //  Allocate space for export buffers.
    if (num_export>0) {
      export_gids = new ZOLTAN_ID_TYPE [num_export * ZOLTAN_GID_SIZE];
      export_lids = new ZOLTAN_ID_TYPE [num_export * ZOLTAN_LID_SIZE];
      export_pids = new int [num_export];
    }
    
    //  Get the GIDs of the export objects. 
    err = Zoltan_Comm_Do(import_plan, 
                         msgtag+1, 
                         (char *) gids,
                         gid_size,
                         (char *) export_gids);

    //  Get the LIDs of the export objects.
    err = Zoltan_Comm_Do(import_plan, 
                         msgtag+2, 
                         (char *) lids,
                         lid_size,
                         (char *) export_lids);

    //  Get the export processors.
    Zoltan_Comm_Info(import_plan, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, export_pids, NULL);

    //  Create plan based on exports so we can set variable sizes.
    err = Zoltan_Comm_Invert_Plan(&import_plan);
    export_plan = import_plan;
    import_plan = NULL;
    break;
  case EXPORT:
    num_export  = num;
    export_gids = gids;
    export_lids = lids;
    export_pids = pids;
    num_import  = 0;
    //  Create zoltan comm plan based on export lists.
    err = Zoltan_Comm_Create(&export_plan, 
                             num_export, 
                             export_pids, 
                             zz->Communicator, 
                             msgtag, 
                             &num_import);
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  
  //  Allocate space for export buffers.
  if (num_export>0) {
    send_sizes   = new int [num_export];
    send_indices = new int [num_export];
  }
  
  //  Allocate space for import buffers.
  if (num_import>0) {
    recv_sizes   = new int [num_import];
    recv_indices = new int [num_import];
    import_gids  = new ZOLTAN_ID_TYPE [num_import * ZOLTAN_GID_SIZE];
  }
  
  zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data, 
                         ZOLTAN_GID_SIZE, 
                         ZOLTAN_LID_SIZE, 
                         num_export,
                         export_gids,
                         export_lids,
                         send_sizes,
                         &err);
                                    
  int tagsize = gid_size + sizeof(int);
  
  int send_size = 0;
  for (int i = 0; i < num_export; ++i) {
    send_size += send_sizes[i]+tagsize;
  }
  if (send_size>0) send_buffer = new char [send_size];

  // Modify sizes[] to contain message sizes, not object sizes
  for (int i=0; i<num_export; ++i) {
    send_sizes[i] += tagsize;
  }
  
  //  Resize the communication plan.
  int recv_size=0;
  err = Zoltan_Comm_Resize(export_plan, send_sizes, msgtag+3, &recv_size);
  if (recv_size>0) recv_buffer = new char [recv_size];
  
  // Modify sizes[] to contain object sizes, not message sizes
  for (int i=0; i<num_export; ++i) {
    send_sizes[i] -= tagsize;
  }
  
}

ContactZoltanCommUtils::ContactZoltanCommUtils(ContactZoltanComm* comm1,
                                               ContactZoltanComm* comm2,
                                               Zoltan_Struct* z_ptr,
                                               ContactParOStream&  os)
{
  type         = 1;
  int err      = 0;
  int msgtag   = 15;
  zz           = z_ptr;
  direction    = BOTH;
  import_gids  = NULL;
  import_lids  = NULL;
  import_pids  = NULL;
  export_gids  = NULL;
  export_lids  = NULL;
  export_pids  = NULL;
  send_buffer  = NULL;
  recv_buffer  = NULL;
  send_sizes   = NULL;
  recv_sizes   = NULL;
  send_indices = NULL;
  recv_indices = NULL;
  import_plan  = NULL;
  export_plan  = NULL;
  int gid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_GID_SIZE);
  int lid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_LID_SIZE);
  
  ZOLTAN_COMM_OBJ* import_plan1;
  int       import_num1  = 0;
  LB_ID_PTR import_gids1 = NULL;
  LB_ID_PTR import_lids1 = NULL;
  int*      import_pids1 = NULL;
  int       export_num1  = 0;
  LB_ID_PTR export_gids1 = NULL;
  LB_ID_PTR export_lids1 = NULL;
  int*      export_pids1 = NULL;
  
  ZOLTAN_COMM_OBJ* import_plan2;
  int       import_num2  = 0;
  LB_ID_PTR import_gids2 = NULL;
  LB_ID_PTR import_lids2 = NULL;
  int*      import_pids2 = NULL;
  int       export_num2  = 0;
  LB_ID_PTR export_gids2 = NULL;
  LB_ID_PTR export_lids2 = NULL;
  int*      export_pids2 = NULL;

  bool export_alloc = false;
  
  switch (comm1->Direction()) {
  case ContactZoltanComm::ZOLTAN_IMPORT:
    import_num1  = comm1->Num_Import();
    import_gids1 = comm1->Import_GIDS();
    import_lids1 = comm1->Import_LIDS();
    import_pids1 = comm1->Import_Procs();
    //  Create zoltan comm plan based on import lists.
    err = Zoltan_Comm_Create(&import_plan1, 
                             import_num1, 
                             import_pids1, 
                             zz->Communicator, 
                             msgtag, 
                             &export_num1);
    //  Allocate space for export buffers.
    if (export_num1>0) {
      export_alloc = true;
      export_gids1 = new ZOLTAN_ID_TYPE [export_num1 * ZOLTAN_GID_SIZE];
      export_lids1 = new ZOLTAN_ID_TYPE [export_num1 * ZOLTAN_LID_SIZE];
      export_pids1 = new int [export_num1];
    }
    
    //  Get the GIDs of the export objects. 
    err = Zoltan_Comm_Do(import_plan1, 
                         msgtag+1, 
                         (char *) import_gids1,
                         gid_size,
                         (char *) export_gids1);

    //  Get the LIDs of the export objects.
    err = Zoltan_Comm_Do(import_plan1, 
                         msgtag+2, 
                         (char *) import_lids1,
                         lid_size,
                         (char *) export_lids1);

    //  Get the export processors.
    Zoltan_Comm_Info(import_plan1, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, export_pids1, NULL);

    Zoltan_Comm_Destroy(&import_plan1);
    import_plan1 = NULL;
    import_num1  = 0;
    import_gids1 = NULL;
    import_lids1 = NULL;
    import_pids1 = NULL;
    break;
  case ContactZoltanComm::ZOLTAN_EXPORT:
    export_num1  = comm1->Num_Export();
    export_gids1 = comm1->Export_GIDS();
    export_lids1 = comm1->Export_LIDS();
    export_pids1 = comm1->Export_Procs();
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  
  switch (comm2->Direction()) {
  case ContactZoltanComm::ZOLTAN_IMPORT:
    import_num2  = comm2->Num_Import();
    import_gids2 = comm2->Import_GIDS();
    import_lids2 = comm2->Import_LIDS();
    import_pids2 = comm2->Import_Procs();
    //  Create zoltan comm plan based on import lists.
    err = Zoltan_Comm_Create(&import_plan2, 
                             import_num2, 
                             import_pids2, 
                             zz->Communicator, 
                             msgtag, 
                             &export_num2);
    //  Allocate space for export buffers.
    if (export_num2>0) {
      export_gids2 = new ZOLTAN_ID_TYPE [export_num2 * ZOLTAN_GID_SIZE];
      export_lids2 = new ZOLTAN_ID_TYPE [export_num2 * ZOLTAN_LID_SIZE];
      export_pids2 = new int [export_num2];
    }
    
    //  Get the GIDs of the export objects. 
    err = Zoltan_Comm_Do(import_plan2, 
                         msgtag+1, 
                         (char *) import_gids2,
                         gid_size,
                         (char *) export_gids2);

    //  Get the LIDs of the export objects.
    err = Zoltan_Comm_Do(import_plan2, 
                         msgtag+2, 
                         (char *) import_lids2,
                         lid_size,
                         (char *) export_lids2);

    //  Get the export processors.
    Zoltan_Comm_Info(import_plan2, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, export_pids2, NULL);

    Zoltan_Comm_Destroy(&import_plan2);
    import_plan2 = NULL;
    import_num2  = 0;
    import_gids2 = NULL;
    import_lids2 = NULL;
    import_pids2 = NULL;
    break;
  case ContactZoltanComm::ZOLTAN_EXPORT:
    export_num2  = comm2->Num_Export();
    export_gids2 = comm2->Export_GIDS();
    export_lids2 = comm2->Export_LIDS();
    export_pids2 = comm2->Export_Procs();
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  
  num_export = export_num1+export_num2;
  if (num_export>0) {
    export_gids = new ZOLTAN_ID_TYPE [num_export * ZOLTAN_GID_SIZE];
    export_lids = new ZOLTAN_ID_TYPE [num_export * ZOLTAN_LID_SIZE];
    export_pids = new int [num_export];
    std::memcpy(export_gids,export_gids1,export_num1*gid_size);
    std::memcpy(&export_gids[export_num1*ZOLTAN_GID_SIZE],export_gids2,export_num2*gid_size);
    std::memcpy(export_lids,export_lids1,export_num1*lid_size);
    std::memcpy(&export_lids[export_num1*ZOLTAN_LID_SIZE],export_lids2,export_num2*lid_size);
    std::memcpy(export_pids,export_pids1,export_num1*sizeof(int));
    std::memcpy(&export_pids[export_num1],export_pids2,export_num2*sizeof(int));
  }
  
  if(export_alloc) {
    if(export_gids1) delete [] export_gids1;
    if(export_lids1) delete [] export_lids1;
    if(export_pids1) delete [] export_pids1;
  }

  //os<<"\nnum_exports = "<<num_export<<"\n";
  //for (int n=0; n<num_export; ++n) {
  //  ZOLTAN_ID_PTR gid = &export_gids[n*ZOLTAN_GID_SIZE];
  //  int entity_type = ContactZoltanGID::Type(gid);
  //  ContactHostGlobalID GID(gid);
  //  os<<"  "<<n<<":  type = "<<entity_type<<",  gid = "<<GID<<",  pid = "<<export_pids[n]<<"\n";
  //}
  //os.flush();
  
  err = Zoltan_Comm_Create(&export_plan, 
                           num_export, 
                           export_pids, 
                           zz->Communicator, 
                           msgtag, 
                           &num_import);
                           
  if (num_import>0) {
    import_gids = new ZOLTAN_ID_TYPE [num_import * ZOLTAN_GID_SIZE];
    import_lids = new ZOLTAN_ID_TYPE [num_import * ZOLTAN_LID_SIZE];
    import_pids = new int [num_import];
  }
                         
  //  Get the GIDs of the import objects. 
  err = Zoltan_Comm_Do(export_plan, 
                       msgtag+1, 
                       (char *) export_gids,
                       gid_size,
                       (char *) import_gids);

  //  Get the LIDs of the import objects.
  err = Zoltan_Comm_Do(export_plan, 
                       msgtag+2, 
                       (char *) export_lids,
                       lid_size,
                       (char *) import_lids);

  //  Get the import processors.
  Zoltan_Comm_Info(export_plan, NULL, NULL, NULL, NULL, NULL, NULL,
                   NULL, NULL, NULL, NULL, NULL, import_pids, NULL);
  
  //os<<"\nnum_imports = "<<num_import<<"\n";
  //for (int n=0; n<num_import; ++n) {
  //  ZOLTAN_ID_PTR gid = &import_gids[n*ZOLTAN_GID_SIZE];
  //  int entity_type = ContactZoltanGID::Type(gid);
  //  ContactHostGlobalID GID(gid);
  //  os<<"  "<<n<<":  type = "<<entity_type<<",  gid = "<<GID<<",  pid = "<<import_pids[n]<<"\n";
  //}
  //os.flush();
  
  //  Allocate space for export buffers.
  if (num_export>0) {
    send_sizes   = new int [num_export];
    send_indices = new int [num_export];
  }
  
  //  Allocate space for import buffers.
  if (num_import>0) {
    recv_sizes   = new int [num_import];
    recv_indices = new int [num_import];
  }
  
  zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data, 
                         ZOLTAN_GID_SIZE, 
                         ZOLTAN_LID_SIZE, 
                         num_export,
                         export_gids,
                         export_lids,
                         send_sizes,
                         &err);
                                    
  int tagsize = gid_size + sizeof(int);
  
  int send_size = 0;
  for (int i = 0; i < num_export; ++i) {
    send_size += send_sizes[i]+tagsize;
  }
  if (send_size>0) send_buffer = new char [send_size];

  // Modify sizes[] to contain message sizes, not object sizes
  for (int i=0; i<num_export; ++i) {
    send_sizes[i] += tagsize;
  }
  
  //  Resize the communication plan.
  int recv_size=0;
  err = Zoltan_Comm_Resize(export_plan, send_sizes, msgtag+3, &recv_size);
  if (recv_size>0) recv_buffer = new char [recv_size];
  
  // Modify sizes[] to contain object sizes, not message sizes
  for (int i=0; i<num_export; ++i) { 
    send_sizes[i] -= tagsize;
  }
  
}

ContactZoltanCommUtils::~ContactZoltanCommUtils()
{
  switch (type) {
  case 0:
    switch (direction) {
    case IMPORT:
      if (export_gids  != NULL) delete [] export_gids;
      if (export_lids  != NULL) delete [] export_lids;
      if (export_pids  != NULL) delete [] export_pids;
      if (import_gids  != NULL) delete [] import_gids;
      break;
    case EXPORT:
      if (import_gids  != NULL) delete [] import_gids;
      if (import_lids  != NULL) delete [] import_lids;
      if (import_pids  != NULL) delete [] import_pids;
      break;
    case BOTH:
      if (import_gids  != NULL) delete [] import_gids;
      if (import_lids  != NULL) delete [] import_lids;
      if (import_pids  != NULL) delete [] import_pids;
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    break;
  case 1:
    if (import_gids  != NULL) delete [] import_gids;
    if (import_lids  != NULL) delete [] import_lids;
    if (import_pids  != NULL) delete [] import_pids;
    if (export_gids  != NULL) delete [] export_gids;
    if (export_lids  != NULL) delete [] export_lids;
    if (export_pids  != NULL) delete [] export_pids;
    break;
  }
  if (send_buffer  != NULL) delete [] send_buffer;
  if (recv_buffer  != NULL) delete [] recv_buffer;
  if (send_sizes   != NULL) delete [] send_sizes;
  if (recv_sizes   != NULL) delete [] recv_sizes;
  if (send_indices != NULL) delete [] send_indices;
  if (recv_indices != NULL) delete [] recv_indices;
  if (import_plan  != NULL) Zoltan_Comm_Destroy(&import_plan);
  if (export_plan  != NULL) Zoltan_Comm_Destroy(&export_plan);
  import_gids  = NULL;
  import_lids  = NULL;
  import_pids  = NULL;
  export_gids  = NULL;
  export_lids  = NULL;
  export_pids  = NULL;
  send_buffer  = NULL;
  recv_buffer  = NULL;
  send_sizes   = NULL;
  recv_sizes   = NULL;
  send_indices = NULL;
  recv_indices = NULL;
  import_plan  = NULL;
  export_plan  = NULL;
}

void ContactZoltanCommUtils::Migrate(ContactParOStream&  os)
{
  int err      = 0;
  int msgtag   = 15;
  int gid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_GID_SIZE);
  int tagsize  = gid_size + sizeof(int);
  
  if (num_export) {
    POSTCONDITION(send_buffer);
  }
  if (num_import) {
    POSTCONDITION(recv_buffer);
  }
  
  int   send_cnt = 0;
  char* tmp = send_buffer;
  for (int i=0; i<num_export; ++i) {
    // Pack the object's global ID
    memcpy(tmp, &(export_gids[i*ZOLTAN_GID_SIZE]), gid_size);
    tmp += gid_size;
  
    /* Pack the object's size */
    *((int *)tmp) = send_sizes[i];
    tmp += sizeof(int);
    
    tmp += send_sizes[i];
    
    send_cnt += tagsize;
    
    send_indices[i] = send_cnt;
    
    send_cnt += send_sizes[i];
  }

  zz->Pack_Obj_Multi(zz->Pack_Obj_Multi_Data,
                     ZOLTAN_GID_SIZE, 
                     ZOLTAN_LID_SIZE,
                     num_export,
                     export_gids,
                     export_lids,
                     export_pids,
                     send_sizes,
                     send_indices,
                     send_buffer, 
                     &err);

  err = Zoltan_Comm_Do(export_plan, msgtag+4, send_buffer, 1, recv_buffer);
  
  int recv_cnt = 0;
  tmp = recv_buffer;
  for (int i=0; i<num_import; ++i) {
    memcpy(&import_gids[i*ZOLTAN_GID_SIZE], tmp, gid_size);
    tmp += gid_size;
    
    recv_sizes[i] = *((int *)tmp);
    tmp += sizeof(int);
    
    recv_cnt += tagsize;
    
    recv_indices[i] = recv_cnt;
    
    tmp += recv_sizes[i];
    recv_cnt += recv_sizes[i];
  }
  zz->Unpack_Obj_Multi(zz->Unpack_Obj_Multi_Data,
                       ZOLTAN_GID_SIZE, 
                       num_import,
                       import_gids,
                       recv_sizes,
                       recv_indices,
                       recv_buffer, 
                       &err);
}

void ContactZoltanCommUtils::Resize(ContactParOStream&  os)
{
  int err      = 0;
  int msgtag   = 15;
  int gid_size = sizeof(ZOLTAN_ID_TYPE)*(ZOLTAN_GID_SIZE);
  int tagsize  = gid_size + sizeof(int);
  
  zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data, 
                         ZOLTAN_GID_SIZE, 
                         ZOLTAN_LID_SIZE, 
                         num_export,
                         export_gids,
                         export_lids,
                         send_sizes,
                         &err); 
  
  int send_size = 0;
  for (int i = 0; i < num_export; ++i) {
    send_size += send_sizes[i]+tagsize;
  }
  if (send_buffer!=NULL) delete [] send_buffer;
  if (send_size>0) 
    send_buffer = new char [send_size];
  else
    send_buffer = NULL;
    
  // Modify sizes[] to contain message sizes, not object sizes
  for (int i=0; i<num_export; ++i) {
    send_sizes[i] += tagsize;
  }
  
  //  Resize the communication plan.
  int recv_size=0;
  err = Zoltan_Comm_Resize(export_plan, send_sizes, msgtag+3, &recv_size);
  if (recv_buffer!=NULL) delete [] recv_buffer;
  if (recv_size>0) 
    recv_buffer = new char [recv_size];
  else
    recv_buffer = NULL;
  
  // Modify sizes[] to contain object sizes, not message sizes
  for (int i=0; i<num_export; ++i) {
    send_sizes[i] -= tagsize;
  }
}

#endif
