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


#ifndef _ContactZoltanComm_h_
#define _ContactZoltanComm_h_

#ifndef CONTACT_NO_MPI

#include "zoltan.h"
#include "Contact_Defines.h"
#include "ContactZoltanID.h"
#include <vector>
#include <map>

class ContactParOStream;

class ContactZoltanComm {

 public:
 
  enum Zoltan_Dir { ZOLTAN_EXPORT, ZOLTAN_IMPORT, ZOLTAN_UNKNOWN };

  ContactZoltanComm();
  ContactZoltanComm( Zoltan_Dir);
  ~ContactZoltanComm();

  void CleanUp();
  void Initialize(Zoltan_Dir);

  void Add_Import( LB_ID_PTR lid, LB_ID_PTR gid, int proc, int check=1 );
  void Add_Export( LB_ID_PTR lid, LB_ID_PTR gid, int proc, int check=1 );
  
  inline int       Num_Import() {return num_import;};
  inline int       Num_Export() {return num_export;};
  inline LB_ID_PTR Import_GIDS() {return import_gids;};
  inline LB_ID_PTR Export_GIDS() {return export_gids;};
  inline LB_ID_PTR Import_LIDS() {return import_lids;};
  inline LB_ID_PTR Export_LIDS() {return export_lids;};
  inline int*      Import_Procs() {return import_procs;};
  inline int*      Export_Procs() {return export_procs;};
  
  void Set_Export(int num, LB_ID_PTR gids, LB_ID_PTR lids, int* procs);

  
  inline Zoltan_Dir Direction() {return direction;};
  
  void Print(ContactParOStream& os);

 private:
  Zoltan_Dir direction;

  int        num_import;
  int        cur_import_capacity;
  LB_ID_PTR  import_gids;
  LB_ID_PTR  import_lids;
  int*       import_procs;

  int        num_export;
  int        cur_export_capacity;
  LB_ID_PTR  export_gids;
  LB_ID_PTR  export_lids;
  int*       export_procs;


  std::map<std::pair<ContactZoltanGID, int>, int> lookup;
  void add_hash_entry( LB_ID_PTR, int );
  int  find_hash_entry( LB_ID_PTR, int );
};


#endif
#endif
