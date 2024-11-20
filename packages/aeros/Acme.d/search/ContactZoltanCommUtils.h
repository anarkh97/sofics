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


#ifndef ContactZoltanCommUtils_h_
#define ContactZoltanCommUtils_h_

#ifndef CONTACT_NO_MPI

#include <mpi.h>
#include <zoltan.h>
#include <zz_const.h>
#include "Contact_Defines.h"
#include "contact_assert.h"
#include "ContactZoltanComm.h"
#include "ContactParOStream.h"


class ContactZoltanCommUtils {

 public:
 
    enum CommDirection { IMPORT, EXPORT, BOTH };
 
    ContactZoltanCommUtils(CommDirection dir,
                           Zoltan_Struct* z_ptr,
                           int num, LB_ID_TYPE* gids, 
                           LB_ID_TYPE* lids, int* pids,
                           ContactParOStream& os);
                           
    ContactZoltanCommUtils(ContactZoltanComm* comm1,
                           ContactZoltanComm* comm2,
                           Zoltan_Struct*     z_ptr,
                           ContactParOStream& os);

    ~ContactZoltanCommUtils();
    
    void Resize(ContactParOStream& os);
    
    void Migrate(ContactParOStream& os);
    
    inline int       Num_Import()   {return num_import;};
    inline int       Num_Export()   {return num_export;};
    inline LB_ID_PTR Import_GIDS()  {return import_gids;};
    inline LB_ID_PTR Export_GIDS()  {return export_gids;};
    inline LB_ID_PTR Import_LIDS()  {return import_lids;};
    inline LB_ID_PTR Export_LIDS()  {return export_lids;};
    inline int*      Import_Procs() {return import_pids;};
    inline int*      Export_Procs() {return export_pids;};
  
  private:
    int              type;
    
    CommDirection    direction;
  
    Zoltan_Struct*   zz;
    
    ZOLTAN_COMM_OBJ* import_plan;
    ZOLTAN_COMM_OBJ* export_plan;
    
    int              num_import;
    LB_ID_TYPE*      import_gids;
    LB_ID_TYPE*      import_lids;
    int*             import_pids;
    
    int              num_export;
    LB_ID_TYPE*      export_gids;
    LB_ID_TYPE*      export_lids;
    int*             export_pids;
    
    char*            send_buffer;
    char*            recv_buffer;
    
    int*             send_sizes;
    int*             recv_sizes;
    int*             send_indices;
    int*             recv_indices;

};

#endif

#endif
