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


#ifndef ContactZoltan_h_
#define ContactZoltan_h_

#ifndef CONTACT_NO_MPI

#include <zoltan.h>
#include "Contact_Defines.h"
#include "ContactTopology.h"

class ContactBoundingBox;

class ContactZoltan
{

public:

 ContactZoltan( MPI_Comm communicator, int& ierr );

 ~ContactZoltan();

 //Load Balance Calls

 void Set_GeomCallBacks( ContactSearch* );

 void Set_ToplevelExportCallBacks( ContactSearch* );

 void Set_InteractionCallBacks( ContactSearch* );

 void Set_HostidQueryCallBacks( ContactSearch* );
 
 void Set_DynamicLoadBalanceCallBacks( ContactSearch* );
 
 void Set_GhostingExportCallBacks( ContactSearch* );
 
 void Set_GhostingImportCallBacks( ContactSearch* );
 
 void Set_UpdateGhostingExportCallBacks( ContactSearch* );
 
 void Set_UpdateGhostingImportCallBacks( ContactSearch* );
 
 void Set_UpdateTiedImportCallBacks( ContactSearch* );

 int Set_Method( char* string );

 int Set_Param( char* param, char * val );

 int Balance();

 void Evaluate( int    print_stats,
		int*   num_objects,
		float* object_weights,
                int*   num_cuts,
    		float* cut_weights,
		int*   num_boundary_objects,
		int*   num_adj_procs,
		int*   ierr );

 int Free_Data();
  
 //Decomposition Augmentation
 
 int Point_Assign ( Real* coords,
		    int* proc );

 int Box_Assign ( Real xmin,
		  Real ymin,
		  Real zmin,
		  Real xmax,
		  Real ymax,
		  Real zmax,
		  int* procs,
		  int* numprocs );
                  
 int Box_Assign (const ContactBoundingBox &box,
                 int* procs,
                 int* numprocs);

 //Migration Functions
 
  int Compute_Destinations ( int        num_import,
			     ZOLTAN_ID_PTR  import_global_ids,
			     ZOLTAN_ID_PTR  import_local_ids,
			     int*       import_procs,
			     int*       num_export,
			     ZOLTAN_ID_PTR* export_global_ids,
			     ZOLTAN_ID_PTR* export_local_ids,
			     int**      export_procs );

  int Help_Migrate ( int       num_import,
		     ZOLTAN_ID_PTR import_global_ids,
		     ZOLTAN_ID_PTR import_local_ids,
		     int*      import_procs,
		     int       num_export,
		     ZOLTAN_ID_PTR export_global_ids,
		     ZOLTAN_ID_PTR export_local_ids,
		     int*      export_procs );
  
  inline float     Version() { return version; };
  inline int       New_Decomp() { return new_decomp; };
  inline int       Num_Import() { return num_import; };
  inline ZOLTAN_ID_PTR Import_Global_IDs() { return import_global_ids; };
  inline ZOLTAN_ID_PTR Import_Local_IDs() { return import_local_ids; };
  inline int*      Import_Procs() { return import_procs; };
  inline int       Num_Export() { return num_export; };
  inline ZOLTAN_ID_PTR Export_Global_IDs() { return export_global_ids; };
  inline ZOLTAN_ID_PTR Export_Local_IDs() { return export_local_ids; };
  inline int*      Export_Procs() { return export_procs; };
  inline void      Position(VariableHandle p) { position=p; };
  inline VariableHandle Position() { return position; };
  inline Zoltan_Struct* Get_ZoltanPtr() { return Zoltan_Ptr_; };

private:

 VariableHandle position;
 
 Zoltan_Struct* Zoltan_Ptr_; 
 float          version;
 int            new_decomp;
 int            num_gid_entries;
 int            num_lid_entries;
 int            num_import;
 ZOLTAN_ID_PTR  import_global_ids;
 ZOLTAN_ID_PTR  import_local_ids;
 int*           import_procs;
 int            num_export;
 ZOLTAN_ID_PTR  export_global_ids;
 ZOLTAN_ID_PTR  export_local_ids;
 int*           export_procs;

};

#endif
#endif
