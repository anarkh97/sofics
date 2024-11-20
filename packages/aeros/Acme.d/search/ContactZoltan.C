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

#include <iostream>
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactZoltan.h"
#include "Zoltan_Interface.h"
#include "ContactBoundingBox.h"

ContactZoltan::ContactZoltan(MPI_Comm communicator, int& error ) :
  version(0.0), 
  new_decomp(0),
  num_gid_entries(0),
  num_lid_entries(0),
  num_import(0),
  import_global_ids(NULL),
  import_local_ids(NULL),
  import_procs(NULL),
  num_export(0),
  export_global_ids(NULL),
  export_local_ids(NULL),
  export_procs(NULL)
{
  Zoltan_Initialize( 0, NULL, &version );
  Zoltan_Ptr_ = Zoltan_Create( communicator );
  POSTCONDITION( Zoltan_Ptr_ );
  if( Zoltan_Ptr_ )
    error = ZOLTAN_OK;
  else
    error = ZOLTAN_FATAL;
}

ContactZoltan::~ContactZoltan()
{
  Zoltan_Destroy( &Zoltan_Ptr_ );
}

//Load Balance Calls

void ContactZoltan::Set_GeomCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Num_Obj_Fn( Zoltan_Ptr_, 
			 ContactQueryNumObjects, 
			 reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Obj_List_Fn( Zoltan_Ptr_, 
			  ContactQueryObjectList, 
			  reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Num_Geom_Fn( Zoltan_Ptr_, 
			  ContactQueryNumGeomObjects, 
			  reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Geom_Multi_Fn( Zoltan_Ptr_, 
                            ContactQueryGeomMultiValues, 
                            reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_ToplevelExportCallBacks( ContactSearch* search )
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
                                ContactMigrateEntityExportSizes, 
                                reinterpret_cast<void*>(search) );

  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
				ContactMigrateEntityExportPack,
				reinterpret_cast<void*>(search) );

  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
				  ContactMigrateEntityUnpack,
				  reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_InteractionCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
				ContactMigrateInteractionSizes, 
				reinterpret_cast<void*>(search) );

  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
				ContactMigratePackInteractions, 
				reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
				  ContactMigrateUnpackInteractions, 
				  reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_HostidQueryCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Fn( Zoltan_Ptr_, 
			  ContactHostidQuerySize, 
			  reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Pack_Obj_Fn( Zoltan_Ptr_, 
			  ContactHostidQueryPack, 
			  reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Fn( Zoltan_Ptr_, 
			    ContactHostidQueryUnpack, 
			    reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_DynamicLoadBalanceCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
			        ContactDynamicLoadBalanceSize, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
			        ContactDynamicLoadBalancePack, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
			          ContactDynamicLoadBalanceUnpack, 
			          reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_GhostingExportCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
			        ContactGhostingExportSizes, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
			        ContactGhostingExportPack, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
			          ContactGhostingUnpack, 
			          reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_GhostingImportCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
			        ContactGhostingImportSizes, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
			        ContactGhostingImportPack, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
			          ContactGhostingUnpack, 
			          reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_UpdateGhostingExportCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
			        ContactUpdateGhostingExportSizes, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
			        ContactUpdateGhostingExportPack, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
			          ContactUpdateGhostingUnpack, 
			          reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_UpdateGhostingImportCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
			        ContactUpdateGhostingImportSizes, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
			        ContactUpdateGhostingImportPack, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
			          ContactUpdateGhostingUnpack, 
			          reinterpret_cast<void*>(search) );
}

void ContactZoltan::Set_UpdateTiedImportCallBacks ( ContactSearch* search)
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
			        ContactUpdateTiedImportSizes, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
			        ContactUpdateTiedImportPack, 
			        reinterpret_cast<void*>(search) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
			          ContactUpdateTiedUnpack, 
			          reinterpret_cast<void*>(search) );
}

int ContactZoltan::Set_Method( char * string )
{
  char* method = (char*) "LB_METHOD";
  return Zoltan_Set_Param( Zoltan_Ptr_, method, string );
}

int ContactZoltan::Set_Param( char * param, char * val )
{
  return Zoltan_Set_Param( Zoltan_Ptr_, param, val );
}

int ContactZoltan::Balance()
{
  //Zoltan_Generate_Files(Zoltan_Ptr_, "heaphy", 1);
  return Zoltan_LB_Balance( Zoltan_Ptr_, &new_decomp, 
			    &num_gid_entries, &num_lid_entries,
			    &num_import, &import_global_ids, 
			    &import_local_ids, &import_procs, 
			    &num_export, &export_global_ids,
			    &export_local_ids, &export_procs );
}

/* doesn't build with zoltan in trilinos version 10.8.3
void ContactZoltan::Evaluate ( int    print_stats,
			       int*   num_objects,
			       float* object_weights,
                               int*   num_cuts,
			       float* cut_weights,
			       int*   num_boundary_objects,
			       int*   num_adj_procs,
			       int*   ierr )
{
  *ierr = Zoltan_LB_Eval( Zoltan_Ptr_, print_stats, num_objects, 
			  object_weights, num_cuts, cut_weights, 
			  num_boundary_objects, num_adj_procs );
}
*/

int ContactZoltan::Free_Data ()
{
  return Zoltan_LB_Free_Data( &import_global_ids, &import_local_ids, 
			      &import_procs, &export_global_ids, 
			      &export_local_ids, &export_procs );
}

//Decomposition Augmentation

int ContactZoltan::Point_Assign ( Real* coords, int* proc )
{
  PRECONDITION( sizeof(Real) == sizeof(double) );
  return Zoltan_LB_Point_Assign( Zoltan_Ptr_, coords, proc );
}

int ContactZoltan::Box_Assign ( Real xmin,
				Real ymin,
				Real zmin,
				Real xmax,
				Real ymax,
				Real zmax,
				int* procs,
				int* numprocs )
{
  PRECONDITION( sizeof(Real) == sizeof(double) );
  return Zoltan_LB_Box_Assign( Zoltan_Ptr_, xmin, ymin, zmin, xmax, 
			       ymax, zmax, procs, numprocs );
}

int ContactZoltan::Box_Assign (const ContactBoundingBox &box,
                               int *procs,
                               int *numprocs) {
  PRECONDITION( sizeof(Real) == sizeof(double) );
  return Zoltan_LB_Box_Assign( Zoltan_Ptr_, 
                               box.get_x_min(), 
                               box.get_y_min(), 
                               box.get_z_min(), 
                               box.get_x_max(),
                               box.get_y_max(), 
                               box.get_z_max(), 
                               procs, numprocs );  
}


int ContactZoltan::Compute_Destinations ( int        Num_import,
					  ZOLTAN_ID_PTR  Import_global_ids,
					  ZOLTAN_ID_PTR  Import_local_ids,
					  int*       Import_procs,
					  int*       Num_export,
					  ZOLTAN_ID_PTR* Export_global_ids,
					  ZOLTAN_ID_PTR* Export_local_ids,
					  int**      Export_procs )
{
  return Zoltan_Compute_Destinations( Zoltan_Ptr_, 
				      Num_import, Import_global_ids, 
				      Import_local_ids, Import_procs, 
				      Num_export, Export_global_ids, 
				      Export_local_ids, Export_procs );
}

int ContactZoltan::Help_Migrate ( int       Num_import,
				  ZOLTAN_ID_PTR Import_global_ids,
				  ZOLTAN_ID_PTR Import_local_ids,
				  int*      Import_procs,
				  int       Num_export,
				  ZOLTAN_ID_PTR Export_global_ids,
				  ZOLTAN_ID_PTR Export_local_ids,
				  int*      Export_procs )
{
  return Zoltan_Help_Migrate( Zoltan_Ptr_, 
			      Num_import, Import_global_ids, 
			      Import_local_ids, Import_procs,
			      Num_export, Export_global_ids, 
			      Export_local_ids, Export_procs );
}


#endif
