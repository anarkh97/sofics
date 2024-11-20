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


#ifndef ContactFaceFaceInteraction_h_
#define ContactFaceFaceInteraction_h_

#include "Contact_Defines.h"
#include "ContactInteractionEntity.h"
#include "ContactFace.h"
#include "ContactSearch.h"
 
template<typename DataType>
struct ContactFaceFaceVertex {
  DataType slave_x;
  DataType slave_y;
  DataType master_x;
  DataType master_y;
  int  slave_edge_id;
  int  master_edge_flag;
  DataType *slave_x_derivatives;
  DataType *slave_y_derivatives;
  DataType *master_x_derivatives;
  DataType *master_y_derivatives;
  DataType *slave_x_second_derivatives;
  DataType *slave_y_second_derivatives;
  DataType *master_x_second_derivatives;
  DataType *master_y_second_derivatives;
};

class CString;
class ContactTopologyEntityList;
class ContactTopologyEntityHash;
class ContactHostGlobalID;
class ContactFixedSizeAllocator;

template <typename DataType>
class ContactFaceFaceInteraction : public ContactInteractionEntity<DataType> {
  
 public:
  
  enum InteractionSource { UNKNOWN_SOURCE=-1,CLOSEST_POINT_PROJECTION_1=1, 
                           CLOSEST_POINT_PROJECTION_2, MOVING_INTERSECTION };

  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#include "contact_variables.def"
#include "contact_variables.undefine"
                     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactFaceFaceInteraction();
  ContactFaceFaceInteraction( ContactFace<DataType>*, ContactFace<DataType>*, 
			      int, int*, int*, DataType*, DataType*, DataType (*)[42] = NULL,
                              DataType (*)[42] = NULL, DataType (*)[42] = NULL, DataType (*)[42] = NULL );
  ContactFaceFaceInteraction( ContactFaceFaceInteraction& );
  static ContactFaceFaceInteraction* new_ContactFaceFaceInteraction(
            ContactFixedSizeAllocator&, ContactFace<DataType>*, ContactFace<DataType>*, 
	    int, int*, int*, DataType*, DataType*, DataType (*)[42] = NULL, DataType (*)[42] = NULL,
            DataType (*)[42] = NULL, DataType (*)[42] = NULL );
  static ContactFaceFaceInteraction* new_ContactFaceFaceInteraction(
	     ContactFixedSizeAllocator& );
  static ContactFaceFaceInteraction* new_ContactFaceFaceInteraction(
	     ContactFixedSizeAllocator&, ContactFaceFaceInteraction& );
  ~ContactFaceFaceInteraction();
  
#ifndef CONTACT_NO_MPI
  inline ContactZoltanLID& Zoltan_LID() { return zoltan_lid; };
  inline ContactZoltanGID& Zoltan_GID() { return zoltan_gid; };
  inline void ZoltanFaceLID(LB_ID_PTR lid, int flag=0) 
    {//if (flag) zoltan_lid.ZoltanLID(CT_FACE, 
     //                               master_face_entity_data.index_in_owner_proc_array, 
     //                               lid);
     //else      zoltan_lid.ZoltanLID(CT_FACE, 
     //                               master_face_entity_data.index_in_proc_array, 
     //                               lid);
     zoltan_lid.ZoltanLID(CT_FACE, 
                                    master_face_entity_data.index_in_owner_proc_array, 
                                    lid);};
  inline void ZoltanFaceGID(LB_ID_PTR gid) 
    {zoltan_gid.ZoltanGID(CT_FACE, 
                          master_face_entity_data.host_gid[0],
                          master_face_entity_data.host_gid[1],
                          gid);};
#endif

  inline ContactFace<DataType>* SlaveFace() {return slave_face;};
  inline typename ContactInteractionEntity<DataType>::entity_data* SlaveFaceEntityData()  {return &slave_face_entity_data;};
  inline ContactFace<DataType>* MasterFace() {return master_face;};
  inline typename ContactInteractionEntity<DataType>::entity_data* MasterFaceEntityData() {return &master_face_entity_data;};
  inline void Set_SlaveFace(ContactFace<DataType>* cf) {slave_face=cf;};
  inline void Set_MasterFace(ContactFace<DataType>* cf) {master_face=cf;};
  int Set_SlaveFaceEntityData();
  int Set_MasterFaceEntityData();
  inline int NumEdges() {return num_edges;};
  inline void NumEdges(int n) {num_edges=n;vertices=new ContactFaceFaceVertex<DataType>[n+1];};
  inline int NumDerivatives() {return num_derivatives;};
  inline int NumSecondDerivatives() {return num_second_derivatives;};
  inline ContactFaceFaceVertex<DataType>* Get_Vertices() { return vertices; };
  inline ContactFaceFaceVertex<DataType>* Get_Vertex(int n) { return &vertices[n]; };
#if (MAX_FFI_DERIVATIVES > 0)
  void Set_Derivatives(ContactFaceFaceInteraction<ActiveScalar> *active_cffi);
#endif

  inline DataType& Scalar_Var( VariableHandle vh ) {return DataArray[vh];};
  inline DataType* Vector_Var( VariableHandle vh ) 
    { return (DataArray+NUMBER_SCALAR_VARS+3*vh); };

  inline int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline void Initialize_Memory() {std::memset(this->DataArray_Buffer(), 0, DataArray_Length()*sizeof(DataType));};
  void Connect_SlaveFace ( ContactTopologyEntityList& );
  void Connect_MasterFace( ContactTopologyEntityList& );
  void Connect_SlaveFace ( ContactTopologyEntityHash& );
  void Connect_MasterFace( ContactTopologyEntityHash& );
  void Connect_SlaveFace ( ContactTopology* );
  void Connect_MasterFace( ContactTopology* );
  void Connect_SlaveFace ( ContactFace<DataType>* );
  void Connect_MasterFace( ContactFace<DataType>* );

  // Parallel packing/unpacking functions
  int  Size();
  void Pack( char* buffer );
  void Unpack( char* buffer );
  void Copy( ContactFaceFaceInteraction* src );

  // Restart Pack/Unpack functions
  int  Restart_Size();
  void Restart_Pack( DataType* buffer );
  void Restart_Unpack( DataType* buffer );

  int Data_Size();
  
 protected:

 private:
  
  DataType DataArray[NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS+1];
  ContactFace<DataType>* slave_face;
  typename ContactInteractionEntity<DataType>::entity_data slave_face_entity_data;
  ContactFace<DataType>* master_face;
  typename ContactInteractionEntity<DataType>::entity_data master_face_entity_data;
  int num_edges;
  int num_derivatives;
  int num_second_derivatives;
  ContactFaceFaceVertex<DataType>* vertices;
#ifndef CONTACT_NO_MPI
  ContactZoltanLID zoltan_lid;
  ContactZoltanGID zoltan_gid;
#endif
};

#endif // ContactFaceFaceInteraction_h_
