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


#ifndef ContactTopologyEntity_H_
#define ContactTopologyEntity_H_

#include "ContactEntity.h"
#include "Contact_Defines.h"
#include "ContactHostGlobalID.h"
#include "ContactZoltanID.h"
#include "contact_assert.h"
#include "ContactBoundingBox.h"
#ifndef CONTACT_NO_MPI
#include "zoltan.h"
#endif

template<typename DataType>
class ContactTopologyEntity {

 public:
 
  typedef struct {
    int owner;
    int owner_proc_array_index;
    int host_gid[2];
  } connection_data;
  
  enum ProcessorOwnership { OWNED, NOT_OWNED };
  
  enum SearchContext      { ACTIVE=1, MASTER=2, SLAVE=4, 
                            TIED=8, GLUED=16, SLIDING=32,
                            IN_PROXIMITY=64, 
                            ACTIVE_IN_GLOBAL_SEARCH=128, 
                            ACTIVE_IN_TRACK_SEARCH=256,
                            ACTIVE_FOR_GLOBAL_GHOSTING=512,
                            ACTIVE_FOR_TRACK_GHOSTING=1024,
			    GEOMETRY_UPDATE=2048,
			    GLOBAL_SEARCH_SLAVE=4096,
			    TRACK_SEARCH_SLAVE=8192,
			    ACTIVE_FOR_GHOSTING_RCB=16384,
                            GHOSTED_FOR_SEARCH=32768 };
                            
  //enum ActiveStatus       { INACTIVE = 0,
  //                          GLOBAL_SEARCH,
  //                          TRACKED_SEARCH,
  //                          ALL_SEARCH };
  enum Packed_Variables   { BASE_TYPE = 0, DERIVED_TYPE, 
			    GID_HI, GID_LO, BLOCK_ID, CONTEXT,
                            PRIM_OWNER, SEC_OWNER, 
			    OWNER_PROC_INDEX, HOST_INDEX,
                            PROC_INDEX, PRIM_PROC_INDEX,
			    OWNERSHIP, 
                            NUMBER_PACKED_VARS };
  enum Hierarchy_Status   { UNDETERMINED = -1, TOPLEVEL = 0, DERIVED = 1 };
  
  enum Migration_Direction   { PRIMARY_TO_PRIMARY, 
                               PRIMARY_TO_SECONDARY, 
                               SECONDARY_TO_PRIMARY };
  //int  search_status;
  int  temp_tag;
  int  temp_tag1;
  int  in_proximity;
  int  fcs_index;

  // Constructors/Destructors
  ContactTopologyEntity(DataType *data_array_, const ContactType base_type_);
  ContactTopologyEntity( int, int, DataType *data_array_, const ContactType base_type_);
  virtual ~ContactTopologyEntity();

  // Access Functions
  inline bool Shared() { return shared; };
  inline void Shared( bool s) { shared = s;};
  inline ProcessorOwnership Ownership() { return ownership; };
  inline void Ownership( ProcessorOwnership  PO) { ownership = PO;};
  inline int  Owner() { return owner; };
  inline void Owner(int proc) { owner = proc; };
  inline int  Secondary_Owner() { return secondary_owner; };
  inline void Secondary_Owner(int proc) { secondary_owner = proc; };
  inline ContactHostGlobalID& Global_ID() { return global_id; };
  inline void Global_ID(int hi, int lo) 
         { global_id.HiInt(hi); global_id.LoInt(lo); };
#ifndef CONTACT_NO_MPI
  inline void ZoltanLID(const int input_base_type, LB_ID_PTR lid, int flag=0) {
    if (flag) {
      lid[0] = input_base_type;
      lid[1] = primary_proc_array_index;
    } else {
      lid[0] = input_base_type;
      lid[1] = proc_array_index;
    }
  };
  inline void ZoltanGID(const int input_base_type, LB_ID_PTR gid) {
    gid[0] = input_base_type;
    gid[1] = global_id.HiInt();
    gid[2] = global_id.LoInt();
  };
#endif
  inline int  BlockID() { return block_id; };
  inline void BlockID(int i) { block_id=i; };
  inline int  PrimaryProcArrayIndex() { return primary_proc_array_index; };
  inline void PrimaryProcArrayIndex(int i) { primary_proc_array_index = i; };
  inline int  ProcArrayIndex() { return proc_array_index; };
  inline void ProcArrayIndex(int i) { proc_array_index = i; };
  inline int  EnfArrayIndex() const { return enf_array_index; };
  inline void EnfArrayIndex(int i) { enf_array_index = i; };
  inline void OwnerProcArrayIndex(int i) { owner_proc_array_index = i; };
  inline int  OwnerProcArrayIndex() { return owner_proc_array_index; };
  inline int  HostArrayIndex() { return host_array_index; };
  inline void HostArrayIndex(int i) { host_array_index = i; };
  inline int  HostGlobalArrayIndex() { return host_global_array_index; };
  inline void HostGlobalArrayIndex(int i) { host_global_array_index = i; };
  //inline ActiveStatus StatusFlag() { return status_flag; };
  //inline void StatusFlag(ActiveStatus status) { status_flag = status; };

  inline void Entity_Key(int key) { entity_key=key; };
  inline int  Entity_Key()        { return entity_key; };
  
  // Access Functions
  inline ContactType Base_Type() {return base_type;};
  inline DataType*  DataArray_Buffer() {return data_array;};
  inline DataType* Variable( const VariableHandle & vh) {return DataArray_Buffer()+vh;};
  
  void Display(ContactParOStream&);
  
  int BaseType(char *buffer);
  int DerivedType(char* buffer);
  
  void ClearContext()                    { context=0; };
  void ClearNonTiedContext()             { context&=(TIED|GHOSTED_FOR_SEARCH); };
  void ClearContextBit(unsigned int mask){ context&=(~mask); };
  void SetContext(unsigned int value)    { context=value; };
  void SetContextBit(unsigned int value) { context|=value; };
  bool CheckContext(unsigned int mask)   { return context&mask; };
  unsigned int GetContext() { return context; };
  
  bool IsActive()                  { return context&ACTIVE; };
  void SetActive()                 { context|=ACTIVE; };
  bool IsMaster()                  { return context&MASTER; };
  void SetMaster()                 { context|=MASTER; };
  bool IsSlave()                   { return context&SLAVE; };
  void SetSlave()                  { context|=SLAVE; };
  bool IsInProximity()             { return context&IN_PROXIMITY; };
  void SetInProximity()            { context|=IN_PROXIMITY; };
  bool IsActiveInGlobalSearch()    { return context&ACTIVE_IN_GLOBAL_SEARCH; };
  void SetActiveInGlobalSearch()   { context|=ACTIVE_IN_GLOBAL_SEARCH; };
  bool IsActiveInTrackSearch()     { return context&ACTIVE_IN_TRACK_SEARCH; };
  void SetActiveInTrackSearch()    { context|=ACTIVE_IN_TRACK_SEARCH; };
  bool IsActiveInGlobalGhosting()  { return context&ACTIVE_FOR_GLOBAL_GHOSTING; };
  void SetActiveInGlobalGhosting() { context|=ACTIVE_FOR_GLOBAL_GHOSTING; };
  bool IsActiveInTrackGhosting()   { return context&ACTIVE_FOR_TRACK_GHOSTING; };
  void SetActiveInTrackGhosting()  { context|=ACTIVE_FOR_TRACK_GHOSTING; };
  bool IsGhostRCB()                { return context&ACTIVE_FOR_GHOSTING_RCB; };
  void SetGhostRCB()               { context|=ACTIVE_FOR_GHOSTING_RCB; };



 protected:
  // Parallel Packing/Unpacking Functions
  // These are protected because the derived class function generally needs
  // to add data and should call this function.  It shouldn't be called
  // directly
  inline int  Size(const int data_array_length);
  inline void Pack( char*, const int data_array_length );
  inline void Unpack( char*, const int data_array_length );
  inline void Copy( ContactTopologyEntity*, const int data_array_length );
  
  inline int  Size_ForSecondary(const int data_array_length);
  inline void Pack_ForSecondary( char*, const int data_array_length );
  inline void Unpack_ForSecondary( char*, const int data_array_length );
  inline void Pack_ForSecondary( char*, DataType* DataArray, const int data_array_length );
  inline void Unpack_ForSecondary( char*, DataType* DataArray, const int data_array_length );
  inline void Copy_ForSecondary( ContactTopologyEntity*, const int data_array_length );
  
  inline int  Size_ForDataUpdate(const int data_array_length);
  inline void Pack_ForDataUpdate( char*, const int data_array_length );
  inline void Unpack_ForDataUpdate( char*, const int data_array_length );

  inline int  PackConnection(ContactTopologyEntity*, int*);
  inline void PackConnection(ContactTopologyEntity*, connection_data*);
  
  inline int  PackConnection(connection_data*, int*);
  inline int  UnPackConnection(connection_data*, int*);
  
  inline int Unpack( Packed_Variables v, char* buffer)
    {
      int* i_buffer = reinterpret_cast<int*>(buffer);
      return (i_buffer[v]);
    };

  int host_global_array_index;  // This is an object's 0..N-1 index into
                                // the host array on a processor
  int host_array_index;         // This is an object's 0..N-1 index into
                                // the host array (on a per block basis)
  int proc_array_index;         // This is an object's 0..N-1 index 
                                // on a processor
  int enf_array_index;          // This is an object's 0..N-1 index 
                                // for use in enforcement
                                // NKC, note, almost certainly possible to 
                                // combine some of these!
  int primary_proc_array_index; // This is an object's 0..N-1 index 
                                // on a processor in the primary topology
  int owner_proc_array_index;   // This is an object's 0..N-1 index 
                                // on a processor in the primary topology
  int block_id;                 // This is the index of the block 
                                // this entity belongs to
  
  unsigned int context;         // This is searching status of the entity.
                                // Ie if it is in a tracked interaction or if
                                // it should be searched normally.


  ContactType base_type;
  DataType *data_array;

  ContactHostGlobalID global_id;
  bool shared;
  ProcessorOwnership ownership;
  int owner;
  int secondary_owner;
  int entity_key;

 private:
  // not defined; all ContactEntities are not copyable or assignable
  ContactTopologyEntity( const ContactTopologyEntity& );
  ContactTopologyEntity& operator=( const ContactTopologyEntity& );
  ContactTopologyEntity( ContactTopologyEntity&);

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
template<typename DataType>
inline int ContactTopologyEntity<DataType>::Size(const int data_array_length)
{
  return( NUMBER_PACKED_VARS*sizeof(int)+data_array_length*sizeof(DataType) );
}

template<typename DataType>
inline void ContactTopologyEntity<DataType>::Pack( char* buffer, const int data_array_length )
{
  // Note that the first location (0) is reserved for the entity type.
  // and the second location (1) is reserved for derived type.
  // The later will be filled in by the derived classes.
  int* i_buffer = reinterpret_cast<int*>(buffer);
  i_buffer[BASE_TYPE]        = base_type;
  i_buffer[GID_HI]           = global_id.HiInt();
  i_buffer[GID_LO]           = global_id.LoInt();
  i_buffer[BLOCK_ID]         = block_id;
  i_buffer[CONTEXT]          = context;
  i_buffer[PRIM_OWNER]       = owner;
  i_buffer[SEC_OWNER]        = secondary_owner;
  i_buffer[OWNER_PROC_INDEX] = owner_proc_array_index;
  i_buffer[PROC_INDEX]       = proc_array_index;
  i_buffer[PRIM_PROC_INDEX]  = primary_proc_array_index;
  i_buffer[HOST_INDEX]       = host_array_index;
  i_buffer[OWNERSHIP]        = ownership;
  DataType* buf_loc = reinterpret_cast<DataType*> 
    (buffer+NUMBER_PACKED_VARS*sizeof(int));
  std::memcpy(buf_loc,DataArray_Buffer(),data_array_length*sizeof(DataType));
}

template<typename DataType>
inline void ContactTopologyEntity<DataType>::Unpack( char* buffer, const int data_array_length )
{
  //We're unpacking an int into an enum. This will fail if the compiler doesn't
  //use ints to store enums. So preconditioning that here.
  //PRECONDITION(sizeof(ActiveStatus) == sizeof(int));
  
  int* i_buffer = reinterpret_cast<int*>(buffer);
  PRECONDITION( i_buffer[BASE_TYPE] == base_type );
  global_id.HiInt( i_buffer[GID_HI] );
  global_id.LoInt( i_buffer[GID_LO] );
  block_id                 = i_buffer[BLOCK_ID];
  context                  = i_buffer[CONTEXT];
  owner                    = i_buffer[PRIM_OWNER];
  secondary_owner          = i_buffer[SEC_OWNER];
  ownership                = (ProcessorOwnership)(i_buffer[OWNERSHIP]);
  host_array_index         = i_buffer[HOST_INDEX];
  proc_array_index         = i_buffer[PROC_INDEX];
  primary_proc_array_index = i_buffer[PRIM_PROC_INDEX];
  owner_proc_array_index   = i_buffer[OWNER_PROC_INDEX];
  DataType* buf_loc            = reinterpret_cast<DataType*> 
    (buffer+NUMBER_PACKED_VARS*sizeof(int));
  std::memcpy(DataArray_Buffer(),buf_loc,data_array_length*sizeof(DataType));

  //
  //  NKC note, should investigate making Pack / Unpack consitently handle
  //  owner and ownership attributes
  //
  //  int my_proc;
  //MPI_Comm_rank( MPI_COMM_WORLD, &my_proc );
  //if(my_proc != owner) ownership = NOT_OWNED;
}

template<typename DataType>
inline void ContactTopologyEntity<DataType>::Copy( ContactTopologyEntity* src, const int data_array_length )
{
  global_id                = src->global_id;
  block_id                 = src->block_id;
  context                  = src->context;
  owner                    = src->owner;
  secondary_owner          = src->secondary_owner;
  ownership                = src->ownership;
  host_array_index         = src->host_array_index;
  proc_array_index         = src->proc_array_index;
  primary_proc_array_index = src->primary_proc_array_index;
  owner_proc_array_index   = src->owner_proc_array_index;
  std::memcpy(DataArray_Buffer(),src->DataArray_Buffer(),data_array_length*sizeof(DataType));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
template<typename DataType>
inline int ContactTopologyEntity<DataType>::Size_ForSecondary(const int data_array_length)
{
  int num_packed_vars = PROC_INDEX;
  return( num_packed_vars*sizeof(int)+data_array_length*sizeof(DataType) );
}

template<typename DataType>
inline void 
ContactTopologyEntity<DataType>::Pack_ForSecondary( char* buffer, 
                                          const int data_array_length )
{
  int num_packed_vars = PROC_INDEX;
  // Note that the first location (0) is reserved for the entity type.
  // and the second location (1) is reserved for derived type.
  // The later will be filled in by the derived classes.
  int* i_buffer = reinterpret_cast<int*> (buffer);
  i_buffer[BASE_TYPE]        = base_type;
  i_buffer[GID_HI]           = global_id.HiInt();
  i_buffer[GID_LO]           = global_id.LoInt();
  i_buffer[BLOCK_ID]         = block_id;
  i_buffer[CONTEXT]          = context;
  i_buffer[PRIM_OWNER]       = owner;
  i_buffer[SEC_OWNER]        = secondary_owner;
  i_buffer[OWNER_PROC_INDEX] = owner_proc_array_index;
  i_buffer[HOST_INDEX]       = host_array_index;
  DataType* buf_loc = reinterpret_cast<DataType*> 
    (buffer+num_packed_vars*sizeof(int));
  std::memcpy(buf_loc,DataArray_Buffer(),data_array_length*sizeof(DataType));
}

template<typename DataType>
inline void 
ContactTopologyEntity<DataType>::Unpack_ForSecondary( char* buffer,  
                                            const int data_array_length )
{
  int num_packed_vars = PROC_INDEX;
  int* i_buffer = reinterpret_cast<int*> (buffer);
  PRECONDITION( i_buffer[BASE_TYPE] == base_type );
  global_id.HiInt( i_buffer[GID_HI] );
  global_id.LoInt( i_buffer[GID_LO] );
  block_id               = i_buffer[BLOCK_ID];
  context                = i_buffer[CONTEXT];
  owner                  = i_buffer[PRIM_OWNER];
  secondary_owner        = i_buffer[SEC_OWNER];
  owner_proc_array_index = i_buffer[OWNER_PROC_INDEX];
  host_array_index       = i_buffer[HOST_INDEX];
  DataType* buf_loc   = reinterpret_cast<DataType*> 
    (buffer+num_packed_vars*sizeof(int));
  std::memcpy(DataArray_Buffer(),buf_loc,data_array_length*sizeof(DataType));
}

template<typename DataType>
inline void 
ContactTopologyEntity<DataType>::Pack_ForSecondary( char* buffer, 
                                          DataType* DataArray, 
                                          const int data_array_length )
{
  int num_packed_vars = PROC_INDEX;
  // Note that the first location (0) is reserved for the entity type.
  // and the second location (1) is reserved for derived type.
  // The later will be filled in by the derived classes.
  int* i_buffer = reinterpret_cast<int*> (buffer);
  i_buffer[BASE_TYPE]        = base_type;
  i_buffer[GID_HI]           = global_id.HiInt();
  i_buffer[GID_LO]           = global_id.LoInt();
  i_buffer[BLOCK_ID]         = block_id;
  i_buffer[CONTEXT]          = context;
  i_buffer[PRIM_OWNER]       = owner;
  i_buffer[SEC_OWNER]        = secondary_owner;
  i_buffer[OWNER_PROC_INDEX] = owner_proc_array_index;
  i_buffer[HOST_INDEX]       = host_array_index;
  DataType* buf_loc = reinterpret_cast<DataType*> 
    (buffer+num_packed_vars*sizeof(int));
  std::memcpy(buf_loc,DataArray,data_array_length*sizeof(DataType));
}

template<typename DataType>
inline void 
ContactTopologyEntity<DataType>::Unpack_ForSecondary( char* buffer, 
                                            DataType* DataArray, 
                                            const int data_array_length )
{
  int num_packed_vars = PROC_INDEX;
  int* i_buffer = reinterpret_cast<int*> (buffer);
  PRECONDITION( i_buffer[BASE_TYPE] == base_type );
  global_id.HiInt( i_buffer[GID_HI] );
  global_id.LoInt( i_buffer[GID_LO] );
  block_id               = i_buffer[BLOCK_ID];
  context                = i_buffer[CONTEXT];
  owner                  = i_buffer[PRIM_OWNER];
  secondary_owner        = i_buffer[SEC_OWNER];
  owner_proc_array_index = i_buffer[OWNER_PROC_INDEX];
  host_array_index       = i_buffer[HOST_INDEX];
  DataType* buf_loc   = reinterpret_cast<DataType*> 
    (buffer+num_packed_vars*sizeof(int));
  std::memcpy(DataArray,buf_loc,data_array_length*sizeof(DataType));
}

template<typename DataType>
inline void ContactTopologyEntity<DataType>::Copy_ForSecondary( ContactTopologyEntity* src, const int data_array_length )
{
  base_type                = src->base_type;
  global_id                = src->global_id;
  block_id                 = src->block_id;
  context                  = src->context;
  owner                    = src->owner;
  secondary_owner          = src->secondary_owner;
  owner_proc_array_index   = src->owner_proc_array_index;
  host_array_index         = src->host_array_index;
  std::memcpy(DataArray_Buffer(),src->DataArray_Buffer(),data_array_length*sizeof(DataType));
}

template<typename DataType>
inline int ContactTopologyEntity<DataType>::Size_ForDataUpdate(const int data_array_length)
{
  return( data_array_length*sizeof(DataType)+sizeof(int) );
}

template<typename DataType>
inline void ContactTopologyEntity<DataType>::Pack_ForDataUpdate( char* buffer, const int data_array_length )
{
  int* i_buffer = reinterpret_cast<int*>(buffer);
  //i_buffer[0]   = Base_Type();
  i_buffer[0]   = block_id;
  DataType* buf_loc = reinterpret_cast<DataType*>(buffer+sizeof(int));
  std::memcpy(buf_loc, DataArray_Buffer(), data_array_length*sizeof(DataType));
}

template<typename DataType>
inline void ContactTopologyEntity<DataType>::Unpack_ForDataUpdate( char* buffer, const int data_array_length )
{
  DataType* buf_loc = reinterpret_cast<DataType*>(buffer+sizeof(int));
  std::memcpy(DataArray_Buffer(), buf_loc, data_array_length*sizeof(DataType));
}

template<typename DataType>
inline int ContactTopologyEntity<DataType>::PackConnection( ContactTopologyEntity* entity,
                                                  int* buffer)
{
  int* buf = buffer;
  if (entity) {
    *buf++ = entity->owner; 
    *buf++ = entity->owner_proc_array_index;
    *buf++ = entity->global_id.HiInt(); 
    *buf++ = entity->global_id.LoInt(); 
  } else {
    buf += 4;
  }
  return 4;
}

template<typename DataType>
inline int ContactTopologyEntity<DataType>::PackConnection( connection_data* data,
                                                  int* buffer)
{
  int* buf = buffer;
  *buf++   = data->owner; 
  *buf++   = data->owner_proc_array_index;
  *buf++   = data->host_gid[0]; 
  *buf++   = data->host_gid[1]; 
  return 4;
}

template<typename DataType>
inline void ContactTopologyEntity<DataType>::PackConnection( ContactTopologyEntity* entity,
                                                   connection_data* data)
{
  if (entity) {
    data->owner                  = entity->owner;
    data->owner_proc_array_index = entity->owner_proc_array_index;
    data->host_gid[0]            = entity->global_id.HiInt();
    data->host_gid[1]            = entity->global_id.LoInt();
  }
}

template<typename DataType>
inline int ContactTopologyEntity<DataType>::UnPackConnection( connection_data* data, 
                                                    int* buffer)
{
  int* buf                     = buffer;
  data->owner                  = *buf++;
  data->owner_proc_array_index = *buf++;
  data->host_gid[0]            = *buf++;
  data->host_gid[1]            = *buf++;
  return 4;
}



#endif  // ifdef ContactTopologyEntity_H_
