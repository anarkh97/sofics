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


#ifndef ContactNodeNodeInteraction_h_
#define ContactNodeNodeInteraction_h_

#include "Contact_Defines.h"
#include "ContactInteractionEntity.h"
#include "ContactNode.h"

class CString;
class ContactTopologyEntityList;
class ContactFixedSizeAllocator;
class ContactTopology;

class ContactNodeNodeInteraction : public ContactInteractionEntity<Real> {

 public:

  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  NNI_SCALAR_VAR
#define NNI_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
                     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactNodeNodeInteraction( ContactNode<Real>*, ContactNode<Real>*,
                              int, Real );
  ContactNodeNodeInteraction( ContactNodeNodeInteraction& nni );
  ContactNodeNodeInteraction();
  static ContactNodeNodeInteraction* new_ContactNodeNodeInteraction(
            ContactFixedSizeAllocator&, 
            ContactNode<Real>*, ContactNode<Real>*,
            int, Real );
  static ContactNodeNodeInteraction* new_ContactNodeNodeInteraction(
             ContactFixedSizeAllocator& );     
  ~ContactNodeNodeInteraction();

  inline ContactNode<Real>* SlaveNode() { return slave_node; };
  inline entity_data* SlaveNodeEntityData() {return &slave_node_entity_data;};
  int Set_SlaveNodeEntityData();
  inline ContactNode<Real>* MasterNode() { return master_node; };
  inline entity_data* MasterNodeEntityData() {return &master_node_entity_data;};
  int Set_MasterNodeEntityData();

  inline Real& Scalar_Var( VariableHandle vh ) {return DataArray[vh];};
  inline Real* Vector_Var( VariableHandle vh ) 
    { return (DataArray+NUMBER_SCALAR_VARS+3*vh); };

  // Parallel pack/unpack functions
  inline int   Size();
  inline char* Pack( char* buffer );
  inline char* Unpack( char* buffer );
  inline void  Copy( ContactNodeNodeInteraction* src );

  // Restart Pack/Unpack functions
  int  Restart_Size();
  void Restart_Pack( Real* buffer );
  void Restart_Unpack( Real* buffer );

  void Connect_SlaveNode( ContactTopologyEntityList& );
  void Connect_SlaveNode( ContactNode<Real>* );
  void Connect_MasterNode( ContactTopologyEntityList& );
  void Connect_MasterNode( ContactNode<Real>* );

  inline bool     Is_Tied() { return is_tied; };
  inline void     Is_Tied( bool IS_TIED ) { is_tied = IS_TIED; };

  inline int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline void Initialize_Memory() {std::memset(DataArray_Buffer(), 0, DataArray_Length()*sizeof(Real));};
 private:
  
  ContactNode<Real>* slave_node;
  entity_data  slave_node_entity_data;
  ContactNode<Real>* master_node;
  entity_data  master_node_entity_data;
  
  Real DataArray[NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS];

  bool is_tied;
};

inline int ContactNodeNodeInteraction::Size()
{
  return(ContactInteractionEntity<Real>::Size()+2*sizeof(entity_data)+
         (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real));
}

inline char* ContactNodeNodeInteraction::Pack( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Pack( buffer );
  buff += ContactInteractionEntity<Real>::Size();
  
  int  cnt    = 0;
  int* i_buf  = reinterpret_cast<int*>(buff);
  cnt        += PackEntityData(&slave_node_entity_data, &i_buf[cnt]);
  cnt        += PackEntityData(&master_node_entity_data, &i_buf[cnt]);
  
  buff += cnt*sizeof(int);
  cnt   = (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real);
  std::memcpy( buff,DataArray,cnt);
  buff += cnt;
  return buff;
}

inline char* ContactNodeNodeInteraction::Unpack( char* buffer )
{
  char* buff = buffer;
  ContactInteractionEntity<Real>::Unpack( buffer );
  buff += ContactInteractionEntity<Real>::Size();
  
  int cnt    = 0;
  int* i_buf = reinterpret_cast<int*>(buff);
  cnt       += UnPackEntityData(&slave_node_entity_data, &i_buf[cnt]);
  cnt       += UnPackEntityData(&master_node_entity_data, &i_buf[cnt]);
  
  buff += cnt*sizeof(int);
  cnt   = (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real);
  std::memcpy( DataArray,buff,cnt);
  buff += cnt;
  return buff;
}

inline void ContactNodeNodeInteraction::Copy( ContactNodeNodeInteraction* src )
{
  ContactInteractionEntity<Real>::Copy( src );
  slave_node_entity_data  = src->slave_node_entity_data;
  master_node_entity_data = src->master_node_entity_data;
  std::memcpy( DataArray,src->DataArray_Buffer(),
          (NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS)*sizeof(Real) );
}


#endif
