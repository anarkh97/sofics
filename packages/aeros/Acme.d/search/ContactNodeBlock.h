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


#ifndef ContactNodeBlock_h_
#define ContactNodeBlock_h_

#include "ContactNode.h"
#include "ContactSearch.h"
#include "ContactBlockEntityList.h"
#include "ContactBoundingBox.h"
#include "ContactParOStream.h"

class ContactTopology;

class ContactNodeBlock {

 public:

   ContactNodeBlock(ContactSearch::ContactNode_Type, int, int, int&, int*, 
                    int*, ContactType*, ContactTopology* );
   ContactNodeBlock(ContactSearch::ContactNode_Type, int, ContactTopology* );
  ~ContactNodeBlock();

  void UpdateNodeBlock( int, int, int&, ContactType* );
  void Allocate_DataArray( int number_scalar_vars, int number_vector_vars );

  inline ContactSearch::ContactNode_Type Type() { return type;};
  inline void Has_Attributes(int i) { has_attributes=i;};
  inline int Has_Attributes() { return has_attributes;};
  inline void Has_Radius_Attributes(int i) { has_radius_attributes=i;};
  inline int Has_Radius_Attributes() { return has_radius_attributes;};
  inline void Has_Normal_Attributes(int i) { has_normal_attributes=i;};
  inline int Has_Normal_Attributes() { return has_normal_attributes;};
  inline int ID() { return id; };
  inline int Number_of_Nodes() { return number_of_nodes; };
  ContactBlockEntityList* NodeList() { return node_list; };
  ContactBoundingBox* LocalBoundingBox() { return &local_bounding_box; };
  ContactBoundingBox* GlobalBoundingBox() { return &global_bounding_box; };
  void ComputeBoundingBox(int, VariableHandle, VariableHandle, VariableHandle, MPI_Comm&);
  
  void Delete_Node_List( );
  void Delete_Nodes( );
  void Add_Nodes( int, int*, int*, ContactType*  );
  void Insert_Node( ContactNode<Real>* );
  void Delete_Node( ContactNode<Real>* );
#ifndef CONTACT_NO_MPI
  void Insert_Node( char* );
  void Insert_Node_ForSecondary( char* );
#endif
  inline void Rmax(Real r) { rmax=r; };
  inline Real Rmax() { return rmax; };
  inline void Entity_Key(int key) {entity_key=key;};
  inline int  Entity_Key() {return entity_key;};
  
  void Display(ContactParOStream&);

  void SetMaster() {master=true;};
  void SetSlave()  {slave=true;};
  void SetMaster(bool val) {master=val;};
  void SetSlave(bool val)  {slave=val;};
  bool IsMaster( )         {return master;};
  bool IsSlave( )          {return slave;};

 private:
  
  ContactNodeBlock(ContactNodeBlock&);
  ContactNodeBlock& operator=(ContactNodeBlock&);
  
  bool master;
  bool slave;
  int entity_key;

  ContactTopology* topology;  // owner topology object

  int number_of_nodes;
  int num_nodes_added; // This is used to place new nodes in Add_Node
  ContactSearch::ContactNode_Type type;
  int id;
  int has_attributes;
  int has_radius_attributes;
  int has_normal_attributes;
  Real rmax;
  
  ContactBlockEntityList* node_list;
  
  ContactBoundingBox local_bounding_box;
  ContactBoundingBox global_bounding_box;

};

#endif  // #ifdef ContactNodeBlock_h_
