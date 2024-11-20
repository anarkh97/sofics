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


// We currently have two kinds of "elements" in ACME std::right now. 
// The first is the "extrusion" of a face to an element (i.e., a T3 to Wedge6
//     or a Q4 to Hex8) used for face-face searches
// The second is an true element that only has node connections (nodes do not
//     have back pointers to the elements std::right now).  These are used for the
//     volume overlap searches.

#ifndef ContactElement_h_
#define ContactElement_h_

#include "ContactEntity.h"
#include "contact_assert.h"
#include "ContactEdge.h"
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactDoublyLinkedList.h"

template<typename DataType> class ContactNode;
template<typename DataType> class ContactEdge;
template<typename DataType> class ContactFace;
class ContactElementElementInteraction;

template<typename DataType>
class ContactElem : public ContactTopologyEntity<DataType> {

 public:

  using ContactTopologyEntity<DataType>::DataArray_Buffer;
  
  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  ELEM_SCALAR_VAR
#define ELEM_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  ELEM_VECTOR_VAR
#define ELEM_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactElem( ContactSearch::ContactElem_Type, 
	       int Block_Index, int Host_Index_in_Block, int key);
  virtual ~ContactElem();
  
  virtual void BuildTopology(int, int, int, ContactFixedSizeAllocator*) = 0;
  virtual void DeleteTopology(ContactFixedSizeAllocator*) = 0;
  virtual void UpdateTopology(ContactFace<DataType>*, VariableHandle, VariableHandle,
                              VariableHandle, Real, bool use_node_normals=false) = 0;

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  virtual int Nodes_Per_Element() = 0;
  virtual int Edges_Per_Element() = 0;
  virtual int Faces_Per_Element() = 0;
  virtual void Evaluate_Shape_Functions( DataType*, DataType* ) = 0;
  virtual void Compute_Global_Coordinates( VariableHandle, DataType*, DataType* ) = 0;
  virtual bool Compute_Local_Coordinates( VariableHandle, DataType*, DataType* ) = 0;
  virtual bool Is_Local_Coordinates_Inside_Element( DataType* ) = 0;
  virtual bool Is_Local_Coordinates_Near_Element( DataType*, DataType ) = 0;
  virtual void Compute_Partial_Face_Normal( int, VariableHandle, VariableHandle,
                                            DataType (*)[3], Real, DataType (*)[3], DataType * )
               { std::cerr << "ContactElem::Compute_Partial_Face_Normal is not implemented\n"; }
  virtual void Compute_Second_Partial_Face_Normal( int, VariableHandle, VariableHandle, DataType (*)[3],
                                                   DataType (*)[3], Real, DataType (*)[3], DataType * )
               { std::cerr << "ContactElem::Compute_Second_Partial_Face_Normal is not implemented\n"; }
  
  virtual ContactSearch::ContactNode_Type Node_Type() = 0;
  virtual ContactSearch::ContactEdge_Type Edge_Type() = 0;
  virtual ContactSearch::ContactFace_Type Face_Type(int) = 0;
  inline  ContactSearch::ContactElem_Type Elem_Type() {return type;};

  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(DataType));};

  virtual ContactNode<DataType>** Nodes() = 0;
  virtual ContactEdge<DataType>** Edges() = 0;
  virtual ContactFace<DataType>** Faces() = 0;

  ContactNode<DataType>* Node( int i );
  ContactEdge<DataType>* Edge( int i );
  ContactFace<DataType>* Face( int i );
  
  void ConnectNode( int,ContactNode<DataType>* );
  void ConnectEdge( int,ContactEdge<DataType>* );
  void ConnectFace( int,ContactFace<DataType>* );

  DataType MaxSize(VariableHandle POSITION);
  // Packing/Unpacking Functions
  int  Size();
  void Pack( char* );
  void Unpack( char* );

  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( DataType* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(DataType) );
  }
  inline void Restart_Unpack( DataType* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(DataType) );
  }

  virtual int* Node_Ids() = 0;
  virtual int* Edge_Ids() = 0;
  virtual int* Face_Ids() = 0;

 protected:
  using ContactTopologyEntity<DataType>::entity_key;
  using ContactTopologyEntity<DataType>::block_id;
  ContactSearch::ContactElem_Type type;

 private:
  DataType DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
};

template<typename DataType>
inline ContactNode<DataType>* ContactElem<DataType>::Node( int i )
{ 
  PRECONDITION( i>=0 && i<Nodes_Per_Element() );
  return( Nodes()[i] );
}

template<typename DataType>
inline ContactEdge<DataType>* ContactElem<DataType>::Edge( int i )
{ 
  PRECONDITION( i>=0 && i<Edges_Per_Element() );
  return( Edges()[i] );
}

template<typename DataType>
inline ContactFace<DataType>* ContactElem<DataType>::Face( int i )
{ 
  PRECONDITION( i>=0 && i<Faces_Per_Element() );
  return( Faces()[i] );
}




class ContactElement : public ContactTopologyEntity<Real> {

 public:

  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  ELEMENT_SCALAR_VAR
#define ELEMENT_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  ELEMENT_VECTOR_VAR
#define ELEMENT_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactElement( ContactFixedSizeAllocator*,
                  ContactSearch::ContactElement_Type, 
		  int Block_Index, int Host_Index_in_Block, int key);
  virtual ~ContactElement();
  
  virtual int Nodes_Per_Element() = 0;
  virtual void Compute_Volume(VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Centroid(VariableHandle, VariableHandle ) = 0;
  virtual void TetDice(int &ntets, Real thex[][4][3], VariableHandle) = 0;

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};

  virtual void Evaluate_Shape_Functions( Real*, Real* ) = 0;
  virtual void Compute_Global_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual void Compute_Local_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual bool Is_Local_Coordinates_Inside_Element( Real* ) = 0;
  virtual bool Is_Local_Coordinates_Near_Element( Real*, Real ) = 0;
  virtual void Interpolate_Scalar_Value( Real *, Real *, Real& ) { return; };

  inline ContactSearch::ContactElement_Type ElementType() {return type;};
  virtual ContactNode<Real>** Nodes() = 0;
  ContactNode<Real>* Node( int i );  
  void ConnectNode( int, ContactNode<Real>* );

  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( Real* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(Real) );
  }
  inline void Restart_Unpack( Real* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(Real) );
  }

#ifndef CONTACT_NO_MPI

  Real MaxSize(VariableHandle POSITION);
  
  // Packing/Unpacking Functions
  inline int  Size();
  inline void Pack( char* );
  inline void Unpack( char* );
  inline void Copy( ContactElement* );
  
  inline int  Size_ForSecondary();
  inline void Pack_ForSecondary( char* );
  inline void Unpack_ForSecondary( char* );
  inline void Copy_ForSecondary( ContactElement* );
  
  inline int  Size_ForDataUpdate();
  inline void Pack_ForDataUpdate( char* );
  inline void Unpack_ForDataUpdate( char* );
  
  virtual int  Size_Interactions( int state=0 );
  virtual void Pack_Interactions( char*, int state=0 );
  virtual void Unpack_Interactions( char*, int state=0 );
  virtual void Copy_Interactions( ContactElement*, int state=0 );
  
  virtual int  Size_Interactions_ForSecondary( int state=0 );
  virtual void Pack_Interactions_ForSecondary( char*, int state=0 );
  virtual void Unpack_Interactions_ForSecondary( char*, int state=0 );
  virtual void Copy_Interactions_ForSecondary( ContactElement*, int state=0 );

#endif

  virtual connection_data* NodeInfo() = 0;
  
  inline int Number_Interactions(int state = 0) 
      { int n = ElementElementInteractions[state]->NumEntities();
        return n; };
  
  inline int Number_ElementElement_Interactions (int state = 0) 
      { return ElementElementInteractions[state]->NumEntities(); };
                                         
  inline ContactInteractionDLL<Real>* Get_ElementElement_Interactions( int state = 0 )
       { return ElementElementInteractions[state]; };
      
  ContactElementElementInteraction* Get_ElementElement_Interaction(int interaction_number,
						                   int state = 0 );
       
  void Store_ElementElement_Interaction( ContactElementElementInteraction*, 
					 int state = 0 );
                                         
  void Delete_ElementElement_Interaction( ContactElementElementInteraction*, 
			 		  int state = 0 );
                                         
  void Display_ElementElement_Interactions( ContactParOStream&, int state = 0 );
  
  void Update_Interactions( );

  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(Real));};
  
  void ComputeBoundingBoxForSearch(const int num_configs,
                                   const VariableHandle &NODE_COORD_START,
                                   const VariableHandle &NODE_COORD_END,
                                   const int  auto_tol,
                                   const Real box_inflation,
                                   const Real user_tol,
                                   ContactBoundingBox &box);

 protected:
  ContactSearch::ContactElement_Type type;
  ContactFixedSizeAllocator* allocators;

 private:
  Real DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
  
  int number_of_states;
  ContactInteractionDLL<Real>** ElementElementInteractions;
};

inline ContactNode<Real>* ContactElement::Node( int i )
{ 
  PRECONDITION( i>=0 && i<Nodes_Per_Element() );
  return( Nodes()[i] );
}

inline void ContactElement::ConnectNode( int i, ContactNode<Real>* node ) {
  PRECONDITION( i>=0 && i<Nodes_Per_Element() );
  Nodes()[i] = node;
  PackConnection((ContactTopologyEntity<Real>*)node, &NodeInfo()[i]);
}

#ifndef CONTACT_NO_MPI
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used 
// to just update the data for ghosted elements
//--------------------------------------------------------------------
inline int ContactElement::Size_ForDataUpdate()
{
  return( ContactTopologyEntity<Real>::Size_ForDataUpdate(DataArray_Length()) );
}

inline void ContactElement::Pack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity<Real>::Pack_ForDataUpdate( buffer, DataArray_Length() );
}

inline void ContactElement::Unpack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity<Real>::Unpack_ForDataUpdate( buffer, DataArray_Length() );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
inline int ContactElement::Size()
{
  return( ContactTopologyEntity<Real>::Size(DataArray_Length()) + 
          Nodes_Per_Element()*sizeof(connection_data) );
}

inline void ContactElement::Pack( char* buffer )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  // ContactTopologyEntity<Real> packs in location 0 as ContactElement 
  // and here we pack in the derived type in location 1.
  i_buf[1] = type;
  ContactTopologyEntity<Real>::Pack( buffer, DataArray_Length() );
  // Add the global ids of the nodes
  i_buf = reinterpret_cast<int*> (buffer+ContactTopologyEntity<Real>::Size(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    cnt += PackConnection(reinterpret_cast<ContactTopologyEntity<Real>*>(Node(i)), &i_buf[cnt]);
  }
}

inline void ContactElement::Unpack( char* buffer )
{
  ContactTopologyEntity<Real>::Unpack( buffer, DataArray_Length() );
  entity_key = block_id;

  PRECONDITION( ((int*)buffer)[1] == type );
  // Store off the global ids of the nodes
  int* i_buf = reinterpret_cast<int*> ( buffer + ContactTopologyEntity<Real>::Size(DataArray_Length()) );
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    cnt += UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
}

inline void ContactElement::Copy( ContactElement* src )
{
  ContactTopologyEntity<Real>::Copy( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    PackConnection(reinterpret_cast<ContactTopologyEntity<Real>*>(src->Node(i)), &NodeInfo()[i]);
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
inline int ContactElement::Size_ForSecondary()
{
  return( ContactTopologyEntity<Real>::Size_ForSecondary(DataArray_Length()) + 
          Nodes_Per_Element()*sizeof(connection_data) );
}

inline void ContactElement::Pack_ForSecondary( char* buffer )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  // ContactTopologyEntity<Real> packs in location 0 as ContactElement 
  // and here we pack in the derived type in location 1.
  i_buf[1] = type;
  ContactTopologyEntity<Real>::Pack_ForSecondary( buffer, DataArray_Length() );
  // Add the global ids of the nodes
  i_buf = reinterpret_cast<int*> (buffer+ContactTopologyEntity<Real>::Size_ForSecondary(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    cnt += PackConnection(reinterpret_cast<ContactTopologyEntity<Real>*>(Node(i)), &i_buf[cnt]);
  }
}

inline void ContactElement::Unpack_ForSecondary( char* buffer )
{
  ContactTopologyEntity<Real>::Unpack_ForSecondary( buffer, DataArray_Length() );
  entity_key = block_id;

  PRECONDITION( ((int*)buffer)[1] == type );
  // Store off the global ids of the nodes
  int* i_buf = reinterpret_cast<int*> ( buffer + ContactTopologyEntity<Real>::Size_ForSecondary(DataArray_Length()) );
  int cnt = 0;
  for( int i=0; i<Nodes_Per_Element(); ++i ){
    cnt += UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
}

inline void ContactElement::Copy_ForSecondary( ContactElement* src )
{
  ContactTopologyEntity<Real>::Copy_ForSecondary( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0; i<Nodes_Per_Element(); ++i ){
    PackConnection(reinterpret_cast<ContactTopologyEntity<Real>*>(src->Node(i)), &NodeInfo()[i]);
  }
}
#endif

#endif // #ifdef ContactElement_h_
