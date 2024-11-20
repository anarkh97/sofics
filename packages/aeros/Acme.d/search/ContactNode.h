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


#ifndef ContactNode_h_
#define ContactNode_h_

/* A Node is a Node is a Node so this is a concrete class unlike
   ContactEdge and ContactFace which are themselves abstract classes */

#include "ContactTopologyEntity.h"
#include "contact_assert.h"
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactEdge.h"
#include "ContactBoundingBox.h"

template<typename DataType> class ContactFace;
class ContactNodeNodeInteraction;
class ContactNodeEntityInteraction;
class ContactNodeFaceInteraction;
class ContactNodeSurfaceInteraction;

template<typename DataType>
class ContactNode : public ContactTopologyEntity<DataType> {

 public:

  using ContactTopologyEntity<DataType>::DataArray_Buffer;
  using ContactTopologyEntity<DataType>::Base_Type;
  using ContactTopologyEntity<DataType>::Variable;

  enum PHYSICAL_TYPE { CONTINUUM_NODE, MIXED_NODE, SHELL_NODE, SHELL_TAB_NODE,
                       MIXED_TAB_NODE };

  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  NODE_SCALAR_VAR
#define NODE_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  NODE_VECTOR_VAR
#define NODE_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactNode( ContactFixedSizeAllocator*,
               ContactSearch::ContactNode_Type, 
               int Block_Index=-1,
	       int Index_in_Block=-1, 
	       ContactType cetype=CT_NODE);
               
  static ContactNode* new_ContactNode( ContactFixedSizeAllocator*,
				       ContactSearch::ContactNode_Type, 
				       int Block_Index=-1,
				       int Index_in_Block=-1,
				       ContactType 
				       cetype=CT_NODE); 
                                       
  virtual ~ContactNode();

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};

  inline ContactFace<DataType>* GetFace( int i );
  inline int          GetFacePFIndex(int i);
  int                 GetFacePFEntityKey(int physical_face_num);
  bool                ConnectedToFace(const ContactHostGlobalID &id);
  void                SortConnectedFaces();
  inline std::vector<std::pair<ContactFace<DataType>*, int> >& ConnectedFaces() { return faces; };



  // Functions to connect up the topology
  void Delete_Face_Connections( );
  void Connect_Face(ContactFace<DataType>* );

  inline int Number_Face_Connections() {return faces.size();};
  inline ContactSearch::ContactNode_Type NodeType() {return node_type;};
  
  //
  //  Determine the entity key for the current node.  Return the entity number
  //  (greater than or equal to zero) if the node belongs to one and only
  //  one entity.  Return -1 if the node belongs to multiple entities.
  //
  int Get_Owning_Entity();
         
  // Packing/Unpacking Functions
  // state = -2: don't include any interactions
  // state = -1: include all interactions
  // state >= 0: include interactions for state #
  inline virtual int  Size() { return Size(-2); }
  inline virtual int  Size(int state );
  inline virtual void Pack( char* buffer ) { Pack(buffer, -2); }
  inline virtual void Pack( char*, int state );
  inline virtual void Unpack( char* );
  inline virtual void Copy( ContactNode* node ) { Copy(node, -2); }
  inline virtual void Copy( ContactNode*, int state );
  
  inline virtual int  Size_ForSecondary() 
           { return Size_ForSecondary(-2); }
  inline virtual int  Size_ForSecondary(int state);
  inline virtual void Pack_ForSecondary( char* buffer )
           { Pack_ForSecondary(buffer, -2); }
  inline virtual void Pack_ForSecondary( char*, int state );
  inline virtual void Unpack_ForSecondary( char* );

  inline virtual void Copy_ForSecondary( ContactNode* node) 
           { Copy_ForSecondary(node, -2); }
  inline virtual void Copy_ForSecondary( ContactNode*, int state );
  
  inline virtual int  Size_ForDataUpdate( );
  inline virtual void Pack_ForDataUpdate( char* );
  inline virtual void Unpack_ForDataUpdate( char* );

  virtual int  Size_Interactions( int state );
  virtual void Pack_Interactions( char*, int state );
  virtual void Unpack_Interactions( char*, int state );
  virtual void Copy_Interactions( ContactNode*, int state );
  
  virtual int  Size_Interactions_ForSecondary( int state );
  virtual void Pack_Interactions_ForSecondary( char*, int state );
  virtual void Unpack_Interactions_ForSecondary( char*, int state );
  virtual void Copy_Interactions_ForSecondary( ContactNode*, int state );
  
  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( DataType* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(DataType) );
  }
  inline void Restart_Unpack( DataType* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(DataType) );
  }

  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(DataType));};

  inline bool Is_a_Shell_Node() { 
    if( Base_Type() == CT_NODE ) return false;
    return true;
  };

  inline void Exodus_ID( int e_id ) {exodus_id=e_id;};
  inline int Exodus_ID(){ return exodus_id; };

  PHYSICAL_TYPE Physical_Type() { return physical_type; };
  void Physical_Type( PHYSICAL_TYPE pt ) { physical_type = pt; };
  
  inline int Number_Interactions(int state = 0) {
    int num_nni = 0;
    if((int)NodeNodeInteractions.size() > state) num_nni = NodeNodeInteractions[state].NumEntities();
    int num_nei = 0;
    if((int)NodeEntityInteractions.size() > state) num_nei = NodeEntityInteractions[state].size();
    return num_nni + num_nei;
  }
      
  inline int Number_NodeNode_Interactions(int state = 0) { 
    if((int)NodeNodeInteractions.size() <= state) return 0;
    return NodeNodeInteractions[state].NumEntities();
  };
      
  inline ContactInteractionDLL<Real>* Get_NodeNode_Interactions( int state = 0 ) {
    if((int)NodeNodeInteractions.size() <= state) return NULL;
    return &(NodeNodeInteractions[state]);
   };
       
  ContactNodeNodeInteraction* 
    Get_NodeNode_Interaction(int interaction_number, int state = 0 );

  void Add_NodeNode_Interaction( ContactNodeNodeInteraction* cnni,
				 int state = 0 );
                                   
  void Delete_NodeNode_Interaction( ContactNodeNodeInteraction* cnni,
				    int state = 0 );
                                   
  void Delete_NodeNode_Interactions( int state = 0 );

  void Print_NodeNode_Interactions( int state = 0 );

  void Display_NodeNode_Interactions( ContactParOStream&, int state = 0 );
                    
  int Number_NodeFace_Interactions(const int state = 0) const;
  int Number_NodeSurface_Interactions(const int state = 0) const;
                    
  inline int Number_NodeEntity_Interactions(int state = 0) { 
    if((int)NodeEntityInteractions.size() <= state) return 0;
    return NodeEntityInteractions[state].size();
  }

  inline ContactNodeEntityInteraction** Get_NodeEntity_Interactions( int state = 0 ) { 
    if((int)NodeEntityInteractions.size() <= state) return NULL;
    if(NodeEntityInteractions[state].size() == 0) return NULL;
    return &(NodeEntityInteractions[state][0]); 
  };
       
  inline ContactNodeEntityInteraction*  Get_NodeEntity_Interaction(int interaction_number, int state = 0 ) { 
    if((int)NodeEntityInteractions.size() <= state) return NULL;
    if((int)NodeEntityInteractions[state].size() <= interaction_number) return NULL;
    return NodeEntityInteractions[state][interaction_number]; 
  };

  void Add_NodeEntity_Interaction( ContactNodeEntityInteraction* cnei,
				   int state = 0 );

  void Store_NodeEntity_Interaction( int number, 
				   ContactNodeEntityInteraction* cnei,
				   int state = 0 );
                                   
  int  Num_Tied_Interactions(int state = 0);
  int  Num_Tracked_Interactions(int state = 0);

  void Delete_NodeEntity_Interaction( ContactNodeEntityInteraction* cnei,
				      int state = 0 );
                                   
  void Delete_NodeEntity_Interactions( int state = 0 );

  void Compact_NodeEntity_Interactions( int state = 0 );

  void Print_NodeEntity_Interactions( int state = 0 );

  void Display_NodeEntity_Interactions( ContactParOStream&, int state = 0 );
 
  void Update_Interactions( );

  //These methods are used to update the tied/tracking status.
  inline void ComputeBoundingBoxForSearch( const int num_configs,
                                           const VariableHandle &CURRENT_POSITION,
                                           const VariableHandle &PREDICTED_POSITION,
                                           const VariableHandle &REMAINING_GAP,
                                           const int  auto_tol,
                                           const DataType box_inflation,
                                           ContactBoundingBox &box);
                                           
  inline void ComputeBoundingBoxForSearch( const int num_configs,
                                           const VariableHandle &CURRENT_POSITION,
                                           const VariableHandle &PREDICTED_POSITION,
                                           const VariableHandle &REMAINING_GAP,
                                           const int  auto_tol,
                                           const DataType box_inflation,
                                           const DataType box_expand,
                                           ContactBoundingBox &box);

 protected:
  using ContactTopologyEntity<DataType>::entity_key;
  using ContactTopologyEntity<DataType>::owner_proc_array_index;
  ContactSearch::ContactNode_Type node_type;
  ContactFixedSizeAllocator* allocators;
  int shell_node_base_id;

 private:
  int exodus_id;

  //
  //  List of faces along with the physical face of this node they belong to.  For 
  //  the physical face portion:
  //   -1   = unknown
  //   0->2 = physical face 0,1, or 2
  //   >3   = Extra face set that is part of no physcial face
  //
  std::vector<std::pair<ContactFace<DataType>*,int> > faces;
  PHYSICAL_TYPE physical_type;

  inline int NumberOfStates() const {return 2;};
  inline int MaxInteractionsPerNode() const {return 3;};



  DataType DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
  
  std::vector<ContactInteractionDLL<Real> > NodeNodeInteractions;
  std::vector< std::vector<ContactNodeEntityInteraction*> >NodeEntityInteractions;
};



template<typename DataType>
inline ContactFace<DataType>* ContactNode<DataType>::GetFace( int i ) {
  PRECONDITION( i>=0 && i<faces.size() );
  return( faces[i].first );
}
template<typename DataType>
inline int ContactNode<DataType>::GetFacePFIndex( int i ) {
  PRECONDITION( i>=0 && i<faces.size() );
  return( faces[i].second );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
template<typename DataType>
inline int ContactNode<DataType>::Size(int state)
{
  int size = ContactTopologyEntity<DataType>::Size(DataArray_Length()) + 6*sizeof(int);
  if (state>=0) {
    size += Size_Interactions(state);
  } else if (state==-1) { 
    for (int i=0; i<NumberOfStates(); i++) {
      size += Size_Interactions(i);
    }
  }
  return size;
}

template<typename DataType>
inline void ContactNode<DataType>::Pack(char* buffer, int state )
{
  char* buff = buffer;
  ContactTopologyEntity<DataType>::Pack( buff, DataArray_Length() );
  buff += ContactTopologyEntity<DataType>::Size(DataArray_Length());
  
  // NOTE: If the location of the 'state' variable in the buffer array
  // is hard coded into the ContactShellNode Pack() function so if anything 
  // with following changes here then need to make a change to the Pack()_
  // function there.
  int cnt      = 0;
  int* i_buf   = reinterpret_cast<int*>(buff);
  i_buf[cnt++] = entity_key;
  i_buf[cnt++] = exodus_id;
  i_buf[cnt++] = shell_node_base_id;
  i_buf[cnt++] = node_type;
  i_buf[cnt++] = physical_type;
  i_buf[cnt++] = state;
  
  if (state>=0) {
    buff += cnt*sizeof(int);
    Pack_Interactions(buff, state);
  } else if (state==-1) {
    buff += cnt*sizeof(int);
    for (int i=0; i<NumberOfStates(); i++) {
      Pack_Interactions(buff, i);
      buff += Size_Interactions(i);
    }
  }
}

template<typename DataType>
inline void ContactNode<DataType>::Unpack( char* buffer )
{
  char* buff = buffer;
  ContactTopologyEntity<DataType>::Unpack( buff, DataArray_Length() );
  buff += ContactTopologyEntity<DataType>::Size(DataArray_Length());
  
  int cnt                       = 0;
  int* i_buf                    = reinterpret_cast<int*>( buff );
  entity_key                    = i_buf[cnt++];
  exodus_id                     = i_buf[cnt++];
  shell_node_base_id            = i_buf[cnt++];
  node_type                     = (ContactSearch::ContactNode_Type) i_buf[cnt++];
  physical_type                 = (typename ContactNode<DataType>::PHYSICAL_TYPE) i_buf[cnt++];
  int state                     = i_buf[cnt++];
  faces.clear();
  
  if (state>=0) {
    buff += cnt*sizeof(int);
    Unpack_Interactions(buff, state);
  } else if (state==-1) {
    buff += cnt*sizeof(int);
    for (int i=0; i<NumberOfStates(); i++) {
      Unpack_Interactions(buff, i);
      buff += Size_Interactions(i);
    }
  }
}

template<typename DataType>
inline void ContactNode<DataType>::Copy( ContactNode* src, int state )
{
  ContactTopologyEntity<DataType>::Copy( src, DataArray_Length() );
  entity_key                    = src->entity_key;
  exodus_id                     = src->exodus_id;
  shell_node_base_id            = src->shell_node_base_id;
  node_type                     = src->node_type;
  physical_type                 = src->physical_type;
  faces.clear();
  
  if (state>=0) {
    Copy_Interactions(src, state);
  } else if (state==-1) {
    for (int i=0; i<NumberOfStates(); i++) {
      Copy_Interactions(src, i);
    }
  }
}

//for data xfer for search, don't send NODE_GHOST_GAP, PROJ_DIR, or LOFTING_VECTOR
#define DATA_LENGTH (DataArray_Length()-3*3)

template<typename DataType>
inline int ContactNode<DataType>::Size_ForDataUpdate()
{
  return( ContactTopologyEntity<DataType>::Size_ForDataUpdate(DATA_LENGTH) );
}

template<typename DataType>
inline void ContactNode<DataType>::Pack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity<DataType>::Pack_ForDataUpdate( buffer, DATA_LENGTH );
}

template<typename DataType>
inline void ContactNode<DataType>::Unpack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity<DataType>::Unpack_ForDataUpdate( buffer, DATA_LENGTH );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
template<typename DataType>
inline int ContactNode<DataType>::Size_ForSecondary(int state)
{
  int size = ContactTopologyEntity<DataType>::Size_ForSecondary(DATA_LENGTH) + 2*sizeof(int);
  if (state>=0) {
    size += Size_Interactions_ForSecondary(state);
  } else if (state==-1) {
    for (int i=0; i<NumberOfStates(); i++) {
      size += Size_Interactions_ForSecondary(i);
    }
  }
  return size;
}

template<typename DataType>
void ContactNode<DataType>::Pack_ForSecondary(char* buffer, int state )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  ContactTopologyEntity<DataType>::Pack_ForSecondary( buffer, DATA_LENGTH );
  
  // ContactTopologyEntity<DataType> packs in location 0 as BaseType
  // and here we packin some misc stuff in location 1.
  i_buf[1] = physical_type |  node_type<<4 | (state+2)<<12;
  
  // NOTE: If the location of the 'state' variable in the buffer array
  // is hard coded into the ContactShellNode Pack() function so if anything 
  // with following changes here then need to make a change to the Pack()_
  // function there.
  int cnt      = 0;
  char* buff   = buffer;
  buff        += ContactTopologyEntity<DataType>::Size_ForSecondary(DATA_LENGTH);
  i_buf        = reinterpret_cast<int*> (buff);
  i_buf[cnt++] = entity_key;
  i_buf[cnt++] = exodus_id;
  
  if (state>=0) {
    buff += cnt*sizeof(int);
    Pack_Interactions_ForSecondary(buff, state);
  } else if (state==-1) {
    buff += cnt*sizeof(int);
    for (int i=0; i<NumberOfStates(); i++) {
      Pack_Interactions_ForSecondary(buff, i);
      buff += Size_Interactions_ForSecondary(i);
    }
  }
}

template<typename DataType>
inline void ContactNode<DataType>::Unpack_ForSecondary( char* buffer )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  ContactTopologyEntity<DataType>::Unpack_ForSecondary( buffer, DATA_LENGTH );
  physical_type      = (typename ContactNode<DataType>::PHYSICAL_TYPE       )((i_buf[1]   ) & 0xf);
  node_type          = (ContactSearch::ContactNode_Type    )((i_buf[1]>>4) & 0xf);
  int state          = ((i_buf[1]>>12) & 0xf)-2;
  
  char* buff         = buffer;
  buff              += ContactTopologyEntity<DataType>::Size_ForSecondary(DATA_LENGTH);
  int cnt            = 0;
  i_buf              = reinterpret_cast<int*> ( buff );
  entity_key         = i_buf[cnt++];
  exodus_id          = i_buf[cnt++];
  
  faces.clear();
  
  if (state>=0) {
    buff += cnt*sizeof(int);
    Unpack_Interactions_ForSecondary(buff, state);
  } else if (state==-1) {
    buff += cnt*sizeof(int);
    for (int i=0; i<NumberOfStates(); i++) {
      Unpack_Interactions_ForSecondary(buff, i);
      buff += Size_Interactions_ForSecondary(i);
    }
  }
}

template<typename DataType>
inline void ContactNode<DataType>::Copy_ForSecondary( ContactNode* src, int state )
{
  ContactTopologyEntity<DataType>::Copy_ForSecondary( src, DATA_LENGTH );
  entity_key    = src->entity_key;
  exodus_id     = src->exodus_id;
  node_type     = src->node_type;
  physical_type = src->physical_type;
  
  faces.clear();
  
  if (state>=0) {
    Copy_Interactions_ForSecondary(src, state);
  } else if (state==-1) {
    for (int i=0; i<NumberOfStates(); i++) {
      Copy_Interactions_ForSecondary(src, i);
    }
  }
}

template<typename DataType>
inline void 
ContactNode<DataType>::ComputeBoundingBoxForSearch( const int num_configs,
                                          const VariableHandle &REMAINING_GAP,
                                          const VariableHandle &CURRENT_POSITION,
                                          const VariableHandle &PREDICTED_POSITION,
                                          const int  auto_tol,
                                          const DataType box_inflation,
                                          ContactBoundingBox &box)
{
  box.set_point(Variable(CURRENT_POSITION));
  if (num_configs>1) {
    box.add_point(Variable(PREDICTED_POSITION));
    box.add_tolerance(Variable(REMAINING_GAP));
  }
  if (auto_tol) {
    DataType box_tol[3] = {0.0, 0.0, 0.0};
    box.get_dimensions(box_tol);
    for (int i=0; i<3; ++i) {
      box_tol[i] *= box_inflation;
    }
    box.add_tolerance(box_tol);
  }
}

template<typename DataType>
inline void 
ContactNode<DataType>::ComputeBoundingBoxForSearch( const int num_configs,
                                          const VariableHandle &REMAINING_GAP,
                                          const VariableHandle &CURRENT_POSITION,
                                          const VariableHandle &PREDICTED_POSITION,
                                          const int  auto_tol,
                                          const DataType box_inflation,
                                          const DataType box_expand,
                                          ContactBoundingBox &box)
{
  box.set_point(Variable(CURRENT_POSITION));
  if (num_configs>1) {
    box.add_point(Variable(PREDICTED_POSITION));
    box.add_tolerance(Variable(REMAINING_GAP));
    box.add_tolerance(box_expand);
  }
  if (auto_tol) {
    DataType box_tol[3] = {0.0, 0.0, 0.0};
    box.add_tolerance(box_expand);
    box.get_dimensions(box_tol);
    for (int i=0; i<3; ++i) {
      box_tol[i] *= box_inflation;
    }
    box.add_tolerance(box_tol);
  }
}

#endif  // #ifdef ContactNode_h_
