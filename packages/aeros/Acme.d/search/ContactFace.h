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


#ifndef ContactFace_h_
#define ContactFace_h_

#include "ContactTopologyEntity.h"
#include "contact_assert.h"
#include "ContactEdge.h"
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactDoublyLinkedList.h"
#include "ContactNode.h"

template<typename DataType> class ContactNode;
template<typename DataType> class ContactEdge;
template<typename DataType> class ContactFaceFaceInteraction;
class ContactFaceCoverageInteraction;

template<typename DataType>
class ContactFace : public ContactTopologyEntity<DataType> {

 public:

  using ContactTopologyEntity<DataType>::DataArray_Buffer;

  enum Ctrcl_Index{MSPARAM = 0,  ICPOINTX = 1,  ICPOINTY = 2,  ICPOINTZ  = 3,
                   IPENMAG = 4,  IPUSHX   = 5,  IPUSHY   = 6,  IPUSHZ    = 7,
                   INORMX  = 8,  INORMY   = 9,  INORMZ   = 10, ILOCATION = 11,
                   ICTIMC  = 12, LENGTH  = 13};
  
  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  FACE_SCALAR_VAR
#define FACE_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  FACE_VECTOR_VAR
#define FACE_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactFace( ContactFixedSizeAllocator*,
               ContactSearch::ContactFace_Type, 
               int Block_Index,
	       int Host_Index_in_Block, 
               int key,
               ContactNode<DataType> **node_list_,
               ContactEdge<DataType> **edge_list_,
               typename ContactTopologyEntity<DataType>::connection_data *node_info_list_,
               typename ContactTopologyEntity<DataType>::connection_data *edge_info_list_);

  virtual ~ContactFace();

  static inline int Nodes_Per_Face(ContactSearch::ContactFace_Type check_face_type) {
    PRECONDITION(array_init);
    return NODES_PER_FACE[check_face_type];
  }



  inline int Nodes_Per_Face() {
    PRECONDITION(array_init);
    return NODES_PER_FACE[face_type];
  };

  inline int Edges_Per_Face() {
    PRECONDITION(array_init);
    return EDGES_PER_FACE[face_type];
  };

  inline ContactNode<DataType>**     Nodes()    {return node_list;};
  inline typename ContactTopologyEntity<DataType>::connection_data*  NodeInfo() {return node_info_list;};
  inline ContactEdge<DataType>**     Edges()    {return edge_list;};
  inline typename ContactTopologyEntity<DataType>::connection_data*  EdgeInfo() {return edge_info_list;};


  virtual ContactSearch::ContactEdge_Type Edge_Type() = 0;
  inline ContactNode<DataType>* Node( const int i );
  static void Initialize_Lookup_Arrays();


  // This function get the two nodes that terminate edge i
  virtual void Get_Edge_Nodes( int, ContactNode<DataType>**) = 0;
  virtual int  Get_Edge_Number( ContactNode<DataType>** ) = 0;
  virtual int  Get_Edge_Number( DataType* ) = 0;
  virtual void Compute_Normal(VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Partial_Face_Normal(VariableHandle, DataType (*)[3] )
               { std::cerr << "ContactFace::Compute_Partial_Face_Normal is not implemented\n"; }
  virtual void Compute_Second_Partial_Face_Normal(VariableHandle, DataType (*)[3] )
               { std::cerr << "ContactFace::Compute_Second_Partial_Face_Normal is not implemented\n"; }
  virtual void Compute_Normal(VariableHandle, DataType*, DataType* ) = 0;
  virtual void Compute_Normal(DataType**, DataType*, DataType* ) = 0;
  virtual void Compute_CharacteristicLength(VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Centroid( VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Edge_Normal( VariableHandle, VariableHandle,
				    int , DataType*) = 0;
  virtual void Compute_Local_Coordinates( DataType, VariableHandle, 
                                          VariableHandle, VariableHandle, 
		                          DataType*, DataType* ) = 0;

  virtual void Compute_Local_Coords(DataType node_positions[MAX_NODES_PER_FACE][3], 
				    DataType global_coords[3],
			            DataType local_coords[3]) = 0;


  virtual void Compute_Local_Coordinates( VariableHandle, DataType*, DataType* ) = 0;
  virtual void Compute_Partial_Local_Coordinates_1(VariableHandle, DataType*, DataType (*)[2] )
               { std::cerr << "ContactFace::Compute_Partial_Local_Coordinates_1 is not implemented\n"; }
  virtual void Compute_Partial_Local_Coordinates_2(VariableHandle, DataType*, DataType[2], DataType[2], DataType[2] )
               { std::cerr << "ContactFace::Compute_Partial_Local_Coordinates_2 is not implemented\n"; }
  virtual void Compute_Second_Partial_Local_Coordinates_1(VariableHandle, DataType*, DataType (*)[2] )
               { std::cerr << "ContactFace::Compute_Second_Partial_Local_Coordinates_1 is not implemented\n"; }
  virtual void Compute_Second_Partial_Local_Coordinates_2(VariableHandle, DataType*, DataType[2], DataType[2], DataType[2],
                                                          DataType[2], DataType[2], DataType[2] )
               { std::cerr << "ContactFace::Compute_Second_Partial_Local_Coordinates_2 is not implemented\n"; }
  virtual void Compute_Second_Partial_Local_Coordinates_12(VariableHandle, DataType*, DataType (*)[2], DataType (*)[2], DataType (*)[2] )
               { std::cerr << "ContactFace::Compute_Second_Partial_Local_Coordinates_12 is not implemented\n"; }
  virtual void Compute_Global_Coordinates( VariableHandle, DataType*, DataType* ) = 0;
  virtual void Evaluate_Shape_Functions( DataType* local_coord, 
					 DataType* shape_fnc) = 0;
  virtual bool Is_Inside_Face( DataType* Local_Coordinates ) = 0;
  virtual ContactFace* Neighbor( DataType* Local_Coordinates ) = 0;

  // this is provided for the curvature modification 
  virtual void Get_Close_Edges( DataType*, int&, int&, int& ) = 0;
  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(DataType));};

  virtual void FacetDecomposition(int& nfacets, 
                                  DataType* coords0       , DataType* normal0       , 
				  VariableHandle POSITION0,
                                  DataType* coords1 = NULL, DataType* normal1 = NULL, 
				  VariableHandle POSITION1 = 0,
                                  DataType* coords2 = NULL, DataType* normal2 = NULL, 
				  VariableHandle POSITION2 = 0) = 0;
  virtual void FacetStaticRestriction(int, DataType*, DataType*, DataType*, DataType*) = 0;
  virtual void FacetDynamicRestriction(int, DataType*, DataType*) = 0;
  virtual int  FaceEdge_Intersection(VariableHandle, ContactEdge<DataType>*, DataType*) = 0;
  virtual bool IsPlanar(VariableHandle) = 0;
 
#ifndef CONTACT_NO_MPI
  DataType MaxSize(VariableHandle POSITION);
#endif

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline ContactEdge<DataType>* Edge( const int i );
  void ConnectNode( const int, ContactNode<DataType>* );
  void ConnectEdge( const int, ContactEdge<DataType>* );

  inline ContactEdge<DataType>* Clockwise_Edge( ContactEdge<DataType>* );
  inline void Clockwise_EdgeNode( int, ContactNode<DataType>** );

  inline ContactSearch::ContactFace_Type FaceType() {return face_type;};

  // Packing/Unpacking Functions
  inline int  Size(int flag=0);
  inline void Pack( char*, int flag=0 );
  inline void Unpack( char* );
  inline void Copy( ContactFace* face, int neighbors=0);
  
  inline int  Size_ForSecondary(int include_neighbors=0, int include_edgeinfo=0);
  inline void Pack_ForSecondary( char*, int include_neighbors=0, int include_edgeinfo=0 );
  inline void Unpack_ForSecondary( char* );
  inline void Copy_ForSecondary( ContactFace* face, int include_neighbors=0, int include_edgeinfo=0);
  
  inline int  Size_ForDataUpdate();
  inline void Pack_ForDataUpdate( char* );
  inline void Unpack_ForDataUpdate( char* );
  
  int  Size_Interactions( int state=0 );
  void Pack_Interactions( char*, int state=0 );
  void Unpack_Interactions( char*, int state=0 );
  void Copy_Interactions( ContactFace*, int state=0 );
  
  int  Size_Interactions_ForSecondary( int state=0 );
  void Pack_Interactions_ForSecondary( char*, int state=0 ); 
  void Unpack_Interactions_ForSecondary( char*, int state=0 );
  void Copy_Interactions_ForSecondary( ContactFace*, int state=0 );
  
  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( DataType* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(DataType) );
  }
  inline void Restart_Unpack( DataType* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(DataType) );
  }

  virtual void Smooth_Normal( VariableHandle, VariableHandle, VariableHandle,
			      VariableHandle, 
			      ContactSearch::Smoothing_Resolution,
			      DataType, DataType*, DataType*, DataType ) = 0;
  
  virtual void Compute_Node_Areas( VariableHandle, VariableHandle, DataType* ) = 0;

  inline int Number_Interactions(int state = 0) {
    int num_ffi = 0;
    if((int)FaceFaceInteractions.size() > state) num_ffi = FaceFaceInteractions[state].NumEntities();
    int num_fci = 0;
    if((int)FaceCoverageInteractions.size() > state) num_fci = FaceCoverageInteractions[state].NumEntities();
    return num_ffi + num_fci;
  };
  
  inline int Number_FaceFace_Interactions (int state = 0) { 
    if((int)FaceFaceInteractions.size() <= state) return 0;
    return FaceFaceInteractions[state].NumEntities(); 
  };
                                         
  inline ContactInteractionDLL<DataType>* Get_FaceFace_Interactions( int state = 0 ) { 
    if((int)FaceFaceInteractions.size() <= state) return NULL;
    return &(FaceFaceInteractions[state]); 
  };
      
  ContactFaceFaceInteraction<DataType>* Get_FaceFace_Interaction(int interaction_number,
						       int state = 0 );
                                         
  void Store_FaceFace_Interaction( ContactFaceFaceInteraction<DataType>*, 
                                   int state = 0 );
                                         
  void Delete_FaceFace_Interaction( ContactFaceFaceInteraction<DataType>*, 
                                    int state = 0 );

  void Clear_FaceFace_Interactions() { FaceFaceInteractions.clear(); }
                                   
  void Display_FaceFace_Interactions( ContactParOStream&, int state = 0 );
  
  inline int Number_FaceCoverage_Interactions (int state = 0) { 
    if((int)FaceCoverageInteractions.size() <= state) return 0;
    return FaceCoverageInteractions[state].NumEntities(); 
  };
                                         
  inline ContactInteractionDLL<DataType>* Get_FaceCoverage_Interactions( int state = 0 ) { 
    if((int)FaceCoverageInteractions.size() <= state) return NULL;
    return &(FaceCoverageInteractions[state]); 
  };
                                         
  void Store_FaceCoverage_Interaction( ContactFaceCoverageInteraction*, 
                                       int state = 0 );
                                       
  void Display_FaceCoverage_Interactions( ContactParOStream&, int state = 0 );
  
  void Update_Interactions( );
  
  void SetNeighborFacesInfo( );
  
  int NumberOfNeighbors() { return number_of_neighbors; };
  typename ContactTopologyEntity<DataType>::connection_data* NeighborInfo() { return neighbor_face_info; };

  void SetEdgeCurvature(VariableHandle);

  void SetEdgeCurvature(VariableHandle &var, ContactEdge<DataType> *edge);


  DataType GetEdgeCurvature(int);
  
  void GetEdgeInfo(ContactNode<DataType>* node, ContactNode<DataType>** edge_nodes, 
                   int* edge_nums);

  void SetEdgeSmoothedNormal(VariableHandle);
  void GetEdgeSmoothedNormal(int, DataType*);
  
  void ComputeBoundingBoxForSearch(const int num_configs,
                                   const VariableHandle &NODE_COORD_START,
                                   const VariableHandle &NODE_COORD_END,
                                   const int  auto_tol,
                                   const DataType box_inflation,
                                   const DataType user_tol,
                                   ContactBoundingBox &box_c,
                                   ContactBoundingBox &box_p,
                                   ContactBoundingBox &box_s);
  
  void ComputeBoundingBoxForSearch(const int num_configs,
                                   const VariableHandle &NODE_COORD_START,
                                   const VariableHandle &NODE_COORD_END,
                                   const int  auto_tol,
                                   const DataType box_inflation,
                                   const DataType* max_node_motion,
                                   const DataType max_remaining_gap_mag,
                                   const DataType user_search_tol,
                                   ContactBoundingBox &box_c,
                                   ContactBoundingBox &box_p,
                                   ContactBoundingBox &box_s);

  

 protected:
  using ContactTopologyEntity<DataType>::entity_key;
  using ContactTopologyEntity<DataType>::block_id;
  int number_of_neighbors;
  ContactSearch::ContactFace_Type face_type;
  ContactFixedSizeAllocator* allocators;

 private:
  //
  //  Arrays for fast lookup of face info based on face type
  //
  static bool array_init;
  static int  NODES_PER_FACE[ContactSearch::NFACE_TYPES];
  static int  EDGES_PER_FACE[ContactSearch::NFACE_TYPES];

  // The edges and node arrays are actually owned in the derived class
  // but we hold a pointer to them to provide access through the base class.
  
  ContactNode<DataType> **node_list;
  ContactEdge<DataType> **edge_list;
  typename ContactTopologyEntity<DataType>::connection_data *node_info_list;   
  typename ContactTopologyEntity<DataType>::connection_data *edge_info_list;

  typename ContactTopologyEntity<DataType>::connection_data* neighbor_face_info;
  
//  int entity_key;

  DataType DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
  
  inline int NumberOfStates() const {return 2;};

  std::vector<ContactInteractionDLL<DataType> > FaceFaceInteractions;
  std::vector<ContactInteractionDLL<Real> > FaceCoverageInteractions;
};

template<typename DataType>
inline ContactNode<DataType>* ContactFace<DataType>::Node( const int i ) { 
  PRECONDITION( i>=0 && i<Nodes_Per_Face() );
  return( Nodes()[i] );
}

template<typename DataType>
inline ContactEdge<DataType>* ContactFace<DataType>::Edge( const int i ) { 
  PRECONDITION( i>=0 && i<Edges_Per_Face() );
  return( Edges()[i] );
}

template<typename DataType>
inline ContactEdge<DataType>* ContactFace<DataType>::Clockwise_Edge( ContactEdge<DataType>* edge ) {
  int i;
  // Find this edge in the edge list
  for( i=0 ; i<Edges_Per_Face() ; ++i){
    if( Edges()[i] == edge ){
      if( i == Edges_Per_Face()-1 )
	return Edges()[0];
      else 
	return Edges()[i+1];
    }
  }
  // We didn't find this edge which is an error
  POSTCONDITION( 0 );
  return NULL;
}

template<typename DataType>
inline void ContactFace<DataType>::Clockwise_EdgeNode( int edge, ContactNode<DataType>** edge_nodes)
{
  int n1 = edge==Edges_Per_Face()-1?0:edge+1;
  int n2 =   n1==Edges_Per_Face()-1?0:n1+1;
  edge_nodes[0] = node_list[n1];
  edge_nodes[1] = node_list[n2];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
template<typename DataType>
inline int ContactFace<DataType>::Size(int include_neighbors)
{
  if (include_neighbors) {
    return( ContactTopologyEntity<DataType>::Size(DataArray_Length()) + sizeof(int) +
	    Nodes_Per_Face()*sizeof(typename ContactTopologyEntity<DataType>::connection_data) +
	    Edges_Per_Face()*sizeof(typename ContactTopologyEntity<DataType>::connection_data) +
            Edges_Per_Face()*sizeof(typename ContactTopologyEntity<DataType>::connection_data));
  } else {
    return( ContactTopologyEntity<DataType>::Size(DataArray_Length()) + sizeof(int) +
            Nodes_Per_Face()*sizeof(typename ContactTopologyEntity<DataType>::connection_data) +
	    Edges_Per_Face()*sizeof(typename ContactTopologyEntity<DataType>::connection_data));
  }
}

template<typename DataType>
inline void ContactFace<DataType>::Pack( char* buffer, int include_neighbors )
{
  int* i_buf = reinterpret_cast<int*>(buffer);
  // ContactTopologyEntity<DataType> packs in location 0 as ContactFace and here we pack
  // in the derived type in location 1.
  i_buf[1] = face_type;
  ContactTopologyEntity<DataType>::Pack( buffer, DataArray_Length() );
  // Add the entity data for the nodes and edges
  i_buf = reinterpret_cast<int*>(buffer+ContactTopologyEntity<DataType>::Size(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    cnt += this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(Node(i)), &i_buf[cnt]);
  }
  for( int i=0 ; i<Edges_Per_Face() ; ++i ){
    cnt += this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(Edge(i)), &i_buf[cnt]);
  }
  if (include_neighbors) {
    i_buf[cnt++] = number_of_neighbors;
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      cnt += this->PackConnection(&neighbor_face_info[i], &i_buf[cnt]);
    }
  } else {
    i_buf[cnt] = 0;
  }
}

template<typename DataType>
inline void ContactFace<DataType>::Unpack( char* buffer )
{
  ContactTopologyEntity<DataType>::Unpack( buffer, DataArray_Length() );
  entity_key = block_id;

  PRECONDITION( ((int*)buffer)[1] == face_type );
  // Store off the entity data for the nodes and edge
  int* i_buf = reinterpret_cast<int*> ( buffer + ContactTopologyEntity<DataType>::Size(DataArray_Length()) );
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    cnt += this->UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
  for( int i=0 ; i<Edges_Per_Face() ; ++i ){
    cnt += this->UnPackConnection(&EdgeInfo()[i], &i_buf[cnt]);
  }
  number_of_neighbors = i_buf[cnt++];
  if (number_of_neighbors>0) {
    neighbor_face_info = new typename ContactTopologyEntity<DataType>::connection_data[Edges_Per_Face()];
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      cnt += this->UnPackConnection(&neighbor_face_info[i], &i_buf[cnt]);
    }
  }
}

template<typename DataType>
inline void ContactFace<DataType>::Copy( ContactFace* src, int include_neighbors )
{
  ContactTopologyEntity<DataType>::Copy( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(src->Node(i)), &NodeInfo()[i]);
  }
  for( int i=0 ; i<Edges_Per_Face() ; ++i ){
    this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(src->Edge(i)), &EdgeInfo()[i]);
  }
  if (include_neighbors) {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info  = NULL;
    number_of_neighbors = src->number_of_neighbors;
    if (number_of_neighbors>0) {
      neighbor_face_info = new typename ContactTopologyEntity<DataType>::connection_data[Edges_Per_Face()];
      for (int i=0; i<Edges_Per_Face(); ++i) {
        neighbor_face_info[i] = src->neighbor_face_info[i];
      }
    }
  } else {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info  = NULL;
    number_of_neighbors = 0;
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
template<typename DataType>
inline int ContactFace<DataType>::Size_ForSecondary(int include_neighbors, 
                                          int include_edgeinfo)
{
  int cnt = sizeof(int) + Nodes_Per_Face()*sizeof(typename ContactTopologyEntity<DataType>::connection_data);
  if (include_edgeinfo) {
    cnt += ContactTopologyEntity<DataType>::Size_ForSecondary(DataArray_Length());
    cnt += Edges_Per_Face()*sizeof(typename ContactTopologyEntity<DataType>::connection_data);
  } else {
    cnt += ContactTopologyEntity<DataType>::Size_ForSecondary(DataArray_Length());
  }
  if (include_neighbors) {
    cnt += number_of_neighbors*(sizeof(typename ContactTopologyEntity<DataType>::connection_data)+sizeof(int));
  }
  return cnt;
}

template<typename DataType>
inline void ContactFace<DataType>::Pack_ForSecondary( char* buffer, 
                                            int include_neighbors, 
                                            int include_edgeinfo )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  // ContactTopologyEntity<DataType> packs in location 0 as ContactFace and here we pack
  // in the derived type in location 1.
  i_buf[1] = face_type;
  
  if (include_edgeinfo) {
    ContactTopologyEntity<DataType>::Pack_ForSecondary( buffer, DataArray_Length() );
  } else {
    ContactTopologyEntity<DataType>::Pack_ForSecondary( buffer, &DataArray[Edge0_Curvature], DataArray_Length() );
  }
  i_buf[ContactTopologyEntity<DataType>::SEC_OWNER] |= include_edgeinfo<<24;
  i_buf = reinterpret_cast<int*>(buffer+ContactTopologyEntity<DataType>::Size_ForSecondary(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    cnt += this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(Node(i)), &i_buf[cnt]);
  }
  if (include_edgeinfo) {
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      cnt += this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(Edge(i)), &i_buf[cnt]);
    }
  }
  if (include_neighbors) {
    int nn=0;
    i_buf[cnt++] = number_of_neighbors;
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      if (neighbor_face_info[i].owner>=0) {
        ++nn;
        i_buf[cnt++] = i;
        cnt += this->PackConnection(&neighbor_face_info[i], &i_buf[cnt]);
      }
    }
    POSTCONDITION(nn==number_of_neighbors);
  } else {
    i_buf[cnt++] = 0;
  }
}

template<typename DataType>
inline void ContactFace<DataType>::Unpack_ForSecondary( char* buffer )
{
  int* i_buf = reinterpret_cast<int*>( buffer );
  PRECONDITION( i_buf[1] == face_type );
  int include_edgeinfo = i_buf[ContactTopologyEntity<DataType>::SEC_OWNER]>>24;
  i_buf[ContactTopologyEntity<DataType>::SEC_OWNER] &= 0xFFFFFF;
  if (include_edgeinfo) {
    ContactTopologyEntity<DataType>::Unpack_ForSecondary( buffer, DataArray_Length() );
  } else {
    ContactTopologyEntity<DataType>::Unpack_ForSecondary( buffer, &DataArray[Edge0_Curvature], DataArray_Length() );
  }
  entity_key = block_id;

  // Store off the entity data for the nodes and edge
  i_buf = reinterpret_cast<int*>( buffer + ContactTopologyEntity<DataType>::Size_ForSecondary(DataArray_Length()) );
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i){
    cnt += this->UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
  if (include_edgeinfo) {
    for( int i=0 ; i<Edges_Per_Face() ; ++i){
      cnt += this->UnPackConnection(&EdgeInfo()[i], &i_buf[cnt]);
    }
  }
  neighbor_face_info = new typename ContactTopologyEntity<DataType>::connection_data[Edges_Per_Face()];
  for (int i=0; i<Edges_Per_Face(); ++i) {
    neighbor_face_info[i].owner = -1;;
  }
  number_of_neighbors = i_buf[cnt++];
  if (number_of_neighbors>0) {
    int nn = 0;
    while (nn<number_of_neighbors) {
      int ii = i_buf[cnt++];
      cnt   += this->UnPackConnection(&neighbor_face_info[ii], &i_buf[cnt]);
      ++nn;
    }
  }
}

template<typename DataType>
inline void ContactFace<DataType>::Copy_ForSecondary( ContactFace* src, 
                                            int include_neighbors, 
                                            int include_edgeinfo )
{
  ContactTopologyEntity<DataType>::Copy_ForSecondary( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(src->Node(i)), &NodeInfo()[i]);
  }
  if (include_edgeinfo) {
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      this->PackConnection(reinterpret_cast<ContactTopologyEntity<DataType>*>(src->Edge(i)), &EdgeInfo()[i]);
    }
  }
  if (include_neighbors) {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info  = NULL;
    number_of_neighbors = src->number_of_neighbors;
    if (number_of_neighbors>0) {
      neighbor_face_info = new typename ContactTopologyEntity<DataType>::connection_data[Edges_Per_Face()];
      for (int i=0; i<Edges_Per_Face(); ++i) {
        neighbor_face_info[i] = src->neighbor_face_info[i];
      }
    }
  } else {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info = new typename ContactTopologyEntity<DataType>::connection_data[Edges_Per_Face()];
    for (int i=0; i<Edges_Per_Face(); ++i) {
      neighbor_face_info[i].owner = -1;;
    }
    number_of_neighbors = 0;
  }
}

template<typename DataType>
inline int ContactFace<DataType>::Size_ForDataUpdate()
{
  return( ContactTopologyEntity<DataType>::Size_ForDataUpdate(DataArray_Length()) );
}

template<typename DataType>
inline void ContactFace<DataType>::Pack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity<DataType>::Pack_ForDataUpdate( buffer, DataArray_Length() );
}

template<typename DataType>
inline void ContactFace<DataType>::Unpack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity<DataType>::Unpack_ForDataUpdate( buffer, DataArray_Length() );
}

#endif // #ifdef ContactFace_h_

