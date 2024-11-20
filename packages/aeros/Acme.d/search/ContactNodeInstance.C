// $Id$

#include "ContactNode.C"

template
ContactNode<Real>::ContactNode(ContactFixedSizeAllocator* alloc,
                         ContactSearch::ContactNode_Type Type, 
                         int Block_Index, 
                         int Host_Index_in_Block,
			 ContactType cetype ); 

template
ContactNode<Real>* ContactNode<Real>::new_ContactNode( ContactFixedSizeAllocator* alloc,
					   ContactSearch::ContactNode_Type Type,
					   int Block_Index, 
					   int Host_Index_in_Block,
					   ContactType cetype );

template
void ContactNode_SizeAllocator<Real>( ContactFixedSizeAllocator& alloc);

template
ContactNode<Real>::~ContactNode();

template
void ContactNode<Real>::Delete_Face_Connections( );

template
void ContactNode<Real>::Connect_Face(ContactFace<Real>* Face );

template
int ContactNode<Real>::Size_Interactions(int state);

template
void ContactNode<Real>::Pack_Interactions( char* Buffer, int state );

template
void ContactNode<Real>::Unpack_Interactions( char* buffer, int state );

template
void ContactNode<Real>::Copy_Interactions( ContactNode<Real>* src, int state );

template
int ContactNode<Real>::Size_Interactions_ForSecondary(int state);

template
void ContactNode<Real>::Pack_Interactions_ForSecondary( char* Buffer, int state );

template
void ContactNode<Real>::Unpack_Interactions_ForSecondary( char* buffer, int state );

template
void ContactNode<Real>::Copy_Interactions_ForSecondary( ContactNode<Real>* src, int state );

template
ContactNodeNodeInteraction* 
ContactNode<Real>::Get_NodeNode_Interaction(int interaction_number, int state );

template
void 
ContactNode<Real>::Add_NodeNode_Interaction( 
                ContactNodeNodeInteraction* nni, int state );

template
void 
ContactNode<Real>::Delete_NodeNode_Interaction( 
                ContactNodeNodeInteraction* nni, int state );

template
void 
ContactNode<Real>::Delete_NodeNode_Interactions( int state );

template
void 
ContactNode<Real>::Display_NodeNode_Interactions( ContactParOStream& postream,
					    int state );

template
void
ContactNode<Real>::Add_NodeEntity_Interaction(ContactNodeEntityInteraction* nei,
				        int state );

template
void
ContactNode<Real>::Store_NodeEntity_Interaction(int interaction_number,
					  ContactNodeEntityInteraction* nei,
					  int state );
  
template
void 
ContactNode<Real>::Delete_NodeEntity_Interaction(ContactNodeEntityInteraction* nei, 
                                           int state );

template
void 
ContactNode<Real>::Delete_NodeEntity_Interactions( int state );
  
template
void 
ContactNode<Real>::Display_NodeEntity_Interactions( ContactParOStream& postream,
					      int state );

template
void 
ContactNode<Real>::Print_NodeEntity_Interactions( int state );

template
void
ContactNode<Real>::Update_Interactions( );

template
int ContactNode<Real>::Number_NodeFace_Interactions(const int state) const;

template
int ContactNode<Real>::Number_NodeSurface_Interactions(const int state) const;

template
int ContactNode<Real>::Get_Owning_Entity();

template
int ContactNode<Real>::Num_Tied_Interactions(int istate);

template
int ContactNode<Real>::Num_Tracked_Interactions(int istate);

template
int ContactNode<Real>::GetFacePFEntityKey(int physical_face_num);

template
bool ContactNode<Real>::ConnectedToFace(const ContactHostGlobalID &id);

template
bool sort_face_pair_by_id<Real>(const pair<ContactFace<Real>*,int> &d1, const pair<ContactFace<Real>*,int> &d2);

template
void ContactNode<Real>::SortConnectedFaces();


#if (MAX_FFI_DERIVATIVES > 0)
template
ContactNode<ActiveScalar>::ContactNode(ContactFixedSizeAllocator* alloc,
                         ContactSearch::ContactNode_Type Type,
                         int Block_Index,
                         int Host_Index_in_Block,
                         ContactType cetype );

template
ContactNode<ActiveScalar>* ContactNode<ActiveScalar>::new_ContactNode( ContactFixedSizeAllocator* alloc,
                                           ContactSearch::ContactNode_Type Type,
                                           int Block_Index,
                                           int Host_Index_in_Block,
                                           ContactType cetype );

template
void ContactNode_SizeAllocator<ActiveScalar>( ContactFixedSizeAllocator& alloc);

template
ContactNode<ActiveScalar>::~ContactNode();

template
void ContactNode<ActiveScalar>::Delete_Face_Connections( );

template
void ContactNode<ActiveScalar>::Connect_Face(ContactFace<ActiveScalar>* Face );

template
int ContactNode<ActiveScalar>::Size_Interactions(int state);

template<>
void ContactNode<ActiveScalar>::Unpack_Interactions( char* buffer, int state )
{
  std::cerr << "ContactNode<ActiveScalar>::Unpack_Interactions is not implemented\n";
}

template<>
void ContactNode<ActiveScalar>::Copy_Interactions( ContactNode* src, int state )
{
  std::cerr << "ContactNode<ActiveScalar>::Copy_Interactions is not implemented\n";
}

template<>
void ContactNode<ActiveScalar>::Unpack_Interactions_ForSecondary( char* buffer, int state )
{
  std::cerr << "ContactNode<ActiveScalar>::Unpack_Interactions_ForSecondary is not implemented\n";
}

template<>
void ContactNode<ActiveScalar>::Copy_Interactions_ForSecondary( ContactNode* src, int state )
{
  std::cerr << "ContactNode<ActiveScalar>::Copy_Interactions_ForSecondary is not implemented\n";
}
#endif
