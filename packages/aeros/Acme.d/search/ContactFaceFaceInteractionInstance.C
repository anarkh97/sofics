// $Id$

#include "ContactFaceFaceInteraction.C"

template
ContactFaceFaceInteraction<Real>::ContactFaceFaceInteraction( );

template
ContactFaceFaceInteraction<Real>::ContactFaceFaceInteraction( ContactFace<Real>* Sface,
							ContactFace<Real>* Mface,
							int Nedges, 
                                                        int* FaceEdge,
                                                        int* EdgeMaster,
                                                        Real* Sarea, 
                                                        Real* Marea,
                                                        Real (*Sarea_derivatives)[42],
                                                        Real (*Marea_derivatives)[42],
                                                        Real (*Sarea_second_derivatives)[42],
                                                        Real (*Marea_second_derivatives)[42] );

template
ContactFaceFaceInteraction<Real>::ContactFaceFaceInteraction( 
                                       ContactFaceFaceInteraction& ffi );

template
ContactFaceFaceInteraction<Real>* 
ContactFaceFaceInteraction<Real>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFace<Real>* Sface,
				     ContactFace<Real>* Mface,
				     int Nedges, 
                                     int* FaceEdge,
                                     int* EdgeMaster,
                                     Real* Sarea, Real* Marea,
                                     Real (*Sarea_derivatives)[42], Real (*Marea_derivatives)[42],
                                     Real (*Sarea_second_derivatives)[42], Real (*Marea_second_derivatives)[42] );

template
ContactFaceFaceInteraction<Real>* 
ContactFaceFaceInteraction<Real>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc );

template
ContactFaceFaceInteraction<Real>* 
ContactFaceFaceInteraction<Real>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFaceFaceInteraction& cffi );

template
void ContactFaceFaceInteraction_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactFaceFaceInteraction<Real>::~ContactFaceFaceInteraction();

template
int ContactFaceFaceInteraction<Real>::Size();

template
void ContactFaceFaceInteraction<Real>::Pack( char* buffer );

template
void ContactFaceFaceInteraction<Real>::Unpack( char* buffer );

template
void ContactFaceFaceInteraction<Real>::Copy( ContactFaceFaceInteraction<Real>* src );

template
void ContactFaceFaceInteraction<Real>::Connect_SlaveFace( ContactTopologyEntityList& hash_table );

template
void ContactFaceFaceInteraction<Real>::Connect_MasterFace( ContactTopologyEntityList& hash_table );

template
void ContactFaceFaceInteraction<Real>::Connect_SlaveFace( ContactTopologyEntityHash& hash_table );

template
void ContactFaceFaceInteraction<Real>::Connect_MasterFace( ContactTopologyEntityHash& hash_table );

template
void ContactFaceFaceInteraction<Real>::Connect_SlaveFace( ContactTopology* topology );

template
void ContactFaceFaceInteraction<Real>::Connect_MasterFace( ContactTopology* topology );

template
void ContactFaceFaceInteraction<Real>::Connect_SlaveFace( ContactFace<Real>* Face );

template
void ContactFaceFaceInteraction<Real>::Connect_MasterFace( ContactFace<Real>* Face );

template
int ContactFaceFaceInteraction<Real>::Data_Size();

template
int ContactFaceFaceInteraction<Real>::Restart_Size();

template
void ContactFaceFaceInteraction<Real>::Restart_Pack( Real* buffer );

template
void ContactFaceFaceInteraction<Real>::Restart_Unpack( Real* buffer );

template
int ContactFaceFaceInteraction<Real>::Set_SlaveFaceEntityData(); 

template
int ContactFaceFaceInteraction<Real>::Set_MasterFaceEntityData();

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactFaceFaceInteraction<ActiveScalar>::ContactFaceFaceInteraction( );

template
ContactFaceFaceInteraction<ActiveScalar>::ContactFaceFaceInteraction( ContactFace<ActiveScalar>* Sface,
							ContactFace<ActiveScalar>* Mface,
							int Nedges, 
                                                        int* FaceEdge,
                                                        int* EdgeMaster,
                                                        ActiveScalar* Sarea, 
                                                        ActiveScalar* Marea,
                                                        ActiveScalar (*Sarea_derivatives)[42],
                                                        ActiveScalar (*Marea_derivatives)[42],
                                                        ActiveScalar (*Sarea_second_derivatives)[42],
                                                        ActiveScalar (*Marea_second_derivatives)[42] );

template
ContactFaceFaceInteraction<ActiveScalar>::ContactFaceFaceInteraction( 
                                       ContactFaceFaceInteraction& ffi );

template
ContactFaceFaceInteraction<ActiveScalar>* 
ContactFaceFaceInteraction<ActiveScalar>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFace<ActiveScalar>* Sface,
				     ContactFace<ActiveScalar>* Mface,
				     int Nedges, 
                                     int* FaceEdge,
                                     int* EdgeMaster,
                                     ActiveScalar* Sarea, ActiveScalar* Marea,
                                     ActiveScalar (*Sarea_derivatives)[42], ActiveScalar (*Marea_derivatives)[42],
                                     ActiveScalar (*Sarea_second_derivatives)[42], ActiveScalar (*Marea_second_derivatives)[42] );

template
ContactFaceFaceInteraction<ActiveScalar>* 
ContactFaceFaceInteraction<ActiveScalar>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc );

template
ContactFaceFaceInteraction<ActiveScalar>* 
ContactFaceFaceInteraction<ActiveScalar>::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFaceFaceInteraction& cffi );

template
void ContactFaceFaceInteraction_SizeAllocator<ActiveScalar>(ContactFixedSizeAllocator& alloc);

template
ContactFaceFaceInteraction<ActiveScalar>::~ContactFaceFaceInteraction();

template
int ContactFaceFaceInteraction<ActiveScalar>::Size();

template
void ContactFaceFaceInteraction<ActiveScalar>::Pack( char* buffer );

template
void ContactFaceFaceInteraction<ActiveScalar>::Unpack( char* buffer );

template
void ContactFaceFaceInteraction<ActiveScalar>::Copy( ContactFaceFaceInteraction<ActiveScalar>* src );
#endif
