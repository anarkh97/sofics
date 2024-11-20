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


#ifndef ContactFaceCoverageInteraction_h_
#define ContactFaceCoverageInteraction_h_

#include "Contact_Defines.h"
#include "ContactInteractionEntity.h"
#include "ContactFace.h"
#include "ContactSearch.h"

class CString;
class ContactTopologyEntityList;
class ContactHostGlobalID;
class ContactFixedSizeAllocator;

struct ContactFaceCoverageVertex {
  Real slave_x;
  Real slave_y;
  struct ContactFaceCoverageVertex* next;
};

class ContactFaceCoverageInteraction : public ContactInteractionEntity<Real> {
  
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

  ContactFaceCoverageInteraction();
  ContactFaceCoverageInteraction( ContactFace<Real>* );
  ContactFaceCoverageInteraction( ContactFaceCoverageInteraction& );
  static ContactFaceCoverageInteraction* new_ContactFaceCoverageInteraction(
            ContactFixedSizeAllocator&, ContactFace<Real>* );
  static ContactFaceCoverageInteraction* new_ContactFaceCoverageInteraction(
	     ContactFixedSizeAllocator& );
  static ContactFaceCoverageInteraction* new_ContactFaceCoverageInteraction(
	     ContactFixedSizeAllocator&, ContactFaceCoverageInteraction& );
  ~ContactFaceCoverageInteraction();

  inline ContactFace<Real>* SlaveFace() {return slave_face;};
  inline ContactHostGlobalID* SlaveFace_Global_ID() {return &slave_face_global_id;};
  inline int NumVertices() {return num_vertices;};
  inline ContactFaceCoverageVertex* Head() { return head; };
  void AddVertex(Real x, Real y);
  void AddVertex(ContactFaceCoverageVertex* v);

  inline Real& Scalar_Var( VariableHandle vh ) {return DataArray[vh];};
  inline Real* Vector_Var( VariableHandle vh ) 
    { return (DataArray+NUMBER_SCALAR_VARS+3*vh); };

  void Connect_SlaveFace( ContactTopologyEntityList& );
  void Connect_SlaveFace( ContactFace<Real>* f) {slave_face=f;};
  int Set_SlaveFaceEntityData();
  inline entity_data* SlaveFaceEntityData() {return &slave_face_entity_data;};

  inline int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline void Initialize_Memory() {std::memset(DataArray_Buffer(), 0, DataArray_Length()*sizeof(Real));};

  // Parallel packing/unpacking functions
#ifndef CONTACT_NO_MPI
  int  Size();
  void Pack( char* buffer );
  void Unpack( char* buffer );
  void Copy( ContactFaceCoverageInteraction* src );
#endif

  // Restart Pack/Unpack functions
  int  Restart_Size();
  void Restart_Pack( Real* buffer );
  void Restart_Unpack( Real* buffer );
  
  int Data_Size();
  
 protected:

 private:
  Real DataArray[NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS+1];
  ContactFace<Real>* slave_face;
  ContactHostGlobalID slave_face_global_id;
  entity_data slave_face_entity_data;
  int num_vertices;
  ContactFaceCoverageVertex* head;
  ContactFaceCoverageVertex* tail;
};

#endif // ContactFaceCoverageInteraction_h_
