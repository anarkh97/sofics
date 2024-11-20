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


#ifndef ContactFaceBlock_h_
#define ContactFaceBlock_h_

#include "ContactSearch.h"
#include "ContactFace.h"
#include "ContactEdge.h"
#include "ContactBlockEntityList.h"
#include "ContactBoundingBox.h"

class ContactTopology;

class ContactFaceBlock {

 public:

  ContactFaceBlock( ContactSearch::ContactFace_Type, int, int, int&, int*,
                    ContactTopology* );
  ContactFaceBlock( ContactSearch::ContactFace_Type, int, ContactTopology*);
  ~ContactFaceBlock();
  
  void Delete_Face_List();
  
  void Delete_Faces();
  void Add_Faces( int, int* );

  inline ContactSearch::ContactFace_Type Type() { return type;};
  inline ContactSearch::ContactEdge_Type EdgeType() { return edge_type;};
  inline int Entity_Key() { return entity_key; };
  inline int ID() { return id; };
  inline int Number_of_Faces() { return number_of_faces; };
  ContactBlockEntityList* FaceList() { return face_list; };
  ContactBoundingBox* LocalBoundingBox() { return &local_bounding_box; };
  ContactBoundingBox* GlobalBoundingBox() { return &global_bounding_box; };
  void ComputeBoundingBox(int, VariableHandle, VariableHandle, MPI_Comm&);

  void Insert_Face( ContactFace<Real>* );
  void Delete_Face( ContactFace<Real>* );
#ifndef CONTACT_NO_MPI
  void Insert_Face( char* );
  void Insert_Face_ForSecondary( char* );
#endif

 void SetMaster() {master=true;};
 void SetSlave()  {slave=true;};
 void SetMaster(bool val) {master=val;};
 void SetSlave(bool val)  {slave=val;};
 bool IsMaster()          {return master;};
 bool IsSlave()           {return slave;};

 private:
  
  ContactFaceBlock(ContactFaceBlock&);
  ContactFaceBlock& operator=(ContactFaceBlock&);
  
  ContactTopology* topology;  // owner topology object
  
  int number_of_faces;
  int num_faces_added;
  ContactSearch::ContactFace_Type type;
  ContactSearch::ContactEdge_Type edge_type;
  bool master;
  bool slave;
  int entity_key;
  int id;

  ContactBlockEntityList* face_list;
  
  ContactBoundingBox local_bounding_box;
  ContactBoundingBox global_bounding_box;

};

#endif //ContactFaceBlock_h_
