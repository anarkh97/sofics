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


#ifndef ContactHexOverlap_h_
#define ContactHexOverlap_h_


#include <iostream>


#include <cmath>
typedef double Real;
#define SMALL_NUMBER 1.0e-04
#define MAX_POINTS 100

class face;
class face_collection;

bool equal(face*,face*);
bool lower_than(face*,face*);
void normalize_vector(Real vec[3]);
void vector_difference(Real v1[3],Real v2[3],Real v3[3]);
void cross_vectors(Real vec_0_1[3], Real vec_0_2[3], Real cross[3]);
Real volume(face_collection * fc,Real pts[MAX_POINTS][3]);
void three_point_normal(Real v1 [3],Real v2 [3], Real v3 [3],Real vout [3]);
face_collection * hull_point_collection(Real pc [MAX_POINTS][3],int npts);
void intersect_for_points(Real t_hex1[4][3], Real t_hex2[4][3],Real pts[MAX_POINTS][3], int & ptcnt);
Real min_hex_length(Real hex[8][3]);
Real dot_vectors(Real v1[3],Real v2[3]);
Real norm_vector(Real vec[3]);
Real squared_vector_difference(Real v1[3],Real v2[3]);
Real intersection_volume(Real thex1[24][4][3],Real thex2[24][4][3]);
void write_tet(char * aname,Real t_hex[4][3]);

class face {
public:
  face(){
    next = NULL;
    center[0] = 0.;
    center[1] = 0.;
    center[2] = 0.;
    normal[0] = 0.;
    normal[1] = 0.;
    normal[2] = 0.;
    active = true;
  };
  face(Real v[MAX_POINTS][3],int i1, int i2, int i3);


  void add(face * incoming_face){
    face * last_face = this;
    face * tface = this;
    while(tface){
      if(tface->active){
	if(equal(tface,incoming_face)){
	  delete incoming_face;
	  tface->active = false;
	  return;
	}
      }
      last_face = tface;
      tface = tface->next;
    }
    last_face->next = incoming_face;
    return;
  };

  virtual ~face(){};
  face * next;
  Real center[3];
  Real normal[3];
  int sorted_indices[3];
  int indices[3];
  bool active;
};

class face_collection {
public:
  face_collection(){
  faces = NULL;
  npts = 0;
  }
  face_collection(int n_points, face * the_faces){
    npts = n_points;
    faces = the_faces;
  }

  void write(char * base_name,Real coords[MAX_POINTS][3]);
  void write (std::ostream * os ,Real coords[MAX_POINTS][3]);


  virtual ~face_collection(){      
    face * first_face = faces;
    while(first_face){
      face *next_face = first_face->next;
      delete first_face;
      first_face = next_face;
    }
  }
  face * faces;
  int npts;
};

#endif
