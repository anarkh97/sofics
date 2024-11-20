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


#include "ContactHexOverlap.h"
#include "contact_tolerances.h"

// the hex intersection proceeds by accepting tetted hexes and then
// intersecting all the tets, the tets are intersected by trimming all
// of their edges agains each other to get a set of points that define
// the intersection between two tets. These points are then convex hulled
// and the resulting volume accumulated.

// the convex hull is slow and would benefit from a judicious re-write or 
// the application of the q-hull library from the University of Illinois

//entry point that runs tets agains each other and accumulates volume ovelap
Real intersection_volume(Real thex1[24][4][3],Real thex2[24][4][3]){
  Real points[MAX_POINTS][3];
  int ptnct = 0;
  
  Real volume_sum = 0.;
  for(int i = 0; i < 24; i ++){
    for(int j = 0; j < 24; j++){
      ptnct = 0;
      intersect_for_points(thex2[i], thex1[j],points, ptnct);
      intersect_for_points(thex1[j], thex2[i],points, ptnct);
      if(ptnct>3){
	face_collection * fc = hull_point_collection(points,ptnct);
	if(fc)volume_sum += volume(fc,points);
	if(fc)  delete fc;
      }
    }
  }
  return volume_sum;
}

// a slow convex hull routine
face_collection * hull_point_collection(Real the_points[MAX_POINTS][3],int n_pts){
    bool used[12];
    static int tet_faces[4][3] = {{0,1,3},  
				  {1,2,3},	
				  {0,3,2},
				  {0,2,1}};
    for(int i = 0; i < n_pts; ++i){
      used[i] = false;
    }
  // create first three faces using any 3 points plus a non co-planar point
  // add its faces to the free face list which is kept sorted. by their ordered verices
  //faces should have normals, center points, and ordered pointers to vertices

  // loop over remaining points
    // loop over free faces
      // if point is positive from face
        // create four new faces 
        // (if they already exist when added, do not add and remove existing) 

  // create a first tet using any 3 points plus a non co-planar point
    
    Real norm[3];
    three_point_normal(the_points[0],the_points[1],the_points[2],norm);
    used[0] = true;
    used[1] = true;
    used[2] = true;


  int  ti[4];
  bool found = false;
  Real tnorm[3];
  
  for(int i = 3; i < n_pts && found == false; ++i){
    vector_difference(the_points[i], the_points[0], tnorm);
    normalize_vector(tnorm);

    Real dot_prod = dot_vectors(norm,tnorm);
    if((dot_prod) > SMALL_NUMBER){
      ti[0] = 0; ti[1] = 1; ti[2] = 2; ti[3] = i;
      found = true;
      used[i] = true;
    }
    if((dot_prod) < - SMALL_NUMBER){
      ti[0] = 1; ti[1] = 0; ti[2] = 2; ti[3] = i;
      found = true;
      used[i] = true;
    }
  }
  if(!found){
    //coplanar points commonly occur in file made by humans.
    return NULL;
  }
  face * first_face = new face(the_points,
			       ti[tet_faces[0][0]],
			       ti[tet_faces[0][1]],
			       ti[tet_faces[0][2]]);

  for(int i = 1; i < 4; ++i){
    face * a_face = new face(the_points,
			     ti[tet_faces[i][0]],
			     ti[tet_faces[i][1]],
			     ti[tet_faces[i][2]]);
    
    first_face->add(a_face);
  }

  for(int i = 3; i < n_pts; ++i){
    if(!used[i]){
      face * tface = first_face;
      while(tface){
	if(tface->active){
	  vector_difference(the_points[i],tface->center,tnorm);
	  normalize_vector(tnorm);
	  
	  Real dot_prod = dot_vectors(tnorm,tface->normal);
	  
	  if((dot_prod) > SMALL_NUMBER /** 1.0e-2*/){
	    for(int k = 0; k < 3; ++k){
	      ti[k] = tface->indices[k];
	    }
	    ti[3] = i;
	    for(int lct = 0; lct < 4; ++lct){
	      face * a_face = new face(the_points,
				       ti[tet_faces[lct][0]],
				       ti[tet_faces[lct][1]],
				       ti[tet_faces[lct][2]]);
	      
	      first_face->add(a_face);
	    }
	    tface = first_face;
	  }
	  else{
	    //starting over here because tface get messed up in the add function
	    // should not really have to do this.
	    tface = tface->next;
	  }
	}
	else{
	  tface = tface->next;
	}
      } 
    }
  }
  face_collection * fc = new face_collection(n_pts, first_face);
  return fc;
}

bool lower_than(face*f1,face*f2){
  for(int i = 0 ; i < 3; i ++){
    if(f1->sorted_indices[i] > f2->sorted_indices[i]){
      return false;
    }
    else if(f1->sorted_indices[i] < f2->sorted_indices[i]){
      return true;
    }
  }
  return false;
}


bool equal(face*f1,face*f2){
  for(int i = 0; i < 3; i ++){
    if(f1->sorted_indices[i] != f2->sorted_indices[i])return false; 
  }
  return true;
}

//Trims the edges of every tet agains the planes defined by the faces of the other tet.
//Each tet operated on by the other to produce a set of points that bound the intersection
// of the tets.
void intersect_for_points(Real t_hex1[4][3], Real t_hex2[4][3], Real result [MAX_POINTS][3], int & ptcnt){
  
  static int tet_edges[6][2] = {{0,1},
				{1,2},
				{2,0},
				{0,3},
				{1,3},
				{2,3}};
  
  static int tet_faces[4][3] = {{0,1,3},  
				{1,2,3},	
				{0,3,2},
				{0,2,1}};
  for( int i = 0; i < 6; ++i ){
    Real Rd [3];//directional vector
    Rd[0] = t_hex2[tet_edges[i][1]][0] - t_hex2[tet_edges[i][0]][0];
    Rd[1] = t_hex2[tet_edges[i][1]][1] - t_hex2[tet_edges[i][0]][1];
    Rd[2] = t_hex2[tet_edges[i][1]][2] - t_hex2[tet_edges[i][0]][2];
    Real tfar = norm_vector(Rd);
    Real tfar_orig = tfar;
    Real small_length = tfar_orig * REL_TANG_TOL;
    Real ls2 = small_length*small_length;
    Real tnear = 0.;
    Rd[0] /= tfar;
    Rd[1] /= tfar;
    Rd[2] /= tfar;
    Real Ro[3];//origin
    Ro[0] = t_hex2[tet_edges[i][0]][0];
    Ro[1] = t_hex2[tet_edges[i][0]][1];
    Ro[2] = t_hex2[tet_edges[i][0]][2];
    bool intersecting = true;
    for(int nfaces = 0; (nfaces < 4) && intersecting /*&& (tfar > tnear)*/; nfaces++){
      Real Pn[3];//plane normal
      Real d1[3];
      Real d2[3];
      d1[0] = t_hex1[tet_faces[nfaces][1]][0] - t_hex1[tet_faces[nfaces][0]][0];
      d1[1] = t_hex1[tet_faces[nfaces][1]][1] - t_hex1[tet_faces[nfaces][0]][1];
      d1[2] = t_hex1[tet_faces[nfaces][1]][2] - t_hex1[tet_faces[nfaces][0]][2];
      
      d2[0] = t_hex1[tet_faces[nfaces][2]][0] - t_hex1[tet_faces[nfaces][0]][0];
      d2[1] = t_hex1[tet_faces[nfaces][2]][1] - t_hex1[tet_faces[nfaces][0]][1];
      d2[2] = t_hex1[tet_faces[nfaces][2]][2] - t_hex1[tet_faces[nfaces][0]][2];
      
      cross_vectors(d1,d2,Pn);
      normalize_vector(Pn);
      
      Real d = -dot_vectors(Pn,t_hex1[tet_faces[nfaces][0]]);//distance from plane to origin
      Real vn = d+dot_vectors(Pn,Ro);// normal distance between plane and line origin
      Real vd = dot_vectors(Pn,Rd); // cosine of angle between plane normal and line direction
      if(fabs(vd) <= REL_TANG_TOL){// parallel to plane
	if(vn > small_length){
	  intersecting = false;
	  tnear = tfar;
	}	
      }
      else if (vd > 0.){
	Real t = -vn/vd;
	if(t < tfar) tfar = t;
      }
      else{
	Real t = -vn/vd;
	if(t > tnear) tnear = t;
      }
    }
    if(intersecting){
      if(tfar > tnear || (fabs(tfar-tnear)/tfar_orig < REL_TANG_TOL)){
	
	result[ptcnt][0] = Ro[0] + Rd[0]*tnear; 
	result[ptcnt][1] = Ro[1] + Rd[1]*tnear; 
	result[ptcnt][2] = Ro[2] + Rd[2]*tnear;
	bool found = false;
	for(int j = 0; j < ptcnt  && !found; ++j){
	  if(squared_vector_difference(result[j],result[ptcnt]) < ls2)found = true;
	}
	if (!found) {ptcnt ++;}
	
	result[ptcnt][0] = Ro[0] + Rd[0]*tfar; 
	result[ptcnt][1] = Ro[1] + Rd[1]*tfar; 
	result[ptcnt][2] = Ro[2] + Rd[2]*tfar; 
	found = false;
	for(int j = 0; j < ptcnt  && !found; ++j){
	  if(squared_vector_difference(result[j],result[ptcnt]) < ls2)found = true;
	}
	if (!found) {ptcnt ++;}
      }
    }
  }
}

void three_point_normal(Real v1 [3],Real v2 [3], Real v3 [3],Real vout [3]){
  Real d1[3];
  Real d2[3];
  vector_difference(v2,v1,d1);
  vector_difference(v3,v1,d2);
  cross_vectors(d1,d2,vout);
  normalize_vector(vout);
}

//Volume calculation for a volume bounded by a collection of triangular faces all
// oriented outward
Real volume(face_collection * fc, Real points[MAX_POINTS][3]){
  Real volume = 0.;
  face * tface = fc->faces;
  while(tface){
    if(tface->active){
      Real cross_prod[3];
      Real sum[3];
      cross_vectors(points[tface->indices[0]],points[tface->indices[1]],cross_prod);
      sum[0] = cross_prod[0];
      sum[1] = cross_prod[1];
      sum[2] = cross_prod[2];
      cross_vectors(points[tface->indices[1]],points[tface->indices[2]],cross_prod);
      sum[0] += cross_prod[0];
      sum[1] += cross_prod[1];
      sum[2] += cross_prod[2];
      cross_vectors(points[tface->indices[2]],points[tface->indices[0]],cross_prod);
      sum[0] += cross_prod[0];
      sum[1] += cross_prod[1];
      sum[2] += cross_prod[2];
      Real mag = fabs(dot_vectors(tface->normal,sum));
      volume += mag*dot_vectors(points[tface->indices[0]],tface->normal);
    }
    tface = tface->next;
  }
  volume /=6.0;
  return volume;
}

void vector_difference(Real v1[3],Real v2[3],Real v3[3]){
  for(int i = 0; i < 3; i ++){
    v3[i] = v1[i] - v2[i];
  }
}

Real squared_vector_difference(Real v1[3],Real v2[3]){
  Real res = 0;
  for(int i = 0; i < 3; i ++){
    Real t = v1[i] - v2[i];
    res += t*t;
  }
  return res;
}

void scale_hex(Real arr[8][3] , Real val){
  for(int i = 0; i < 8; i ++){
    arr[i][0]*=val;
    arr[i][1]*=val;
    arr[i][2]*=val;
  }
}

Real min_hex_length(Real hex[8][3]){
  Real res = squared_vector_difference(hex[0],hex[1]);
  for(int i = 0; i < 7 ; i++){
    for (int j = i+1; j < 8; j++){
      Real tdif = squared_vector_difference(hex[i],hex[j]);
      if(tdif < res)res = tdif;
    }
  }
  if(res)return sqrt(res);
  return 0.;
}

void cross_vectors(Real vec_0_1[3], Real vec_0_2[3], Real cross[3]){
  cross[0] = vec_0_1[1]*vec_0_2[2] - vec_0_1[2]*vec_0_2[1];
  cross[1] = vec_0_1[2]*vec_0_2[0] - vec_0_1[0]*vec_0_2[2];
  cross[2] = vec_0_1[0]*vec_0_2[1] - vec_0_1[1]*vec_0_2[0];
}

Real dot_vectors(Real v1[3],Real v2[3]){
  return(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

Real norm_vector(Real vec[3]){
  Real tn = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
  if(tn){
    return sqrt(tn);
  }
  return 0.;
}

void normalize_vector(Real vec[3]){
  Real tn = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
  if(tn>0.){
   tn =  sqrt(tn);
   vec[0]/=tn;
   vec[1]/=tn;
   vec[2]/=tn;
  }
}

face::face(Real v[MAX_POINTS][3],int i1, int i2, int i3){
    if(i1 == i2 || i2 == i3 || i1 == i3){
      std::cout << "AAAACCCKKKK" << std::endl;
    }
    active = true;
    indices[0] = i1;
    indices[1] = i2;
    indices[2] = i3;
    
    center[0] = 0.;
    center[1] = 0.;
    center[2] = 0.;
    for(int i = 0; i < 3; i ++){
      center[i]+=v[i1][i];
      center[i]+=v[i2][i];
      center[i]+=v[i3][i];      
    }

    center[0]/=3.;
    center[1]/=3.;
    center[2]/=3.;
    
    three_point_normal(v[i1],v[i2],v[i3],normal);

    sorted_indices[0] = i1;
    sorted_indices[1] = i2;
    sorted_indices[2] = i3;
    for(int j = 0; j < 2; j++){
      for(int k = 0; k < 2; k++){
	if(sorted_indices[k+1] < sorted_indices[k]){
	  int t = sorted_indices[k];
	  sorted_indices[k] = sorted_indices[k+1];
	  sorted_indices[k+1] = t;
	}
      }
    }
    next = NULL;
  };
