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


#include "ContactFaceBlock.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSearch.h"
#include "ContactSearchData.h"
#include "ContactNode.h"
#include "ContactEdge.h"
#include "ContactFace.h"
#include "ContactElement.h"
#include "ContactTriFaceL3.h"
#include "ContactTopology.h"
#include "ContactFaceFaceInteraction.h"
#include "Contact_Defines.h"
#include "contact_tolerances.h"

#include <new>
#include <cmath>
#include <limits>

template<typename DataType>
struct PointStruct {
  DataType x, y, z;
  DataType xm, ym, zm;
};

template<typename DataType>
struct PolyStruct {
  int np;
  PointStruct<DataType> p[20];
};

template<typename DataType>
struct PlaneStruct {
  PointStruct<DataType> normal;  
  DataType d;
};

#define V3_Dot(a,b) \
        ((a)->x * (b)->x + (a)->y * (b)->y + (a)->z * (b)->z)

template<typename DataType>
ContactFaceFaceInteraction<DataType>*
ContactSearch::Face_Face_Search(ContactFace<DataType>* slave_face, 
                                ContactFace<DataType>* master_face, 
                                ContactElem<DataType>* element,
                                VariableHandle POSITION,
                                ContactFixedSizeAllocator *allocators)
{
  ContactFaceFaceInteraction<DataType>* cffi = NULL;
  typedef PointStruct<DataType> Point;
  typedef PolyStruct<DataType> Poly;
  typedef PlaneStruct<DataType> Plane;

  using std::abs;
  using std::sqrt;

  int   i, j, k, in_cnt;
  int   in[6], out[6];
  int   num_area=0;
  int   ifaceedge[20], iedge_m[20];
  DataType  area_s[42], area_m[42];
  Poly  poly0, poly1, poly2;
  Poly* p=&poly0;
  Poly* q=&poly1;
  Point  p0;
  Plane  face_plane;
  Plane* plane=&face_plane;
  ContactFace<DataType>* face = slave_face;
  VariableHandle FACE_NORMAL = 
    search_topology->Variable_Handle( ContactTopology::Face_Normal );
  
  for (i=0; i<20; ++i) {
    ifaceedge[i] = 0;
    iedge_m[i]   = 0;
  }
  
  p->np = face->Nodes_Per_Face();
  for (i=0; i<face->Nodes_Per_Face(); ++i) {
    ContactNode<DataType>* node = face->Node(i);
    DataType* position = node->Variable(POSITION);
    p->p[i].x = position[0];
    p->p[i].y = position[1];
    p->p[i].z = position[2];
  }

  //============================================
  // FIRST, CHECK TO SEE IF ALL THE VERTICES OF
  // THE FACE ARE INSIDE OR OUTSIDE THE ELEMENT
  //============================================
  for (i=0; i<6; ++i) {
    in[i]  = 0;
    out[i] = 0;
  }
  for (i=0; i<element->Faces_Per_Element(); ++i) {
    ContactNode<DataType>* node = element->Face(i)->Node(0);
    DataType* position = node->Variable(POSITION);
    p0.x               = position[0];
    p0.y               = position[1];
    p0.z               = position[2];
    DataType* normal   = element->Face(i)->Variable(FACE_NORMAL);
    plane->normal.x    = normal[0];
    plane->normal.y    = normal[1];
    plane->normal.z    = normal[2];
    plane->d           = -V3_Dot(&(plane->normal), &p0);
    for (j=0; j<face->Nodes_Per_Face(); ++j) {
      DataType value = V3_Dot(&(p->p[j]), &plane->normal) + plane->d;
      if (value > THICK_PLANE_TOL) {
        out[i]++;
      } else {
        in[i]++;
      }
    }
  }

  for (i=0; i<element->Faces_Per_Element(); ++i) {
    if (out[i]==face->Nodes_Per_Face()) return NULL;
  }

  int all_planar = 0;
  for (in_cnt=0, i=0; i<element->Faces_Per_Element(); ++i) {
    all_planar += element->Face(i)->IsPlanar(POSITION);
    in_cnt     += in[i];
  }

  if (all_planar==element->Faces_Per_Element()) {

    //=========================================================
    // If all the faces of the element are planar, then use
    // them as clipping planes to get the intersecting polygon
    //=========================================================
    if (in_cnt!=element->Faces_Per_Element()*face->Nodes_Per_Face()) {
      p = &poly0;
      q = &poly1;
      Poly*  r;
      Point* u;
      Point* v;
      Point* w;
      DataType t, tu, tv, tw;
      for (i=0; i<element->Faces_Per_Element(); ++i) {
        ContactNode<DataType>* node = element->Face(i)->Node(0);
        DataType* position = node->Variable(POSITION);
        p0.x               = position[0];
        p0.y               = position[1];
        p0.z               = position[2];
        DataType* normal   = element->Face(i)->Variable(FACE_NORMAL);
        plane->normal.x    = normal[0];
        plane->normal.y    = normal[1];
        plane->normal.z    = normal[2];
        plane->d           = -V3_Dot(&(plane->normal), &p0);
        //========================
        // START WITH U=VERT[N-1]
        //========================
        q->np = 0;
        u     = &(p->p[p->np-1]);
        tu    = V3_Dot(u, &plane->normal)+plane->d;
        for (j=0, k=p->np; k>0; k--, u=v, tu=tv, ++j) {
          //===================================================
          // ON OLD POLYGON (P), U IS PREVIOUS VERTEX, V IS
          // CURRENT VERTEX, TV IS POSITIVE IF VERTEX V IS OUT
          //===================================================
          v  = &(p->p[j]);
          tv = V3_Dot(v, &plane->normal)+plane->d;
          if ((tu>THICK_PLANE_TOL && tv<-THICK_PLANE_TOL) ||  // V is in; U is out
              (tu<-THICK_PLANE_TOL && tv>THICK_PLANE_TOL)) {  // U is in; V is out
            //=================================================
            // EDGE CROSSES PLANE; ADD INTERSECTION POINT TO Q
            //=================================================
            w = &(q->p[q->np]);
            if(abs(tu) < abs(tv)) {
              t    = tu/(tu-tv);
              w->x = u->x+t*(v->x-u->x);
              w->y = u->y+t*(v->y-u->y);
              w->z = u->z+t*(v->z-u->z);
            }
            else {
              t    = tv/(tv-tu);
              w->x = v->x+t*(u->x-v->x);
              w->y = v->y+t*(u->y-v->y);
              w->z = v->z+t*(u->z-v->z);
            }
            q->np++;
          }
          if (tv<-THICK_PLANE_TOL) { // V is in
            //==============================
            // VERTEX V IS IN, COPY IT TO Q
            //==============================
            w    = &(q->p[q->np]);
            w->x = v->x;
            w->y = v->y;
            w->z = v->z;
            q->np++;
          }
          if(tv >= -THICK_PLANE_TOL && tv <= THICK_PLANE_TOL) { // V is "on the plane"
            //=============================================================
            // VERTEX V IS ON THE "THICK PLANE", COPY IT TO Q IF NECESSARY
            //=============================================================
            w  = (j==p->np-1) ? &(p->p[0]) : &(p->p[j+1]);
            tw = V3_Dot(w, &plane->normal)+plane->d;
            if ((tu < -THICK_PLANE_TOL) ||                                                    // U is in 
                (tu >= -THICK_PLANE_TOL && tu <= THICK_PLANE_TOL && tw < -THICK_PLANE_TOL) || // U is "on the plane"; W is in
                (tu>THICK_PLANE_TOL && tw<-THICK_PLANE_TOL)) {                                // U is out; W is in
              w    = &(q->p[q->np]);
              w->x = v->x;
              w->y = v->y;
              w->z = v->z;
              q->np++;
            }
          }
        }
        r = p;
        p = q;
        q = r;
        if (p->np<3) {
          p->np = 0;
          break;
        }
      }
    }
 
    if (p->np>0) {
#ifdef COMPUTE_CENTROID_AND_LOCAL_EDGE_COORDS
      //=========================================================
      // If there is an intersecting polygon, it is stored in the
      // global coordinate system so calculate the centroid and
      // convert to the face and element local coordinate system
      //=========================================================
      DataType xc = 0.0;
      DataType yc = 0.0;
      DataType zc = 0.0;
      for (i=0; i<p->np; ++i) {
        xc += p->p[i].x;
        yc += p->p[i].y;
        zc += p->p[i].z;
      }
      xc /= p->np;
      yc /= p->np;
      zc /= p->np;
      DataType xbar  = 0.0;
      DataType ybar  = 0.0;
      DataType zbar  = 0.0;
      DataType tarea = 0.0;
      for (i=0; i<p->np; ++i) {
        int  i1    = i;
        int  i2    = (i1+1)%p->np;
        DataType x1    = p->p[i1].x;
        DataType y1    = p->p[i1].y;
        DataType z1    = p->p[i1].z;
        DataType x2    = p->p[i2].x;
        DataType y2    = p->p[i2].y;
        DataType z2    = p->p[i2].z;
        DataType area  = 0.5*sqrt(((y2-y1)*(zc-z1)-(yc-y1)*(z2-z1))*
                              ((y2-y1)*(zc-z1)-(yc-y1)*(z2-z1))+
                              ((z2-z1)*(xc-x1)-(zc-z1)*(x2-x1))*
                              ((z2-z1)*(xc-x1)-(zc-z1)*(x2-x1))+
                              ((x2-x1)*(yc-y1)-(xc-x1)*(y2-y1))*
                              ((x2-x1)*(yc-y1)-(xc-x1)*(y2-y1)));
        DataType area3 = area/3.0;
        xbar      += x1*area3 + x2*area3 + xc*area3;
        ybar      += y1*area3 + y2*area3 + yc*area3;
        zbar      += z1*area3 + z2*area3 + zc*area3;
        tarea     += area;
      }

      if(tarea > SMALL_AREA_TOL) { // avoid ctc polygon of null (small) area

        xbar /= tarea;
        ybar /= tarea;
        zbar /= tarea;
#else
      {
#endif
        num_area = p->np;
        DataType local_coords[6];
        DataType global_coords[6];
        for (i=0; i<p->np; ++i) {
          global_coords[0] = p->p[i].x;
          global_coords[1] = p->p[i].y;
          global_coords[2] = p->p[i].z;
          face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
          area_s[2*i+0] = local_coords[0];
          area_s[2*i+1] = local_coords[1];
          // If the element is a prism, then we can just compute the area_m local_coords directly on the master_face.
          // This is particularly useful for tri faces because Compute_Local_Coordinates is implemented more efficiently for the tri
          // face than for the general wedge element. The element is a prism if (a) the master_face is planar, and (b) the element
          // is formed by extruding the master_face in the direction of the face normal (i.e. the normal at the center of the face).
          master_face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
          area_m[2*i+0] = local_coords[0];
          area_m[2*i+1] = local_coords[1];
        }
        i = p->np;
#ifdef COMPUTE_CENTROID_AND_LOCAL_EDGE_COORDS
        global_coords[0] = xbar;
        global_coords[1] = ybar;
        global_coords[2] = zbar;
        face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
        area_s[2*i+0] = local_coords[0];
        area_s[2*i+1] = local_coords[1];
        master_face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
        area_m[2*i+0] = local_coords[0];
        area_m[2*i+1] = local_coords[1];
 
        for (i=0; i<num_area; ++i) {
          int  i1 = i;
          int  i2 = (i1+1)%num_area;
          DataType local_edge_coords[4];
          local_edge_coords[0] = area_s[i1*2];
          local_edge_coords[1] = area_s[i1*2+1];
          local_edge_coords[2] = area_s[i2*2];
          local_edge_coords[3] = area_s[i2*2+1];
          ifaceedge[i]         = face->Get_Edge_Number(local_edge_coords)+1;
          local_edge_coords[0] = area_m[i1*2];
          local_edge_coords[1] = area_m[i1*2+1];
          local_edge_coords[2] = area_m[i2*2];
          local_edge_coords[3] = area_m[i2*2+1];
          iedge_m[i]           = master_face->Get_Edge_Number(local_edge_coords)+1;
        }
#else
        area_s[2*i+0] = 0;
        area_s[2*i+1] = 0;
        area_m[2*i+0] = 0;
        area_m[2*i+1] = 0;
        for (i=0; i<num_area; ++i) {
          ifaceedge[i]         = 0;
          iedge_m[i]           = 0;
        }
#endif
        cffi = ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction(
              allocators[ALLOC_ContactFaceFaceInteraction],
              slave_face, master_face, num_area,
              ifaceedge, iedge_m, area_s, area_m );
      }
    }

  } else {
 
    //=========================================================
    // If all the faces of the element are not planar, then 
    // find all the edge-face intersections in the slave face
    // local coordinate system.  Then order these intersections
    // in a ccw direction and remove duplicate entries.
    //=========================================================
    DataType local_coords[4];
    DataType global_coords[3];

    q     = &poly1;
    q->np = 0;
 
    //=======================================================
    // Find all the vertices of the face that are completely
    // inside the element and add them to the poly
    //=======================================================
    for (i=0; i<face->Nodes_Per_Face(); ++i) {
      DataType* position = face->Node(i)->Variable(POSITION);
      global_coords[0] = position[0];
      global_coords[1] = position[1];
      global_coords[2] = position[2];
      element->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
      if (element->Is_Local_Coordinates_Inside_Element(local_coords)) {
        q->p[q->np].xm = local_coords[0];
        q->p[q->np].ym = local_coords[1];
        q->p[q->np].zm = local_coords[2];
        face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
        q->p[q->np].x = local_coords[0];
        q->p[q->np].y = local_coords[1];
        q->p[q->np].z = 1-local_coords[0]-local_coords[1];
        q->np++;
      }
    }

    //========================================
    // Find all the face-edge/element-face
    // intersections and add them to the poly
    //========================================
    for (i=0; i<element->Faces_Per_Element(); ++i) {
      ContactFace<DataType>* Face = element->Face(i);
      for (j=0; j<face->Edges_Per_Face(); ++j) {
        ContactEdge<DataType>* Edge = face->Edge(j);
        if (Face->FaceEdge_Intersection(POSITION, Edge, global_coords)) {
          element->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
          q->p[q->np].xm = local_coords[0];
          q->p[q->np].ym = local_coords[1];
          q->p[q->np].zm = local_coords[2];
          face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
          q->p[q->np].x = local_coords[0];
          q->p[q->np].y = local_coords[1];
          q->p[q->np].z = 1-local_coords[0]-local_coords[1];
          q->np++;
        }
      }
    }
 
    //========================================
    // Find all the element-edge/face
    // intersections and add them to the poly
    //========================================
    ContactFace<DataType>* Face = face;
    for (i=0; i<element->Edges_Per_Element(); ++i) {
      ContactEdge<DataType>* Edge = element->Edge(i);
      if (Face->FaceEdge_Intersection(POSITION, Edge, global_coords)) {
        element->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
        q->p[q->np].xm = local_coords[0];
        q->p[q->np].ym = local_coords[1];
        q->p[q->np].zm = local_coords[2];
        face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
        q->p[q->np].x = local_coords[0];
        q->p[q->np].y = local_coords[1];
        q->p[q->np].z = 1-local_coords[0]-local_coords[1];
        q->np++;
      }
    }
 
    if(q->np < 3) p->np = 0; // PJSA 1/27/2015: previously was if(q->np==0) ...
    else {

      //============================================
      // Order the points of the poly in ccw order
      // based on their positions on the slave face
      //============================================
      int  imin=-1;
      DataType ymin=2.0;
      int*  hit = new int [q->np];
      DataType* vec = new DataType [q->np];
      for (i=0; i<q->np; ++i) hit[i] = 0;
 
      switch (face->FaceType()) {
      case SHELLTRIFACEL3:
      case TRIFACEL3: {
        DataType t3_node_positions[3][3];
        t3_node_positions[0][0] = 0.0;
        t3_node_positions[0][1] = 0.0;
        t3_node_positions[0][2] = 0.0;
        t3_node_positions[1][0] = 1.0;
        t3_node_positions[1][1] = 0.0;
        t3_node_positions[1][2] = 0.0;
        t3_node_positions[2][0] = 0.0;
        t3_node_positions[2][1] = 1.0;
        t3_node_positions[2][2] = 0.0;
        ContactTriFaceL3<DataType>* t3face = static_cast<ContactTriFaceL3<DataType>*>(face);
        for (i=0; i<q->np; ++i) {
          local_coords[0] = q->p[i].x;
          local_coords[1] = q->p[i].y;
          local_coords[2] = q->p[i].z;
          t3face->Compute_Global_Coords(t3_node_positions, local_coords, global_coords);
          if (global_coords[1] < ymin) {
            imin = i;
            ymin = global_coords[1];
          }
        }
        POSTCONDITION(imin>=0);
        local_coords[0] = q->p[imin].x;
        local_coords[1] = q->p[imin].y;
        local_coords[2] = q->p[imin].z;
        t3face->Compute_Global_Coords(t3_node_positions, local_coords, global_coords);
        DataType dx0 = global_coords[0];
        DataType dy0 = global_coords[1];
        for (i=0; i<q->np; ++i) {
          if (i==imin) continue;
          local_coords[0] = q->p[i].x;
          local_coords[1] = q->p[i].y;
          local_coords[2] = q->p[i].z;
          t3face->Compute_Global_Coords(t3_node_positions, local_coords, global_coords);
          DataType dx  = global_coords[0] - dx0;
          DataType dy  = global_coords[1] - dy0;
          DataType mag = dx*dx + dy*dy;
          if( mag > 0.0 ) {
            mag = sqrt(mag);
            dx /= mag;
            dy /= mag;
          } else {
            mag = 1.0;
            dx  = 0.0;
            dy  = 0.0;
          }
          vec[i] = dx;
        }
        } break;
      case SHELLQUADFACEL4:
      case QUADFACEL4:
        for (i=0; i<q->np; ++i) {
          if (q->p[i].y < ymin) {
            imin = i;
            ymin = q->p[i].y;
          }
        }
        POSTCONDITION(imin>=0);
        for (i=0; i<q->np; ++i) {
          if (i==imin) continue;
          DataType dx  = q->p[i].x - q->p[imin].x;
          DataType dy  = q->p[i].y - q->p[imin].y;
          DataType mag = dx*dx + dy*dy;
          if( mag > 0.0 ) {
            mag = sqrt(mag);
            dx /= mag;
            dy /= mag;
          } else {
            mag = 1.0;
            dx  = 0.0;
            dy  = 0.0;
          }
          vec[i] = dx;
        }
        break;
      default:
        PRECONDITION(0);
        break;
      }

      if(imin < 0) p->np = 0; // PJSA 27/1/2015: this can happen for highly distorted meshes,
                              //                 e.g. in a divergent Newton iterate.
      else {
        for (i=0; i<q->np; ++i) {
          if( i==imin ) continue;
          DataType dx = q->p[i].x - q->p[imin].x;
          DataType dy = q->p[i].y - q->p[imin].y;
          DataType d  = dx*dx + dy*dy;
          if (d <= 1.0e-10) hit[i] = 1;
        }
        hit[imin] = 1;
        int kmin  = imin;

        p     = &poly2;
        p->np = 1;
        p->p[0] = q->p[imin];
        for (i=0; i<q->np; ++i) {
          int  jmin = -1;
          DataType vmax = -2.0;
          for (j=0; j<q->np; ++j) {
            if( !hit[j] && vec[j]>vmax ) {
               vmax = vec[j];
               jmin = j;
            }
          }
          if (jmin>=0) {
            hit[jmin] = 1;
            int equivalent=0;
            DataType dx = q->p[jmin].x - q->p[kmin].x;
            DataType dy = q->p[jmin].y - q->p[kmin].y;
            DataType d  = dx*dx + dy*dy;
            if (d <= 1.0e-10) equivalent=1;
            if (!equivalent) {
              kmin = jmin;
              p->p[p->np] = q->p[jmin];
              p->np++;
            }
          }
        }
      }
      delete [] hit;
      delete [] vec;
      if (p->np<3) p->np = 0;
    }

    if (p->np>0) {
#ifdef COMPUTE_CENTROID_AND_LOCAL_EDGE_COORDS
      //=========================================================
      // If there is an intersecting polygon, it is stored in the
      // local coordinate system so calculate the centroid
      //=========================================================
      DataType xc = 0.0;
      DataType yc = 0.0;
      DataType zc = 0.0;
      for (i=0; i<p->np; ++i) {
        local_coords[0] = p->p[i].x;
        local_coords[1] = p->p[i].y;
        local_coords[2] = p->p[i].z;
        face->Compute_Global_Coordinates(POSITION, local_coords, global_coords);
        xc += global_coords[0];
        yc += global_coords[1];
        zc += global_coords[2];
      }
      xc /= p->np;
      yc /= p->np;
      zc /= p->np;
      DataType xbar  = 0.0;
      DataType ybar  = 0.0;
      DataType zbar  = 0.0;
      DataType tarea = 0.0;
      for (i=0; i<p->np; ++i) {
        int  i1    = i;
        int  i2    = (i1+1)%p->np;
        local_coords[0] = p->p[i1].x;
        local_coords[1] = p->p[i1].y;
        local_coords[2] = p->p[i1].z;
        face->Compute_Global_Coordinates(POSITION, local_coords, global_coords);
        DataType x1    = global_coords[0];
        DataType y1    = global_coords[1];
        DataType z1    = global_coords[2];
        local_coords[0] = p->p[i2].x;
        local_coords[1] = p->p[i2].y;
        local_coords[2] = p->p[i2].z;
        face->Compute_Global_Coordinates(POSITION, local_coords, global_coords);
        DataType x2    = global_coords[0];
        DataType y2    = global_coords[1];
        DataType z2    = global_coords[2];
        DataType area  = 0.5*sqrt(((y2-y1)*(zc-z1)-(yc-y1)*(z2-z1))*
                              ((y2-y1)*(zc-z1)-(yc-y1)*(z2-z1))+
                              ((z2-z1)*(xc-x1)-(zc-z1)*(x2-x1))*
                              ((z2-z1)*(xc-x1)-(zc-z1)*(x2-x1))+
                              ((x2-x1)*(yc-y1)-(xc-x1)*(y2-y1))*
                              ((x2-x1)*(yc-y1)-(xc-x1)*(y2-y1)));
        DataType area3 = area/3.0;
        xbar      += x1*area3 + x2*area3 + xc*area3;
        ybar      += y1*area3 + y2*area3 + yc*area3;
        zbar      += z1*area3 + z2*area3 + zc*area3;
        tarea     += area;
      }

      if(tarea > SMALL_AREA_TOL) { // avoid ctc polygon of null (small) area

        xbar /= tarea;
        ybar /= tarea;
        zbar /= tarea;
#else
      {
#endif
        num_area = p->np;
        for (i=0; i<p->np; ++i) {
          area_s[i*2+0] = p->p[i].x;
          area_s[i*2+1] = p->p[i].y;
          area_m[i*2+0] = p->p[i].xm;
          area_m[i*2+1] = p->p[i].ym;
        }

#ifdef COMPUTE_CENTROID_AND_LOCAL_EDGE_COORDS
        global_coords[0] = xbar;
        global_coords[1] = ybar;
        global_coords[2] = zbar;
        face->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
        area_s[i*2+0] = local_coords[0];
        area_s[i*2+1] = local_coords[1];
        element->Compute_Local_Coordinates(POSITION, global_coords, local_coords);
        area_m[i*2+0] = local_coords[0];
        area_m[i*2+1] = local_coords[1];
 
        for (i=0; i<num_area; ++i) {
          int  i1 = i;
          int  i2 = (i1+1)%num_area;
          DataType local_edge_coords[4];
          local_edge_coords[0] = area_s[i1*2];
          local_edge_coords[1] = area_s[i1*2+1];
          local_edge_coords[2] = area_s[i2*2];
          local_edge_coords[3] = area_s[i2*2+1];
          ifaceedge[i]         = face->Get_Edge_Number(local_edge_coords)+1;
          local_edge_coords[0] = area_m[i1*2];
          local_edge_coords[1] = area_m[i1*2+1];
          local_edge_coords[2] = area_m[i2*2];
          local_edge_coords[3] = area_m[i2*2+1];
          iedge_m[i]           = master_face->Get_Edge_Number(local_edge_coords)+1;
        }
#else
        area_s[i*2+0] = 0;
        area_s[i*2+1] = 0;
        area_m[i*2+0] = 0;
        area_m[i*2+1] = 0;
        for (i=0; i<num_area; ++i) {
          ifaceedge[i] = 0;
          iedge_m[i] = 0;
        }
#endif
        cffi = ContactFaceFaceInteraction<DataType>::new_ContactFaceFaceInteraction(
              allocators[ALLOC_ContactFaceFaceInteraction],
              slave_face, master_face, num_area,
              ifaceedge, iedge_m, area_s, area_m );
      }
    }
  }
  return cffi;
}

