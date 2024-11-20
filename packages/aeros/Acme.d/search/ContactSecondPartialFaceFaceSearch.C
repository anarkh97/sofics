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
struct PointStructExp2 {
  DataType x, y, z;
  DataType xm, ym, zm;

  DataType *partial_x, *partial_y, *partial_z;
  DataType *partial_t;
  DataType *second_partial_x, *second_partial_y, *second_partial_z;
  DataType *second_partial_t;
  PointStructExp2() : partial_x(NULL), partial_y(NULL), partial_z(NULL), partial_t(NULL),
                      second_partial_x(NULL), second_partial_y(NULL), second_partial_z(NULL), second_partial_t(NULL) {}
  ~PointStructExp2() {
    if(partial_x) delete [] partial_x;
    if(partial_y) delete [] partial_y;
    if(partial_z) delete [] partial_z;
    if(partial_t) delete [] partial_t;
    if(second_partial_x) delete [] second_partial_x;
    if(second_partial_y) delete [] second_partial_y;
    if(second_partial_z) delete [] second_partial_z;
    if(second_partial_t) delete [] second_partial_t;
  }
};

template<typename DataType>
struct PolyStructExp2 {
  int np;
  PointStructExp2<DataType> p[20];
};

template<typename DataType>
struct PlaneStructExp2 {
  PointStructExp2<DataType> normal;  
  DataType d;
};

#define V3_Dot(a,b) \
        ((a)->x * (b)->x + (a)->y * (b)->y + (a)->z * (b)->z)

template<typename DataType>
ContactFaceFaceInteraction<DataType>*
ContactSearch::Second_Partial_Face_Face_Search(ContactFace<DataType>* slave_face, 
                                        ContactFace<DataType>* master_face, 
                                        ContactElem<DataType>* element,
                                        VariableHandle POSITION, Real tol,
                                        ContactFixedSizeAllocator *allocators)
{
  ContactFaceFaceInteraction<DataType>* cffi = NULL;
  typedef PointStructExp2<DataType> Point;
  typedef PolyStructExp2<DataType> Poly;
  typedef PlaneStructExp2<DataType> Plane;

  using std::abs;
  using std::sqrt;

  const int nbDer = 3*(slave_face->Nodes_Per_Face() + master_face->Nodes_Per_Face());

  int   i, j, k, l, m, in_cnt;
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

  // initialize the partial derivatives of slave_face nodes' x,y,z coordinates
  for (i=0; i<face->Nodes_Per_Face(); ++i) {
    p->p[i].partial_x = new DataType[nbDer];
    p->p[i].partial_y = new DataType[nbDer];
    p->p[i].partial_z = new DataType[nbDer];
    for (j=0; j<nbDer; ++j) p->p[i].partial_x[j] = p->p[i].partial_y[j] = p->p[i].partial_z[j] = 0;
    p->p[i].partial_x[3*i  ] = 1;
    p->p[i].partial_y[3*i+1] = 1;
    p->p[i].partial_z[3*i+2] = 1;
  }

  DataType (*partial_area_s)[42] = new DataType[nbDer][42];
  DataType (*partial_area_m)[42] = new DataType[nbDer][42];
  DataType (*second_partial_area_s)[42] = new DataType[nbDer*nbDer][42];
  DataType (*second_partial_area_m)[42] = new DataType[nbDer*nbDer][42];

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
      DataType (*master_face_dnormal)[3] = NULL;
      DataType (*master_face_d2normal)[3] = NULL;
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
        DataType (*plane_dnormal)[3] = NULL;
        DataType *plane_dd           = NULL;
        DataType (*plane_d2normal)[3] = NULL;
        DataType *plane_d2d           = NULL;
        //========================
        // START WITH U=VERT[N-1]
        //========================
        q->np = 0;
        u     = &(p->p[p->np-1]);
        tu    = V3_Dot(u, &plane->normal)+plane->d;
        for (j=0; j<20; ++j) {
          if(p->p[j].partial_t) { delete [] p->p[j].partial_t; p->p[j].partial_t = NULL; }
          if(q->p[j].partial_t) { delete [] q->p[j].partial_t; q->p[j].partial_t = NULL; }
          if(p->p[j].second_partial_t) { delete [] p->p[j].second_partial_t; p->p[j].second_partial_t = NULL; }
          if(q->p[j].second_partial_t) { delete [] q->p[j].second_partial_t; q->p[j].second_partial_t = NULL; }
        }
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
            if(!plane_dnormal) {
              if(!master_face_dnormal) {
                // compute derivatives of master_face normal w.r.t. global X-, Y- and Z-coordinates of master_face nodes
                master_face_dnormal = new DataType[3*master_face->Nodes_Per_Face()][3];
                master_face->Compute_Partial_Face_Normal(POSITION, master_face_dnormal);
                master_face_d2normal = new DataType[3*master_face->Nodes_Per_Face()*3*master_face->Nodes_Per_Face()][3];
                master_face->Compute_Second_Partial_Face_Normal(POSITION, master_face_d2normal);
              }
              // compute derivatives of plane->normal and plane->d w.r.t. global X-, Y- and Z-coordinates of master_face nodes
              plane_dnormal = new DataType[3*master_face->Nodes_Per_Face()][3];
              plane_dd = new DataType[3*master_face->Nodes_Per_Face()];
              element->Compute_Partial_Face_Normal(i, POSITION, FACE_NORMAL, master_face_dnormal, tol, plane_dnormal, plane_dd);
              plane_d2normal = new DataType[3*master_face->Nodes_Per_Face()*3*master_face->Nodes_Per_Face()][3];
              plane_d2d = new DataType[3*master_face->Nodes_Per_Face()*3*master_face->Nodes_Per_Face()];
              element->Compute_Second_Partial_Face_Normal(i, POSITION, FACE_NORMAL, master_face_dnormal, master_face_d2normal,
                                                          tol, plane_d2normal, plane_d2d); 
            }
            if(!u->partial_t) {
              // compute derivatives of tu w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              // due to symmetry of the second derivatives, only the upper triangular part is computed
              u->partial_t = new DataType[nbDer];
              for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                u->partial_t[l] = u->partial_x[l]*plane->normal.x + u->partial_y[l]*plane->normal.y + u->partial_z[l]*plane->normal.z;
              }
              u->second_partial_t = new DataType[nbDer*nbDer];
              if(u->second_partial_x != NULL) {
                for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                  for (m=l; m<3*slave_face->Nodes_Per_Face(); ++m) {
                    u->second_partial_t[nbDer*l+m] = u->second_partial_x[nbDer*l+m]*plane->normal.x
                                                   + u->second_partial_y[nbDer*l+m]*plane->normal.y
                                                   + u->second_partial_z[nbDer*l+m]*plane->normal.z;
                  }
                  for (m=0; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    u->second_partial_t[nbDer*l+M] = u->second_partial_x[nbDer*l+M]*plane->normal.x + u->partial_x[l]*plane_dnormal[m][0]
                                                   + u->second_partial_y[nbDer*l+M]*plane->normal.y + u->partial_y[l]*plane_dnormal[m][1]
                                                   + u->second_partial_z[nbDer*l+M]*plane->normal.z + u->partial_z[l]*plane_dnormal[m][2];
                  }
                }
                for (l=0; l<3*master_face->Nodes_Per_Face(); ++l) {
                  int L = 3*slave_face->Nodes_Per_Face()+l;
                  u->partial_t[L] = u->partial_x[L]*plane->normal.x + u->x*plane_dnormal[l][0]
                                  + u->partial_y[L]*plane->normal.y + u->y*plane_dnormal[l][1]
                                  + u->partial_z[L]*plane->normal.z + u->z*plane_dnormal[l][2]
                                  + plane_dd[l];
                  //Note: this block is entirely in the strictly lower triangular part of the matrix
                  //for (m=0; m<3*slave_face->Nodes_Per_Face(); ++m) {
                  //  u->second_partial_t[nbDer*L+m] = u->second_partial_x[nbDer*L+m]*plane->normal.x + u->partial_x[m]*plane_dnormal[l][0]
                  //                                 + u->second_partial_y[nbDer*L+m]*plane->normal.y + u->partial_y[m]*plane_dnormal[l][1]
                  //                                 + u->second_partial_z[nbDer*L+m]*plane->normal.z + u->partial_z[m]*plane_dnormal[l][2];
                  //}
                  for (m=l; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    u->second_partial_t[nbDer*L+M] = u->second_partial_x[nbDer*L+M]*plane->normal.x + u->partial_x[L]*plane_dnormal[m][0]
                                                   + u->partial_x[M]*plane_dnormal[l][0] + u->x*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][0]
                                                   + u->second_partial_y[nbDer*L+M]*plane->normal.y + u->partial_y[L]*plane_dnormal[m][1]
                                                   + u->partial_y[M]*plane_dnormal[l][1] + u->y*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][1]
                                                   + u->second_partial_z[nbDer*L+M]*plane->normal.z + u->partial_z[L]*plane_dnormal[m][2]
                                                   + u->partial_z[M]*plane_dnormal[l][2] + u->z*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][2]
                                                   + plane_d2d[3*master_face->Nodes_Per_Face()*l+m];
                  }
                }
              }
              else {
                for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                  for (m=l; m<3*slave_face->Nodes_Per_Face(); ++m) {
                    u->second_partial_t[nbDer*l+m] = 0;
                  }
                  for (m=0; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    u->second_partial_t[nbDer*l+M] = u->partial_x[l]*plane_dnormal[m][0]
                                                   + u->partial_y[l]*plane_dnormal[m][1]
                                                   + u->partial_z[l]*plane_dnormal[m][2];
                  }
                }
                for (l=0; l<3*master_face->Nodes_Per_Face(); ++l) {
                  int L = 3*slave_face->Nodes_Per_Face()+l;
                  u->partial_t[L] = u->partial_x[L]*plane->normal.x + u->x*plane_dnormal[l][0]
                                  + u->partial_y[L]*plane->normal.y + u->y*plane_dnormal[l][1]
                                  + u->partial_z[L]*plane->normal.z + u->z*plane_dnormal[l][2]
                                  + plane_dd[l];
                  //Note: this block is entirely in the strictly lower triangular part of the matrix
                  //for (m=0; m<3*slave_face->Nodes_Per_Face(); ++m) {
                  //  u->second_partial_t[nbDer*L+m] = u->partial_x[m]*plane_dnormal[l][0]
                  //                                 + u->partial_y[m]*plane_dnormal[l][1]
                  //                                 + u->partial_z[m]*plane_dnormal[l][2];
                  //}
                  for (m=l; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    u->second_partial_t[nbDer*L+M] = u->partial_x[L]*plane_dnormal[m][0]
                                                   + u->partial_x[M]*plane_dnormal[l][0] + u->x*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][0]
                                                   + u->partial_y[L]*plane_dnormal[m][1]
                                                   + u->partial_y[M]*plane_dnormal[l][1] + u->y*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][1]
                                                   + u->partial_z[L]*plane_dnormal[m][2]
                                                   + u->partial_z[M]*plane_dnormal[l][2] + u->z*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][2]
                                                   + plane_d2d[3*master_face->Nodes_Per_Face()*l+m];
                  }
                }
              }
              // fill in the strictly lower triangular part
              //for (l=0; l<nbDer; ++l)
              //  for (m=0; m<l; ++m) u->second_partial_t[nbDer*l+m] = u->second_partial_t[nbDer*m+l];
            }
            if(!v->partial_t) {
              // compute derivatives of tv w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              // due to symmetry of the second derivatives, only the upper triangular part is computed
              v->partial_t = new DataType[nbDer];
              for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                v->partial_t[l] = v->partial_x[l]*plane->normal.x + v->partial_y[l]*plane->normal.y + v->partial_z[l]*plane->normal.z;
              }
              v->second_partial_t = new DataType[nbDer*nbDer];
              if(v->second_partial_x != NULL) {
                for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                  for (m=l; m<3*slave_face->Nodes_Per_Face(); ++m) {
                    v->second_partial_t[nbDer*l+m] = v->second_partial_x[nbDer*l+m]*plane->normal.x
                                                   + v->second_partial_y[nbDer*l+m]*plane->normal.y
                                                   + v->second_partial_z[nbDer*l+m]*plane->normal.z;
                  }
                  for (m=0; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    v->second_partial_t[nbDer*l+M] = v->second_partial_x[nbDer*l+M]*plane->normal.x + v->partial_x[l]*plane_dnormal[m][0]
                                                   + v->second_partial_y[nbDer*l+M]*plane->normal.y + v->partial_y[l]*plane_dnormal[m][1]
                                                   + v->second_partial_z[nbDer*l+M]*plane->normal.z + v->partial_z[l]*plane_dnormal[m][2];
                  }
                }
                for (l=0; l<3*master_face->Nodes_Per_Face(); ++l) {
                  int L = 3*slave_face->Nodes_Per_Face()+l;
                  v->partial_t[L] = v->partial_x[L]*plane->normal.x + v->x*plane_dnormal[l][0]
                                  + v->partial_y[L]*plane->normal.y + v->y*plane_dnormal[l][1]
                                  + v->partial_z[L]*plane->normal.z + v->z*plane_dnormal[l][2]
                                  + plane_dd[l];
                  //Note: this block is entirely in the strictly lower triangular part of the matrix
                  //for (m=0; m<3*slave_face->Nodes_Per_Face(); ++m) {
                  //  v->second_partial_t[nbDer*L+m] = v->second_partial_x[nbDer*L+m]*plane->normal.x + v->partial_x[m]*plane_dnormal[l][0]
                  //                                 + v->second_partial_y[nbDer*L+m]*plane->normal.y + v->partial_y[m]*plane_dnormal[l][1]
                  //                                 + v->second_partial_z[nbDer*L+m]*plane->normal.z + v->partial_z[m]*plane_dnormal[l][2];
                  //}
                  for (m=l; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    v->second_partial_t[nbDer*L+M] = v->second_partial_x[nbDer*L+M]*plane->normal.x + v->partial_x[L]*plane_dnormal[m][0]
                                                   + v->partial_x[M]*plane_dnormal[l][0] + v->x*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][0]
                                                   + v->second_partial_y[nbDer*L+M]*plane->normal.y + v->partial_y[L]*plane_dnormal[m][1]
                                                   + v->partial_y[M]*plane_dnormal[l][1] + v->y*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][1]
                                                   + v->second_partial_z[nbDer*L+M]*plane->normal.z + v->partial_z[L]*plane_dnormal[m][2] 
                                                   + v->partial_z[M]*plane_dnormal[l][2] + v->z*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][2]
                                                   + plane_d2d[3*master_face->Nodes_Per_Face()*l+m];
                  }
                }
              }
              else {
                for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                  for (m=l; m<3*slave_face->Nodes_Per_Face(); ++m) {
                    v->second_partial_t[nbDer*l+m] = 0;
                  }
                  for (m=0; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    v->second_partial_t[nbDer*l+M] = v->partial_x[l]*plane_dnormal[m][0]
                                                   + v->partial_y[l]*plane_dnormal[m][1]
                                                   + v->partial_z[l]*plane_dnormal[m][2];
                  }
                }
                for (l=0; l<3*master_face->Nodes_Per_Face(); ++l) {
                  int L = 3*slave_face->Nodes_Per_Face()+l;
                  v->partial_t[L] = v->partial_x[L]*plane->normal.x + v->x*plane_dnormal[l][0]
                                  + v->partial_y[L]*plane->normal.y + v->y*plane_dnormal[l][1]
                                  + v->partial_z[L]*plane->normal.z + v->z*plane_dnormal[l][2]
                                  + plane_dd[l];
                  //Note: this block is entirely in the strictly lower triangular part of the matrix
                  //for (m=0; m<3*slave_face->Nodes_Per_Face(); ++m) {
                  //  v->second_partial_t[nbDer*L+m] = v->partial_x[m]*plane_dnormal[l][0]
                  //                                 + v->partial_y[m]*plane_dnormal[l][1]
                  //                                 + v->partial_z[m]*plane_dnormal[l][2];
                  //}
                  for (m=l; m<3*master_face->Nodes_Per_Face(); ++m) {
                    int M = 3*slave_face->Nodes_Per_Face()+m;
                    v->second_partial_t[nbDer*L+M] = v->partial_x[L]*plane_dnormal[m][0]
                                                   + v->partial_x[M]*plane_dnormal[l][0] + v->x*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][0]
                                                   + v->partial_y[L]*plane_dnormal[m][1]
                                                   + v->partial_y[M]*plane_dnormal[l][1] + v->y*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][1]
                                                   + v->partial_z[L]*plane_dnormal[m][2]
                                                   + v->partial_z[M]*plane_dnormal[l][2] + v->z*plane_d2normal[3*master_face->Nodes_Per_Face()*l+m][2]
                                                   + plane_d2d[3*master_face->Nodes_Per_Face()*l+m];
                  }
                }
              }
              // fill in the strictly lower triangular part
              //for (l=0; l<nbDer; ++l)
              //  for (m=0; m<l; ++m) v->second_partial_t[nbDer*l+m] = v->second_partial_t[nbDer*m+l];
            }
            w = &(q->p[q->np]);
            if(!w->partial_x) {
              w->partial_x = new DataType[nbDer];
              w->partial_y = new DataType[nbDer];
              w->partial_z = new DataType[nbDer];
            }
            if(!w->second_partial_x) {
              w->second_partial_x = new DataType[nbDer*nbDer];
              w->second_partial_y = new DataType[nbDer*nbDer];
              w->second_partial_z = new DataType[nbDer*nbDer];
            }
            if(abs(tu) < abs(tv)) {
              t    = tu/(tu-tv);
              w->x = u->x+t*(v->x-u->x);
              w->y = u->y+t*(v->y-u->y);
              w->z = u->z+t*(v->z-u->z);
              // compute derivatives of w->x, w->y, and w->z w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              // due to symmetry of the second derivatives, only the upper triangular part is computed
              DataType s = (tu-tv)*(tu-tv);
              DataType s2 = s*s;
              DataType *partial_t = new DataType[nbDer];
              DataType *partial_s = new DataType[nbDer];
              for (l=0; l<nbDer; ++l) {
                partial_t[l] = (tu*v->partial_t[l] - tv*u->partial_t[l])/s;
                partial_s[l] = 2*(u->partial_t[l]-v->partial_t[l])*(tu-tv);
              }
              for (l=0; l<nbDer; ++l) {
                w->partial_x[l] = u->partial_x[l]+t*(v->partial_x[l]-u->partial_x[l]) + partial_t[l]*(v->x-u->x);
                w->partial_y[l] = u->partial_y[l]+t*(v->partial_y[l]-u->partial_y[l]) + partial_t[l]*(v->y-u->y);
                w->partial_z[l] = u->partial_z[l]+t*(v->partial_z[l]-u->partial_z[l]) + partial_t[l]*(v->z-u->z);
              }
              if(u->second_partial_x != NULL && v->second_partial_x != NULL) {
                for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((u->partial_t[m]*v->partial_t[l]+tu*v->second_partial_t[nbDer*l+m]
                                                - v->partial_t[m]*u->partial_t[l]-tv*u->second_partial_t[nbDer*l+m])*s
                                                - (tu*v->partial_t[l] - tv*u->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = u->second_partial_x[nbDer*l+m]
                                           + t*(v->second_partial_x[nbDer*l+m]-u->second_partial_x[nbDer*l+m]) + partial_t[m]*(v->partial_x[l]-u->partial_x[l])
                                           + second_partial_t*(v->x-u->x) + partial_t[l]*(v->partial_x[m]-u->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = u->second_partial_y[nbDer*l+m]
                                           + t*(v->second_partial_y[nbDer*l+m]-u->second_partial_y[nbDer*l+m]) + partial_t[m]*(v->partial_y[l]-u->partial_y[l])
                                           + second_partial_t*(v->y-u->y) + partial_t[l]*(v->partial_y[m]-u->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = u->second_partial_z[nbDer*l+m]
                                           + t*(v->second_partial_z[nbDer*l+m]-u->second_partial_z[nbDer*l+m]) + partial_t[m]*(v->partial_z[l]-u->partial_z[l])
                                           + second_partial_t*(v->z-u->z) + partial_t[l]*(v->partial_z[m]-u->partial_z[m]);
                  }
                }
              }
              else if(u->second_partial_x != NULL && v->second_partial_x == NULL) {
                for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((u->partial_t[m]*v->partial_t[l]+tu*v->second_partial_t[nbDer*l+m]
                                                - v->partial_t[m]*u->partial_t[l]-tv*u->second_partial_t[nbDer*l+m])*s
                                                - (tu*v->partial_t[l] - tv*u->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = u->second_partial_x[nbDer*l+m]
                                           - t*u->second_partial_x[nbDer*l+m] + partial_t[m]*(v->partial_x[l]-u->partial_x[l])
                                           + second_partial_t*(v->x-u->x) + partial_t[l]*(v->partial_x[m]-u->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = u->second_partial_y[nbDer*l+m]
                                           - t*u->second_partial_y[nbDer*l+m] + partial_t[m]*(v->partial_y[l]-u->partial_y[l])
                                           + second_partial_t*(v->y-u->y) + partial_t[l]*(v->partial_y[m]-u->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = u->second_partial_z[nbDer*l+m]
                                           - t*u->second_partial_z[nbDer*l+m] + partial_t[m]*(v->partial_z[l]-u->partial_z[l])
                                           + second_partial_t*(v->z-u->z) + partial_t[l]*(v->partial_z[m]-u->partial_z[m]);
                  }
                }
              }
              else if(u->second_partial_x == NULL && v->second_partial_x != NULL) {
                for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((u->partial_t[m]*v->partial_t[l]+tu*v->second_partial_t[nbDer*l+m]
                                                - v->partial_t[m]*u->partial_t[l]-tv*u->second_partial_t[nbDer*l+m])*s
                                                - (tu*v->partial_t[l] - tv*u->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = t*v->second_partial_x[nbDer*l+m] + partial_t[m]*(v->partial_x[l]-u->partial_x[l])
                                           + second_partial_t*(v->x-u->x) + partial_t[l]*(v->partial_x[m]-u->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = t*v->second_partial_y[nbDer*l+m] + partial_t[m]*(v->partial_y[l]-u->partial_y[l])
                                           + second_partial_t*(v->y-u->y) + partial_t[l]*(v->partial_y[m]-u->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = t*v->second_partial_z[nbDer*l+m] + partial_t[m]*(v->partial_z[l]-u->partial_z[l])
                                           + second_partial_t*(v->z-u->z) + partial_t[l]*(v->partial_z[m]-u->partial_z[m]);
                  }
                }
              }
              else {
                for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((u->partial_t[m]*v->partial_t[l]+tu*v->second_partial_t[nbDer*l+m]
                                                - v->partial_t[m]*u->partial_t[l]-tv*u->second_partial_t[nbDer*l+m])*s
                                                - (tu*v->partial_t[l] - tv*u->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = partial_t[m]*(v->partial_x[l]-u->partial_x[l])
                                           + second_partial_t*(v->x-u->x) + partial_t[l]*(v->partial_x[m]-u->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = partial_t[m]*(v->partial_y[l]-u->partial_y[l])
                                           + second_partial_t*(v->y-u->y) + partial_t[l]*(v->partial_y[m]-u->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = partial_t[m]*(v->partial_z[l]-u->partial_z[l])
                                           + second_partial_t*(v->z-u->z) + partial_t[l]*(v->partial_z[m]-u->partial_z[m]);
                  }
                }
              }
              delete [] partial_t;
              delete [] partial_s;
            }
            else {
              t    = tv/(tv-tu);
              w->x = v->x+t*(u->x-v->x);
              w->y = v->y+t*(u->y-v->y);
              w->z = v->z+t*(u->z-v->z);
              // compute derivatives of w->x, w->y, and w->z w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              // due to symmetry of the second derivatives, only the upper triangular part is computed
              DataType s = (tv-tu)*(tv-tu);
              DataType s2 = s*s;
              DataType *partial_t = new DataType[nbDer];
              DataType *partial_s = new DataType[nbDer];
              for (l=0; l<nbDer; ++l) {
                partial_t[l] = (tv*u->partial_t[l] - tu*v->partial_t[l])/s;
                partial_s[l] = 2*(v->partial_t[l]-u->partial_t[l])*(tv-tu);
              }
              for (l=0; l<nbDer; ++l) {
                w->partial_x[l] = v->partial_x[l]+t*(u->partial_x[l]-v->partial_x[l]) + partial_t[l]*(u->x-v->x);
                w->partial_y[l] = v->partial_y[l]+t*(u->partial_y[l]-v->partial_y[l]) + partial_t[l]*(u->y-v->y);
                w->partial_z[l] = v->partial_z[l]+t*(u->partial_z[l]-v->partial_z[l]) + partial_t[l]*(u->z-v->z);
              }
              if(u->second_partial_x != NULL && v->second_partial_x != NULL) {
                for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((v->partial_t[m]*u->partial_t[l]+tv*u->second_partial_t[nbDer*l+m]
                                                - u->partial_t[m]*v->partial_t[l]-tu*v->second_partial_t[nbDer*l+m])*s
                                                - (tv*u->partial_t[l] - tu*v->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = v->second_partial_x[nbDer*l+m]
                                           + t*(u->second_partial_x[nbDer*l+m]-v->second_partial_x[nbDer*l+m]) + partial_t[m]*(u->partial_x[l]-v->partial_x[l])
                                           + second_partial_t*(u->x-v->x) + partial_t[l]*(u->partial_x[m]-v->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = v->second_partial_y[nbDer*l+m]
                                           + t*(u->second_partial_y[nbDer*l+m]-v->second_partial_y[nbDer*l+m]) + partial_t[m]*(u->partial_y[l]-v->partial_y[l])
                                           + second_partial_t*(u->y-v->y) + partial_t[l]*(u->partial_y[m]-v->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = v->second_partial_z[nbDer*l+m]
                                           + t*(u->second_partial_z[nbDer*l+m]-v->second_partial_z[nbDer*l+m]) + partial_t[m]*(u->partial_z[l]-v->partial_z[l])
                                           + second_partial_t*(u->z-v->z) + partial_t[l]*(u->partial_z[m]-v->partial_z[m]);
                  }
                }
              }
              else if(u->second_partial_x != NULL && v->second_partial_x == NULL) {
                for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((v->partial_t[m]*u->partial_t[l]+tv*u->second_partial_t[nbDer*l+m]
                                                - u->partial_t[m]*v->partial_t[l]-tu*v->second_partial_t[nbDer*l+m])*s
                                                - (tv*u->partial_t[l] - tu*v->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = t*u->second_partial_x[nbDer*l+m] + partial_t[m]*(u->partial_x[l]-v->partial_x[l])
                                           + second_partial_t*(u->x-v->x) + partial_t[l]*(u->partial_x[m]-v->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = t*u->second_partial_y[nbDer*l+m] + partial_t[m]*(u->partial_y[l]-v->partial_y[l])
                                           + second_partial_t*(u->y-v->y) + partial_t[l]*(u->partial_y[m]-v->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = t*u->second_partial_z[nbDer*l+m] + partial_t[m]*(u->partial_z[l]-v->partial_z[l])
                                           + second_partial_t*(u->z-v->z) + partial_t[l]*(u->partial_z[m]-v->partial_z[m]);
                  }
                }
              }
              else if(u->second_partial_x == NULL && v->second_partial_x != NULL) {
                for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((v->partial_t[m]*u->partial_t[l]+tv*u->second_partial_t[nbDer*l+m]
                                                - u->partial_t[m]*v->partial_t[l]-tu*v->second_partial_t[nbDer*l+m])*s
                                                - (tv*u->partial_t[l] - tu*v->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = v->second_partial_x[nbDer*l+m]
                                           - t*v->second_partial_x[nbDer*l+m] + partial_t[m]*(u->partial_x[l]-v->partial_x[l])
                                           + second_partial_t*(u->x-v->x) + partial_t[l]*(u->partial_x[m]-v->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = v->second_partial_y[nbDer*l+m]
                                           - t*v->second_partial_y[nbDer*l+m] + partial_t[m]*(u->partial_y[l]-v->partial_y[l])
                                           + second_partial_t*(u->y-v->y) + partial_t[l]*(u->partial_y[m]-v->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = v->second_partial_z[nbDer*l+m]
                                           - t*v->second_partial_z[nbDer*l+m] + partial_t[m]*(u->partial_z[l]-v->partial_z[l])
                                           + second_partial_t*(u->z-v->z) + partial_t[l]*(u->partial_z[m]-v->partial_z[m]);
                  }
                }
              }
              else {
               for (l=0; l<nbDer; ++l) {
                  for (m=l; m<nbDer; ++m) {
                    DataType second_partial_t = ((v->partial_t[m]*u->partial_t[l]+tv*u->second_partial_t[nbDer*l+m]
                                                - u->partial_t[m]*v->partial_t[l]-tu*v->second_partial_t[nbDer*l+m])*s
                                                - (tv*u->partial_t[l] - tu*v->partial_t[l])*partial_s[m])/s2;
                    w->second_partial_x[nbDer*l+m] = partial_t[m]*(u->partial_x[l]-v->partial_x[l])
                                           + second_partial_t*(u->x-v->x) + partial_t[l]*(u->partial_x[m]-v->partial_x[m]);
                    w->second_partial_y[nbDer*l+m] = partial_t[m]*(u->partial_y[l]-v->partial_y[l])
                                           + second_partial_t*(u->y-v->y) + partial_t[l]*(u->partial_y[m]-v->partial_y[m]);
                    w->second_partial_z[nbDer*l+m] = partial_t[m]*(u->partial_z[l]-v->partial_z[l])
                                           + second_partial_t*(u->z-v->z) + partial_t[l]*(u->partial_z[m]-v->partial_z[m]);
                  }
                }
              }
              delete [] partial_t;
              delete [] partial_s;
            }
            // fill in the strictly lower triangular part of w->second_partial_x/y/z
            //for (l=0; l<nbDer; ++l)
            //  for (m=0; m<l; ++m) {
            //    w->second_partial_x[nbDer*l+m] = w->second_partial_x[nbDer*m+l];
            //    w->second_partial_y[nbDer*l+m] = w->second_partial_y[nbDer*m+l];
            //    w->second_partial_z[nbDer*l+m] = w->second_partial_z[nbDer*m+l];
            //  }
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
            if(!w->partial_x) {
              w->partial_x = new DataType[nbDer];
              w->partial_y = new DataType[nbDer];
              w->partial_z = new DataType[nbDer];
            }
            for (l=0; l<nbDer; ++l) {
              w->partial_x[l] = v->partial_x[l];
              w->partial_y[l] = v->partial_y[l];
              w->partial_z[l] = v->partial_z[l];
            }
            if(v->second_partial_x) {
              if(!w->second_partial_x) {
                w->second_partial_x = new DataType[nbDer*nbDer];
                w->second_partial_y = new DataType[nbDer*nbDer];
                w->second_partial_z = new DataType[nbDer*nbDer];
              }
              for (l=0; l<nbDer*nbDer; ++l) {
                w->second_partial_x[l] = v->second_partial_x[l];
                w->second_partial_y[l] = v->second_partial_y[l];
                w->second_partial_z[l] = v->second_partial_z[l];
              }
            }
            else {
              if(w->second_partial_x) {
                delete [] w->second_partial_x; w->second_partial_x = NULL;
                delete [] w->second_partial_y; w->second_partial_y = NULL;
                delete [] w->second_partial_z; w->second_partial_z = NULL;
              }
            }
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
              if(!w->partial_x) {
                w->partial_x = new DataType[nbDer];
                w->partial_y = new DataType[nbDer];
                w->partial_z = new DataType[nbDer];
              }
              for (l=0; l<nbDer; ++l) {
                w->partial_x[l] = v->partial_x[l];
                w->partial_y[l] = v->partial_y[l];
                w->partial_z[l] = v->partial_z[l];
              }
              if(v->second_partial_x) {
                if(!w->second_partial_x) {
                  w->second_partial_x = new DataType[nbDer*nbDer];
                  w->second_partial_y = new DataType[nbDer*nbDer];
                  w->second_partial_z = new DataType[nbDer*nbDer];
                }
                for (l=0; l<nbDer*nbDer; ++l) {
                  w->second_partial_x[l] = v->second_partial_x[l];
                  w->second_partial_y[l] = v->second_partial_y[l];
                  w->second_partial_z[l] = v->second_partial_z[l];
                }
              }
              else {
                if(w->second_partial_x) {
                  delete [] w->second_partial_x; w->second_partial_x = NULL;
                  delete [] w->second_partial_y; w->second_partial_y = NULL;
                  delete [] w->second_partial_z; w->second_partial_z = NULL;
                }
              }
              q->np++;
            }
          }
        }
        r = p;
        p = q;
        q = r;
        if(plane_dnormal) delete [] plane_dnormal;
        if(plane_dd) delete [] plane_dd;
        if(plane_d2normal) delete [] plane_d2normal;
        if(plane_d2d) delete [] plane_d2d;
        if (p->np<3) {
          p->np = 0;
          break;
        }
      }
      if(master_face_dnormal) delete [] master_face_dnormal;
      if(master_face_d2normal) delete [] master_face_d2normal;
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
        DataType (*face_dlocal_coords)[2] = new DataType[3*face->Nodes_Per_Face()][2];
        DataType (*master_face_dlocal_coords)[2] = new DataType[3*master_face->Nodes_Per_Face()][2];
        DataType dmdX[2], dmdY[2], dmdZ[2];
        DataType (*face_d2local_coords)[2] = new DataType[3*face->Nodes_Per_Face()*3*face->Nodes_Per_Face()][2];
        DataType (*master_face_d2local_coords)[2] = new DataType[3*master_face->Nodes_Per_Face()*3*master_face->Nodes_Per_Face()][2];
        DataType d2mdX2[2], d2mdY2[2], d2mdZ2[2], d2mdXdY[2], d2mdYdZ[2], d2mdXdZ[2];
        DataType (*face_ddlocal_coords_dX)[2] = new DataType[3*face->Nodes_Per_Face()][2],
                 (*face_ddlocal_coords_dY)[2] = new DataType[3*face->Nodes_Per_Face()][2],
                 (*face_ddlocal_coords_dZ)[2] = new DataType[3*face->Nodes_Per_Face()][2];
        DataType (*master_face_ddlocal_coords_dX)[2] = new DataType[3*master_face->Nodes_Per_Face()][2],
                 (*master_face_ddlocal_coords_dY)[2] = new DataType[3*master_face->Nodes_Per_Face()][2],
                 (*master_face_ddlocal_coords_dZ)[2] = new DataType[3*master_face->Nodes_Per_Face()][2];
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

          // compute partial_area_s[2*i+0], partial_area_s[2*i+1] 
          // and second_partial_area_s[[2*i+0], second_partial_area_s[2*i+1] 
          // due to symmetry of the second derivatives, only the upper triangular part is computed
          face->Compute_Partial_Local_Coordinates_1(POSITION, global_coords, face_dlocal_coords);
          face->Compute_Partial_Local_Coordinates_2(POSITION, global_coords, dmdX, dmdY, dmdZ);
          face->Compute_Second_Partial_Local_Coordinates_1(POSITION, global_coords, face_d2local_coords);
          face->Compute_Second_Partial_Local_Coordinates_2(POSITION, global_coords, d2mdX2, d2mdY2, d2mdZ2, d2mdXdY, d2mdYdZ, d2mdXdZ);
          face->Compute_Second_Partial_Local_Coordinates_12(POSITION, global_coords, face_ddlocal_coords_dX, face_ddlocal_coords_dY,
                                                            face_ddlocal_coords_dZ);
          for (j=0; j<nbDer; ++j) {
            partial_area_s[j][2*i+0] = dmdX[0]*p->p[i].partial_x[j] + dmdY[0]*p->p[i].partial_y[j] + dmdZ[0]*p->p[i].partial_z[j];
            partial_area_s[j][2*i+1] = dmdX[1]*p->p[i].partial_x[j] + dmdY[1]*p->p[i].partial_y[j] + dmdZ[1]*p->p[i].partial_z[j];
            for (k=j; k<nbDer; ++k) { // add the term d2area_s/dglobal_coords2*(dglobal_coords/dQ)^2
              second_partial_area_s[nbDer*j+k][2*i+0] = d2mdX2[0]*p->p[i].partial_x[k]*p->p[i].partial_x[j]
                                               + d2mdY2[0]*p->p[i].partial_y[k]*p->p[i].partial_y[j] + d2mdZ2[0]*p->p[i].partial_z[k]*p->p[i].partial_z[j]
                                               + d2mdXdY[0]*(p->p[i].partial_y[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_y[j])
                                               + d2mdYdZ[0]*(p->p[i].partial_z[k]*p->p[i].partial_y[j] + p->p[i].partial_y[k]*p->p[i].partial_z[j])
                                               + d2mdXdZ[0]*(p->p[i].partial_z[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_z[j]);
              second_partial_area_s[nbDer*j+k][2*i+1] = d2mdX2[1]*p->p[i].partial_x[k]*p->p[i].partial_x[j]
                                               + d2mdY2[1]*p->p[i].partial_y[k]*p->p[i].partial_y[j] + d2mdZ2[1]*p->p[i].partial_z[k]*p->p[i].partial_z[j]
                                               + d2mdXdY[1]*(p->p[i].partial_y[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_y[j])
                                               + d2mdYdZ[1]*(p->p[i].partial_z[k]*p->p[i].partial_y[j] + p->p[i].partial_y[k]*p->p[i].partial_z[j])
                                               + d2mdXdZ[1]*(p->p[i].partial_z[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_z[j]);
            }
            if(p->p[i].second_partial_x != NULL) {
              for (k=j; k<nbDer; ++k) { // add the term darea_s/dglobal_coords*d2global_coords/dQ2
                second_partial_area_s[nbDer*j+k][2*i+0] += dmdX[0]*p->p[i].second_partial_x[nbDer*j+k] + dmdY[0]*p->p[i].second_partial_y[nbDer*j+k]
                                                         + dmdZ[0]*p->p[i].second_partial_z[nbDer*j+k];
                second_partial_area_s[nbDer*j+k][2*i+1] += dmdX[1]*p->p[i].second_partial_x[nbDer*j+k] + dmdY[1]*p->p[i].second_partial_y[nbDer*j+k]
                                                         + dmdZ[1]*p->p[i].second_partial_z[nbDer*j+k];
              }
            }
            for (k=j; k<3*face->Nodes_Per_Face(); ++k) { // add the term d2area_s/dglobal_coords/dQ*dglobal_coords/dQ
              second_partial_area_s[nbDer*j+k][2*i+0] += face_ddlocal_coords_dX[k][0]*p->p[i].partial_x[j]
                                                       + face_ddlocal_coords_dY[k][0]*p->p[i].partial_y[j]
                                                       + face_ddlocal_coords_dZ[k][0]*p->p[i].partial_z[j];
              second_partial_area_s[nbDer*j+k][2*i+1] += face_ddlocal_coords_dX[k][1]*p->p[i].partial_x[j]
                                                       + face_ddlocal_coords_dY[k][1]*p->p[i].partial_y[j]
                                                       + face_ddlocal_coords_dZ[k][1]*p->p[i].partial_z[j];
            }
          }
          for (j=0; j<3*face->Nodes_Per_Face(); ++j) {
            partial_area_s[j][2*i+0] += face_dlocal_coords[j][0];
            partial_area_s[j][2*i+1] += face_dlocal_coords[j][1];

            for (k=j; k<nbDer; ++k) { // add the term d2area_s/dglobal_coords/dQ*dglobal_coords/dQ
              second_partial_area_s[nbDer*j+k][2*i+0] += face_ddlocal_coords_dX[j][0]*p->p[i].partial_x[k]
                                                       + face_ddlocal_coords_dY[j][0]*p->p[i].partial_y[k]
                                                       + face_ddlocal_coords_dZ[j][0]*p->p[i].partial_z[k];
              second_partial_area_s[nbDer*j+k][2*i+1] += face_ddlocal_coords_dX[j][1]*p->p[i].partial_x[k]
                                                       + face_ddlocal_coords_dY[j][1]*p->p[i].partial_y[k]
                                                       + face_ddlocal_coords_dZ[j][1]*p->p[i].partial_z[k];
            }

            for (k=j; k<3*face->Nodes_Per_Face(); ++k) { // add the term: d2area_s/dQ2
              second_partial_area_s[nbDer*j+k][2*i+0] += face_d2local_coords[3*face->Nodes_Per_Face()*j+k][0];
              second_partial_area_s[nbDer*j+k][2*i+1] += face_d2local_coords[3*face->Nodes_Per_Face()*j+k][1];
            }
          }
          // compute partial_area_m[2*i+0], partial_area_m[2*i+1]
          // and second_partial_area_m[[2*i+0], second_partial_area_m[2*i+1]
          master_face->Compute_Partial_Local_Coordinates_1(POSITION, global_coords, master_face_dlocal_coords);
          master_face->Compute_Partial_Local_Coordinates_2(POSITION, global_coords, dmdX, dmdY, dmdZ);
          master_face->Compute_Second_Partial_Local_Coordinates_1(POSITION, global_coords, master_face_d2local_coords);
          master_face->Compute_Second_Partial_Local_Coordinates_2(POSITION, global_coords, d2mdX2, d2mdY2, d2mdZ2, d2mdXdY, d2mdYdZ, d2mdXdZ);
          master_face->Compute_Second_Partial_Local_Coordinates_12(POSITION, global_coords, master_face_ddlocal_coords_dX,
                                                                   master_face_ddlocal_coords_dY, master_face_ddlocal_coords_dZ);

          for (j=0; j<nbDer; ++j) {
            partial_area_m[j][2*i+0] = dmdX[0]*p->p[i].partial_x[j] + dmdY[0]*p->p[i].partial_y[j] + dmdZ[0]*p->p[i].partial_z[j];
            partial_area_m[j][2*i+1] = dmdX[1]*p->p[i].partial_x[j] + dmdY[1]*p->p[i].partial_y[j] + dmdZ[1]*p->p[i].partial_z[j];
            for (k=j; k<nbDer; ++k) { // add the term d2area_m/dglobal_coords2*(dglobal_coords/dQ)^2
              second_partial_area_m[nbDer*j+k][2*i+0] = d2mdX2[0]*p->p[i].partial_x[k]*p->p[i].partial_x[j]
                                               + d2mdY2[0]*p->p[i].partial_y[k]*p->p[i].partial_y[j] + d2mdZ2[0]*p->p[i].partial_z[k]*p->p[i].partial_z[j]
                                               + d2mdXdY[0]*(p->p[i].partial_y[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_y[j])
                                               + d2mdYdZ[0]*(p->p[i].partial_z[k]*p->p[i].partial_y[j] + p->p[i].partial_y[k]*p->p[i].partial_z[j])
                                               + d2mdXdZ[0]*(p->p[i].partial_z[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_z[j]);
              second_partial_area_m[nbDer*j+k][2*i+1] = d2mdX2[1]*p->p[i].partial_x[k]*p->p[i].partial_x[j]
                                               + d2mdY2[1]*p->p[i].partial_y[k]*p->p[i].partial_y[j] + d2mdZ2[1]*p->p[i].partial_z[k]*p->p[i].partial_z[j]
                                               + d2mdXdY[1]*(p->p[i].partial_y[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_y[j])
                                               + d2mdYdZ[1]*(p->p[i].partial_z[k]*p->p[i].partial_y[j] + p->p[i].partial_y[k]*p->p[i].partial_z[j])
                                               + d2mdXdZ[1]*(p->p[i].partial_z[k]*p->p[i].partial_x[j] + p->p[i].partial_x[k]*p->p[i].partial_z[j]);
            }
            if(p->p[i].second_partial_x != NULL) {
              for (k=j; k<nbDer; ++k) { // add the term darea_m/dglobal_coords*d2global_coords/dQ2
                second_partial_area_m[nbDer*j+k][2*i+0] += dmdX[0]*p->p[i].second_partial_x[nbDer*j+k] + dmdY[0]*p->p[i].second_partial_y[nbDer*j+k]
                                                         + dmdZ[0]*p->p[i].second_partial_z[nbDer*j+k];
                second_partial_area_m[nbDer*j+k][2*i+1] += dmdX[1]*p->p[i].second_partial_x[nbDer*j+k] + dmdY[1]*p->p[i].second_partial_y[nbDer*j+k]
                                                         + dmdZ[1]*p->p[i].second_partial_z[nbDer*j+k];
              }
            }
            for (k=std::max(0,j-3*face->Nodes_Per_Face()); k<3*master_face->Nodes_Per_Face(); ++k) { // add the term d2area_m/dglobal_coords/dQ*dglobal_coords/dQ
              int K = 3*face->Nodes_Per_Face()+k;
              second_partial_area_m[nbDer*j+K][2*i+0] += master_face_ddlocal_coords_dX[k][0]*p->p[i].partial_x[j]
                                                       + master_face_ddlocal_coords_dY[k][0]*p->p[i].partial_y[j]
                                                       + master_face_ddlocal_coords_dZ[k][0]*p->p[i].partial_z[j];
              second_partial_area_m[nbDer*j+K][2*i+1] += master_face_ddlocal_coords_dX[k][1]*p->p[i].partial_x[j]
                                                       + master_face_ddlocal_coords_dY[k][1]*p->p[i].partial_y[j]
                                                       + master_face_ddlocal_coords_dZ[k][1]*p->p[i].partial_z[j];
            }
          }
          for (j=0; j<3*master_face->Nodes_Per_Face(); ++j) {
            int J = 3*face->Nodes_Per_Face()+j;
            partial_area_m[J][2*i+0] += master_face_dlocal_coords[j][0];
            partial_area_m[J][2*i+1] += master_face_dlocal_coords[j][1];

            for (k=3*face->Nodes_Per_Face()+j; k<nbDer; ++k) { // add the term d2area_m/dglobal_coords/dQ*dglobal_coords/dQ
              second_partial_area_m[nbDer*J+k][2*i+0] += master_face_ddlocal_coords_dX[j][0]*p->p[i].partial_x[k]
                                                       + master_face_ddlocal_coords_dY[j][0]*p->p[i].partial_y[k]
                                                       + master_face_ddlocal_coords_dZ[j][0]*p->p[i].partial_z[k];
              second_partial_area_m[nbDer*J+k][2*i+1] += master_face_ddlocal_coords_dX[j][1]*p->p[i].partial_x[k]
                                                       + master_face_ddlocal_coords_dY[j][1]*p->p[i].partial_y[k]
                                                       + master_face_ddlocal_coords_dZ[j][1]*p->p[i].partial_z[k];
            }

            for (k=j; k<3*master_face->Nodes_Per_Face(); ++k) { // add the term: d2area_m/dQ2
              int K = 3*face->Nodes_Per_Face()+k;
              second_partial_area_m[nbDer*J+K][2*i+0] += master_face_d2local_coords[3*master_face->Nodes_Per_Face()*j+k][0];
              second_partial_area_m[nbDer*J+K][2*i+1] += master_face_d2local_coords[3*master_face->Nodes_Per_Face()*j+k][1];
            }
          }
          // fill in the strictly lower triangular part
          for (j=0; j<nbDer; ++j) {
            for (k=0; k<j; ++k) {
              second_partial_area_s[nbDer*j+k][2*i+0] = second_partial_area_s[nbDer*k+j][2*i+0];
              second_partial_area_s[nbDer*j+k][2*i+1] = second_partial_area_s[nbDer*k+j][2*i+1];
              second_partial_area_m[nbDer*j+k][2*i+0] = second_partial_area_m[nbDer*k+j][2*i+0];
              second_partial_area_m[nbDer*j+k][2*i+1] = second_partial_area_m[nbDer*k+j][2*i+1];
            }
          }
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
              ifaceedge, iedge_m, area_s, area_m, partial_area_s, partial_area_m,
              second_partial_area_s, second_partial_area_m );
        delete [] face_dlocal_coords;
        delete [] master_face_dlocal_coords;
        delete [] face_d2local_coords;
        delete [] master_face_d2local_coords;
        delete [] face_ddlocal_coords_dX;
        delete [] face_ddlocal_coords_dY;
        delete [] face_ddlocal_coords_dZ;
        delete [] master_face_ddlocal_coords_dX;
        delete [] master_face_ddlocal_coords_dY;
        delete [] master_face_ddlocal_coords_dZ;
      }
    }

  } else {
    std::cerr << " *** ERROR: ContactSearch::Second_Partial_Face_Face_Search is not implemented for non-planar case\n";
    exit(-1);
  }
  delete [] partial_area_s;
  delete [] partial_area_m;
  delete [] second_partial_area_s;
  delete [] second_partial_area_m;
  return cffi;
}

