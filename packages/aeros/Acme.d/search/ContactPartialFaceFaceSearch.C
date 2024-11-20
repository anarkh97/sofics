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
struct PointStructExp {
  DataType x, y, z;
  DataType xm, ym, zm;

  DataType *partial_x, *partial_y, *partial_z;
  DataType *partial_t;
  PointStructExp() : partial_x(NULL), partial_y(NULL), partial_z(NULL), partial_t(NULL) {}
  ~PointStructExp() {
    if(partial_x) delete [] partial_x;
    if(partial_y) delete [] partial_y;
    if(partial_z) delete [] partial_z;
    if(partial_t) delete [] partial_t;
  }
};

template<typename DataType>
struct PolyStructExp {
  int np;
  PointStructExp<DataType> p[20];
};

template<typename DataType>
struct PlaneStructExp {
  PointStructExp<DataType> normal;  
  DataType d;
};

#define V3_Dot(a,b) \
        ((a)->x * (b)->x + (a)->y * (b)->y + (a)->z * (b)->z)

template<typename DataType>
ContactFaceFaceInteraction<DataType>*
ContactSearch::Partial_Face_Face_Search(ContactFace<DataType>* slave_face, 
                                        ContactFace<DataType>* master_face, 
                                        ContactElem<DataType>* element,
                                        VariableHandle POSITION, Real tol,
                                        ContactFixedSizeAllocator *allocators)
{
  ContactFaceFaceInteraction<DataType>* cffi = NULL;
  typedef PointStructExp<DataType> Point;
  typedef PolyStructExp<DataType> Poly;
  typedef PlaneStructExp<DataType> Plane;

  using std::abs;
  using std::sqrt;

  const int nbDer = 3*(slave_face->Nodes_Per_Face() + master_face->Nodes_Per_Face());

  int   i, j, k, l, in_cnt;
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
        //========================
        // START WITH U=VERT[N-1]
        //========================
        q->np = 0;
        u     = &(p->p[p->np-1]);
        tu    = V3_Dot(u, &plane->normal)+plane->d;
        for (j=0; j<20; ++j) {
          if(p->p[j].partial_t) { delete [] p->p[j].partial_t; p->p[j].partial_t = NULL; }
          if(q->p[j].partial_t) { delete [] q->p[j].partial_t; q->p[j].partial_t = NULL; }
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
              }
              // compute derivatives of plane->normal and plane->d w.r.t. global X-, Y- and Z-coordinates of master_face nodes
              plane_dnormal = new DataType[3*master_face->Nodes_Per_Face()][3];
              plane_dd = new DataType[3*master_face->Nodes_Per_Face()];
              element->Compute_Partial_Face_Normal(i, POSITION, FACE_NORMAL, master_face_dnormal, tol, plane_dnormal, plane_dd);
            }
            if(!u->partial_t) {
              // compute derivatives of tu w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              u->partial_t = new DataType[nbDer];
              for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                u->partial_t[l] = u->partial_x[l]*plane->normal.x + u->partial_y[l]*plane->normal.y + u->partial_z[l]*plane->normal.z;
              }
              for (l=0; l<3*master_face->Nodes_Per_Face(); ++l) {
                int L = 3*slave_face->Nodes_Per_Face()+l;
                u->partial_t[L] = u->partial_x[L]*plane->normal.x + u->x*plane_dnormal[l][0]
                                + u->partial_y[L]*plane->normal.y + u->y*plane_dnormal[l][1]
                                + u->partial_z[L]*plane->normal.z + u->z*plane_dnormal[l][2]
                                + plane_dd[l];
              }
            }
            if(!v->partial_t) {
              // compute derivatives of tv w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              v->partial_t = new DataType[nbDer];
              for (l=0; l<3*slave_face->Nodes_Per_Face(); ++l) {
                v->partial_t[l] = v->partial_x[l]*plane->normal.x + v->partial_y[l]*plane->normal.y + v->partial_z[l]*plane->normal.z;
              }
              for (l=0; l<3*master_face->Nodes_Per_Face(); ++l) {
                int L = 3*slave_face->Nodes_Per_Face()+l;
                v->partial_t[L] = v->partial_x[L]*plane->normal.x + v->x*plane_dnormal[l][0]
                                + v->partial_y[L]*plane->normal.y + v->y*plane_dnormal[l][1]
                                + v->partial_z[L]*plane->normal.z + v->z*plane_dnormal[l][2]
                                + plane_dd[l];
              }
            }
            w = &(q->p[q->np]);
            if(!w->partial_x) {
              w->partial_x = new DataType[nbDer];
              w->partial_y = new DataType[nbDer];
              w->partial_z = new DataType[nbDer];
            }
            if(abs(tu) < abs(tv)) {
              t    = tu/(tu-tv);
              w->x = u->x+t*(v->x-u->x);
              w->y = u->y+t*(v->y-u->y);
              w->z = u->z+t*(v->z-u->z);
              // compute derivatives of w->x, w->y, and w->z w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              DataType s = (tu-tv)*(tu-tv);
              for (l=0; l<nbDer; ++l) {
                DataType partial_t = (tu*v->partial_t[l] - tv*u->partial_t[l])/s;
                w->partial_x[l] = u->partial_x[l]+t*(v->partial_x[l]-u->partial_x[l]) + partial_t*(v->x-u->x);
                w->partial_y[l] = u->partial_y[l]+t*(v->partial_y[l]-u->partial_y[l]) + partial_t*(v->y-u->y);
                w->partial_z[l] = u->partial_z[l]+t*(v->partial_z[l]-u->partial_z[l]) + partial_t*(v->z-u->z);
              }
            }
            else {
              t    = tv/(tv-tu);
              w->x = v->x+t*(u->x-v->x);
              w->y = v->y+t*(u->y-v->y);
              w->z = v->z+t*(u->z-v->z);
              // compute derivatives of w->x, w->y, and w->z w.r.t. global X-, Y- and Z-coordinates of slave_face and master_face nodes
              DataType s = (tv-tu)*(tv-tu);
              for (l=0; l<nbDer; ++l) {
                DataType partial_t = (tv*u->partial_t[l] - tu*v->partial_t[l])/s;
                w->partial_x[l] = v->partial_x[l]+t*(u->partial_x[l]-v->partial_x[l]) + partial_t*(u->x-v->x);
                w->partial_y[l] = v->partial_y[l]+t*(u->partial_y[l]-v->partial_y[l]) + partial_t*(u->y-v->y);
                w->partial_z[l] = v->partial_z[l]+t*(u->partial_z[l]-v->partial_z[l]) + partial_t*(u->z-v->z);
              }
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
              q->np++;
            }
          }
        }
        r = p;
        p = q;
        q = r;
        if(plane_dnormal) delete [] plane_dnormal;
        if(plane_dd) delete [] plane_dd;
        if (p->np<3) {
          p->np = 0;
          break;
        }
      }
      if(master_face_dnormal) delete [] master_face_dnormal;
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

          // compute partial_area_s[2*i+0], partial_area_s[2*i+1] = darea_s(Q,global_coords(Q))/dQ = darea_s/dQ + darea_s/dglobal_coords*dglobal_coords/dQ
          face->Compute_Partial_Local_Coordinates_1(POSITION, global_coords, face_dlocal_coords);
          face->Compute_Partial_Local_Coordinates_2(POSITION, global_coords, dmdX, dmdY, dmdZ);
          for (j=0; j<nbDer; ++j) {
            partial_area_s[j][2*i+0] = dmdX[0]*p->p[i].partial_x[j] + dmdY[0]*p->p[i].partial_y[j] + dmdZ[0]*p->p[i].partial_z[j];
            partial_area_s[j][2*i+1] = dmdX[1]*p->p[i].partial_x[j] + dmdY[1]*p->p[i].partial_y[j] + dmdZ[1]*p->p[i].partial_z[j];
          }
          for (j=0; j<3*face->Nodes_Per_Face(); ++j) {
            partial_area_s[j][2*i+0] += face_dlocal_coords[j][0];
            partial_area_s[j][2*i+1] += face_dlocal_coords[j][1];
          }
          // compute partial_area_m[2*i+0], partial_area_m[2*i+1]
          master_face->Compute_Partial_Local_Coordinates_1(POSITION, global_coords, master_face_dlocal_coords);
          master_face->Compute_Partial_Local_Coordinates_2(POSITION, global_coords, dmdX, dmdY, dmdZ);
          for (j=0; j<nbDer; ++j) {
            partial_area_m[j][2*i+0] = dmdX[0]*p->p[i].partial_x[j] + dmdY[0]*p->p[i].partial_y[j] + dmdZ[0]*p->p[i].partial_z[j];
            partial_area_m[j][2*i+1] = dmdX[1]*p->p[i].partial_x[j] + dmdY[1]*p->p[i].partial_y[j] + dmdZ[1]*p->p[i].partial_z[j];
          }
          for (j=0; j<3*master_face->Nodes_Per_Face(); ++j) {
            int J = 3*face->Nodes_Per_Face()+j;
            partial_area_m[J][2*i+0] += master_face_dlocal_coords[j][0];
            partial_area_m[J][2*i+1] += master_face_dlocal_coords[j][1];
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
              ifaceedge, iedge_m, area_s, area_m, partial_area_s, partial_area_m );
        delete [] face_dlocal_coords;
        delete [] master_face_dlocal_coords;
      }
    }

  } else {
    std::cerr << " *** ERROR: ContactSearch::Partial_Face_Face_Search is not implemented for non-planar case\n";
    exit(-1);
  }
  delete [] partial_area_s;
  delete [] partial_area_m;
  return cffi;
}

