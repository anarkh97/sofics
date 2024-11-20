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


#include "ContactAnalyticPlane.h"
#include "ContactErrors.h"
#include "ContactParOStream.h"
#include <iostream>
#include <cmath>

ContactAnalyticPlane::ContactAnalyticPlane( int ID, int Key, const Real* data )
  : ContactAnalyticSurface( ContactSearch::PLANE, ID, Key )
{
  point[0]         = data[0];
  point[1]         = data[1];
  point[2]         = data[2];
  normal_vector[0] = data[3];
  normal_vector[1] = data[4];
  normal_vector[2] = data[5];

  // Ensure that vector is a unit vector
  Real vmag = std::sqrt( normal_vector[0]*normal_vector[0] + 
		         normal_vector[1]*normal_vector[1] +
		         normal_vector[2]*normal_vector[2] );
  if( vmag > 0.0){
    vmag = 1.0/vmag;
    normal_vector[0] *= vmag;
    normal_vector[1] *= vmag;
    normal_vector[2] *= vmag;
  } else {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "ContactAnalyticPlane: The normal vector is zero_length" << std::endl;
    POSTCONDITION(0);
#endif
  }
}

ContactAnalyticPlane::ContactAnalyticPlane( const ContactAnalyticPlane& plane )
  : ContactAnalyticSurface( ContactSearch::PLANE, 
			    (ContactAnalyticSurface*)&plane )
{
  point[0]         = plane.point[0];
  point[1]         = plane.point[1];
  point[2]         = plane.point[2];
  normal_vector[0] = plane.normal_vector[0];
  normal_vector[1] = plane.normal_vector[1];
  normal_vector[2] = plane.normal_vector[2];
}

ContactAnalyticPlane::ContactAnalyticPlane( int ID )
  : ContactAnalyticSurface( ContactSearch::PLANE, ID )
{
  point[0]         = 0.0;
  point[1]         = 0.0;
  point[2]         = 0.0;
  normal_vector[0] = 0.0;
  normal_vector[1] = 0.0;
  normal_vector[2] = 0.0;
}


ContactSearch::ContactErrorCode 
ContactAnalyticPlane::Check_for_Errors( ContactErrors* errors )
{
  // Check the normal_vector is not a zero length vector
  Real vmag = std::sqrt( normal_vector[0]*normal_vector[0] + 
		         normal_vector[1]*normal_vector[1] +
		         normal_vector[2]*normal_vector[2] );
  if( vmag ) return ContactSearch::NO_ERROR;
  errors->Add_Error_Message("AnalyticPlane Normal Vector is zero length");
  return ContactSearch::INVALID_DATA;
}

ContactAnalyticPlane::~ContactAnalyticPlane() {}

void ContactAnalyticPlane::ComputeBoundingBox(ContactBoundingBox* bb)
{
  // since a plane is essentially unbounded, we artificially bound it 
  // by finding it's intersection with the topology's global bounding box.
  
  Real p[3], p0[3], p1[3];
  
  // loop thru all the edges of the topology's global bounding box 
  // and find it's intersection with the plane and add it to the 
  // plane's bounding box
  for (int i=0; i<10; ++i) {
    switch (i) {
    case 0:
      p0[0] = bb->get_x_min();
      p0[1] = bb->get_y_min();
      p0[2] = bb->get_z_min();
      p1[0] = bb->get_x_max();
      p1[1] = bb->get_y_min();
      p1[2] = bb->get_z_min();
      break;
    case 1:
      p0[0] = bb->get_x_min();
      p0[1] = bb->get_y_min();
      p0[2] = bb->get_z_min();
      p1[0] = bb->get_x_min();
      p1[1] = bb->get_y_max();
      p1[2] = bb->get_z_min();
      break;
    case 2:
      p0[0] = bb->get_x_min();
      p0[1] = bb->get_y_min();
      p0[2] = bb->get_z_min();
      p1[0] = bb->get_x_min();
      p1[1] = bb->get_y_min();
      p1[2] = bb->get_z_max();
      break;
      
    case 3:
      p0[0] = bb->get_x_max();
      p0[1] = bb->get_y_min();
      p0[2] = bb->get_z_max();
      p1[0] = bb->get_x_max();
      p1[1] = bb->get_y_min();
      p1[2] = bb->get_z_min();
      break;
    case 4:
      p0[0] = bb->get_x_max();
      p0[1] = bb->get_y_min();
      p0[2] = bb->get_z_max();
      p1[0] = bb->get_x_max();
      p1[1] = bb->get_y_max();
      p1[2] = bb->get_z_max();
      break;
    case 5:
      p0[0] = bb->get_x_max();
      p0[1] = bb->get_y_min();
      p0[2] = bb->get_z_max();
      p1[0] = bb->get_x_min();
      p1[1] = bb->get_y_min();
      p1[2] = bb->get_z_max();
      break;
      
    case 6:
      p0[0] = bb->get_x_max();
      p0[1] = bb->get_y_max();
      p0[2] = bb->get_z_min();
      p1[0] = bb->get_x_min();
      p1[1] = bb->get_y_max();
      p1[2] = bb->get_z_min();
      break;
    case 7:
      p0[0] = bb->get_x_max();
      p0[1] = bb->get_y_max();
      p0[2] = bb->get_z_min();
      p1[0] = bb->get_x_max();
      p1[1] = bb->get_y_min();
      p1[2] = bb->get_z_min();
      break;
    case 8:
      p0[0] = bb->get_x_max();
      p0[1] = bb->get_y_max();
      p0[2] = bb->get_z_min();
      p1[0] = bb->get_x_max();
      p1[1] = bb->get_y_max();
      p1[2] = bb->get_z_max();
      break;
      
    case 9:
      p0[0] = bb->get_x_min();
      p0[1] = bb->get_y_max();
      p0[2] = bb->get_z_max();
      p1[0] = bb->get_x_min();
      p1[1] = bb->get_y_max();
      p1[2] = bb->get_z_min();
      break;
    case 10:
      p0[0] = bb->get_x_min();
      p0[1] = bb->get_y_max();
      p0[2] = bb->get_z_max();
      p1[0] = bb->get_x_min();
      p1[1] = bb->get_y_min();
      p1[2] = bb->get_z_max();
      break;
    case 11:
      p0[0] = bb->get_x_min();
      p0[1] = bb->get_y_max();
      p0[2] = bb->get_z_max();
      p1[0] = bb->get_x_max();
      p1[1] = bb->get_y_max();
      p1[2] = bb->get_z_max();
      break;
    }
    Real t0 = (p0[0]-point[0])*normal_vector[0]+
              (p0[1]-point[1])*normal_vector[1]+
              (p0[2]-point[2])*normal_vector[2];
    Real t1 = (p1[0]-point[0])*normal_vector[0]+
              (p1[1]-point[1])*normal_vector[1]+
              (p1[2]-point[2])*normal_vector[2];
    if (t0==0.0) {
      bounding_box.add_point(p0);
    } else if (t1==0.0) {
      bounding_box.add_point(p1);
    } else  if ((t0>0.0 && t1<0.0) ||
                (t0<0.0 && t1>0.0)) {
      p[0]   = p0[0]+(p1[0]-p0[0])*t0/(t0-t1);
      p[1]   = p0[1]+(p1[1]-p0[1])*t0/(t0-t1);
      p[2]   = p0[2]+(p1[2]-p0[2])*t0/(t0-t1);
      bounding_box.add_point(p);
    }
  }   
}

ContactSearch::ContactErrorCode 
ContactAnalyticPlane::Set_Configuration( const Real* data )
{
  point[0]         = data[0];
  point[1]         = data[1];
  point[2]         = data[2];
  normal_vector[0] = data[3];
  normal_vector[1] = data[4];
  normal_vector[2] = data[5];

  // Ensure that vector is a unit vector
  Real vmag = std::sqrt( normal_vector[0]*normal_vector[0] + 
		         normal_vector[1]*normal_vector[1] +
		         normal_vector[2]*normal_vector[2] );
  if( vmag > 0.0){
    vmag = 1.0/vmag;
    normal_vector[0] *= vmag;
    normal_vector[1] *= vmag;
    normal_vector[2] *= vmag;
  } else {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "ContactAnalyticPlane: The normal vector is zero_length" << std::endl;
#endif
    return ContactSearch::INVALID_DATA;
  }
  return ContactSearch::NO_ERROR;
}

int ContactAnalyticPlane::Restart_Size()
{
  return 6;
}

int ContactAnalyticPlane::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;

  restart_data[words_added++] = point[0];
  restart_data[words_added++] = point[1];
  restart_data[words_added++] = point[2];
  restart_data[words_added++] = normal_vector[0];
  restart_data[words_added++] = normal_vector[1];
  restart_data[words_added++] = normal_vector[2];

  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactAnalyticPlane::Implant_Restart_Data( Real* restart_data )
{
  int words_read = 0;

  point[0]         = restart_data[words_read++];
  point[1]         = restart_data[words_read++];
  point[2]         = restart_data[words_read++];
  normal_vector[0] = restart_data[words_read++];
  normal_vector[1] = restart_data[words_read++];
  normal_vector[2] = restart_data[words_read++];

  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

void ContactAnalyticPlane::Bounding_Box( Real* min, Real* max )
{
  // Just set min and max to the entire space since a plane is unbounded in
  // two dimensions.  A cartesian box can be limited only in the special case
  // of the plane being normal to one of the cartesian directions.
  min[0] = -1.E30;
  min[1] = -1.E30;
  min[2] = -1.E30;
  max[0] = +1.E30;
  max[1] = +1.E30;
  max[2] = +1.E30;
}

bool ContactAnalyticPlane::Process( Real* node_position,
                                    Real& gap, 
                                    Real* contact_point,
                                    Real* surface_normal,
                                    Real* pushback_dir,
                                    Real& time_to_contact,
                                    int&  location,
                                    bool PRINT_THIS_NODE,
                                    ContactParOStream* postream )
{
  // construct a vector for the point defining the plane to the node
  Real vec_p_n_x = node_position[0] - point[0];
  Real vec_p_n_y = node_position[1] - point[1];
  Real vec_p_n_z = node_position[2] - point[2];

  // Dot this vector onto the normal to the plane
  Real dot = vec_p_n_x*normal_vector[0] + 
             vec_p_n_y*normal_vector[1] +
             vec_p_n_z*normal_vector[2] ;

  if( dot > 0.0 ) return false; // The point is outside the plane

  gap = dot;

  surface_normal[0] = normal_vector[0];
  surface_normal[1] = normal_vector[1];
  surface_normal[2] = normal_vector[2];

  contact_point[0] = node_position[0] - gap*normal_vector[0];
  contact_point[1] = node_position[1] - gap*normal_vector[1];
  contact_point[2] = node_position[2] - gap*normal_vector[2];
  
  if (gap==0.0) {
    pushback_dir[0] = surface_normal[0];
    pushback_dir[1] = surface_normal[1];
    pushback_dir[2] = surface_normal[2];
  } else {
    pushback_dir[0] = (contact_point[0] - node_position[0]) / std::fabs(gap);
    pushback_dir[1] = (contact_point[1] - node_position[1]) / std::fabs(gap);
    pushback_dir[2] = (contact_point[2] - node_position[2]) / std::fabs(gap);
  }
  
  time_to_contact = 0.0;
  location        = 1;

  return true;
}

bool ContactAnalyticPlane::Process( Real* node_position,    //1st configuration
              Real* node_position_2,  //2nd configuration (aug or predicted); 
              Real& penetration_mag, 
              Real* contact_point,
              Real* surface_normal,
              Real* pushback_dir,
              Real& time_to_contact,
              int&  location,
              bool PRINT_THIS_NODE,
              ContactParOStream* postream)
{    
  //setup the other values like normal
  bool valid = Process(node_position, penetration_mag, contact_point, surface_normal, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, postream);
  
  //if interaction is valid then calculate the normalized time to contact
  if(valid == true){
    // construct a vector for the point defining the plane to the node
    Real vec_p_n_x = node_position_2[0] - point[0];
    Real vec_p_n_y = node_position_2[1] - point[1];
    Real vec_p_n_z = node_position_2[2] - point[2];

    // Dot this vector onto the normal to the plane
    Real pgap = vec_p_n_x*normal_vector[0] + 
                vec_p_n_y*normal_vector[1] +
                vec_p_n_z*normal_vector[2] ;

    if((pgap-penetration_mag) == 0.0){
      time_to_contact = 0.0;
    } else{
      time_to_contact = -penetration_mag/(pgap - penetration_mag);
    }
  }
  
  return valid;
}

void ContactAnalyticPlane::Display(ContactParOStream& postream)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream << "\n  Analytic Plane" << "\n";
  postream << "\tID = " << Global_ID() << "\n";
  postream << "\tPoint:" << point[0] << " " << point[1] << " " 
         << point[2] << "\n";
  postream << "\tNormal: " << normal_vector[0] << " " << normal_vector[1] << " "
         << normal_vector[2] << "\n";
#endif
}
