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


#include "ContactAnalyticSphere.h"
#include "ContactErrors.h"
#include "ContactParOStream.h"
#include <iostream>
#include <cmath>
#include <cstdio>

ContactAnalyticSphere::ContactAnalyticSphere( int ID, int Key, const Real* data )
  : ContactAnalyticSurface( ContactSearch::SPHERE, ID, Key )
{
  center[0] = data[0];
  center[1] = data[1];
  center[2] = data[2];
  radius    = data[3];

  if( radius <= 0 ){ 
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "ContactAnalyticPlane: The radius must be positive" << std::endl;
    POSTCONDITION(0);
#endif
  }
}

ContactAnalyticSphere::ContactAnalyticSphere(const ContactAnalyticSphere& s)
  : ContactAnalyticSurface( ContactSearch::SPHERE, 
			    (ContactAnalyticSurface*)&s )
{
  center[0] = s.center[0];
  center[1] = s.center[1];
  center[2] = s.center[2];
  radius    = s.radius;
}

ContactAnalyticSphere::ContactAnalyticSphere( int ID )
  : ContactAnalyticSurface( ContactSearch::SPHERE, ID )
{
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;
  radius    = 0.0;
}

ContactAnalyticSphere::~ContactAnalyticSphere() {}


void ContactAnalyticSphere::ComputeBoundingBox(ContactBoundingBox* bb)
{
  bounding_box.Reset();
  bounding_box.add_sphere(center, radius);
}

ContactSearch::ContactErrorCode
ContactAnalyticSphere::Set_Configuration( const Real* data )
{
  center[0] = data[0];
  center[1] = data[1];
  center[2] = data[2];
  radius    = data[3];

  if( radius <= 0 ){ 
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "ContactAnalyticPlane: The radius must be positive" << std::endl;
#endif
    return ContactSearch::INVALID_DATA;
  }
  return ContactSearch::NO_ERROR;
}


int ContactAnalyticSphere::Restart_Size()
{
  return 4;
}

int ContactAnalyticSphere::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;

  restart_data[words_added++] = center[0];
  restart_data[words_added++] = center[1];
  restart_data[words_added++] = center[2];
  restart_data[words_added++] = radius;

  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactAnalyticSphere::Implant_Restart_Data( Real* restart_data )
{
  int words_read = 0;

  center[0] = restart_data[words_read++];
  center[1] = restart_data[words_read++];
  center[2] = restart_data[words_read++];
  radius    = restart_data[words_read++];

  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode 
ContactAnalyticSphere::Check_for_Errors( ContactErrors* errors )
{
  // Check the radius is positive
  if( radius > 0 ) return ContactSearch::NO_ERROR;
  char* message = new char[81];
  std::sprintf(message,"Error: Sphere Radius non-positive: Radius = %g",radius );
  errors->Add_Error_Message( message );
  delete [] message;
  return ContactSearch::INVALID_DATA;
}

void ContactAnalyticSphere::Bounding_Box( Real* min, Real* max )
{
  min[0] = center[0] - radius;
  min[1] = center[1] - radius;
  min[2] = center[2] - radius;
  max[0] = center[0] + radius;
  max[1] = center[1] + radius;
  max[2] = center[2] + radius;
}

bool ContactAnalyticSphere::Process( Real* node_position,
                                     Real& gap, 
				     Real* contact_center,
                                     Real* surface_normal,
                                     Real* pushback_dir,
                                     Real& time_to_contact,
                                     int&  location,
                                     bool PRINT_THIS_NODE,
                                     ContactParOStream* postream )
{
  // Compute the distance from the center of the sphere to the center
  Real xdis = node_position[0] - center[0];
  Real ydis = node_position[1] - center[1];
  Real zdis = node_position[2] - center[2];
  Real distance = std::sqrt( xdis*xdis + ydis*ydis + zdis*zdis );

  if( distance > radius ) return false;

  if( distance != 0 ){
    // Create a unit vector in the direction from the center to the node
    Real d = 1.0/distance;
    xdis  *= d;
    ydis  *= d;
    zdis  *= d;
    gap    = distance-radius;
    contact_center[0] = center[0] + radius*xdis;
    contact_center[1] = center[1] + radius*ydis;
    contact_center[2] = center[2] + radius*zdis;
    surface_normal[0] = xdis;
    surface_normal[1] = ydis;
    surface_normal[2] = zdis;
  } else {
    // the node position is exactly at the center of the sphere
    // the center along the x-axis
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Infinite # of closest points in ContactAnalyticSphere::Process"
	 << std::endl;
#endif
    gap = -radius;
    contact_center[0] = center[0] + radius;
    contact_center[1] = center[1];
    contact_center[2] = center[2];
    surface_normal[0] = 1.0;
    surface_normal[1] = 0.0;
    surface_normal[2] = 0.0;
  }
  
  if (gap==0.0) {
    pushback_dir[0] = surface_normal[0];
    pushback_dir[1] = surface_normal[1];
    pushback_dir[2] = surface_normal[2];
  } else {
    pushback_dir[0] = (contact_center[0] - node_position[0]) / std::fabs(gap);
    pushback_dir[1] = (contact_center[1] - node_position[1]) / std::fabs(gap);
    pushback_dir[2] = (contact_center[2] - node_position[2]) / std::fabs(gap);
  }
  
  time_to_contact = 0.0;
  location        = 1;

  return true;
}

bool ContactAnalyticSphere::Process( Real* node_position,    //1st configuration
              Real* node_position_2,  //2nd configuration (aug or predicted); 
              Real& penetration_mag, 
              Real* contact_point,
              Real* surface_normal,
              Real* pushback_dir,
              Real& time_to_contact,
              int&  location,
              bool PRINT_THIS_NODE,
              ContactParOStream* postream )
{
  //setup the other values like normal
  bool valid = Process(node_position, penetration_mag, contact_point, surface_normal, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, postream);
  
  //if interaction is valid then calculate the normalized time to contact
  if(valid == true){
    // Compute the distance from the center of the sphere to the center
    Real xdis = node_position_2[0] - center[0];
    Real ydis = node_position_2[1] - center[1];
    Real zdis = node_position_2[2] - center[2];
    Real distance = std::sqrt( xdis*xdis + ydis*ydis + zdis*zdis );
    Real pgap;
    if( distance != 0 ){
      pgap = distance-radius;
    } else{
      pgap = -radius;
    }
    
    if(pgap == penetration_mag){
      time_to_contact = 0.0;
    } else{
      time_to_contact = -penetration_mag/(pgap - penetration_mag);
    }
  }
  return valid;
}

void ContactAnalyticSphere::Display(ContactParOStream& postream)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream << "\n  Analytic Sphere" << "\n";
  postream << "\tID = " << Global_ID() << "\n";
  postream << "\tCenter: " << center[0] << " " << center[1] << " " 
         << center[2] << "\n";
  postream << "\tRadius: " << radius << "\n";
#endif
}
