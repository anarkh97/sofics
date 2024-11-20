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


#include "ContactAnalyticCylinderInside.h"
#include "ContactErrors.h"
#include "ContactParOStream.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

ContactAnalyticCylinderInside::ContactAnalyticCylinderInside( int ID, int Key, 
						              const Real* data )
  : ContactAnalyticSurface( ContactSearch::CYLINDER_INSIDE, ID, Key )
{
  center[0]    = data[0];
  center[1]    = data[1];
  center[2]    = data[2];
  axial_dir[0] = data[3];
  axial_dir[1] = data[4];
  axial_dir[2] = data[5];
  radius       = data[6];
  length       = data[7]*0.5;  // We want 1/2 the length

  // Ensure that vector is a unit vector
  Real vmag = std::sqrt( axial_dir[0]*axial_dir[0] + 
		         axial_dir[1]*axial_dir[1] +
		         axial_dir[2]*axial_dir[2] );
  if( vmag > 0.0){
    vmag = 1.0/vmag;
    axial_dir[0] *= vmag;
    axial_dir[1] *= vmag;
    axial_dir[2] *= vmag;
  } else {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "AnalyticCylinderInside: The axial direction vector is 0 length" 
	      << std::endl;
    POSTCONDITION(0);
#endif
  }
}

ContactAnalyticCylinderInside::ContactAnalyticCylinderInside( 
		           const ContactAnalyticCylinderInside& cylinder )
  : ContactAnalyticSurface( ContactSearch::CYLINDER_INSIDE,
                            (ContactAnalyticSurface*)&cylinder )
{
  center[0]    = cylinder.center[0];
  center[1]    = cylinder.center[1];
  center[2]    = cylinder.center[2];
  axial_dir[0] = cylinder.axial_dir[0];
  axial_dir[1] = cylinder.axial_dir[1];
  axial_dir[2] = cylinder.axial_dir[2];
  radius       = cylinder.radius;
  length       = cylinder.length;
  bounding_box.set_box(cylinder.bounding_box);

}

ContactAnalyticCylinderInside::ContactAnalyticCylinderInside( int ID )
  : ContactAnalyticSurface( ContactSearch::CYLINDER_INSIDE, ID )
{
  center[0]    = 0.0;
  center[1]    = 0.0;
  center[2]    = 0.0;
  axial_dir[0] = 0.0;
  axial_dir[1] = 0.0;
  axial_dir[2] = 0.0;
  radius       = 0.0;
  length       = 0.0;
}

void ContactAnalyticCylinderInside::ComputeBoundingBox(ContactBoundingBox* bb)
{
  int i;
  Real u[3];       // orthoganal vector to the cylinder axis;
  Real v[3];       // orthoganal vector to the cylinder axis and u vector 
  Real d0[3];      // endpoint of radius along u vector
  Real d1[3];      // endpoint of radius along v vector
  Real c0[3];      // center of 1st endcap
  Real c1[3];      // center of 2nd endcap
  Real point0[3];  // radius endpoints on the 1st endcap at +/- u and +/i v vectors
  Real point1[3];  // radius endpoints on the 2nd endcap at +/- u and +/i v vectors
  
  bounding_box.Reset();
  
  // determine cross axes to the axial direction
  if (std::fabs(axial_dir[2])<std::fabs(axial_dir[0])) {
    if (std::fabs(axial_dir[2])<std::fabs(axial_dir[1])) {
      u[0] = -axial_dir[1];
      u[1] =  axial_dir[0];
      u[2] =  0.0;
    } else {
      u[0] = -axial_dir[2];
      u[1] =  0.0;
      u[2] =  axial_dir[0];
    }
  } else {
    u[0] =  0.0;
    u[1] = -axial_dir[2];
    u[2] =  axial_dir[1];
  }
  v[0] = (u[1]*axial_dir[2]) - (u[2]*axial_dir[1]);
  v[1] = (u[2]*axial_dir[0]) - (u[0]*axial_dir[2]);
  v[2] = (u[0]*axial_dir[1]) - (u[1]*axial_dir[0]);
  
  // determine radial distances along cross axes
  for (i=0; i<3; ++i) {
    d0[i] = radius*u[i];
    d1[i] = radius*v[i];
  }
  for (i=0; i<3; ++i) {
    c0[i] = center[i]-length*axial_dir[i];
    c1[i] = center[i]+length*axial_dir[i];
  }
  for (i=0; i<3; ++i) {
    point0[i] = c0[i]+d0[i]+d1[i];
    point1[i] = c1[i]+d0[i]+d1[i];
  }
  bounding_box.add_point(point0);
  bounding_box.add_point(point1);
  
  for (i=0; i<3; ++i) {
    point0[i] = c0[i]-d0[i]+d1[i];
    point1[i] = c1[i]-d0[i]+d1[i];
  }
  bounding_box.add_point(point0);
  bounding_box.add_point(point1);
  
  for (i=0; i<3; ++i) {
    point0[i] = c0[i]+d0[i]-d1[i];
    point1[i] = c1[i]+d0[i]-d1[i];
  }
  bounding_box.add_point(point0);
  bounding_box.add_point(point1);
  
  for (i=0; i<3; ++i) {
    point0[i] = c0[i]-d0[i]-d1[i];
    point1[i] = c1[i]-d0[i]-d1[i];
  }
  bounding_box.add_point(point0);
  bounding_box.add_point(point1);
}


int ContactAnalyticCylinderInside::Restart_Size()
{
  return 8;
}

int ContactAnalyticCylinderInside::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;

  restart_data[words_added++] = center[0];
  restart_data[words_added++] = center[1];
  restart_data[words_added++] = center[2];
  restart_data[words_added++] = axial_dir[0];
  restart_data[words_added++] = axial_dir[1];
  restart_data[words_added++] = axial_dir[2];
  restart_data[words_added++] = radius;
  restart_data[words_added++] = length;

  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactAnalyticCylinderInside::Implant_Restart_Data( Real* restart_data )
{
  int words_read = 0;

  center[0]    = restart_data[words_read++];
  center[1]    = restart_data[words_read++];
  center[2]    = restart_data[words_read++];
  axial_dir[0] = restart_data[words_read++];
  axial_dir[1] = restart_data[words_read++];
  axial_dir[2] = restart_data[words_read++];
  radius       = restart_data[words_read++];
  length       = restart_data[words_read++];

  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode 
ContactAnalyticCylinderInside::Check_for_Errors( ContactErrors* errors )
{
  ContactSearch::ContactErrorCode error = ContactSearch::NO_ERROR;

  // Check the axial_dir is not a zero length vector
  Real vmag = std::sqrt( axial_dir[0]*axial_dir[0] + 
		         axial_dir[1]*axial_dir[1] +
		         axial_dir[2]*axial_dir[2] );
  if( vmag==0 ){ 
    error = ContactSearch::INVALID_DATA;
    errors->Add_Error_Message("AnalyticCylinderInside Normal Vector is zero length.");
  }

  // Check the radius
  if( radius <= 0 ){
    error = ContactSearch::INVALID_DATA;
    errors->Add_Error_Message("AnalyticCylinder Radius is non-positive.");
  }

  // Check the length
  if( length <= 0 ){
    error = ContactSearch::INVALID_DATA;
    errors->Add_Error_Message("AnalyticCylinder Lenght is non-positive.");
  }

  return error;
}

ContactAnalyticCylinderInside::~ContactAnalyticCylinderInside() {}

ContactSearch::ContactErrorCode
ContactAnalyticCylinderInside::Set_Configuration( const Real* data )
{
  center[0]    = data[0];
  center[1]    = data[1];
  center[2]    = data[2];
  axial_dir[0] = data[3];
  axial_dir[1] = data[4];
  axial_dir[2] = data[5];
  radius       = data[6];
  length       = data[7]*0.5;  // We want 1/2 the length

  // Ensure that vector is a unit vector
  Real vmag = std::sqrt( axial_dir[0]*axial_dir[0] + 
		         axial_dir[1]*axial_dir[1] +
		         axial_dir[2]*axial_dir[2] );
  if( vmag > 0.0){
    vmag = 1.0/vmag;
    axial_dir[0] *= vmag;
    axial_dir[1] *= vmag;
    axial_dir[2] *= vmag;
  } else {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "AnalyticCylinderInside: The axial direction vector is 0 length" 
	 << std::endl;
#endif
    return ContactSearch::INVALID_DATA;
  }
  return ContactSearch::NO_ERROR;
}

void ContactAnalyticCylinderInside::Bounding_Box( Real* min, Real* max )
{
  // Just set min and max to the entire space.  For efficiency, this should
  // be changed.
  min[0] = -1.E30;
  min[1] = -1.E30;
  min[2] = -1.E30;
  max[0] = +1.E30;
  max[1] = +1.E30;
  max[2] = +1.E30;
}

bool ContactAnalyticCylinderInside::Process( Real* node_position,
                                             Real& gap, 
					     Real* contact_point,
					     Real* surface_normal,
                                             Real* pushback_dir,
                                             Real& time_to_contact,
                                             int&  location,
                                             bool PRINT_THIS_NODE,
                                             ContactParOStream* postream )
{
  // construct a vector from the center to the node
  Real vec_p_n_x = node_position[0] - center[0];
  Real vec_p_n_y = node_position[1] - center[1];
  Real vec_p_n_z = node_position[2] - center[2];

  // Dot this vector onto the axial direction of the cylinder
  Real dot = vec_p_n_x*axial_dir[0] + 
             vec_p_n_y*axial_dir[1] +
             vec_p_n_z*axial_dir[2] ;

  // Compute the radial distance from the centerline
  Real vec_rad_x = vec_p_n_x - dot*axial_dir[0];
  Real vec_rad_y = vec_p_n_y - dot*axial_dir[1];
  Real vec_rad_z = vec_p_n_z - dot*axial_dir[2];
  Real rad_distance = std::sqrt( vec_rad_x*vec_rad_x + 
                                 vec_rad_y*vec_rad_y +
			         vec_rad_z*vec_rad_z );

  // If The point is inside the cap and the radius, its not interacting.
  if( rad_distance < radius && std::fabs(dot) < length ) return false;

#ifdef CONTACT_DEBUG_NODE
  if(PRINT_THIS_NODE && postream != NULL){
    *postream << "Interaction with Analytic Surface " << ProcArrayIndex() << " found.\n";
  }
#endif
  Real radial_gap = radius - rad_distance;
  Real axial_gap  = length - std::fabs(dot);

  // Choose the "minimum" gap (this corresponds to the least negative)
#define RADIAL 0
#define AXIAL  1
#define EDGE   2
  int method;
  if( radial_gap < 0.0 && axial_gap < 0.0 ){
    method = EDGE;
  } else if( radial_gap < 0 ) {
    method = RADIAL;
  } else if( axial_gap < 0 ){
    method = AXIAL;
  } else {
    // The case where both gaps are zero. Take the axial
    method = AXIAL;
  }
  
#ifdef CONTACT_DEBUG_NODE
  if(PRINT_THIS_NODE && postream != NULL){
    *postream << "  Interaction method is: ";
    switch(method){
      case RADIAL:
        *postream << "RADIAL\n";
        *postream << "  Radial gap = " << radial_gap << "\n";
        break;
      case AXIAL:
        *postream << "AXIAL\n";
        *postream << "  Axial gap  = " << axial_gap  << "\n";
        break;
      case EDGE:
        *postream << "EDGE\n";
        *postream << "  Radial gap = " << radial_gap << "\n"
                   << "  Axial gap  = " << axial_gap  << "\n";
        break;
    }
  }
#endif
  if( method == AXIAL ){
    gap = axial_gap;
    int sign = dot > 0 ? -1 : 1;
    surface_normal[0] = axial_dir[0]*sign;
    surface_normal[1] = axial_dir[1]*sign;
    surface_normal[2] = axial_dir[2]*sign;
    contact_point[0]  = axial_dir[0]*length+vec_rad_x;
    contact_point[1]  = axial_dir[1]*length+vec_rad_y;
    contact_point[2]  = axial_dir[2]*length+vec_rad_z;
  } else if( method == RADIAL ){
    gap = radial_gap;
    surface_normal[0] = -vec_rad_x/rad_distance;
    surface_normal[1] = -vec_rad_y/rad_distance;
    surface_normal[2] = -vec_rad_z/rad_distance;
    contact_point[0]  = center[0]+dot*axial_dir[0]+(radius*-surface_normal[0]);
    contact_point[1]  = center[1]+dot*axial_dir[1]+(radius*-surface_normal[1]);
    contact_point[2]  = center[2]+dot*axial_dir[2]+(radius*-surface_normal[2]);
  } else {
    //this is a little weird and probably wrong. Needs to be looked at more.
    int sgn = dot > 0 ? 1 : -1;
    Real cap_x     = center[0]+sgn*length*axial_dir[0];
    Real cap_y     = center[1]+sgn*length*axial_dir[1];
    Real cap_z     = center[2]+sgn*length*axial_dir[2];
    Real vec_c_p_x = node_position[0] - cap_x;
    Real vec_c_p_y = node_position[1] - cap_y;
    Real vec_c_p_z = node_position[2] - cap_z;
    Real dot_prod  = vec_c_p_x*axial_dir[0] + 
                     vec_c_p_y*axial_dir[1] +
                     vec_c_p_z*axial_dir[2];
    Real rad_dir_x = vec_c_p_x-dot_prod*axial_dir[0];
    Real rad_dir_y = vec_c_p_y-dot_prod*axial_dir[1];
    Real rad_dir_z = vec_c_p_z-dot_prod*axial_dir[2];
    Real rdmag     = std::sqrt( rad_dir_x*rad_dir_x + 
                                rad_dir_y*rad_dir_y +
		                rad_dir_z*rad_dir_z );
    contact_point[0] = cap_x + radius*rad_dir_x/rdmag;
    contact_point[1] = cap_y + radius*rad_dir_y/rdmag;
    contact_point[2] = cap_z + radius*rad_dir_z/rdmag;
    // Make the surface normal the direction from the contact point
    // to the node position.  This will vary this normal continuously.
    surface_normal[0] = node_position[0] - contact_point[0];
    surface_normal[1] = node_position[1] - contact_point[1];
    surface_normal[2] = node_position[2] - contact_point[2];
    gap = std::sqrt( surface_normal[0]*surface_normal[0] +
		     surface_normal[1]*surface_normal[1] +
		     surface_normal[2]*surface_normal[2] );
    surface_normal[0] /= gap; 
    surface_normal[1] /= gap; 
    surface_normal[2] /= gap; 
  }
  
  if (gap==0.0) {
    pushback_dir[0] = surface_normal[0];
    pushback_dir[1] = surface_normal[1];
    pushback_dir[2] = surface_normal[2];
  } else {
    pushback_dir[0] = (contact_point[0] - node_position[0]) / std::fabs(gap);
    pushback_dir[1] = (contact_point[1] - node_position[1]) / std::fabs(gap);
    pushback_dir[2] = (contact_point[2] - node_position[2]) / std::fabs(gap);
  }

#ifdef CONTACT_DEBUG_NODE
  if(PRINT_THIS_NODE && postream != NULL){
     *postream << " Surface Normal = " << surface_normal[0] << ","
                                      << surface_normal[1] << ","
                                      << surface_normal[2] << "\n"
              << "  Contact Point = " << contact_point[0]  << ","
                                      << contact_point[1]  << ","
                                      << contact_point[1]  << "\n"
              << "  Pushback Dir  = " << pushback_dir[0]   << ","
                                      << pushback_dir[1]   << ","
                                      << pushback_dir[1]   << "\n";
  }
#endif
  
  time_to_contact = 0.0;
  location        = 1;

  return true;
}

bool ContactAnalyticCylinderInside::Process( Real* node_position,    //1st configuration
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
    // construct a vector from the center to the node
    Real vec_p_n_x = node_position_2[0] - center[0];
    Real vec_p_n_y = node_position_2[1] - center[1];
    Real vec_p_n_z = node_position_2[2] - center[2];

    // Dot this vector onto the axial direction of the cylinder
    Real dot = vec_p_n_x*axial_dir[0] + 
               vec_p_n_y*axial_dir[1] +
               vec_p_n_z*axial_dir[2] ;

    // Compute the radial distance from the centerline
    Real vec_rad_x = vec_p_n_x - dot*axial_dir[0];
    Real vec_rad_y = vec_p_n_y - dot*axial_dir[1];
    Real vec_rad_z = vec_p_n_z - dot*axial_dir[2];
    Real rad_distance = std::sqrt( vec_rad_x*vec_rad_x + 
                                   vec_rad_y*vec_rad_y +
			           vec_rad_z*vec_rad_z );

    Real radial_gap = radius - rad_distance;
    Real axial_gap  = length - std::fabs(dot);
    Real pgap;

    if( radial_gap < 0.0 && axial_gap < 0.0 ){
      Real normal[3];
      normal[0] = node_position_2[0] - contact_point[0];
      normal[1] = node_position_2[1] - contact_point[1];
      normal[2] = node_position_2[2] - contact_point[2];
      pgap = std::sqrt(normal[0]*normal[0]+
                       normal[1]*normal[1]+
                       normal[2]*normal[2]);
    } else if( radial_gap < 0 ) {
      pgap = radial_gap;
    } else if( axial_gap < 0 ){
      pgap = axial_gap;
    } else {
      // The case where both gaps are zero. Take the axial
      pgap = axial_gap;
    } 
    
    if(pgap == penetration_mag){
      time_to_contact = 0.0;
    } else{
      time_to_contact = -penetration_mag/(pgap - penetration_mag);
    }
  }
  
  return valid;
}

void ContactAnalyticCylinderInside::Display(ContactParOStream& postream)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream << "\n  Analytic Cylinder_Inside" << "\n";
  postream << "\tID = " << Global_ID() << "\n";
  postream << "\tCenter:" << center[0] << " " << center[1] << " " 
         << center[2] << "\n";
  postream << "\tAxial Dir: " << axial_dir[0] << " " << axial_dir[1] << " "
         << axial_dir[2] << "\n";
  postream << "\tRadius: " << radius << "\n";
  postream << "\tLength: " << length*2.0 << "\n";
#endif
}
