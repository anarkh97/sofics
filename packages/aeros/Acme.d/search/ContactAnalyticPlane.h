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


#ifndef ContactAnalyticPlane_h_
#define ContactAnalyticPlane_h_

#include "ContactAnalyticSurface.h"
#include "ContactSearch.h"
#include "Contact_Defines.h"

class ContactParOStream;

class ContactAnalyticPlane : public ContactAnalyticSurface {

 public:

  ContactAnalyticPlane( int id, int key, const Real* data );
  ContactAnalyticPlane( const ContactAnalyticPlane& );
  ContactAnalyticPlane( int id );
  ~ContactAnalyticPlane();

  ContactSearch::ContactErrorCode Check_for_Errors( ContactErrors* );
  ContactSearch::ContactErrorCode Set_Configuration( const Real* data );

  void Bounding_Box( Real* min, Real* max );
  void ComputeBoundingBox(ContactBoundingBox* bb);
  virtual bool Process( Real* node_position,
                        Real& penetration_mag,
                        Real* contact_point,
                        Real* surface_normal,
                        Real* pushback_dir,
                        Real& time_to_contact,
                        int&  location,
                        bool PRINT_THIS_NODE,
                        ContactParOStream* postream );
  virtual bool Process( Real* node_position,    //1st configuration
                        Real* node_position_2,  //2nd configuration (aug or predicted); 
                        Real& penetration_mag, 
			Real* contact_point,
                        Real* surface_normal,
                        Real* pushback_dir,
                        Real& time_to_contact,
                        int&  location,
                        bool PRINT_THIS_NODE,
                        ContactParOStream* postream );
  void Display( ContactParOStream& );

  virtual int Restart_Size();
  virtual int Extract_Restart_Data( Real* );
  virtual int Implant_Restart_Data( Real* );

 private:
  
  Real point[3];
  Real normal_vector[3];

};

#endif // ContactAnalyticPlane_h_
