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


#ifndef ContactTDFaceFaceEnf_h_
#define ContactTDFaceFaceEnf_h_

#ifdef CONTACT_TD_FACE_FACE_ENF

#include "ContactSearch.h"
#include "ContactEnforcement.h"

class ContactTDFaceFaceEnf : public ContactEnforcement {

 public:

  enum Enforcement_Data_Index{ KINEMATIC_PARTITION=0,
			       NSIZED };

  ContactTDFaceFaceEnf( const double*, ContactSearch*,
			ContactSearch::ContactErrorCode& error );
  ContactTDFaceFaceEnf( ContactSearch*, const double* restart_data,
			ContactSearch::ContactErrorCode& error );
  ~ContactTDFaceFaceEnf();

  ContactSearch::ContactErrorCode Compute_Forces(const Real* mass, Real* force);
  
  // regression test output functions
  virtual int Number_of_Nodal_Plot_Variables() { return number_nodal_plot_vars; };
  virtual int Number_of_Element_Plot_Variables() { return number_element_plot_vars; };
  virtual Real Get_Global_Plot_Variable( int ) { return 0.0; };
  virtual void Get_Nodal_Plot_Variable( int, Real* ) { return; };
  virtual void Get_Element_Plot_Variable( int, Real* ) { return; };

  // Restart Functions
  virtual int Number_General_Restart_Variables() {return 0;};
  virtual int Number_Nodal_Restart_Variables() {return 0;};
  virtual int Number_Edge_Restart_Variables() {return 0;};
  virtual int Number_Face_Restart_Variables() {return 0;};
  virtual int Number_Element_Restart_Variables() {return 0;};
  virtual ContactSearch::ContactErrorCode
    Extract_General_Restart_Variable( Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable( int n, Real* data, int* node_ids ) 
       {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Extract_Edge_Restart_Variable( int n, Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Extract_Face_Restart_Variable( int n, Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Extract_Element_Restart_Variable( int n, Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Implant_General_Restart_Variable( Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable( int n, Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Implant_Edge_Restart_Variable( int n, Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Implant_Face_Restart_Variable( int n, Real* data ) 
      {return ContactSearch::NO_ERROR;};
  virtual ContactSearch::ContactErrorCode
    Implant_Element_Restart_Variable( int n, Real* data ) 
      {return ContactSearch::NO_ERROR;};
 
  virtual ContactSearch::ContactErrorCode
    Update_For_Topology_Change(int new_number_of_nodes,int* old_to_new_map)
      {return ContactSearch::NO_ERROR;};

  // restart functions
  virtual int Restart_Size() {return 0;};
  virtual ContactSearch::ContactErrorCode 
    Extract_Restart_Data(Real* restart_data ) {return ContactSearch::NO_ERROR;};


 private:
};

#endif
#endif
