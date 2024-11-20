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


#ifndef ContactMPCs_h_
#define ContactMPCs_h_

#include "ContactSearch.h"
#include "ContactEnforcement.h"

template<typename DataType> class ContactNode;
class ContactTopology;
class ContactAsymComm;
struct Contact_MPC_Equation;

class ContactMPCs : public ContactEnforcement {

 public:

  enum Enforcement_Data_Index { KINEMATIC_PARTITION=0, NSIZED };

  ContactMPCs( const double*, ContactSearch*,
			 ContactSearch::ContactErrorCode& );

  ContactMPCs( ContactSearch*, const double* restart_data,
			 ContactSearch::ContactErrorCode& );

  ~ContactMPCs();

  ContactSearch::ContactErrorCode Compute_MPCs();


  // Function to extract the MPCs
  ContactSearch::ContactErrorCode Number_of_MPC_Equations(int&);
  ContactSearch::ContactErrorCode Get_MPC_Equations(int, int*, int*, int*, int*,
						    int*, int*, int*, Real*);

  // regression test output functions
  virtual Real Get_Global_Plot_Variable( int ) { return -1.0; };
  virtual void Get_Nodal_Plot_Variable( int, Real* );
  virtual void Get_Element_Plot_Variable( int, Real* ) { return; };

  // Restart Functions
  virtual int Number_General_Restart_Variables() 
    {return Number_Base_General_Restart_Variables();};
  virtual int Number_Nodal_Restart_Variables() {return 0;};
  virtual int Number_Edge_Restart_Variables() {return 0;};
  virtual int Number_Face_Restart_Variables() {return 0;};
  virtual int Number_Element_Restart_Variables() {return 0;};
  virtual ContactSearch::ContactErrorCode
    Extract_General_Restart_Variable( double* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable( int n, double* data, int* node_ids );
  virtual ContactSearch::ContactErrorCode
    Extract_Edge_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Face_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Element_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_General_Restart_Variable( double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Edge_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Face_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Element_Restart_Variable( int n, double* data );

  virtual ContactSearch::ContactErrorCode
    Update_For_Topology_Change( int new_number_of_nodes, int* old_to_new_map );

 private:
  
  ContactMPCs(ContactMPCs&);
  ContactMPCs& operator=(ContactMPCs&);
  
  void Remove_MPCs();
  void Communicate_MPCs( );
  void Display_MPCs();

  ContactAsymComm* node_asymcomm;
  Contact_MPC_Equation** slave_mpcs;
  int num_face_mpcs;
  Contact_MPC_Equation** face_mpcs;

  //regression test
  void Plot_MPC_Vars( );
  void Set_NVARS( int num_elem_vars, Real* nodal_data );
  
};

struct Contact_MPC_Equation {
  int slave_node_proc_id;
  int slave_node_local_id;
  int master_face_proc_id;
  int master_face_local_id;
  int num_face_nodes;
  int face_node_proc_id[8];
  int face_node_local_id[8];
  double face_node_coeff[8];
};

#endif
