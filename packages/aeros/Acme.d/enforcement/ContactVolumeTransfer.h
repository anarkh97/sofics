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


#ifndef ContactVolumeTransfer_h_
#define ContactVolumeTransfer_h_

#include "ContactSearch.h"
#include "ContactEnforcement.h"

template<typename DataType> class ContactNode;
class ContactTopology;

class ContactVolumeTransfer : public ContactEnforcement {
  
 public:
  
  enum Enforcement_Data_Index { NEAR_TOLERANCE=0,
                                NSIZED };
  enum Contact_Status { };
  enum Contact_VolTran_Plot_Vars { };
  
  ContactVolumeTransfer( const Real*, ContactSearch*,
			 ContactSearch::ContactErrorCode& error );
  ContactVolumeTransfer( ContactSearch*, const Real* restart_data,
			 ContactSearch::ContactErrorCode& error );
  ~ContactVolumeTransfer();
  
  // regression test output functions
  virtual int Number_of_Global_Plot_Variables() { return 0; };
  virtual int Number_of_Nodal_Plot_Variables() { return number_nodal_plot_vars; };
  virtual int Number_of_Element_Plot_Variables() { return number_element_plot_vars; };
  virtual Real Get_Global_Plot_Variable( int ) { return -1.0; };
  virtual void Get_Nodal_Plot_Variable( int, Real* );
  virtual void Get_Element_Plot_Variable( int, Real* );

  ContactSearch::ContactErrorCode Get_Plot_Variable( Contact_VolTran_Plot_Vars, 
						     Real* buffer );
  
  ContactSearch::ContactErrorCode
    Compute_Volume_Transfer( int num_node_vars,
			     int num_elem_vars,
			     const Real* donar_node_vars,
			     const Real* donar_elem_vars,
			     Real* receiver_node_vars,
			     Real* receiver_elem_vars,
			     Real* volume_fraction );
  
  // Restart Functions
  virtual int Number_General_Restart_Variables() 
    {return Number_Base_General_Restart_Variables();};
  virtual int Number_Nodal_Restart_Variables();
  virtual int Number_Edge_Restart_Variables() {return 0;};
  virtual int Number_Face_Restart_Variables() {return 0;};
  virtual int Number_Element_Restart_Variables() {return 0;};
  virtual ContactSearch::ContactErrorCode
    Extract_General_Restart_Variable( Real* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable( int n, Real* data, int* node_ids );
  virtual ContactSearch::ContactErrorCode
    Extract_Edge_Restart_Variable( int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Face_Restart_Variable( int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Element_Restart_Variable( int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_General_Restart_Variable( Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable( int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Edge_Restart_Variable( int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Face_Restart_Variable( int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Element_Restart_Variable( int n, Real* data );
  
  // restart functions
  virtual int Restart_Size();
  virtual ContactSearch::ContactErrorCode Extract_Restart_Data(Real* restart_data );
  
  virtual ContactSearch::ContactErrorCode
    Update_For_Topology_Change( int new_number_of_nodes, int* old_to_new_map );

  virtual void Output_Timings() {};

 private:
  
  ContactVolumeTransfer(ContactVolumeTransfer&);
  ContactVolumeTransfer& operator=(ContactVolumeTransfer&);
  
  Real* final_position;
#ifndef CONTACT_NO_MPI
  int* recv_node_var_set;
#endif

  ContactSearch::ContactErrorCode Volume_Transfer_Set_Up();
  void Volume_Transfer_Clean_Up();
  void Set_CVARS();

  // These are place holders to keep from passing the host code data around
  
  // Internal memory
  
  //inlines
  
  // variable plotting
  bool SAVE_CVARS;
  // regression test plot tests
  bool SAVE_NVARS;
  bool SAVE_EVARS;
  void Set_EVARS( int num_elem_vars, Real* element_data );
  void Set_NVARS( int num_elem_vars, Real* nodal_data );
  
  // Enforcement Data (near search tolerance)
  double Get_Near_Tolerance ( ContactElementElementInteraction* );
  void   Display_Enforcement_Data();
  double Enforcement_Data( Enforcement_Data_Index, int, int );
};

#endif
