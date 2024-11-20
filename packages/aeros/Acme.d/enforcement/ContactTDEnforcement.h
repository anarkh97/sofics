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


#ifndef ContactTDEnforcement_h_
#define ContactTDEnforcement_h_
 
#include "ContactSearch.h"
#include "ContactEnforcement.h"
#include "ContactTDEnfModel.h"
#include "ContactScratchManager.h"

template<typename DataType> class ContactNode;
class ContactNodeFaceInteraction;
class ContactNodeSurfaceInteraction;
class ContactTopology;

class ContactTDEnforcement : public ContactEnforcement {
  
 public:

  enum Enforcement_Data_Index { KINEMATIC_PARTITION=0, 
				FRICTION_MODEL_ID,
                                NSIZED };
  enum Contact_Status { INACTIVE=0, ACTIVE };
  enum Contact_TDEnf_Plot_Vars { CONFACE=1,
                                 NORMAL_FORCE_MAG,
                                 NORMAL_TRACTION_MAG,
				 TANGENTIAL_FORCE_MAG,
                                 TANGENTIAL_TRACTION_MAG,
				 CDIRNORX, CDIRNORY, CDIRNORZ,
                                 CDIRTANX, CDIRTANY, CDIRTANZ,
                                 SLIP_MAG,
                                 NODAL_DISSIPATION,
                                 NODAL_DISSIPATION_DENSITY,
                                 CONTACT_AREA,
                                 GAP_CUR,
                                 GAP_OLD,
                                 KINEMATIC_PARTITION_VALUE};
  enum Contact_TDEnf_Global_Vars {
                                 FORCE_X=1,
                                 FORCE_Y,
                                 FORCE_Z,
                                 FORCE_NORM,
                                 DISSIPATION,
                                 CONSTRAINT_NORM,
                                 INC_FORCE_NORM };

  ContactTDEnforcement( const double*, ContactSearch*,
			ContactSearch::ContactErrorCode& error,
                        bool get_cvars,
                        bool calc_plot_force );
  ContactTDEnforcement( ContactSearch*, const double* restart_data,
			ContactSearch::ContactErrorCode& error,
                        bool get_cvars,
                        bool calc_plot_force );
  virtual ~ContactTDEnforcement();

  virtual ContactSearch::ContactErrorCode Compute_Contact_Force( double dt_old, 
                                                                 double dt,
                                                                 double* mass,
                                                                 double* density,
                                                                 double* wavespeed,
                                                                 double* force );


  ContactSearch::ContactErrorCode 
  Enforce_Symmetry_on_Nodes( int node1_exodus_id, int node_2_exodus_id );

  // internal debugging and regression test output functions
  virtual Real Get_Global_Plot_Variable( int );
  virtual void Get_Nodal_Plot_Variable( int, Real* );
  virtual void Get_Element_Plot_Variable( int, Real* ) { return; };

  // external output of contact information
  ContactSearch::ContactErrorCode Get_Plot_Variable( Contact_TDEnf_Plot_Vars, 
						     Real* );
  ContactSearch::ContactErrorCode Get_Global_Variable( Contact_TDEnf_Global_Vars,
						     Real&); 
 
  // external interface functions
  ContactSearch::ContactErrorCode Set_Number_of_Iterations( int iter );
  ContactSearch::ContactErrorCode Set_Convergence_Tolerance(Real tol);
  void Set_Save_Nodal_Vars (bool flag);
  void Get_Old_Displacement_Ptrs( int, Real** );
  void Get_Old_Displacement( int, Real* );
  void Get_Old_Disp_Restart( int, Real* );
  void Set_Old_Disp_Restart( int, Real* );
  Real Get_Old_Time_Multiplier() { return time_mult_old_aug; }; 

  // Restart Functions
  virtual int Number_General_Restart_Variables() 
    {return 1+Number_Base_General_Restart_Variables();};
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

  virtual ContactSearch::ContactErrorCode
    Update_For_Topology_Change( int new_number_of_nodes,  int* old_to_new_map );

  // restart functions
  virtual int Restart_Size();
  virtual ContactSearch::ContactErrorCode Extract_Restart_Data(Real* restart_data );

  static Real Invert_3x3Matrix(double a[3][3], double i[3][3]);
  
  void Display_Enforcement_Data();
  void Display0_Enforcement_Data();

 protected:
  // These are place holders to keep from passing the host code data around
  double dt;
  double dt_old;
  double dt2;
  double time_mult_old_aug;

  bool only_frictionless_or_tied;

  Real* plot_force;
  Real* disp_old_aug;

  // variable plotting
  bool SAVE_CVARS;
  bool SAVE_CGVARS;
  const bool CALC_PLOT_FORCE;
  Real* conface;

  Real* normal_force_mag;
  Real* normal_traction_mag;
  Real* tangential_force_mag;
  Real* tangential_traction_mag;

  Real* normal_dir_x;
  Real* normal_dir_y;
  Real* normal_dir_z;

  Real* tangential_dir_x;
  Real* tangential_dir_y;
  Real* tangential_dir_z;

  Real* slipmag;
  Real* nodal_dissipation;
  Real* nodal_dissipation_density;
  Real* contact_area;
  Real* plot_gap_cur;
  Real* plot_gap_old;
  Real* kinematic_partition;

  // timing variables
  int enforcement_time;
  int initialization_time;
  int scratch_allocation_time;
  int copy_nodal_to_scratch;
  int init_host_code_vars;
  int setup_interaction_list;
  int td_enforcement_setup;

  int iteration_time;
  int compute_gap_time;
  int penalty_time;
  int assemble_forces_time;

  int accel_cor_time;
  int response_prediction_time;
  int response_correction_time;
  int swap_add_time;
  int response_constitutive_correction_time;
  int contact_force_time;
  int force_adjust_time;


  int kinematic_consistent_time;
  int update_config_time;
  int extra_ghost_len_time;
  int cleanup_time;

  Real  global_force_x, global_force_y, global_force_z, global_force_norm;
  Real  global_dissipation;
  Real  global_constraint_norm;
  Real  global_inc_force_norm, initial_global_inc_force_norm;

  int max_interactions;
  int num_iterations;
  Real convergence_tolerance;

  ContactTDEnfModel* TDEnfModel( ContactNodeEntityInteraction* cnei);
  ContactTDEnfModel* TDEnfModel( int i ){
    return (ContactTDEnfModel*) enforcement_models[i]; };
  bool Only_Frictionless_or_Tied(void);
  void Make_Scratch_Vector_Consistent_with_BCs( ScratchVariable &var );
  void Compute_Extra_Ghosting_Length();
  //
  //  Nodal Scratch Variables
  //
  ScratchVariable NEW_POSITION;  

  void Debug_Dump(int iteration);

 private:
  void Apply_Face_Forces(Real face_forces[3], 
                         Node_Constraint_Group &cnei_group,
                         int inter_num,
			 Real m);                         

  void Verify_Enforcement_Data();
  
  ContactSearch::ContactErrorCode TD_Enforcement_Set_Up();
  ContactSearch::ContactErrorCode Symmetry_Node_Correction();
  void TD_Enforcement_Clean_Up();
  void Compute_Penalty_Forces();
  void Assemble_Nodal_Forces_and_Masses();
  void Assemble_Nodal_Forces_No_Multiple();
  void Assemble_Masses();
  void Response_Prediction(Real mult);
  void Response_Correction(); 
  void Response_Constitutive_Correction();
  double Auto_Kinematic_Partition( ContactNodeFaceInteraction*);
  double Enforcement_Data( Enforcement_Data_Index, int, int );
  double Enforcement_Data( Enforcement_Data_Index, ContactNodeEntityInteraction*);
  void Get_Normals_and_Gaps_at_Node(const Node_Constraint_Group cnei_list, 
                                    double**, double*);
  void Unify(double**, double*, double*);
  void Partition_Mass(int, double**, double*, double*);
  void Partition_Force(int, double**, double*, double*, double[3][3]);
  void Partition_Gap(double**, double*, double* );
  void Compute_Kinematic_Quantities( int iteration );
 protected:
  void Initialize_CVARS();
  void Set_CVARS();
  void Set_CGVARS();
  void Store_Shell_Final_Lofted_Positions();
  void Assemble_Shell_Forces();
 private:

  Real third_vector[3];

  //
  //  Nodal Scratch Variables
  //
 protected: // PJSA: make these available to derived class ContactTDEnfPenalty
  ScratchVariable ASSEMBLED_MASS;
  ScratchVariable NODAL_MASS;
  ScratchVariable NODAL_WAVESPEED;
  ScratchVariable NODAL_DENSITY;
  ScratchVariable ASSEMBLED_FORCE;
  ScratchVariable SN_ACCELERATION;
  ScratchVariable INC_FORCE;
  ScratchVariable A_PREDICTION;
  ScratchVariable A_CORRECTION;
  ScratchVariable A_CONSTITUTIVE;  
  ScratchVariable TOTAL_FORCE;  
 protected:
  void Compact_Node_Entity_List();
  int* number_of_constraints;
  ContactNodeEntityInteraction** node_entity_constraints;

  std::vector<ContactNodeEntityInteraction*> extra_constraints;


  ScratchVariable PREDICTED_POS;
  ScratchVariable CURRENT_POS;
  ScratchVariable NUM_KIN_CONSTR;
  ScratchVariable KIN_CONSTR_VECTOR;
 private:
  //
  //  Common Node-Face / Node-Surface Variables
  //
  ScratchVariable NEI_TOTAL_SLIP;
  ScratchVariable NEI_GAP_TO_ENFORCE;
  ScratchVariable NEI_REL_DISP;
  ScratchVariable NEI_FORCE;
  ScratchVariable NEI_TAN_SLIP;
  ScratchVariable NEI_TARGET_GAP;
  ScratchVariable NEI_KINEMATIC_PARTITION;


  bool have_symmetry_nodes;
  int symmetry_node_1_exodus_id;
  int symmetry_node_2_exodus_id;

  int have_multiple;
  int have_auto_kin_part;

  
  ContactSearch::Search_Option_Status new_tied_enforcement;
};

#endif

