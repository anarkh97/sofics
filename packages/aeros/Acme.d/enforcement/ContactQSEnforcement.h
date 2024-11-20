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


#ifndef ContactQSEnforcement_h_
#define ContactQSEnforcement_h_

#include "ContactSearch.h"
#include "ContactEnforcement.h"
#include "ContactTDEnfModel.h"

template<typename DataType> class ContactNode;
class ContactNodeFaceInteraction;
class ContactNodeSurfaceInteraction;
class ContactTopology;

class ContactQSEnforcement : public ContactEnforcement {
  
 public:

  enum Enforcement_Data_Index { KINEMATIC_PARTITION=0, 
				FRICTION_MODEL_ID,
                                NSIZED };

  ContactQSEnforcement( const double*, ContactSearch*,
			ContactSearch::ContactErrorCode& error );
  ContactQSEnforcement( ContactSearch*, const double* restart_data,
			ContactSearch::ContactErrorCode& error );
  ~ContactQSEnforcement();

// NEED : capture_tolerance, tension_release
//        penalty/s

  /** \name Modify_Predictor_Velocity
   * a.Retrieves the active constraints
   * b.Applies the rate constraint for the normal constraints
   * c.Allows an amount of slip:
   * slip_predictor * slip_implied_in_the_predictor ??
   *  usually zero
  */
  void Modify_Predictor_Velocity(Real *PredictedVelocity, const Real delta_time);


  void Modify_Preconditioner(Real *Preconditioner);

  void Remove_Gaps(Real *Position);

  
  /** \name Compute_Contact_Force
   * a.Retrieves the active constraints
   * b.Given the current residual vector, compute contact force
   * c.Determines if a node is sticking or slipping.
   *   if slipping it leaves a tangential residual on the slave node
   *   (1-slip_penalty)*(tangential_residual)
   */
  void Compute_Contact_Force(Real *Residual);

  /** \name Compute_Consistent_Tied
   * a. retreives the consistent tying constraints
   * b. given current configuration, computes contact force and updates 
   *    multipliers
   * c. updates consistent tying constraints i.e. the operator
   */
  void Compute_Consistent_Tied(void);

  /** \name Update_Search_Direction
   * a. retrieves the active contact constraints
   * b. applies normal contact rate constraint
   * c. applies tangential contact rate constraint if a node is sticking.
   *    adds tangential motion on slave node if a node is slipping. 
   */
  void Update_Search_Direction(Real *SearchDirection);

  /** \name Compute_Slip
   * a. retrieves the active contact constraints
   * b. computes and stores relative slip at slave nodes
   */
  void Compute_Slip(void);

  /** \name Compute_Slip_Residual
   * takes active constraints and computes L2 normal of residual...
   */
  void Compute_Slip_Residual(void);

  /** \name Advance_Contact_State
   * a. retrieves the active contact constraints
   * b. makes permanent slip increment and state variables for load step
   * c. gathers data stored with contact constraints and scattes them to 
   *    nodal arrays
   */
  void Advance_Contact_State(void);

  void Compute_Consistent_Tied_Control(void);
  void Compute_FC_Controlled(void);
  void Update_Search_Direction_Control(void);


  // external interface functions 
  void Get_Old_Displacement( int, Real* );
  void Get_Old_Disp_Restart( int, Real* );
  void Set_Old_Disp_Restart( int, Real* );

  // Restart Functions (DISABLED presently)
  virtual int Number_General_Restart_Variables() {return 0;}; //
  virtual int Number_Nodal_Restart_Variables(){return 0;}; //
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

  // restart functions
  virtual int Restart_Size() {return 0;}; // 
  virtual ContactSearch::ContactErrorCode 
    Extract_Restart_Data(Real* restart_data )
    {return ContactSearch::NO_ERROR;};

  // internal debugging & regression test output functions (DISABLED presently)
  virtual Real Get_Global_Plot_Variable( int ) {return 0.0;};
  virtual void Get_Nodal_Plot_Variable( int, Real* ) {return;};
  virtual void Get_Element_Plot_Variable( int, Real* ) { return; };

 private:
  int max_interactions;
  int num_iterations;
  Real convergence_tolerance;

};

#endif

