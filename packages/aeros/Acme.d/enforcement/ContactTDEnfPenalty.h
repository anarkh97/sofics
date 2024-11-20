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


#ifndef ContactTDEnfPenalty_h_
#define ContactTDEnfPenalty_h_

#include "ContactSearch.h"
#include "ContactEnforcement.h"
#include "ContactTDEnforcement.h"
#include "ContactTDEnfModel.h"

template<typename DataType> class ContactNode;
class ContactNodeFaceInteraction;
class ContactNodeSurfaceInteraction;
class ContactNodeEntityInteraction;
class ContactTopology;

class ContactTDEnfPenalty : public ContactTDEnforcement {
  
 public:
  ContactTDEnfPenalty( const double*, ContactSearch*,
 			     ContactSearch::ContactErrorCode& error,
                             bool get_cvars, 
                             bool get_plot_force );
  ContactTDEnfPenalty( ContactSearch*, const double* restart_data,
			     ContactSearch::ContactErrorCode& error,
                             bool get_cvars,
                             bool get_plot_force );
  ~ContactTDEnfPenalty();

  ContactSearch::ContactErrorCode Compute_Contact_Force( double dt_old, 
							 double dt,
							 double* mass,
							 double* density,
							 double* wavespeed,
							 double* force );
   ContactSearch::ContactErrorCode Set_Penalty_Scale( double scale );

 private:
  void TDPenalty_Release_Scratch(void);

  // Nodal Scratch Variable Handles
/* PJSA use base class
  ScratchVariable NODAL_MASS;
  ScratchVariable INC_FORCE;
  ScratchVariable TOTAL_FORCE;
*/
//  ScratchVariable NEW_POSITION;

  ContactNodeFaceInteraction* cnfi;
  ContactNodeEntityInteraction* cnei;
  double penalty_scale;

};

#endif

