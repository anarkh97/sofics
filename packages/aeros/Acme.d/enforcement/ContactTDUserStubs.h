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



#include "Contact_Defines.h"
#include "ContactTDUserSubTypes.h"

extern "C" {
  
  CONTACT_INIT_MODEL_FN                 FORTRAN(user_initialize_model);
  CONTACT_INIT_TIME_STEP_FN             FORTRAN(user_initialize_time_step);
  CONTACT_INIT_NODE_STATE_DATA_FN       FORTRAN(user_init_node_state_data);
  CONTACT_LIMIT_FORCE_FN                FORTRAN(user_limit_force);
  CONTACT_INTERACTION_ACTIVE_FN         FORTRAN(user_is_active);
  CONTACT_INTERACTION_TYPE_FN           FORTRAN(user_interaction_type);
  
}

