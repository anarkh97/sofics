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


#ifndef ContactTDUserSubTypes_h_
#define ContactTDUserSubTypes_h_

typedef int CONTACT_INIT_MODEL_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &istat    // return code
);

typedef int CONTACT_INIT_TIME_STEP_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &istat    // return code
);

typedef void CONTACT_INIT_NODE_STATE_DATA_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  Real *sdata,   // node state data
  int  &istat    // return code
);

typedef void CONTACT_NUM_NODE_STATE_VARS_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &istat    // return code
);

typedef int CONTACT_LIMIT_FORCE_FN(
  void *enf,     // enforcement object
  void *ni,      // interaction object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &type,    // type of interaction, nfi or nsi
  int  &index,   // node index
  Real &gap, 
  Real *rel_disp, 
  Real *slip, 
  Real *normal,
  Real &dt, 
  Real &area, 
  Real *force,
  int  &istat    // return code
);

typedef int CONTACT_INTERACTION_ACTIVE_FN(
  void *enf,     // enforcement object
  void *ni,      // interaction object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &type,    // type of interaction, nfi or nsi
  int  &index,   // node index
  Real &gap,     // gap
  int  &istat    // return code
);

typedef int CONTACT_INTERACTION_TYPE_FN(
  void *enf,     // enforcement object
  void *ni,      // interaction object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &type,    // type of interaction, nfi or nsi
  int  &index,   // node index
  int  &istat    // return code
);

#endif
