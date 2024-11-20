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


#ifndef ContactEntity_H_
#define ContactEntity_H_

#include "Contact_Defines.h"
#include "contact_assert.h"
#include <cstring>

  enum ContactType { CT_UNKNOWN, CT_NODE, CT_SHELL_NODE,
		     CT_EDGE, CT_FACE, CT_ELEM, 
		     CT_ELEMENT, CT_ANALYTIC_SURFACE, 
                     CT_NFI, CT_NSI, CT_NNI,
                     CT_FFI, CT_FCI, CT_EEI,
		     CT_NUM_ENTITY_TYPES };
  
#endif  // ifdef ContactEntity_H_
