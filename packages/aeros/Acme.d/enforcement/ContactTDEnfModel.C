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


#include "ContactTDEnfModel.h"
#include "ContactEnforcement.h"
#include "ContactEnfModel.h"

ContactTDEnfModel::ContactTDEnfModel( int ID,
         ContactEnforcement::Enforcement_Model_Types Type,
				      ContactTopology* Topology )
  : ContactEnfModel( ID, Type, Topology )
{
}

ContactTDEnfModel::ContactTDEnfModel( 
	 ContactEnforcement::Enforcement_Model_Types Type,
	 ContactTopology* Topology )
  : ContactEnfModel( Type, Topology )
{
}

ContactTDEnfModel::~ContactTDEnfModel()
{}

bool 
ContactTDEnfModel::Active_Interaction
(ContactNodeEntityInteraction* cnfi,Real gap)
{
  return (gap < 0 ? true : false);
}

void 
ContactTDEnfModel::Set_Node_State_Data(Real value, int node, 
                                       int offset, int state)
{
  Node_State_Data(node,state)[offset] = value;
}

Real 
ContactTDEnfModel::Get_Node_State_Data(int node, int offset, int state)
{
  return Node_State_Data(node,state)[offset];
}

#if 0
ContactParOStream& 
ContactTDEnfModel::ParOStream()
{
	  return search->postream;
}
#endif

