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


#include "ContactAnalyticSurface.h"

ContactAnalyticSurface::ContactAnalyticSurface( 
    ContactSearch::AnalyticSurface_Type Type, int Host_Index, int Entity_Key )
  : ContactTopologyEntity<Real>( 0, Host_Index, NULL, CT_ANALYTIC_SURFACE)
{
  type = Type;
  entity_key = Entity_Key;
  ProcArrayIndex(Host_Index);
  Global_ID(0,Host_Index);
}

ContactAnalyticSurface::ContactAnalyticSurface( 
    ContactSearch::AnalyticSurface_Type Type, int Host_Index)
  : ContactTopologyEntity<Real>( 0, Host_Index, NULL, CT_ANALYTIC_SURFACE)
{
  type = Type;
  entity_key = -1;
  ProcArrayIndex(Host_Index);
  Global_ID(0,Host_Index);
}

ContactAnalyticSurface::ContactAnalyticSurface( 
    ContactSearch::AnalyticSurface_Type Type, 
    ContactAnalyticSurface* surf)
  : ContactTopologyEntity<Real>( 0, surf->host_array_index, NULL, CT_ANALYTIC_SURFACE)
{
  type = Type;
  entity_key = surf->entity_key;
  
  BlockID(surf->block_id);
  HostArrayIndex(surf->host_array_index);
  ProcArrayIndex(surf->proc_array_index);
  Global_ID(surf->global_id.HiInt(),surf->global_id.LoInt());
}

ContactAnalyticSurface::~ContactAnalyticSurface() {}

void
ContactAnalyticSurface::Display(ContactParOStream& postream)
{
  postream<<"ContactEntity: "<<global_id<<"\n";
  postream<<"               entity type:       CT_ANALYTIC_SURFACE\n";
  postream<<"               surface type:      ";
  switch (type) {
  case ContactSearch::PLANE:
    postream<<"PLANE\n";
    break;
  case ContactSearch::SPHERE:
    postream<<"SPHERE\n";
    break;
  case ContactSearch::CYLINDER_INSIDE:
    postream<<"CYLINDER_INSIDE\n";
    break;
  case ContactSearch::CYLINDER_OUTSIDE:
    postream<<"CYLINDER_OUTSIDE\n";
    break;
  default:
    postream<<"\n";
    break;
  }
  postream<<"               surface_id:        "<<proc_array_index<<"\n";
}

