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
#include "zoltan.h"

extern "C" {
  ZOLTAN_NUM_OBJ_FN            ContactQueryNumObjects;
  ZOLTAN_OBJ_LIST_FN           ContactQueryObjectList;
  ZOLTAN_NUM_GEOM_FN           ContactQueryNumGeomObjects;
  ZOLTAN_GEOM_MULTI_FN         ContactQueryGeomMultiValues;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactMigrateEntityExportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactMigrateEntityExportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactMigrateEntityUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactMigrateInteractionSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactMigratePackInteractions;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactMigrateUnpackInteractions;
  
  ZOLTAN_OBJ_SIZE_FN           ContactHostidQuerySize;
  ZOLTAN_PACK_OBJ_FN           ContactHostidQueryPack;
  ZOLTAN_UNPACK_OBJ_FN         ContactHostidQueryUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactDynamicLoadBalanceSize;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactDynamicLoadBalancePack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactDynamicLoadBalanceUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactGhostingExportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactGhostingExportPack;
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactGhostingImportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactGhostingImportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactGhostingUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactUpdateGhostingExportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactUpdateGhostingExportPack;
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactUpdateGhostingImportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactUpdateGhostingImportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactUpdateGhostingUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactUpdateTiedImportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactUpdateTiedImportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactUpdateTiedUnpack;

}

