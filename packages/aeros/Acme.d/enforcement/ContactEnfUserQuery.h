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


extern "C" {
  
  void FORTRAN(userquery_table_last_abscissa)(void *obj, int *id, Real* value);
  void FORTRAN(userquery_table_interpolate)(void *obj, int *id, 
                                            Real* abscissa, Real* ordinate);
  void FORTRAN(userquery_number_of_nodes)(void *obj, int *nnodes);
  void FORTRAN(userquery_node_state_data)(void *obj, int* id, int* node, 
                                          int* offset, Real* value);
  void FORTRAN(userset_node_state_data)(void *obj, int* id, int* node,
                                        int* offset, Real* value);
  void FORTRAN(userset_nfi_failure_tied)(void *enf_obj, int* model_id, 
                                         void *nfi_obj);

}
