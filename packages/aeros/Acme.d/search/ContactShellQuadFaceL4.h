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


#ifndef ContactShellQuadFaceL4_h_
#define ContactShellQuadFaceL4_h_

#include "ContactQuadFaceL4.h"
#include <cstring>

class ContactFixedSizeAllocator;

/* This class represents the bilinear four node quadrilateral shell face with 
   the following fortran numbering convention.

                             E2
                       3-------------2
                       |             |
                       |             |
                    E3 |             | E1
                       |             |
                       0-------------1
                             E0
*/

template<typename DataType>
class ContactShellQuadFaceL4 : public ContactQuadFaceL4<DataType> {

  public :
    ContactShellQuadFaceL4( ContactFixedSizeAllocator*, 
                            int blk_indx=-1, int indx_in_block=-1, int key=-1);
  static ContactShellQuadFaceL4<DataType>* new_ContactShellQuadFaceL4(
		    ContactFixedSizeAllocator*,
                    int blk_indx=-1, int indx_in_block=-1, int key=-1);

  ~ContactShellQuadFaceL4();

  private :
    
#include "shell_functions.h"

};

#endif // ContactShellQuadFaceL4_h_
