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


#ifndef ContactShellTriFaceL3_h_
#define ContactShellTriFaceL3_h_

#include "ContactTriFaceL3.h"
#include <cstring>

/* This class represents the linear three node triangle shell face with the
   following fortran numbering convention.

                              2
                              /\
                             /  \
                            /    \
                        E2 /      \ E1
                          /        \
                         /          \
                        0------------1
                              E0
*/


template<typename DataType>
class ContactShellTriFaceL3 : public ContactTriFaceL3<DataType> {

  public :
    ContactShellTriFaceL3( ContactFixedSizeAllocator*, 
                           int blk_indx=-1, int indx_in_block=-1, int key=-1);
  static 
    ContactShellTriFaceL3<DataType>* new_ContactShellTriFaceL3(ContactFixedSizeAllocator*,
                                                     int blk_indx=-1, 
                                                     int indx_in_block=-1, 
                                                     int key=-1);
  ~ContactShellTriFaceL3();

 private:

#include "shell_functions.h"

};

#endif // ContactShellTriFaceL3_h_
