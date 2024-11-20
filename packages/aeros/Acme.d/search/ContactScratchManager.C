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


#include "contact_assert.h"
#include "Contact_Defines.h"
#include "ContactScratchManager.h"
#include "ContactTDEnfModel.h"
#include "ContactErrors.h"
#include "ContactTDEnforcement.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "ContactShellNode.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"

//
//  Scratch variable constructor creates an empty scratch variable
//
ScratchVariable::ScratchVariable() :
  is_allocated(false),
  num_objects(0),
  size_per_object(0),
  data_array(NULL)
{}


//
//  Clear operation frees memory and resets all meta-data about the scratch array
//
void ScratchVariable::Clear_Scratch() {
  if(data_array) delete [] data_array;
  data_array = NULL;
  is_allocated = false;
  num_objects = 0;
  size_per_object = 0;
}

//
//  Scratch variable destructor frees all allocated memory
//
ScratchVariable::~ScratchVariable() {
  Clear_Scratch();
}

//
//  Allocate routine creates the actual data array based on the input number of objects and size per object
//
void ScratchVariable::Allocate_Scratch(const int &num_objects_, const int &size_per_object_, const int zero_scratch) {
  PRECONDITION(!is_allocated);
  num_objects = num_objects_;
  size_per_object = size_per_object_;
  PRECONDITION(num_objects >= 0);
  PRECONDITION(size_per_object > 0);
  is_allocated = true;
  data_array = new Real[num_objects * size_per_object];
  if(zero_scratch == ZERO_SCRATCH) Zero_Scratch();
}

//
//  Zero scratch routine fills entire data array with zeros
//
void ScratchVariable::Zero_Scratch() {
  PRECONDITION(is_allocated);
  std::memset(data_array, 0, num_objects * size_per_object * sizeof(Real));
}
