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

#ifndef ContactScratchManager_h_ 
#define ContactScratchManager_h_

#include "ContactTopologyEntity.h"
#include "ContactInteractionEntity.h"
//
//  An extended array class for dealing with scratch variables.  All data is contiguous.  
//  Stores the data size and stride for ease of data extraction and bounds checking.
//

class ScratchVariable {
 public:
  enum {DO_NOT_ZERO_SCRATCH,
        ZERO_SCRATCH};

  ScratchVariable();
  ~ScratchVariable();

  void Allocate_Scratch(const int &num_objects_, const int &size_per_object_, const int zero_scratch = DO_NOT_ZERO_SCRATCH);

  void Clear_Scratch();

  void Zero_Scratch();

  void Duplicate_Scratch(ScratchVariable &source) {
    PRECONDITION(is_allocated);
    PRECONDITION(num_objects == source.num_objects);
    PRECONDITION(size_per_object = source.size_per_object);
    std::memcpy(data_array, source.data_array, num_objects * size_per_object * sizeof(Real));
  }

  inline Real* Get_Scratch(const int &object_num) const {
    PRECONDITION(is_allocated);
    PRECONDITION(object_num >= 0 && object_num < num_objects);
    return data_array + (object_num * size_per_object);
  }

  inline Real* Get_Scratch(const ContactTopologyEntity<Real> *entity) {
    return Get_Scratch(entity->EnfArrayIndex());
  }

  inline Real* Get_Scratch(const ContactInteractionEntity<Real> *entity) {
    return Get_Scratch(entity->EnfArrayIndex());
  }


  inline int Get_Size() const {return size_per_object;}

  
  inline Real *Get_Data() const {return data_array;}

 private:
  bool is_allocated;
  int num_objects;
  int size_per_object;
  Real *data_array;
};



#endif
