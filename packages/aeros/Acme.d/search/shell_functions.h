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


// These functions override base class virtual functions for shells.
// I chose to do this with an include file instead of having a 
// diamond multiple inheritance.


 public: 

#ifndef CONTACT_NO_MPI
  virtual int Size() { return (ContactFace<DataType>::Size() + sizeof(Real)); };
  virtual void Pack( char* buffer ){ 
    ContactFace<DataType>::Pack(buffer);
    std::memcpy( buffer+ContactFace<DataType>::Size(), &thickness, sizeof(Real) );
  };
  virtual void Unpack( char* buffer ){
    ContactFace<DataType>::Unpack(buffer);
    std::memcpy( &thickness, buffer+ContactFace<DataType>::Size(), sizeof(Real) );
  };
#endif

  Real Thickness() { return thickness; };
  void Thickness(Real t) {thickness = t;};

  Real Lofting_Factor() { return lofting_factor; };
  void Lofting_Factor( Real lf ) { lofting_factor = lf; };

 private:

  Real thickness;
  Real lofting_factor;

