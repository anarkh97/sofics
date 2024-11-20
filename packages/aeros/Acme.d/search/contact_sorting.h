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

#ifdef WIN32

  void __stdcall FORTRAN(contact_make_rank)
                        ( int&, int&, int&, Real*, int*, int*, int* );

  void __stdcall FORTRAN(contact_get_bound)
                        ( int&, Real*, int*, int&, Real*, Real*, 
                          int&, int*, int*, int&, Real* );

  void __stdcall FORTRAN(contact_make_list)
                        ( int&, int*, int*, int*, int*, 
                          int&, int*, int&, int&, int& );

  void __stdcall FORTRAN(contact_make_rank_t)
                        ( int&, int&, int&, Real*, Real*, int*, int*, int* ,
                          int* ,int& );

  void __stdcall FORTRAN(contact_get_bound_t)
                        ( int&, Real*, int*, int&, Real*, Real*, 
                          int&, int*, int*, int&, Real*, int& );

 void __stdcall FORTRAN(contact_indexx)
                        ( int&, Real*, int*, int& );

 void __stdcall FORTRAN(contact_indexx_float)
                        ( int&, float*, int*, int& );

  void __stdcall FORTRAN(contact_rank)
                        ( int&, int*, int*, int& );

#else

  void FORTRAN(contact_make_rank)( int&, int&, int&, Real*, int*, int*, int*);

  void FORTRAN(contact_get_bound)( int&, Real*, int*, int&, Real*, Real*, int&, 
				   int*, int*, int&, Real*);

  void FORTRAN(contact_make_list)( int&, int*, int*, int*, int*, int&, int*, 
				   int&, int&, int& );
				    
  void FORTRAN(contact_make_rank_t)( int&, int&, int&, Real*, Real*, int*, int*, int*,
                                     int*,int&  );

  void FORTRAN(contact_get_bound_t)( int&, Real*, int*, int&, Real*, Real*, int&, 
				     int*, int*, int&, Real*, int& );

  void FORTRAN(contact_indexx)
                        ( const int&, Real*, int*, const int& );

  void FORTRAN(contact_indexx_float)
                        ( const int&, float*, int*, const int& );

  void FORTRAN(contact_rank)
                        ( const int&, int*, int*, const int& );

#endif

}
