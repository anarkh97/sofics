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

#ifdef _WIN32

  void __stdcall FORTRAN(cnodetriangle_cpproj)
                        ( int&, Real*, Real*, Real*, int*,
                          Real*, Real*, Real&, Real&, Real& );

  void __stdcall FORTRAN(cnodetriangle_movsrch)
                        ( int&,  Real*, Real*, Real*, 
                          Real*, Real*, Real*, Real*, Real& );

  void __stdcall FORTRAN(cnodetriangle_cpproj_aug)
                        ( int&, 
                          Real*, Real*, Real*,
                          Real*, Real*, Real*,
                          int*, Real*, Real*, Real&, Real&, Real& );

  void __stdcall FORTRAN(cnodetriangle_movsrch_aug)
                        ( int&,  
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*,
                          Real*, Real&, Real&, Real&);

  void __stdcall FORTRAN(cnodetriangle_cpproj_aug_noauto)
                        ( int&, 
                          Real*, Real*, Real*,
                          Real*, Real*, Real*,
                          int*, Real*, Real*, Real&, Real&, Real& );

  void __stdcall FORTRAN(cnodetri_movsrch_aug_noauto)
                        ( int&,  
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*,
                          Real*, Real&, Real& );

  void __stdcall FORTRANNO_(oq4h8)
                           ( Real*, Real*, int*, int*, 
		             int*, int*, Real*, Real* );

#else

  void FORTRAN(cnodetriangle_cpproj)( int&, Real*, Real*, Real*, int*,
				      Real*, Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetriangle_movsrch)( int&, Real*, Real*, Real*, 
				       Real*, Real*, Real*, Real*, Real& );

  void FORTRAN(cnodetriangle_cpproj_aug)( int&, 
                                          Real*, Real*, Real*,
                                          Real*, Real*, Real*,
                                          int*, Real*, Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetriangle_movsrch_aug)( int&,  
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*,
                                           Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetriangle_cpproj_aug_noauto)( int&, 
                                          Real*, Real*, Real*,
                                          Real*, Real*, Real*,
                                          int*, Real*, Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetri_movsrch_aug_noauto)( int&,  
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*,
                                           Real*, Real&, Real& );

  void FORTRAN(cnodeline_cpproj)( int&, Real*, Real*, Real*, int*,
				  Real*, Real*, Real&, Real&, Real& );

  void FORTRANNO_(oq4h8)( Real*, Real*, int*, int*, int*, int*, Real*, Real* );

#endif

}

