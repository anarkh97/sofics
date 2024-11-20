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


#ifndef contact_assert_h_
#define contact_assert_h_

#ifdef CONTACT_DEBUG

#include <cstdio>
#include <cstdlib>

#ifndef PRECONDITION
#define PRECONDITION(e) \
  (!(e)? std::fprintf(stderr,"Precondition not satisfied: file \"%s\",\
  line %d\n%s\n", __FILE__, __LINE__, #e), (void)std::abort() : \
  (void)0)
#endif

#ifndef POSTCONDITION
#define POSTCONDITION(e) \
  (!(e)? std::fprintf(stderr,"Postcondition not satisfied: file \"%s\",\
  line %d\n%s\n", __FILE__, __LINE__, #e), (void)std::abort() : \
  (void)0)
#endif

#ifndef INVARIANT
#define INVARIANT(e) \
  (!(e)? std::fprintf(stderr,"Postcondition not satisfied: file \"%s\",\
  line %d\n%s\n", __FILE__, __LINE__, #e), (void)std::abort() : \
  (void)0)
#endif

#ifndef REMEMBER
#define REMEMBER(e) e
#endif

#else

#ifndef PRECONDITION
#define PRECONDITION(e) 
#endif

#ifndef POSTCONDITION
#define POSTCONDITION(e) 
#endif

#ifndef INVARIANT
#define INVARIANT(e) 
#endif

#ifndef REMEMBER
#define REMEMBER(e)
#endif

#endif // #ifdef CONTACT_DEBUG

#endif // #ifndef contact_asserts_h_


