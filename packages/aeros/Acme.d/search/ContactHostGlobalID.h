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


#ifndef ContactHostGlobalID_h_
#define ContactHostGlobalID_h_
#ifndef CONTACT_NO_MPI
#include <zoltan.h>
#endif
#ifndef CONTACT_NO_MPI
#include "ContactZoltanID.h"
#endif
#include "Contact_Defines.h"

#include <iostream>

class ContactParOStream;

class ContactHostGlobalID {

 public:
  ContactHostGlobalID( int hi, int lo);
#ifndef CONTACT_NO_MPI
  ContactHostGlobalID( LB_ID_PTR );
#endif
  ContactHostGlobalID();
  ~ContactHostGlobalID();
  
  bool operator> ( const ContactHostGlobalID& ) const;
  bool operator< ( const ContactHostGlobalID& ) const;
  bool operator==( const ContactHostGlobalID& ) const;
  bool operator<=( const ContactHostGlobalID& ) const;
  void operator= ( const ContactHostGlobalID& );
  int  operator% ( const int& ) const;

  inline void Set( int hi, int lo ) { hi_int = hi; lo_int = lo; };

  inline int LoInt() const { return lo_int; } ;
  inline int HiInt() const { return hi_int; } ;

  inline void LoInt( int lo ) { lo_int = lo; };
  inline void HiInt( int hi ) { hi_int = hi; };

  friend std::ostream& operator<<( std::ostream& os, const ContactHostGlobalID& gid );
  friend ContactParOStream& operator<<( ContactParOStream& pos, 
					const ContactHostGlobalID& gid );

 private:

  int lo_int;
  int hi_int;
};

class GIDSortClass {
 public:
  ContactHostGlobalID key;
  int old_index;
};

inline ContactHostGlobalID::ContactHostGlobalID() 
  : lo_int(-1), hi_int(-1) {}

inline ContactHostGlobalID::ContactHostGlobalID( int hi, int lo )
  : lo_int(lo), hi_int(hi) {}

#ifndef CONTACT_NO_MPI
inline ContactHostGlobalID::ContactHostGlobalID( LB_ID_PTR zid )
{
  ContactZoltanGID zoltan_gid;
  lo_int = zoltan_gid.Lo(zid);
  hi_int = zoltan_gid.Hi(zid);
}
#endif

inline ContactHostGlobalID::~ContactHostGlobalID(){}

inline bool 
ContactHostGlobalID::operator>( const ContactHostGlobalID& gid ) const
{
  return lo_int != gid.lo_int ? lo_int > gid.lo_int : hi_int > gid.hi_int;
  
  // Sort first by lo_int and then by hi_int
  //if( lo_int > gid.lo_int )
  //  return true;
  //else if( lo_int < gid.lo_int )
  //  return false;
  //else if( hi_int > gid.hi_int )
  //  return true;
  //return false;
}

inline bool 
ContactHostGlobalID::operator<( const ContactHostGlobalID& gid ) const
{
  return lo_int != gid.lo_int ? lo_int < gid.lo_int : hi_int < gid.hi_int;
  
  // Sort first by lo_int and then by hi_int
  //if( lo_int < gid.lo_int )
  //  return true;
  //else if( lo_int > gid.lo_int )
  //  return false;
  //else if( hi_int < gid.hi_int )
  //  return true;
  //return false;
}

inline bool 
ContactHostGlobalID::operator<=( const ContactHostGlobalID& gid ) const
{
  return lo_int != gid.lo_int ? lo_int < gid.lo_int : hi_int <= gid.hi_int;
  
  // Sort first by lo_int and then by hi_int
  //if( lo_int < gid.lo_int ) {
  //  return true;
  //} else if(lo_int > gid.lo_int) {
  //  return false;
  //} else {
  //  if(hi_int <= gid.hi_int) {
  //    return true;
  //  } else {
  //    return false;
  //  }
  //}
}

inline bool 
ContactHostGlobalID::operator==( const ContactHostGlobalID& gid ) const
{
  if( gid.lo_int == lo_int && gid.hi_int == hi_int ) return true;
  return false;
}

inline void 
ContactHostGlobalID::operator=( const ContactHostGlobalID& gid )
{
  lo_int = gid.lo_int;
  hi_int = gid.hi_int;
}

inline int 
ContactHostGlobalID::operator%( const int& divisor ) const
{
  if (lo_int<0) return(0);
  else return ( lo_int % divisor );
}

#endif
