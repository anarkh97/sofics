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


#ifndef _ContactEntityHash_
#define _ContactEntityHash_

#include "ContactTopologyEntity.h"
#include "ContactParOStream.h"
#include "ContactBlockEntityList.h"
#include <iostream>
#include <vector>

namespace topology_hash_0
{
  struct implementation {
    implementation();
    implementation(int, ContactTopologyEntity<Real>**);
    implementation(std::vector<ContactTopologyEntity<Real>*> *);
    ~implementation();

    std::ostream& stream_data(std::ostream& os) const;
    ContactParOStream& stream_data(ContactParOStream& os) const;

    void ReHash(int, ContactTopologyEntity<Real> **);
    void SetupHash(int, ContactTopologyEntity<Real> **);
    void ClearHash();
    ContactTopologyEntity<Real> *find( ContactHostGlobalID & );
    ContactTopologyEntity<Real>* find( ContactInteractionEntity<Real>::entity_data *data );
    void insert(ContactHostGlobalID &Global_ID, ContactTopologyEntity<Real>* Entity);

    int nbins;
    int nbins_orig;
    int number_of_entities;

    struct hash {
      hash() : entity(0), next(0) {}
      ContactHostGlobalID global_id;
      ContactTopologyEntity<Real>* entity;
      struct hash *next;
    };

    hash* hash_space;
    hash** bins;
    hash* next_free;



#ifdef CONTACT_ANALYZE_HASH
    int hash_collisions;
#endif

    int  hash_func(int);
    void create_space();
    void ComputeNbins(int);
  };
}

#ifdef CONTACT_HAVE_COMPILER_HASH
#include <ext/hash_map>
namespace topology_hash_1
{
  using __gnu_cxx::hash_map;
  struct implementation {
    implementation() {}
    implementation(int, ContactTopologyEntity<Real> **);
    implementation(std::vector<ContactTopologyEntity<Real>*> *);
    ~implementation() {}

    int hash_func(int) {return 0;}
    void create_space() {}
    void ComputeNbins(int) {}
    void ReHash(int, ContactTopologyEntity<Real> **);
    void SetupHash(int, ContactTopologyEntity<Real> **);
    void AddEntities(ContactBlockEntityList *);
    void ClearHash();
    ContactTopologyEntity<Real> *find( ContactHostGlobalID & );
    ContactTopologyEntity<Real>* find( ContactInteractionEntity<Real>::entity_data *data );
    void insert(ContactHostGlobalID &, ContactTopologyEntity<Real> *);

    std::ostream& stream_data(std::ostream &) const;
    ContactParOStream& stream_data(ContactParOStream &) const;

  private:
    hash_map<int, ContactTopologyEntity<Real> *> container;
  };
}
#endif

class ContactTopologyEntityHash : public topology_hash_0::implementation {

public:
  ContactTopologyEntityHash( );
  ContactTopologyEntityHash( int, ContactTopologyEntity<Real>** );
  ContactTopologyEntityHash( std::vector<ContactTopologyEntity<Real>*> *list);
  ~ContactTopologyEntityHash();

  friend std::ostream& operator<<( std::ostream &, const ContactTopologyEntityHash &);
  friend ContactParOStream& operator<<( ContactParOStream &, const ContactTopologyEntityHash &);

private:
  std::ostream& stream_data(std::ostream& os) const;
  ContactParOStream& stream_data(ContactParOStream& os) const;

};

#endif // _ContactEntityHash_
