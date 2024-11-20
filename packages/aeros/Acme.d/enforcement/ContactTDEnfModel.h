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


#ifndef ContactTDEnfModel_h_
#define ContactTDEnfModel_h_

#include "ContactEnfModel.h"
#include <cmath>

#define NDIM 3
#define ZERO_TOL 1.0e-10
class ContactParOStream;

class ContactTDEnfModel : public ContactEnfModel {

 public:
  ContactTDEnfModel( int ID, 
		     ContactEnforcement::Enforcement_Model_Types Type,
		     ContactTopology* Topology );
  ContactTDEnfModel(ContactEnforcement::Enforcement_Model_Types Type,
		    ContactTopology* Topology );
  virtual ~ContactTDEnfModel();


  void Set_Node_State_Data(Real value, int node, int offset, int state=0);
  Real Get_Node_State_Data(int node, int offset, int state=0);

  // Virtual functions required by the base class
  virtual int Num_Interaction_State_Variables() = 0;
  virtual int Num_Node_State_Variables() = 0;
  virtual ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
		      int num_tables, ContactTable** tables ) = 0;
  // virtual void Initialize_Interaction_State_Data( Real* data ) = 0;
  virtual void Initialize_Node_State_Data( Real* data ) = 0;
  virtual int Restart_Size() = 0;
  virtual int Extract_Restart_Data( Real* restart_data ) = 0;
  virtual int Implant_Restart_Data( const Real* restart_data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable( int n, const Real* data ) = 0;

  // Virtual functions specific to TD EnfModel Models
  virtual int  Interaction_Type(ContactNodeEntityInteraction*)=0 ;
  virtual bool Active_Interaction(ContactNodeEntityInteraction*, Real);
  /** typically the direction of the tangential component of the contact
   * force is not changed by these models e.g. Coulomb friction, since
   * the incoming predicted force reflects the contributions of many
   * constraints not just the one at the node */
  virtual bool Limit_Force(ContactNodeEntityInteraction*,
	Real, Real*, Real*,Real*,Real, Real, Real*) = 0;

  virtual void Initialize_for_Time_Step(){};

  bool Needs_Glued_Search(int itype)
    { if (itype == TDEM_TIED || itype == TDEM_SPRINGY) return true;
	  else return false;}

  bool Needs_CPP_Search(int itype)
    { if (itype == TDEM_ADHESIVE || itype == TDEM_ADHESIVE_FRICTION)return true;
	  else return false;}

  // NOTE : if this changes SIMOD needs to know
  enum TDEnfModel_Type {
		  TDEM_FRICTIONLESS=0,
		  TDEM_FRICTIONAL,
		  TDEM_TIED,
		  TDEM_ADHESIVE,
		  TDEM_SPRINGY,
		  TDEM_ADHESIVE_FRICTION,
		  TDEM_NUM_MODEL_TYPES
  };


 protected:
  /* inline */
  inline void   Zero(double*);
  inline double Magnitude(double*);
  inline double Dot(double*,double*);
  inline void   Scale(double*,double,double*);
  inline void   Normalize(double*);
  inline void   Remove_Component(double*,double*);
  inline void   Remove_Component(double*,double*,double*);
  inline void   Add(double*,double,double*,double*);
  inline void   Linear_Combination(double,double*,double,double*,double*);

#if 0
  ContactParOStream& ParOStream();
#endif


 private:

};

inline void
ContactTDEnfModel::Zero
(double* v)
{ int i;
  for (i=0;i<NDIM; ++i) {v[i]=0;}
}


inline double
ContactTDEnfModel::Magnitude
(double* v)
{ int i;
  double mag= 0.0;
  for (i=0;i<NDIM; ++i) {mag += v[i]*v[i];}
  mag = std::sqrt(mag);
  return mag;
}

inline double
ContactTDEnfModel::Dot
(double* v1, double* v2)
{ int i;
  double dot= 0.0;
  for (i=0;i<NDIM; ++i) {dot += v1[i]*v2[i];}
  return dot;
}

inline void
ContactTDEnfModel::Scale
(double* v1, double scale, double* v2)
{ int i;
  for (i=0;i<NDIM; ++i) {v2[i] = scale*v1[i];}
}

inline void
ContactTDEnfModel::Normalize
( double* v)
{  int j;
   const Real tol = 1.E-30;
   double mag = 0.0;
   for (j=0;j<NDIM; ++j){ mag +=  v[j]*v[j]; }
   if (mag > tol) {
     mag = 1.0/std::sqrt(mag);
     for (j=0;j<NDIM; ++j){ v[j] =  v[j]*mag; }
   }
}

inline void
ContactTDEnfModel::Remove_Component
( double* v, double* nm )
{  int j;
   double dot_n = 0.0;
   for (j=0;j<NDIM; ++j){ dot_n +=  v[j]*nm[j]; }
   for (j=0;j<NDIM; ++j){ v[j]  -= dot_n*nm[j]; }
}

inline void
ContactTDEnfModel::Remove_Component
( double* v, double* nm, double* v2 )
{  int j;
   double dot_n = 0.0;
   for (j=0;j<NDIM; ++j){ dot_n +=  v[j]*nm[j]; }
   for (j=0;j<NDIM; ++j){ v2[j] =   v[j] - dot_n*nm[j]; }
}

inline void
ContactTDEnfModel::Add
(double* v1, double a2, double* v2, double* v )
{  int j;
   for (j=0;j<NDIM; ++j){ v[j] = v1[j] + a2*v2[j]; }
}

inline void
ContactTDEnfModel::Linear_Combination
(double a1, double* v1, double a2, double* v2, double* v )
{  int j;
   for (j=0;j<NDIM; ++j){ v[j] = a1*v1[j] + a2*v2[j]; }
}

#endif
