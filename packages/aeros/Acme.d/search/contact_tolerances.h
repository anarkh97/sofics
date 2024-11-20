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


#ifndef contact_tolerances_h_
#define contact_tolerances_h_

//---------------------------------------------------------------
//
//        ACME's G L O B A L   T O L E R A N C E S
//
//---------------------------------------------------------------
// see also: search_parameters.par & Contact_Defines.h
// (as well as the individual source files)

// These are the MAX & MIN tolerances for time to contact
// They are used in Interaction_Definition() to determine
// if a moving interaction is valid.
//static Real time_to_contact_tolerance_min = -0.55;
//static Real time_to_contact_tolerance_max =  1.05;
#define time_to_contact_tolerance_min  -0.55
#define time_to_contact_tolerance_max   1.05

// inflation tolerance on bounding boxes: related to JAS's tol_msl & tol_snm
//static Real BOX_INFLATION_FACTOR = 0.01;
#define BOX_INFLATION_FACTOR  0.01

// tolerance that inflates the max remaining gap to get a "reasonable gap" 
//static Real GAP_INFLATION_FACTOR = 0.01;
#define GAP_INFLATION_FACTOR  0.01

// The computed areas in cnodetriange_cpproj for a point exactly on a
// corner (or possibly on an edge) may yield a very small negative number.
// REL_TANG_TOL is used in an if-test to avoid bad consequences from this
// possibility.  A better approach to handle this is desirable. (RMS)
// the tolerance on subtriangle to triangle area in cnodetriangle_cpproj, etc
// NOTE: this was "CTOLT1"
//static Real REL_TANG_TOL  = 1.e-12;
#define REL_TANG_TOL   1.e-12

// The following two tolerances are used in ContactFaceFaceSearch.C
#include <limits>
#define THICK_PLANE_TOL std::numeric_limits<double>::epsilon()
#define SMALL_AREA_TOL 10*std::numeric_limits<double>::epsilon()

// The following tolerances are used in Compute_Local_Coordinates for some ContactElems
#define spatial_tolerance_pre 0
#define spatial_tolerance_post 0
#define newton_tolerance 1.0e-12
#define abs_newton_tolerance 1.0e-29

// The following tolerance is used to define a planar face
#define is_planar_tol 0.999999

// The following tolerance is used to define when a ray is parallel to a plane
#define is_parallel_tol 1e-10

#endif
