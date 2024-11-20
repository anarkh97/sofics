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



#ifndef ContactRangeSearch_h_
#define ContactRangeSearch_h_

#include <vector>
#include <Contact_Defines.h>
#include <ContactVector.h>
#include <ContactBoundingBoxHierarchy.h>
#include <ContactNodeFaceInteraction.h>
#include <ContactSearch.h>
#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif

class ContactTopology;
template<typename DataType> class ContactNode;
class ContactSearch;
//
//  This range search class holds data necessary for doing a node on face range searches.
//  The search can be initialized with a variety of options which are then used to create
//  build and perform the search.
//
class ContactNodeFaceRangeSearch {
 public:
  //
  //  Specifies the type of range search to perform.  All should get the same awnser.
  //  However the performance of these options will vary.
  //
  enum SearchType{BASIC_SEARCH, KEY_AWARE_SEARCH, ENTITY_SEARCH};

  ContactNodeFaceRangeSearch(SearchType my_type, 
                             ContactTopology *my_topo,
                             const int num_active_nodes_,
                             int *node_search_status,
                             ContactSearch *my_search,
                             VariableHandle &NODE_COORD_START_,
                             VariableHandle &NODE_COORD_END);
                             
  ContactNodeFaceRangeSearch(SearchType my_type, 
                             ContactTopology *my_topo,
                             const int num_active_nodes_,
                             int *node_search_status,
                             ContactSearch *my_search,
                             VariableHandle &NODE_COORD_START_,
                             VariableHandle &NODE_COORD_END,
                             Real capture_tol);

  ~ContactNodeFaceRangeSearch();
  //
  //  Determine if a search using this search object can be sucessefully performed
  //
  bool valid_search();

  void search_for_overlap(ContactBoundingBox &face_search_box,
                          ACME::Int_Vector &node_keys,
                          ACME::ContactNode_Vector &node_list,
                          std::vector<bool> &valid_inter);

  void process_analytic_surfaces(ContactNodeFaceInteraction::InteractionSource Process_Method,
                                 VariableHandle POSITION, VariableHandle POSITION_2 = -1);

 private:
  //
  //  This data is used to perform a variety of range searches.  Not all of it is used for any given search type
  // 
  ObjectBoundingBoxHierarchy *node_hierarchy;
  ObjectBoundingBoxHierarchy *node_block_BB_hierarchy;
  ObjectBoundingBoxHierarchy **node_hierarchy_ptrs;
  int *node_key_array;
  Real *scratch;
  int *index;
  int *rank;
  int *rank2;
  Real *position;
  int *node_map;  

  Real *max_remaining_gap;
  Real *max_node_motion;
  ContactNode<Real> **Nodes;
  std::vector<int> node_entity_keys;
  ContactTopology *search_topology; 
  ContactSearch *search;
  SearchType search_type;
  int num_active_nodes;
  VariableHandle NODE_COORD_START; 
  VariableHandle NODE_COORD_END;
  VariableHandle NODE_NORMAL;
  VariableHandle FACE_NORMAL;
  VariableHandle REMAINING_GAP;
};


//
// These classes and functions provide an interface to a number of Generic range search functions.  
// These functions can be called by the host code directly to peform geometric searches without the overhead
// of creating an ACME search object or topology.  Additionally these routines may be used as workhorse routines 
// by the ACME search objects to find likely interactions based on purely geometric information.
//
namespace ACME {
  //
  //  Pick a range search type based on various input parameters.  All searches yield identical results,
  //  but will take different ammounts of time for different problems.
  //
  ContactNodeFaceRangeSearch::SearchType determine_range_search_type(ContactSearchData *search_data);
  //
  //  Overlapping sphere ghost takes as input a number of spheres owned by the current processor.
  //  It returns the spheres on other processors that need to be ghosted.
  //
  //  Input:
  //    int         num_sphere        The number of spheres owned by the current processor
  //
  //    Real*       sphere_data       Array of the sphere data.  That data layout is {c1_x, c1_y, c1_z, r1,
  //                                  c2_x, c2_y, c2_z, r2, etc.}
  // 
  //    MPI_Comm    mpi_communicator  Communicator for the processor group where at least each processor that
  //                                  owns a sphere is included in the communicator.
  //
  //    Int_Vector  ghost_indexes     The index numbers for the ghosts to send to another processor.  The index
  //                                  numbers refer to the indexes of the input particles.
  //
  //    Int_Vector  ghost_procs       The processor numbers to send ghosts to.  This array will be the same length
  //                                  as ghost_indexes.
  //
  int Overlapping_Sphere_Ghost(const int num_sphere, 
                               const Real* sphere_data,
		               MPI_Comm &mpi_communicator,
			       Int_Vector &ghost_indexes,
                               Int_Vector &ghost_procs);
  //
  //  Overlapping sphere search takes as input a number of spheres defined by center points and radii.  
  //  The spheres are divided into two groups.  Group1 and Group2.  the search finds all spheres in Group1
  //  one that overlap spheres in Group2.
  //
  //  Input:
  //    int         num_sphere1       The number of spheres in group1
  //
  //    Real*       sphere_data1      Array of the sphere data.  The data layout is {c1_x, c1_y, c1_z, r1, 
  //                                  c2_x, c2_y, c2_z, re, etc.}
  // 
  //    int         num_sphere2       The number of spheres in group2
  //
  //    Real*       sphere_data2      Array of the sphere data.  The data layout is {c1_x, c1_y, c1_z, r1,
  //                                  c2_x, c2_y, c2_z, r2, etc.}
  // 
  //  Output:
  //    Int_Vector  interaction_list  This array contains the interactions.  A value of '5' would indicate
  //                                  that the sphere in group1 owning this section of the list would interact with the
  //                                  5th sphere in group2.
  //
  //    Int_Vector  first_interaction  This array contains the locations of the first interaction in the 
  //                                   interaction list for each sphere of group1.
  //
  //    Int_Vector  last_interaction  This array contains the locations of the last interaction in the 
  //                                  interaction list for each sphere group1.
  //
  //  Returns:
  //    An integer error code.  Zero = success.  Non zero = some failure error code.
  //
  int Overlapping_Sphere_Search(const int num_sphere1, 
                                const Real* sphere_data1,
                                const int num_sphere2, 
                                const Real* sphere_data2,
			        Int_Vector& interaction_list,
                                Int_Vector& first_interaction,
                                Int_Vector& last_interaction);
  //-------------------------------------------------------------------------------------------------------------
  //


  //
  //  Sphere_Point_Ghost takes as input a number of spheres and points owned by the current processor.
  //  It returns the points that need to be ghosted to other processors
  //
  //  Input:
  //    int         num_sphere        The number of spheres owned by the current processor
  //
  //    Real*       sphere_data       Array of the sphere data.  That data layout is {c1_x, c1_y, c1_z, r1,
  //                                  c2_x, c2_y, c2_z, r2, etc.}
  //
  //    int         num_point         The number of points owned by the current processor
  //
  //    Real*       point_data        Array of point data.  The data layout is {p1_x, p1_y, p1_z,
  //                                  p2_x, p2_y, p2_z, etc.};
  // 
  //    MPI_Comm    mpi_communicator  Communicator for the processor group where at least each processor that
  //                                  owns a sphere is included in the communicator.
  //
  //    Int_Vector  ghost_indexes     The index numbers for the ghosts to send to another processor.  The index
  //                                  numbers refer to the indexes of the input points.
  //
  //    Int_Vector  ghost_procs       The processor numbers to send ghosts to.  This array will be the same length
  //                                  as ghost_indexes.
  //
  int Sphere_Point_Ghost(const int num_sphere, 
                         const Real* sphere_data,
                         const int num_point,
                         const Real* point_data,
		         MPI_Comm &mpi_communicator,
		         Int_Vector &ghost_indexes,
                         Int_Vector &ghost_procs);
  //
  //  Sphere_Point_Search takes as input a number of spheres defined by center points and radii.  
  //  And a series of points.  The search finds intersections between the points and the spheres
  //
  //  Input:
  //    int         num_sphere1       The number of spheres in group1
  //
  //    Real*       sphere_data1      Array of the sphere data.  The data layout is {c1_x, c1_y, c1_z, r1, 
  //                                  c2_x, c2_y, c2_z, re, etc.}
  // 
  //    int         num_sphere2       The number of spheres in group2
  //
  //    Real*       sphere_data2      Array of the sphere data.  The data layout is {c1_x, c1_y, c1_z, r1,
  //                                  c2_x, c2_y, c2_z, r2, etc.}
  // 
  //  Output:
  //    Int_Vector  interaction_list  This array contains the interactions.  A value of '5' would indicate
  //                                  that the sphere owning this section of the list would interact with the
  //                                  5th point passed in.
  //
  //    Int_Vector  first_interaction  This array contains the locations of the first interaction in the 
  //                                   interaction list for each sphere
  //
  //    Int_Vector  last_interaction  This array contains the locations of the last interaction in the 
  //                                  interaction list for each sphere.
  //
  //  Returns:
  //    An integer error code.  Zero = success.  Non zero = some failure error code.
  //
  int Sphere_Point_Search(const int num_sphere, 
                          const Real* sphere_data,
                          const int num_point, 
                          const Real* point_data,
		          Int_Vector& interaction_list,
                          Int_Vector& first_interaction,
                          Int_Vector& last_interaction);




  //
  //  BoxA_BoxB_Ghost takes as input a number of boxes (set A) and boxes (set B) owned by the current processor.
  //  It returns the set of boxes in set B that overlap the boxes in set A.
  //
  //  Input:
  //    int         num_boxA          The number of boxes in set A owned by the current processor
  //
  //    Real*       boxA_data         Array of the box A data.  That data layout is {a1_x_lo, a1_y_lo, a1_z_lo, 
  //                                                                                 a1_x_hi, a1_y_hi, a1_z_hi,
  //                                                                                 a2_x_lo, a2_y_lo, a2_z_lo, 
  //                                                                                 a2_x_hi, a2_y_hi, a2_z_hi}
  //
  //    int         num_boxB          The number of boxes in set B owned by the current processor
  //
  //    Real*       boxB_data         Array of box B data.  The data layout is {b1_x_lo, b1_y_lo, b1_z_lo, 
  //                                                                            b1_x_hi, b1_y_hi, b1_z_hi,
  //                                                                            b2_x_lo, b2_y_lo, b2_z_lo, 
  //                                                                            b2_x_hi, b2_y_hi, b2_z_hi}
  // 
  //    MPI_Comm    mpi_communicator  Communicator for the processor group where at least each processor that
  //                                  owns a sphere is included in the communicator.
  //
  //    Int_Vector  ghost_indexes     The index numbers for the ghosts to send to another processor.  The index
  //                                  numbers refer to the indexes of the input objects of set B.
  //
  //    Int_Vector  ghost_procs       The processor numbers to send ghosts to.  This array will be the same length
  //                                  as ghost_indexes.
  //
  int BoxA_BoxB_Ghost(const int num_boxA, 
                      const Real* boxA_data,
                      const int num_boxB,
                      const Real* boxB_data,
		      MPI_Comm &mpi_communicator,
		      Int_Vector &ghost_indexes,
                      Int_Vector &ghost_procs);
  //
  //  Sphere_Point_Search takes as input a number of spheres defined by center points and radii.  
  //  And a series of points.  The search finds intersections between the points and the spheres
  //
  //  Input:
  //    int         num_boxA          The number of boxes in set A owned by the current processor
  //
  //    Real*       boxA_data         Array of the box A data.  That data layout is {a1_x_lo, a1_y_lo, a1_z_lo, 
  //                                                                                 a1_x_hi, a1_y_hi, a1_z_hi,
  //                                                                                 a2_x_lo, a2_y_lo, a2_z_lo, 
  //                                                                                 a2_x_hi, a2_y_hi, a2_z_hi}
  //
  //    int         num_boxB          The number of boxes in set B owned by the current processor
  //
  //    Real*       boxB_data         Array of box B data.  The data layout is {b1_x_lo, b1_y_lo, b1_z_lo, 
  //                                                                            b1_x_hi, b1_y_hi, b1_z_hi,
  //                                                                            b2_x_lo, b2_y_lo, b2_z_lo, 
  //                                                                            b2_x_hi, b2_y_hi, b2_z_hi}
  //
  //  Output:
  //    Int_Vector  interaction_list  This array contains the interactions.  A value of '5' would indicate
  //                                  that the boxA owning this section of the list would interact with the
  //                                  5th boxB passed in.
  //
  //    Int_Vector  first_interaction  This array contains the locations of the first interaction in the 
  //                                   interaction list for each boxA
  //
  //    Int_Vector  last_interaction  This array contains the locations of the last interaction in the 
  //                                  interaction list for each boxA.
  //
  //  Returns:
  //    An integer error code.  Zero = success.  Non zero = some failure error code.
  //
  int BoxA_BoxB_Search(const int num_boxA, 
                       const Real* boxA_data,
                       const int num_boxB, 
                       const Real* boxB_data,
		       Int_Vector& interaction_list,
                       Int_Vector& first_interaction,
                       Int_Vector& last_interaction);


  //
  //-------------------------------------------------------------------------------------------------------------
  //
  //  ##### Helper routines. 

  //
  //  Sphere_Sphere_Overlap: Determine if two spheres overlap.
  //
  //  Input:
  //    Real *   sphere_data1   Array of data describing sphere 1, in c_x, c_y, c_z, r
  //
  //    Real *   sphere_data2   Array of data describing sphere 2, in c_x, c_y, c_z, r
  //
  //  Returns:
  //    True if the spheres overlap.  False otherwise.
  //
  bool Sphere_Sphere_Overlap(const Real* sphere_data1,
                             const Real* sphere_data2);
  bool Sphere_Point_Overlap(const Real* sphere_data,
                            const Real* point_data);

} // end namespace ACME


#endif // #ifdef ContactRangeSearch_h_
