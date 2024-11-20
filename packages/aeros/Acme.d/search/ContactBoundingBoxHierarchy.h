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

 
#ifndef ContactBoundingBoxHierarchy_h_
#define ContactBoundingBoxHierarchy_h_

#include "ContactBoundingBox.h"
#include "ContactVector.h"
#include <vector>
#include <map>

//
//  Class for general hierarchy of bounding boxes.  This class is used for fast overlap searches.  Note, the 
//  calculations for creating the hierarchy and searching the hierachy have been HEAVILY optimized for both
//  reduction of operation count as well as for cache efficency.
//
//  NKC NOTE, possible further optimizations include
//    2)  Implement a more efficent search when search for multple objects.  Currently each object is searched for
//        independiantly.  Might be able to get some better performance by taking into account search is searching
//        some set of objects against another set of objects.
//
class ObjectBoundingBoxHierarchy : public ContactBoundingBox{
 public:
  ObjectBoundingBoxHierarchy();
  virtual ~ObjectBoundingBoxHierarchy();
  //
  //  Create a hierarchy from a list of object bounding boxes.  The top level of the hierarchy will contain all boxes.  In the next level
  //  each leave of the tree will contain approximiatly half of the boxes.  The hierarchy continues down until the tree node contains only 
  //  a single box
  //
  static void create_hierarchy(ObjectBoundingBoxHierarchy *hierarchy_data, 
                               ObjectBoundingBox *boxes,
                               const int &num_boxes);
  //
  //  Routine to find overlap between the input box and all boxes in the hierarchy.  The object indexes that the input object overlaps will
  //  be stored in the obj_list array which must be pre-allocated.  The number of objects found is stored in the num_objs data item.
  //
  static void search_for_overlap_loop(const ObjectBoundingBoxHierarchy *const hierarchy_start_ptr,
                                const ContactBoundingBox &search_box,
                                int *obj_list,
                                int  &num_objs);

  static bool find_any_overlap_loop(const ObjectBoundingBoxHierarchy *const hierarchy_start_ptr,
                                    const ContactBoundingBox &search_box);

  static void search_for_overlap(const ObjectBoundingBoxHierarchy *const hierarchy_start_ptr,
                                 const ContactBoundingBox &search_box,
                                 ContactNode<Real> **Nodes,
                                 ACME::ContactNode_Vector &node_list);
  //
  //  Search for overlaps in a hybrid tree, tree one contains the roots nodes for sets of trees in root_tree
  //
  static void search_for_overlap_hybrid_tree(const ObjectBoundingBoxHierarchy *const tree1,
                                             const ObjectBoundingBoxHierarchy *const *const tree2_ptrs,
				             const ContactBoundingBox &search_box, 
					     ContactNode<Real> **Nodes, 
					     ACME::Int_Vector &node_keys,
                                             ACME::ContactNode_Vector &node_list,
                                             std::vector<bool> &valid_inter);
  //
  //  Perform a symmetric tree on itself search.  Can make various optimizations for this search such as, if A overlaps B,
  //  B must overlap A.
  //
  void search_for_overlap_recurse_sym(const ObjectBoundingBoxHierarchy &object_box,
                                      int *obj_list,
                                      int &num_objs,
                                      const ObjectBoundingBoxHierarchy *max_owned_ptr
                                      );

  //
  //  Extract the object number for the current tree node.  If the node is not a terminal node, return -1
  //
  inline int get_object_number() const {
    if(right_child_offset <= 0) {
      return -right_child_offset;
    } else {
      return -1;
    }
  };


  int get_right_child_offset() const {return right_child_offset;};
  inline void set_right_child_offset(const int right_child_offset_) {right_child_offset = right_child_offset_;};
  //
  //  Update the bounding boxes for all objects in the tree based on new leaf node bounding boxes
  //
  void update_hierarchy_recurse();
  //
  //  Search for overlap taking into account node object keys.  The tree overlap calculations uses both the spatial 
  //  overlap of boxes, and checks if the box can contain valid interactions based on the nodal entity keys.
  //
  static void search_for_overlap_key_aware(const ObjectBoundingBoxHierarchy *const tree1,
				           const ContactBoundingBox &search_box, 
					   ContactNode<Real> **Nodes, 
				           ACME::Int_Vector &node_keys,
                                           ACME::ContactNode_Vector &node_list,
                                           std::vector<bool> &valid_inter,
                                           int* node_key_array);
  //
  //  Set the correct node_key on all nodes
  //
  void calc_keys_recurse(int *hierarchy_node_keys, std::vector<int> &node_entity_keys);
 protected:
  //
  //  Use a loop control structure to create the bounding box hierarchy
  //
  static void create_hierarchy_loop(ObjectBoundingBoxHierarchy *const hierarchy_start_ptr,
                                    ObjectBoundingBox *const boxes, 
                                    ObjectBoundingBox *const scratch_boxes,
                                    const int num_boxes,
                                    int *const index_array0,
                                    int *const index_array1,
                                    int *const index_array2);
  //
  //  These static data items are used when constructing the tree.  They will be allocated and then deleted by the create_hierarchy method.
  //  They are made static to allow all the recurise tree creation routine to quickley access them
  //
  static int *index_array;
  static int *index_array0;
  static int *index_array1;
  static int *index_array2;
  //
  //  Right child offset stores one of two things.
  //    If the offset is <= 0 the current object is a terminal node of the tree.  The value is the negative of the 
  //  index back to the object that the terminal node represents.  If the value is positive it is the offset from
  //  the current object to the objects std::right child.  Note that the offset to the std::left child is always one.  Thus
  //  for a given object the std::left child can be found at this[1] and the std::right child at this[right_child_offset]
  //
  int right_child_offset;
};

#endif  // #ifdef ContactBoundingBoxHierarchy_h_
